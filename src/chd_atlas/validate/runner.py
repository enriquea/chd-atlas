# src/chd_atlas/validate/runner.py
"""Runs every validator and aggregates the result into one deterministic report."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path

from chd_atlas.corpus import load_curation
from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.tables import (
    TABLE_SCHEMAS,
    mirror_paths,
    read_table,
    unexpected_mirror_entries,
    validate_table,
)
from chd_atlas.validate.ids import load_id_registry, validate_ids
from chd_atlas.validate.ontology import OntologyRegistry, validate_terms
from chd_atlas.validate.referential import validate_mirror_references, validate_references
from chd_atlas.validate.sort_order import validate_sort_order
from chd_atlas.validate.sources import load_sources, validate_source_references


@dataclass
class ValidationReport:
    issues: list[ValidationIssue] = field(default_factory=list)

    def __post_init__(self) -> None:
        # Sorted by location first so everything wrong with one file reads together.
        # ValidationIssue's own ordering is code-first, which groups by error type and
        # scatters a single file's problems across the report — fine as a stable
        # default, wrong for a human fixing one file at a time.
        self.issues = sorted(
            set(self.issues),
            key=lambda issue: (issue.location, issue.code, issue.severity, issue.message),
        )

    @property
    def error_count(self) -> int:
        return sum(1 for issue in self.issues if issue.severity is Severity.ERROR)

    @property
    def warning_count(self) -> int:
        return sum(1 for issue in self.issues if issue.severity is Severity.WARNING)

    @property
    def ok(self) -> bool:
        return self.error_count == 0

    def render(self) -> str:
        lines = [issue.format() for issue in self.issues]
        lines.append(f"{self.error_count} error(s), {self.warning_count} warning(s)")
        return "\n".join(lines)


def _used_sources(root: Path) -> dict[str, set[str]]:
    """Collect the distinct ``source`` values each mirror table relies on."""
    used: dict[str, set[str]] = {}
    for path, schema_name in mirror_paths(root):
        schema = TABLE_SCHEMAS[schema_name]
        if not any(column.name == "source" for column in schema.columns):
            continue
        frame, _ = read_table(path, schema)
        if frame is None or "source" not in frame.columns:
            continue
        values = {str(value) for value in frame["source"].to_list() if value is not None}
        if values:
            used[str(path)] = values
    return used


def _known_genes(root: Path) -> set[str] | None:
    """The gene registry, or None when it could not be read.

    None is deliberately distinct from an empty set. An empty registry makes
    every assertion and every functional record report REF001, burying the one
    real problem — the registry itself — under one error per record, with
    nothing anywhere naming the cause.
    """
    path = root / "mirrors" / "genes.tsv"
    if not path.is_file():
        return None
    frame, _ = read_table(path, TABLE_SCHEMAS["genes"])
    if frame is None or "hgnc_id" not in frame.columns:
        return None
    return {value for value in frame["hgnc_id"].to_list() if value is not None}


def _relative_to_root(issues: list[ValidationIssue], root: Path) -> list[ValidationIssue]:
    """Rewrite absolute locations as repo-relative.

    Locations are embedded verbatim in the rendered report, so absolute paths
    would make output differ between a curator's machine and CI, and would sort
    a hardcoded relative literal away from the same file's other issues.
    """
    prefix = f"{root}/"
    return [
        replace(issue, location=issue.location.removeprefix(prefix))
        for issue in issues
    ]


def validate_repository(root: Path) -> ValidationReport:
    issues: list[ValidationIssue] = []

    corpus, corpus_issues = load_curation(root)
    issues.extend(corpus_issues)

    registry, source_issues = load_sources(root)
    issues.extend(source_issues)

    ontologies = OntologyRegistry.from_files(root, registry.ontology_files())
    issues.extend(ontologies.load_issues)

    for path, schema_name in mirror_paths(root):
        schema = TABLE_SCHEMAS[schema_name]
        issues.extend(validate_table(path, schema))
        issues.extend(validate_sort_order(path, schema))
    # Unconditional: the loop above iterates what the schemas expect, so an entry
    # nothing claims — a shard directory lost to a typo — is invisible to it.
    issues.extend(unexpected_mirror_entries(root))

    # Referential checks over a corpus that failed to load would compare against a
    # knowingly incomplete set. The file wrappers validate their whole record list
    # as one model, so a single malformed publication empties `corpus.publications`
    # entirely — and every assertion citing a perfectly valid PMID would then be
    # reported as a dangling reference. Hundreds of false positives would bury the
    # one real error. Report the skip instead.
    if corpus_issues:
        issues.append(
            ValidationIssue(
                "REF000",
                Severity.WARNING,
                str(root / "curation"),
                f"skipped referential checks: {len(corpus_issues)} record(s) failed to load",
            )
        )
    else:
        known_genes = _known_genes(root)
        if known_genes is None:
            issues.append(
                ValidationIssue(
                    "TBL008",
                    Severity.ERROR,
                    str(root / "mirrors" / "genes.tsv"),
                    "gene registry is missing or unreadable, so gene references "
                    "cannot be checked",
                )
            )
        issues.extend(validate_references(corpus, known_genes=known_genes))
        issues.extend(validate_mirror_references(root, corpus))

    # Same reasoning as the referential skip above: on a failed registry load
    # `registry` is empty, so every source every mirror table uses would report
    # SRC003 and bury the one real error.
    if source_issues:
        issues.append(
            ValidationIssue(
                "SRC000",
                Severity.WARNING,
                str(root / "mirrors" / "sources.yaml"),
                "skipped source-reference checks: the source registry did not load",
            )
        )
    else:
        issues.extend(validate_source_references(_used_sources(root), registry))

    # The ontology pins live in sources.yaml, so a registry that failed to load
    # pins nothing and every term reports ONT003 "no pinned ontology" — one per
    # term, burying the single SRC001 that caused them. Same cascade as above.
    if source_issues:
        issues.append(
            ValidationIssue(
                "ONT000",
                Severity.WARNING,
                str(root / "mirrors" / "sources.yaml"),
                "skipped ontology checks: the source registry did not load, "
                "so no releases are pinned",
            )
        )
    else:
        for assertion in corpus.assertions:
            issues.extend(
                validate_terms(
                    [*assertion.phenotypes, *assertion.extracardiac_features],
                    ontologies,
                    f"assertion {assertion.id}",
                )
            )
        for record in corpus.functional:
            issues.extend(
                validate_terms(
                    [record.organism, *record.cardiac_phenotype],
                    ontologies,
                    f"functional evidence {record.id}",
                )
            )
        for term in corpus.phenotypes:
            issues.extend(
                validate_terms([term.id], ontologies, str(root / "curation" / "phenotypes.yaml"))
            )

    id_registry, id_issues = load_id_registry(root)
    issues.extend(id_issues)
    # None means the counter file is missing or malformed, which ID001/ID004
    # already named. Checking against an empty registry would additionally make
    # every atlas ID exceed its ceiling, burying that one error under one per record.
    if id_registry is not None:
        issues.extend(
            validate_ids(
                [assertion.id for assertion in corpus.assertions]
                + [record.id for record in corpus.functional],
                id_registry,
            )
        )

    # Relativised before the report is built, so `__post_init__` sorts the
    # locations the reader will actually see.
    return ValidationReport(issues=_relative_to_root(issues, root))
