# src/chd_atlas/validate/runner.py
"""Runs every validator and aggregates the result into one deterministic report."""

from __future__ import annotations

import os
from dataclasses import dataclass, field, replace
from pathlib import Path

from chd_atlas.corpus import load_curation, unexpected_curation_entries
from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.tables import (
    TABLE_SCHEMAS,
    mirror_paths,
    read_table,
    unexpected_mirror_entries,
    validate_table,
)
from chd_atlas.validate.burden import validate_burden, validate_burden_references
from chd_atlas.validate.ids import load_id_registry, validate_ids
from chd_atlas.validate.ontology import (
    OntologyRegistry,
    validate_labels,
    validate_terms,
)
from chd_atlas.validate.referential import (
    validate_mirror_references,
    validate_ptm_evidence_is_reachable,
    validate_references,
)
from chd_atlas.validate.scope import (
    scope_candidates,
    validate_curation_is_in_scope,
    validate_scope_attribution,
    validate_scope_terms,
)
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
        # Warnings are deliberately ignored, so a warnings-only report exits 0.
        # Two kinds of warning exist, and the reasoning differs for each.
        #
        # SKIP warnings (REF000, SRC000, ONT000) mean "a check did not run".
        # Ignoring one would let CI go green on a repository whose checks were
        # silently skipped, so each is emitted in the same branch as the error
        # that caused the skip — a failed corpus load, or a failed source
        # registry load — and a skip can therefore never be the only thing
        # wrong. A skip warning that could arrive alone must be raised to an
        # error. Pinned by
        # test_every_skip_warning_arrives_with_the_error_that_caused_it.
        #
        # GAP warnings (REF013) mean "a check ran and found curated data that
        # will not reach a reader". They arrive alone, deliberately: the record
        # is sound and the site is publishable, and what is missing is one
        # mirror cell. Blocking a whole deploy on that would be the wrong trade
        # for a resource curated incrementally. This is a real weakening — the
        # site can publish with a known gap — and it is accepted because the
        # alternative before REF013 existed was the same gap reported nowhere at
        # all. The safety they rest on is different from the skip warnings': the
        # check ran, so the report names the gene, the assertion and the file.
        # Pinned by test_a_gap_warning_is_reported_without_blocking_the_build.
        return self.error_count == 0

    def render(self) -> str:
        lines = [issue.format() for issue in self.issues]
        # Blank line so the summary reads as a summary rather than as one more
        # issue running straight on from the last.
        lines.append("")
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


def _mirrored_validity(
    root: Path,
) -> (
    tuple[
        dict[str, str],
        dict[str, set[str]],
        dict[str, set[str]],
        dict[str, frozenset[str]],
        dict[str, frozenset[str]],
    ]
    | None
):
    """Disease labels and gene/disease cross-links, from every validity mirror.

    Returns three mappings the scope validators need, built once from both
    mirror tables so each validator stays a pure check over already-loaded
    data:

    - a MONDO term -> the label its authority (ClinGen or GenCC) published
    - a MONDO term -> every HGNC id curated under it
    - an HGNC id -> every MONDO term it is curated under
    - a MONDO term -> every ClinGen expert panel that curated a record naming it
    - a MONDO term -> every GenCC submitter that did

    The last two are kept APART rather than merged into one set of names, and
    that separation is the whole of SCP005's strength. Merged, the check can
    only ask "did some authority use this term", so a curator may take a term
    whose sole warrant is a commercial laboratory's GenCC submission, label it
    `admitted_by: clingen_chd_panel`, and pass -- publishing on 24 pages that
    ClinGen's chartered panel scoped a disease it never curated. Found by
    adversarial review of #30, on the guard added in the same commit.

    None is deliberately distinct from an empty result, for exactly the reason
    `_known_genes` returns None: checking every scope term against an empty
    mapping would report SCP001 once per term — 68 failures naming the
    symptom, none naming the cause — so a mirror that could not be read at all
    is reported once, as SCP000, and nothing else runs. Unlike `_known_genes`,
    which needs one file, this needs *neither* mirror to be readable before
    giving up — one working mirror still gives real, if partial, coverage.
    """
    labels: dict[str, str] = {}
    genes_by_disease: dict[str, set[str]] = {}
    diseases_by_gene: dict[str, set[str]] = {}
    clingen_panels: dict[str, set[str]] = {}
    gencc_submitters: dict[str, set[str]] = {}
    read_any = False
    # GenCC first, so ClinGen's label overwrites it below: ClinGen is the
    # primary source for gene-disease validity and GenCC harmonises rather
    # than authors, so preferring GenCC's label would attribute it to the
    # resource that did not write it.
    for schema_name, filename in (
        ("gencc_submissions", "gencc_submissions.tsv"),
        ("clingen_validity", "clingen_gene_validity.tsv"),
    ):
        path = root / "mirrors" / filename
        frame, _ = read_table(path, TABLE_SCHEMAS[schema_name])
        # All three columns this loop selects, not only `disease`: an upstream
        # rename can drop `disease_label` or `gene` just as easily, and
        # `frame.select` raises `ColumnNotFoundError` on any one of them
        # missing — `validate_table` (called separately, over the same file)
        # already reports that as TBL001, so nothing is lost by skipping the
        # mirror here rather than crashing the whole validation run on it.
        if frame is None or not {"disease", "disease_label", "gene"}.issubset(frame.columns):
            continue
        read_any = True
        # The column naming the body that made each call: ClinGen records carry
        # the expert panel, GenCC rows the submitter. Selected separately from
        # the three above because either mirror may lack it without the other
        # checks being affected -- `validate_table` reports a missing column as
        # TBL001 over the same file, so skipping it here loses nothing.
        clingen = schema_name == "clingen_validity"
        who = "gcep" if clingen else "submitter"
        into = clingen_panels if clingen else gencc_submitters
        curators = frame[who].to_list() if who in frame.columns else [None] * frame.height
        rows = frame.select(["disease", "disease_label", "gene"]).iter_rows()
        for (disease, label, gene), curator in zip(rows, curators, strict=True):
            if disease is None:
                continue
            if label is not None:
                labels[str(disease)] = str(label)
            if gene is not None:
                genes_by_disease.setdefault(str(disease), set()).add(str(gene))
                diseases_by_gene.setdefault(str(gene), set()).add(str(disease))
            if curator is not None:
                into.setdefault(str(disease), set()).add(str(curator))
    if not read_any:
        return None
    return (
        labels,
        genes_by_disease,
        diseases_by_gene,
        {term: frozenset(names) for term, names in clingen_panels.items()},
        {term: frozenset(names) for term, names in gencc_submitters.items()},
    )


def _relative_to_root(issues: list[ValidationIssue], root: Path) -> list[ValidationIssue]:
    """Rewrite absolute locations as repo-relative.

    Locations are embedded verbatim in the rendered report, so absolute paths
    would make output differ between a curator's machine and CI, and would sort
    a hardcoded relative literal away from the same file's other issues.

    The separator comes from `os.sep` rather than a literal `/` because every
    location is built as `str(path)`, which renders `C:\\repo\\mirrors\\genes.tsv`
    on Windows. A hardcoded `/` would match no prefix there, leave every location
    absolute, and silently void the guarantee this function exists to provide.
    """
    prefix = f"{root}{os.sep}"
    return [replace(issue, location=issue.location.removeprefix(prefix)) for issue in issues]


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
    # Also unconditional, and deliberately outside the corpus branch below: the
    # burden table's cross-column rules read nothing but the table itself, so a
    # corpus that failed to load must not silently take these checks with it.
    issues.extend(validate_burden(root))
    # Same reasoning for the interpretive tree, which holds the curator
    # judgement the atlas exists to record.
    issues.extend(unexpected_curation_entries(root))

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
                    "gene registry is missing or unreadable, so gene references cannot be checked",
                )
            )
        # Without the phenotype vocabulary, REF007/REF009/REF010 silently no-op:
        # `lesion_group_of` is empty so every lookup misses. An error, not a
        # warning, because `ok` ignores warnings and a check stopped running.
        if corpus.assertions and not corpus.phenotypes:
            issues.append(
                ValidationIssue(
                    "CUR002",
                    Severity.ERROR,
                    str(root / "curation" / "phenotypes.yaml"),
                    "phenotype vocabulary is missing or empty, so lesion-group "
                    "consistency checks cannot run",
                )
            )
        # Otherwise every citation reports REF002 and nothing names the cause.
        if (corpus.assertions or corpus.functional) and not corpus.publications:
            issues.append(
                ValidationIssue(
                    "CUR003",
                    Severity.ERROR,
                    str(root / "curation" / "publications.yaml"),
                    "publication registry is missing or empty, so every citation "
                    "will be reported as dangling",
                )
            )
        issues.extend(validate_references(corpus, known_genes=known_genes))
        issues.extend(validate_mirror_references(root, corpus))
        # Deliberately not folded into the `validate_burden` call above: these
        # need registries a failed corpus load empties, and the comparator rules
        # must keep running when that happens.
        #
        # **Each registry is passed as None when it is empty, not as an empty
        # set**, and this is the load-bearing part. The first version passed the
        # sets unconditionally and broke the guarantee in §1 of CLAUDE.md that a
        # missing registry must not report as hundreds of dangling references.
        # Measured 2026-08-05 against the committed corpus:
        #
        #   deleting curation/cohorts.yaml   -> 4 BUR009 errors, none naming the cause
        #   header-only mirrors/genes.tsv    -> 146 errors (1 REF001 + 145 BUR011),
        #                                       against 1 error before this table existed
        #
        # `_known_genes` returns None only for a missing or unreadable file; a
        # header-only mirror reads fine and yields an empty set, which is why the
        # emptiness test lives here rather than there.
        burden_path = root / "mirrors" / "burden.tsv"
        cohort_ids = {str(cohort.id) for cohort in corpus.cohorts}
        study_ids = {str(publication.id) for publication in corpus.publications}
        if burden_path.is_file() and not cohort_ids:
            issues.append(
                ValidationIssue(
                    "CUR004",
                    Severity.ERROR,
                    str(root / "curation" / "cohorts.yaml"),
                    "cohort registry is missing or empty, so every cohort a burden "
                    "row cites will be reported as dangling",
                )
            )
        # CUR003 above fires only when an assertion or a functional record needs
        # the registry. A burden row cites a publication too, so the same cause
        # needs the same report when burden is the only thing citing one.
        if burden_path.is_file() and not study_ids and not (corpus.assertions or corpus.functional):
            issues.append(
                ValidationIssue(
                    "CUR003",
                    Severity.ERROR,
                    str(root / "curation" / "publications.yaml"),
                    "publication registry is missing or empty, so every citation "
                    "will be reported as dangling",
                )
            )
        issues.extend(
            validate_burden_references(
                root,
                known_cohorts=cohort_ids or None,
                known_genes=known_genes or None,
                known_studies=study_ids or None,
            )
        )
        # Inside the same branch as the two above, and for the same reason: it
        # reads the gene registry, so on a corpus that failed to load it would
        # report a missing accession for every gene at once.
        issues.extend(validate_ptm_evidence_is_reachable(root, corpus))

        # Same branch, same reasoning again: a corpus that failed to load empties
        # `corpus.chd_scope`, and every scope term would then look absent.
        scope_location = str(root / "curation" / "chd_scope.yaml")
        mirrored = _mirrored_validity(root)
        if mirrored is None:
            # SCP000 alone is a WARNING, and `ValidationReport.ok` ignores
            # warnings — so without an accompanying ERROR, a repository missing
            # both validity mirrors would validate "clean" while every scope
            # check silently could not run, and `build_site` would publish it
            # anyway. Same shape as `TBL008` for a missing gene registry: name
            # the cause as an error, located at the mirrors directory because
            # two files are implicated and neither alone explains it.
            issues.append(
                ValidationIssue(
                    "TBL012",
                    Severity.ERROR,
                    str(root / "mirrors"),
                    "the validity mirrors are missing or unreadable, so no gene "
                    "can be placed in scope",
                )
            )
            issues.extend(validate_scope_terms(corpus.chd_scope, None, scope_location))
            issues.extend(validate_scope_attribution(corpus.chd_scope, None, None, scope_location))
        else:
            labels, genes_by_disease, diseases_by_gene, panels, submitters = mirrored
            issues.extend(validate_scope_terms(corpus.chd_scope, labels, scope_location))
            issues.extend(
                validate_scope_attribution(corpus.chd_scope, panels, submitters, scope_location)
            )
            scope_terms = {str(entry.id) for entry in corpus.chd_scope}
            in_scope_genes = {
                gene for term in scope_terms for gene in genes_by_disease.get(term, set())
            }
            issues.extend(
                scope_candidates(
                    scope_terms, labels, diseases_by_gene, in_scope_genes, scope_location
                )
            )
            # `Corpus` does not record which file each assertion was loaded from
            # (see `corpus.py::Corpus`), so every curated gene is attributed to
            # the assertions directory rather than a specific file within it —
            # coarser than a curator fixing SCP004 would want, but not a
            # precision this module has and is pretending not to.
            curated_genes = {
                str(assertion.gene): str(root / "curation" / "assertions")
                for assertion in corpus.assertions
            }
            issues.extend(validate_curation_is_in_scope(curated_genes, in_scope_genes))

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
    # ...and likewise when a pinned release itself failed to load: ONT003 reads
    # "no pinned ontology for prefix 'HP'", which is false when the prefix IS
    # pinned and the file is simply unreadable. A curator would go add a pin
    # that already exists. The ONT004 naming the real cause is already reported.
    if source_issues or ontologies.load_issues:
        cause = (
            "the source registry did not load, so no releases are pinned"
            if source_issues
            else "a pinned ontology release failed to load"
        )
        issues.append(
            ValidationIssue(
                "ONT000",
                Severity.WARNING,
                str(root / "mirrors" / "sources.yaml"),
                f"skipped ontology checks: {cause}",
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
        phenotypes_file = str(root / "curation" / "phenotypes.yaml")
        for term in corpus.phenotypes:
            issues.extend(validate_terms([term.id], ontologies, phenotypes_file))
        # The file declares its labels transcribed from the pinned release, and
        # those labels are what the atlas displays. Pinning makes the claim
        # checkable, so check it.
        issues.extend(
            validate_labels(
                [(term.id, term.label) for term in corpus.phenotypes],
                ontologies,
                phenotypes_file,
            )
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
