# src/chd_atlas/validate/referential.py
"""Cross-entity reference integrity across the curation corpus and mirror tables."""

from __future__ import annotations

from pathlib import Path
from typing import Final

from chd_atlas.corpus import Corpus
from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.tables import TABLE_SCHEMAS, mirror_paths, read_table

# Mirror tables whose rows point back at a curated dataset (and usually a contrast).
_DATASET_LINKED_TABLES: Final[frozenset[str]] = frozenset(
    {"expression", "profiles", "proteomics", "phospho"}
)


def validate_references(
    corpus: Corpus, known_genes: set[str] | None
) -> list[ValidationIssue]:
    """Check that every reference in the corpus resolves to something real.

    ``known_genes`` comes from ``mirrors/genes.tsv``; it is passed in rather
    than read here so this function stays a pure check over loaded data.

    ``None`` means that registry could not be read, so the gene checks are
    skipped. An empty set is indistinguishable from a missing file, and a
    missing file would otherwise report one REF001 per assertion and per
    functional record — hundreds of failures naming the cause nowhere. The
    caller reports the unreadable registry once instead.
    """
    issues: list[ValidationIssue] = []
    known_publications = {publication.id for publication in corpus.publications}
    functional_by_id = {record.id: record for record in corpus.functional}
    known_datasets = {dataset.id for dataset in corpus.datasets}
    lesion_group_of = {term.id: term.lesion_group for term in corpus.phenotypes}

    def error(code: str, location: str, message: str) -> None:
        issues.append(ValidationIssue(code, Severity.ERROR, location, message))

    for assertion in corpus.assertions:
        location = f"assertion {assertion.id}"

        if known_genes is not None and assertion.gene not in known_genes:
            error("REF001", location, f"gene {assertion.gene} is not in mirrors/genes.tsv")

        # An assertion pairing an atrial septal defect with lesion_groups [conotruncal]
        # is internally inconsistent. The mapping lives in curation/phenotypes.yaml, so
        # the single-record model cannot see it — it has to be checked here.
        for phenotype in assertion.phenotypes:
            expected = lesion_group_of.get(phenotype)
            if expected is not None and expected not in assertion.lesion_groups:
                declared = sorted(group.value for group in assertion.lesion_groups)
                error(
                    "REF007",
                    location,
                    f"phenotype {phenotype} belongs to lesion group "
                    f"'{expected.value}', which is not among {declared}",
                )

        # The converse of REF007: a declared lesion group that no phenotype accounts for.
        # It files the assertion under a browse facet nothing in the record supports.
        # Sound only when every phenotype is mapped — an unmapped one could legitimately
        # justify the extra group, which is why REF007 exempts unmapped phenotypes too.
        if all(phenotype in lesion_group_of for phenotype in assertion.phenotypes):
            justified = {lesion_group_of[p] for p in assertion.phenotypes}
            unjustified = sorted(
                group.value for group in assertion.lesion_groups if group not in justified
            )
            if unjustified:
                error(
                    "REF009",
                    location,
                    f"lesion groups {unjustified} are not justified by any declared phenotype",
                )

        # curation/phenotypes.yaml is the atlas's own register of cardiac lesions, so a
        # term listed there cannot also be what makes the assertion syndromic.
        for feature in assertion.extracardiac_features:
            feature_group = lesion_group_of.get(feature)
            if feature_group is not None:
                error(
                    "REF010",
                    location,
                    f"extracardiac feature {feature} is a cardiac lesion in "
                    f"group '{feature_group.value}'",
                )

        for index, evidence in enumerate(assertion.evidence):
            evidence_location = f"{location} evidence[{index}]"

            if evidence.publication not in known_publications:
                error(
                    "REF002",
                    evidence_location,
                    f"publication {evidence.publication} is not in curation/publications.yaml",
                )
            if evidence.functional_evidence is not None:
                cited = functional_by_id.get(evidence.functional_evidence)
                if cited is None:
                    error(
                        "REF003",
                        evidence_location,
                        f"functional evidence {evidence.functional_evidence} does not exist",
                    )
                # Functional evidence feeds ClinGen-style scoring, so a record about
                # another gene can carry an assertion to 'definitive' on model data
                # that says nothing about the gene being asserted.
                elif cited.gene != assertion.gene:
                    error(
                        "REF008",
                        evidence_location,
                        f"functional evidence {evidence.functional_evidence} is about "
                        f"gene {cited.gene}, not {assertion.gene}",
                    )
            if evidence.dataset is not None and evidence.dataset not in known_datasets:
                error(
                    "REF004",
                    evidence_location,
                    f"dataset {evidence.dataset} is not in curation/datasets/",
                )

    for record in corpus.functional:
        location = f"functional evidence {record.id}"
        if known_genes is not None and record.gene not in known_genes:
            error("REF001", location, f"gene {record.gene} is not in mirrors/genes.tsv")
        if record.publication not in known_publications:
            error(
                "REF002",
                location,
                f"publication {record.publication} is not in curation/publications.yaml",
            )

    for entry in corpus.featured:
        if entry.publication not in known_publications:
            error(
                "REF002",
                f"featured entry {entry.order}",
                f"publication {entry.publication} is not in curation/publications.yaml",
            )

    for dataset in corpus.datasets:
        if dataset.publication is not None and dataset.publication not in known_publications:
            error(
                "REF002",
                f"dataset {dataset.id}",
                f"publication {dataset.publication} is not in curation/publications.yaml",
            )

    return issues


def validate_mirror_references(root: Path, corpus: Corpus) -> list[ValidationIssue]:
    """Check that mirror rows point at datasets and contrasts that actually exist.

    Distinct (dataset, contrast) pairs are checked once rather than per row, so a
    100,000-row table with one bad contrast yields one issue, not 100,000.
    """
    known_datasets = {dataset.id for dataset in corpus.datasets}
    known_contrasts = {
        (dataset.id, contrast.id)
        for dataset in corpus.datasets
        for contrast in dataset.contrasts
    }
    issues: list[ValidationIssue] = []

    for path, schema_name in mirror_paths(root):
        if schema_name not in _DATASET_LINKED_TABLES:
            continue

        frame, _ = read_table(path, TABLE_SCHEMAS[schema_name])
        if frame is None or "dataset" not in frame.columns:
            continue

        has_contrast = "contrast" in frame.columns
        selected = ["dataset", "contrast"] if has_contrast else ["dataset"]
        # A null dataset cell is a schema violation `validate_table` already reports as
        # TBL003 — both columns are non-nullable in every linked schema. Carrying it
        # here would add "dataset None is not in curation/datasets/", naming a dataset
        # that was never written down. A null contrast is skipped below, by the same
        # guard that skips tables having no contrast column at all.
        pairs = {
            (row[0], row[1] if has_contrast else None)
            for row in frame.select(selected).rows()
            if row[0] is not None
        }

        for dataset, contrast in sorted(pairs, key=lambda pair: (str(pair[0]), str(pair[1]))):
            if dataset not in known_datasets:
                issues.append(
                    ValidationIssue(
                        "REF005",
                        Severity.ERROR,
                        str(path),
                        f"dataset {dataset} is not in curation/datasets/",
                    )
                )
            elif contrast is not None and (dataset, contrast) not in known_contrasts:
                issues.append(
                    ValidationIssue(
                        "REF006",
                        Severity.ERROR,
                        str(path),
                        f"contrast '{contrast}' is not defined by dataset {dataset}",
                    )
                )

    return issues
