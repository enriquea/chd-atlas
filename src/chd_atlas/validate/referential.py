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


def validate_references(corpus: Corpus, known_genes: set[str]) -> list[ValidationIssue]:
    """Check that every reference in the corpus resolves to something real.

    ``known_genes`` comes from ``mirrors/genes.tsv``; it is passed in rather
    than read here so this function stays a pure check over loaded data.
    """
    issues: list[ValidationIssue] = []
    known_publications = {publication.id for publication in corpus.publications}
    known_functional = {record.id for record in corpus.functional}
    known_datasets = {dataset.id for dataset in corpus.datasets}
    lesion_group_of = {term.id: term.lesion_group for term in corpus.phenotypes}

    def error(code: str, location: str, message: str) -> None:
        issues.append(ValidationIssue(code, Severity.ERROR, location, message))

    for assertion in corpus.assertions:
        location = f"assertion {assertion.id}"

        if assertion.gene not in known_genes:
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

        for index, evidence in enumerate(assertion.evidence):
            evidence_location = f"{location} evidence[{index}]"

            if evidence.publication not in known_publications:
                error(
                    "REF002",
                    evidence_location,
                    f"publication {evidence.publication} is not in curation/publications.yaml",
                )
            if (
                evidence.functional_evidence is not None
                and evidence.functional_evidence not in known_functional
            ):
                error(
                    "REF003",
                    evidence_location,
                    f"functional evidence {evidence.functional_evidence} does not exist",
                )
            if evidence.dataset is not None and evidence.dataset not in known_datasets:
                error(
                    "REF004",
                    evidence_location,
                    f"dataset {evidence.dataset} is not in curation/datasets/",
                )

    for record in corpus.functional:
        location = f"functional evidence {record.id}"
        if record.gene not in known_genes:
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
        pairs = {
            (row[0], row[1] if has_contrast else None)
            for row in frame.select(selected).rows()
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
