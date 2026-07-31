# tests/unit/test_referential.py
from datetime import date
from pathlib import Path

from chd_atlas.corpus import Corpus
from chd_atlas.models.assertion import (
    Evidence,
    GeneDiseaseAssertion,
    SupplementaryLocator,
)
from chd_atlas.models.functional import FunctionalEvidence
from chd_atlas.models.literature import FeaturedManuscript, PhenotypeTerm, Publication
from chd_atlas.validate.referential import (
    validate_mirror_references,
    validate_references,
)


def _publication(pmid: str = "PMID:8988165") -> Publication:
    return Publication(
        id=pmid,
        title="T",
        journal="J",
        year=1997,
        authors=["A"],
        study_type="case_control",
    )


def _evidence(**overrides: object) -> Evidence:
    payload: dict[str, object] = {
        "publication": "PMID:8988165",
        "evidence_class": "genetic_case",
        "locator": SupplementaryLocator(file="T.xlsx"),
        "summary": "s",
        "strength": "strong",
    }
    payload.update(overrides)
    return Evidence.model_validate(payload)


def _assertion(**overrides: object) -> GeneDiseaseAssertion:
    payload: dict[str, object] = {
        "id": "CHDA:AST:0000001",
        "gene": "HGNC:11604",
        "phenotypes": ["HP:0001631"],
        "lesion_groups": ["septal"],
        "classification": "definitive",
        "inheritance": ["AD"],
        "mechanism": "haploinsufficiency",
        "syndromic": "both",
        "evidence": [_evidence()],
        "source_tier": "own_curation",
        "curator": "c",
        "curated_on": date(2026, 7, 1),
        "last_reviewed": date(2026, 7, 1),
    }
    payload.update(overrides)
    return GeneDiseaseAssertion.model_validate(payload)


def _corpus(**overrides: object) -> Corpus:
    payload: dict[str, object] = {
        "root": Path("."),
        "assertions": (_assertion(),),
        "publications": (_publication(),),
    }
    payload.update(overrides)
    return Corpus(**payload)  # type: ignore[arg-type]


def test_accepts_a_referentially_intact_corpus() -> None:
    assert validate_references(_corpus(), known_genes={"HGNC:11604"}) == []


def test_reports_an_assertion_gene_absent_from_the_gene_registry() -> None:
    issues = validate_references(_corpus(), known_genes=set())
    assert [i.code for i in issues] == ["REF001"]
    assert "HGNC:11604" in issues[0].message


def test_reports_evidence_citing_an_unknown_publication() -> None:
    corpus = _corpus(assertions=(_assertion(evidence=[_evidence(publication="PMID:1")]),))
    issues = validate_references(corpus, known_genes={"HGNC:11604"})
    assert [i.code for i in issues] == ["REF002"]


def test_reports_an_unresolvable_functional_evidence_reference() -> None:
    evidence = _evidence(
        evidence_class="functional_model", functional_evidence="CHDA:FUN:0000009"
    )
    corpus = _corpus(assertions=(_assertion(evidence=[evidence]),))
    issues = validate_references(corpus, known_genes={"HGNC:11604"})
    assert [i.code for i in issues] == ["REF003"]


def test_resolves_a_present_functional_evidence_reference() -> None:
    functional = FunctionalEvidence.model_validate(
        {
            "id": "CHDA:FUN:0000009",
            "gene": "HGNC:11604",
            "organism": "NCBITaxon:10090",
            "perturbation": "knockout",
            "zygosity": "heterozygous",
            "cardiac_phenotype": ["MP:0010402"],
            "phenocopies_human": "partial",
            "rescue_outcome": "not_attempted",
            "publication": "PMID:8988165",
        }
    )
    evidence = _evidence(
        evidence_class="functional_model", functional_evidence="CHDA:FUN:0000009"
    )
    corpus = _corpus(
        assertions=(_assertion(evidence=[evidence]),), functional=(functional,)
    )
    assert validate_references(corpus, known_genes={"HGNC:11604"}) == []


def test_reports_a_featured_manuscript_citing_an_unknown_publication() -> None:
    featured = FeaturedManuscript(
        publication="PMID:999", order=1, blurb="b", topic="genomics"
    )
    issues = validate_references(_corpus(featured=(featured,)), known_genes={"HGNC:11604"})
    assert [i.code for i in issues] == ["REF002"]


def test_reports_evidence_citing_an_unknown_dataset() -> None:
    corpus = _corpus(assertions=(_assertion(evidence=[_evidence(dataset="PXD999999")]),))
    issues = validate_references(corpus, known_genes={"HGNC:11604"})
    assert [i.code for i in issues] == ["REF004"]


def _septal_term() -> PhenotypeTerm:
    return PhenotypeTerm.model_validate(
        {"id": "HP:0001631", "label": "Atrial septal defect", "lesion_group": "septal"}
    )


def test_phenotype_consistent_with_declared_lesion_group_passes() -> None:
    corpus = _corpus(phenotypes=(_septal_term(),))
    assert validate_references(corpus, known_genes={"HGNC:11604"}) == []


def test_reports_a_phenotype_contradicting_its_lesion_group() -> None:
    corpus = _corpus(
        assertions=(_assertion(lesion_groups=["conotruncal"]),),
        phenotypes=(_septal_term(),),
    )
    issues = validate_references(corpus, known_genes={"HGNC:11604"})
    assert [i.code for i in issues] == ["REF007"]
    assert "septal" in issues[0].message


def test_unmapped_phenotype_is_not_reported_as_inconsistent() -> None:
    """A phenotype absent from curation/phenotypes.yaml has no group to contradict."""
    corpus = _corpus(phenotypes=())
    assert validate_references(corpus, known_genes={"HGNC:11604"}) == []


DATASET_YAML = """\
id: PXD012345
archive: pride
technology: phosphoproteomics
tissue: right ventricle
developmental_stage: infant
organism: NCBITaxon:9606
n_samples: 24
licence: CC-BY-4.0
contrasts:
  - id: tof_vs_control
    description: TOF right ventricle versus control
    case_group: TOF RV
    control_group: Control RV
    statistical_method: limma-voom
    software: limma 3.58.1
"""

EXPRESSION_TSV_HEADER = (
    "dataset\tcontrast\tgene\tlog2fc\tpvalue\tfdr\tdirection\t"
    "n_case\tn_control\ttissue\tstage\n"
)


def _mirror_repo(tmp_path: Path, row: str) -> Path:
    from chd_atlas.corpus import load_curation

    (tmp_path / "curation" / "datasets").mkdir(parents=True)
    (tmp_path / "curation" / "datasets" / "PXD012345.yaml").write_text(DATASET_YAML)
    (tmp_path / "curation" / "publications.yaml").write_text(
        "publications:\n"
        "  - id: PMID:8988165\n"
        "    title: T\n"
        "    journal: J\n"
        "    year: 1997\n"
        "    authors: [A]\n"
        "    study_type: case_control\n"
    )
    (tmp_path / "mirrors" / "expression").mkdir(parents=True)
    (tmp_path / "mirrors" / "expression" / "PXD012345.tsv").write_text(
        EXPRESSION_TSV_HEADER + row
    )
    load_curation(tmp_path)
    return tmp_path


def test_mirror_rows_referencing_a_known_contrast_pass(tmp_path: Path) -> None:
    from chd_atlas.corpus import load_curation

    root = _mirror_repo(
        tmp_path,
        "PXD012345\ttof_vs_control\tHGNC:11604\t1.2\t0.001\t0.01\tup\t10\t10\tRV\tinfant\n",
    )
    corpus, _ = load_curation(root)

    assert validate_mirror_references(root, corpus) == []


def test_reports_a_mirror_row_referencing_an_unknown_dataset(tmp_path: Path) -> None:
    from chd_atlas.corpus import load_curation

    root = _mirror_repo(
        tmp_path,
        "PXD999999\ttof_vs_control\tHGNC:11604\t1.2\t0.001\t0.01\tup\t10\t10\tRV\tinfant\n",
    )
    corpus, _ = load_curation(root)

    issues = validate_mirror_references(root, corpus)

    assert [i.code for i in issues] == ["REF005"]


def test_reports_a_mirror_row_referencing_an_unknown_contrast(tmp_path: Path) -> None:
    from chd_atlas.corpus import load_curation

    root = _mirror_repo(
        tmp_path,
        "PXD012345\tmystery_contrast\tHGNC:11604\t1.2\t0.001\t0.01\tup\t10\t10\tRV\tinfant\n",
    )
    corpus, _ = load_curation(root)

    issues = validate_mirror_references(root, corpus)

    assert [i.code for i in issues] == ["REF006"]
    assert "mystery_contrast" in issues[0].message
