import pytest
from pydantic import ValidationError

from chd_atlas.models.functional import FunctionalEvidence, FunctionalFile


def _functional(**overrides: object) -> dict[str, object]:
    base: dict[str, object] = {
        "id": "CHDA:FUN:0000001",
        "gene": "HGNC:11604",
        "organism": "NCBITaxon:10090",
        "perturbation": "knockout",
        "zygosity": "heterozygous",
        "cardiac_phenotype": ["MP:0010402"],
        "phenocopies_human": "partial",
        "rescue_performed": False,
        "publication": "PMID:11572777",
    }
    base.update(overrides)
    return base


def test_parses_a_well_formed_functional_record() -> None:
    record = FunctionalEvidence.model_validate(_functional())
    assert record.organism == "NCBITaxon:10090"
    assert record.rescue_performed is False


def test_rejects_an_organism_outside_the_allowlist() -> None:
    with pytest.raises(ValidationError, match="not an allowed model organism"):
        FunctionalEvidence.model_validate(_functional(organism="NCBITaxon:6239"))


def test_requires_at_least_one_cardiac_phenotype() -> None:
    with pytest.raises(ValidationError):
        FunctionalEvidence.model_validate(_functional(cardiac_phenotype=[]))


def test_rejects_a_malformed_model_phenotype_term() -> None:
    with pytest.raises(ValidationError):
        FunctionalEvidence.model_validate(_functional(cardiac_phenotype=["HP:0001631"]))


def test_functional_file_wraps_a_non_empty_list() -> None:
    parsed = FunctionalFile.model_validate({"functional_evidence": [_functional()]})
    assert len(parsed.functional_evidence) == 1
