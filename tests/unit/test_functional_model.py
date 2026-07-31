import pytest
from pydantic import ValidationError

from chd_atlas.models.functional import FunctionalEvidence, FunctionalFile
from chd_atlas.vocab import MODEL_ORGANISMS


def _functional(**overrides: object) -> dict[str, object]:
    base: dict[str, object] = {
        "id": "CHDA:FUN:0000001",
        "gene": "HGNC:11604",
        "organism": "NCBITaxon:10090",
        "perturbation": "knockout",
        "zygosity": "heterozygous",
        "cardiac_phenotype": ["MP:0010402"],
        "phenocopies_human": "partial",
        "rescue_outcome": "not_attempted",
        "publication": "PMID:11572777",
    }
    base.update(overrides)
    return base


def test_parses_a_well_formed_functional_record() -> None:
    record = FunctionalEvidence.model_validate(_functional())
    assert record.organism == "NCBITaxon:10090"
    assert record.rescue_outcome == "not_attempted"


def test_rejects_an_organism_outside_the_allowlist() -> None:
    with pytest.raises(ValidationError, match="not an allowed model organism"):
        FunctionalEvidence.model_validate(_functional(organism="NCBITaxon:6239"))


def test_requires_at_least_one_cardiac_phenotype() -> None:
    with pytest.raises(ValidationError):
        FunctionalEvidence.model_validate(_functional(cardiac_phenotype=[]))


def test_rejects_a_malformed_model_phenotype_term() -> None:
    with pytest.raises(ValidationError):
        FunctionalEvidence.model_validate(_functional(cardiac_phenotype=["MONDO:0007259"]))


def test_functional_file_wraps_a_non_empty_list() -> None:
    parsed = FunctionalFile.model_validate({"functional_evidence": [_functional()]})
    assert len(parsed.functional_evidence) == 1


def test_morpholino_rejected_for_a_mammalian_organism() -> None:
    with pytest.raises(ValidationError, match="morpholino"):
        FunctionalEvidence.model_validate(
            _functional(perturbation="morpholino", organism="NCBITaxon:10090", zygosity="n/a")
        )


def test_morpholino_accepted_for_zebrafish() -> None:
    record = FunctionalEvidence.model_validate(
        _functional(perturbation="morpholino", organism="NCBITaxon:7955", zygosity="n/a")
    )
    assert record.perturbation == "morpholino"


def test_germline_perturbation_must_state_zygosity() -> None:
    with pytest.raises(ValidationError, match="zygosity"):
        FunctionalEvidence.model_validate(_functional(perturbation="knockout", zygosity="n/a"))


def test_transient_perturbation_may_omit_zygosity() -> None:
    record = FunctionalEvidence.model_validate(
        _functional(perturbation="overexpression", zygosity="n/a")
    )
    assert record.zygosity == "n/a"


def test_human_ipsc_model_can_record_an_hpo_phenotype() -> None:
    """MODEL_ORGANISMS allows human iPSC-CM models, which have no MP or ZP terms."""
    record = FunctionalEvidence.model_validate(
        _functional(organism="NCBITaxon:9606", cardiac_phenotype=["HP:0001631"])
    )
    assert record.cardiac_phenotype == ["HP:0001631"]


@pytest.mark.parametrize("organism", sorted(MODEL_ORGANISMS))
def test_every_allowlisted_organism_is_accepted(organism: str) -> None:
    record = FunctionalEvidence.model_validate(_functional(organism=organism))
    assert record.organism == organism


def test_functional_file_does_not_reject_duplicate_ids() -> None:
    """Intentional: repo-wide atlas ID uniqueness is enforced by validate_ids()."""
    parsed = FunctionalFile.model_validate({"functional_evidence": [_functional(), _functional()]})
    assert len(parsed.functional_evidence) == 2
