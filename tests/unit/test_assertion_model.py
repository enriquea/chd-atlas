# tests/unit/test_assertion_model.py
import pytest
from pydantic import ValidationError

from chd_atlas.models.assertion import (
    AssertionFile,
    Evidence,
    GeneDiseaseAssertion,
    InTextLocator,
    SupplementaryLocator,
)


def _evidence(**overrides: object) -> dict[str, object]:
    base: dict[str, object] = {
        "publication": "PMID:8988165",
        "evidence_class": "genetic_case",
        "locator": {"kind": "supplementary", "file": "TableS3.xlsx", "sheet": "S3", "row": 42},
        "summary": "De novo nonsense variant in a proband with ASD.",
        "strength": "strong",
    }
    base.update(overrides)
    return base


def _assertion(**overrides: object) -> dict[str, object]:
    base: dict[str, object] = {
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
        "curator": "0000-0002-1825-0097",
        "curated_on": "2026-07-01",
        "last_reviewed": "2026-07-15",
    }
    base.update(overrides)
    return base


def test_parses_a_well_formed_assertion() -> None:
    assertion = GeneDiseaseAssertion.model_validate(_assertion())
    assert assertion.gene == "HGNC:11604"
    assert isinstance(assertion.evidence[0].locator, SupplementaryLocator)
    assert assertion.evidence[0].locator.row == 42


def test_locator_discriminates_on_kind() -> None:
    evidence = Evidence.model_validate(
        _evidence(locator={"kind": "in_text", "element": "Figure 2B"})
    )
    assert isinstance(evidence.locator, InTextLocator)
    assert evidence.locator.element == "Figure 2B"


def test_functional_evidence_class_requires_a_functional_reference() -> None:
    with pytest.raises(ValidationError, match="requires 'functional_evidence'"):
        Evidence.model_validate(_evidence(evidence_class="functional_model"))


def test_functional_reference_rejected_for_non_functional_evidence() -> None:
    with pytest.raises(ValidationError, match="only valid for evidence_class"):
        Evidence.model_validate(_evidence(functional_evidence="CHDA:FUN:0000001"))


def test_functional_evidence_class_accepts_a_functional_reference() -> None:
    evidence = Evidence.model_validate(
        _evidence(evidence_class="functional_model", functional_evidence="CHDA:FUN:0000001")
    )
    assert evidence.functional_evidence == "CHDA:FUN:0000001"


def test_extracardiac_features_rejected_for_isolated_assertions() -> None:
    with pytest.raises(ValidationError, match="cannot be set when syndromic is 'isolated'"):
        GeneDiseaseAssertion.model_validate(
            _assertion(syndromic="isolated", extracardiac_features=["HP:0009777"])
        )


def test_last_reviewed_cannot_precede_curated_on() -> None:
    with pytest.raises(ValidationError, match="cannot precede"):
        GeneDiseaseAssertion.model_validate(
            _assertion(curated_on="2026-07-15", last_reviewed="2026-07-01")
        )


def test_assertion_requires_at_least_one_evidence_item() -> None:
    with pytest.raises(ValidationError):
        GeneDiseaseAssertion.model_validate(_assertion(evidence=[]))


def test_unknown_fields_are_rejected() -> None:
    with pytest.raises(ValidationError):
        GeneDiseaseAssertion.model_validate(_assertion(confidence="high"))


def test_assertion_file_wraps_a_non_empty_list() -> None:
    parsed = AssertionFile.model_validate({"assertions": [_assertion()]})
    assert len(parsed.assertions) == 1
