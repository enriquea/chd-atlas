# tests/unit/test_vocab.py
import pytest

from chd_atlas.vocab import (
    CLASSIFICATION_RANK,
    MODEL_ORGANISMS,
    Classification,
    EvidenceClass,
    LesionGroup,
    SourceTier,
    strongest,
)


def test_classification_values_match_clingen_terms() -> None:
    assert Classification.DEFINITIVE == "definitive"
    assert Classification.NO_KNOWN_ASSOCIATION == "no_known_association"


def test_every_classification_has_a_rank() -> None:
    assert set(CLASSIFICATION_RANK) == set(Classification)


def test_strongest_picks_the_highest_ranked_classification() -> None:
    assert (
        strongest([Classification.LIMITED, Classification.DEFINITIVE]) is Classification.DEFINITIVE
    )
    assert strongest([Classification.DISPUTED, Classification.MODERATE]) is Classification.MODERATE


def test_strongest_rejects_an_empty_sequence() -> None:
    with pytest.raises(ValueError, match="at least one classification"):
        strongest([])


def test_lesion_groups_cover_the_curated_chd_vocabulary() -> None:
    assert {g.value for g in LesionGroup} == {
        "conotruncal",
        "septal",
        "lvoto",
        "rvoto",
        "avsd",
        "heterotaxy",
        "anomalous_venous_return",
        "complex",
        "other",
    }


def test_regulatory_is_an_evidence_class() -> None:
    assert EvidenceClass.REGULATORY == "regulatory"


def test_source_tiers_distinguish_own_curation() -> None:
    assert SourceTier.OWN_CURATION == "own_curation"


def test_model_organisms_map_taxon_ids_to_labels() -> None:
    assert MODEL_ORGANISMS["NCBITaxon:10090"] == "Mus musculus"
    assert "NCBITaxon:7955" in MODEL_ORGANISMS
