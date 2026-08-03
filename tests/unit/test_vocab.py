# tests/unit/test_vocab.py
import pytest

from chd_atlas.vocab import (
    CLASSIFICATION_RANK,
    MODEL_ORGANISMS,
    Archive,
    Classification,
    EvidenceClass,
    EvidenceStrength,
    FeaturedTopic,
    Inheritance,
    LesionGroup,
    Mechanism,
    Perturbation,
    PhenocopyAssessment,
    Redistribution,
    RescueOutcome,
    SourceTier,
    StudyType,
    SyndromicStatus,
    Technology,
    Zygosity,
    has_conflicting_evidence,
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


def test_strongest_hides_a_refutation_behind_a_stronger_claim() -> None:
    """Deliberate, and the reason has_conflicting_evidence exists.

    Pinned so that changing CLASSIFICATION_RANK alters this on purpose rather
    than by accident.
    """
    assert (
        strongest([Classification.DEFINITIVE, Classification.REFUTED]) is Classification.DEFINITIVE
    )


def test_conflicting_evidence_detects_a_contested_gene() -> None:
    assert has_conflicting_evidence([Classification.DEFINITIVE, Classification.REFUTED])
    assert has_conflicting_evidence([Classification.MODERATE, Classification.DISPUTED])


def test_uncontested_classifications_are_not_conflicting() -> None:
    assert not has_conflicting_evidence([Classification.DEFINITIVE, Classification.MODERATE])
    assert not has_conflicting_evidence([Classification.REFUTED])
    assert not has_conflicting_evidence([])


def test_enum_values_are_unique_within_each_enum() -> None:
    """StrEnum silently aliases a duplicate value instead of raising.

    Iterating an Enum skips aliases, so __members__ is used: it includes them.
    """
    for enum_cls in (
        Archive,
        Classification,
        EvidenceClass,
        EvidenceStrength,
        FeaturedTopic,
        Inheritance,
        LesionGroup,
        Mechanism,
        PhenocopyAssessment,
        Perturbation,
        Redistribution,
        RescueOutcome,
        SourceTier,
        StudyType,
        SyndromicStatus,
        Technology,
        Zygosity,
    ):
        members = enum_cls.__members__
        assert len({m.value for m in members.values()}) == len(members), enum_cls.__name__
