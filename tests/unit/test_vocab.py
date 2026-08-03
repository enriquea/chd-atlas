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


def test_every_mirrored_classification_maps_and_an_unknown_one_raises() -> None:
    """The map is total over the allowed sets, and has no fallback.

    A fallback would let a term an authority adds in a later release be silently
    coerced to something adjacent -- publishing a confidence nobody asserted.
    Refusing is the whole point: an unmapped term must reach a human.

    `Supportive` maps to None deliberately. It is GenCC's bucket for submitters
    that assert an association without grading evidence, so it has no rung on
    the ladder and must not be given one.

    Binds to the schemas' own `classification` column rather than importing
    `tables._CLINGEN_CLASSIFICATIONS` / `_GENCC_CLASSIFICATIONS` directly: those
    names are private to `tables.py`, and `Column.allowed` is already the
    schema's public way of exposing its vocabulary (see
    `test_tables.py::test_gencc_keys_on_the_submitter_because_it_publishes_no_verdict`,
    which reads the same column the same way). Going through that door instead
    of a private one keeps this test from opening its own back channel into
    another module -- and still fails exactly the same way if the two
    definitions drift: adding a term to one allowed set and not the other, or to
    a map and not the schema, still trips the set-equality assertions below.
    """
    from chd_atlas.tables import CLINGEN_VALIDITY, GENCC_SUBMISSIONS
    from chd_atlas.vocab import (
        CLINGEN_CLASSIFICATIONS,
        GENCC_CLASSIFICATIONS,
        Classification,
    )

    clingen_allowed = next(
        c for c in CLINGEN_VALIDITY.columns if c.name == "classification"
    ).allowed
    gencc_allowed = next(c for c in GENCC_SUBMISSIONS.columns if c.name == "classification").allowed
    assert clingen_allowed is not None
    assert gencc_allowed is not None
    assert set(CLINGEN_CLASSIFICATIONS) == clingen_allowed
    assert set(GENCC_CLASSIFICATIONS) == gencc_allowed
    assert CLINGEN_CLASSIFICATIONS["No Known Disease Relationship"] is (
        Classification.NO_KNOWN_ASSOCIATION
    )
    assert GENCC_CLASSIFICATIONS["Disputed Evidence"] is Classification.DISPUTED
    assert GENCC_CLASSIFICATIONS["Supportive"] is None
    assert "Supportive" not in CLINGEN_CLASSIFICATIONS


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
