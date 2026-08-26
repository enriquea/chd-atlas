# tests/unit/test_models_phases.py
"""The phase vocabulary is transcribed, not authored, so it is shaped like one.

Most fixture values here exercise interval logic only and are not a claim
about embryology -- `_phase()` below fills placeholder Carnegie-stage/HsapDv
fields for exactly that reason. The one exception is
`test_the_committed_phase_vocabulary_loads_and_matches_the_verified_table`,
which reads the real `curation/cardiac_phases.yaml` and pins its actual
embryology, the same way `test_models_scope.py`'s own
`test_the_committed_scope_list_loads_and_names_no_individual` pins the real
`chd_scope.yaml` rather than trusting a unit fixture to stand in for it.

**Phases may overlap, and this file used to forbid that.** Human cardiac
morphogenetic processes run concurrently -- see `models/phases.py`'s module
docstring for the measured concurrency table this vocabulary is built to
represent -- so `phases_for` returns a tuple, possibly with more than one
member, and `CardiacPhaseFile` no longer rejects an overlapping pair at load
time. Only duplicate ids are still rejected.

**A phase may also have no stated end at all, and this file used to invent
one.** An earlier revision capped every open-ended phase at a curator-chosen
boundary (`end_basis: curator_capped`); that concept is retired, because two
of the three phases originally capped turned out to have a real, stated end
the first extraction missed, and the third (`heart_looping`) must never be
asserted as running at a stage nothing gives a basis for.

`EndBasis` has **three** members since 2026-08-24 -- `STATED`, `DERIVED` and
`NOT_STATED`. `DERIVED` was added for `heart_looping`, whose end follows from
GO:0001947's own definition ("ends ... preceding septation") read against a
septation boundary this file already carries as `STATED`. It is deliberately
*not* `curator_capped` returning: that took a sentence about a different
process, this is an entailment from the phase's own ontology term. See
`EndBasis` in `models/phases.py`, where the distinction is argued in full.

**No curated phase is `NOT_STATED` any more**, and the member is kept and
exercised here on fixtures precisely for that reason -- `_open_phase()` below
builds that shape, so the rule survives the corpus no longer containing an
instance of it.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from pydantic import ValidationError
from ruamel.yaml import YAML

from chd_atlas.models.phases import CardiacPhase, CardiacPhaseFile, EndBasis


def _phase(id_: str, start_wpc: float, end_wpc: float, **overrides: object) -> CardiacPhase:
    """A validly-constructed, fully-`STATED` `CardiacPhase` for interval-logic tests.

    The Carnegie-stage/HsapDv/`go_id` values are placeholders, not a claim
    about embryology -- only `id`, `start_wpc` and `end_wpc` vary across the
    tests that use this factory. `label` defaults to `id` verbatim so a
    caller need not repeat it for every fixture phase.
    """
    base: dict[str, object] = {
        "id": id_,
        "go_id": "GO:0000001",
        "label": id_,
        "start_wpc": start_wpc,
        "end_wpc": end_wpc,
        "start_carnegie_stage": "CS1",
        "end_carnegie_stage": "CS2",
        "start_hsapdv_id": "HsapDv:0000001",
        "end_hsapdv_id": "HsapDv:0000002",
        "end_basis": EndBasis.STATED,
    }
    base.update(overrides)
    return CardiacPhase.model_validate(base)


def _open_phase(id_: str, start_wpc: float, **overrides: object) -> CardiacPhase:
    """A validly-constructed `CardiacPhase` with a real start and NO stated end.

    The `heart_looping` shape: every end field is `None` and `end_basis` is
    `NOT_STATED`, matching what `end_fields_match_end_basis` requires.
    """
    base: dict[str, object] = {
        "id": id_,
        "go_id": "GO:0000001",
        "label": id_,
        "start_wpc": start_wpc,
        "end_wpc": None,
        "start_carnegie_stage": "CS1",
        "end_carnegie_stage": None,
        "start_hsapdv_id": "HsapDv:0000001",
        "end_hsapdv_id": None,
        "end_basis": EndBasis.NOT_STATED,
    }
    base.update(overrides)
    return CardiacPhase.model_validate(base)


def _file(**overrides: object) -> dict[str, object]:
    base: dict[str, object] = {
        "attributed_to": "O'Rahilly & Muller 1987",
        "citation": "ISBN:0872796248",
        "phases": [],
    }
    base.update(overrides)
    return base


def test_phases_are_half_open() -> None:
    """Half-open [start, end) so a stage on a boundary falls in exactly one phase.

    Under closed intervals an 8 wpc stage would fall in both a 5-8 and an 8-10
    phase, and a reader could not tell adjacency from an overlap.
    """
    file = CardiacPhaseFile(
        attributed_to="O'Rahilly & Muller 1987",
        citation="ISBN:0872796248",
        phases=[
            _phase("looping", 3.0, 5.0),
            _phase("septation", 5.0, 8.0),
        ],
    )
    matched = file.phases_for(5.0)
    assert [phase.id for phase in matched] == ["septation"]  # start is inclusive
    assert file.phases_for(8.0) == ()  # end is exclusive
    assert file.phases_for(2.0) == ()  # before every phase


def test_overlapping_phases_load_and_both_match_a_shared_wpc() -> None:
    """The disjointness check is gone: two overlapping phases must both load,
    and a wpc inside the overlap must report BOTH, not one arbitrarily chosen.

    This is the unique killer of a reinstated disjointness check: the old
    validator raised `ValidationError` on exactly this fixture (matching
    the old `test_phases_are_half_open_and_may_not_overlap`), so if anyone
    ever restores that check, `CardiacPhaseFile(...)` below raises and this
    test fails at construction, before the `phases_for` assertions even run.
    """
    file = CardiacPhaseFile(
        attributed_to="x",
        citation="PMID:1",
        phases=[
            _phase("a", 3.0, 6.0),
            _phase("b", 5.0, 8.0),
        ],
    )
    assert {phase.id for phase in file.phases_for(5.5)} == {"a", "b"}


def test_a_wpc_covered_by_one_phase_is_distinguished_from_one_covered_by_three() -> None:
    """The fixture-diversity guard this project's own history calls for
    (CLAUDE.md 4.14/15b/30/36): a fixture where every match count is the same
    number cannot prove `phases_for` returns ALL matches rather than
    hard-coding a length. Three phases overlap at 4.5; only one of them
    (`a`) also reaches 1.5, and none reaches 9.0.
    """
    file = CardiacPhaseFile(
        attributed_to="x",
        citation="PMID:1",
        phases=[
            _phase("a", 1.0, 6.0),
            _phase("b", 4.0, 7.0),
            _phase("c", 4.0, 5.0),
        ],
    )
    assert [phase.id for phase in file.phases_for(1.5)] == ["a"]
    assert [phase.id for phase in file.phases_for(4.5)] == ["a", "b", "c"]
    assert file.phases_for(9.0) == ()


def test_phases_for_is_sorted_by_start_then_id_regardless_of_declaration_order() -> None:
    """Declared out of both start-order and id-order, to make a dropped sort
    compare (or return) the wrong sequence rather than merely look
    coincidentally right -- the same reasoning
    `test_prf006_reports_one_issue_per_interior_gap_in_wpc_order` applies on
    the validator side.
    """
    file = CardiacPhaseFile(
        attributed_to="x",
        citation="PMID:1",
        phases=[
            _phase("z_late", 5.0, 9.0),
            _phase("b_tied", 1.0, 9.0),
            _phase("a_tied", 1.0, 9.0),
        ],
    )
    assert [phase.id for phase in file.phases_for(6.0)] == ["a_tied", "b_tied", "z_late"]


def test_a_null_wpc_is_outside_every_phase() -> None:
    """A post-natal stage carries no wpc; that is not evidence against a phase, or for one."""
    file = CardiacPhaseFile(**_file(phases=[_phase("looping", 3.0, 5.0)]))
    assert file.phases_for(None) == ()


def test_a_phase_must_be_attributed() -> None:
    """The atlas transcribes embryology; it does not author it."""
    with pytest.raises(ValidationError):
        CardiacPhaseFile(attributed_to="", citation="PMID:1", phases=[])


def test_duplicate_phase_ids_are_rejected() -> None:
    """Two entries claiming one id is not two disjoint intervals -- something upstream is wrong.

    The two intervals here (1-2, 3-4) do not overlap, so this is the unique
    killer of the duplicate-id guard: an overlap guard (even if one still
    existed) could not also catch it.
    """
    with pytest.raises(ValidationError, match="duplicate"):
        CardiacPhaseFile(
            **_file(
                phases=[
                    _phase("a", 1.0, 2.0),
                    _phase("a", 3.0, 4.0, label="A again"),
                ]
            )
        )


@pytest.mark.parametrize(
    "start_wpc, end_wpc",
    [
        pytest.param(5.0, 5.0, id="equal"),
        pytest.param(5.0, 3.0, id="reversed"),
    ],
)
def test_end_wpc_must_exceed_start_wpc(start_wpc: float, end_wpc: float) -> None:
    """An empty or a reversed interval names no window at all."""
    with pytest.raises(ValidationError, match="end_wpc must exceed start_wpc"):
        _phase("a", start_wpc, end_wpc)


@pytest.mark.parametrize("field", ["attributed_to", "citation"])
def test_whitespace_only_attribution_is_rejected(field: str) -> None:
    """A value of all whitespace passes a bare `min_length` and renders as nothing.

    That is worse than an absent field because it looks answered. Same defect,
    same reasoning, as `Dataset.licence_is_not_blank`: a blank here records no
    provenance for a developmental-biology claim at all.
    """
    with pytest.raises(ValidationError, match="must not be blank"):
        CardiacPhaseFile(**_file(**{field: "   "}))


def test_end_basis_distinguishes_a_stated_end_from_a_not_stated_one() -> None:
    """`end_basis` is not decorative: a stated and an open-ended phase must
    both construct, and must round-trip as genuinely different values --
    never collapsed to one, and never silently defaulted for the other.
    """
    stated = _phase("septation", 3.0, 5.0, end_basis=EndBasis.STATED)
    open_ended = _open_phase("looping", 3.0)
    assert stated.end_basis is EndBasis.STATED
    assert stated.end_wpc == 5.0
    assert open_ended.end_basis is EndBasis.NOT_STATED
    assert open_ended.end_wpc is None
    assert stated.end_basis is not open_ended.end_basis


def test_end_basis_rejects_a_value_outside_the_vocabulary() -> None:
    """A closed vocabulary: 'stated', 'derived' and 'not_stated' are the only
    three facts a boundary can be, and a fourth spelling (a typo, the retired
    'curator_capped', or a future value nobody has designed yet) must not
    silently pass through as a string.

    `curator_capped` is named explicitly rather than left to the generic case.
    It is the one rejected value that was once *valid* here, so it is the one a
    stale YAML file or an old branch could actually carry -- and the one whose
    meaning `DERIVED` is closest to and must not be confused with.
    """
    with pytest.raises(ValidationError):
        _phase("a", 1.0, 2.0, end_basis="estimated")
    with pytest.raises(ValidationError):
        _phase("a", 1.0, 2.0, end_basis="curator_capped")


@pytest.mark.parametrize("basis", [EndBasis.STATED, EndBasis.DERIVED])
@pytest.mark.parametrize("missing_field", ["end_wpc", "end_carnegie_stage", "end_hsapdv_id"])
def test_a_basis_that_claims_an_end_requires_every_end_field(
    missing_field: str, basis: EndBasis
) -> None:
    """An `end_basis` that claims a boundary, with an end field left `None`, is
    a claimed boundary with no number behind it -- caught for each end field
    independently, so a guard checking only `end_wpc` cannot silently let the
    other two drift.

    **Parametrised over the basis, not just `STATED`, and that is the point.**
    `end_fields_match_end_basis` is written as `is not NOT_STATED` rather than
    as a list of the members that do require fields, precisely so a member
    added later inherits the requirement. Measured 2026-08-24: reverting it to
    `is EndBasis.STATED` -- which is what it said before `DERIVED` existed --
    passed the whole suite. The defensive phrasing was itself unguarded
    (CLAUDE.md section 4.40), and only a case that exercises the non-`STATED`
    branch can tell the two spellings apart.

    Built from an explicit dict rather than through `_phase()` -- `_phase()`
    binds `end_wpc` positionally, which collides with also overriding it by
    keyword when `missing_field == "end_wpc"`.
    """
    fields: dict[str, object] = {
        "id": "a",
        "go_id": "GO:0000001",
        "label": "a",
        "start_wpc": 1.0,
        "end_wpc": 2.0,
        "start_carnegie_stage": "CS1",
        "end_carnegie_stage": "CS2",
        "start_hsapdv_id": "HsapDv:0000001",
        "end_hsapdv_id": "HsapDv:0000002",
        "end_basis": basis,
    }
    fields[missing_field] = None
    with pytest.raises(ValidationError, match="end field is missing"):
        CardiacPhase.model_validate(fields)


@pytest.mark.parametrize(
    "extra_field, value",
    [("end_wpc", 9.0), ("end_carnegie_stage", "CS9"), ("end_hsapdv_id", "HsapDv:9")],
)
def test_end_basis_not_stated_forbids_every_end_field(extra_field: str, value: object) -> None:
    """The reverse direction: `end_basis: not_stated` with a real value in any
    end field would be silently treated as sourced by `phases_for` and
    `_coverage_spans` (both key on `end_wpc is None`), reintroducing the
    invented-boundary defect `NOT_STATED` exists to rule out. Checked per
    field for the same reason as the STATED direction above.
    """
    with pytest.raises(ValidationError, match="end field is set"):
        _open_phase("a", 1.0, **{extra_field: value})


def test_a_phase_with_no_stated_end_is_never_returned_by_phases_for() -> None:
    """Unconditional, at every wpc -- not merely 'too far past its start'.

    There is no source-given point at which this atlas could say an
    open-ended phase has stopped, so picking any cutoff (its own start, a
    round number, the vocabulary's last known end) would itself be an
    invented boundary. Checked at the phase's own start, just after it, and
    far past every other phase in the file, so a mutant that excludes it only
    near one of those three regions is still caught.
    """
    file = CardiacPhaseFile(
        attributed_to="x",
        citation="PMID:1",
        phases=[_open_phase("looping", 3.0), _phase("septation", 5.0, 8.0)],
    )
    assert file.phases_for(3.0) == ()
    assert file.phases_for(4.0) == ()
    assert [phase.id for phase in file.phases_for(6.0)] == ["septation"]
    assert file.phases_for(100.0) == ()


def test_the_committed_phase_vocabulary_loads_and_matches_the_verified_table() -> None:
    """A smoke test over the real file, which no unit fixture can stand in for.

    Pins the actual embryology this atlas ships, verified against the
    boundary table supplied for this file (HsapDv day-post-fertilization
    ranges / 7.0; see `curation/cardiac_phases.yaml`'s own header for the
    ordinal-vs-elapsed-week trap this guards against, and for the correction
    history: two of the three phases originally marked `curator_capped` had
    a real stated end the first extraction simply missed; only `heart_
    looping` genuinely has none).

    **`heart_looping` never appears in any `phases_for` result here**, even
    though it is curated (with a real, stated start) and even though it is
    biologically underway throughout 3.14-8.57 wpc -- it has no stated end,
    so it is structurally excluded (`end_wpc is None`) rather than asserted
    past a point the source never gave. A mutant that lets it leak back into
    a match is exactly what this test's absence of `"heart_looping"` from
    every list below is watching for.

    **Every wpc pinned gives a DIFFERENT match set**, not just a different
    count, per the fixture-diversity rule: 1 phase, 2 phases, 3 phases, 1
    phase again (a different one), then 0 -- a mutant that returns the wrong
    *members* while keeping the right *count* cannot survive all five.
    """
    path = Path(__file__).parent.parent.parent / "curation" / "cardiac_phases.yaml"
    raw = path.read_text(encoding="utf-8")
    parsed = CardiacPhaseFile.model_validate(YAML(typ="safe").load(raw))

    assert parsed.attributed_to.strip() != ""
    assert "32048790" in parsed.citation

    ids = {phase.id for phase in parsed.phases}
    assert ids == {
        "embryonic_heart_tube_morphogenesis",
        "heart_looping",
        "atrial_septum_morphogenesis",
        "ventricular_septum_morphogenesis",
        "outflow_tract_septum_morphogenesis",
        "heart_valve_morphogenesis",
    }

    go_ids = {phase.id: phase.go_id for phase in parsed.phases}
    assert go_ids == {
        "embryonic_heart_tube_morphogenesis": "GO:0003143",
        "heart_looping": "GO:0001947",
        "atrial_septum_morphogenesis": "GO:0060413",
        "ventricular_septum_morphogenesis": "GO:0060412",
        "outflow_tract_septum_morphogenesis": "GO:0003148",
        "heart_valve_morphogenesis": "GO:0003179",
    }

    basis = {phase.id: phase.end_basis for phase in parsed.phases}
    assert basis["embryonic_heart_tube_morphogenesis"] is EndBasis.STATED
    assert basis["atrial_septum_morphogenesis"] is EndBasis.STATED
    assert basis["outflow_tract_septum_morphogenesis"] is EndBasis.STATED
    assert basis["ventricular_septum_morphogenesis"] is EndBasis.STATED
    assert basis["heart_valve_morphogenesis"] is EndBasis.STATED
    assert basis["heart_looping"] is EndBasis.DERIVED

    by_id = {phase.id: phase for phase in parsed.phases}
    assert by_id["heart_looping"].start_wpc == pytest.approx(3.14)
    # The derived end, pinned against literals rather than read back from the
    # septation phase it is anchored to. Anchoring the assertion the same way
    # the value is anchored would pass for any pair that happened to agree,
    # including both being wrong together; these three numbers are what
    # `curation/cardiac_phases.yaml` must actually publish.
    assert by_id["heart_looping"].end_wpc == pytest.approx(3.71)
    assert by_id["heart_looping"].end_carnegie_stage == "CS12"
    assert by_id["heart_looping"].end_hsapdv_id == "HsapDv:0000019"
    # ...and they are the same boundary septation starts at, which is the
    # entailment `end_basis: derived` records. Asserted as an equality so the
    # two cannot drift apart silently: if a curator moves septation's start,
    # this fails rather than leaving looping ending where nothing begins.
    assert by_id["heart_looping"].end_wpc == by_id["atrial_septum_morphogenesis"].start_wpc
    assert by_id["ventricular_septum_morphogenesis"].end_wpc == pytest.approx(7.29)
    assert by_id["heart_valve_morphogenesis"].end_wpc == pytest.approx(8.57)

    # heart_looping now *does* match, and this assertion is the one that
    # changed when its end became derived rather than absent. Through
    # 2026-08-24 it read `== ["embryonic_heart_tube_morphogenesis"]`, because a
    # phase with a null `end_wpc` is excluded at every wpc including its own
    # start. With an end it participates, and 3.14 is the wpc where that is
    # most visible: both phases genuinely contain it, since
    # embryonic_heart_tube_morphogenesis runs [2.71, 3.29) and heart_looping
    # runs [3.14, 3.71). The overlap is real biology, not a transcription
    # error -- see the header of `curation/cardiac_phases.yaml`.
    assert [phase.id for phase in parsed.phases_for(3.14)] == [
        "embryonic_heart_tube_morphogenesis",
        "heart_looping",
    ]

    # The two ends of heart_looping's own window, because a half-open interval
    # is exactly where an off-by-one hides and neither bound is checked above.
    # Just inside its end, it is the only match -- the heart tube closed at
    # 3.29. Exactly at its end it is gone and septation has begun, which is
    # the entailment `derived` encodes: the two phases meet at 3.71 with no
    # gap and no overlap.
    assert [phase.id for phase in parsed.phases_for(3.70)] == ["heart_looping"]
    assert [phase.id for phase in parsed.phases_for(3.71)] == [
        "atrial_septum_morphogenesis",
        "ventricular_septum_morphogenesis",
    ]
    assert [phase.id for phase in parsed.phases_for(4.0)] == [
        "atrial_septum_morphogenesis",
        "ventricular_septum_morphogenesis",
    ]
    assert [phase.id for phase in parsed.phases_for(6.0)] == [
        "atrial_septum_morphogenesis",
        "ventricular_septum_morphogenesis",
        "outflow_tract_septum_morphogenesis",
    ]
    assert [phase.id for phase in parsed.phases_for(7.0)] == [
        "ventricular_septum_morphogenesis",
        "heart_valve_morphogenesis",
    ]
    assert [phase.id for phase in parsed.phases_for(8.0)] == ["heart_valve_morphogenesis"]
    assert parsed.phases_for(9.0) == ()
    assert parsed.phases_for(100.0) == ()
