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
"""

from __future__ import annotations

from pathlib import Path

import pytest
from pydantic import ValidationError
from ruamel.yaml import YAML

from chd_atlas.models.phases import CardiacPhase, CardiacPhaseFile, EndBasis


def _phase(id_: str, start_wpc: float, end_wpc: float, **overrides: object) -> CardiacPhase:
    """A validly-constructed `CardiacPhase` for interval-logic tests.

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


def test_end_basis_distinguishes_a_stated_end_from_a_curator_capped_one() -> None:
    """`end_basis` is not decorative: a stated and a capped phase must both
    construct, and must round-trip as genuinely different values -- never
    collapsed to one, and never silently defaulted for the other.
    """
    stated = _phase("septation", 3.0, 5.0, end_basis=EndBasis.STATED)
    capped = _phase("looping", 3.0, 9.0, end_basis=EndBasis.CURATOR_CAPPED)
    assert stated.end_basis is EndBasis.STATED
    assert capped.end_basis is EndBasis.CURATOR_CAPPED
    assert stated.end_basis is not capped.end_basis


def test_end_basis_rejects_a_value_outside_the_two_member_vocabulary() -> None:
    """A closed vocabulary: 'stated' and 'curator_capped' are the only two
    facts a boundary can be, and a third spelling (a typo, or a future
    'estimated' nobody has designed yet) must not silently pass through as a
    string.
    """
    with pytest.raises(ValidationError):
        _phase("a", 1.0, 2.0, end_basis="estimated")


def test_the_committed_phase_vocabulary_loads_and_matches_the_verified_table() -> None:
    """A smoke test over the real file, which no unit fixture can stand in for.

    Pins the actual embryology this atlas ships, verified against the
    boundary table in the task that populated this file (HsapDv day-post-
    fertilization ranges / 7.0; see `curation/cardiac_phases.yaml`'s own
    header for the ordinal-vs-elapsed-week trap this guards against).

    **4 wpc and 7 wpc are both triple matches, but not the same triple.**
    4 wpc: heart looping + atrial septation + ventricular septation (matches
    the concurrency table in `models/phases.py`'s module docstring exactly).
    7 wpc: heart looping + ventricular septation + heart valve morphogenesis
    -- outflow tract septation has already ended by 6.86 wpc, so it does NOT
    appear here even though a same-day illustrative reading might expect it.
    Both are pinned explicitly so a mutant collapsing either count to 1, or
    silently substituting one triple's members for the other's, is caught.

    **Capped vs. stated is pinned by name, not just by count.** A mutant that
    flips which three phases are capped (rather than merely dropping the
    field) would survive a "some are stated, some are capped" assertion; it
    cannot survive naming each side.
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
    assert basis["heart_looping"] is EndBasis.CURATOR_CAPPED
    assert basis["ventricular_septum_morphogenesis"] is EndBasis.CURATOR_CAPPED
    assert basis["heart_valve_morphogenesis"] is EndBasis.CURATOR_CAPPED

    assert [phase.id for phase in parsed.phases_for(4.0)] == [
        "heart_looping",
        "atrial_septum_morphogenesis",
        "ventricular_septum_morphogenesis",
    ]
    assert [phase.id for phase in parsed.phases_for(7.0)] == [
        "heart_looping",
        "ventricular_septum_morphogenesis",
        "heart_valve_morphogenesis",
    ]
