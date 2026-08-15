# tests/unit/test_models_phases.py
"""The phase vocabulary is transcribed, not authored, so it is shaped like one.

Every fixture value here exercises the interval logic; none is a claim about
embryology. The values that ARE a claim about embryology live only in
`curation/cardiac_phases.yaml`, which this task ships empty on purpose -- no
source has been verified yet.
"""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from chd_atlas.models.phases import CardiacPhase, CardiacPhaseFile


def _file(**overrides: object) -> dict[str, object]:
    base: dict[str, object] = {
        "attributed_to": "O'Rahilly & Muller 1987",
        "citation": "ISBN:0872796248",
        "phases": [],
    }
    base.update(overrides)
    return base


def test_phases_are_half_open_and_may_not_overlap() -> None:
    """Half-open [start, end) so a stage on a boundary falls in exactly one phase.

    Under closed intervals an 8 wpc stage would fall in both a 5-8 and an 8-10
    phase, and "overlapping" would then be indistinguishable from "adjacent".
    """
    file = CardiacPhaseFile(
        attributed_to="O'Rahilly & Muller 1987",
        citation="ISBN:0872796248",
        phases=[
            CardiacPhase(id="looping", label="Cardiac looping", start_wpc=3.0, end_wpc=5.0),
            CardiacPhase(id="septation", label="Septation", start_wpc=5.0, end_wpc=8.0),
        ],
    )
    assert file.phase_for(5.0) is not None
    assert file.phase_for(5.0).id == "septation"  # start is inclusive
    assert file.phase_for(8.0) is None  # end is exclusive
    assert file.phase_for(2.0) is None  # before every phase

    with pytest.raises(ValidationError, match="overlap"):
        CardiacPhaseFile(
            attributed_to="x",
            citation="PMID:1",
            phases=[
                CardiacPhase(id="a", label="A", start_wpc=3.0, end_wpc=6.0),
                CardiacPhase(id="b", label="B", start_wpc=5.0, end_wpc=8.0),
            ],
        )


def test_a_phase_must_be_attributed() -> None:
    """The atlas transcribes embryology; it does not author it."""
    with pytest.raises(ValidationError):
        CardiacPhaseFile(attributed_to="", citation="PMID:1", phases=[])


def test_a_null_wpc_is_outside_every_phase() -> None:
    """A post-natal stage carries no wpc; that is not evidence against a phase, or for one."""
    file = CardiacPhaseFile(
        **_file(
            phases=[CardiacPhase(id="looping", label="Cardiac looping", start_wpc=3.0, end_wpc=5.0)]
        )
    )
    assert file.phase_for(None) is None


def test_duplicate_phase_ids_are_rejected() -> None:
    """Two entries claiming one id is not two disjoint intervals -- something upstream is wrong.

    The two intervals here (1-2, 3-4) do not overlap, so this is the unique
    killer of the duplicate-id guard: the overlap guard cannot also catch it.
    """
    with pytest.raises(ValidationError, match="duplicate"):
        CardiacPhaseFile(
            **_file(
                phases=[
                    CardiacPhase(id="a", label="A", start_wpc=1.0, end_wpc=2.0),
                    CardiacPhase(id="a", label="A again", start_wpc=3.0, end_wpc=4.0),
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
        CardiacPhase(id="a", label="A", start_wpc=start_wpc, end_wpc=end_wpc)


@pytest.mark.parametrize("field", ["attributed_to", "citation"])
def test_whitespace_only_attribution_is_rejected(field: str) -> None:
    """A value of all whitespace passes a bare `min_length` and renders as nothing.

    That is worse than an absent field because it looks answered. Same defect,
    same reasoning, as `Dataset.licence_is_not_blank`: a blank here records no
    provenance for a developmental-biology claim at all.
    """
    with pytest.raises(ValidationError, match="must not be blank"):
        CardiacPhaseFile(**_file(**{field: "   "}))
