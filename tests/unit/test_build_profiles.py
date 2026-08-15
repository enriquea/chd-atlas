# tests/unit/test_build_profiles.py
"""Tests for `build/profiles.py`: the percentile band and tau (Tasks 9-10),
phase assignment (Task 11) and quantile shard emission (Task 12).

`percentile_of` and `band` are pure lookups over a hand-built grid, so every
expected value here is read off `GRID` by construction rather than computed by
the code under test. `placement` adds the D41 detection-floor gate and the
published payload shape.

Verified by direct computation before writing these assertions
(`.venv/bin/python`, `bisect.bisect_left`): `GRID[46] == 0.1`,
`GRID[50] == 0.5`, `GRID[100] == 5.5`, `bisect_left(GRID, 5.5) == 100` (an
ordinary match, no clamp needed) and `bisect_left(GRID, 99.0) == 101` (the
`min(..., 100)` clamp is what brings that one back to 100) -- so the plan's
sketch literals for `percentile_of` are correct as written and exercise two
genuinely different code paths, not one value asserted twice.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from chd_atlas.build.emit import Emitter
from chd_atlas.build.omics import build_omics
from chd_atlas.build.profiles import (
    LOOKUP_RULE,
    TAU_METHOD,
    PhaseOutcome,
    assign_phase,
    band,
    build_profile_quantiles,
    gene_expression_profiles,
    percentile_annotations,
    percentile_of,
    placement,
    specificity,
)
from chd_atlas.models.dataset import Dataset, Stage
from chd_atlas.models.phases import CardiacPhase, CardiacPhaseFile

# 101 breakpoints, percentile i at index i. Hand-built so the expected answers
# below are read off the grid rather than computed by the code under test.
GRID = [0.0] * 46 + [round(0.1 * (i - 45), 4) for i in range(46, 101)]
# GRID[46] == 0.1, GRID[50] == 0.5, GRID[100] == 5.5


def test_a_value_reports_the_lowest_percentile_of_a_tie() -> None:
    """The lookup rule is named because a tie is 46 percentile points wide here.

    GRID's first 46 breakpoints are all 0.0 -- roughly half a bulk library
    sits at zero. Reading the *top* of that tie would publish a gene at zero
    as the 45th percentile. The rule is the bottom of the tie, which
    understates rather than overstates, and D41 removes those genes from play
    anyway by refusing a percentile below the detection floor.
    """
    assert percentile_of(0.0, GRID) == 0
    assert percentile_of(0.1, GRID) == 46
    assert percentile_of(0.5, GRID) == 50
    assert percentile_of(5.5, GRID) == 100
    assert percentile_of(99.0, GRID) == 100  # above every breakpoint, clamped


def test_the_band_is_q25_to_q75_and_falls_back_when_they_are_null() -> None:
    """A point estimate at n=2 has an honest range of roughly top 2% to top 15%.

    Publishing the band costs one extra lookup and no sampling model -- the
    numbers are already in the mirror. `None` passes straight through rather
    than falling back to the median's own percentile, which would publish a
    band with no reported width as though it had been measured.
    """
    assert band(median=0.5, q25=0.1, q75=5.5, breakpoints=GRID) == (46, 50, 100)
    assert band(median=0.5, q25=None, q75=None, breakpoints=GRID) == (None, 50, None)


# `unit`/`n_samples`/`n_genes` are required, non-defaulted arguments to
# `placement` (see the module docstring for why) -- the plan's Step 1 sketch
# called `placement` with only five arguments, which no longer matches the
# five-argument-plus-payload-fields signature Step 3 asks for. Threaded
# through every call below rather than defaulted, matching this project's
# stance against silent defaults for a field a page must be able to render.
_UNIT = "rpkm"
_N_SAMPLES = 3
_N_GENES = 19842


@pytest.mark.parametrize(
    ("median", "floor", "expect_placed"),
    [
        (0.0, 1.0, False),  # below the floor, and also zero
        (0.5, 1.0, False),  # below the floor, but NOT zero -- see the docstring
        (1.0, 1.0, True),  # exactly AT the floor: the pinned boundary
        (2.0, 1.0, True),  # clearly above
    ],
)
def test_the_floor_gate_is_strict_less_than(
    median: float, floor: float, expect_placed: bool
) -> None:
    """D41, and the boundary judgement call this module makes: `<`, not `<=`.

    Below the floor, a percentile is an artefact of the tied breakpoints at
    the bottom of the grid, not a measurement -- `None`, not a number a page
    would render as though it meant something.

    A median exactly equal to the declared floor counts as detected, not
    below it: the floor is a source's own declared threshold for "present at
    all", so it is itself a detected value, and `<=` would exclude the one
    point the source chose as its own boundary. Task 10's `specificity` uses
    the same comparison (`peak_raw < floor`) on a gene's peak across organs,
    so one gene cannot read "detected" under one figure and "below the floor"
    under the other from the same value.

    0.5 is a second below-floor case distinct from 0.0: a mutant that checked
    `median == 0` instead of the real comparison would still pass the 0.0 case
    and only fail here.
    """
    result = placement(
        median=median,
        q25=None,
        q75=None,
        breakpoints=GRID,
        floor=floor,
        unit=_UNIT,
        n_samples=_N_SAMPLES,
        n_genes=_N_GENES,
    )
    assert (result is not None) == expect_placed


def test_placement_carries_the_band_and_every_qualifier_needed_to_render_it() -> None:
    """Design of `Placement`: nothing a page needs is left to a default.

    `q25_percentile` (46), `median_percentile` (50) and `q75_percentile` (100)
    are three different numbers on purpose -- a fixture where any two of them
    coincided would let a swap between them survive. `n_samples` and `n_genes`
    are likewise two different numbers so a mutant that dropped or confused
    either is caught: `n_samples` is the *count backing this gene's median*
    (the `count_unit` precedent -- the schema permits `n_samples` = 1, so
    "median" can be a single observation and the page must say so), `n_genes`
    is the *size of the ranked universe* the percentile is a rank within
    (D39(b) -- "top 4%" names no population without its denominator).
    `method` carries `LOOKUP_RULE` (D39(a)) so a consumer re-deriving this
    figure uses the same tie-breaking rule this module used.
    """
    result = placement(
        median=0.5,
        q25=0.1,
        q75=5.5,
        breakpoints=GRID,
        floor=0.05,
        unit=_UNIT,
        n_samples=_N_SAMPLES,
        n_genes=_N_GENES,
    )
    assert result is not None
    assert result["q25_percentile"] == 46
    assert result["median_percentile"] == 50
    assert result["q75_percentile"] == 100
    assert result["median_abundance"] == 0.5
    assert result["unit"] == _UNIT
    assert result["n_samples"] == _N_SAMPLES
    assert result["n_genes"] == _N_GENES
    assert result["method"] == LOOKUP_RULE


def test_a_single_sample_still_publishes_with_its_count_named() -> None:
    """The schema permits `n_samples` = 1, and `Placement` must be able to say so.

    Also confirms `placement` forwards a null band through to the payload
    rather than defaulting it: `q25`/`q75` are null below `n_samples` = 3 in
    the real mirror, and n=1 is the sharpest case of that.
    """
    result = placement(
        median=2.0,
        q25=None,
        q75=None,
        breakpoints=GRID,
        floor=1.0,
        unit=_UNIT,
        n_samples=1,
        n_genes=_N_GENES,
    )
    assert result is not None
    assert result["n_samples"] == 1
    assert result["q25_percentile"] is None
    assert result["q75_percentile"] is None
    assert result["median_percentile"] == 65  # bisect_left(GRID, 2.0), verified by hand


# --- tau (Yanai et al. 2005), plan Task 10 --------------------------------
#
# Every literal below was checked by independent arithmetic (`.venv/bin/
# python`, `math.log2`, a scratch script that does not import this module)
# before being pinned here:
#
#   heart 100, six organs at 10 (10-fold)  linear tau 0.900   log2 tau 0.480426
#   heart 100, six organs at 25 (4-fold)   linear tau 0.750   log2 tau 0.294039
#   heart  50, six organs at  2 (25-fold)  linear tau 0.960   log2 tau 0.720585
#
# The argmax fixture (heart 20 / liver 200 / five organs at 5) measures tau
# 0.9625 on LINEAR and 0.622772 on log2 -- "tau 0.963" is the linear figure,
# not this module's. No test below asserts a tau value for that fixture,
# because the point it makes is scale-*independent*: log2 is a strictly
# increasing transform, so it can move tau's magnitude but never its argmax.
# Recorded here so the 0.963-vs-0.623 gap is not rediscovered as a bug.


def test_tau_is_computed_on_log2_and_the_scale_is_published() -> None:
    """Measured: the scale decides the verdict, so it is part of the method.

    A 10-fold enriched gene (heart 100, six organs at 10) reads
    "tissue-specific" on linear RPKM (tau 0.900) and "broadly expressed" on
    log2 (tau 0.480). Yanai and the Kryuchkova-Mostacci & Robinson-Rechavi
    benchmark both compute tau on log expression -- computing on linear would
    call ordinary genes heart-preferential across the board.
    """
    medians = {"Heart": 100.0} | {
        organ: 10.0 for organ in ("Brain", "Cerebellum", "Kidney", "Liver", "Ovary", "Testis")
    }
    result = specificity(medians, floor=1.0)
    assert result is not None
    assert result["scale"] == "log2(x+1)"
    assert round(result["tau"], 3) == 0.480  # NOT 0.900


def test_tau_names_the_organ_it_peaked_in_because_tau_alone_cannot() -> None:
    """Tau measures concentration, not location, so a "heart-preferential"
    gloss keyed on tau alone would state the opposite of the truth here: the
    gene peaks in liver. Argmax does not depend on the scale -- log2 is
    monotonic -- so this fixture makes the point on whichever scale is used.
    """
    medians = {
        "Heart": 20.0,
        "Liver": 200.0,
        "Brain": 5.0,
        "Cerebellum": 5.0,
        "Kidney": 5.0,
        "Ovary": 5.0,
        "Testis": 5.0,
    }
    result = specificity(medians, floor=1.0)
    assert result is not None
    assert result["highest_in"] == "Liver"
    assert result["highest_in"] not in {"Heart"}


@pytest.mark.parametrize(
    ("medians", "why"),
    [
        ({"Heart": 10.0}, "one organ sampled: tau divides by n-1"),
        ({"Heart": 0.0, "Liver": 0.0, "Brain": 0.0}, "peak below floor: tau divides by max"),
    ],
)
def test_tau_is_null_never_zero_when_it_is_undefined(medians: dict[str, float], why: str) -> None:
    """Measured: a genuinely ubiquitous gene scores exactly 0.000, so a guard
    returning 0.0 for an undefined tau publishes the opposite claim -- and it
    would do so on the earliest stages, where the organ panel is smallest.
    """
    assert specificity(medians, floor=1.0) is None, why

    ubiquitous = {organ: 10.0 for organ in ("Heart", "Brain", "Liver")}
    result = specificity(ubiquitous, floor=1.0)
    assert result is not None and result["tau"] == 0.0


def test_tau_publishes_the_tissue_list_not_only_the_count() -> None:
    """A count cannot distinguish two different four-organ panels, and the
    real panel contains a correlated pair (cerebrum and cerebellum are both
    CNS). Insertion order here (Heart, Brain, Liver) is deliberately not
    alphabetical, so a mutant dropping `sorted()` would return a different
    tuple than the one asserted.
    """
    result = specificity({"Heart": 100.0, "Brain": 10.0, "Liver": 5.0}, floor=1.0)
    assert result is not None
    assert result["tissues"] == ("Brain", "Heart", "Liver")  # sorted, published
    assert result["n_tissues"] == 3


def test_a_below_floor_organ_is_retained_as_zero_not_dropped() -> None:
    """Dropping silently redefines n_tissues from "sampled" to "detected" --
    and it is self-defeating: the most heart-exclusive gene in the atlas,
    detected in heart alone, would fall to n=1 and get no tau at all.
    """
    result = specificity({"Heart": 100.0, "Liver": 0.0, "Brain": 0.0}, floor=1.0)
    assert result is not None
    assert result["n_tissues"] == 3
    assert result["tissues"] == ("Brain", "Heart", "Liver")
    assert result["tau"] > 0.9


def test_a_below_floor_nonzero_value_uses_its_raw_median_not_a_floored_zero() -> None:
    """The test above cannot distinguish "raw" from "floored to zero",
    because its below-floor organs are already raw 0.0 either way. Liver
    here is 0.5 -- below floor=1.0, but not itself zero, so the two readings
    diverge and this is the only test in this file that can tell them apart.

    An earlier version of this function floored a below-floor organ to zero
    before the log transform, reasoning that a sub-floor value is noise.
    That broke D39(b): `medians` is published raw, so a consumer recomputing
    tau from the published inputs got a *different* number than the one
    published beside them. Reverted; see `specificity`'s docstring.

    Verified independently (`math.log2`, not this module): log2(101) =
    6.658211, log2(11) = 3.459432, log2(1.5) = 0.584963 (Liver's raw
    contribution, UNFLOORED):

        tau = ((1 - 6.658211/6.658211) + (1 - 3.459432/6.658211) + (1 - 0.584963/6.658211)) / 2
            = (0 + 0.480401 + 0.912142) / 2 = 0.696285

    against 0.740213 if Liver had been floored to zero instead (its log2(1)
    = 0 contributes the full 1.0 term rather than 0.912142) -- a real
    difference, not a cosmetic one, and floored is the *wrong* answer here.
    """
    result = specificity({"Liver": 0.5, "Heart": 100.0, "Kidney": 10.0}, floor=1.0)
    assert result is not None
    assert round(result["tau"], 3) == 0.696  # NOT 0.740 -- that was the bug
    # `medians` names the true measured value, and tau's own arithmetic used
    # that same raw value -- D39(b) needs both halves of this to agree.
    assert result["medians"] == {"Liver": 0.5, "Heart": 100.0, "Kidney": 10.0}


def test_a_negative_median_is_clamped_to_zero_before_log2() -> None:
    """The one rewrite tau's arithmetic makes to any organ's raw median: a
    negative abundance is not a value, and `log2` of one is not a number.
    Unlike a below-floor value, this is not a judgement call about noise --
    without the clamp, `math.log2(-5.0 + 1.0)` raises `ValueError: math
    domain error`, so the fixture also proves the guard doesn't crash.

    Verified independently: log2(101) = 6.658211, and a clamped Liver
    contributes log2(0+1) = 0, the same as the fully-concentrated case:
    tau = ((1 - 6.658211/6.658211) + (1 - 0/6.658211)) / 1 = 1.0.
    """
    result = specificity({"Heart": 100.0, "Liver": -5.0}, floor=1.0)
    assert result is not None
    assert result["tau"] == 1.0


def test_tau_is_re_derivable_from_its_own_published_medians() -> None:
    """D39(b)'s whole point, made a checked invariant rather than a claim in
    a docstring: a consumer re-derives tau from what the payload publishes,
    not from this module's internals. This is the test that would have
    caught the floor-vs-raw mismatch directly -- publish raw medians while
    computing on floored ones, and recomputing tau from the published
    medians (independently, by hand, using only the published `scale`'s
    formula) disagrees with the published `tau`.

    Deliberately reuses the fixture above, where floored-vs-raw actually
    diverge (0.740 vs 0.696) -- re-deriving from a fixture where they
    coincide would not have caught the original bug either.
    """
    medians = {"Liver": 0.5, "Heart": 100.0, "Kidney": 10.0}
    result = specificity(medians, floor=1.0)
    assert result is not None

    values = [math.log2(v + 1.0) for v in result["medians"].values()]
    x_max = max(values)
    rederived = sum(1.0 - v / x_max for v in values) / (len(values) - 1)
    assert rederived == pytest.approx(result["tau"])


@pytest.mark.parametrize(
    ("peak", "floor", "expect_computed"),
    [
        (0.999, 1.0, False),  # just below the floor
        (1.0, 1.0, True),  # exactly at the floor: the pinned boundary
        (2.0, 1.0, True),  # clearly above
    ],
)
def test_the_floor_gate_on_tau_is_strict_less_than(
    peak: float, floor: float, expect_computed: bool
) -> None:
    """Task 9 pins `median < floor` for `placement`
    (`test_the_floor_gate_is_strict_less_than`); `specificity` gates its peak
    with the same comparison, so one gene cannot read "detected" under one
    figure and "below the floor" under the other from the same value.
    """
    result = specificity({"Heart": peak, "Liver": 0.1}, floor=floor)
    assert (result is not None) == expect_computed


def test_a_tie_at_the_peak_publishes_no_argmax() -> None:
    """ "Highest in X" would be an arbitrary choice between tied organs, so a
    tie publishes `None` instead. n_tissues is 4 here, deliberately different
    from the 3 used by several tests above, so a mutant hardcoding
    `n_tissues = 3` cannot survive both. The tie is between Heart and Liver,
    not the alphabetically-first organ (Brain), so a mutant that fell back to
    `leaders[0]` would return "Heart", not `None` -- still a wrong,
    detectable answer.
    """
    medians = {"Heart": 50.0, "Liver": 50.0, "Brain": 5.0, "Kidney": 5.0}
    result = specificity(medians, floor=1.0)
    assert result is not None
    assert result["highest_in"] is None
    assert result["n_tissues"] == 4
    assert round(result["tau"], 3) == 0.363


def test_tau_names_its_method() -> None:
    """D39(a): the formula travels with the number, the same discipline
    `LOOKUP_RULE` enforces for the percentile side.
    """
    result = specificity({"Heart": 100.0, "Liver": 10.0}, floor=1.0)
    assert result is not None
    assert result["method"] == TAU_METHOD


def test_a_degenerate_zero_floor_refuses_rather_than_dividing_by_zero() -> None:
    """floor=0.0 is not a realistic detection floor -- every real one is a
    positive RPKM/TPM threshold -- but a curated value of exactly 0 must
    still refuse cleanly rather than crash the build over one gene: with
    peak_raw=0.0, `peak_raw < floor` is `0.0 < 0.0`, False under the strict
    boundary convention, so the peak-below-floor gate does not fire, and
    `x_max = log2(0+1) = 0` would otherwise divide every term by zero.
    """
    assert specificity({"Heart": 0.0, "Liver": 0.0}, floor=0.0) is None


# --- phase assignment (Task 11) -------------------------------------------
#
# `_STAGES` deliberately spans all four `PhaseOutcome`s and both directions of
# "outside": "2wpc" is before every declared phase, "7wpc" matches septation,
# "20wpc" is past every declared phase (the common case per the source's own
# series, not an edge case), and "senior" is post-natal. "unknown" is
# deliberately absent from `_STAGES` so `assign_phase` is asked about a token
# with no `Stage` behind it at all.
_STAGES = (
    Stage(token="2wpc", wpc=2.0),
    Stage(token="7wpc", wpc=7.0),
    Stage(token="20wpc", wpc=20.0),
    Stage(token="senior", wpc=None),
)
_PHASES = CardiacPhaseFile(
    attributed_to="O'Rahilly & Muller 1987",
    citation="ISBN:0872796248",
    phases=[
        CardiacPhase(id="looping", label="Cardiac looping", start_wpc=3.0, end_wpc=5.0),
        CardiacPhase(id="septation", label="Septation", start_wpc=5.0, end_wpc=8.0),
    ],
)


def test_a_stage_outside_every_phase_publishes_a_reason_not_a_gap() -> None:
    """A wpc in no phase, and a null wpc, must both reach the page.

    PRF006 catches gaps *between* declared ranges; it cannot catch a stage
    past the last one. The source's series runs well past the morphogenetic
    window, so this is the common case, not an edge case.
    """
    assert assign_phase("7wpc", _STAGES, _PHASES).phase_id == "septation"
    assert assign_phase("20wpc", _STAGES, _PHASES).phase_id is None
    assert assign_phase("20wpc", _STAGES, _PHASES).reason == "after the curated window"
    assert assign_phase("senior", _STAGES, _PHASES).reason == "post-natal"
    assert assign_phase("unknown", _STAGES, _PHASES).reason == "stage not declared by this dataset"


def test_a_matched_phase_carries_no_reason_the_page_would_have_to_suppress() -> None:
    """`reason` explains an ABSENT `phase_id`, and is `None` exactly when
    `phase_id` is not -- a caller must never face both populated at once.
    """
    result = assign_phase("7wpc", _STAGES, _PHASES)
    assert result.outcome is PhaseOutcome.MATCHED
    assert result.phase_id == "septation"
    assert result.reason is None


@pytest.mark.parametrize(
    ("token", "expected_outcome"),
    [
        pytest.param("senior", PhaseOutcome.POST_NATAL, id="post-natal"),
        pytest.param("unknown", PhaseOutcome.UNDECLARED, id="undeclared"),
    ],
)
def test_post_natal_and_undeclared_carry_the_matching_outcome(
    token: str, expected_outcome: PhaseOutcome
) -> None:
    """`reason` is display text; `outcome` is what a caller should branch on
    (see `PhaseOutcome`'s docstring) -- pinned separately from the `reason`
    strings above so a mutant swapping the enum value while leaving the text
    intact is still caught.
    """
    result = assign_phase(token, _STAGES, _PHASES)
    assert result.outcome is expected_outcome
    assert result.phase_id is None


def test_a_stage_before_the_curated_window_is_distinguished_from_one_after_it() -> None:
    """Both are OUTSIDE_WINDOW with a real wpc and a populated vocabulary;
    only the direction of the miss differs, and the pinned test above only
    exercises the 'after' half. A mutant collapsing 'before' into 'after'
    (or the reverse) survives that test and must fail here.
    """
    before = assign_phase("2wpc", _STAGES, _PHASES)
    after = assign_phase("20wpc", _STAGES, _PHASES)
    assert before.outcome is PhaseOutcome.OUTSIDE_WINDOW
    assert after.outcome is PhaseOutcome.OUTSIDE_WINDOW
    assert before.reason == "before the curated window"
    assert after.reason == "after the curated window"


def test_a_stage_exactly_at_the_last_phase_boundary_reads_as_after() -> None:
    """Intervals are half-open, so a wpc exactly at the last phase's own
    `end_wpc` is not inside it (`phase_for(8.0) is None` -- pinned by
    `models/phases.py`'s own boundary test) -- but the fallback comparison
    here must still be `>=`, not `>`, or this exact value falls through to
    the generic "outside" text instead of the more specific "after" one. The
    pinned test's own '20wpc' is far past this boundary and cannot
    distinguish the two operators.
    """
    boundary = assign_phase("8wpc", (Stage(token="8wpc", wpc=8.0),), _PHASES)
    assert boundary.outcome is PhaseOutcome.OUTSIDE_WINDOW
    assert boundary.reason == "after the curated window"


@pytest.mark.parametrize(
    "phases",
    [
        pytest.param(None, id="vocabulary-absent"),
        pytest.param(
            CardiacPhaseFile(attributed_to="x", citation="PMID:1", phases=[]),
            id="vocabulary-empty",
        ),
    ],
)
def test_an_absent_or_empty_phase_vocabulary_reads_as_outside_the_window(
    phases: CardiacPhaseFile | None,
) -> None:
    """The committed `curation/cardiac_phases.yaml` IS the second case today
    -- ships with zero phases, deliberately, until a source is verified -- so
    this is the real corpus's own state, not a hypothetical.

    Both parametrisations must be exercised, not just one: a fix that special
    -cases `phases is None` but forgets `not phases.phases` (or the reverse)
    passes half of this and fails the other. Both would otherwise index into
    `phases.phases` for the before/after boundary comparison, which raises
    `IndexError` on an empty sequence -- so this also proves the guard
    against that crash.
    """
    result = assign_phase("7wpc", _STAGES, phases)
    assert result.outcome is PhaseOutcome.OUTSIDE_WINDOW
    assert result.reason == "outside the curated window"
    assert result.phase_id is None


def test_an_interior_gap_is_outside_the_window_not_an_empty_vocabulary() -> None:
    """Distinguishes 'nothing curated at all' from 'something is curated, and
    this wpc simply is not in it' -- the two cases the test above and this one
    must not share a fixture over, per the recorded fixture-size lesson: a
    fixture whose values all share the quantity under test measures nothing.
    `phases` here is genuinely populated, so a fix that only special-cases an
    empty list must still find the right answer by searching it.
    """
    gapped = CardiacPhaseFile(
        attributed_to="O'Rahilly & Muller 1987",
        citation="ISBN:0872796248",
        phases=[
            CardiacPhase(id="looping", label="Cardiac looping", start_wpc=3.0, end_wpc=5.0),
            CardiacPhase(id="septation", label="Septation", start_wpc=6.0, end_wpc=8.0),
        ],
    )
    result = assign_phase("7wpc", (Stage(token="7wpc", wpc=5.5),), gapped)
    assert result.outcome is PhaseOutcome.OUTSIDE_WINDOW
    assert result.reason == "outside the curated window"


# --- quantile shard emission (Task 12) -------------------------------------


def _emitter(tmp_path: Path) -> Emitter:
    return Emitter(root=tmp_path)


# Deliberately NOT already in percentile order on disk. `read_table` preserves
# a TSV's row order verbatim (`pl.read_csv`, no grouping or sort of its own),
# so if this fixture were pre-sorted, a builder that forgot to sort at all
# would still pass -- the "row sort removed" mutation named in the self-review
# would be invisible. Reversed is the simplest permutation that is
# unambiguously not ascending.
_SHUFFLED_PERCENTILES = tuple(reversed(range(101)))


def _root_with_quantiles(tmp_path: Path, *, dataset: str = "E-MTAB-6814") -> Path:
    """One shard, one (tissue, stage) grid, 101 rows written out of order."""
    directory = tmp_path / "mirrors" / "profile_quantiles"
    directory.mkdir(parents=True, exist_ok=True)
    header = "dataset\ttissue\tstage\tpercentile\tvalue\tunit\tn_genes\n"
    rows = "".join(
        f"{dataset}\tHeart\t7wpc\t{percentile}\t{float(percentile)}\trpkm\t19842\n"
        for percentile in _SHUFFLED_PERCENTILES
    )
    (directory / f"{dataset}.tsv").write_text(header + rows)
    return tmp_path


def test_the_quantile_grid_is_published_so_the_percentile_can_be_checked(tmp_path: Path) -> None:
    """D39(b) is unmet unless this file is fetchable.

    `build_omics` emits a shard only for tables in `_GENE_COLUMN`, and a
    quantile table has no gene column -- so without this the grid is
    mirrored, schema-validated, sorted and checksummed, and reaches no
    published byte. The table whose entire purpose is making the headline
    number auditable would be the one thing a consumer could not fetch, and
    nothing would catch it: the validators check the mirror against itself,
    the arithmetic tests pass on unpublished input, and a `diff -rq` between
    two builds shows a file that appears, never one that should have.
    """
    emitter = _emitter(tmp_path)
    build_profile_quantiles(_root_with_quantiles(tmp_path), emitter)
    payload = json.loads((tmp_path / "omics/profile_quantiles/E-MTAB-6814.json").read_bytes())
    assert payload["table"] == "profile_quantiles"
    assert len(payload["rows"]) == 101
    assert [row["percentile"] for row in payload["rows"]] == list(range(101))


def test_a_shard_with_two_grids_keeps_each_ones_percentiles_contiguous(tmp_path: Path) -> None:
    """A single-grid fixture cannot distinguish 'sort by percentile alone'
    from 'sort by the table's own (dataset, tissue, stage, percentile) key' --
    the two rules agree whenever there is only one (tissue, stage) pair, which
    is the shape of the fixture above and NOT the shape of a real shard (many
    organs, many stages per dataset). Two three-row grids here, interleaved
    and each individually out of order, make the two rules disagree: sorting
    by percentile alone would interleave Heart's and Liver's rows by
    breakpoint value; sorting by the full key keeps each tissue's own rows
    together and internally ascending.
    """
    directory = tmp_path / "mirrors" / "profile_quantiles"
    directory.mkdir(parents=True)
    header = "dataset\ttissue\tstage\tpercentile\tvalue\tunit\tn_genes\n"
    rows = (
        "E-MTAB-6814\tLiver\t7wpc\t2\t9.0\trpkm\t19842\n"
        "E-MTAB-6814\tHeart\t7wpc\t1\t8.0\trpkm\t19842\n"
        "E-MTAB-6814\tLiver\t7wpc\t0\t7.0\trpkm\t19842\n"
        "E-MTAB-6814\tHeart\t7wpc\t2\t6.0\trpkm\t19842\n"
        "E-MTAB-6814\tLiver\t7wpc\t1\t5.0\trpkm\t19842\n"
        "E-MTAB-6814\tHeart\t7wpc\t0\t4.0\trpkm\t19842\n"
    )
    (directory / "E-MTAB-6814.tsv").write_text(header + rows)

    build_profile_quantiles(tmp_path, _emitter(tmp_path))

    payload = json.loads((tmp_path / "omics/profile_quantiles/E-MTAB-6814.json").read_bytes())
    assert [row["tissue"] for row in payload["rows"]] == ["Heart"] * 3 + ["Liver"] * 3
    assert [row["percentile"] for row in payload["rows"]] == [0, 1, 2, 0, 1, 2]


def test_a_shard_filename_needing_escape_is_slugged_before_use(tmp_path: Path) -> None:
    """The stem becomes a URL; `build_omics` slugs its own shard stems for the
    same reason (`paths.slug`'s docstring). Without the `slug` call this would
    publish `omics/profile_quantiles/bad name.json` -- a space that a browser
    or a bare `curl` invocation would need to escape before it could fetch it.

    Also pins the returned mapping's shape: keyed on the raw accession (what a
    caller already has from `profiles.dataset`/`Dataset.id`), valued with the
    exact slugged path this function wrote -- so a future caller reaches the
    file by reading this return rather than by reconstructing the path with
    a second call to `slug` that could drift from this one.
    """
    directory = tmp_path / "mirrors" / "profile_quantiles"
    directory.mkdir(parents=True)
    header = "dataset\ttissue\tstage\tpercentile\tvalue\tunit\tn_genes\n"
    (directory / "bad name.tsv").write_text(header + "bad name\tHeart\t7wpc\t0\t1.0\trpkm\t100\n")

    shards = build_profile_quantiles(tmp_path, _emitter(tmp_path))

    assert (tmp_path / "omics" / "profile_quantiles" / "bad_name.json").is_file()
    assert not (tmp_path / "omics" / "profile_quantiles" / "bad name.json").exists()
    assert shards == {"bad name": "omics/profile_quantiles/bad_name.json"}


def test_an_unreadable_shard_does_not_stop_the_others(tmp_path: Path) -> None:
    """A zero-length shard makes polars raise, and one bad file is not the build.

    `validate_table` reports the same file against the same path; raising
    here would report it twice and abort every shard behind it -- the same
    reasoning as `build_omics`'s own
    `test_an_unreadable_mirror_does_not_stop_the_others`.
    """
    directory = tmp_path / "mirrors" / "profile_quantiles"
    directory.mkdir(parents=True)
    (directory / "E-MTAB-0000001.tsv").write_text("")
    _root_with_quantiles(tmp_path, dataset="E-MTAB-0000002")

    shards = build_profile_quantiles(tmp_path, _emitter(tmp_path))

    assert list(shards) == ["E-MTAB-0000002"]
    assert (tmp_path / "omics" / "profile_quantiles" / "E-MTAB-0000002.json").is_file()
    assert not (tmp_path / "omics" / "profile_quantiles" / "E-MTAB-0000001.json").exists()


# --- gene_expression_profiles: the per-gene assembly (Task 13) -------------
#
# Four recorded recurrences in this project of a fixture whose rows all share
# the value under test, so the mutant it was meant to catch survives. This
# section is built to distinguish, in the sense of separate, dedicated rows
# rather than one fixture asserted on twice: a gene below the detection floor
# from one above it; a gene sampled in one organ from one sampled in several;
# a gene whose argmax is the heart from one whose argmax is not; n_tissues
# from a hardcoded literal (two different panel sizes); a null-stage row from
# a declared-stage row; a gene with no profile data at all.

_ASCENDING_GRID = tuple(float(percentile) for percentile in range(101))


def _profile_row(
    gene: str = "HGNC:11604",
    tissue: str = "Heart",
    stage: str | None = "7wpc",
    median: float = 50.0,
    unit: str = "rpkm",
    q25: float | None = 40.0,
    q75: float | None = 60.0,
    n_samples: int = 3,
) -> dict[str, object]:
    return {
        "gene": gene,
        "tissue": tissue,
        "stage": stage,
        "median_abundance": median,
        "unit": unit,
        "q25": q25,
        "q75": q75,
        "n_samples": n_samples,
    }


def _tsv_cell(value: object) -> str:
    return "" if value is None else str(value)


def _write_profiles(
    root: Path, dataset: str, rows: list[dict[str, object]], *, filename: str | None = None
) -> None:
    """Write one `profiles` shard, naming the file `{dataset}.tsv` by default.

    `filename` decouples the shard's name from the `dataset` *column* value it
    writes: `mirror_paths` sorts shards by filename, so a fixture that always
    lets filename and column agree can never tell "sorted by dataset id" from
    "left in the order `mirror_paths` read the files" -- the two coincide by
    construction whenever the accession names the file that carries it, which
    is the ordinary case this atlas curates but not one this test may lean on.
    """
    directory = root / "mirrors" / "profiles"
    directory.mkdir(parents=True, exist_ok=True)
    header = "dataset\tgene\ttissue\tstage\tmedian_abundance\tunit\tq25\tq75\tn_samples\n"
    body = "".join(
        "\t".join(
            _tsv_cell(value)
            for value in (
                dataset,
                row["gene"],
                row["tissue"],
                row["stage"],
                row["median_abundance"],
                row["unit"],
                row["q25"],
                row["q75"],
                row["n_samples"],
            )
        )
        + "\n"
        for row in rows
    )
    (directory / f"{filename or dataset}.tsv").write_text(header + body)


def _write_quantiles(
    root: Path,
    dataset: str,
    grids: dict[tuple[str, str], tuple[float, ...]],
    *,
    unit: str = "rpkm",
    n_genes: int = 100,
    filename: str | None = None,
) -> None:
    """One shard holding every (tissue, stage) grid supplied, each 101 rows.

    A single `.write_text` per dataset: `_write_quantiles` called twice for
    one dataset would overwrite rather than accumulate, so every grid a test
    needs for one dataset is passed in one call. `filename` -- see
    `_write_profiles`'s docstring for why a fixture needs it decoupled.
    """
    directory = root / "mirrors" / "profile_quantiles"
    directory.mkdir(parents=True, exist_ok=True)
    header = "dataset\ttissue\tstage\tpercentile\tvalue\tunit\tn_genes\n"
    body = "".join(
        f"{dataset}\t{tissue}\t{stage}\t{percentile}\t{value}\t{unit}\t{n_genes}\n"
        for (tissue, stage), values in grids.items()
        for percentile, value in enumerate(values)
    )
    (directory / f"{filename or dataset}.tsv").write_text(header + body)


def _dataset(
    accession: str = "E-MTAB-6814",
    *,
    detection_floor: float = 1.0,
    cardiac_tissues: tuple[str, ...] = ("Heart",),
    stages: tuple[Stage, ...] = (Stage(token="7wpc", wpc=7.0),),
) -> Dataset:
    return Dataset(
        id=accession,
        archive="arrayexpress",
        technology="bulk_rnaseq",
        design="profile",
        tissue="whole embryo",
        developmental_stage="embryonic",
        organism="NCBITaxon:9606",
        n_samples=3,
        licence="CC BY 4.0",
        contrasts=[],
        cardiac_tissues=cardiac_tissues,
        detection_floor=detection_floor,
        floor_source="source methods, section 4",
        quantile_estimator="linear",
        stages=stages,
    )


def _contrast_dataset(accession: str) -> Dataset:
    """A registered dataset with no declared floor.

    `a_profile_dataset_is_fully_declared` refuses to load a `design="profile"`
    dataset with no floor, so the only way to construct a *registered*
    dataset with `detection_floor is None` is one of a different design --
    standing in for a curator's typo pointing a profiles row at the wrong
    accession, which `validate_profile_references`'s own docstring names as
    out of scope for that validator (`_prf004_issues`: "a profiles row naming
    a dataset with no curated record at all is out of scope for this check").
    """
    return Dataset(
        id=accession,
        archive="geo",
        technology="bulk_rnaseq",
        design="contrast",
        tissue="Heart",
        developmental_stage="adult",
        organism="NCBITaxon:9606",
        n_samples=6,
        licence="CC BY 4.0",
        contrasts=[
            {
                "id": "c1",
                "description": "d",
                "case_group": "case",
                "control_group": "control",
                "statistical_method": "m",
                "software": "s",
            }
        ],
    )


def test_a_gene_with_no_rows_in_any_mirror_has_no_entry(tmp_path: Path) -> None:
    """The assembler reports what it read; `bundles.py` supplies the
    always-present empty shape for a gene absent from this mapping.
    """
    _write_profiles(tmp_path, "E-MTAB-6814", [_profile_row(gene="HGNC:11604")])
    _write_quantiles(tmp_path, "E-MTAB-6814", {("Heart", "7wpc"): _ASCENDING_GRID})

    result = gene_expression_profiles(tmp_path, (_dataset(),), None, {})

    assert "HGNC:11604" in result
    assert "HGNC:99999" not in result


def test_a_gene_above_the_floor_is_placed_and_one_below_it_is_not(tmp_path: Path) -> None:
    """D41's refusal must be reachable, which needs a fixture that separates
    the two rather than asserting the same gene twice.
    """
    _write_profiles(
        tmp_path,
        "E-MTAB-6814",
        [
            _profile_row(gene="HGNC:1", median=50.0),  # well above floor=1.0
            _profile_row(gene="HGNC:2", median=0.5),  # below floor=1.0, not zero
        ],
    )
    _write_quantiles(tmp_path, "E-MTAB-6814", {("Heart", "7wpc"): _ASCENDING_GRID})

    result = gene_expression_profiles(tmp_path, (_dataset(detection_floor=1.0),), None, {})

    placed = result["HGNC:1"]["datasets"][0]["stages"][0]["tissues"][0]
    unplaced = result["HGNC:2"]["datasets"][0]["stages"][0]["tissues"][0]
    assert placed["placement"] is not None
    assert placed["not_placed_reason"] is None
    assert unplaced["placement"] is None
    assert unplaced["not_placed_reason"] == "below_detection_floor"
    # The raw measurement still publishes even when it is not placed -- an
    # unplaced gene is not an absent one.
    assert unplaced["median_abundance"] == 0.5
    assert unplaced["unit"] == "rpkm"


def test_a_gene_sampled_in_one_organ_gets_no_specificity_and_one_in_several_does(
    tmp_path: Path,
) -> None:
    """Tau's `None` refusal at n<2 organs must be reachable, not just asserted."""
    _write_profiles(
        tmp_path,
        "E-MTAB-6814",
        [_profile_row(gene="HGNC:1", tissue="Heart", median=50.0)]
        + [
            _profile_row(gene="HGNC:2", tissue=tissue, median=value)
            for tissue, value in (("Heart", 50.0), ("Liver", 5.0), ("Brain", 5.0))
        ],
    )
    _write_quantiles(
        tmp_path,
        "E-MTAB-6814",
        {(tissue, "7wpc"): _ASCENDING_GRID for tissue in ("Heart", "Liver", "Brain")},
    )

    result = gene_expression_profiles(tmp_path, (_dataset(),), None, {})

    one_organ = result["HGNC:1"]["datasets"][0]["stages"][0]
    several = result["HGNC:2"]["datasets"][0]["stages"][0]
    assert one_organ["specificity"] is None
    assert one_organ["specificity_unavailable_reason"] == "one_organ_sampled"
    assert several["specificity"] is not None
    assert several["specificity_unavailable_reason"] is None


def test_the_argmax_organ_is_published_whether_or_not_it_is_the_heart(tmp_path: Path) -> None:
    """Keying a gloss on the wrong organ (or a hardcoded one) must be
    reachable to be caught -- two genes, one peaking in Heart and one in
    Liver, never the same gene asserted on twice.
    """
    _write_profiles(
        tmp_path,
        "E-MTAB-6814",
        [
            _profile_row(gene="HGNC:1", tissue=tissue, median=value)
            for tissue, value in (("Heart", 100.0), ("Liver", 10.0), ("Brain", 10.0))
        ]
        + [
            _profile_row(gene="HGNC:2", tissue=tissue, median=value)
            for tissue, value in (("Heart", 10.0), ("Liver", 200.0), ("Brain", 5.0))
        ],
    )
    _write_quantiles(
        tmp_path,
        "E-MTAB-6814",
        {(tissue, "7wpc"): _ASCENDING_GRID for tissue in ("Heart", "Liver", "Brain")},
    )

    result = gene_expression_profiles(tmp_path, (_dataset(),), None, {})

    heart_peak = result["HGNC:1"]["datasets"][0]["stages"][0]["specificity"]
    liver_peak = result["HGNC:2"]["datasets"][0]["stages"][0]["specificity"]
    assert heart_peak is not None and heart_peak["highest_in"] == "Heart"
    assert liver_peak is not None and liver_peak["highest_in"] == "Liver"


def test_n_tissues_reflects_the_real_panel_size_not_a_hardcoded_literal(tmp_path: Path) -> None:
    """Two different panel sizes, so a hardcoded `n_tissues` cannot survive both."""
    _write_profiles(
        tmp_path,
        "E-MTAB-6814",
        [_profile_row(gene="HGNC:1", tissue=tissue, median=10.0) for tissue in ("Heart", "Liver")]
        + [
            _profile_row(gene="HGNC:2", tissue=tissue, median=10.0)
            for tissue in ("Heart", "Liver", "Brain", "Kidney", "Ovary")
        ],
    )
    _write_quantiles(
        tmp_path,
        "E-MTAB-6814",
        {
            (tissue, "7wpc"): _ASCENDING_GRID
            for tissue in ("Heart", "Liver", "Brain", "Kidney", "Ovary")
        },
    )

    result = gene_expression_profiles(tmp_path, (_dataset(),), None, {})

    small = result["HGNC:1"]["datasets"][0]["stages"][0]["specificity"]
    large = result["HGNC:2"]["datasets"][0]["stages"][0]["specificity"]
    assert small is not None and small["n_tissues"] == 2
    assert large is not None and large["n_tissues"] == 5


def test_a_null_stage_row_publishes_with_a_stated_reason_not_a_dropped_row(
    tmp_path: Path,
) -> None:
    """A null-stage row can have no quantile partner by construction
    (`profile_quantiles.stage` is never null). The join keeps the row; only
    the PRF002 curation check is what skips it -- this is that same design
    fact, reached from the assembly side.
    """
    _write_profiles(
        tmp_path,
        "E-MTAB-6814",
        [
            _profile_row(gene="HGNC:1", tissue="Heart", stage=None, median=50.0),
            _profile_row(gene="HGNC:1", tissue="Heart", stage="7wpc", median=50.0),
        ],
    )
    _write_quantiles(tmp_path, "E-MTAB-6814", {("Heart", "7wpc"): _ASCENDING_GRID})

    result = gene_expression_profiles(tmp_path, (_dataset(),), None, {})

    stages = {entry["stage"]: entry for entry in result["HGNC:1"]["datasets"][0]["stages"]}
    assert set(stages) == {None, "7wpc"}  # the null-stage row is NOT dropped

    null_stage_tissue = stages[None]["tissues"][0]
    assert null_stage_tissue["placement"] is None
    assert null_stage_tissue["not_placed_reason"] == "no_quantile_grid"
    # The raw measurement is still published even with no percentile.
    assert null_stage_tissue["median_abundance"] == 50.0

    declared_stage_tissue = stages["7wpc"]["tissues"][0]
    assert declared_stage_tissue["placement"] is not None

    # A null stage has no token to place in the developmental window at all --
    # distinct from `PhaseOutcome.UNDECLARED`, which means a real token the
    # dataset simply never declared.
    assert stages[None]["phase"] == {
        "outcome": None,
        "phase_id": None,
        "reason": "no developmental stage recorded for this measurement",
    }


def test_a_dataset_absent_from_the_curated_registry_still_publishes_its_rows(
    tmp_path: Path,
) -> None:
    """Must not crash, must not silently vanish -- `_prf004_issues`'s own
    docstring says an unregistered dataset accession is out of scope for
    `validate_profile_references`, so this is reachable on a validated
    repository, not only on a bypassed gate.
    """
    _write_profiles(tmp_path, "E-MTAB-9999", [_profile_row(gene="HGNC:1", median=50.0)])
    _write_quantiles(tmp_path, "E-MTAB-9999", {("Heart", "7wpc"): _ASCENDING_GRID})

    result = gene_expression_profiles(tmp_path, (), None, {})  # no Dataset record at all

    stage = result["HGNC:1"]["datasets"][0]["stages"][0]
    tissue = stage["tissues"][0]
    assert tissue["placement"] is None
    assert tissue["not_placed_reason"] == "dataset_not_registered"
    assert stage["specificity_unavailable_reason"] == "dataset_not_registered"
    # The phase can still be asked about the token with no dataset record at
    # all, and correctly finds nothing declared.
    assert stage["phase"]["outcome"] == "undeclared"


def test_a_registered_dataset_with_no_declared_floor_still_publishes_its_rows(
    tmp_path: Path,
) -> None:
    """A `design="contrast"` dataset naturally has `detection_floor=None` --
    pydantic refuses a `design="profile"` dataset with no floor, so this is
    the only way to construct a *registered* dataset this gap can reach --
    standing in for a curator's typo referencing the wrong accession.
    """
    _write_profiles(tmp_path, "GSE1000", [_profile_row(gene="HGNC:1", median=50.0)])
    _write_quantiles(tmp_path, "GSE1000", {("Heart", "7wpc"): _ASCENDING_GRID})

    result = gene_expression_profiles(tmp_path, (_contrast_dataset("GSE1000"),), None, {})

    tissue = result["HGNC:1"]["datasets"][0]["stages"][0]["tissues"][0]
    assert tissue["not_placed_reason"] == "detection_floor_undeclared"


def test_an_incomplete_quantile_grid_is_treated_as_unavailable(tmp_path: Path) -> None:
    """`percentile_of` assumes position i holds percentile i's own value
    (Task 9's docstring); a grid missing even one index would silently
    misplace a gene rather than merely cost one point of precision, and no
    PRF check pins grid *completeness* today -- PRF001-009 check units, cell
    presence and monotonicity, never the count of rows in one cell.
    """
    _write_profiles(tmp_path, "E-MTAB-6814", [_profile_row(gene="HGNC:1", median=50.0)])
    directory = tmp_path / "mirrors" / "profile_quantiles"
    directory.mkdir(parents=True)
    header = "dataset\ttissue\tstage\tpercentile\tvalue\tunit\tn_genes\n"
    # Percentile 50 is missing from an otherwise complete grid.
    rows = "".join(
        f"E-MTAB-6814\tHeart\t7wpc\t{percentile}\t{float(percentile)}\trpkm\t100\n"
        for percentile in range(101)
        if percentile != 50
    )
    (directory / "E-MTAB-6814.tsv").write_text(header + rows)

    result = gene_expression_profiles(tmp_path, (_dataset(),), None, {})

    tissue = result["HGNC:1"]["datasets"][0]["stages"][0]["tissues"][0]
    assert tissue["placement"] is None
    assert tissue["not_placed_reason"] == "no_quantile_grid"


def test_placement_carries_the_grids_own_n_genes(tmp_path: Path) -> None:
    """`n_genes` must come from the grid actually read, not a hardcoded
    literal shared with every other test in this file.
    """
    _write_profiles(tmp_path, "E-MTAB-6814", [_profile_row(gene="HGNC:1", median=50.0)])
    _write_quantiles(tmp_path, "E-MTAB-6814", {("Heart", "7wpc"): _ASCENDING_GRID}, n_genes=19842)

    result = gene_expression_profiles(tmp_path, (_dataset(),), None, {})

    placed = result["HGNC:1"]["datasets"][0]["stages"][0]["tissues"][0]["placement"]
    assert placed is not None
    assert placed["n_genes"] == 19842


def test_the_dataset_entry_carries_its_own_quantile_shard_link(tmp_path: Path) -> None:
    """`build_profile_quantiles`'s own return threads through unchanged, so a
    consumer reaches the grid its percentile came from without a second call
    to `slug` that could drift from the one that wrote the file.
    """
    _write_profiles(tmp_path, "E-MTAB-6814", [_profile_row(gene="HGNC:1")])
    _write_quantiles(tmp_path, "E-MTAB-6814", {("Heart", "7wpc"): _ASCENDING_GRID})

    linked = gene_expression_profiles(
        tmp_path,
        (_dataset(),),
        None,
        {"E-MTAB-6814": "omics/profile_quantiles/E-MTAB-6814.json"},
    )
    unlinked = gene_expression_profiles(tmp_path, (_dataset(),), None, {})

    assert (
        linked["HGNC:1"]["datasets"][0]["quantile_shard"]
        == "omics/profile_quantiles/E-MTAB-6814.json"
    )
    assert unlinked["HGNC:1"]["datasets"][0]["quantile_shard"] is None


def test_a_real_stage_resolves_its_phase_through_the_datasets_own_wpc_mapping(
    tmp_path: Path,
) -> None:
    """Wires `assign_phase` against a populated vocabulary, not only the
    placeholder empty one every other test in this section uses.
    """
    _write_profiles(tmp_path, "E-MTAB-6814", [_profile_row(gene="HGNC:1", stage="7wpc")])
    _write_quantiles(tmp_path, "E-MTAB-6814", {("Heart", "7wpc"): _ASCENDING_GRID})
    phases = CardiacPhaseFile(
        attributed_to="O'Rahilly & Muller 1987",
        citation="ISBN:0872796248",
        phases=[CardiacPhase(id="septation", label="Septation", start_wpc=5.0, end_wpc=8.0)],
    )

    result = gene_expression_profiles(tmp_path, (_dataset(),), phases, {})

    phase = result["HGNC:1"]["datasets"][0]["stages"][0]["phase"]
    assert phase == {"outcome": "matched", "phase_id": "septation", "reason": None}


def test_tau_is_computed_per_stage_not_merged_across_stages(tmp_path: Path) -> None:
    """A mutation that grouped tau's inputs by tissue instead of by stage
    would either collapse the two stages into one bucket per tissue name or
    mix one stage's organs with another's -- either way, the two stages'
    specificities could not independently reproduce what an outside call to
    `specificity()` gives the same per-stage medians.
    """
    _write_profiles(
        tmp_path,
        "E-MTAB-6814",
        [
            _profile_row(gene="HGNC:1", tissue="Heart", stage="s1", median=100.0),
            _profile_row(gene="HGNC:1", tissue="Liver", stage="s1", median=10.0),
            _profile_row(gene="HGNC:1", tissue="Brain", stage="s1", median=10.0),
            _profile_row(gene="HGNC:1", tissue="Heart", stage="s2", median=10.0),
            _profile_row(gene="HGNC:1", tissue="Liver", stage="s2", median=10.0),
            _profile_row(gene="HGNC:1", tissue="Brain", stage="s2", median=100.0),
        ],
    )
    _write_quantiles(
        tmp_path,
        "E-MTAB-6814",
        {
            ("Heart", "s1"): _ASCENDING_GRID,
            ("Liver", "s1"): _ASCENDING_GRID,
            ("Brain", "s1"): _ASCENDING_GRID,
            ("Heart", "s2"): _ASCENDING_GRID,
            ("Liver", "s2"): _ASCENDING_GRID,
            ("Brain", "s2"): _ASCENDING_GRID,
        },
    )
    dataset = _dataset(
        stages=(Stage(token="s1", wpc=1.0), Stage(token="s2", wpc=2.0)),
    )

    result = gene_expression_profiles(tmp_path, (dataset,), None, {})

    stages = {entry["stage"]: entry for entry in result["HGNC:1"]["datasets"][0]["stages"]}
    s1, s2 = stages["s1"]["specificity"], stages["s2"]["specificity"]
    assert s1 is not None and s2 is not None
    assert s1["highest_in"] == "Heart"
    assert s2["highest_in"] == "Brain"
    # Independently re-derived from the same per-stage medians, using the
    # already-tested `specificity()` directly -- not a magic literal.
    assert s1 == specificity({"Heart": 100.0, "Liver": 10.0, "Brain": 10.0}, floor=1.0)
    assert s2 == specificity({"Heart": 10.0, "Liver": 10.0, "Brain": 100.0}, floor=1.0)


def test_datasets_stages_and_tissues_are_all_sorted(tmp_path: Path) -> None:
    """A dropped sort anywhere in the nesting must be visible.

    Every dataset id, stage token and tissue name below is written to disk
    out of its eventual sort order, matching this project's own precedent
    (`_SHUFFLED_PERCENTILES` above) for why a pre-sorted fixture cannot catch
    a builder that merely preserves encounter order.
    """
    # `E-ZZZZ-9`/`E-AAAA-1`, not literal "Z-LATER"/"A-EARLIER": `Dataset.id` is
    # an `AccessionId`, pattern-checked against ArrayExpress's own grammar.
    #
    # The shard *filenames* are deliberately the reverse of the dataset ids'
    # own alphabetical order: `mirror_paths` sorts shards by filename, so a
    # fixture that lets a dataset's filename and its `dataset` column agree
    # (the ordinary case, and every other fixture in this file) can never
    # distinguish "sorted by dataset id" from "left in the order the mirror
    # was read" -- the two always coincide when the id names the file. Naming
    # the file that carries "E-ZZZZ-9" `aaa-reads-first.tsv` and the file
    # carrying "E-AAAA-1" `zzz-reads-second.tsv` makes read order and sorted
    # order disagree, so a dropped `sorted(by_dataset)` publishes
    # ["E-ZZZZ-9", "E-AAAA-1"] instead of the correct
    # ["E-AAAA-1", "E-ZZZZ-9"].
    _write_profiles(
        tmp_path,
        "E-ZZZZ-9",
        [_profile_row(gene="HGNC:1", tissue="Zebra", stage="s1")],
        filename="aaa-reads-first",
    )
    _write_profiles(
        tmp_path,
        "E-AAAA-1",
        [
            _profile_row(gene="HGNC:1", tissue="Zebra", stage="s3"),
            _profile_row(gene="HGNC:1", tissue="Alpha", stage="s3"),
            _profile_row(gene="HGNC:1", tissue="Zebra", stage="s1"),
            _profile_row(gene="HGNC:1", tissue="Zebra", stage="s2"),
        ],
        filename="zzz-reads-second",
    )
    _write_quantiles(
        tmp_path,
        "E-ZZZZ-9",
        {("Zebra", "s1"): _ASCENDING_GRID},
        filename="aaa-reads-first",
    )
    _write_quantiles(
        tmp_path,
        "E-AAAA-1",
        {
            ("Zebra", "s3"): _ASCENDING_GRID,
            ("Alpha", "s3"): _ASCENDING_GRID,
            ("Zebra", "s1"): _ASCENDING_GRID,
            ("Zebra", "s2"): _ASCENDING_GRID,
        },
        filename="zzz-reads-second",
    )

    result = gene_expression_profiles(
        tmp_path,
        (
            _dataset("E-ZZZZ-9", stages=(Stage(token="s1", wpc=1.0),)),
            _dataset(
                "E-AAAA-1",
                stages=(
                    Stage(token="s3", wpc=3.0),
                    Stage(token="s1", wpc=1.0),
                    Stage(token="s2", wpc=2.0),
                ),
            ),
        ),
        None,
        {},
    )

    datasets = result["HGNC:1"]["datasets"]
    assert [entry["dataset"] for entry in datasets] == ["E-AAAA-1", "E-ZZZZ-9"]

    stages = datasets[0]["stages"]
    assert [entry["stage"] for entry in stages] == ["s1", "s2", "s3"]

    tissues = next(entry for entry in stages if entry["stage"] == "s3")["tissues"]
    assert [entry["tissue"] for entry in tissues] == ["Alpha", "Zebra"]


# --- percentile_annotations: the flat lookup build_omics ranks on (Task 8b) -


def test_percentile_annotations_reads_back_every_placed_cell_and_only_those(
    tmp_path: Path,
) -> None:
    """The one flat read of the `Placement`s `gene_expression_profiles`
    already computed -- never a second derivation of a percentile.

    HGNC:1 is above the floor and gets an entry keyed on the exact cell it was
    measured in; HGNC:2 is below it and gets none, the same way `placement`
    itself returns `None` rather than a number that would mislead a reader.
    Asserting the whole returned dict, rather than checking HGNC:1's entry in
    isolation, is what proves the second gene contributes nothing at all.
    """
    _write_profiles(
        tmp_path,
        "E-MTAB-6814",
        [
            _profile_row(gene="HGNC:1", tissue="Heart", stage="7wpc", median=42.0),
            _profile_row(gene="HGNC:2", tissue="Heart", stage="7wpc", median=0.5),
        ],
    )
    _write_quantiles(tmp_path, "E-MTAB-6814", {("Heart", "7wpc"): _ASCENDING_GRID})

    profiles = gene_expression_profiles(tmp_path, (_dataset(detection_floor=1.0),), None, {})
    annotations = percentile_annotations(profiles)

    assert annotations == {("HGNC:1", "E-MTAB-6814", "Heart", "7wpc"): 42}


def test_build_omics_and_the_published_bundle_rank_on_the_same_percentile(
    tmp_path: Path,
) -> None:
    """`build_omics` must never compute its own percentile.

    Both figures below trace to the one call to `gene_expression_profiles`:
    `percentile_annotations` only reads `median_percentile` back out of the
    `Placement`s that call already returned, the same object `build_genes`
    would go on to publish as the gene's `expression_profile`. Stage tokens
    are chosen so alphabetical order and rank order disagree ("s1" sorts
    first, "s9" ranks first) -- agreeing fixtures pass whether or not the
    annotation is wired up at all, the trap recorded at §4.14/15b/30/36.

    This reproduces `runner.py`'s own call sequence (`build_profile_quantiles`
    then `gene_expression_profiles` then `percentile_annotations` then
    `build_omics`, one shared `Emitter`) rather than `build_site`'s YAML
    corpus, so it is cheap while still exercising the real wiring between the
    two modules -- the one place a key-shape mismatch between the producer and
    the consumer would actually show up.
    """
    dataset = _dataset(stages=(Stage(token="s1", wpc=1.0), Stage(token="s9", wpc=9.0)))
    _write_profiles(
        tmp_path,
        "E-MTAB-6814",
        [
            _profile_row(gene="HGNC:11604", tissue="Heart", stage="s1", median=10.0),
            _profile_row(gene="HGNC:11604", tissue="Heart", stage="s9", median=90.0),
        ],
    )
    _write_quantiles(
        tmp_path,
        "E-MTAB-6814",
        {("Heart", "s1"): _ASCENDING_GRID, ("Heart", "s9"): _ASCENDING_GRID},
    )
    emitter = Emitter(root=tmp_path / "dist")

    quantile_shards = build_profile_quantiles(tmp_path, emitter)
    profiles = gene_expression_profiles(tmp_path, (dataset,), None, quantile_shards)
    percentiles = percentile_annotations(profiles)
    omics = build_omics(
        tmp_path,
        emitter,
        cardiac={"E-MTAB-6814": frozenset({"Heart"})},
        percentiles=percentiles,
    )

    stages = {entry["stage"]: entry for entry in profiles["HGNC:11604"]["datasets"][0]["stages"]}
    s1_placement = stages["s1"]["tissues"][0]["placement"]
    s9_placement = stages["s9"]["tissues"][0]["placement"]
    assert s1_placement is not None and s9_placement is not None
    assert (s1_placement["median_percentile"], s9_placement["median_percentile"]) == (10, 90)

    top = omics["HGNC:11604"]["profiles"]["top"]
    # `build_omics` puts the row carrying the higher `median_percentile`
    # above first -- read from the same computation, not a second one -- even
    # though "s1" sorts first alphabetically.
    assert [row["stage"] for row in top] == ["s9", "s1"]
    # The raw shard carries the same annotation the summary ranked on, not a
    # second, unpublished copy of it.
    shard = json.loads((tmp_path / "dist" / "omics" / "profiles" / "E-MTAB-6814.json").read_text())
    assert {row["stage"]: row["percentile"] for row in shard["rows"]} == {"s1": 10, "s9": 90}
