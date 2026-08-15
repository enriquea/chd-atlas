# tests/unit/test_build_profiles.py
"""Tests for the percentile half of `build/profiles.py` (plan Task 9).

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

import pytest

from chd_atlas.build.profiles import (
    LOOKUP_RULE,
    TAU_METHOD,
    band,
    percentile_of,
    placement,
    specificity,
)

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


def test_a_below_floor_nonzero_value_is_floored_to_zero_not_passed_through() -> None:
    """The test above cannot distinguish "floored to zero" from "left as
    measured", because its below-floor organs are already raw 0.0 either way.
    Liver here is 0.5 -- below floor=1.0, but not itself zero.

    Verified independently (`math.log2`, not this module): log2(101) =
    6.658211, log2(11) = 3.459432, and a *floored* Liver contributes
    log2(0+1) = 0 to the sum:

        tau = ((1 - 6.658211/6.658211) + (1 - 3.459432/6.658211) + (1 - 0/6.658211)) / 2
            = (0 + 0.480401 + 1) / 2 = 0.740213

    against 0.696285 if the raw 0.5 (log2(1.5) = 0.584963) had been used
    unfloored instead -- a real difference, not a cosmetic one.
    """
    result = specificity({"Liver": 0.5, "Heart": 100.0, "Kidney": 10.0}, floor=1.0)
    assert result is not None
    assert round(result["tau"], 3) == 0.740
    # `medians` still names the true measured value: D39(b) requires tau's
    # actual inputs, and 0.5 is what was measured, not what tau computed with.
    assert result["medians"] == {"Liver": 0.5, "Heart": 100.0, "Kidney": 10.0}


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
