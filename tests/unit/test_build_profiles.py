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

from chd_atlas.build.profiles import LOOKUP_RULE, band, percentile_of, placement

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
