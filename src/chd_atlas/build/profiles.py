# src/chd_atlas/build/profiles.py
"""The figures this atlas derives from a developmental transcriptome.

Neither figure is a mirror column. The build derives both from the two
published tables, so a consumer's re-derivation and the build's are the same
arithmetic over the same bytes -- the pattern `omics.py` uses for `genes`,
where `count` is built from the attribution the rows publish so the two cannot
drift. Emitting them as columns instead would mean two computations of each
number with no way to tell which is right.

D39(b) is a one-level guarantee, not an unbounded one: a gene's percentile is
re-derivable from its median and the published breakpoints; the breakpoints
are not, because D32 forbids publishing the matrix they came from. Stated as a
trade rather than claimed as a proof.

This module carries the percentile band only (plan Task 9). Task 10 adds tau
to this same module; Task 11 phase assignment; Task 12 quantile shard
emission.
"""

from __future__ import annotations

import bisect
from collections.abc import Sequence
from typing import Final, TypedDict

# D39(a): the method travels with the number. Published beside every
# percentile so a consumer can re-derive it without reading this module's
# source to learn which side of a tied breakpoint was chosen.
LOOKUP_RULE: Final = "lowest percentile of a tied breakpoint (bisect_left)"


def percentile_of(value: float, breakpoints: Sequence[float]) -> int:
    """Where `value` falls in a 101-point grid, as a percentile 0-100.

    `bisect_left`, so a value tied with a run of identical breakpoints reports
    the *lowest* percentile in that run. The alternative reads the top of the
    tie, and roughly half a bulk library ties at zero -- which would publish an
    unexpressed gene above the median. Understating is the safe direction here,
    and D41 keeps those genes out of the published figure regardless (see
    `placement` below).

    Clamped at 100 rather than allowed to return 101 for a value above every
    breakpoint: the grid's top is the maximum, so exceeding it is a rounding
    artefact, not a 101st percentile. A value that merely *equals* the top
    breakpoint already reaches 100 through `bisect_left` alone -- the clamp is
    only ever exercised by a value that exceeds every breakpoint, which is a
    genuinely different input, not the same case asserted twice.
    """
    return min(bisect.bisect_left(breakpoints, value), 100)


def band(
    median: float,
    q25: float | None,
    q75: float | None,
    breakpoints: Sequence[float],
) -> tuple[int | None, int, int | None]:
    """The percentile band for one gene: (q25's percentile, median's, q75's).

    Published as a band rather than a point. At `n_samples` = 2 a point
    estimate of "top 4%" has an honest range of roughly top 2% to top 15%;
    publishing it costs one extra lookup per side and no sampling model,
    because `q25`/`q75` are already columns in `mirrors/profiles/`.

    `q25`/`q75` are `None` exactly when the mirror's are -- the schema nulls
    them below `n_samples` = 3, where a quartile of two points is not a
    quartile -- and `None` passes straight through rather than falling back to
    the median's own percentile, which would publish a band with no reported
    width as though it had been measured.
    """
    q25_percentile = percentile_of(q25, breakpoints) if q25 is not None else None
    q75_percentile = percentile_of(q75, breakpoints) if q75 is not None else None
    return (q25_percentile, percentile_of(median, breakpoints), q75_percentile)


class Placement(TypedDict):
    """One gene's percentile band in one (dataset, tissue, stage) cell.

    Every field a page needs to render the figure honestly travels with it,
    the same discipline `count_unit` enforces for every burden count:

    - `n_samples` -- the schema permits `n_samples` = 1, so "median" can be a
      single observation, and a page must be able to say so.
    - `n_genes` -- "top 4%" names no population without its denominator
      (D39(b)); this is the source matrix's whole gene universe, not the 92
      genes this atlas publishes.
    - `unit` -- a median with no unit is not a quantity.
    - `method` -- `LOOKUP_RULE` (D39(a)), so a consumer re-deriving this
      figure resolves a tied breakpoint the same way this module did.
    """

    q25_percentile: int | None
    median_percentile: int
    q75_percentile: int | None
    median_abundance: float
    unit: str
    n_samples: int
    n_genes: int
    method: str


def placement(
    median: float,
    q25: float | None,
    q75: float | None,
    breakpoints: Sequence[float],
    floor: float,
    unit: str,
    n_samples: int,
    n_genes: int,
) -> Placement | None:
    """The published figure for one gene in one cell, or `None` below the floor.

    D41: below the detection floor, a percentile is an artefact of the tied
    breakpoints at the bottom of the grid, not a measurement -- so `None` is
    returned rather than a number a page would render as though it meant
    something.

    **Judgement call: the comparison is `median < floor`, strictly.** A
    median exactly equal to the declared floor counts as detected, not below
    it. The floor is a source's own declared threshold for "present at all",
    so the boundary value itself is a detected value under that source's own
    convention; `<=` would exclude the one point the source chose to draw the
    line at, which is stricter than the source's own rule. Task 10's
    `specificity` gates tau's peak the same way (`peak_raw < floor`), so one
    gene cannot read "detected" under one figure and "below the floor" under
    the other from the same value. Pinned at the exact boundary by
    `test_the_floor_gate_is_strict_less_than`.
    """
    if median < floor:
        return None
    q25_percentile, median_percentile, q75_percentile = band(median, q25, q75, breakpoints)
    return Placement(
        q25_percentile=q25_percentile,
        median_percentile=median_percentile,
        q75_percentile=q75_percentile,
        median_abundance=median,
        unit=unit,
        n_samples=n_samples,
        n_genes=n_genes,
        method=LOOKUP_RULE,
    )
