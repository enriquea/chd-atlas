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

This module carries the percentile band (Task 9) and tau (Task 10). Task 11
adds phase assignment; Task 12 quantile shard emission.
"""

from __future__ import annotations

import bisect
import math
from collections.abc import Mapping, Sequence
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


# D39(a): the scale is part of the method, not an implementation detail --
# see `specificity`'s docstring for the measurement that makes this load-
# bearing (a 10-fold enriched gene reads tissue-specific on linear and broad
# on log2).
TAU_SCALE: Final = "log2(x+1)"

# D39(a): named beside every tau so a consumer can re-derive the number from
# `medians` without reading this module's source. The dataset's own detection
# floor is deliberately absent from this description -- it gates whether tau
# is published at all (see `specificity`), but no longer changes what any
# organ contributes once it is, so re-deriving tau needs nothing but
# `medians` itself. Same one-level guarantee `LOOKUP_RULE` documents for the
# percentile side (see the module docstring).
TAU_METHOD: Final = (
    "tau (Yanai et al. 2005): mean over organs of (1 - x_i/x_max), "
    "x = log2(median+1); a negative median is clamped to 0 before the "
    "transform, and every organ's raw median is used even below the "
    "dataset's detection floor"
)


class Specificity(TypedDict):
    """Tau with everything needed to read it, and nothing that reads as a verdict.

    No adjective and no band: tau 0.71 is "intermediate" under the common
    banding (<=0.5 broad, >=0.8-0.9 specific), and choosing a threshold is a
    classification the atlas would author, which D39(c) forbids. `highest_in`
    is what makes a page's wording safe -- tau measures concentration, not
    location, so a "heart-preferential" gloss keyed on tau alone can state
    the opposite of the truth for a gene concentrated elsewhere.

    - `scale` -- names the transform (D39(a)); the same input reads
      tissue-specific on linear and broad on log2, so this is not decoration.
    - `method` -- the formula and citation (D39(a)), `LOOKUP_RULE`'s
      counterpart for tau.
    - `tissues` -- the ordered organ panel tau was computed over. A count
      cannot say *which* organs contributed; a real panel can contain a
      correlated pair (cerebrum/cerebellum), so this matters.
    - `n_tissues` -- organs *sampled*, not organs *detected*: an organ below
      the detection floor is retained here, never dropped.
    - `highest_in` -- the argmax organ, or `None` when two or more organs tie
      at the peak, where "highest in X" would be an arbitrary choice.
    - `medians` -- tau's actual inputs (D39(b)), because the bundle's own
      slice of a gene's organs is truncated elsewhere (`omics.select_top`).
      Published exactly as measured, and tau's own arithmetic uses these
      same raw values -- no internal flooring -- so recomputing tau from
      `medians` alone reproduces `tau` exactly. See `specificity`'s
      docstring for why an earlier version of this function broke that.
    """

    tau: float
    scale: str
    method: str
    tissues: tuple[str, ...]
    n_tissues: int
    highest_in: str | None
    medians: dict[str, float]


def specificity(medians: Mapping[str, float], floor: float) -> Specificity | None:
    """Yanai tau for one gene's per-organ medians at one developmental stage.

    tau = sum(1 - x_i/x_max) / (n - 1), computed on `x = log2(median + 1)`.
    Measured: a 10-fold enriched gene (heart 100, six organs at 10) scores
    0.900 on linear RPKM -- "tissue-specific" -- and 0.480 on log2 --
    "broadly expressed". Yanai et al. 2005 and the Kryuchkova-Mostacci &
    Robinson-Rechavi benchmark both compute tau on log expression; computing
    on linear would call ordinary genes heart-preferential across the board.

    Returns `None`, never `0.0`, when tau is undefined:

    - fewer than two organs were sampled at this stage (the denominator is
      `n - 1`), or
    - the gene's peak across every sampled organ is below `floor` (the
      normaliser is `x_max`, and computing one from an unreliable peak would
      publish a ratio measured against noise).

    A genuinely ubiquitous gene -- every organ equal -- measures *exactly*
    0.000, so a guard that defaults an undefined tau to `0.0` (`n < 2`, or
    `except ZeroDivisionError`) publishes "expressed identically everywhere"
    for a gene measured in a single organ -- the opposite claim, and it lands
    on the earliest stages, where the organ panel is smallest.

    The peak's floor comparison is `peak_raw < floor`, strictly -- the same
    comparison `placement` uses (D41), so one gene cannot read "detected"
    under one figure and "below the floor" under the other from the same
    value; pinned at the boundary by `test_the_floor_gate_on_tau_is_strict_
    less_than`, matching Task 9's own boundary test.

    **The floor gates whether tau is published at all. It does not change
    what any organ contributes once tau is computed.** A non-peak organ
    below the floor is retained in `tissues`/`n_tissues` (so those keep
    meaning "organs sampled", not "organs detected" -- dropping is
    self-defeating, since the most heart-exclusive gene in the atlas,
    detected in heart alone, would fall to n=1 and lose tau entirely), *and
    its raw median is used in the sum exactly as measured* -- not floored to
    zero. An earlier version of this function floored it, reasoning that a
    sub-floor value is unreliable noise. Reverted on review, for two reasons:

    1. D39(b) requires tau to be re-derivable from the `medians` this same
       payload publishes, and `medians` is (and was always) published raw.
       Flooring internally while publishing raw inputs makes those two
       requirements contradict each other: a consumer who recomputes tau
       from the published `medians` gets a different number than the one
       published beside them (measured: 0.696 re-derived vs. 0.740
       published, same fixture). `test_tau_is_re_derivable_from_its_own_
       published_medians` makes this a checked invariant rather than a claim
       nobody watches.
    2. Flooring is not a neutral rounding choice: it measurably *raises*
       tau (0.740 vs. 0.696 on that same fixture -- `test_a_below_floor_
       nonzero_value_uses_its_raw_median_not_a_floored_zero`), pushing
       toward a stronger specificity claim than the data supports, in the
       same direction tau's own single-organ normaliser is already biased
       toward. Where two readings of an under-specified rule differ, this
       atlas takes the one that claims less.

    "Below the detection floor" means the source does not vouch for the
    value as a detection -- it does not mean the source measured zero, and
    substituting zero for it is this atlas authoring a number the source did
    not report (the act D12/D39 exist to forbid). Gating *publication* on
    the floor and *rewriting a value* because of it are different acts; this
    function does the first and, deliberately, not the second.

    The one rewrite tau's arithmetic does make is a negative median, clamped
    to zero before the log transform: a negative abundance is not a value,
    and `log2` of one is not a number. `test_a_negative_median_is_clamped_
    to_zero_before_log2` pins this, including that it does not raise.

    `medians` is published alongside because the bundle's own slice of a
    gene's organs is truncated elsewhere (`omics.select_top`), and D39(b)
    requires tau's actual inputs to be reachable from the same payload that
    carries tau.
    """
    tissues = tuple(sorted(medians))
    if len(tissues) < 2:
        return None

    peak_raw = max(medians[tissue] for tissue in tissues)
    if peak_raw < floor:
        return None

    def _log2p1(value: float) -> float:
        # The only rewrite applied to any organ's raw median: a negative
        # abundance is not a value, and log2 of one is not a number. A
        # value below the detection floor is NOT rewritten here -- see the
        # docstring's "gate publication, don't rewrite the value" argument.
        return math.log2(max(value, 0.0) + 1.0)

    logged = {tissue: _log2p1(value) for tissue, value in medians.items()}
    x_max = _log2p1(peak_raw)
    if x_max <= 0:
        return None

    tau = sum(1.0 - value / x_max for value in logged.values()) / (len(tissues) - 1)

    leaders = [tissue for tissue in tissues if medians[tissue] == peak_raw]
    return Specificity(
        tau=tau,
        scale=TAU_SCALE,
        method=TAU_METHOD,
        tissues=tissues,
        n_tissues=len(tissues),
        highest_in=leaders[0] if len(leaders) == 1 else None,
        medians=dict(medians),
    )
