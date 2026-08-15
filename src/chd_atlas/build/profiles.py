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

This module carries the percentile band (Task 9), tau (Task 10), phase
assignment (Task 11) and quantile shard emission (Task 12). Task 12 is the
other half of D39(b)'s bargain: Task 9's percentile is re-derivable only if
the breakpoints it was read against are themselves fetchable, and
`build_omics` never emits this table -- it skips every schema absent from its
own `_GENE_COLUMN`, and a quantile grid has no gene column at all.
"""

from __future__ import annotations

import bisect
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Final, TypedDict

from chd_atlas.build.emit import Emitter
from chd_atlas.build.paths import slug
from chd_atlas.models.dataset import Stage
from chd_atlas.models.phases import CardiacPhaseFile
from chd_atlas.tables import TABLE_SCHEMAS, mirror_paths, read_table

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


class PhaseOutcome(StrEnum):
    """Why `assign_phase` returned the `phase_id` (or lack of one) that it did.

    Four outcomes, never a bare `None`: a page that only checked
    `phase_id is None` would render one sentence for a post-natal stage, a
    stage past the curated window, and a stage the dataset's own record never
    declared -- three different facts, and `UNDECLARED` is not a biological
    statement at all, it is a curation gap PRF005 also reports from the other
    direction. `Resolution` in `chd_atlas/genes.py` is the precedent: a
    caller branches on the named outcome, never on which fields happen to be
    `None` or on matching the display text in `reason`.
    """

    MATCHED = "matched"
    """`wpc` falls inside a declared cardiac phase; `phase_id` names it."""

    OUTSIDE_WINDOW = "outside_window"
    """A real `wpc` exists, but no curated phase covers it -- before the
    first, after the last, or in an interior gap between two (PRF006 flags
    the gap shape as a probable transcription slip; this is the runtime fact
    a page renders regardless of *why* the gap exists). The committed
    `curation/cardiac_phases.yaml` declares zero phases today -- the
    boundaries have not yet been transcribed from a verified source -- so on
    the real corpus every stage with a real `wpc` resolves here. That is
    correct, not a bug and not a blank: Task 4's placeholder file exists
    precisely so this reads as "outside the curated window" rather than as
    silence."""

    POST_NATAL = "post_natal"
    """The stage's own `wpc` is null. Every `CardiacPhase` is prenatal by
    construction (`start_wpc`/`end_wpc` are both `gt=0`), so this is a
    different fact than OUTSIDE_WINDOW: there is no developmental age to
    place at all, and a wider curated window could never change the
    answer."""

    UNDECLARED = "undeclared"
    """`token` is not in the dataset's own `stages`, so there is no `Stage`
    to read a `wpc` from in the first place. PRF005 reports the reverse
    direction (a `profiles` row naming a stage the dataset record does not
    declare); this is what a caller sees asking about a token with nothing
    behind it, and `assign_phase` must not guess a biological placement for
    a stage it was never told about."""


@dataclass(frozen=True)
class PhaseAssignment:
    """Where one of a dataset's own stage tokens falls in the curated phase vocabulary.

    `phase_id` is set only for `PhaseOutcome.MATCHED`. `reason` is a short,
    renderable phrase explaining an absent `phase_id`, and is `None` exactly
    when `phase_id` is not -- a page renders one field or the other, never
    both, and never neither.
    """

    token: str
    outcome: PhaseOutcome
    phase_id: str | None = None
    reason: str | None = None


def assign_phase(
    token: str, stages: Sequence[Stage], phases: CardiacPhaseFile | None
) -> PhaseAssignment:
    """Place one of a dataset's own stage tokens in the curated phase vocabulary.

    Checked in this order, each a stronger claim than the last is absent:

    1. `token` is not in `stages` at all -- `UNDECLARED`. Nothing about
       development can be said for a token with no `Stage` behind it, so this
       is decided before either of the wpc-based questions below are even
       asked.
    2. The matching `Stage.wpc` is `None` -- `POST_NATAL`. Every `CardiacPhase`
       is prenatal by construction, so a null wpc can never match one
       regardless of what `phases` declares.
    3. `phases` is absent entirely, or declares zero phases -- `OUTSIDE_WINDOW`,
       reason "outside the curated window". Checked before indexing into
       `phases.phases` for the boundary comparison in step 4, which would
       raise `IndexError` on an empty sequence otherwise -- and this is not a
       hypothetical: it is the committed corpus's own state today.
    4. Otherwise `phases.phase_for(wpc)` is tried. A hit is `MATCHED`. A miss
       is `OUTSIDE_WINDOW` again, with a reason naming *where*: "before the
       curated window" (earlier than every phase's `start_wpc`), "after the
       curated window" (at or past the last phase's `end_wpc` -- the common
       case, since a real developmental series runs well past any
       morphogenetic window, not an edge case), or the same generic phrase as
       step 3 for the rarer case of an interior gap between two declared
       phases, which PRF006 already names as a probable transcription slip.
    """
    stage = next((candidate for candidate in stages if candidate.token == token), None)
    if stage is None:
        return PhaseAssignment(
            token=token,
            outcome=PhaseOutcome.UNDECLARED,
            reason="stage not declared by this dataset",
        )
    if stage.wpc is None:
        return PhaseAssignment(token=token, outcome=PhaseOutcome.POST_NATAL, reason="post-natal")
    if phases is None or not phases.phases:
        return PhaseAssignment(
            token=token,
            outcome=PhaseOutcome.OUTSIDE_WINDOW,
            reason="outside the curated window",
        )
    found = phases.phase_for(stage.wpc)
    if found is not None:
        return PhaseAssignment(token=token, outcome=PhaseOutcome.MATCHED, phase_id=found.id)

    starts = [phase.start_wpc for phase in phases.phases]
    ends = [phase.end_wpc for phase in phases.phases]
    if stage.wpc < min(starts):
        reason = "before the curated window"
    elif stage.wpc >= max(ends):
        reason = "after the curated window"
    else:
        reason = "outside the curated window"
    return PhaseAssignment(token=token, outcome=PhaseOutcome.OUTSIDE_WINDOW, reason=reason)


def build_profile_quantiles(root: Path, emitter: Emitter) -> dict[str, str]:
    """Publish every quantile grid, so a percentile this atlas derives can be checked.

    D39(b) is unmet without this. A gene page states a percentile whose input
    is the whole source transcriptome (tens of thousands of genes), and D32
    forbids re-hosting that matrix -- so `mirrors/profile_quantiles/<accession>.tsv`
    exists solely to publish the 101 breakpoints a consumer needs to place any
    abundance and reproduce the number themselves. `build_omics` never emits
    it: that function skips every schema absent from its own `_GENE_COLUMN`
    (`omics.py`'s `if schema_name not in _GENE_COLUMN: continue`), and a
    quantile grid has no gene column -- there is no gene to attribute a
    breakpoint row to. Without this function the table is mirrored,
    schema-validated, sort-checked and sha256'd, and reaches no published
    byte: the one file whose entire purpose is auditability would be the one
    thing a consumer could not fetch, and nothing would catch it -- the
    validators check the mirror against itself, `percentile_of`'s arithmetic
    is tested against unpublished input either way, and a `diff -rq` between
    two builds shows a file that appears, never one that should have.

    Emits `omics/profile_quantiles/<accession>.json` as
    `{"table": "profile_quantiles", "rows": [...]}` -- the same envelope
    `build_omics` writes for its own shards, so a consumer reads one shape
    regardless of which table it fetched. The accession is the shard's own
    filename stem, exactly as `build_omics` reads a dataset accession from a
    `profiles`/`expression`/`proteomics`/`phospho` shard's name; it is put
    through `paths.slug` for the same reason every shard stem is -- the stem
    becomes a URL, and a space or a colon in one would need escaping before it
    could be fetched.

    Rows are sorted by the table's own canonical `sort_key`
    (`dataset, tissue, stage, percentile`) before being written, never trusted
    to already be in that order: `read_table` preserves a TSV's file order
    verbatim (`pl.read_csv`, no grouping or sort of its own), so an unsorted
    row order on disk would otherwise publish unsorted, silently, on a build
    that never re-validates a mirror it has already read. Sorting on the full
    key rather than on `percentile` alone matters the moment one shard holds
    more than one `(tissue, stage)` grid, which is the ordinary shape for a
    real dataset (many organs, many stages) rather than the single-grid shape
    of the test fixture -- sorting on `percentile` alone would interleave
    every grid in the shard by breakpoint value instead of keeping each one
    contiguous and internally ordered.

    Shards are emitted in `mirror_paths`' own order, which is sorted by
    filename -- so two builds of one commit write these files in the same
    order, though (unlike row order within one file) that has no bearing on
    any file's *content*.

    Returns `{dataset accession: relative shard path}`. Not yet read by any
    caller: Task 13 is what threads a gene's `expression_profile` bundle key
    to the grid its percentile came from, the way `ModalitySummary.shards`
    lets a gene bundle reach an omics shard today. Returning the exact path
    this function wrote -- rather than leaving Task 13 to reconstruct it via
    a second call to `slug` -- is what keeps that link from becoming a second
    computation that could drift from the first, the same discipline
    `omics.py`'s own module docstring states for `count`.
    """
    shards: dict[str, str] = {}
    for path, schema_name in mirror_paths(root):
        if schema_name != "profile_quantiles":
            continue
        frame, _ = read_table(path, TABLE_SCHEMAS[schema_name])
        if frame is None:
            # Unreadable is `validate_table`'s to report against this same
            # path; failing the build here would say it a second time with
            # less context. Unreachable behind `build_site`'s validation
            # gate in the same way `build_omics`'s equivalent branch is --
            # kept as a guard on that bypassed gate rather than trusted to
            # stay true as the two modules evolve apart.
            continue
        rows = frame.sort(list(TABLE_SCHEMAS[schema_name].sort_key)).to_dicts()
        accession = path.stem
        relative = f"omics/profile_quantiles/{slug(accession)}.json"
        emitter.write_json(relative, {"table": schema_name, "rows": rows})
        shards[accession] = relative
    return shards
