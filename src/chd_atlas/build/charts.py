# src/chd_atlas/build/charts.py
"""Inline-SVG primitives. This module does not know what a gene is.

`render.py` turns text into a page and knows nothing about genes; this turns
numbers into geometry on the same terms. `pages.py` is the only module that
joins either to data, which is what keeps 2,000 lines of copy and a coordinate
system from growing into each other.

**No primitive here emits a colour.** Every stroke and fill comes from a class
defined once in `render.STYLESHEET`, against the custom properties that already
carry the page's `prefers-color-scheme: dark` palette. A hardcoded hex would
render a light chart on a dark page for every reader who has that preference
set, and nothing in the suite would fail.

Deterministic for the same reason every other builder here is: nothing reads a
clock, and every coordinate is formatted rather than `repr`'d, so two builds of
one commit are byte-identical.

**There is deliberately no primitive, and no `chart-diamond` class, for a
forest plot's summary diamond.** A diamond draws a *pooled* estimate, and D33
forbids one: `pages._POOLING_NOTICE` states the reason in full and
`build/burden.py::shared_cohorts` is what acts on it. Defining the class would
be the first half of that defect, and the forest plot's tests assert it never
appears -- so a caller wanting one has a data decision to make, not a missing
shape to add. Recorded here rather than in `render.STYLESHEET`, because that
string is inlined verbatim into every page: a comment there is published bytes,
and the token a page-wide check greps for must not be among them.
"""

from __future__ import annotations

import html
import math
from collections.abc import Sequence
from dataclasses import dataclass


def coordinate(value: float) -> str:
    """One axis coordinate, fixed to a tenth of a pixel.

    A signed zero is normalised **after** formatting, not by adding zero
    before it, and the difference is load-bearing. Measured 2026-08-20:
    `-0.0 + 0.0` is `+0.0`, so adding zero handles an exact negative zero --
    but `-1e-9 + 0.0` is still `-1e-9`, and `f"{-1e-9:.1f}"` is `"-0.0"`. Any
    coordinate a hair below the axis origin formats to a signed zero, so the
    pre-format trick fixes the case nobody hits and misses the case a
    subtraction actually produces. `.1f` emits one decimal place, so `"-0.0"`
    is the only string form a negative zero can take here.

    It matters because `-0.0` and `0.0` are the same point: a subtraction
    landing on one rather than the other would change published bytes on the
    strength of operand order alone, and two builds of one commit must be
    byte-identical.
    """
    text = f"{value:.1f}"
    return "0.0" if text == "-0.0" else text


@dataclass(frozen=True)
class LogScale:
    """Maps a positive value onto a pixel axis, base 10.

    Log rather than linear because both axes this serves span orders of
    magnitude: one gene's heart median runs 5 to 275 tpm across development,
    and published odds ratios run 0.0 to unbounded. On a linear axis the
    post-natal half of a trajectory collapses onto the floor.
    """

    low: float
    high: float
    left: float
    width: float

    def __post_init__(self) -> None:
        """Refuses an axis that cannot place anything, at construction.

        Each of the three costs something different, and one of them costs
        silence. Measured 2026-08-20 against the arithmetic in `x` below:

        - `low <= 0` raises `ValueError: math domain error` from `log10`, but
          only once a caller draws its first point, and the message names
          neither the axis nor the argument.
        - `high == low` gives a zero span and raises `ZeroDivisionError:
          float division by zero`, again at draw time and again opaquely.
          `high < low` raises nothing at all: the span goes negative and the
          axis silently inverts, so a rising trajectory is drawn falling.
        - `width <= 0` **raises nothing and returns `left` for every value.**
          Every measurement stacks on one pixel, which reads as "all these
          numbers are the same" rather than as a broken axis. That is the one
          this guard exists for; the other two would surface eventually.

        `raise`, never `assert` -- `-O` strips `assert`.
        """
        if self.low <= 0:
            raise ValueError(f"a log axis needs a positive low bound, got {self.low!r}")
        if self.high <= self.low:
            raise ValueError(
                f"a log axis needs high above low, got low={self.low!r} high={self.high!r}"
            )
        if self.width <= 0:
            raise ValueError(
                f"a log axis needs a positive width, got {self.width!r}; a zero-width axis "
                f"places every value at the same pixel rather than failing"
            )

    def x(self, value: float) -> float:
        """Where `value` sits. Refuses a non-positive value rather than clamping.

        A clamped zero would land on the axis origin and read as "the lowest
        measurement here", which is a different claim from "below the detection
        floor" and from "not sampled". The caller owns that decision because
        only the caller knows which of the three it has.
        """
        if value <= 0:
            raise ValueError(f"a log axis needs a positive value, got {value!r}")
        span = math.log10(self.high) - math.log10(self.low)
        offset = (math.log10(value) - math.log10(self.low)) / span
        return self.left + self.width * offset


def svg_figure(*, width: int, height: int, title: str, body: str) -> str:
    """One `<svg>`, labelled for assistive technology.

    `role="img"` with a `<title>` rather than a bare `<svg>`: a screen reader
    otherwise announces the individual shapes or nothing at all. The `<details>`
    table each chart sits above is the full text alternative; this is the
    one-sentence summary that says whether opening it is worth it.
    """
    return (
        f'<svg class="chart" viewBox="0 0 {width} {height}" role="img" '
        f'preserveAspectRatio="xMinYMin meet">'
        f"<title>{html.escape(title)}</title>{body}</svg>"
    )


def polyline(points: Sequence[tuple[float, float]], *, css_class: str) -> str:
    """A connected run of measurements.

    Refuses fewer than two points. One point is not a line, and a renderer that
    quietly emitted a zero-length polyline for a gene measured at a single
    stage would draw a trend from one measurement. Live case, not
    hypothetical: measured 2026-08-20 on the committed corpus, FOXH1 and
    TFAP2B are each placed in heart at 1 of the 19 stages this dataset samples
    heart at, and 18 of the 581 (gene, organ) series with anything placed at
    all have exactly one placed point.
    """
    if len(points) < 2:
        raise ValueError(f"a polyline needs at least two points, got {len(points)}")
    drawn = " ".join(f"{coordinate(x)},{coordinate(y)}" for x, y in points)
    return f'<polyline class="{html.escape(css_class)}" points="{drawn}"/>'


def marker(x: float, y: float, *, css_class: str) -> str:
    """One measured point."""
    return (
        f'<circle class="{html.escape(css_class)}" '
        f'cx="{coordinate(x)}" cy="{coordinate(y)}" r="2.5"/>'
    )
