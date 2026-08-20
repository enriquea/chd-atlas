# tests/unit/test_build_charts.py
"""Primitives for inline SVG. Nothing here knows what a gene is."""

from __future__ import annotations

import re
from typing import Final

import pytest

from chd_atlas.build.charts import LogScale, coordinate, marker, polyline, svg_figure
from chd_atlas.build.render import STYLESHEET


def test_a_log_scale_places_its_endpoints_on_the_axis_ends() -> None:
    scale = LogScale(low=1.0, high=1000.0, left=10.0, width=300.0)
    assert scale.x(1.0) == 10.0
    assert scale.x(1000.0) == 310.0
    assert scale.x(31.62277) == pytest.approx(160.0, abs=0.1)


def test_a_value_at_or_below_zero_is_refused_rather_than_placed() -> None:
    """`log10(0)` is `-inf` and `log10(-1)` raises.

    A silently clamped zero would place a below-floor measurement at the axis
    origin, which reads as "lowest measured" rather than "not measured".
    Callers must decide what an unplaceable value renders as; this refuses to
    decide for them. `raise`, never `assert` -- `-O` strips `assert`.
    """
    scale = LogScale(low=1.0, high=100.0, left=0.0, width=100.0)
    with pytest.raises(ValueError, match="positive"):
        scale.x(0.0)
    with pytest.raises(ValueError, match="positive"):
        scale.x(-1.0)


@pytest.mark.parametrize(
    ("low", "high", "width", "message"),
    [
        (0.0, 100.0, 100.0, "positive low"),
        (-1.0, 100.0, 100.0, "positive low"),
        (10.0, 10.0, 100.0, "high above low"),
        (100.0, 1.0, 100.0, "high above low"),
        (1.0, 100.0, 0.0, "positive width"),
        (1.0, 100.0, -5.0, "positive width"),
    ],
)
def test_an_axis_that_cannot_place_anything_is_refused_at_construction(
    low: float, high: float, width: float, message: str
) -> None:
    """Three degenerate axes, and only one of them announces itself.

    Measured 2026-08-20 against the arithmetic in `LogScale.x` with the
    constructor guard removed: `low <= 0` raises `math domain error` and
    `high == low` raises `float division by zero`, both at draw time and
    neither naming the axis -- but `width <= 0` **raises nothing and returns
    `left` for every value**, stacking every measurement on one pixel so the
    chart reads as "these numbers are all equal". `high < low` is the same
    shape: no error, a negative span, and a rising series drawn falling.

    Refused at construction rather than on the first point, so a caller that
    built a scale and drew nothing still fails.
    """
    with pytest.raises(ValueError, match=message):
        LogScale(low=low, high=high, left=0.0, width=width)


def test_coordinates_are_fixed_precision_so_two_builds_agree() -> None:
    """Formatted, not `repr`'d. A fixed format makes the byte-for-byte
    reproducibility guarantee hold for chart output too."""
    assert coordinate(1.0 / 3.0) == "0.3"
    assert coordinate(2.0) == "2.0"


def test_negative_zero_formats_as_zero() -> None:
    """`-0.0` and `0.0` are the same point and must be the same bytes.

    A subtraction landing on `-0.0` would otherwise emit `-0.0` on one build
    and `0.0` on another depending only on operand order.

    The second case is why this is normalised *after* formatting rather than by
    adding zero first. Measured 2026-08-20: `-0.0 + 0.0` is `+0.0`, so adding
    zero does fix the first case -- but `-1e-9 + 0.0` is still `-1e-9`, and
    `f"{-1e-9:.1f}"` is `"-0.0"`. Any coordinate a hair below the axis origin
    formats to a signed zero, and a rounding difference of one part in a
    billion must not change a published byte.
    """
    assert coordinate(-0.0) == "0.0"
    assert coordinate(-1e-9) == "0.0"


def test_a_figure_is_labelled_for_a_screen_reader() -> None:
    """`role="img"` plus a `<title>`. A chart a screen reader cannot convey is
    a chart that loses evidence for the reader who most needs the text."""
    out = svg_figure(width=100, height=50, title="TBX5 in heart", body="")
    assert 'role="img"' in out
    assert "<title>TBX5 in heart</title>" in out


def test_a_figure_escapes_its_title() -> None:
    """Titles carry mirror-supplied tissue and stage tokens."""
    out = svg_figure(width=10, height=10, title='x<script>"&', body="")
    assert "<script>" not in out
    assert "&lt;script&gt;" in out


def test_no_primitive_emits_a_colour() -> None:
    """Every colour comes from `render.STYLESHEET` via a class.

    The page defines a full dark palette under `prefers-color-scheme: dark`.
    A baked-in hex renders a light-mode chart on a dark page.
    """
    rendered = (
        svg_figure(width=10, height=10, title="t", body="")
        + polyline([(0.0, 0.0), (1.0, 1.0)], css_class="chart-line")
        + marker(1.0, 2.0, css_class="chart-point")
    )
    assert "#" not in rendered
    assert "rgb" not in rendered


def test_a_polyline_needs_at_least_two_points() -> None:
    """One point is not a line.

    A later guard depends on this: one published gene is placed at 1 of 19
    stages, and a polyline through it would draw a trend the data does not
    contain.
    """
    with pytest.raises(ValueError, match="two"):
        polyline([(0.0, 0.0)], css_class="chart-line")


def test_the_chart_stylesheet_names_no_colour_of_its_own() -> None:
    """Every chart rule resolves its colour through a custom property.

    `test_no_primitive_emits_a_colour` above guards `charts.py`, which is only
    half the surface: the geometry carries class names and the *stylesheet*
    carries the colours. Measured while adding this -- swapping
    `stroke: var(--link)` for `stroke: #0b5fa5` in the chart block passed the
    entire suite, because nothing read the stylesheet at all. Dark mode would
    have broken for every reader with that preference set, silently, under a
    green build.

    Scoped to chart rules rather than the whole sheet: the palette blocks at
    the top of `STYLESHEET` define the hex literals every `var()` here resolves
    to, and must keep them.
    """
    chart_rules = [
        rule for rule in re.findall(r"(\.(?:chart|spark|forest)[^{}]*\{[^}]*\})", STYLESHEET)
    ]
    assert chart_rules, "no chart rules found; this guard would pass vacuously"
    offenders = [rule for rule in chart_rules if re.search(r"#[0-9a-fA-F]{3,8}\b", rule)]
    assert offenders == [], f"chart CSS hardcodes a colour instead of a token: {offenders}"


# Which declarations carry an encoding a reader can see. `stroke-width` and
# `stroke-dasharray` are deliberately absent: they are shape, and a chart may
# distinguish two marks by weight alone without either resolving a colour.
_PAINT: Final = ("fill", "stroke")
_COLOURED: Final = (*_PAINT, "color", "background")

_RULE = re.compile(r"(\.(?:chart|spark|forest)[^{}]*)\{([^}]*)\}")
_TOKEN = re.compile(r"var\((--[a-z-]+)\)")


def _chart_rules() -> dict[str, dict[str, str]]:
    """Every chart rule in `STYLESHEET`, as `{selector: {property: value}}`."""
    rules: dict[str, dict[str, str]] = {}
    for selector, body in _RULE.findall(STYLESHEET):
        declarations = {}
        for declaration in body.split(";"):
            name, sep, value = declaration.partition(":")
            if sep:
                declarations[name.strip()] = value.strip()
        rules[selector.strip()] = declarations
    return rules


def test_every_chart_mark_resolves_a_paint_and_no_two_meanings_share_one() -> None:
    """The stylesheet *is* the encoding, and only its hex-freedom was guarded.

    `test_the_chart_stylesheet_names_no_colour_of_its_own` above greps for a
    literal `#`, so every rule below it can be gutted or collapsed into its
    neighbour and the suite stays green. Two measurements from 2026-08-20, on
    a build of the committed corpus:

    * `.chart-point-open { fill: var(--bg) }` changed to `var(--link)` erases
      the "did not survive its own study's correction" encoding for **all 604
      hollow circles on the site**, leaving them identical to the 250 filled
      ones. Its sibling `.chart-arrow-open` *is* pinned against a literal by
      `test_an_arrow_carries_the_same_correction_fill_as_a_circle` -- guarded
      on one of a pair, in the release that cites the other (CLAUDE.md
      section 4.31).
    * `.chart-label` losing its `fill` gives black text on a dark page for all
      **2,289 labels**, which is precisely the defect the comment beside that
      rule names: SVG `<text>` takes no `color`, only `fill`.

    Four rules, structural rather than per-declaration, so a chart class added
    later is covered without anyone remembering this test:

    1. every paint resolves through a custom property or is `none` -- the
       positive form of the hex grep, which a `fill: black` would satisfy;
    2. every token named is one the palette actually declares, since a
       misspelt `var(--lnk)` renders as nothing at all and no hex check sees
       it;
    3. every mark class declares a paint, which is what `.chart-label` lost;
    4. the two hollow/solid pairs differ from their solid partners, which is
       what makes the encoding impossible to collapse.
    """
    rules = _chart_rules()
    assert rules, "no chart rules found; this guard would pass vacuously"
    declared = set(re.findall(r"(--[a-z-]+):", STYLESHEET))

    for selector, declarations in rules.items():
        for name, value in declarations.items():
            if name in _COLOURED and value != "none":
                assert "var(" in value, f"{selector} sets {name} to a literal colour: {value}"
            for token in _TOKEN.findall(value):
                assert token in declared, f"{selector} resolves {token}, which no palette declares"

    # `.chart` itself is layout only -- width, height, overflow -- and paints
    # nothing; every other chart class names a mark and must colour it.
    marks = {
        selector: declarations
        for selector, declarations in rules.items()
        if selector.startswith(".chart-")
    }
    assert len(marks) > 8, f"only {len(marks)} mark classes found; the pattern has drifted"
    for selector, declarations in marks.items():
        assert any(name in declarations for name in _PAINT), (
            f"{selector} declares no fill and no stroke, so it inherits the browser's black"
        )

    # A hollow mark is hollow: filled with the page's own background, never
    # with the colour its solid partner uses. Both pairs, because a guard on
    # one of a pair is evidence about one.
    for solid, hollow in (
        (".chart-point", ".chart-point-open"),
        (".chart-arrow", ".chart-arrow-open"),
    ):
        assert rules[hollow]["fill"] == "var(--bg)", f"{hollow} is not drawn hollow"
        assert rules[solid]["fill"] != rules[hollow]["fill"], (
            f"{solid} and {hollow} are the same paint; the correction encoding is gone"
        )
