# tests/unit/test_convert_cardoso_moreira.py
"""The Cardoso-Moreira converter's contracts, checked on synthetic records.

Deliberately NOT written against the committed mirror. A test that reads a
committed artifact cannot test the code that wrote it -- measured elsewhere in
this repository, where mutants relabelling a burden arm's `count_unit`, its
`origin` and its denominator all survived the full suite because every
assertion read the mirror rather than calling the builder. These call the
builders directly.

The two contracts that carry real risk here:

- **An absent measurement must not become a zero.** The Expression Atlas emits
  each cell as a `min,q25,median,q75,max` quintuple; a blank or `NA` cell means
  the gene was not quantified in that assay group, which is a different fact
  from a measured zero. Collapsing them would publish a floor-level abundance
  the source never reported, and D41 would then refuse a percentile for a gene
  that was never measured rather than one measured as absent.

- **The quantile estimator is this file's own.** It is written out rather than
  taken from a library so the rule that produced a published percentile is
  auditable from the converter alone, which means it needs pinning against
  worked literals -- a dependency's default could otherwise drift underneath it
  and move every breakpoint in the atlas without a line of this repo changing.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from convert_cardoso_moreira import (  # noqa: E402
    UNIT,
    AssayGroup,
    parse_cell,
    profile_row,
    quantile_linear,
    quantile_rows,
)


@pytest.mark.parametrize(
    ("cell", "expected"),
    [
        ("8,15,18,20,21", (15.0, 18.0, 20.0)),
        ("13,13,13,13,13", (13.0, 13.0, 13.0)),
        ("7", (7.0, 7.0, 7.0)),
    ],
)
def test_a_quintuple_yields_q25_median_q75(cell: str, expected: tuple[float, float, float]) -> None:
    """The Atlas hands over the quartiles; nothing here recomputes them.

    A single value is the degenerate quintuple the Atlas emits for a one-assay
    group, and is accepted as all five being equal.
    """
    assert parse_cell(cell) == expected


@pytest.mark.parametrize("cell", ["", "   ", "NA", "nan", "1,2,3"])
def test_an_unquantified_cell_is_none_and_never_zero(cell: str) -> None:
    """ "Not measured" and "measured as zero" are different facts.

    Returning 0.0 here would put a fabricated floor-level value into the
    quantile grid *and* give the gene a `median_abundance` the source never
    reported. A malformed cell (`1,2,3` -- neither one value nor five) is
    refused for the same reason rather than guessed at.
    """
    assert parse_cell(cell) is None


def test_quartiles_are_dropped_below_three_replicates_but_the_median_survives() -> None:
    """Quartiles of two points are not quartiles.

    The median of two is a mean of two with none of the robustness the word
    implies, so it still publishes -- with `n_samples` beside it saying how
    little it rests on -- while `q25`/`q75` go null rather than assert a spread
    two points cannot support.
    """
    three = AssayGroup("g1", "7 week post conception", "heart", 3)
    two = AssayGroup("g2", "4 week post conception", "heart", 2)

    rich = profile_row("E-MTAB-6814", "HGNC:11604", three, 15.0, 18.0, 20.0)
    assert (rich["q25"], rich["q75"], rich["n_samples"]) == (15.0, 20.0, 3)

    thin = profile_row("E-MTAB-6814", "HGNC:11604", two, 15.0, 18.0, 20.0)
    assert (thin["q25"], thin["q75"]) == (None, None)
    assert thin["median_abundance"] == 18.0
    assert thin["n_samples"] == 2
    assert thin["unit"] == UNIT == "tpm"


@pytest.mark.parametrize(
    ("percentile", "expected"),
    [(0, 0.0), (25, 24.75), (50, 49.5), (75, 74.25), (100, 99.0)],
)
def test_the_quantile_estimator_is_r_type_7(percentile: int, expected: float) -> None:
    """Pinned against worked literals, not against whatever a library returns.

    With n = 100 and h = (n - 1) * p, p25 falls at position 24.75, i.e. three
    quarters of the way from 24 to 25. These are the values numpy's default
    `method="linear"` and R's `type=7` produce, and they are asserted here so
    the atlas's published percentiles cannot move because a dependency changed
    its default.
    """
    assert quantile_linear(list(range(100)), percentile) == pytest.approx(expected)


def test_a_grid_is_101_rows_non_decreasing_and_states_its_denominator() -> None:
    """PRF003 requires non-decreasing, and ties at the floor are the normal case.

    Roughly half a real bulk library sits at exactly zero, so several dozen
    consecutive breakpoints are legitimately identical -- a strictly increasing
    grid would be the anomaly. `n_genes` is the whole gene universe the rank was
    taken over, without which "top 4%" names no population.
    """
    group = AssayGroup("g1", "7 week post conception", "heart", 3)
    medians = [0.0] * 46 + [float(i) for i in range(54)]
    rows = list(quantile_rows("E-MTAB-6814", group, medians))

    assert len(rows) == 101
    assert [row["percentile"] for row in rows] == list(range(101))
    values = [row["value"] for row in rows]
    # Not `strict=True`: the two sequences differ in length by one by
    # design -- this walks adjacent pairs, it does not zip equals.
    assert all(a <= b for a, b in zip(values, values[1:]))
    assert values[0] == 0.0
    assert {row["n_genes"] for row in rows} == {100}
    assert {row["unit"] for row in rows} == {"tpm"}
