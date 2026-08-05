# tests/unit/test_convert_burden.py
"""The converter's two published contracts, checked against the real mirror.

There is no synthetic workbook fixture here, deliberately. The xlsx reader's one
realistic failure -- assuming the transitional OOXML namespace when this
workbook uses the strict one -- yields *zero* rows, and `_report` prints
"0 row(s) for 0 of the 154 registered genes", so it cannot pass unnoticed. What
can pass unnoticed is a wrong number in a right-looking row, and the committed
mirror is where that would live.
"""

from __future__ import annotations

import sys
from pathlib import Path

import polars as pl
import pytest

from chd_atlas.tables import BURDEN, TABLE_SCHEMAS, read_table

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from convert_burden import _number  # noqa: E402


@pytest.fixture(scope="module")
def burden() -> pl.DataFrame:
    frame, issues = read_table(REPO_ROOT / "mirrors" / "burden.tsv", BURDEN)
    assert frame is not None, issues
    return frame


@pytest.mark.parametrize(
    ("published", "expected"),
    [
        # What the workbook actually stores, and what it means. Both strings
        # denote the identical IEEE-754 double, so this is a rendering choice
        # and not a rounding one -- `float(expected) == float(published)`.
        ("0.61099999999999999", "0.611"),
        ("2.4500000000000002", "2.45"),
        ("0.29799999999999999", "0.298"),
        ("1E-3", "0.001"),
        # Exponents survive: a p-value of 1.21e-10 must not become 0.0.
        ("1.21e-10", "1.21e-10"),
        ("0", "0.0"),
    ],
)
def test_a_published_float_is_written_at_its_shortest_exact_spelling(
    published: str, expected: str
) -> None:
    """17 significant digits on a page is false precision the study never claimed."""
    assert _number(published) == expected
    assert float(_number(published)) == float(published)


def test_the_stratum_each_row_names_is_the_one_its_odds_ratio_was_computed_from(
    burden: pl.DataFrame,
) -> None:
    """The converter's single most dangerous claim, re-derived from the mirror.

    Supplementary Data 3 gives every sheet all four carrier counts and differs
    only in which comparison its `fet.*` columns describe. Nothing in the file
    says which -- the mapping was established by recomputing the odds ratio as
    (a/(N-a)) / (c/(M-c)) and seeing which pairing agreed. Getting it wrong
    would publish the syndromic-only odds ratio against the whole-cohort
    numerator: a plausible number, a wrong claim, and invisible in any per-column
    check.

    Checked here on the committed mirror rather than on the workbook, so it
    guards the artifact consumers actually read. Rows with no case carrier or no
    control carrier are skipped: their odds ratio is 0 or unbounded, which
    carries no information about which pairing produced it. 860 of the 1,192
    rows remain, and every one agrees to within 2%, the precision the paper
    published at.
    """
    checked = 0
    for row in burden.iter_rows(named=True):
        a, n = row["n_case_carriers"], row["n_cases"]
        c, m = row["n_control_carriers"], row["n_controls"]
        if row["effect"] is None or a == 0 or c == 0:
            continue
        expected = (a / (n - a)) / (c / (m - c))
        assert abs(expected - row["effect"]) <= max(0.02, 0.02 * expected), row
        checked += 1

    assert checked == 860


def test_an_unbounded_odds_ratio_publishes_no_number_at_all(burden: pl.DataFrame) -> None:
    """`encode_json` sets `allow_nan=False`, so `Infinity` cannot reach a page.

    It is accepted by Python's `json.loads` and rejected by `JSON.parse`, so a
    row carrying one would break every page while the manifest checksum still
    verified. These 34 rows also hold the study's strongest results, so dropping
    them was not an option either.

    The biconditional is what matters: `effect_bound` is set on exactly the rows
    with no effect, and both directions are asserted. A converter that flagged
    the bound but wrote a placeholder number, or wrote null without the flag,
    would satisfy one half and publish either a fabricated ceiling or a blank
    cell a reader would read as "not tested".
    """
    unbounded = burden.filter(pl.col("effect_bound").is_not_null())
    assert unbounded.height == 34
    assert unbounded["effect_bound"].unique().to_list() == ["unbounded_above"]
    assert unbounded["effect"].null_count() == 34
    assert unbounded["ci_high"].null_count() == 34
    # ...and no other row has a null effect.
    assert burden["effect"].null_count() == 34
    # The lower bound survives: it is finite, and it is the whole result. TAB2's
    # syndromic row says the true odds ratio is at least 28.1.
    assert unbounded["ci_low"].null_count() == 0


def test_the_committed_mirror_covers_every_gene_the_site_publishes(
    burden: pl.DataFrame,
) -> None:
    """Coverage, and the shape of it.

    Table 1 of the paper lists 21 genes, which is what made the first reading of
    this supplement conclude 4 of the 23 published genes were covered. The full
    table covers all 23.

    **The matrix is sparse, and this docstring said otherwise until measured.**
    It claimed nine rows per gene -- three strata by three consequence classes.
    Measured: of the 145 genes, 103 have 9 rows, 13 have 8, 24 have 6, 1 has 5
    and 4 have 3; of the 23 published genes, 11 have 9, 8 have 8 and 4 have 6.

    What a missing cell means is pinned by
    `test_a_missing_cell_means_no_carrier_rather_than_no_test` below, and the
    renderer depends on it.

    The synonymous count is asserted because it is the one a well-meaning filter
    would remove: synonymous variants are not a result, they are the study's own
    negative control, and a gene whose synonymous row is significant is a gene
    whose comparison is poorly calibrated. That they are the *most* numerous
    class -- 435 against 345 loss-of-function -- is itself the point: a gene is
    likelier to carry a synonymous variant than a truncating one, so the control
    is available almost everywhere the result is.
    """
    genes, _ = read_table(REPO_ROOT / "mirrors" / "genes.tsv", TABLE_SCHEMAS["genes"])
    assert genes is not None

    assert burden.height == 1192
    assert burden["gene"].n_unique() == 145
    assert set(burden["gene"].to_list()) <= set(genes["hgnc_id"].to_list())
    assert burden["consequence_class"].value_counts(sort=True).to_dicts() == [
        {"consequence_class": "synonymous", "count": 435},
        {"consequence_class": "missense_damaging", "count": 412},
        {"consequence_class": "lof", "count": 345},
    ]
    assert sorted(burden["cohort_stratum"].unique().to_list()) == [
        "all",
        "nonsyndromic",
        "syndromic",
    ]


def test_a_missing_cell_means_no_carrier_rather_than_no_test(burden: pl.DataFrame) -> None:
    """The semantic the gene page's empty cells depend on.

    The burden matrix is sparse: 42 of the 145 genes are missing at least one
    (stratum, consequence) cell. A reader meeting a blank will read it as
    "not tested" unless the page says otherwise, and "not tested" is a different
    and much weaker claim than "nobody carried one".

    Zero of the 1,192 rows have no carrier in the cases *and* none in the
    controls. A 2x2 table of all zeros supports no test, so the study emitted no
    row for it -- which makes an absent cell mean exactly "no qualifying variant
    was seen in either group", and the renderer may say so.

    The two one-sided cases are real rows and stay: 298 rows have no case
    carrier (evidence against enrichment) and 34 have no control carrier (the
    unbounded odds ratios, the strongest signal in the study).
    """
    assert burden.filter(
        (pl.col("n_case_carriers") == 0) & (pl.col("n_control_carriers") == 0)
    ).is_empty()
    assert burden.filter(pl.col("n_case_carriers") == 0).height == 298
    assert burden.filter(pl.col("n_control_carriers") == 0).height == 34


def test_the_denominators_are_the_studys_post_qc_set_not_its_abstract(
    burden: pl.DataFrame,
) -> None:
    """3,876 / 45,082, and the two strata that partition the cases.

    The abstract reports 4,747 cases and 52,881 controls. Those are recruitment
    figures; every odds ratio in this mirror was computed against the post-QC
    numbers below. Publishing the abstract's figures beside these statistics
    would put a denominator on the page that no statistic beside it used -- the
    exact over-claim this project ranks worse than publishing nothing.
    """
    by_stratum = {
        stratum: sorted(
            burden.filter(pl.col("cohort_stratum") == stratum)["n_cases"].unique().to_list()
        )
        for stratum in ("all", "syndromic", "nonsyndromic")
    }
    assert by_stratum == {"all": [3876], "syndromic": [1471], "nonsyndromic": [2405]}
    assert by_stratum["syndromic"][0] + by_stratum["nonsyndromic"][0] == by_stratum["all"][0]
    assert burden["n_controls"].unique().to_list() == [45082]
