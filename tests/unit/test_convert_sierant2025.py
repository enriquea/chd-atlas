# tests/unit/test_convert_sierant2025.py
"""The Sierant 2025 converter's contracts, checked against the real mirror.

Same shape as `test_convert_burden.py` and for the same reason: what can pass
unnoticed here is not a crash but a right-looking row holding the wrong number.

This converter has one failure mode the others do not. Dataset S6 has a
three-level merged header, so every column is addressed by a bare integer index
into a positional map. A supplement re-uploaded with one column inserted would
transcribe a D-Mis p-value into a loss-of-function row, and every downstream gate
would pass, because the numbers are individually plausible. `_check_layout` is
the guard, and `test_a_shifted_workbook_is_refused_rather_than_transcribed` is
what proves it fires.
"""

from __future__ import annotations

import sys
from pathlib import Path

import polars as pl
import pytest

from chd_atlas.tables import BURDEN, read_table

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from convert_sierant2025 import (  # noqa: E402
    _BLOCKS,
    _EXPECTED_HEADERS,
    _HEADER_ROW,
    GENE_COLUMN,
    _case_control_row,
    _check_layout,
    _de_novo_row,
)

STUDY = "PMID:40127276"


@pytest.fixture(scope="module")
def burden() -> pl.DataFrame:
    """Only this converter's rows; the mirror holds three studies."""
    frame, issues = read_table(REPO_ROOT / "mirrors" / "burden.tsv", BURDEN)
    assert frame is not None, issues
    mine = frame.filter(pl.col("study") == STUDY)
    assert mine.height, f"no rows for {STUDY}; did the merge drop them?"
    return mine


def test_a_shifted_workbook_is_refused_rather_than_transcribed() -> None:
    """The guard that stops a re-uploaded supplement being silently relabelled.

    `_BLOCKS` maps consequence and arm to bare column integers. Insert one column
    into the sheet and every index after it reads its neighbour: the numbers stay
    plausible, the schema stays satisfied, `validate_burden` stays clean, and the
    page publishes a D-Mis statistic under a loss-of-function heading.

    Nothing downstream can catch that, because there is nothing wrong with the
    values -- only with which cell they came from. So the converter refuses at
    the door, and this is the test that proves the door is shut.
    """
    # A header row exactly as the real workbook has it, then the same row with a
    # single column spliced in at the front of the de novo block.
    good = [""] * 70
    for index, label in _REAL_HEADER.items():
        good[index] = label
    _check_layout(good)  # the real layout: no exception

    shifted = good[:7] + ["Obs"] + good[7:]
    with pytest.raises(SystemExit) as raised:
        _check_layout(shifted)
    # Names the offending columns rather than failing with a bare index error,
    # so a curator meeting this knows which supplement changed and where.
    assert "does not match the layout" in str(raised.value)
    assert "column 8" in str(raised.value)


def test_every_gene_carries_both_arms_for_all_three_consequence_classes(
    burden: pl.DataFrame,
) -> None:
    """Six rows per gene, and the partition that makes them distinct.

    Dataset S6 reports a de novo arm and a case-control arm for each of
    `damaging`, `lof` and `missense_damaging`. The two arms of one consequence
    differ in `comparator`, in `origin` and in `count_unit` but share a gene, a
    study and a stratum -- so if `origin` were dropped from the mirror's sort key
    the pair would collide, and one of the two would be the row that survived.

    Measured 2026-08-05: 30 registered genes x 3 consequences x 2 arms = 180.
    """
    assert burden.height == 180
    assert burden["gene"].n_unique() == 30
    assert set(burden["consequence_class"].unique()) == set(_BLOCKS)

    per_gene = burden.group_by("gene").len()["len"].unique().to_list()
    assert per_gene == [6]

    de_novo = burden.filter(pl.col("comparator") == "mutation_model")
    case_control = burden.filter(pl.col("comparator") == "control_cohort")
    assert de_novo.height == case_control.height == 90


def test_the_two_arms_count_different_things_and_say_so(burden: pl.DataFrame) -> None:
    """The distinction `count_unit` was added for, on the study that forced it.

    The de novo arm counts mutations against a denominator of trios, constant at
    3,887. The case-control arm counts alleles against per-gene allele totals
    that run above 20,000 for a cohort of 11,555 probands -- which is the
    arithmetic that proves they cannot both be people.
    """
    de_novo = burden.filter(pl.col("comparator") == "mutation_model")
    assert set(de_novo["count_unit"].unique()) == {"de_novo_mutations"}
    assert set(de_novo["origin"].unique()) == {"de_novo"}
    assert set(de_novo["n_cases"].unique()) == {3887}
    # A modelled expectation, and never an observed control count.
    assert de_novo["expected_count"].null_count() == 0
    assert de_novo["n_controls"].null_count() == de_novo.height

    case_control = burden.filter(pl.col("comparator") == "control_cohort")
    assert set(case_control["count_unit"].unique()) == {"alleles"}
    assert set(case_control["origin"].unique()) == {"transmitted_or_unphased"}
    assert case_control["n_cases"].max() > 2 * 10_000  # alleles, not the 11,555 probands
    assert case_control["expected_count"].null_count() == case_control.height


def test_the_composite_row_is_the_sum_of_its_two_components(burden: pl.DataFrame) -> None:
    """`damaging` is the union of `lof` and `missense_damaging`, not a third result.

    This is the check that the three blocks were read from the columns they were
    meant to be read from: the arithmetic only holds if `_BLOCKS` points at the
    right band of the sheet for each consequence. A transposed pair of indices
    would satisfy every schema rule and fail here.

    Verified on the counts, not the p-values -- the composite's p-value is *not*
    a function of its components', which is why the composite is mirrored at all
    rather than left to be reconstructed.
    """
    counts = {
        (row["gene"], row["consequence_class"], row["comparator"]): row["n_case_carriers"]
        for row in burden.iter_rows(named=True)
    }
    genes = sorted({gene for gene, _, _ in counts})
    assert len(genes) == 30

    for gene in genes:
        for comparator in ("mutation_model", "control_cohort"):
            composite = counts[(gene, "damaging", comparator)]
            parts = (
                counts[(gene, "lof", comparator)] + counts[(gene, "missense_damaging", comparator)]
            )
            assert composite == parts, (gene, comparator, composite, parts)


def test_an_infinite_odds_ratio_publishes_no_number_and_still_names_its_measure() -> None:
    """Dataset S6 stores a literal `Inf` where no control allele carries.

    `encode_json` sets `allow_nan=False`, so infinity cannot be published; the
    row must carry a null effect plus `effect_bound`, and `pages._effect` then
    renders `OR ∞`. The measure survives because "unbounded" is as meaningless
    without it as "3.1" is.

    **Driven directly rather than through the mirror, because no committed row
    reaches this branch.** Measured 2026-08-05 over the 248 panel genes: exactly
    two `Inf` cells, RPL5 (LOF, 1 case allele of 21,901 against 0 of 270,942) and
    KLF4 (D-Mis, 1 of 20,558 against 0 of 234,502) -- and neither gene is in
    `mirrors/genes.tsv`, so the converter skips both and the mirror holds zero
    unbounded rows from this study. The assertion below said 2 until the test
    was run; the mirror said 0, and the mirror was right.

    A guard nobody has seen fire is a guard nobody knows works, so the branch is
    exercised on a synthetic record shaped like RPL5's real one.
    """
    record = [""] * 70
    record[GENE_COLUMN] = "RPL5"
    ix = _BLOCKS["lof"]
    record[ix["cc_obs"]], record[ix["cc_total"]] = "1", "21901"
    record[ix["cc_ctrl"]], record[ix["cc_ctrl_total"]] = "0", "270942"
    record[ix["cc_effect"]] = "Inf"
    record[ix["cc_p"]], record[ix["cc_fdr"]] = "0.0747875141287311", "0.3"

    row = _case_control_row("HGNC:10360", "lof", record, ix)

    assert row["effect"] == ""
    assert row["effect_bound"] == "unbounded_above"
    assert row["effect_measure"] == "odds_ratio"
    # This study publishes no confidence intervals at all, so the unbounded row
    # carries neither bound. `validate_burden` permits that only because
    # `effect_bound` explains the absence -- BUR014 is suspended exactly there.
    assert row["ci_low"] == row["ci_high"] == ""

    # A finite effect keeps its number and takes no bound, so the branch above is
    # a branch rather than the only path. The spelling shortens on the way
    # through -- the workbook stores KMT2D's odds ratio as
    # "3.5919514463315401" and the mirror writes "3.59195144633154", which is
    # the same IEEE-754 double at its shortest exact spelling. That is
    # `_number`'s contract, shared with `convert_burden`: a rendering choice, not
    # a rounding one.
    record[ix["cc_effect"]] = "3.5919514463315401"
    finite = _case_control_row("HGNC:10360", "lof", record, ix)
    assert finite["effect"] == "3.59195144633154"
    assert float(finite["effect"]) == float("3.5919514463315401")
    assert finite["effect_bound"] == ""


def test_no_committed_row_from_this_study_is_unbounded_or_carries_an_interval(
    burden: pl.DataFrame,
) -> None:
    """The mirror-side half of the case above, and a contract of its own.

    Both `Inf` genes are unregistered, so every committed row from this study has
    a finite effect. And the study publishes no confidence intervals for either
    arm, which is why `ci_low` and `ci_high` are null on all 180 -- not because
    they were dropped in transcription. BUR014 would report a half-interval, and
    BUR015 an interval qualifying nothing; neither fires, and this says why.
    """
    assert burden.filter(pl.col("effect_bound").is_not_null()).height == 0
    assert burden["ci_low"].null_count() == burden.height
    assert burden["ci_high"].null_count() == burden.height
    assert burden["effect"].null_count() == 0
    assert set(burden["effect_measure"].unique()) == {"odds_ratio", "enrichment_ratio"}


def test_the_corrected_p_is_the_studys_benjamini_hochberg_and_never_its_local_fdr(
    burden: pl.DataFrame,
) -> None:
    """Which of Dataset S6's two corrections is publishable, and why only one.

    Each arm carries a Benjamini-Hochberg FDR, which is a step-up adjusted p and
    is never smaller than the p it corrects -- BUR018 asserts exactly that, and
    it holds on all 180 rows.

    The sheet also carries a JL-FDR combining the two arms. **That one is a
    local false discovery rate, not an adjusted p**, and measured over the 248
    panel genes it is smaller than its own p-value for 118 damaging, 135 LOF and
    129 D-Mis genes. Putting it here would trip BUR018, correctly: it is not a
    corrected p and must not be rendered in a column headed "corrected p".
    """
    assert set(burden["pvalue_adjustment"].unique()) == {"benjamini_hochberg"}
    assert burden["pvalue_adjusted"].null_count() == 0
    assert burden.filter(pl.col("pvalue_adjusted") < pl.col("pvalue")).height == 0
    assert set(burden["pvalue_test"].unique()) == {"poisson", "fisher_exact"}


# Row 7 of the real workbook, as the converter's own expectation table describes
# it. Written out here rather than imported so that a change to `_EXPECTED_HEADERS`
# has to be made twice deliberately, instead of the guard silently agreeing with
# whatever it was changed to.
_REAL_HEADER = {
    2: "pLI",
    3: "MisZ",
    5: "LOF",
    6: "MisD",
    7: "Obs",
    8: "Exp",
    9: "Enrich",
    10: "P-Val",
    11: "FDR",
    13: "LOF",
    14: "MisD",
    15: "Obs",
    16: "Total",
    17: "AF",
    18: "LOF",
    19: "MisD",
    20: "Obs",
    21: "Total",
    22: "AF",
    23: "Enrich",
    24: "P-Val",
    25: "FDR",
    27: "P-Val",
    28: "FDR",
    30: "Obs",
    31: "Exp",
    32: "Enrich",
    33: "P-Val",
    34: "FDR",
    36: "Obs",
    37: "Total",
    38: "AF",
    39: "Obs",
    40: "Total",
    41: "AF",
    42: "Enrich",
    43: "P-Val",
    44: "FDR",
    46: "P-Val",
    47: "FDR",
    49: "Obs",
    50: "Exp",
    51: "Enrich",
    52: "P-Val",
    53: "FDR",
    55: "Obs",
    56: "Total",
    57: "AF",
    58: "Obs",
    59: "Total",
    60: "AF",
    61: "Enrich",
    62: "P-Val",
    63: "FDR",
    65: "P-Val",
    66: "FDR",
}


def test_the_header_row_index_is_the_one_the_expectation_table_describes() -> None:
    """`_HEADER_ROW` and `_EXPECTED_HEADERS` have to agree about which row is which.

    Trivial, and here because the two are separate constants: pointing
    `_HEADER_ROW` at row 6 would compare the leaf labels against the *band*
    labels, `_check_layout` would reject the real workbook, and the failure would
    read as a changed supplement rather than as a changed constant.
    """
    assert _HEADER_ROW == 7
    good = [""] * 70
    for index, label in _REAL_HEADER.items():
        good[index] = label
    _check_layout(good)


def _synthetic(**cells: str) -> list[str]:
    """A DS6 data row with only the cells a test names, positioned as in the sheet."""
    record = [""] * 70
    record[GENE_COLUMN] = "KMT2D"
    for key, value in cells.items():
        record[_BLOCKS["lof"][key]] = value
    return record


def test_the_two_arms_are_built_with_the_units_and_denominators_they_measured() -> None:
    """Driven through the row builders, not read back off the mirror.

    **The mirror-reading tests above cannot catch a converter defect**, and that
    is not a hypothetical: `mirrors/burden.tsv` is a committed artifact, so
    mutating `convert_sierant2025.py` changes nothing any of them observes. A
    mutation matrix on 2026-08-05 relabelled the de novo arm `individuals` and
    every one of those tests still passed.

    So the labelling contract is asserted where it is decided. Each field here is
    a claim the row makes about what its numbers mean:

    * the de novo arm counts mutations against 3,887 *trios*, carries a modelled
      expectation and no control column, and its effect is an enrichment;
    * the case-control arm counts alleles against a per-gene allele total,
      carries observed controls and no expectation, and its effect is an odds
      ratio.

    Swap either `count_unit`, either `origin`, or the trio denominator and this
    fails.
    """
    record = _synthetic(
        dnm_obs="21",
        dnm_exp="0.16558",
        dnm_effect="126.82",
        dnm_p="6.65e-37",
        dnm_fdr="1.6491e-34",
        cc_obs="12",
        cc_total="21768",
        cc_ctrl="35",
        cc_ctrl_total="227965",
        cc_effect="3.5919",
        cc_p="5.4935e-4",
        cc_fdr="1.2385e-2",
    )
    ix = _BLOCKS["lof"]

    de_novo = _de_novo_row("HGNC:7133", "lof", record, ix)
    assert de_novo["count_unit"] == "de_novo_mutations"
    assert de_novo["origin"] == "de_novo"
    assert de_novo["comparator"] == "mutation_model"
    assert (de_novo["n_case_carriers"], de_novo["n_cases"]) == ("21", "3887")
    assert de_novo["expected_count"] == "0.16558"
    assert de_novo["effect_measure"] == "enrichment_ratio"
    assert de_novo["pvalue_test"] == "poisson"
    # A mutation model observes no controls, and BUR001 refuses a row that
    # claims otherwise.
    assert de_novo["n_control_carriers"] == de_novo["n_controls"] == ""
    assert de_novo["control_cohorts"] == ""
    # No frequency filter: a de novo mutation is defined by absence from both
    # parents. Empty rather than 1.0 -- "no filter" and "filtered at fixation"
    # are different claims.
    assert de_novo["maf_max"] == ""

    case_control = _case_control_row("HGNC:7133", "lof", record, ix)
    assert case_control["count_unit"] == "alleles"
    assert case_control["origin"] == "transmitted_or_unphased"
    assert case_control["comparator"] == "control_cohort"
    assert (case_control["n_case_carriers"], case_control["n_cases"]) == ("12", "21768")
    assert (case_control["n_control_carriers"], case_control["n_controls"]) == ("35", "227965")
    assert case_control["effect_measure"] == "odds_ratio"
    assert case_control["pvalue_test"] == "fisher_exact"
    assert case_control["expected_count"] == ""
    assert case_control["maf_max"] == "1e-05"

    # Both arms publish the study's own Benjamini-Hochberg FDR, never the JL-FDR
    # local false discovery rate, which is smaller than its own p on more than
    # half the panel and is not a corrected p at all.
    for row in (de_novo, case_control):
        assert row["pvalue_adjustment"] == "benjamini_hochberg"
        assert row["case_cohorts"] == "pcgc"
        assert row["study"] == STUDY


def test_a_header_row_too_short_to_hold_the_layout_is_refused_by_name() -> None:
    """The short-row branch, which used to raise `IndexError` from its own message.

    `_check_layout` built the complaint string with `header[index]` in the same
    comprehension that tested `index >= len(header)`, so a truncated header --
    the shape a re-uploaded sheet with columns *removed* has -- died inside the
    error path rather than in it. A traceback naming neither the column nor the
    workbook is not a refusal a curator can act on.
    """
    with pytest.raises(SystemExit) as raised:
        _check_layout(["", ""])
    message = str(raised.value)
    assert "past the end of the row" in message
    assert "column 7" in message


def test_the_three_blocks_read_three_disjoint_bands_of_the_sheet() -> None:
    """`_BLOCKS` is 36 bare integers, and a transposition among them is silent.

    Point the damaging block's `dnm_obs` at column 30 -- the loss-of-function
    column -- and every gate stays green: the value is a real observed count, the
    schema accepts it, `validate_burden` accepts it, and the page publishes a
    loss-of-function count under a damaging heading. A mutation matrix on
    2026-08-05 confirmed it: that mutant SURVIVED the whole suite, because every
    other test here reads the committed mirror, which a converter mutation does
    not touch.

    Two structural properties close it, and neither needs the workbook:

    * **Every index is used exactly once.** A block reading another block's
      column necessarily duplicates one.
    * **Every index points at a header whose label matches its role.** `*_obs`
      must land on an `Obs` column, `*_exp` on `Exp`, `*_total` on `Total`,
      `*_effect` on `Enrich`, `*_p` on `P-Val` and `*_fdr` on `FDR`. This is what
      catches a swap *within* a block, which the uniqueness check alone would
      miss -- reading the expectation where the observation belongs keeps every
      index distinct.

    `_EXPECTED_HEADERS` is the shared reference, and `_check_layout` verifies it
    against the real sheet on every run, so the two halves meet.
    """
    used = [index for block in _BLOCKS.values() for index in block.values()]
    assert len(used) == 36
    assert len(set(used)) == 36, "two blocks read the same column"

    role_label = {
        "obs": "Obs",
        "exp": "Exp",
        "total": "Total",
        "ctrl": "Obs",
        "ctrl_total": "Total",
        "effect": "Enrich",
        "p": "P-Val",
        "fdr": "FDR",
    }
    for consequence, block in _BLOCKS.items():
        for key, index in block.items():
            role = key.split("_", 1)[1]
            assert _EXPECTED_HEADERS[index] == role_label[role], (
                f"{consequence}.{key} points at column {index}, whose header is "
                f"{_EXPECTED_HEADERS[index]!r} and not {role_label[role]!r}"
            )
