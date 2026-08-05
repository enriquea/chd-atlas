# scripts/convert_sierant2025.py
"""Convert the Sierant 2025 PCGC panel supplement into `mirrors/burden.tsv`.

Run locally by a curator, never by CI or the build -- design decision D6.

    .venv/bin/python scripts/convert_sierant2025.py \\
        local/data/sierant2025/pnas.2420343122.sd06.xlsx mirrors/burden.tsv

Merges beside the rows other studies already contributed rather than replacing
the file; see `burden_mirror.merge_study`.

    Sierant MC, Jin SC, Bilguvar K, et al. (2025) Genomic analysis of 11,555
    probands identifies 60 dominant congenital heart disease genes.
    PNAS 122(13): e2420343122.

Licence CC BY 4.0. Confirmed three ways on 2026-08-05, because PNAS uses
CC BY-NC-ND on some content: the article's own permissions statement
("This open access article is distributed under Creative Commons Attribution
License 4.0 (CC BY)"), Crossref's licence record for the DOI, and NCBI's OA
service, which returns `license="CC BY"` for PMC12002227.

WHY THIS STUDY AND NOT JIN 2017
-------------------------------
Its cases are the Pediatric Cardiac Genomics Consortium's CHD GENES cohort
(NCT01196182) -- the same consortium as Jin et al. 2017 (PMID:28991257), at
11,555 probands and 3,887 trios against Jin's ~2,871. **Jin's trios are a subset
of these.** Curating both would publish the same children twice under two study
ids, and `shared_cohorts` could not show it, because the overlap sits inside one
study's cohort list rather than between two. `curation/cohorts.yaml`'s `pcgc`
entry records that constraint so a later curator meets it before writing a row.

WHICH TABLE, AND WHY IT CARRIES THE NULLS
-----------------------------------------
Dataset S6: all 248 prespecified MIPseq panel genes, with both analytic arms for
each of three consequence classes. Counted 2026-08-05: 248 gene rows, no
duplicates. It includes the genes that showed nothing, which is what keeps this
mirror from being a list of hits -- the same reason `convert_audain2021.py` takes
S6 over S8 there.

Of the 248, 30 are in `mirrors/genes.tsv` and 15 are among the 23 genes the site
publishes: CHD7, FLT4, GATA4, GATA6, KDM6A, KMT2D, MYH6, NKX2-5, NODAL, NOTCH2,
SMAD2, TAB2, TBX5, ZEB2, ZIC3. The other eight published genes -- EVC, EVC2,
ISL1, NR2F2, PLD1, RBM10, TBX20, TFAP2B -- are not on the panel and will show
nothing from this study.

TWO ARMS, TWO COMPARATORS, TWO COUNT UNITS
------------------------------------------
The schema's `comparator` discriminator carries both arms without a migration,
and this is the first study to exercise `mutation_model` with real data:

* **De novo, `comparator: mutation_model`.** Observed de novo mutations in the
  3,887 trios against a mutability-based expectation, tested by Poisson.
  `count_unit: de_novo_mutations` -- the numerator counts mutations and the
  denominator counts *trios*.
* **Case-control, `comparator: control_cohort`.** Ultra-rare transmitted/unphased
  variants in all 11,555 probands against 133,743 gnomAD samples, one-tailed
  Fisher. `count_unit: alleles`.

**The case-control arm counts alleles, and that was measured rather than
assumed.** Its `Obs` equals the number of qualifying variant rows in Dataset S4
for 245 of 248 genes, and the number of distinct probands for only 235 --
CACNA1A 80 variants across 79 probands, TSC1 36 across 34, LRP1 135 across 132 --
and `Obs` follows the variant every time. The denominators are allele counts too:
they run to 23,022 against 11,555 probands, vary gene by gene with coverage, and
are not all even. See `CountUnit`.

`origin` is `transmitted_or_unphased`, not `inherited` and not `any`. The arm
removes de novo mutations, then pools variants whose transmission is known with
variants that were never phased; in Dataset S4, 27,429 of 36,054 variant rows are
`Unphased`, so the unphased majority is the whole difficulty.

THE THIRD CONSEQUENCE ROW IS THE UNION OF THE OTHER TWO
-------------------------------------------------------
Dataset S6 reports each arm for `Damaging (LOF + D-Mis)`, for `Loss of Function`
alone, and for `D-Mis` alone. All three are mirrored, because the composite is
the study's *primary* analysis -- its 60 genes are defined by it -- and its
p-value is not a function of its components'. Verified on KMT2D: the damaging
de novo `Obs` of 27 is the LOF 21 plus the D-Mis 6, and the damaging
case-control `Obs` of 108 is 12 plus 96.

`ConsequenceClass.DAMAGING` sorts first so the components read as a breakdown of
the headline rather than as a third independent finding, and `pages.py` states
the relationship outright.

WHAT IS NOT MIRRORED, AND WHY
-----------------------------
* **The JL-FDR meta-analysis columns.** They combine the two arms, and are a
  *local* false discovery rate rather than a step-up adjusted p: measured over
  the 248 genes, the JL-FDR is **smaller** than its own p-value for 118 damaging,
  135 LOF and 129 D-Mis genes. Putting it in `pvalue_adjusted` would trip BUR018,
  which is right to refuse it -- "corrected p" it is not.
* **`pLI` and `MisZ`** are gnomAD constraint scores, a property of the gene
  rather than a burden statistic of this study.
* **Dataset S7's TDT** is a third design again (transmission disequilibrium in
  trios) and needs no new columns, but it is a separate analysis with separate
  denominators and belongs in its own pass.

THE EFFECT COLUMN IS THE STUDY'S NUMBER, NOT A RECOMPUTED ONE
-------------------------------------------------------------
`Enrich` in the case-control arm is labelled `odds_ratio`, and the reading was
checked rather than assumed: against the published `Obs`/`Total` cells, the
sample odds ratio agrees with it to within a relative 5.2e-4 over 169 LOF genes
and 1.2e-4 over 205 D-Mis genes, while the allele-frequency rate ratio agrees far
less well and neither reproduces it exactly. The methods say only "one-tailed
Fisher's exact test", whose effect estimate is an odds ratio, so that is the
label; the residual is consistent with a conditional maximum-likelihood estimate
or with the displayed counts differing slightly from the tested 2x2 after the
study's platform harmonisation. **The mirror carries the study's printed value
unchanged** -- the atlas computes no effect size -- and at the three significant
figures `pages._fmt` renders, a 5e-4 relative difference is invisible anyway.

Two genes carry a literal `Inf`: RPL5 (LOF, 1 case allele of 21,901 against 0 of
270,942) and KLF4 (D-Mis, 1 of 20,558 against 0 of 234,502). Neither is
published. They take `effect_bound: unbounded_above` with a null effect, for the
reason `EffectBound` exists -- `allow_nan=False` refuses to publish infinity.
Neither has a confidence interval, because this study publishes none at all.

DENOMINATORS
------------
`n_cases` differs between the arms and is taken from each arm's own column, never
from the abstract: the de novo arm is 3,887 trios (a constant), and the
case-control arm's allele totals are per gene.
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path
from typing import Final

from burden_mirror import merge_study
from convert_burden import gene_symbols, read_sheet_rows

STUDY: Final = "PMID:40127276"
SOURCE: Final = "sierant2025_ds6"
SHEET: Final = "DS6. MetaAnalysis"

# `curation/cohorts.yaml`. One case cohort, because the study is one consortium.
CASE_COHORTS: Final = "pcgc"
CONTROL_COHORTS: Final = "gnomad_controls"

# "Additional WES increased totals to 3,887 CHD trios and 1,739 singletons, which
# combined with MIPseq probands yielded 11,555 CHD probands." The de novo arm
# uses the trios; the case-control arm's denominator is per gene and read from
# the sheet.
N_TRIOS: Final = 3887

# "We filtered for variants with a MAF <= 10-5 in both BRAVO and EVS databases."
TUV_MAF_MAX: Final = "1e-05"

# Row 8 onward is data; rows 0-7 are the title, the caption and a three-level
# header. Column indices are read off row 7 (the leaf header) against row 4 (the
# consequence block) and row 5 (the arm), and verified in `_check_layout` on
# every run rather than trusted -- a supplement re-uploaded with one column
# inserted would otherwise silently transcribe the wrong statistic.
FIRST_DATA_ROW: Final = 8
GENE_COLUMN: Final = 1

# consequence_class -> column indices for each arm.
#   dnm_*: Poisson test of de novo mutations in trios
#   cc_*:  Fisher case-control of ultra-rare TUVs against gnomAD
#
# `fmt: off` through `_HEADER_ROW`: these are a positional map into a
# spreadsheet, and the line grouping is the audit trail -- one line per arm,
# ascending, so a reader can check it against the workbook's own header bands by
# eye. Exploded one key per line, 60 lines of `"cc_obs": 15,` are unreviewable
# and a transposed pair would hide in them.
# fmt: off
_BLOCKS: Final[dict[str, dict[str, int]]] = {
    "damaging": {
        "dnm_obs": 7, "dnm_exp": 8, "dnm_effect": 9, "dnm_p": 10, "dnm_fdr": 11,
        "cc_obs": 15, "cc_total": 16, "cc_ctrl": 20, "cc_ctrl_total": 21,
        "cc_effect": 23, "cc_p": 24, "cc_fdr": 25,
    },
    "lof": {
        "dnm_obs": 30, "dnm_exp": 31, "dnm_effect": 32, "dnm_p": 33, "dnm_fdr": 34,
        "cc_obs": 36, "cc_total": 37, "cc_ctrl": 39, "cc_ctrl_total": 40,
        "cc_effect": 42, "cc_p": 43, "cc_fdr": 44,
    },
    "missense_damaging": {
        "dnm_obs": 49, "dnm_exp": 50, "dnm_effect": 51, "dnm_p": 52, "dnm_fdr": 53,
        "cc_obs": 55, "cc_total": 56, "cc_ctrl": 58, "cc_ctrl_total": 59,
        "cc_effect": 61, "cc_p": 62, "cc_fdr": 63,
    },
}

# What row 7 must read at each index for the mapping above to be the mapping the
# workbook actually has. Checked on every run; see `_check_layout`.
_EXPECTED_HEADERS: Final[dict[int, str]] = {
    GENE_COLUMN: "",  # row 7 is blank under "Gene", which sits on row 6
    7: "Obs", 8: "Exp", 9: "Enrich", 10: "P-Val", 11: "FDR",
    15: "Obs", 16: "Total", 20: "Obs", 21: "Total", 23: "Enrich", 24: "P-Val", 25: "FDR",
    30: "Obs", 31: "Exp", 32: "Enrich", 33: "P-Val", 34: "FDR",
    36: "Obs", 37: "Total", 39: "Obs", 40: "Total", 42: "Enrich", 43: "P-Val", 44: "FDR",
    49: "Obs", 50: "Exp", 51: "Enrich", 52: "P-Val", 53: "FDR",
    55: "Obs", 56: "Total", 58: "Obs", 59: "Total", 61: "Enrich", 62: "P-Val", 63: "FDR",
}
# fmt: on
_HEADER_ROW: Final = 7

METHOD_NOTE: Final = (
    "Cases are the Pediatric Cardiac Genomics Consortium CHD GENES cohort "
    "(NCT01196182), which subsumes the trios of the earlier PCGC exome papers. "
    "The de novo and case-control arms count different things: de novo rows "
    "count mutations against a denominator of trios, and case-control rows count "
    "alleles, so a proband carrying two qualifying variants is counted twice and "
    "the denominator is an allele count that varies per gene with coverage. The "
    "damaging row is the union of the loss-of-function and damaging-missense "
    "rows, not an independent third result."
)


def _check_layout(header: list[str]) -> None:
    """Refuse a workbook whose columns are not where `_BLOCKS` says they are.

    Every index in `_BLOCKS` is a bare integer into a three-level header, which
    is the most brittle thing in this converter: a supplement re-uploaded with
    one column inserted would transcribe a D-Mis p-value into a LOF row, and
    every downstream gate would pass, because the numbers are individually
    plausible. That is this project's characteristic failure with the sign
    flipped -- not evidence lost, but evidence relabelled.
    """
    # `found` is computed before the comparison rather than inside it: a header
    # short enough to be missing the index would otherwise raise `IndexError`
    # from the message-building expression itself, turning "this workbook is the
    # wrong shape" into a traceback that names neither the column nor the file.
    wrong = []
    for index, expected in sorted(_EXPECTED_HEADERS.items()):
        found = header[index].strip() if index < len(header) else "<past the end of the row>"
        if found != expected:
            wrong.append(f"column {index}: expected {expected!r}, found {found!r}")
    if wrong:
        raise SystemExit(
            "the workbook's header row does not match the layout this converter "
            "was written against, so every column index below is suspect:\n  " + "\n  ".join(wrong)
        )


def _cell(row: list[str], index: int) -> str:
    return row[index].strip() if index < len(row) else ""


def _number(row: list[str], index: int) -> str:
    """A published float at its shortest exact spelling; see `convert_burden`.

    `Inf` is returned as the empty string: the caller pairs that with
    `effect_bound: unbounded_above`, because `encode_json` sets `allow_nan=False`
    and infinity cannot be published as a number.
    """
    value = _cell(row, index)
    if not value or value.lower() in {"inf", "-inf", "nan", "na"}:
        return ""
    return repr(float(value))


def _integer(row: list[str], index: int) -> str:
    value = _cell(row, index)
    return str(int(float(value))) if value else ""


def _de_novo_row(
    gene: str, consequence: str, record: list[str], ix: dict[str, int]
) -> dict[str, str]:
    """The Poisson arm: observed de novo mutations against a modelled expectation."""
    return {
        "study": STUDY,
        "gene": gene,
        "cohort_stratum": "all",
        "lesion_group": "",
        "variant_class": "snv_indel",
        "consequence_class": consequence,
        "origin": "de_novo",
        # No frequency filter: a de novo mutation is defined by its absence from
        # both parents, not by a population threshold. Empty rather than 1.0 --
        # "no filter" and "filtered at fixation" are different claims.
        "maf_max": "",
        "count_unit": "de_novo_mutations",
        "n_case_carriers": _integer(record, ix["dnm_obs"]),
        "n_cases": str(N_TRIOS),
        "comparator": "mutation_model",
        "n_control_carriers": "",
        "n_controls": "",
        "expected_count": _number(record, ix["dnm_exp"]),
        "effect": _number(record, ix["dnm_effect"]),
        "effect_measure": "enrichment_ratio",
        "effect_bound": "",
        "ci_low": "",
        "ci_high": "",
        "pvalue": _number(record, ix["dnm_p"]),
        "pvalue_test": "poisson",
        "pvalue_adjusted": _number(record, ix["dnm_fdr"]),
        "pvalue_adjustment": "benjamini_hochberg",
        "case_cohorts": CASE_COHORTS,
        "control_cohorts": "",
        "method_note": METHOD_NOTE,
        "source": SOURCE,
    }


def _case_control_row(
    gene: str, consequence: str, record: list[str], ix: dict[str, int]
) -> dict[str, str]:
    """The Fisher arm: ultra-rare TUV alleles in probands against gnomAD."""
    effect = _number(record, ix["cc_effect"])
    unbounded = not effect and _cell(record, ix["cc_effect"]).lower() == "inf"
    return {
        "study": STUDY,
        "gene": gene,
        "cohort_stratum": "all",
        "lesion_group": "",
        "variant_class": "snv_indel",
        "consequence_class": consequence,
        "origin": "transmitted_or_unphased",
        "maf_max": TUV_MAF_MAX,
        "count_unit": "alleles",
        "n_case_carriers": _integer(record, ix["cc_obs"]),
        "n_cases": _integer(record, ix["cc_total"]),
        "comparator": "control_cohort",
        "n_control_carriers": _integer(record, ix["cc_ctrl"]),
        "n_controls": _integer(record, ix["cc_ctrl_total"]),
        "expected_count": "",
        "effect": effect,
        # Present even where the effect is not: an unbounded odds ratio still
        # has to say it is an odds ratio. `validate_burden` requires the measure
        # for either an effect or a bound, and `pages._effect` has no branch that
        # omits it.
        "effect_measure": "odds_ratio",
        "effect_bound": "unbounded_above" if unbounded else "",
        # This study publishes no confidence intervals, for either arm.
        "ci_low": "",
        "ci_high": "",
        "pvalue": _number(record, ix["cc_p"]),
        "pvalue_test": "fisher_exact",
        "pvalue_adjusted": _number(record, ix["cc_fdr"]),
        "pvalue_adjustment": "benjamini_hochberg",
        "case_cohorts": CASE_COHORTS,
        "control_cohorts": CONTROL_COHORTS,
        "method_note": METHOD_NOTE,
        "source": SOURCE,
    }


def convert(source: Path, symbols: dict[str, str]) -> tuple[list[dict[str, str]], list[str]]:
    """Build the mirror rows, and the report lines describing what was dropped."""
    sheet = read_sheet_rows(source, SHEET)
    _check_layout(sheet[_HEADER_ROW])

    records = [row for row in sheet[FIRST_DATA_ROW:] if _cell(row, GENE_COLUMN)]
    rows: list[dict[str, str]] = []
    seen: set[str] = set()

    for record in records:
        symbol = _cell(record, GENE_COLUMN)
        seen.add(symbol)
        hgnc = symbols.get(symbol)
        if hgnc is None:
            continue
        for consequence, ix in _BLOCKS.items():
            rows.append(_de_novo_row(hgnc, consequence, record, ix))
            rows.append(_case_control_row(hgnc, consequence, record, ix))

    kept = {row["gene"] for row in rows}
    missing = sorted(symbol for symbol, hgnc in symbols.items() if hgnc not in kept)
    report = [
        f"{len(rows)} row(s) for {len(kept)} of the {len(symbols)} registered genes",
        f"{len(seen)} genes on the panel; {len(seen) - len(kept)} not in mirrors/genes.tsv",
    ]
    if missing:
        report.append(
            f"WARNING: {len(missing)} registered gene(s) are not on this study's 248-gene "
            f"panel, so they will show no burden from it: {', '.join(missing)}"
        )
    return rows, report


def main(argv: list[str]) -> int:
    if not 2 <= len(argv) <= 3:
        print(__doc__, file=sys.stderr)
        return 2
    source, out = Path(argv[0]), Path(argv[1])
    root = Path(argv[2]) if len(argv) == 3 else Path.cwd()

    rows, report = convert(source, gene_symbols(root))
    kept, written = merge_study(out, STUDY, rows)

    for line in report:
        print(line)
    print(f"wrote {written} row(s) for {STUDY}, keeping {kept} row(s) from other studies")
    print(f"sha256 of {source.name}: {hashlib.sha256(source.read_bytes()).hexdigest()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
