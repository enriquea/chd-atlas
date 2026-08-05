# scripts/convert_audain2021.py
"""Convert the Audain 2021 CNV-deletion burden supplement into `mirrors/burden.tsv`.

Run locally by a curator, never by CI or the build -- design decision D6.

    .venv/bin/python scripts/convert_audain2021.py \\
        local/data/audain2021/s011.xlsx mirrors/burden.tsv

Merges beside the rows other studies already contributed rather than replacing
the file; see `burden_mirror.merge_study`.

    Audain E, Wilsdon A, Breckpot J, et al. (2021) Integrative analysis of
    genomic variants reveals new associations of candidate haploinsufficient
    genes with congenital heart disease. PLoS Genet 17(7): e1009679.

Licence CC0 -- read off the article, which states it is "free of all copyright
... made available under the Creative Commons CC0 public domain dedication".
Not assumed from the publisher: PLOS uses CC BY on other titles.

THE ERRATUM WAS READ FIRST
--------------------------
PMID:34547032 (PLoS Genet 17(9): e1009809) corrects this paper. It reads, in
full: "The thirteenth author's name is spelled incorrectly. The correct name is:
D. Woodrow Benson." No data is corrected, so every number below stands. Checked
before a single figure was transcribed, because an erratum on a burden paper
could as easily have corrected the counts.

WHICH TABLE, AND WHY IT IS NOT THE OBVIOUS ONE
----------------------------------------------
S6 Table Sheet A: all 11,515 genes PLINK tested, with `case`, `control`, EMP1
and EMP2. S8 Table has something better in one respect -- it splits the case set
into CHD-only and TAA-only -- and something disqualifying in another: **it
contains only the 528 genes that reached significance.**

Measured 2026-08-05 against the registry: S6 covers 103 of the 154 in-scope
genes and 13 of the 23 published ones; S8 covers 12 and 4. Every one of S8's
four (GATA4, NKX2-5, NR2F2, TAB2) is significant by construction, while S6 also
carries the nulls that qualify them -- CHD7 at EMP1 0.0068 and EMP2 0.99, EVC2
at 0.22, ISL1 at 0.34. Mirroring S8 alone would have put a significant CNV
result on every gene page that showed one, which is publication bias built into
the atlas rather than merely reported by it.

So S6 it is, and with it the cost below.

THE CASE SET IS NOT PURELY CHD
------------------------------
S6's `case` column counts CHD **and** sporadic thoracic aortic aneurysm
probands together. The study says so in its abstract -- "probands with CHD
including cases with sporadic thoracic aortic aneurysm" -- and S2 Table puts
the TAA series (PMID 21092924) at 777 of the 7,958 contributing individuals,
about 10%. `curation/chd_scope.yaml` places aortopathy outside this atlas.

Every row therefore names `taa_cases` in `case_cohorts`, whose registry entry
states the exclusion in as many words, and carries a `method_note` saying it on
the row itself. Nothing here silently relabels a CHD+TAA analysis as CHD.

WHY THE DE NOVO TABLE IS NOT MIRRORED
-------------------------------------
S10 Table (gene-based DNV analysis, 3,258 genes, 2,489 trios) fits this schema
technically -- it even carries `n_syn`, so the synonymous negative control would
work. It is deliberately left out.

**Its de novo variants are Jin 2017's and Sifrim 2016's.** The paper describes
itself as "a genome-wide integrative meta-analysis of published and publicly
available datasets of CNVs and DNVs", says "the data used in this study have
already been published", and cites Sifrim 2016 as reference [6] and Jin 2017 as
reference [7]. Those are the next two studies on this atlas's roadmap. Mirroring
S10 now and their own tables later would publish the same de novo variants twice
under two study ids, with the cohort columns unable to show it -- exactly the
double counting design decision D33 forbids, arriving through the one route the
overlap machinery cannot see.

The CNV part does not have that problem, and the difference is the point. S2
Table's ~190 sources are mostly single-case reports plus four collections this
atlas will never curate as studies in their own right, so aggregating them
creates no future double count. The de novo part aggregates two studies this
atlas intends to carry individually.

When Jin 2017 and Sifrim 2016 land, their own de novo tables are the right
representation, and their cohorts will name the overlap between them properly.

DENOMINATORS
------------
S1 Table gives the per-CNV-type counts: 7,107 cases and 14,040 controls have
deletion data, against 7,958 and 14,082 for CNVs of any type. The deletion
figures are the ones this analysis was computed on, so they are the ones the
mirror carries -- the same rule applied to PMID:42230622, whose abstract
likewise quotes a larger set than its statistics used.
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path
from typing import Final

from burden_mirror import merge_study
from convert_burden import gene_symbols, read_sheet

STUDY: Final = "PMID:34324492"
SOURCE: Final = "audain2021_s6"

# S2/S3 Tables. The ~190 single-case PMID sources are not registered: they are
# individual published case reports, not collections, and a registry entry per
# one would bury the four that carry the cohort-overlap signal.
CASE_COHORTS: Final = "decipher;isca;signature;taa_cases"
CONTROL_COHORTS: Final = "gain_controls;ottawa_controls;wtccc_controls"

# S1 Table, the deletion column.
N_CASES: Final = 7107
N_CONTROLS: Final = 14040

METHOD_NOTE: Final = (
    "Cases combine congenital heart disease with 777 sporadic thoracic aortic "
    "aneurysm probands (of 7,958 contributing individuals); the study reports "
    "the CHD-only contribution separately for significant genes. A carrier is "
    "an individual with a rare deletion (MAF below 0.01) overlapping the gene."
)


def convert(source: Path, symbols: dict[str, str]) -> tuple[list[dict[str, str]], list[str]]:
    """Build the mirror rows, and the report lines describing what was dropped."""
    records = read_sheet(source, "Sheet A")
    rows: list[dict[str, str]] = []
    seen: set[str] = set()

    for record in records:
        symbol = record["symbol"]
        seen.add(symbol)
        hgnc = symbols.get(symbol)
        if hgnc is None:
            continue
        rows.append(
            {
                "study": STUDY,
                "gene": hgnc,
                # Not `all`: `CohortStratum.ALL` means a study's undivided case
                # set, and this one is undivided in a way that includes a
                # disease the atlas excludes. The stratum stays `all` because
                # that is what the study analysed; `taa_cases` and the method
                # note are what say the set is not purely CHD.
                "cohort_stratum": "all",
                "lesion_group": "",
                "variant_class": "cnv_deletion",
                # A deletion removes the gene; the study frames it as
                # haploinsufficiency throughout, and PLINK counts individuals
                # with an overlapping deletion rather than a coding consequence.
                "consequence_class": "lof",
                "origin": "any",
                # "Rare CNV deletions (MAF 0.01)" -- S15 Table's own title.
                "maf_max": "0.01",
                "n_case_carriers": record["case"],
                "n_cases": str(N_CASES),
                "comparator": "control_cohort",
                "n_control_carriers": record["control"],
                "n_controls": str(N_CONTROLS),
                "expected_count": "",
                # PLINK's CNV association test reports no effect size and no
                # confidence interval, so none is invented. `_effect` renders an
                # em dash, and the counts beside it carry the finding.
                "effect": "",
                "effect_measure": "",
                "effect_bound": "",
                "ci_low": "",
                "ci_high": "",
                "pvalue": _number(record["p_value"]),
                "pvalue_test": "permutation",
                # EMP2, PLINK's max(T) family-wise corrected p over the same
                # permutations that produced EMP1. This study publishes its own
                # correction, which PMID:42230622 does not.
                "pvalue_adjusted": _number(record["p_value_adj"]),
                "pvalue_adjustment": "familywise_permutation",
                "case_cohorts": CASE_COHORTS,
                "control_cohorts": CONTROL_COHORTS,
                "method_note": METHOD_NOTE,
                "source": SOURCE,
            }
        )

    kept = {row["gene"] for row in rows}
    missing = sorted(symbol for symbol, hgnc in symbols.items() if hgnc not in kept)
    report = [
        f"{len(rows)} row(s) for {len(kept)} of the {len(symbols)} registered genes",
        f"{len(seen)} genes tested in the supplement; {len(seen) - len(kept)} out of scope",
    ]
    if missing:
        report.append(
            f"WARNING: {len(missing)} registered gene(s) were not tested by this study, "
            f"so they will show no CNV burden from it: {', '.join(missing)}"
        )
    return rows, report


def _number(value: str) -> str:
    """A published float at its shortest exact spelling; see `convert_burden`."""
    return repr(float(value))


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
