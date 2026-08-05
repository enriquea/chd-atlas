# src/chd_atlas/build/burden.py
"""Per-gene rare-variant burden, loaded from `mirrors/burden.tsv`.

Loading and grouping only. The bundle shape is `burden_payload`'s and the page
markup is `pages._burden_section`'s, the same split `omics.py` has with
`bundles.py`: this module knows the mirror, and nothing here knows HTML.

**Nothing filters a row out.** Not the synonymous rows, which are the study's own
negative control and the reason this layer is worth publishing at all; not the
rows with no case carrier, which are evidence *against* enrichment; not the rows
whose odds ratio came back unbounded, which carry the strongest signals in the
data. A gene's rows reach the page exactly as the study computed them, and the
page's job is to label them, not to choose between them.

The one thing that is *not* published is anything this atlas computed. There is
no pooled p-value across studies and no per-gene summary verdict, because the
CHD literature reuses cohorts -- DDD contributed to several of the papers this
atlas plans to curate, though only one cited today names it -- so a pooled
statistic would count the same children twice. `shared_cohorts` is the
alternative: it names the overlap instead of pretending it away, and it
derives nothing beyond a set intersection over what each study already declared.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Final

from chd_atlas.build.emit import Json
from chd_atlas.models.cohort import Cohort
from chd_atlas.tables import BURDEN, read_table

# The order strata are presented in, everywhere. Not alphabetical, which would
# put `nonsyndromic` before `syndromic` and both before the whole-cohort figure a
# reader looks at first; not the mirror's sort order either, which is
# alphabetical for the same reason. This is the reading order.
STRATUM_ORDER: Final[tuple[str, ...]] = ("all", "syndromic", "nonsyndromic")

# Within a stratum: the result, then the other result, then the control that
# calibrates them. `synonymous` last is a deliberate reading order rather than an
# ordering of importance -- it is the row that tells you whether to believe the
# two above it.
CONSEQUENCE_ORDER: Final[tuple[str, ...]] = (
    "lof",
    "missense_damaging",
    "missense_all",
    "all_coding",
    "synonymous",
)


@dataclass(frozen=True)
class BurdenRow:
    """One published burden statistic, as `mirrors/burden.tsv` records it.

    A flat mirror of the schema rather than a reshaped view of it. The columns
    that constrain each other -- `comparator` against the control counts,
    `effect` against `effect_measure` and `effect_bound` -- are already checked
    by `validate/burden.py`, and reshaping them here into, say, a union type
    would move that check into the build, where a violation raises mid-write
    instead of being reported by the gate.
    """

    study: str
    gene: str
    cohort_stratum: str
    lesion_group: str | None
    variant_class: str
    consequence_class: str
    origin: str
    maf_max: float | None
    n_case_carriers: int
    n_cases: int
    comparator: str
    n_control_carriers: int | None
    n_controls: int | None
    expected_count: float | None
    effect: float | None
    effect_measure: str | None
    effect_bound: str | None
    ci_low: float | None
    ci_high: float | None
    pvalue: float | None
    pvalue_test: str | None
    pvalue_adjusted: float | None
    pvalue_adjustment: str | None
    case_cohorts: tuple[str, ...]
    control_cohorts: tuple[str, ...]
    method_note: str | None
    source: str


def _cohorts(value: object) -> tuple[str, ...]:
    """The `;`-joined cohort column as a tuple, in the order it was written.

    Not sorted: the converter writes them in the order the study's own cohort
    description names them, and a reader meeting "cnchd, ddd, nottingham"
    should meet the study's ordering rather than the alphabet's. The mirror cell
    is a single string, so this order is stable across builds by construction.
    """
    if value is None:
        return ()
    return tuple(part for part in str(value).split(";") if part)


def _sort_key(row: BurdenRow) -> tuple[object, ...]:
    """Reading order within one gene: study, then stratum, then consequence.

    `STRATUM_ORDER` and `CONSEQUENCE_ORDER` are looked up by index rather than
    compared as strings, so the sequence a reader sees is the one declared above
    and not the alphabet's. A value outside either tuple sorts last rather than
    raising: it is a vocabulary member `validate_table` accepted and this module
    has not been taught to order, and dropping it or dying on it would be worse
    than putting it at the end.
    """
    stratum = (
        STRATUM_ORDER.index(row.cohort_stratum)
        if row.cohort_stratum in STRATUM_ORDER
        else len(STRATUM_ORDER)
    )
    consequence = (
        CONSEQUENCE_ORDER.index(row.consequence_class)
        if row.consequence_class in CONSEQUENCE_ORDER
        else len(CONSEQUENCE_ORDER)
    )
    return (row.study, stratum, row.cohort_stratum, consequence, row.consequence_class, row.origin)


def load_burden(root: Path) -> dict[str, list[BurdenRow]]:
    """Every burden row, grouped by HGNC id and ordered for reading.

    Returns an empty mapping when the mirror is absent -- there is none before
    any is curated -- and when it is unreadable or missing a column, which
    `validate_table` reports as TBL000/TBL001 and `build_site` refuses on, so
    neither is reachable behind the gate.

    The tie-break in `_sort_key` runs to `origin` and stops there. That is not
    the mirror's full key: `lesion_group`, `variant_class` and `maf_max` are
    excluded, and today every row shares one value of each, so the order is
    total. When a second frequency threshold or a CNV row arrives, two rows can
    compare equal here and `sorted` is stable, so they will hold the mirror's
    order -- which `validate_sort_order` pins. Order stays deterministic either
    way; it is only the *reading* order that would need extending.
    """
    path = root / "mirrors" / "burden.tsv"
    if not path.is_file():
        return {}
    frame, _ = read_table(path, BURDEN)
    if frame is None or any(column not in frame.columns for column in BURDEN.column_names):
        return {}

    grouped: dict[str, list[BurdenRow]] = {}
    for record in frame.iter_rows(named=True):
        row = BurdenRow(
            study=str(record["study"]),
            gene=str(record["gene"]),
            cohort_stratum=str(record["cohort_stratum"]),
            lesion_group=record["lesion_group"],
            variant_class=str(record["variant_class"]),
            consequence_class=str(record["consequence_class"]),
            origin=str(record["origin"]),
            maf_max=record["maf_max"],
            n_case_carriers=int(record["n_case_carriers"]),
            n_cases=int(record["n_cases"]),
            comparator=str(record["comparator"]),
            n_control_carriers=record["n_control_carriers"],
            n_controls=record["n_controls"],
            expected_count=record["expected_count"],
            effect=record["effect"],
            effect_measure=record["effect_measure"],
            effect_bound=record["effect_bound"],
            ci_low=record["ci_low"],
            ci_high=record["ci_high"],
            pvalue=record["pvalue"],
            pvalue_test=record["pvalue_test"],
            pvalue_adjusted=record["pvalue_adjusted"],
            pvalue_adjustment=record["pvalue_adjustment"],
            case_cohorts=_cohorts(record["case_cohorts"]),
            control_cohorts=_cohorts(record["control_cohorts"]),
            method_note=record["method_note"],
            source=str(record["source"]),
        )
        grouped.setdefault(row.gene, []).append(row)

    return {gene: sorted(rows, key=_sort_key) for gene, rows in grouped.items()}


def shared_cohorts(rows: Iterable[BurdenRow]) -> dict[tuple[str, str], tuple[str, ...]]:
    """Which sample collections each pair of studies has in common.

    The reason `curation/cohorts.yaml` exists. The atlas publishes per-study rows
    and no pooled statistic, because pooling studies that share samples inflates
    the result -- but a reader looking at two rows will compare them anyway, and
    nothing on the page would tell them the two are not independent. This names
    the overlap so the page can.

    Keyed on the ordered pair with the smaller study id first, so one pair yields
    one entry however the rows are encountered, and the mapping iterates in a
    stable order. Pairs sharing nothing are omitted: the interesting fact is the
    overlap, and an entry per independent pair would bury it.

    **Latent on the committed corpus**, which carries one study. Implemented and
    tested against a fixture rather than left until the second study lands,
    because a feature nobody exercised is a feature that breaks on the day it is
    first needed.
    """
    by_study: dict[str, set[str]] = {}
    for row in rows:
        by_study.setdefault(row.study, set()).update(row.case_cohorts, row.control_cohorts)

    overlaps: dict[tuple[str, str], tuple[str, ...]] = {}
    studies = sorted(by_study)
    for index, left in enumerate(studies):
        for right in studies[index + 1 :]:
            common = by_study[left] & by_study[right]
            if common:
                overlaps[(left, right)] = tuple(sorted(common))
    return overlaps


def burden_payload(rows: Iterable[BurdenRow]) -> list[Json]:
    """One gene's burden rows, JSON-ready, in the order they will be rendered.

    Every key is present on every row, `null` where the comparator does not
    populate it, for the reason `bundles._validity_record` publishes ClinGen-only
    and GenCC-only fields as `null` rather than omitting them: an array whose
    objects differ in shape is a trap for a consumer that reads a field off
    `rows[0]` and expects the same key on `rows[1]`.

    Cohorts publish as arrays rather than as the mirror's `;`-joined string. The
    separator is a storage detail of a TSV that keeps one value per cell; a
    consumer should not have to know it, and a client splitting on `;` itself is
    a client that will one day meet a cohort id containing one.
    """
    return [
        {
            "study": row.study,
            "cohort_stratum": row.cohort_stratum,
            "lesion_group": row.lesion_group,
            "variant_class": row.variant_class,
            "consequence_class": row.consequence_class,
            "origin": row.origin,
            "maf_max": row.maf_max,
            "n_case_carriers": row.n_case_carriers,
            "n_cases": row.n_cases,
            "comparator": row.comparator,
            "n_control_carriers": row.n_control_carriers,
            "n_controls": row.n_controls,
            "expected_count": row.expected_count,
            "effect": row.effect,
            "effect_measure": row.effect_measure,
            "effect_bound": row.effect_bound,
            "ci_low": row.ci_low,
            "ci_high": row.ci_high,
            "pvalue": row.pvalue,
            "pvalue_test": row.pvalue_test,
            "pvalue_adjusted": row.pvalue_adjusted,
            "pvalue_adjustment": row.pvalue_adjustment,
            "case_cohorts": list(row.case_cohorts),
            "control_cohorts": list(row.control_cohorts),
            "method_note": row.method_note,
            "source": row.source,
        }
        for row in rows
    ]


def cohort_registry(cohorts: Iterable[Cohort]) -> dict[str, Cohort]:
    """Cohort id -> the curated record, from `corpus.cohorts`.

    Returns the whole record rather than just the name. The first version
    returned `{id: name}` and dropped `description` -- which is where
    `models/cohort.py` says the caveats qualifying every number drawn from a
    collection belong, and where a curator had written that UK Biobank's
    controls are adults recruited at 40-69 while the CHD cases were largely
    enrolled in childhood. Measured 2026-08-05: no cohort description reached a
    single published byte. The curator wrote the caveat, the model documented
    where it goes, and no reader could reach it.

    A cohort absent from the registry keeps its id as its label at the render
    site -- the same fallback `build_genes` applies to a gene missing from
    `mirrors/genes.tsv`. `validate_burden_references` reports the absence as
    BUR009 and the gate refuses on it, so that fallback is unreachable behind
    `build_site`.
    """
    return {str(cohort.id): cohort for cohort in cohorts}
