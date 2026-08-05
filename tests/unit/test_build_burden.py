# tests/unit/test_build_burden.py
"""Loading, reading order, and the overlap the atlas refuses to pool away."""

from __future__ import annotations

from pathlib import Path

from chd_atlas.build.burden import (
    BurdenRow,
    burden_payload,
    cohort_labels,
    load_burden,
    shared_cohorts,
)
from chd_atlas.models.cohort import Cohort
from chd_atlas.tables import BURDEN

_BASE: dict[str, str] = {
    "study": "PMID:42230622",
    "gene": "HGNC:17075",
    "cohort_stratum": "all",
    "lesion_group": "",
    "variant_class": "snv_indel",
    "consequence_class": "lof",
    "origin": "any",
    "maf_max": "0.001",
    "n_case_carriers": "6",
    "n_cases": "3876",
    "comparator": "control_cohort",
    "n_control_carriers": "0",
    "n_controls": "45082",
    "expected_count": "",
    "effect": "",
    "effect_measure": "odds_ratio",
    "effect_bound": "unbounded_above",
    "ci_low": "13.7",
    "ci_high": "",
    "pvalue": "2.45e-07",
    "pvalue_test": "fisher_exact",
    "case_cohorts": "cnchd;ddd;nottingham",
    "control_cohorts": "ukbb",
    "method_note": "",
    "source": "audain2026_sd3",
}


def _write(tmp_path: Path, *rows: dict[str, str]) -> Path:
    mirrors = tmp_path / "mirrors"
    mirrors.mkdir(parents=True, exist_ok=True)
    header = "\t".join(BURDEN.column_names)
    body = ["\t".join(row[name] for name in BURDEN.column_names) for row in rows]
    (mirrors / "burden.tsv").write_text("\n".join([header, *body]) + "\n", encoding="utf-8")
    return tmp_path


def test_rows_read_in_the_order_a_reader_needs_not_the_order_they_sort_in(
    tmp_path: Path,
) -> None:
    """Two orderings that both differ from the alphabet, for different reasons.

    Strata: the whole cohort first, then its two partitions. Alphabetically
    `nonsyndromic` would come before `syndromic` and both before `all`, which
    puts the headline figure last.

    Consequence: result, other result, then the control.

    **`all_coding` is in this fixture because without it the test proved
    nothing.** Measured 2026-08-05: with only `lof`, `missense_damaging` and
    `synonymous` -- the three classes the committed mirror carries -- a mutant
    collapsing `CONSEQUENCE_ORDER` to a constant survived, because those three
    happen to sort into reading order alphabetically anyway. `all_coding` sorts
    first as a string and fourth in reading order, so it is what makes the
    ordering rule observable at all.
    """
    written = [
        {**_BASE, "cohort_stratum": stratum, "consequence_class": consequence}
        for stratum in ("syndromic", "nonsyndromic", "all")
        for consequence in ("synonymous", "all_coding", "missense_damaging", "lof")
    ]
    rows = load_burden(_write(tmp_path, *written))["HGNC:17075"]

    assert [(row.cohort_stratum, row.consequence_class) for row in rows] == [
        ("all", "lof"),
        ("all", "missense_damaging"),
        ("all", "all_coding"),
        ("all", "synonymous"),
        ("syndromic", "lof"),
        ("syndromic", "missense_damaging"),
        ("syndromic", "all_coding"),
        ("syndromic", "synonymous"),
        ("nonsyndromic", "lof"),
        ("nonsyndromic", "missense_damaging"),
        ("nonsyndromic", "all_coding"),
        ("nonsyndromic", "synonymous"),
    ]


def test_the_cohort_columns_are_split_but_not_reordered(tmp_path: Path) -> None:
    """`;` is a TSV storage detail; the order within it is the study's own.

    Sorting them would put `cnchd` before `ddd` before `nottingham` here by
    coincidence, so the fixture writes them in an order the alphabet disagrees
    with. A consumer reads arrays and never has to know the separator.
    """
    row = load_burden(
        _write(tmp_path, {**_BASE, "case_cohorts": "nottingham;cnchd;ddd"})
    )["HGNC:17075"][0]

    assert row.case_cohorts == ("nottingham", "cnchd", "ddd")
    assert row.control_cohorts == ("ukbb",)
    assert burden_payload([row])[0]["case_cohorts"] == ["nottingham", "cnchd", "ddd"]  # type: ignore[index,call-overload]


def test_two_studies_drawing_on_one_collection_are_named_as_not_independent() -> None:
    """The whole reason `curation/cohorts.yaml` exists.

    **Latent on the committed corpus**, which carries one study, so this is the
    only place the behaviour is exercised until a second lands. DDD contributes
    cases to more than one of the papers this atlas will cite, and the atlas
    computes no pooled statistic precisely because of that -- but a reader
    comparing two rows by eye needs to be told, and nothing else on the page
    would tell them.

    The third study shares nothing and is deliberately absent from the result:
    the interesting fact is the overlap, and an entry per independent pair would
    bury it. The key is the ordered pair, so one pair yields one entry however
    the rows are encountered.
    """
    rows = [
        _row(study="PMID:1", case_cohorts=("cnchd", "ddd"), control_cohorts=("ukbb",)),
        _row(study="PMID:2", case_cohorts=("ddd", "pcgc"), control_cohorts=()),
        _row(study="PMID:3", case_cohorts=("gel",), control_cohorts=()),
    ]

    assert shared_cohorts(rows) == {("PMID:1", "PMID:2"): ("ddd",)}


def _row(**overrides: object) -> BurdenRow:
    payload: dict[str, object] = {
        "study": "PMID:1",
        "gene": "HGNC:17075",
        "cohort_stratum": "all",
        "lesion_group": None,
        "variant_class": "snv_indel",
        "consequence_class": "lof",
        "origin": "any",
        "maf_max": 0.001,
        "n_case_carriers": 6,
        "n_cases": 3876,
        "comparator": "control_cohort",
        "n_control_carriers": 0,
        "n_controls": 45082,
        "expected_count": None,
        "effect": None,
        "effect_measure": "odds_ratio",
        "effect_bound": "unbounded_above",
        "ci_low": 13.7,
        "ci_high": None,
        "pvalue": 2.45e-07,
        "pvalue_test": "fisher_exact",
        "case_cohorts": ("cnchd",),
        "control_cohorts": ("ukbb",),
        "method_note": None,
        "source": "audain2026_sd3",
    }
    payload.update(overrides)
    return BurdenRow(**payload)  # type: ignore[arg-type]


def test_every_published_row_carries_every_key_whatever_its_comparator() -> None:
    """An array whose objects differ in shape is a trap for a consumer reading a
    field off `rows[0]` and expecting it on `rows[1]`.

    A case-control row and a de novo row populate disjoint halves of the
    comparator block, and both publish all 24 keys with `null` where the
    comparator does not reach. Same rule `bundles._validity_record` follows for
    ClinGen-only and GenCC-only fields.
    """
    case_control = _row()
    de_novo = _row(
        comparator="mutation_model",
        n_control_carriers=None,
        n_controls=None,
        control_cohorts=(),
        expected_count=0.42,
        effect=9.5,
        effect_measure="enrichment_ratio",
        effect_bound=None,
        ci_high=24.1,
        pvalue_test="poisson",
    )
    payload = burden_payload([case_control, de_novo])

    assert set(payload[0]) == set(payload[1])  # type: ignore[arg-type]
    assert len(payload[0]) == 24  # type: ignore[arg-type]
    assert payload[1]["expected_count"] == 0.42  # type: ignore[index,call-overload]
    assert payload[1]["n_controls"] is None  # type: ignore[index,call-overload]


def test_an_absent_mirror_loads_as_no_burden_rather_than_failing(tmp_path: Path) -> None:
    """There is no burden mirror before any is curated, and every gene page has
    to render without one. `_burden_section` returns "" for such a gene rather
    than an empty table, which would read as "studied, nothing found".
    """
    assert load_burden(tmp_path) == {}


def test_a_cohort_missing_from_the_registry_keeps_its_id_as_its_label() -> None:
    """The same fallback `build_genes` applies to a gene missing from
    `mirrors/genes.tsv`, and for the same reason: an id is something a reader can
    look up, and a blank is a broken page. Unreachable behind `build_site`, which
    refuses on the BUR009 the absence reports.
    """
    labels = cohort_labels(
        [Cohort(id="ddd", name="Deciphering Developmental Disorders (DDD)", description="d")]
    )

    assert labels == {"ddd": "Deciphering Developmental Disorders (DDD)"}
    assert labels.get("pcgc", "pcgc") == "pcgc"
