# tests/unit/test_validate_profiles.py
"""The rules that stop a gene's percentile being read against another
dataset's ruler.

Every fixture here is synthetic rather than transcribed from a real study --
unlike `test_validate_burden.py`, which builds on one real row -- because no
`profiles`/`profile_quantiles` mirror has ever been committed (measured:
`git log --all` returns nothing for either path). `_p` and `_grid` build the
smallest valid row and the smallest valid 101-row quantile grid respectively,
so each test only has to state what it is varying.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

from chd_atlas.tables import PROFILE_QUANTILES, PROFILES
from chd_atlas.validate.profiles import validate_profiles

_DATASET = "GSE999999"


def _p(
    *,
    dataset: str = _DATASET,
    gene: str = "HGNC:11604",
    tissue: str = "Heart",
    stage: str | None = "7wpc",
    median_abundance: float = 10.0,
    unit: str = "tpm",
    q25: float | None = None,
    q75: float | None = None,
    n_samples: int = 5,
) -> dict[str, str]:
    """One `mirrors/profiles/<accession>.tsv` row, as TSV cell strings.

    `stage=None` writes an empty cell -- `profiles.stage` is the one nullable
    column these tests exercise -- everything else defaults to a plain,
    always-valid row so a test only overrides what it cares about.
    """
    return {
        "dataset": dataset,
        "gene": gene,
        "tissue": tissue,
        "stage": "" if stage is None else stage,
        "median_abundance": str(median_abundance),
        "unit": unit,
        "q25": "" if q25 is None else str(q25),
        "q75": "" if q75 is None else str(q75),
        "n_samples": str(n_samples),
    }


def _grid(
    *,
    dataset: str = _DATASET,
    tissue: str = "Heart",
    stage: str = "7wpc",
    unit: str = "tpm",
    values: list[float] | None = None,
    n_genes: int = 20_000,
) -> list[dict[str, str]]:
    """101 `mirrors/profile_quantiles/<accession>.tsv` rows, percentile 0..100.

    Defaults to non-decreasing values (`value == percentile`), so a test only
    has to state what it is varying. `values[i]` is the value placed at
    percentile `i`. The length is asserted rather than silently zipped short
    or padded: a mismatched fixture would cover fewer than the 101 percentiles
    a real quantile grid always has, and every test built on a short grid
    would be exercising a shape production never produces.
    """
    if values is None:
        values = [float(percentile) for percentile in range(101)]
    if len(values) != 101:
        raise ValueError(f"_grid needs exactly 101 values, one per percentile; got {len(values)}")
    return [
        {
            "dataset": dataset,
            "tissue": tissue,
            "stage": stage,
            "percentile": str(percentile),
            "value": str(value),
            "unit": unit,
            "n_genes": str(n_genes),
        }
        for percentile, value in enumerate(values)
    ]


def _write_shard(
    root: Path, directory: str, column_names: tuple[str, ...], rows: list[dict[str, str]]
) -> None:
    """Write `rows` as shard TSVs under `mirrors/<directory>/`.

    One file per distinct `dataset` value, named `<dataset>.tsv` -- matching
    the real convention (`mirrors/profiles/<accession>.tsv`) closely enough
    that `dataset` and the shard's filename agree, as they will for a
    genuinely curated mirror.
    """
    shard_dir = root / "mirrors" / directory
    shard_dir.mkdir(parents=True, exist_ok=True)
    by_dataset: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_dataset.setdefault(row["dataset"], []).append(row)
    header = "\t".join(column_names)
    for dataset, dataset_rows in by_dataset.items():
        body = ["\t".join(row[name] for name in column_names) for row in dataset_rows]
        (shard_dir / f"{dataset}.tsv").write_text(
            "\n".join([header, *body]) + "\n", encoding="utf-8"
        )


def _root_with(
    *,
    profiles: list[dict[str, str]] | None = None,
    quantiles: list[dict[str, str]] | None = None,
) -> Path:
    """A fresh repository root carrying the given profiles/profile_quantiles rows.

    `None` (the default for each) omits that mirror directory entirely --
    standing in for a dataset whose quantiles, or profiles, have not been
    curated yet. That is the ordinary state mid-curation
    `validate_profiles` must not read as "everything is missing"; see its
    module docstring.

    Each call gets its own directory via `tempfile.mkdtemp()` rather than a
    shared `tmp_path` fixture: some tests below build two independent
    repositories in one test body (a fresh root per call), and none of the
    test functions in this file take a fixture argument.
    """
    root = Path(tempfile.mkdtemp())
    if profiles:
        _write_shard(root, "profiles", PROFILES.column_names, profiles)
    if quantiles:
        _write_shard(root, "profile_quantiles", PROFILE_QUANTILES.column_names, quantiles)
    return root


def test_a_clean_matching_repository_reports_nothing() -> None:
    """The smoke test every check below rests on.

    Without it, a rule that fired unconditionally would still make every
    negative case below go green, and the whole suite would agree the
    validator worked.
    """
    root = _root_with(
        profiles=[
            _p(tissue="Heart", stage="7wpc"),
            _p(tissue="Liver", stage="8wpc", gene="HGNC:2020"),
        ],
        quantiles=_grid(tissue="Heart", stage="7wpc") + _grid(tissue="Liver", stage="8wpc"),
    )
    assert validate_profiles(root) == []


def test_prf002_skips_a_null_stage_but_reports_a_real_gap() -> None:
    """A null-stage row can have no quantile partner by construction.

    `profiles.stage` is nullable and `profile_quantiles.stage` is not, so
    firing PRF002 on a null would report a design fact as a curation error,
    once per tissue, on every dataset that has one -- and a curator who learns
    to ignore PRF002 is a curator who misses the real one below.

    The row is not dropped either: a later join is expected to be a left join,
    so it publishes with a null percentile. Both halves are needed.
    """
    issues = validate_profiles(
        _root_with(
            profiles=[
                _p(tissue="Heart", stage=None),
                _p(tissue="Heart", stage="7wpc"),
                _p(tissue="Liver", stage="7wpc"),
            ],
            quantiles=_grid(tissue="Heart", stage="7wpc"),
        )
    )
    codes = sorted(issue.code for issue in issues)
    assert codes == ["PRF002"]
    assert "Liver" in issues[0].location or "Liver" in issues[0].message
    assert "None" not in issues[0].message


def test_prf003_allows_ties_and_refuses_a_decrease() -> None:
    """Ties at the floor are expected on real data.

    Roughly half a bulk library sits at zero, so several dozen consecutive
    breakpoints are legitimately identical. Strict monotonicity would fire on
    every real dataset and the check would be turned off.
    """
    flat = _grid(tissue="Heart", stage="7wpc", values=[0.0] * 46 + [float(i) for i in range(55)])
    assert [i.code for i in validate_profiles(_root_with(quantiles=flat))] == []

    backwards = _grid(tissue="Heart", stage="7wpc", values=[float(i) for i in range(100)] + [3.0])
    assert [i.code for i in validate_profiles(_root_with(quantiles=backwards))] == ["PRF003"]


def test_prf001_and_prf008_catch_the_two_ways_units_can_disagree() -> None:
    """One column cannot merge two quantities if nothing checks it per dataset.

    `allowed` is per cell, so `rpkm` in one organ and `tpm` in another inside
    one dataset passes TBL005 -- which is exactly the merge the column exists
    to prevent.
    """
    mixed = validate_profiles(
        _root_with(
            profiles=[_p(tissue="Heart", unit="rpkm"), _p(tissue="Liver", unit="tpm")],
        )
    )
    assert [i.code for i in mixed] == ["PRF008"]

    crossed = validate_profiles(
        _root_with(
            profiles=[_p(unit="rpkm")],
            quantiles=_grid(unit="tpm"),
        )
    )
    assert "PRF001" in [i.code for i in crossed]


def test_multiple_datasets_flagged_by_prf001_are_reported_in_sorted_order() -> None:
    """PRF001 is sorted by dataset name, and each message's unit list is
    sorted too -- the PRF001 counterpart of the two PRF008 order tests above.

    `set(profile_units) & set(quantile_units)` is a plain `set` built fresh
    by the intersection; its iteration order can follow `PYTHONHASHSEED`
    regardless of file or insertion order, the same risk as PRF002's set
    difference. **Measured**: with three datasets, a dropped `sorted()` around
    this exact loop survived this test in 2 of 11 repeated runs, each a fresh
    interpreter with its own random hash seed -- consistent with a roughly
    1-in-6 chance of a 3-element set coincidentally iterating already
    alphabetical. Five datasets (1-in-120), and all five known units on one of
    them for the message check, measured 0 survivors in 10 further runs.
    """
    datasets = ["Alpha", "Bravo", "Middle", "Yankee", "Zebra"]
    root = _root_with(
        profiles=[_p(dataset=dataset, tissue="Heart", unit="rpkm") for dataset in datasets],
        quantiles=(
            # Every non-Zebra dataset: a plain single-unit mismatch.
            _grid(dataset="Alpha", tissue="Heart", unit="tpm")
            + _grid(dataset="Bravo", tissue="Heart", unit="tpm")
            + _grid(dataset="Middle", tissue="Heart", unit="tpm")
            + _grid(dataset="Yankee", tissue="Heart", unit="tpm")
            # Zebra: all five known units, over five tissues, for the
            # message-content sort check.
            + _grid(dataset="Zebra", tissue="Heart", unit="tpm")
            + _grid(dataset="Zebra", tissue="Liver", unit="nx")
            + _grid(dataset="Zebra", tissue="Brain", unit="cpm")
            + _grid(dataset="Zebra", tissue="Kidney", unit="lfq")
            + _grid(dataset="Zebra", tissue="Lung", unit="rpkm")
        ),
    )

    issues = validate_profiles(root)
    assert [issue.code for issue in issues] == ["PRF001"] * 5
    assert "Alpha" in issues[0].message
    assert "Bravo" in issues[1].message
    assert "Middle" in issues[2].message
    assert "Yankee" in issues[3].message
    assert "Zebra" in issues[4].message
    assert "['cpm', 'lfq', 'nx', 'rpkm', 'tpm']" in issues[4].message


def test_prf008_message_lists_every_offending_unit_in_sorted_order() -> None:
    """PRF008's unit list is sorted, not whatever order a `set` iterates in.

    All five known units land on one dataset, deliberately more than two or
    three: with a smaller set, an unsorted `set` has a non-negligible chance
    of iterating in alphabetical order by pure luck on any given interpreter,
    which would let a dropped `sorted()` pass undetected on an unlucky run.
    Five elements have a 1-in-120 chance of a coincidentally-sorted iteration,
    a small enough risk to trust the literal below.
    """
    profiles = [
        _p(tissue=tissue, unit=unit)
        for tissue, unit in zip(
            ["Heart", "Liver", "Brain", "Kidney", "Lung"],
            ["tpm", "rpkm", "nx", "cpm", "lfq"],
            strict=True,
        )
    ]
    issues = validate_profiles(_root_with(profiles=profiles))
    assert [issue.code for issue in issues] == ["PRF008"]
    assert "['cpm', 'lfq', 'nx', 'rpkm', 'tpm']" in issues[0].message


def test_multiple_missing_triples_are_reported_in_sorted_order() -> None:
    """PRF002 issues are sorted by (dataset, tissue, stage), not set/hash order.

    `profile_triples - quantile_triples` is a freshly built `set`; iterating
    it unsorted would let the report's order follow `PYTHONHASHSEED`, so two
    validate runs over one identical repository could disagree on which issue
    comes first. **Measured**: with three tissues, a dropped `sorted()` around
    this exact loop survived this test in 1 of 11 repeated runs (each a fresh
    interpreter, its own random hash seed) -- consistent with a roughly 1-in-6
    chance of a 3-element set coincidentally iterating already alphabetical.
    Five tissues (1-in-120 by the same reasoning) measured 0 survivors in 10
    further runs.
    """
    tissues = ["Alpha", "Bravo", "Middle", "Yankee", "Zebra"]
    profiles = [_p(tissue=tissue, stage="7wpc") for tissue in tissues]
    # A grid for a tissue none of the profiles rows use, so profile_quantiles
    # is present (PRF002 is not skipped wholesale) without satisfying any of
    # the five triples above.
    quantiles = _grid(tissue="Nowhere", stage="7wpc")
    issues = validate_profiles(_root_with(profiles=profiles, quantiles=quantiles))

    assert [issue.code for issue in issues] == ["PRF002"] * 5
    assert "Alpha" in issues[0].message
    assert "Bravo" in issues[1].message
    assert "Middle" in issues[2].message
    assert "Yankee" in issues[3].message
    assert "Zebra" in issues[4].message


def test_multiple_datasets_flagged_by_prf008_are_reported_in_sorted_order() -> None:
    """PRF008 is sorted by dataset name, not row/insertion order.

    `profile_units` is a `dict` keyed by dataset; a `dict` iterates in
    insertion order, which here follows file row order -- not dataset name.
    Three datasets are written to *one* shard in Zebra/Alpha/Middle row order,
    deliberately not alphabetical, so a dropped `sorted()` around the dataset
    loop would report them in that same row order instead of alphabetically.
    Unlike the `set`-based checks above, a `dict`'s order does not depend on
    `PYTHONHASHSEED`, so this is the one place that risk cannot be leaned on
    to catch a missing sort -- the row order has to disagree with the
    alphabetical one directly, by construction of the fixture.
    """
    root = _root_with()
    rows = [
        _p(dataset="Zebra", unit="tpm"),
        _p(dataset="Zebra", unit="rpkm"),
        _p(dataset="Alpha", unit="tpm"),
        _p(dataset="Alpha", unit="rpkm"),
        _p(dataset="Middle", unit="tpm"),
        _p(dataset="Middle", unit="rpkm"),
    ]
    shard_dir = root / "mirrors" / "profiles"
    shard_dir.mkdir(parents=True, exist_ok=True)
    header = "\t".join(PROFILES.column_names)
    body = ["\t".join(row[name] for name in PROFILES.column_names) for row in rows]
    (shard_dir / "combined.tsv").write_text("\n".join([header, *body]) + "\n", encoding="utf-8")

    issues = validate_profiles(root)
    assert [issue.code for issue in issues] == ["PRF008", "PRF008", "PRF008"]
    assert "Alpha" in issues[0].message
    assert "Middle" in issues[1].message
    assert "Zebra" in issues[2].message


def test_prf003_multiple_triples_are_reported_in_sorted_order() -> None:
    """PRF003 issues are sorted by (dataset, tissue, stage), not dict/row order.

    `breakpoints` is a `dict` keyed by triple; three tissues, each with its
    own single-pair inversion, are written in Zebra/Alpha/Middle row order (one
    dataset, so `_write_shard` preserves that order in the one file it
    writes). A dropped `sorted()` around the triple loop would report these in
    file order instead of alphabetically -- the `dict`-order counterpart of
    `test_multiple_missing_triples_are_reported_in_sorted_order`, which tests
    the `set`-based PRF002 collection instead.
    """
    rows: list[dict[str, str]] = []
    for tissue in ("Zebra", "Alpha", "Middle"):
        rows.extend(
            _grid(tissue=tissue, stage="7wpc", values=[float(i) for i in range(100)] + [3.0])
        )

    issues = validate_profiles(_root_with(quantiles=rows))
    assert [issue.code for issue in issues] == ["PRF003", "PRF003", "PRF003"]
    assert "Alpha" in issues[0].message
    assert "Middle" in issues[1].message
    assert "Zebra" in issues[2].message


def test_prf003_sorts_by_percentile_before_checking_monotonicity() -> None:
    """A defensive re-sort, not trust in file/row order.

    `read_table` returns rows in whatever order the TSV has them; nothing
    guarantees a `profile_quantiles` shard arrives sorted by `percentile`
    before `validate_profiles` sees it (`validate_sort_order` is a separate
    check, and this function is exercised standalone throughout this file).
    Two adjacent rows swapped in the file, with values that are perfectly
    non-decreasing once put back in percentile order, must not read as a
    decrease.
    """
    rows = _grid()
    rows[1], rows[2] = rows[2], rows[1]  # percentiles 1 and 2 swap file position
    assert validate_profiles(_root_with(quantiles=rows)) == []


def test_an_absent_or_unreadable_mirror_reports_nothing_and_does_not_raise() -> None:
    """Three distinct absences, none of which is this validator's to report.

    No `profiles` or `profile_quantiles` mirror exists before either is
    curated. An unreadable shard is TBL000's to report and a shard missing a
    column is TBL001's, both ERRORs against the same path, so the build
    refuses either way and a second report here would only say it twice.
    Raising instead would abort validation of every other shard in the
    repository.
    """
    root = _root_with()
    assert validate_profiles(root) == []

    (root / "mirrors" / "profiles").mkdir(parents=True, exist_ok=True)
    (root / "mirrors" / "profiles" / "bad.tsv").write_bytes(b"\xff\xfe not utf-8")
    assert validate_profiles(root) == []

    (root / "mirrors" / "profile_quantiles").mkdir(parents=True, exist_ok=True)
    (root / "mirrors" / "profile_quantiles" / "bad.tsv").write_text(
        f"dataset\ttissue\n{_DATASET}\tHeart\n", encoding="utf-8"
    )
    assert validate_profiles(root) == []


def test_an_unreadable_shard_does_not_stop_the_others() -> None:
    """One corrupt shard costs one TBL000 elsewhere, not every other dataset's
    checks.

    `_shard_rows` skips a bad shard with `continue` inside the loop over
    `mirror_paths`, never an early `return` -- the same "one bad file, one
    issue" discipline as `validate_burden`.
    """
    root = _root_with(
        profiles=[_p(dataset="GSE_GOOD", tissue="Liver", stage="7wpc")],
        quantiles=_grid(dataset="GSE_GOOD", tissue="Heart", stage="7wpc"),
    )
    (root / "mirrors" / "profiles" / "GSE_BAD.tsv").write_bytes(b"\xff\xfe not utf-8")

    issues = validate_profiles(root)
    assert [issue.code for issue in issues] == ["PRF002"]
    assert "GSE_GOOD" in issues[0].message
