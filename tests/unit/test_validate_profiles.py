# tests/unit/test_validate_profiles.py
"""The rules that stop a gene's percentile being read against another
dataset's ruler.

Every fixture here is synthetic rather than transcribed from a real study,
unlike `test_validate_burden.py`, which builds on one real row -- these checks
are about shapes two tables can disagree on (a unit, a cell, a monotone
breakpoint run), not about any one dataset's actual numbers, so a synthetic
row states what each test varies without pulling in E-MTAB-6814's own
1.2M-row mirror. `_p` and `_grid` build the smallest valid row and the
smallest valid 101-row quantile grid respectively, so each test only has to
state what it is varying.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

from chd_atlas.issues import Severity
from chd_atlas.models.dataset import Dataset
from chd_atlas.models.phases import CardiacPhase, CardiacPhaseFile, EndBasis
from chd_atlas.tables import PROFILE_QUANTILES, PROFILES
from chd_atlas.validate.profiles import validate_profile_references, validate_profiles

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


def _profile_dataset(**overrides: object) -> Dataset:
    """A minimal, valid profile-design `Dataset`, for `validate_profile_references`.

    Defaults to `id=_DATASET` -- matching `_p`/`_grid`'s own default dataset --
    so a fixture combining this with `_root_with` addresses the same dataset
    unless a test overrides both. `tests/unit/test_dataset_model.py::
    _profile_dataset` is the model-level precedent this mirrors; its `id`
    defaults to a real ArrayExpress accession that has no reason to agree
    with this file's `_DATASET`, so it is not reused directly.
    """
    base: dict[str, object] = {
        "id": _DATASET,
        "archive": "geo",
        "technology": "bulk_rnaseq",
        "design": "profile",
        "tissue": "Heart",
        "developmental_stage": "embryonic",
        "organism": "NCBITaxon:9606",
        "n_samples": 12,
        "licence": "CC BY 4.0",
        "contrasts": [],
        "cardiac_tissues": ["Heart"],
        "detection_floor": 1.0,
        "floor_source": "test fixture",
        "quantile_estimator": "linear",
        "stages": [{"token": "7wpc", "wpc": 7.0, "order": 1}],
    }
    base.update(overrides)
    return Dataset(**base)


def _phase(phase_id: str, start_wpc: float, end_wpc: float) -> CardiacPhase:
    """A validly-constructed `CardiacPhase` for interval-logic tests.

    The Carnegie-stage/HsapDv/`go_id`/`end_basis` values are placeholders,
    not a claim about embryology -- only `id`/`start_wpc`/`end_wpc` vary
    across the tests that use this factory, matching `test_models_phases.py`'s
    own `_phase()`.
    """
    return CardiacPhase(
        id=phase_id,
        go_id="GO:0000001",
        label=phase_id,
        start_wpc=start_wpc,
        end_wpc=end_wpc,
        start_carnegie_stage="CS1",
        end_carnegie_stage="CS2",
        start_hsapdv_id="HsapDv:0000001",
        end_hsapdv_id="HsapDv:0000002",
        end_basis=EndBasis.STATED,
    )


def _open_phase(phase_id: str, start_wpc: float) -> CardiacPhase:
    """A validly-constructed `CardiacPhase` with a real start and NO stated end.

    The `heart_looping` shape, matching `test_models_phases.py`'s and
    `test_build_profiles.py`'s own `_open_phase()`.
    """
    return CardiacPhase(
        id=phase_id,
        go_id="GO:0000001",
        label=phase_id,
        start_wpc=start_wpc,
        end_wpc=None,
        start_carnegie_stage="CS1",
        end_carnegie_stage=None,
        start_hsapdv_id="HsapDv:0000001",
        end_hsapdv_id=None,
        end_basis=EndBasis.NOT_STATED,
    )


def _phase_file(*phases: CardiacPhase) -> CardiacPhaseFile:
    """A minimal, valid `CardiacPhaseFile` wrapping the given phases."""
    return CardiacPhaseFile(
        attributed_to="Test Reference Atlas of Human Heart Development",
        citation="Test et al., 2026",
        phases=list(phases),
    )


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
    difference.

    **The outer, dataset-name sort measured clean.** A fast probe
    constructing `{"Alpha", "Bravo", "Middle", "Yankee", "Zebra"}` directly
    under 300 explicit `PYTHONHASHSEED` values (not however many random ones
    a CI run happens to sample) measured 0/300 coincidentally-sorted
    iterations, and a slow pytest-based run of the equivalent mutation
    (dropping this exact `sorted()`) measured 20/20 killed across two
    fresh-process batches.

    **The inner, message-content sort could not be fixture-widened the same
    way, because the real abundance-unit vocabulary has exactly five
    members** (`tpm`, `nx`, `cpm`, `lfq`, `rpkm` -- `tables.py`'s
    `_ABUNDANCE_UNITS`), and this docstring's previous claim that five
    elements were "1-in-120," a small enough risk to trust, was never
    checked past 10 runs and does not hold: the same probe method measured
    4/500 (0.8%) for that exact five-member set. There is no sixth real unit
    to widen into. `_prf001_issues` treats `unit` as an opaque string for
    this comparison and never checks it against `_ABUNDANCE_UNITS` (that is
    TBL004's job, a different validator, exercised in `test_tables.py`), so
    Zebra's eight quantile rows below use eight *synthetic* unit tokens
    instead of real ones. The synthetic set measured 0/500 in the same probe
    and 20/20 killed in the slow pytest-based check.
    """
    datasets = ["Alpha", "Bravo", "Middle", "Yankee", "Zebra"]
    synthetic_units = ["alpha", "bravo", "charlie", "delta", "echo", "foxtrot", "golf", "hotel"]
    organs = ["Heart", "Liver", "Brain", "Kidney", "Lung", "Spleen", "Pancreas", "Thymus"]
    zebra_quantiles: list[dict[str, str]] = []
    for organ, unit in zip(organs, synthetic_units, strict=True):
        zebra_quantiles += _grid(dataset="Zebra", tissue=organ, unit=unit)

    root = _root_with(
        profiles=[_p(dataset=dataset, tissue="Heart", unit="rpkm") for dataset in datasets],
        quantiles=(
            # Every non-Zebra dataset: a plain single-unit mismatch.
            _grid(dataset="Alpha", tissue="Heart", unit="tpm")
            + _grid(dataset="Bravo", tissue="Heart", unit="tpm")
            + _grid(dataset="Middle", tissue="Heart", unit="tpm")
            + _grid(dataset="Yankee", tissue="Heart", unit="tpm")
            # Zebra: eight distinct (synthetic) units, for the
            # message-content sort check.
            + zebra_quantiles
        ),
    )

    issues = validate_profiles(root)
    assert [issue.code for issue in issues] == ["PRF001"] * 5
    assert "Alpha" in issues[0].message
    assert "Bravo" in issues[1].message
    assert "Middle" in issues[2].message
    assert "Yankee" in issues[3].message
    assert "Zebra" in issues[4].message
    assert (
        "['alpha', 'bravo', 'charlie', 'delta', 'echo', 'foxtrot', 'golf', 'hotel']"
        in issues[4].message
    )


def test_prf008_message_lists_every_offending_unit_in_sorted_order() -> None:
    """PRF008's unit list is sorted, not whatever order a `set` iterates in.

    **Measured, and the real vocabulary is too small to fixture safely.** The
    abundance-unit vocabulary has exactly five members (`tpm`, `nx`, `cpm`,
    `lfq`, `rpkm`). This docstring previously claimed five elements have "a
    1-in-120 chance of a coincidentally-sorted iteration, a small enough risk
    to trust" -- reasoned from permutation counting, never measured. A fast
    probe constructing that exact five-member `set` under 500 explicit
    `PYTHONHASHSEED` values measured 4/500 (0.8%) coincidentally-sorted
    iterations: real, and about 100x the claimed rate. There is no sixth real
    unit to widen into. `_prf008_issues` treats `unit` as an opaque string --
    membership in the real vocabulary is TBL004's check, not this one's -- so
    the eight rows below use eight synthetic unit tokens instead, matching
    `test_multiple_datasets_flagged_by_prf001_are_reported_in_sorted_order`'s
    fix for the same constraint. The synthetic set measured 0/500 in the same
    probe and 20/20 killed in a slow pytest-based fresh-process check. See
    `test_multiple_missing_triples_are_reported_in_sorted_order` for the
    fuller correction of the permutation-counting reasoning itself.
    """
    organs = ["Heart", "Liver", "Brain", "Kidney", "Lung", "Spleen", "Pancreas", "Thymus"]
    synthetic_units = ["alpha", "bravo", "charlie", "delta", "echo", "foxtrot", "golf", "hotel"]
    profiles = [
        _p(tissue=tissue, unit=unit) for tissue, unit in zip(organs, synthetic_units, strict=True)
    ]
    issues = validate_profiles(_root_with(profiles=profiles))
    assert [issue.code for issue in issues] == ["PRF008"]
    assert (
        "['alpha', 'bravo', 'charlie', 'delta', 'echo', 'foxtrot', 'golf', 'hotel']"
        in issues[0].message
    )


def test_prf001_message_also_sorts_the_profiles_side_unit_list() -> None:
    """PRF001's message has two independent sorted lists -- `p_units` and
    `q_units` -- and only `q_units` had a fixture that could tell a dropped
    `sorted()` apart from doing nothing.

    **A real, pre-existing gap, found by checking rather than assuming.**
    Every PRF001 fixture in this file gives the *profiles* side of a
    mismatched dataset exactly one unit (`p_units` is always a 1-element
    set), so `sorted(p_units)` and `list(p_units)` render identically
    regardless of whether the sort actually runs. Measured directly: applying
    that exact mutation and re-running the full file reported all 27 tests
    passing, unchanged -- the mutant survived every test that existed before
    this one.

    A dataset whose own profiles rows carry more than one unit also trips
    PRF008 (a separate, unconditional check on the profiles side alone), so
    this cannot share the fixtures above without conflating two findings --
    hence a dedicated dataset and an assertion that expects both codes.
    """
    synthetic_units = ["alpha", "bravo", "charlie", "delta", "echo", "foxtrot", "golf", "hotel"]
    profiles = [_p(tissue="Heart", unit=unit) for unit in synthetic_units]
    quantiles = _grid(tissue="Heart", unit="india")

    issues = validate_profiles(_root_with(profiles=profiles, quantiles=quantiles))
    assert sorted(issue.code for issue in issues) == ["PRF001", "PRF008"]
    prf001 = next(issue for issue in issues if issue.code == "PRF001")
    assert (
        "['alpha', 'bravo', 'charlie', 'delta', 'echo', 'foxtrot', 'golf', 'hotel']"
        in prf001.message
    )


def test_multiple_missing_triples_are_reported_in_sorted_order() -> None:
    """PRF002 issues are sorted by (dataset, tissue, stage), not set/hash order.

    `profile_triples - quantile_triples` is a freshly built `set`; iterating
    it unsorted would let the report's order follow `PYTHONHASHSEED`, so two
    validate runs over one identical repository could disagree on which issue
    comes first.

    **This docstring previously reasoned from permutation counting -- "five
    tissues (1-in-120 by the same reasoning) measured 0 survivors in 10
    further runs" -- and that reasoning does not hold. It described a
    fixture that was never independently checked past those 10 samples.**
    Measured directly instead: constructing the actual object this function
    sorts -- a `set` of `("GSE999999", tissue, "7wpc")` triples, not a bare
    tissue string -- under 500 explicit `PYTHONHASHSEED` values (not 10
    random ones), the original five words (`Alpha, Bravo, Middle, Yankee,
    Zebra`) measured **9/500 (1.8%)** coincidentally-sorted iterations once
    wrapped in that triple. That is consistent with the original 10-run
    measurement's 0 survivors -- at a true rate of 1.8%, seeing 0 survivors
    in 10 runs happens roughly 84% of the time -- but the 0 was never proof
    of safety, only a plausible outcome of an underpowered sample; a
    slow pytest-based re-run of the equivalent mutation measured the real
    exposure directly: 1/10 fresh-process survivors.

    The same five words as *bare strings*, with no triple, had separately
    measured 0/300 in an earlier probe -- so the coincidence is a property of
    the **container shape** as much as the string content: wrapping a value
    in a tuple changes which hash-table slot it lands in relative to the
    others. Permutation counting (1-in-N! for N elements) never modelled
    either shape; it counts orderings, not hash-table layouts.

    Eight qualitatively distinct words (`Alpha`..`Hotel`, the NATO alphabet)
    measured 0/500 in the same tuple-wrapped probe, and 20/20 killed in a
    slow pytest-based fresh-process re-run. That is the fixture below.
    """
    tissues = ["Alpha", "Bravo", "Charlie", "Delta", "Echo", "Foxtrot", "Golf", "Hotel"]
    profiles = [_p(tissue=tissue, stage="7wpc") for tissue in tissues]
    # A grid for a tissue none of the profiles rows use, so profile_quantiles
    # is present (PRF002 is not skipped wholesale) without satisfying any of
    # the eight triples above.
    quantiles = _grid(tissue="Nowhere", stage="7wpc")
    issues = validate_profiles(_root_with(profiles=profiles, quantiles=quantiles))

    assert [issue.code for issue in issues] == ["PRF002"] * len(tissues)
    for index, tissue in enumerate(tissues):
        assert tissue in issues[index].message


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


# ---------------------------------------------------------------------------
# validate_profile_references -- PRF004/005/006/007/009
# ---------------------------------------------------------------------------


def test_a_clean_repository_reports_nothing_via_references() -> None:
    """The references-function smoke test every check below rests on.

    Mirrors `test_a_clean_matching_repository_reports_nothing` above: without
    it, a rule that fired unconditionally would still make every negative
    case below look like it passed for the right reason.
    """
    root = _root_with(
        profiles=[_p(gene="HGNC:11604", tissue="Heart", stage="7wpc")],
        quantiles=_grid(tissue="Heart", stage="7wpc"),
    )
    issues = validate_profile_references(
        root,
        datasets=(_profile_dataset(),),
        known_genes={"HGNC:11604"},
        published_genes={"HGNC:11604"},
        phases=_phase_file(_phase("septation", 5.0, 8.0)),
    )
    assert issues == []


def test_prf004_names_a_declared_cardiac_tissue_missing_from_profiles() -> None:
    """D40: a tau computation needs the declared token to have real data behind it.

    `cardiac_tissues` is checked as a declared set, never matched against the
    literal string "heart" -- so a labelling drift between the curated
    dataset record and the converted profiles rows (a different tissue
    string, or no row at all) must be caught, not silently computed against
    zero rows.
    """
    issues = validate_profile_references(
        _root_with(profiles=[_p(tissue="Heart", stage="7wpc")]),
        datasets=(_profile_dataset(cardiac_tissues=["Heart", "Atrium"]),),
        known_genes={"HGNC:11604"},
        published_genes=set(),
        phases=None,
    )
    assert [i.code for i in issues] == ["PRF004"]
    assert issues[0].severity is Severity.ERROR
    assert _DATASET in issues[0].message
    assert "['Atrium']" in issues[0].message


def test_prf004_fires_even_when_the_dataset_has_no_profiles_rows_at_all() -> None:
    """A curated Dataset record with nothing behind it yet is still an error.

    Unlike an absent `profile_quantiles` mirror (Task 7's PRF000/PRF010),
    which the runner treats as the ordinary mid-curation state because
    profiles genuinely precedes it, here the Dataset record is the newer
    artefact and its declared `cardiac_tissues` token already asserts a claim
    the mirror has not backed up. Also confirms PRF005 does not spuriously
    fire in the mirror-image empty case: zero observed stage tokens is zero
    undeclared ones.
    """
    issues = validate_profile_references(
        _root_with(),
        datasets=(_profile_dataset(),),
        known_genes=None,
        published_genes=set(),
        phases=None,
    )
    assert [i.code for i in issues] == ["PRF004"]
    assert _DATASET in issues[0].location


def test_prf005_names_a_stage_token_the_dataset_record_never_declared() -> None:
    """The reverse direction from PRF004: the mirror has data the record does
    not vouch for, rather than a declaration the mirror does not back up.
    """
    issues = validate_profile_references(
        _root_with(profiles=[_p(tissue="Heart", stage="99wpc")]),
        datasets=(_profile_dataset(),),  # declares stages=[{"token": "7wpc", ...}]
        known_genes={"HGNC:11604"},
        published_genes=set(),
        phases=None,
    )
    assert [i.code for i in issues] == ["PRF005"]
    assert issues[0].severity is Severity.ERROR
    assert "['99wpc']" in issues[0].message


def test_prf004_and_prf005_are_reported_in_sorted_order_by_dataset_and_within_message() -> None:
    """Two independent dropped-`sorted()` mutants per code, guarded at once --
    sizes measured directly, not carried over from Task 5's docstring unchecked.

    `_prf004_issues` and `_prf005_issues` each sort twice: once across
    datasets, and once within one dataset's missing-token message. The first
    is a list sort over a tuple the function re-sorts itself rather than
    trusting caller order (not hash-seed dependent), so five datasets are
    constructed here in a scrambled (Zebra, Alfa, Yankee, Midl, Brav)
    construction order -- not the alphabetical order the report must come
    back in.

    **The message-list sorts needed more than five elements, and a different
    string family, measured.** Task 5 documented five as enough (1-in-120 by
    permutation counting) and measured 0/10 survivors for its own fixtures.
    This test originally used five elements too (`Atrium, Heart, OFT,
    SinusVenosus, Ventricle` for the tissue set; `10wpc`..`14wpc` for the
    stage set) and one slow pytest-based fresh-process run measured 10/10
    killed for both -- which turned out to be a lucky sample, not a safe
    fixture. A follow-up probe under 500 explicit PYTHONHASHSEED values (not
    10 random ones) found a real, nonzero coincidence rate behind that one
    clean run: 4/500 for the original tissue set, 10/500 for the original
    stage set (both differed only in a short numeric or near-numeric suffix,
    the same shape that made PRF007/PRF009's gene-id fixtures measurably
    flaky below). Eight elements, and words that do not share a short
    numeric suffix, measured 0/500 for both replacement sets in the same
    probe.
    """
    ids_alphabetical = ["E-ALFA-1", "E-BRAV-1", "E-MIDL-1", "E-YANK-1", "E-ZEBR-1"]
    cardiac_structures = [
        "Atrium",
        "Bulbus",
        "Conotruncus",
        "Ductus",
        "Endocardium",
        "Foramen",
        "GreatVessel",
        "Horn",
    ]
    developmental_stages = [
        "cleavage",
        "gastrula",
        "neurula",
        "somite",
        "pharyngula",
        "limbbud",
        "fetal",
        "neonatal",
    ]
    datasets = (
        _profile_dataset(id="E-ZEBR-1", archive="arrayexpress", cardiac_tissues=cardiac_structures),
        _profile_dataset(id="E-ALFA-1", archive="arrayexpress", cardiac_tissues=["Heart"]),
        _profile_dataset(id="E-YANK-1", archive="arrayexpress", cardiac_tissues=["Heart"]),
        _profile_dataset(id="E-MIDL-1", archive="arrayexpress", cardiac_tissues=["Heart"]),
        _profile_dataset(id="E-BRAV-1", archive="arrayexpress", cardiac_tissues=["Heart"]),
    )
    rows = [
        _p(dataset=dataset_id, tissue="Liver", stage="99wpc")
        for dataset_id in ["E-ALFA-1", "E-BRAV-1", "E-MIDL-1", "E-YANK-1"]
    ]
    rows += [_p(dataset="E-ZEBR-1", tissue="Liver", stage=stage) for stage in developmental_stages]

    issues = validate_profile_references(
        _root_with(profiles=rows),
        datasets=datasets,
        known_genes=None,
        published_genes=set(),
        phases=None,
    )
    prf004 = [i for i in issues if i.code == "PRF004"]
    prf005 = [i for i in issues if i.code == "PRF005"]
    assert len(prf004) == 5
    assert len(prf005) == 5
    for index, dataset_id in enumerate(ids_alphabetical):
        assert dataset_id in prf004[index].message
        assert dataset_id in prf005[index].message
    assert (
        "['Atrium', 'Bulbus', 'Conotruncus', 'Ductus', 'Endocardium', "
        "'Foramen', 'GreatVessel', 'Horn']" in prf004[4].message
    )
    assert (
        "['cleavage', 'fetal', 'gastrula', 'limbbud', 'neonatal', "
        "'neurula', 'pharyngula', 'somite']" in prf005[4].message
    )


def test_prf011_reports_two_stages_claiming_one_position() -> None:
    """Duplicate `order` within one dataset.

    Two stages at position 3 hand the sequence back to the sort's tie-break --
    the token -- which is exactly the alphabetical ordering this field was
    added to remove, reappearing silently inside one dataset.

    **Numbered PRF011, not PRF010.** `validate/runner.py` already issues
    PRF010 for a different failure (profiles rows present, no
    `profile_quantiles` data), where it is PRF000's partner error. One code
    naming two unrelated failures would let a lone skip warning pair with the
    wrong error, which is the whole point of section 4.41.

    **The colliding pair is listed `d` before `b`**, and that is the whole
    reason the message is asserted as a sequence rather than as two
    independent `in` checks. This check scans `sorted(stages, key=(order,
    token))` precisely so that which of a pair is named first does not depend
    on which record a curator happened to type first -- `Stage`'s own
    docstring says the list's order is deliberately not load-bearing -- and
    with the pair listed in already-sorted order the sort is unobservable:
    dropping it left both `in` checks green. A list sort over a tuple, so it
    is comparison-based and owes no seed probe (section 4.42); one
    deliberately unsorted fixture is the whole guard.

    **`wpc` still ascends in the order the pair is listed**, and that is not
    incidental. `_prf012_issues` sorts on `order` alone, and Python's sort is
    stable, so a tie inherits the curated list's order -- a fixture that
    listed the later `wpc` first would provoke a real PRF012 as well and turn
    a test about one check's scan order into a test about two.
    """
    issues = validate_profile_references(
        _root_with(profiles=[_p(tissue="Heart", stage="a")]),
        datasets=(
            _profile_dataset(
                stages=[
                    {"token": "a", "wpc": 1.0, "order": 1},
                    {"token": "d", "wpc": 2.0, "order": 3},
                    {"token": "b", "wpc": 3.0, "order": 3},
                ]
            ),
        ),
        known_genes={"HGNC:11604"},
        published_genes=set(),
        phases=None,
    )
    assert [issue.code for issue in issues] == ["PRF011"]
    assert issues[0].severity is Severity.ERROR
    assert "3" in issues[0].message
    assert "stages 'b' and 'd'" in issues[0].message


def test_prf012_reports_an_order_that_contradicts_wpc() -> None:
    """`order` and `wpc` disagree about which of two stages is earlier.

    The failure `Stage.order` introduces by existing: each record is
    individually valid -- an integer and a positive float -- so pydantic
    cannot see it, and a stage transcribed into the wrong slot would publish
    a trajectory with two of its points swapped.
    """
    issues = validate_profile_references(
        _root_with(profiles=[_p(tissue="Heart", stage="early")]),
        datasets=(
            _profile_dataset(
                stages=[
                    {"token": "early", "wpc": 4.0, "order": 2},
                    {"token": "late", "wpc": 9.0, "order": 1},
                ]
            ),
        ),
        known_genes={"HGNC:11604"},
        published_genes=set(),
        phases=None,
    )
    assert [issue.code for issue in issues] == ["PRF012"]
    assert issues[0].severity is Severity.ERROR
    assert "'late'" in issues[0].message
    assert "'early'" in issues[0].message


def test_prf012_refuses_two_stages_that_claim_the_same_wpc() -> None:
    """The boundary is `>=`, not `>`, and that is a decision rather than an
    off-by-one.

    Two tokens at one `wpc` are two names for one point in time, so `order`
    is claiming a sequence its own evidence does not support -- and `wpc` is
    the key `assign_phase` joins on, so both tokens would resolve to
    identical phases while the axis drew them apart. The likelier cause is a
    transcription slip that duplicated a number.

    Without this case the `>` mutant survives: every other PRF012 fixture
    uses strictly decreasing values, which both operators catch.
    """
    issues = validate_profile_references(
        _root_with(profiles=[_p(tissue="Heart", stage="a")]),
        datasets=(
            _profile_dataset(
                stages=[
                    {"token": "a", "wpc": 4.0, "order": 1},
                    {"token": "b", "wpc": 4.0, "order": 2},
                ]
            ),
        ),
        known_genes={"HGNC:11604"},
        published_genes=set(),
        phases=None,
    )
    assert [issue.code for issue in issues] == ["PRF012"]


def test_prf012_is_silent_when_only_post_natal_stages_lack_wpc() -> None:
    """A null `wpc` is never evidence of a contradiction.

    Every real profile dataset ends in post-natal stages that carry no `wpc`
    at all, so a check that compared against the null would report every
    dataset this atlas will ever curate. It is also why the post-natal block
    of `curation/datasets/E-MTAB-6814.yaml` is pinned by a literal in
    `tests/test_repository_validates.py` instead: PRF012 cannot reach it.
    """
    issues = validate_profile_references(
        _root_with(profiles=[_p(tissue="Heart", stage="19wpc")]),
        datasets=(
            _profile_dataset(
                stages=[
                    {"token": "19wpc", "wpc": 19.0, "order": 1},
                    {"token": "neonate", "wpc": None, "order": 2},
                    {"token": "elderly", "wpc": None, "order": 3},
                ]
            ),
        ),
        known_genes={"HGNC:11604"},
        published_genes=set(),
        phases=None,
    )
    assert issues == []


def test_prf011_and_prf012_are_reported_in_dataset_id_order() -> None:
    """Both checks re-sort `datasets` rather than trusting caller order.

    `corpus.datasets` follows `curation/datasets/`'s directory-listing order,
    not id, and a direct caller -- this test, and `validate_repository`'s own
    single call -- reads what these functions return, not a report that has
    been re-sorted afterwards. Constructed Zebra-before-Alfa so a dropped
    `sorted(datasets, ...)` in either function reverses the pair.

    A list sort over a tuple, so it is comparison-based and carries no
    PYTHONHASHSEED dependence (section 4.42) -- two datasets are enough, and
    no seed probe is owed here.

    PRF004 also fires twice (neither dataset has a profiles row behind its
    declared tissue) and is filtered out: this test is about two orderings,
    not about how many codes one fixture can provoke.
    """
    stages = [
        {"token": "a", "wpc": 1.0, "order": 1},
        {"token": "b", "wpc": 5.0, "order": 2},
        {"token": "c", "wpc": 3.0, "order": 3},
        {"token": "d", "wpc": 4.0, "order": 3},
    ]
    issues = validate_profile_references(
        _root_with(),
        datasets=(
            _profile_dataset(id="E-ZEBR-1", archive="arrayexpress", stages=stages),
            _profile_dataset(id="E-ALFA-1", archive="arrayexpress", stages=stages),
        ),
        known_genes=None,
        published_genes=set(),
        phases=None,
    )

    for code in ("PRF011", "PRF012"):
        reported = [issue for issue in issues if issue.code == code]
        assert len(reported) == 2, code
        assert "E-ALFA-1" in reported[0].message, code
        assert "E-ZEBR-1" in reported[1].message, code
        # The location names the curated file a curator would open, not the
        # mirror shard PRF004/PRF005 blame.
        assert reported[0].location.endswith("curation/datasets/E-ALFA-1.yaml"), code


def test_prf006_is_silent_on_contiguous_or_absent_phases() -> None:
    """No interior gap, and no curated vocabulary at all, both report nothing.

    The region past the last phase (here, wpc >= 8) is the ordinary case --
    an embryology window that legitimately stops before the post-natal
    stages -- and must never fire: a check that fires on every correct
    vocabulary is a check a curator learns to ignore.
    """
    touching = _phase_file(_phase("early", 3.0, 5.0), _phase("late", 5.0, 8.0))
    issues = validate_profile_references(
        _root_with(), datasets=(), known_genes=None, published_genes=set(), phases=touching
    )
    assert issues == []

    issues_no_vocabulary = validate_profile_references(
        _root_with(), datasets=(), known_genes=None, published_genes=set(), phases=None
    )
    assert issues_no_vocabulary == []


def test_prf006_reports_one_issue_per_interior_gap_in_wpc_order() -> None:
    """Phases declared out of wpc order in the file, to prove the re-sort is load-bearing.

    Declared D, A, C, B: adjacent-in-list pairs (D,A), (A,C), (C,B) are not
    adjacent in developmental time, so a dropped `sorted()` here does not
    merely reorder the report -- it compares the wrong phases, missing real
    gaps and (D ends at 25, A starts at 1) inventing none where the raw
    pairing happens not to overlap either. Four phases sorted by wpc are
    A(1-3), B(5-8), C(10-15), D(20-25): three interior gaps, and neither the
    region before A (wpc < 1) nor after D (wpc >= 25) may be reported.
    """
    phases = _phase_file(
        _phase("D", 20.0, 25.0),
        _phase("A", 1.0, 3.0),
        _phase("C", 10.0, 15.0),
        _phase("B", 5.0, 8.0),
    )
    issues = validate_profile_references(
        _root_with(), datasets=(), known_genes=None, published_genes=set(), phases=phases
    )
    assert [i.code for i in issues] == ["PRF006"] * 3
    assert all(i.severity is Severity.WARNING for i in issues)
    assert "'A'" in issues[0].message and "'B'" in issues[0].message
    assert "3.0" in issues[0].message and "5.0" in issues[0].message
    assert "'B'" in issues[1].message and "'C'" in issues[1].message
    assert "8.0" in issues[1].message and "10.0" in issues[1].message
    assert "'C'" in issues[2].message and "'D'" in issues[2].message
    assert "15.0" in issues[2].message and "20.0" in issues[2].message


def test_prf006_does_not_report_a_gap_a_wider_overlapping_phase_already_covers() -> None:
    """The unique killer of a pairwise-adjacent PRF006 that never learned to merge.

    Phases may overlap by design (`models/phases.py`'s module docstring), and
    once they do, "sort by start, compare each phase only to its immediate
    successor" stops being a correct gap test. `A=[1,10)` fully covers
    `B=[2,3)` and `C=[8,9)`; sorted by start this is A, B, C, and a check that
    compares B directly to C would see `C.start (8) > B.end (3)` and report a
    false gap in `[3,8)` -- even though A already covers every wpc in that
    range. The real vocabulary this test guards is exactly this shape: `heart
    _looping` spans almost the entire curated window and several shorter,
    later phases sit inside it.
    """
    nested = _phase_file(
        _phase("A", 1.0, 10.0),
        _phase("B", 2.0, 3.0),
        _phase("C", 8.0, 9.0),
    )
    issues = validate_profile_references(
        _root_with(), datasets=(), known_genes=None, published_genes=set(), phases=nested
    )
    assert issues == []


def test_prf006_reports_a_real_gap_bordered_by_merged_overlapping_phases() -> None:
    """A genuine hole must still be found once its borders are themselves
    merged from more than one overlapping phase -- naming every phase on
    each side of the gap, not an arbitrary one of them.

    `A=[1,3)` and `B=[2,4)` overlap and merge into one coverage span
    `[1,4)`; `C=[6,9)` is a separate span. The gap `[4,6)` is real and must
    be reported exactly once, crediting both `A` and `B` on its near side.
    """
    phases = _phase_file(
        _phase("A", 1.0, 3.0),
        _phase("B", 2.0, 4.0),
        _phase("C", 6.0, 9.0),
    )
    issues = validate_profile_references(
        _root_with(), datasets=(), known_genes=None, published_genes=set(), phases=phases
    )
    assert [i.code for i in issues] == ["PRF006"]
    assert "'A'" in issues[0].message and "'B'" in issues[0].message and "'C'" in issues[0].message
    assert "4.0" in issues[0].message and "6.0" in issues[0].message


def test_prf006_does_not_let_an_open_ended_phase_close_a_gap_it_was_never_shown_to_close() -> None:
    """The real vocabulary's own shape: `heart_looping` starts at 3.14, well
    before the true gap [3.29, 3.71), and has no stated end -- an
    implementation that let it contribute coverage regardless (treating its
    missing end as "runs forever", or defaulting it to some large number)
    would paper over the gap this test's fixture reproduces exactly.

    `looping=[1.0, no end)`, `A=[0.5, 2.0)`, `B=[3.0, 5.0)`: without
    `looping`'s exclusion, a naive implementation might read `looping` as
    covering from 1.0 onward and merge straight through to `B`, reporting no
    gap at all. With the exclusion applied correctly, `A`'s span is `[0.5,
    2.0)`, `looping` contributes nothing, and `B`'s span is `[3.0, 5.0)`:
    exactly one real gap, `[2.0, 3.0)`.
    """
    phases = _phase_file(
        _phase("A", 0.5, 2.0),
        _open_phase("looping", 1.0),
        _phase("B", 3.0, 5.0),
    )
    issues = validate_profile_references(
        _root_with(), datasets=(), known_genes=None, published_genes=set(), phases=phases
    )
    assert [i.code for i in issues] == ["PRF006"]
    assert "'A'" in issues[0].message and "'B'" in issues[0].message
    assert "looping" not in issues[0].message
    assert "2.0" in issues[0].message and "3.0" in issues[0].message


def test_prf007_names_a_gene_id_no_registry_knows() -> None:
    """Measured 2026-08-14: without this, such a row costs zero issues.

    A repository carrying a profiles row for HGNC:99999 -- in no registry, no
    mirror, nothing -- validated at 0 errors, 3 warnings, identical to
    baseline. There is no REF001 equivalent for a mirror table's gene column:
    REF001 covers curated records, `validate_mirror_references` checks
    `dataset` and `contrast`, and TBL005 checks the id's *shape*.

    HGNC withdraws and merges ids routinely, so a matrix mapped through an
    older release yields pattern-valid ids whose whole expression evidence
    publishes into a shard and reaches no bundle, under a green build.
    BUR011 is the in-repo precedent.
    """
    issues = validate_profile_references(
        _root_with(
            profiles=[
                _p(gene="HGNC:11604", tissue="Heart"),
                _p(gene="HGNC:99999", tissue="Liver"),
            ]
        ),
        datasets=(),
        known_genes={"HGNC:11604"},
        published_genes=set(),
        phases=None,
    )
    assert [i.code for i in issues] == ["PRF007"]
    assert issues[0].severity is Severity.ERROR
    assert "HGNC:99999" in issues[0].message


def test_prf007_multiple_unknown_genes_are_reported_in_sorted_order() -> None:
    """PRF007 issues are sorted by gene id, not by row or set/hash order.

    `profile_genes - known_genes` is a freshly built `set`; iterating it
    unsorted would let the report's order follow PYTHONHASHSEED, so two
    validate runs over one identical repository could disagree on which gene
    is named first.

    **The string family matters, and this was measured, not assumed.** Task
    5 documented five elements as enough (1-in-120 by permutation counting)
    and measured 0/10 survivors for its own dataset-name and unit fixtures.
    Eight HGNC-id-shaped strings differing only in a trailing digit
    (`HGNC:100001`..`HGNC:100008`) measured *worse* than Task 5's
    five-element fixtures: dropping this exact `sorted()` call survived 1/10
    twice over (2/20), and PRF009's sibling test below (same mutation shape)
    survived 1/10 once (1/20) at the same size. A follow-up probe outside
    pytest -- constructing each candidate `set` directly under 300 explicit
    `PYTHONHASHSEED` values rather than 10 random ones -- measured why:
    numeric-suffix HGNC ids coincidentally iterate in sorted order at
    18-21/300 (6-7%) even at eight elements, where qualitatively different
    words measured **0/300 at both five and eight elements**. Permutation
    counting (1-in-120 for five elements) is the wrong model for how a small
    `set`'s hash-table slots land; it does not hold for strings that differ
    in only a short numeric suffix. The fix here is the string family, not
    the count: eight distinct words, each still prefixed `HGNC:` to read as
    a gene id (nothing downstream of this function checks the shape).
    """
    genes = [
        "HGNC:AMBER",
        "HGNC:BLUE",
        "HGNC:CORAL",
        "HGNC:DENIM",
        "HGNC:EBONY",
        "HGNC:FUCHSIA",
        "HGNC:GOLD",
        "HGNC:HAZEL",
    ]
    root = _root_with(profiles=[_p(gene=gene, tissue="Heart", stage="7wpc") for gene in genes])
    issues = validate_profile_references(
        root, datasets=(), known_genes=set(), published_genes=set(), phases=None
    )
    assert [i.code for i in issues] == ["PRF007"] * len(genes)
    assert [i.message for i in issues] == [
        f"gene {gene} is not in mirrors/genes.tsv" for gene in sorted(genes)
    ]


def test_an_unreadable_registry_skips_prf007_rather_than_reporting_every_row() -> None:
    """`known_genes=None` means "the registry did not load", not "no gene exists".

    Reporting hundreds of dangling references for one unreadable file is the
    cascade REF000/SRC000/ONT000 exist to prevent.
    """
    issues = validate_profile_references(
        _root_with(profiles=[_p(gene="HGNC:99999")]),
        datasets=(),
        known_genes=None,
        published_genes=set(),
        phases=None,
    )
    assert [i.code for i in issues] == []


def test_an_empty_but_loaded_registry_reports_every_gene_unlike_none() -> None:
    """`known_genes=set()` -- a header-only mirrors/genes.tsv -- is loaded and
    genuinely empty, which is different from `None` (could not be read at
    all): every profiles gene is then unknown, correctly. This is also what
    pins the guard as `is not None` rather than a truthy check -- `set()` is
    falsy in Python, and a truthy check would wrongly treat it like `None`.
    """
    issues = validate_profile_references(
        _root_with(profiles=[_p(gene="HGNC:11604")]),
        datasets=(),
        known_genes=set(),
        published_genes=set(),
        phases=None,
    )
    assert [i.code for i in issues] == ["PRF007"]
    assert "HGNC:11604" in issues[0].message


def test_prf009_distinguishes_not_sampled_from_dropped() -> None:
    """PRF004 catches a token matching nothing; the realistic failure is partial.

    A converter that drops NA rows leaves `cardiac_tissues` matching thousands
    of rows while some published genes have no heart row at all -- PRF004
    passes, the cardiac sentence has nothing to render, and tau silently runs
    over six organs. `profile_quantiles` records which cells were assayed, so
    the question can be inverted.

    `known_genes` and `published_genes` are deliberately different sets here
    (registered-but-unpublished HGNC:3030 is in the first, not the second):
    PRF009 must read only `published_genes`, so HGNC:3030's absence from this
    one assayed cell must not be named, while published HGNC:2020's must.
    """
    issues = validate_profile_references(
        _root_with(
            profiles=[_p(gene="HGNC:11604", tissue="Heart", stage="7wpc")],
            quantiles=_grid(tissue="Heart", stage="7wpc"),
        ),
        datasets=(),
        known_genes={"HGNC:11604", "HGNC:2020", "HGNC:3030"},
        published_genes={"HGNC:11604", "HGNC:2020"},
        phases=None,
    )
    assert [i.code for i in issues] == ["PRF009"]
    assert issues[0].severity is Severity.WARNING
    assert "HGNC:2020" in issues[0].message
    assert "HGNC:3030" not in issues[0].message


def test_prf009_groups_by_assayed_cell_not_by_gene_or_by_row() -> None:
    """Aggregation choice, pinned: one PRF009 issue per assayed cell.

    The realistic failure PRF009 guards against -- a converter that drops a
    scattered subset of rows within one tissue/stage slice (dropped NA rows,
    a bad range in a copy-paste) -- leaves every dropped gene sharing the
    *same* cell, so one issue per cell names the whole defect once. Grouping
    by gene instead would split that single root cause into as many issues as
    genes it dropped; grouping by (gene, cell) pair would make the report's
    size the product of published genes and assayed cells (92 x ~150 today).
    Bounded instead by the number of assayed cells.

    Two cells: Heart is missing three published genes, Liver is missing one.
    Grouping by cell costs exactly two issues -- by gene it would cost three
    (HGNC:2 and HGNC:3 only at Heart, HGNC:4 at both, HGNC:1 nowhere), and by
    row it would cost four.
    """
    quantiles = _grid(tissue="Heart", stage="7wpc") + _grid(tissue="Liver", stage="7wpc")
    profiles = [
        _p(gene="HGNC:1", tissue="Heart", stage="7wpc"),
        _p(gene="HGNC:1", tissue="Liver", stage="7wpc"),
        _p(gene="HGNC:2", tissue="Liver", stage="7wpc"),
        _p(gene="HGNC:3", tissue="Liver", stage="7wpc"),
    ]
    published = {"HGNC:1", "HGNC:2", "HGNC:3", "HGNC:4"}

    issues = validate_profile_references(
        _root_with(profiles=profiles, quantiles=quantiles),
        datasets=(),
        known_genes=published,
        published_genes=published,
        phases=None,
    )
    assert [i.code for i in issues] == ["PRF009", "PRF009"]
    heart_issue = next(i for i in issues if "Heart" in i.message)
    liver_issue = next(i for i in issues if "Liver" in i.message)
    assert heart_issue is not liver_issue
    assert "['HGNC:2', 'HGNC:3', 'HGNC:4']" in heart_issue.message
    assert "['HGNC:4']" in liver_issue.message
    assert "HGNC:2" not in liver_issue.message
    assert "HGNC:3" not in liver_issue.message


def test_prf009_issues_and_their_message_lists_are_sorted() -> None:
    """Two independent dropped-`sorted()` mutants, guarded at once -- both
    fixture sizes measured directly rather than reused from Task 5's docstring.

    `_prf009_issues` sorts twice: once across assayed (dataset, tissue,
    stage) triples -- a `set` scan whose unsorted order can follow
    PYTHONHASHSEED -- and once within one triple's missing-gene list, a
    second `set` difference with the same risk.

    **The triple-order guard is a `set` of 3-tuples, not of bare strings, and
    that changes the number needed.** Task 5's own five-tissue fixture
    (`Alpha, Bravo, Middle, Yankee, Zebra`) is what this test used first,
    matching its `PRF002` precedent, and a slow pytest-based run measured it
    surviving 1/10. A fast probe outside pytest -- constructing
    `{("GSE999999", tissue, "7wpc") for tissue in tissues}` directly (the
    actual object `_prf009_issues` sorts, not a bare string) under 500
    explicit `PYTHONHASHSEED` values -- found why: wrapped in a triple with a
    constant dataset and stage, those five words coincidentally iterate
    already sorted at 9/500 (1.8%), where the *same five words as bare
    strings* had measured 0/300 -- the tuple wrapping itself changes the
    hash-table layout enough to matter. Eight NATO-alphabet words
    (`Alpha`..`Hotel`) measured 0/500 in the same tuple-wrapped probe, and
    that is the fixture below.

    **The message-list guard needed a different string family, not just more
    elements.** Task 5 documented five elements as enough for this mutation
    shape (1-in-120 by permutation counting) and measured 0/10 survivors for
    its own fixtures. Genes shaped like `HGNC:1`..`HGNC:9` (differing only in
    a trailing digit) measured worse here at both five *and* eight elements:
    dropping the missing-gene `sorted()` call survived 1/10 once at five and
    once more at eight (and PRF007's sibling test above, same mutation shape,
    survived 1/10 twice at eight). A bare-string probe under 300 seeds found
    why: numeric-suffix HGNC ids coincidentally iterate already sorted at
    18-21/300 (6-7%) regardless of element count in the 8-12 range, where
    qualitatively different words measured 0/300 at both five and eight.
    Permutation counting does not model how a small `set`'s hash-table slots
    land for strings that differ only in a digit, or for a tuple built from
    them. The missing-gene set below is therefore eight distinct *words*,
    still `HGNC:`-prefixed to read as gene ids (nothing downstream of this
    function checks the shape) -- the same fix, and the same eight words, as
    the PRF007 test above.

    Both fixed fixtures measured 10/10 killed in a slow pytest-based
    fresh-process re-run after this rewrite.
    """
    tissues = ["Alpha", "Bravo", "Charlie", "Delta", "Echo", "Foxtrot", "Golf", "Hotel"]
    quantiles = [row for tissue in tissues for row in _grid(tissue=tissue, stage="7wpc")]
    candidates = [
        "HGNC:AMBER",
        "HGNC:BLUE",
        "HGNC:CORAL",
        "HGNC:DENIM",
        "HGNC:EBONY",
        "HGNC:FUCHSIA",
        "HGNC:GOLD",
        "HGNC:HAZEL",
    ]
    published = {"HGNC:PRESENT", *candidates}
    profiles = [_p(gene="HGNC:PRESENT", tissue=tissue, stage="7wpc") for tissue in tissues]
    # Every tissue but Hotel also carries every candidate gene but one
    # ("HGNC:HAZEL"), so each still produces its own single-gene PRF009 --
    # Hotel alone carries none of the eight, guarding the eight-element
    # missing-gene list sort.
    profiles += [
        _p(gene=gene, tissue=tissue, stage="7wpc")
        for tissue in tissues[:-1]
        for gene in candidates[:-1]
    ]

    issues = validate_profile_references(
        _root_with(profiles=profiles, quantiles=quantiles),
        datasets=(),
        known_genes=published,
        published_genes=published,
        phases=None,
    )
    assert [i.code for i in issues] == ["PRF009"] * len(tissues)
    for index, tissue in enumerate(tissues):
        assert tissue in issues[index].message
    assert "['HGNC:HAZEL']" in issues[0].message
    assert (
        "['HGNC:AMBER', 'HGNC:BLUE', 'HGNC:CORAL', 'HGNC:DENIM', 'HGNC:EBONY', "
        "'HGNC:FUCHSIA', 'HGNC:GOLD', 'HGNC:HAZEL']" in issues[-1].message
    )


def test_prf009_reports_nothing_when_no_quantile_shard_can_be_read() -> None:
    """No assayed-cell evidence exists at all, so there is nothing to invert.

    Unlike PRF001/PRF002 in `validate_profiles` -- which subtract *from* the
    profiles side, so treating an unread quantile mirror as merely empty
    would report every profiles triple as missing -- PRF009 iterates the
    quantile side directly. An unread `profile_quantiles` mirror makes
    `quantile_triples` empty by construction, and the loop over it is then
    simply a no-op: no separate `quantile_available` guard is needed here,
    unlike Task 5's flag for PRF001/PRF002. Measured directly rather than
    reasoned about: a published gene with a plausible profiles gap, and no
    profile_quantiles mirror at all, must report nothing.
    """
    issues = validate_profile_references(
        _root_with(profiles=[_p(gene="HGNC:1", tissue="Heart", stage="7wpc")]),
        datasets=(),
        known_genes={"HGNC:1", "HGNC:2"},
        published_genes={"HGNC:1", "HGNC:2"},
        phases=None,
    )
    assert [i.code for i in issues] == []
