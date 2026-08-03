# tests/unit/test_tables.py
import re
from pathlib import Path

import polars as pl
import pytest

from chd_atlas.tables import (
    CHROMOSOMES,
    CLINGEN_VALIDITY,
    GENCC_SUBMISSIONS,
    TABLE_SCHEMAS,
    Column,
    TableSchema,
    mirror_paths,
    read_table,
    unexpected_mirror_entries,
    validate_table,
)

SCHEMA = TableSchema(
    name="demo",
    columns=(
        Column("gene", pl.String, pattern=r"^HGNC:\d+$"),
        Column("pos", pl.Int64),
        Column("tier", pl.String, allowed=frozenset({"a", "b"})),
        Column("note", pl.String, nullable=True),
    ),
    sort_key=("gene", "pos"),
)


def _write(path: Path, rows: str) -> Path:
    path.write_text(rows)
    return path


def test_accepts_a_conforming_table(tmp_path: Path) -> None:
    path = _write(
        tmp_path / "demo.tsv",
        "gene\tpos\ttier\tnote\nHGNC:1\t100\ta\t\nHGNC:2\t200\tb\thello\n",
    )
    assert validate_table(path, SCHEMA) == []


def test_reports_a_missing_column(tmp_path: Path) -> None:
    path = _write(tmp_path / "demo.tsv", "gene\tpos\ttier\nHGNC:1\t100\ta\n")
    issues = validate_table(path, SCHEMA)
    assert [i.code for i in issues] == ["TBL001"]
    assert "note" in issues[0].message


def test_reports_an_unexpected_column(tmp_path: Path) -> None:
    path = _write(
        tmp_path / "demo.tsv",
        "gene\tpos\ttier\tnote\textra\nHGNC:1\t100\ta\t\tx\n",
    )
    issues = validate_table(path, SCHEMA)
    assert [i.code for i in issues] == ["TBL002"]


def test_reports_a_null_in_a_non_nullable_column(tmp_path: Path) -> None:
    path = _write(
        tmp_path / "demo.tsv",
        "gene\tpos\ttier\tnote\nHGNC:1\t\ta\t\n",
    )
    issues = validate_table(path, SCHEMA)
    assert [i.code for i in issues] == ["TBL003"]
    assert "row 2" in issues[0].location


def test_reports_a_value_outside_the_allowed_set(tmp_path: Path) -> None:
    path = _write(
        tmp_path / "demo.tsv",
        "gene\tpos\ttier\tnote\nHGNC:1\t100\tz\t\n",
    )
    issues = validate_table(path, SCHEMA)
    assert [i.code for i in issues] == ["TBL004"]


def test_reports_a_value_failing_its_pattern(tmp_path: Path) -> None:
    path = _write(
        tmp_path / "demo.tsv",
        "gene\tpos\ttier\tnote\nENSG001\t100\ta\t\n",
    )
    issues = validate_table(path, SCHEMA)
    assert [i.code for i in issues] == ["TBL005"]


def test_integer_columns_stay_integers(tmp_path: Path) -> None:
    """Guards the reason polars was chosen over pandas: nullable ints must not become floats."""
    path = _write(
        tmp_path / "demo.tsv",
        "gene\tpos\ttier\tnote\nHGNC:1\t100\ta\t\nHGNC:2\t200\tb\t\n",
    )
    frame, issues = read_table(path, SCHEMA)
    assert issues == []
    assert frame is not None
    assert frame["pos"].dtype == pl.Int64
    assert frame["pos"].to_list() == [100, 200]


def test_registry_covers_every_mirror_table() -> None:
    assert set(TABLE_SCHEMAS) == {
        "genes",
        "ptm_sites",
        "variants",
        "expression",
        "profiles",
        "proteomics",
        "phospho",
        "clingen_validity",
        "gencc_submissions",
    }


def test_clingen_validity_is_keyed_on_the_triple_not_the_pair() -> None:
    """(gene, disease, moi) is the key; (gene, disease) is not.

    Measured 2026-08-03 against the ClinGen bulk CSV: the triple is unique at
    3,653/3,653 rows, while (gene, disease) collapses 59 of them -- pairs curated
    twice under different modes of inheritance with different classifications.
    Sorting on the pair would make those rows compare equal and let a real
    ordering violation pass; more importantly it documents which columns a
    consumer must key on to avoid dropping one of each duplicate.
    """
    assert CLINGEN_VALIDITY.sort_key == ("gene", "disease", "moi")
    assert "sop" in CLINGEN_VALIDITY.column_names
    sop = next(c for c in CLINGEN_VALIDITY.columns if c.name == "sop")
    # Not nullable: 22% of the file predates SOP8 and ClinGen publishes no
    # crosswalk between framework versions, so a classification without the
    # framework that produced it is an unqualified claim.
    assert sop.nullable is False


def test_gencc_keys_on_the_submitter_because_it_publishes_no_verdict() -> None:
    """GenCC aggregates; it does not adjudicate.

    Its own terms: "The GenCC does not independently verify the submitted
    information." Two submitters may disagree about the same gene-disease pair
    and GenCC publishes both. Omitting `submitter` from the key would make those
    rows compare equal and silently keep whichever sorted last.

    `Supportive` is in the allowed set but is NOT a rung on the evidence ladder:
    it is a mapping-exception bucket for submitters that do not grade evidence.
    `vocab.GENCC_CLASSIFICATIONS` (Task 5) maps it to None.
    """
    assert GENCC_SUBMISSIONS.sort_key == ("gene", "disease", "moi", "submitter")
    classification = next(c for c in GENCC_SUBMISSIONS.columns if c.name == "classification")
    assert classification.allowed is not None
    assert "Supportive" in classification.allowed


def test_phospho_protein_normalized_is_mandatory_and_non_nullable() -> None:
    column = next(c for c in TABLE_SCHEMAS["phospho"].columns if c.name == "protein_normalized")
    assert column.nullable is False
    assert column.dtype == pl.Boolean


def test_unreadable_file_is_reported_not_raised(tmp_path: Path) -> None:
    """One bad byte must not abort validation of every other shard."""
    path = tmp_path / "demo.tsv"
    path.write_bytes("gene\tpos\ttier\tnote\nHGNC:1\t100\ta\tFran\xe7ois\n".encode("latin-1"))

    issues = validate_table(path, SCHEMA)

    assert [i.code for i in issues] == ["TBL000"]


def test_empty_file_is_reported_not_raised(tmp_path: Path) -> None:
    path = tmp_path / "demo.tsv"
    path.write_bytes(b"")

    issues = validate_table(path, SCHEMA)

    assert [i.code for i in issues] == ["TBL000"]


def test_column_wide_failure_collapses_to_one_issue(tmp_path: Path) -> None:
    rows = "".join(f"HGNC:{n}\t{n}\tz\t\n" for n in range(1, 51))
    path = tmp_path / "demo.tsv"
    path.write_text("gene\tpos\ttier\tnote\n" + rows)

    issues = validate_table(path, SCHEMA)

    assert [i.code for i in issues] == ["TBL004"]
    assert "50 rows" in issues[0].message


def test_few_failures_are_reported_per_row(tmp_path: Path) -> None:
    path = tmp_path / "demo.tsv"
    path.write_text("gene\tpos\ttier\tnote\nHGNC:1\t1\tz\t\nHGNC:2\t2\tz\t\n")

    issues = validate_table(path, SCHEMA)

    assert [i.code for i in issues] == ["TBL004", "TBL004"]
    assert all("row" in i.location for i in issues)


def test_reports_a_value_outside_its_numeric_bounds(tmp_path: Path) -> None:
    schema = TableSchema(
        name="bounded",
        columns=(Column("p", pl.Float64, minimum=0, maximum=1),),
        sort_key=("p",),
    )
    path = tmp_path / "bounded.tsv"
    path.write_text("p\n5.0\n")

    issues = validate_table(path, schema)

    assert [i.code for i in issues] == ["TBL006"]


@pytest.mark.parametrize("literal", ["NaN", "nan", "inf", "-inf", "Infinity"])
def test_reports_a_non_finite_value_in_a_bounded_column(tmp_path: Path, literal: str) -> None:
    """A missing statistic written as NaN must not validate as a present one.

    R's `write.table` emits `NaN` and numpy's `savetxt` emits `nan`; both parse
    into a Float64 column without a read error. Every comparison against NaN is
    False, so the bounds check cannot fire for it, and polars models NaN and null
    as distinct states, so the non-null check does not see it either.
    """
    schema = TableSchema(
        name="bounded",
        columns=(Column("p", pl.Float64, minimum=0, maximum=1),),
        sort_key=("p",),
    )
    path = tmp_path / "bounded.tsv"
    path.write_text(f"p\n{literal}\n")

    issues = validate_table(path, schema)

    assert [i.code for i in issues] == ["TBL010"]
    assert "row 2" in issues[0].location


def test_reports_a_non_finite_value_in_an_unbounded_column(tmp_path: Path) -> None:
    """The check belongs to the dtype, not to the bounds.

    A log2 fold change has no natural range, so it declares no minimum or
    maximum. It is still a statistic that must be present to be reported.
    """
    schema = TableSchema(
        name="unbounded",
        columns=(Column("log2fc", pl.Float64),),
        sort_key=("log2fc",),
    )
    path = tmp_path / "unbounded.tsv"
    path.write_text("log2fc\nNaN\n")

    issues = validate_table(path, schema)

    assert [i.code for i in issues] == ["TBL010"]


def test_a_nullable_float_column_still_rejects_nan(tmp_path: Path) -> None:
    """Nullable means "may be absent", which the empty string already says.

    NaN is a third state that reads as present, so it stays an error even where
    a null is allowed.
    """
    schema = TableSchema(
        name="optional",
        columns=(Column("score", pl.Float64, nullable=True),),
        sort_key=("score",),
    )
    path = tmp_path / "optional.tsv"
    path.write_text("score\n\nNaN\n")

    issues = validate_table(path, schema)

    assert [i.code for i in issues] == ["TBL010"]
    assert "row 3" in issues[0].location


def test_variants_chrom_rejects_a_chr_prefix() -> None:
    column = next(c for c in TABLE_SCHEMAS["variants"].columns if c.name == "chrom")
    assert column.allowed is not None
    assert "1" in column.allowed
    assert "chr1" not in column.allowed


def test_mirror_paths_finds_flat_and_sharded_tables(tmp_path: Path) -> None:
    (tmp_path / "mirrors" / "variants").mkdir(parents=True)
    (tmp_path / "mirrors" / "genes.tsv").write_text("hgnc_id\n")
    (tmp_path / "mirrors" / "variants" / "12.tsv").write_text("vrs_id\n")

    found = mirror_paths(tmp_path)

    assert [(p.name, schema) for p, schema in found] == [
        ("genes.tsv", "genes"),
        ("12.tsv", "variants"),
    ]


@pytest.mark.parametrize("name", sorted(TABLE_SCHEMAS))
def test_every_real_schema_accepts_a_header_only_file(name: str, tmp_path: Path) -> None:
    """A header-only shard must validate cleanly, and pins each schema's column list."""
    schema = TABLE_SCHEMAS[name]
    path = tmp_path / f"{name}.tsv"
    path.write_text("\t".join(schema.column_names) + "\n")

    assert validate_table(path, schema) == []


def test_reports_a_duplicate_sort_key(tmp_path: Path) -> None:
    path = tmp_path / "demo.tsv"
    path.write_text("gene\tpos\ttier\tnote\nHGNC:1\t100\ta\t\nHGNC:1\t100\tb\t\nHGNC:2\t200\ta\t\n")

    issues = validate_table(path, SCHEMA)

    assert [i.code for i in issues] == ["TBL007"]
    assert "row 3" in issues[0].location


def test_distinct_sort_keys_are_not_reported_as_duplicates(tmp_path: Path) -> None:
    path = tmp_path / "demo.tsv"
    path.write_text("gene\tpos\ttier\tnote\nHGNC:1\t100\ta\t\nHGNC:1\t200\ta\t\n")

    assert validate_table(path, SCHEMA) == []


def test_missing_file_is_reported_not_raised(tmp_path: Path) -> None:
    issues = validate_table(tmp_path / "absent.tsv", SCHEMA)

    assert [i.code for i in issues] == ["TBL000"]


def test_unexpected_mirror_entries_reports_a_directory_no_schema_claims(tmp_path: Path) -> None:
    """`mirrors/variants` renamed to `mirrors/varaints` glob-matches nothing at all."""
    (tmp_path / "mirrors" / "varaints").mkdir(parents=True)

    issues = unexpected_mirror_entries(tmp_path)

    assert [i.code for i in issues] == ["TBL009"]
    assert "unexpected directory" in issues[0].message
    assert issues[0].location.endswith("varaints")


def test_unexpected_mirror_entries_reports_a_shard_directory_that_is_a_file(
    tmp_path: Path,
) -> None:
    (tmp_path / "mirrors").mkdir(parents=True)
    (tmp_path / "mirrors" / "variants").write_text("oops\n")

    issues = unexpected_mirror_entries(tmp_path)

    assert [i.code for i in issues] == ["TBL009"]
    assert "should be a directory of shards" in issues[0].message


def test_unexpected_mirror_entries_reports_a_file_no_schema_claims(tmp_path: Path) -> None:
    (tmp_path / "mirrors").mkdir(parents=True)
    (tmp_path / "mirrors" / "notes.txt").write_text("stray\n")

    issues = unexpected_mirror_entries(tmp_path)

    assert [i.code for i in issues] == ["TBL009"]
    assert "unexpected file" in issues[0].message


def test_unexpected_mirror_entries_accepts_the_declared_layout(tmp_path: Path) -> None:
    """An absent shard directory is legitimate: nothing curated yet."""
    (tmp_path / "mirrors" / "variants").mkdir(parents=True)
    (tmp_path / "mirrors" / "genes.tsv").write_text("hgnc_id\n")
    (tmp_path / "mirrors" / "ptm_sites.tsv").write_text("site_id\n")
    (tmp_path / "mirrors" / "sources.yaml").write_text("sources: []\n")

    assert unexpected_mirror_entries(tmp_path) == []


def test_unexpected_mirror_entries_ignores_a_missing_mirrors_directory(tmp_path: Path) -> None:
    assert unexpected_mirror_entries(tmp_path) == []


def test_a_shard_with_the_wrong_extension_is_reported(tmp_path: Path) -> None:
    """mirror_paths globs *.tsv, so 12.txt is invisible to every table check."""
    (tmp_path / "mirrors" / "variants").mkdir(parents=True)
    (tmp_path / "mirrors" / "variants" / "12.txt").write_text("vrs_id\n")

    issues = unexpected_mirror_entries(tmp_path)

    assert [i.code for i in issues] == ["TBL009"]
    assert "12.txt" in issues[0].location


def test_a_correctly_named_shard_is_not_reported(tmp_path: Path) -> None:
    (tmp_path / "mirrors" / "variants").mkdir(parents=True)
    (tmp_path / "mirrors" / "variants" / "12.tsv").write_text("vrs_id\n")

    assert unexpected_mirror_entries(tmp_path) == []


@pytest.mark.parametrize(
    "schema_name,column_name",
    [
        ("genes", "uniprot"),
        ("ptm_sites", "protein"),
        ("proteomics", "protein"),
        ("phospho", "protein"),
    ],
)
def test_every_protein_column_constrains_the_accession(schema_name: str, column_name: str) -> None:
    """A protein identifier column must not accept free text.

    `UniprotAccession` constrains the curated YAML side, but the mirror tables
    are where the bulk protein data lives, and a column with no pattern took
    anything at all — a gene symbol, a RefSeq accession, an empty placeholder.
    """
    column = next(c for c in TABLE_SCHEMAS[schema_name].columns if c.name == column_name)
    assert column.pattern is not None
    compiled = re.compile(column.pattern)
    assert compiled.fullmatch("Q99593")
    assert compiled.fullmatch("P12345-2")
    assert not compiled.fullmatch("TBX5")
    assert not compiled.fullmatch("NP_852259.1")


def test_the_chromosome_vocabulary_is_free_of_duplicates() -> None:
    """What deriving the schema from this tuple does not make impossible.

    `VARIANTS`'s `chrom` column takes its `allowed` set from `CHROMOSOMES`, so
    the two can no longer disagree about membership — the test that used to
    check that was retired rather than kept as decoration.

    A duplicate entry survives that refactor and fails quietly in two directions
    at once: `frozenset` silently collapses it, so the schema looks right, while
    `build/variants.py`'s `_ORDER` keeps the *last* index, so one chromosome
    sorts where its twin sat. Neither surfaces until a curator files a shard for
    that chromosome.
    """
    assert len(CHROMOSOMES) == len(set(CHROMOSOMES))
    # Karyotype order is the published contract: the variant index drives a
    # chromosome picker, and it is the one property no other test would notice
    # losing, since every shard-emitting test uses a handful of chromosomes.
    assert CHROMOSOMES[:3] == ("1", "2", "3")
    assert CHROMOSOMES[-3:] == ("X", "Y", "MT")


def test_a_variant_shard_named_for_no_chromosome_is_reported(tmp_path: Path) -> None:
    """The filename is the whole addressing scheme, and nothing else checked it.

    `mirrors/variants/<chrom>.tsv` becomes `variants/<chrom>.json.gz`, and a
    consumer turns a chromosome into a URL by that rule alone. `validate_table`
    checks the `chrom` column; the sweep here checked only the extension. So
    `chr12.tsv` validated at 0 errors and 0 warnings and then raised inside the
    build, giving a curator a traceback where a report belongs — measured before
    this rule existed.

    The vocabulary is the schema's own, derived from the same `CHROMOSOMES`
    tuple as the `chrom` column's `allowed` set, so this cannot come to disagree
    with the column it addresses.

    A correctly named shard is asserted alongside, because a rule that fired on
    every shard would be worse than none.
    """
    shards = tmp_path / "mirrors" / "variants"
    shards.mkdir(parents=True)
    (shards / "12.tsv").write_text("")
    (shards / "MT.tsv").write_text("")
    (shards / "chr12.tsv").write_text("")
    (shards / "12_part2.tsv").write_text("")

    issues = unexpected_mirror_entries(tmp_path)

    reported = sorted(Path(issue.location).name for issue in issues if issue.code == "TBL011")
    assert reported == ["12_part2.tsv", "chr12.tsv"]
    assert all("is not one of" in issue.message for issue in issues if issue.code == "TBL011")
