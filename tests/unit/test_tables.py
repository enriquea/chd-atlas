# tests/unit/test_tables.py
from pathlib import Path

import polars as pl
import pytest

from chd_atlas.tables import (
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
    }


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
    path.write_text(
        "gene\tpos\ttier\tnote\nHGNC:1\t100\ta\t\nHGNC:1\t100\tb\t\nHGNC:2\t200\ta\t\n"
    )

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
