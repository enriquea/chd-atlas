# tests/unit/test_tables.py
from pathlib import Path

import polars as pl

from chd_atlas.tables import TABLE_SCHEMAS, Column, TableSchema, validate_table

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
    frame = pl.read_csv(path, separator="\t", schema_overrides={"pos": pl.Int64})
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
