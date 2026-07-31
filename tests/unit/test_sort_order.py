# tests/unit/test_sort_order.py
from pathlib import Path

import polars as pl

from chd_atlas.tables import Column, TableSchema
from chd_atlas.validate.sort_order import validate_sort_order

SCHEMA = TableSchema(
    name="demo",
    columns=(Column("chrom", pl.String), Column("pos", pl.Int64)),
    sort_key=("chrom", "pos"),
)


def test_accepts_a_correctly_sorted_table(tmp_path: Path) -> None:
    path = tmp_path / "demo.tsv"
    path.write_text("chrom\tpos\n1\t100\n1\t200\n2\t50\n")
    assert validate_sort_order(path, SCHEMA) == []


def test_reports_the_first_out_of_order_row(tmp_path: Path) -> None:
    path = tmp_path / "demo.tsv"
    path.write_text("chrom\tpos\n1\t200\n1\t100\n2\t50\n")

    issues = validate_sort_order(path, SCHEMA)

    assert [i.code for i in issues] == ["SORT001"]
    assert "row 3" in issues[0].location


def test_sorts_numerically_not_lexically(tmp_path: Path) -> None:
    """'100' < '90' lexically but 100 > 90 numerically; the typed read must win."""
    path = tmp_path / "demo.tsv"
    path.write_text("chrom\tpos\n1\t90\n1\t100\n")
    assert validate_sort_order(path, SCHEMA) == []


def test_empty_table_is_trivially_sorted(tmp_path: Path) -> None:
    path = tmp_path / "demo.tsv"
    path.write_text("chrom\tpos\n")
    assert validate_sort_order(path, SCHEMA) == []


def test_missing_sort_column_is_skipped_not_crashed(tmp_path: Path) -> None:
    path = tmp_path / "demo.tsv"
    path.write_text("chrom\n1\n")
    assert validate_sort_order(path, SCHEMA) == []
