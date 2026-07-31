# tests/unit/test_sort_order.py
from pathlib import Path

import polars as pl
import pytest

from chd_atlas.tables import TABLE_SCHEMAS, Column, TableSchema
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


def _write(path: Path, header: list[str], rows: list[list[str]]) -> Path:
    lines = ["\t".join(header)] + ["\t".join(row) for row in rows]
    path.write_text("\n".join(lines) + "\n")
    return path


def test_real_variants_schema_accepts_a_correctly_sorted_shard(tmp_path: Path) -> None:
    schema = TABLE_SCHEMAS["variants"]
    rows = [
        ["1", "100", "A", "T", "ga4gh:VA.aaa"],
        ["1", "200", "C", "G", "ga4gh:VA.bbb"],
        ["2", "50", "G", "A", "ga4gh:VA.ccc"],
    ]
    path = _write(tmp_path / "variants.tsv", list(schema.sort_key), rows)

    assert validate_sort_order(path, schema) == []


def test_real_variants_schema_rejects_a_misordered_shard(tmp_path: Path) -> None:
    schema = TABLE_SCHEMAS["variants"]
    rows = [
        ["1", "200", "C", "G", "ga4gh:VA.bbb"],
        ["1", "100", "A", "T", "ga4gh:VA.aaa"],
    ]
    path = _write(tmp_path / "variants.tsv", list(schema.sort_key), rows)

    issues = validate_sort_order(path, schema)

    assert [i.code for i in issues] == ["SORT001"]
    assert "row 2" in issues[0].message


@pytest.mark.parametrize("name", sorted(TABLE_SCHEMAS))
def test_every_schema_sort_key_uses_declared_columns(name: str) -> None:
    """A sort key naming a column the schema does not declare is silently unenforceable."""
    schema = TABLE_SCHEMAS[name]
    assert set(schema.sort_key) <= set(schema.column_names)
