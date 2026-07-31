# src/chd_atlas/tables.py
"""Schema definition and validation for the sorted TSV tables under ``mirrors/``.

Polars rather than pandas: pandas coerces a nullable integer column to float,
which would silently turn a genomic coordinate of 12345 into 12345.0.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Final

import polars as pl

# Imported from polars.datatypes, NOT as pl.DataTypeClass: the top-level alias was
# deprecated in polars 1.0 and survives only as a runtime shim, so `pl.DataTypeClass`
# fails `mypy --strict` with "Name is not defined" on every version we support.
from polars.datatypes import DataTypeClass

from chd_atlas.identifiers import (
    HGNC_PATTERN,
    MODIFICATION_PATTERN,
    SEQUENCE_ONTOLOGY_PATTERN,
)
from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.vocab import SourceTier


@dataclass(frozen=True)
class Column:
    name: str
    # DataTypeClass, not DataType: we always pass the class (pl.Int64), never an instance.
    dtype: DataTypeClass
    nullable: bool = False
    allowed: frozenset[str] | None = None
    pattern: str | None = None


@dataclass(frozen=True)
class TableSchema:
    name: str
    columns: tuple[Column, ...]
    sort_key: tuple[str, ...]

    @property
    def column_names(self) -> tuple[str, ...]:
        return tuple(column.name for column in self.columns)

    def polars_overrides(self) -> dict[str, DataTypeClass]:
        return {column.name: column.dtype for column in self.columns}


def read_table(path: Path, schema: TableSchema) -> pl.DataFrame:
    """Read a TSV with declared dtypes, treating the empty string as null."""
    present = pl.read_csv(path, separator="\t", n_rows=0).columns
    overrides = {
        name: dtype for name, dtype in schema.polars_overrides().items() if name in present
    }
    return pl.read_csv(
        path,
        separator="\t",
        schema_overrides=overrides,
        null_values=[""],
        infer_schema_length=0 if not overrides else None,
    )


def validate_table(path: Path, schema: TableSchema) -> list[ValidationIssue]:
    """Check one TSV file against its schema, returning every problem found."""
    issues: list[ValidationIssue] = []

    def error(code: str, location: str, message: str) -> None:
        issues.append(ValidationIssue(code, Severity.ERROR, location, message))

    frame = read_table(path, schema)
    present = set(frame.columns)
    expected = set(schema.column_names)

    for missing in sorted(expected - present):
        error("TBL001", str(path), f"missing required column '{missing}'")
    for unexpected in sorted(present - expected):
        error("TBL002", str(path), f"unexpected column '{unexpected}'")

    for column in schema.columns:
        if column.name not in present:
            continue
        series = frame[column.name]

        if not column.nullable:
            for index in _null_row_indices(series):
                error(
                    "TBL003",
                    f"{path}:row {index}",
                    f"column '{column.name}' must not be empty",
                )

        if column.allowed is not None:
            for index, value in _string_values(series):
                if value not in column.allowed:
                    error(
                        "TBL004",
                        f"{path}:row {index}",
                        f"column '{column.name}' value '{value}' is not one of "
                        f"{sorted(column.allowed)}",
                    )

        if column.pattern is not None:
            compiled = re.compile(column.pattern)
            for index, value in _string_values(series):
                if not compiled.fullmatch(value):
                    error(
                        "TBL005",
                        f"{path}:row {index}",
                        f"column '{column.name}' value '{value}' does not match "
                        f"{column.pattern}",
                    )

    return issues


def _null_row_indices(series: pl.Series) -> list[int]:
    """Row numbers (1-based, header counted as row 1) holding a null."""
    return [index + 2 for index, value in enumerate(series.to_list()) if value is None]


def _string_values(series: pl.Series) -> list[tuple[int, str]]:
    return [
        (index + 2, str(value))
        for index, value in enumerate(series.to_list())
        if value is not None
    ]


_SOURCE_TIERS = frozenset(tier.value for tier in SourceTier)
_DIRECTIONS = frozenset({"up", "down", "unchanged"})

GENES = TableSchema(
    name="genes",
    columns=(
        Column("hgnc_id", pl.String, pattern=HGNC_PATTERN),
        Column("symbol", pl.String),
        Column("name", pl.String),
        Column("aliases", pl.String, nullable=True),
        Column("ensembl_gene", pl.String, nullable=True),
        Column("ncbi_gene", pl.Int64, nullable=True),
        Column("locus", pl.String, nullable=True),
        Column("uniprot", pl.String, nullable=True),
        Column("mane_select", pl.String, nullable=True),
    ),
    sort_key=("hgnc_id",),
)

PTM_SITES = TableSchema(
    name="ptm_sites",
    columns=(
        Column("site_id", pl.String),
        Column("protein", pl.String),
        Column("residue", pl.String, allowed=frozenset({"S", "T", "Y", "K", "R", "C", "N", "Q"})),
        Column("position", pl.Int64),
        Column("mod_type", pl.String, pattern=MODIFICATION_PATTERN),
        Column("flanking_sequence", pl.String, nullable=True),
        Column("known_kinases", pl.String, nullable=True),
        Column("source", pl.String),
    ),
    sort_key=("protein", "position", "mod_type"),
)

VARIANTS = TableSchema(
    name="variants",
    columns=(
        Column("vrs_id", pl.String),
        Column("assembly", pl.String, allowed=frozenset({"GRCh38"})),
        Column("chrom", pl.String),
        Column("pos", pl.Int64),
        Column("ref", pl.String),
        Column("alt", pl.String),
        Column("hgvs_g", pl.String, nullable=True),
        Column("hgvs_c", pl.String, nullable=True),
        Column("hgvs_p", pl.String, nullable=True),
        Column("mane_transcript", pl.String, nullable=True),
        Column("gene", pl.String, pattern=HGNC_PATTERN),
        Column("consequence", pl.String, pattern=SEQUENCE_ONTOLOGY_PATTERN),
        Column("clinvar_vcv", pl.String, nullable=True),
        Column("clinvar_significance", pl.String, nullable=True),
        Column("clinvar_review_status", pl.String, nullable=True),
        Column("gnomad_af_popmax", pl.Float64, nullable=True),
        Column("gnomad_version", pl.String, nullable=True),
        Column("acmg_criteria", pl.String, nullable=True),
        Column("inheritance", pl.String, nullable=True),
        Column("de_novo", pl.Boolean, nullable=True),
        Column("zygosity", pl.String, nullable=True),
        Column("phenotypes", pl.String, nullable=True),
        Column("publication", pl.String, nullable=True),
        Column("locator", pl.String, nullable=True),
        Column("source_tier", pl.String, allowed=_SOURCE_TIERS),
        Column("source", pl.String),
        # Non-coding extension (design spec D10); null for coding variants.
        Column("regulatory_element", pl.String, nullable=True),
        Column("target_gene", pl.String, nullable=True),
        Column("assay_evidence", pl.String, nullable=True),
    ),
    sort_key=("chrom", "pos", "ref", "alt", "vrs_id"),
)

EXPRESSION = TableSchema(
    name="expression",
    columns=(
        Column("dataset", pl.String),
        Column("contrast", pl.String),
        Column("gene", pl.String, pattern=HGNC_PATTERN),
        Column("log2fc", pl.Float64),
        Column("pvalue", pl.Float64),
        Column("fdr", pl.Float64),
        Column("direction", pl.String, allowed=_DIRECTIONS),
        Column("n_case", pl.Int64, nullable=True),
        Column("n_control", pl.Int64, nullable=True),
        Column("tissue", pl.String),
        Column("stage", pl.String, nullable=True),
    ),
    sort_key=("contrast", "gene"),
)

PROFILES = TableSchema(
    name="profiles",
    columns=(
        Column("dataset", pl.String),
        Column("gene", pl.String, pattern=HGNC_PATTERN),
        Column("tissue", pl.String),
        Column("stage", pl.String, nullable=True),
        Column("median_abundance", pl.Float64),
        Column("unit", pl.String, allowed=frozenset({"tpm", "nx", "cpm", "lfq"})),
        Column("q25", pl.Float64, nullable=True),
        Column("q75", pl.Float64, nullable=True),
        Column("n_samples", pl.Int64),
    ),
    sort_key=("dataset", "gene", "tissue", "stage"),
)

PROTEOMICS = TableSchema(
    name="proteomics",
    columns=(
        Column("dataset", pl.String),
        Column("contrast", pl.String),
        Column("protein", pl.String),
        Column("gene", pl.String, pattern=HGNC_PATTERN, nullable=True),
        Column("log2fc", pl.Float64),
        Column("pvalue", pl.Float64),
        Column("fdr", pl.Float64),
        Column("direction", pl.String, allowed=_DIRECTIONS),
        Column("n_peptides", pl.Int64, nullable=True),
        Column("sequence_coverage", pl.Float64, nullable=True),
        Column("quant_method", pl.String),
    ),
    sort_key=("contrast", "protein"),
)

PHOSPHO = TableSchema(
    name="phospho",
    columns=(
        Column("dataset", pl.String),
        Column("contrast", pl.String),
        Column("site_id", pl.String),
        Column("protein", pl.String),
        Column("residue", pl.String, allowed=frozenset({"S", "T", "Y"})),
        Column("position", pl.Int64),
        Column("mod_type", pl.String, pattern=MODIFICATION_PATTERN),
        Column("flanking_sequence", pl.String, nullable=True),
        Column("log2fc", pl.Float64),
        Column("pvalue", pl.Float64),
        Column("fdr", pl.Float64),
        # Mandatory: whether the change was corrected for protein abundance decides
        # whether the site is genuinely regulated or merely tracking total protein.
        Column("protein_normalized", pl.Boolean),
        Column("occupancy", pl.Float64, nullable=True),
        Column("known_kinases", pl.String, nullable=True),
    ),
    sort_key=("contrast", "protein", "position", "mod_type"),
)

TABLE_SCHEMAS: Final[dict[str, TableSchema]] = {
    schema.name: schema
    for schema in (GENES, PTM_SITES, VARIANTS, EXPRESSION, PROFILES, PROTEOMICS, PHOSPHO)
}

# Which directory layout each schema lives under, relative to ``mirrors/``.
SHARDED_TABLES: Final[dict[str, str]] = {
    "variants": "variants",
    "expression": "expression",
    "profiles": "profiles",
    "proteomics": "proteomics",
    "phospho": "phospho",
}
FLAT_TABLES: Final[dict[str, str]] = {
    "genes": "genes.tsv",
    "ptm_sites": "ptm_sites.tsv",
}


def mirror_paths(root: Path) -> list[tuple[Path, str]]:
    """Every mirror TSV under ``root/mirrors``, paired with its schema name."""
    paths: list[tuple[Path, str]] = []
    mirrors = root / "mirrors"
    for schema_name, filename in sorted(FLAT_TABLES.items()):
        path = mirrors / filename
        if path.is_file():
            paths.append((path, schema_name))
    for schema_name, directory in sorted(SHARDED_TABLES.items()):
        for path in sorted((mirrors / directory).glob("*.tsv")):
            paths.append((path, schema_name))
    return paths
