# src/chd_atlas/tables.py
"""Schema definition and validation for the sorted TSV tables under ``mirrors/``.

Polars rather than pandas: pandas coerces a nullable integer column to float,
which would silently turn a genomic coordinate of 12345 into 12345.0.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Final

import polars as pl

# Imported from polars.datatypes, NOT as pl.DataTypeClass: the top-level alias was
# deprecated in polars 1.0 and survives only as a runtime shim, so `pl.DataTypeClass`
# fails `mypy --strict` with "Name is not defined" on every version we support.
from polars.datatypes import DataTypeClass
from polars.exceptions import PolarsError

from chd_atlas.duplicates import duplicates
from chd_atlas.fs import list_dir
from chd_atlas.identifiers import (
    HGNC_PATTERN,
    MODIFICATION_PATTERN,
    SEQUENCE_ONTOLOGY_PATTERN,
    UNIPROT_PATTERN,
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
    minimum: float | None = None
    maximum: float | None = None


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


def read_table(
    path: Path, schema: TableSchema
) -> tuple[pl.DataFrame | None, list[ValidationIssue]]:
    """Read a TSV with declared dtypes, treating the empty string as null.

    Returns ``(None, issues)`` when the file cannot be read at all — a non-UTF-8
    byte, a zero-length file, or a missing path. Polars raises for the first two
    and a missing path raises ``FileNotFoundError``; an uncaught raise here would
    abort validation of every other shard in the repository rather than
    reporting one unreadable file.
    """
    try:
        present = pl.read_csv(path, separator="\t", n_rows=0).columns
        overrides = {
            name: dtype
            for name, dtype in schema.polars_overrides().items()
            if name in present
        }
        frame = pl.read_csv(
            path,
            separator="\t",
            schema_overrides=overrides,
            null_values=[""],
            infer_schema_length=0 if not overrides else None,
        )
    except (PolarsError, OSError) as exc:
        issue = ValidationIssue(
            "TBL000", Severity.ERROR, str(path), f"could not read TSV: {exc}"
        )
        return None, [issue]
    return frame, []


# Past this many offending rows in one column, report a count instead. A column
# renamed upstream makes every row fail, and mirrors are reviewed by digest.
_MAX_ROW_ISSUES: Final = 20


def _row_issues(
    code: str, path: Path, offenders: list[tuple[int, str]], summary: str
) -> list[ValidationIssue]:
    """One issue per offending row, collapsed to a count past the threshold."""
    if not offenders:
        return []
    if len(offenders) > _MAX_ROW_ISSUES:
        first_row, first_detail = offenders[0]
        return [
            ValidationIssue(
                code,
                Severity.ERROR,
                str(path),
                f"{len(offenders)} rows: {summary} (first at row {first_row}: {first_detail})",
            )
        ]
    return [
        ValidationIssue(code, Severity.ERROR, f"{path}:row {row}", f"{summary}: {detail}")
        for row, detail in offenders
    ]


def validate_table(path: Path, schema: TableSchema) -> list[ValidationIssue]:
    """Check one TSV file against its schema, returning every problem found."""
    issues: list[ValidationIssue] = []

    def error(code: str, location: str, message: str) -> None:
        issues.append(ValidationIssue(code, Severity.ERROR, location, message))

    frame, read_issues = read_table(path, schema)
    if frame is None:
        return read_issues
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
            offenders_null = [(row, "empty") for row in _null_row_indices(series)]
            issues.extend(
                _row_issues(
                    "TBL003", path, offenders_null, f"column '{column.name}' must not be empty"
                )
            )

        # Keyed on the dtype rather than on the bounds, because a column with no
        # natural range (a log2 fold change) needs this exactly as much as a
        # bounded one. Null and NaN are distinct states in polars: the non-null
        # check above tests `value is None` and so never sees a NaN, while every
        # comparison against NaN is False by IEEE-754 and so the bounds check
        # below cannot fire for one either. Without this, a p-value written as
        # `NaN` — what R's `write.table` and numpy's `savetxt` emit for a missing
        # statistic — validates clean and ships as a real measurement.
        if column.dtype is pl.Float64:
            offenders_finite = [
                (row, repr(value))
                for row, value in _numeric_values(series)
                if not math.isfinite(value)
            ]
            issues.extend(
                _row_issues(
                    "TBL010",
                    path,
                    offenders_finite,
                    f"column '{column.name}' has non-finite values",
                )
            )

        if column.allowed is not None:
            allowed_sorted = sorted(column.allowed)
            offenders_allowed = [
                (row, repr(value))
                for row, value in _string_values(series)
                if value not in column.allowed
            ]
            issues.extend(
                _row_issues(
                    "TBL004",
                    path,
                    offenders_allowed,
                    f"column '{column.name}' has values outside {allowed_sorted}",
                )
            )

        if column.pattern is not None:
            compiled = re.compile(column.pattern)
            offenders_pattern = [
                (row, repr(value))
                for row, value in _string_values(series)
                if not compiled.fullmatch(value)
            ]
            issues.extend(
                _row_issues(
                    "TBL005",
                    path,
                    offenders_pattern,
                    f"column '{column.name}' has values not matching {column.pattern}",
                )
            )

        if column.minimum is not None or column.maximum is not None:
            low = column.minimum if column.minimum is not None else float("-inf")
            high = column.maximum if column.maximum is not None else float("inf")
            # Non-finite values are excluded so exactly one code fires per
            # problem: an infinity is out of every finite range, and reporting it
            # as both TBL010 and TBL006 would say the same thing twice while
            # TBL010 names the cause more precisely.
            offenders_bounds = [
                (row, repr(value))
                for row, value in _numeric_values(series)
                if math.isfinite(value) and (value < low or value > high)
            ]
            issues.extend(
                _row_issues(
                    "TBL006",
                    path,
                    offenders_bounds,
                    f"column '{column.name}' has values outside [{low}, {high}]",
                )
            )

    if all(column in present for column in schema.sort_key):
        keys = frame.select(schema.sort_key).rows()
        repeated = set(duplicates(keys))
        offenders_dupe: list[tuple[int, str]] = []
        seen: set[tuple[object, ...]] = set()
        for index, key in enumerate(keys):
            if key not in repeated:
                continue
            if key in seen:
                offenders_dupe.append((index + 2, repr(key)))
            else:
                seen.add(key)
        issues.extend(
            _row_issues(
                "TBL007", path, offenders_dupe, f"duplicate {list(schema.sort_key)} key"
            )
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


def _numeric_values(series: pl.Series) -> list[tuple[int, float]]:
    return [
        (index + 2, float(value))
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
        Column("uniprot", pl.String, nullable=True, pattern=UNIPROT_PATTERN),
        Column("mane_select", pl.String, nullable=True),
    ),
    sort_key=("hgnc_id",),
)

PTM_SITES = TableSchema(
    name="ptm_sites",
    columns=(
        Column("site_id", pl.String),
        Column("protein", pl.String, pattern=UNIPROT_PATTERN),
        Column("residue", pl.String, allowed=frozenset({"S", "T", "Y", "K", "R", "C", "N", "Q"})),
        Column("position", pl.Int64, minimum=1),
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
        Column(
            "chrom",
            pl.String,
            allowed=frozenset({str(n) for n in range(1, 23)} | {"X", "Y", "MT"}),
        ),
        Column("pos", pl.Int64, minimum=1),
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
        Column("pvalue", pl.Float64, minimum=0, maximum=1),
        Column("fdr", pl.Float64, minimum=0, maximum=1),
        Column("direction", pl.String, allowed=_DIRECTIONS),
        Column("n_case", pl.Int64, nullable=True, minimum=0),
        Column("n_control", pl.Int64, nullable=True, minimum=0),
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
        Column("n_samples", pl.Int64, minimum=1),
    ),
    sort_key=("dataset", "gene", "tissue", "stage"),
)

PROTEOMICS = TableSchema(
    name="proteomics",
    columns=(
        Column("dataset", pl.String),
        Column("contrast", pl.String),
        Column("protein", pl.String, pattern=UNIPROT_PATTERN),
        Column("gene", pl.String, pattern=HGNC_PATTERN, nullable=True),
        Column("log2fc", pl.Float64),
        Column("pvalue", pl.Float64, minimum=0, maximum=1),
        Column("fdr", pl.Float64, minimum=0, maximum=1),
        Column("direction", pl.String, allowed=_DIRECTIONS),
        Column("n_peptides", pl.Int64, nullable=True, minimum=1),
        Column("sequence_coverage", pl.Float64, nullable=True, minimum=0, maximum=100),
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
        Column("protein", pl.String, pattern=UNIPROT_PATTERN),
        Column("residue", pl.String, allowed=frozenset({"S", "T", "Y"})),
        Column("position", pl.Int64, minimum=1),
        Column("mod_type", pl.String, pattern=MODIFICATION_PATTERN),
        Column("flanking_sequence", pl.String, nullable=True),
        Column("log2fc", pl.Float64),
        Column("pvalue", pl.Float64, minimum=0, maximum=1),
        Column("fdr", pl.Float64, minimum=0, maximum=1),
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


def unexpected_mirror_entries(root: Path) -> list[ValidationIssue]:
    """Anything under ``mirrors/`` that no schema claims.

    ``Path.glob`` yields nothing for a missing or non-directory path, so a shard
    directory lost to a typo — ``mirrors/varaints`` — simply vanishes from
    validation and the repository reports clean. Naming the stray entry is what
    catches that. An *absent* shard directory stays legitimate: before any
    variants are curated there is no ``mirrors/variants/``.
    """
    mirrors = root / "mirrors"
    if not mirrors.is_dir():
        return []

    expected_dirs = set(SHARDED_TABLES.values())
    expected_files = set(FLAT_TABLES.values()) | {"sources.yaml"}
    issues: list[ValidationIssue] = []

    entries, listing_issues = list_dir(mirrors, "TBL009")
    issues.extend(listing_issues)

    for entry in entries:
        if entry.is_dir():
            if entry.name not in expected_dirs:
                issues.append(
                    ValidationIssue(
                        "TBL009",
                        Severity.ERROR,
                        str(entry),
                        f"unexpected directory under mirrors/; expected one of "
                        f"{sorted(expected_dirs)}",
                    )
                )
                continue
            # mirror_paths globs *.tsv, so a shard given the wrong extension is
            # invisible to every table check while the gate still passes.
            shards, shard_listing_issues = list_dir(entry, "TBL009")
            issues.extend(shard_listing_issues)
            for shard in shards:
                if not shard.is_file():
                    issues.append(
                        ValidationIssue(
                            "TBL009",
                            Severity.ERROR,
                            str(shard),
                            f"'{shard.name}' should be a .tsv shard, not a directory",
                        )
                    )
                elif shard.suffix != ".tsv":
                    issues.append(
                        ValidationIssue(
                            "TBL009",
                            Severity.ERROR,
                            str(shard),
                            f"shard files must end .tsv; '{shard.name}' does not",
                        )
                    )
        # A shard directory replaced by a regular file is checked before the
        # unexpected-file case: the two name sets are disjoint, so testing
        # membership of `expected_files` first would report `mirrors/variants`
        # as merely unnamed rather than as the wrong kind of entry.
        elif entry.name in expected_dirs:
            issues.append(
                ValidationIssue(
                    "TBL009",
                    Severity.ERROR,
                    str(entry),
                    f"'{entry.name}' should be a directory of shards, not a file",
                )
            )
        elif entry.name not in expected_files:
            issues.append(
                ValidationIssue(
                    "TBL009",
                    Severity.ERROR,
                    str(entry),
                    f"unexpected file under mirrors/; expected one of "
                    f"{sorted(expected_files)}",
                )
            )
    return issues


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
