# scripts/convert_clingen.py
"""Convert the ClinGen gene-validity bulk CSV into the mirror TSV.

Run locally by a curator, never by CI or the build -- design decision D6: CI
validates, builds and deploys; it never authors data. The output is committed,
so a refresh is a reviewable diff rather than a silent change under a build.

    curl -sSL -o /tmp/clingen.csv \\
        https://search.clinicalgenome.org/kb/gene-validity/download
    .venv/bin/python scripts/convert_clingen.py /tmp/clingen.csv \\
        mirrors/clingen_gene_validity.tsv

The source file carries no version tag -- only a same-day `FILE CREATED` banner
stamp -- so `mirrors/sources.yaml` records that stamp as the version. Fetching at
build time instead would drift the corpus between builds with nothing to diff.
"""

from __future__ import annotations

import sys
from pathlib import Path

import polars as pl

# Lines 1-3 are a title block; line 4 is a separator; line 5 is the header.
_PREAMBLE_ROWS = 4
_SEPARATOR = "+++"

_COLUMNS: dict[str, str] = {
    "GENE ID (HGNC)": "gene",
    "GENE SYMBOL": "gene_symbol",
    "DISEASE ID (MONDO)": "disease",
    "DISEASE LABEL": "disease_label",
    "MOI": "moi",
    "SOP": "sop",
    "CLASSIFICATION": "classification",
    "CLASSIFICATION DATE": "classification_date",
    "GCEP": "gcep",
    "ONLINE REPORT": "report_url",
}

_SORT_KEY = ("gene", "disease", "moi")


def convert(source: Path, out: Path) -> int:
    """Write the mirror TSV. Returns the number of data rows written."""
    frame = pl.read_csv(
        source,
        skip_rows=_PREAMBLE_ROWS,
        has_header=True,
        infer_schema_length=0,
    )
    missing = sorted(set(_COLUMNS) - set(frame.columns))
    if missing:
        raise ValueError(f"{source}: expected columns absent: {missing}")

    frame = frame.filter(~pl.col("GENE SYMBOL").str.starts_with(_SEPARATOR))
    frame = frame.select(
        [pl.col(source_name).alias(target) for source_name, target in _COLUMNS.items()]
    ).sort(_SORT_KEY)

    out.parent.mkdir(parents=True, exist_ok=True)
    frame.write_csv(out, separator="\t")
    return frame.height


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print(f"usage: {argv[0]} <clingen.csv> <out.tsv>", file=sys.stderr)
        return 2
    rows = convert(Path(argv[1]), Path(argv[2]))
    print(f"wrote {rows} rows to {argv[2]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
