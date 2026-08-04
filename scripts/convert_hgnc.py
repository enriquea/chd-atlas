# scripts/convert_hgnc.py
"""Convert the HGNC complete set into `mirrors/genes.tsv`.

Run locally by a curator, never by CI or the build -- design decision D6: CI
validates, builds and deploys; it never authors data. The output is committed,
so a refresh is a reviewable diff rather than a silent change under a build.

    curl -sSL -o /tmp/hgnc.tsv \\
        https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt
    .venv/bin/python scripts/convert_hgnc.py /tmp/hgnc.tsv mirrors/genes.tsv

With no third argument the keep set is every gene named by a row of
`mirrors/clingen_gene_validity.tsv` or `mirrors/gencc_submissions.tsv` whose
disease is listed in `curation/chd_scope.yaml` -- 154 genes as committed. The
registry is deliberately wider than what the site publishes (23 genes under
D21), because it is a lookup table rather than published content: widening the
gate later must not require re-mirroring.

Two columns need a rule, and both were measured against the live file rather
than assumed:

`mane_select` is published by HGNC as "ENST00000405440.7|NM_181486.4" for all
154 in-scope genes. The mirror carries the RefSeq half, which is what the
committed TBX5 row already holds and what a variant's `mane_transcript` is
written against.

`uniprot_ids` carries more than one accession for exactly one of the 154 --
RBM10, "P0DW28|P98175" -- and RBM10 is one of the 23 genes D21 publishes, so
this is live. `UNIPROT_PATTERN` rejects a pipe, so the pair cannot be one cell.
The first is taken. Publishing null instead would drop the gene out of every
accession-keyed join with nothing raised, which is worse: a wrong accession
fails to match and can be found, a null looks like a gene with no protein.
`tests/unit/test_convert_hgnc.py` pins the choice.
"""

from __future__ import annotations

import re
import sys
from collections.abc import Collection
from pathlib import Path

import polars as pl

# The source-to-mirror column mapping, for the reader. `convert` spells the
# same mapping out as expressions; a column HGNC renames raises out of polars
# there rather than mirroring a silently absent field.
_COLUMNS: dict[str, str] = {
    "hgnc_id": "hgnc_id",
    "symbol": "symbol",
    "name": "name",
    "alias_symbol": "aliases",
    "ensembl_gene_id": "ensembl_gene",
    "entrez_id": "ncbi_gene",
    "location": "locus",
    "uniprot_ids": "uniprot",
    "mane_select": "mane_select",
}

# The order `TABLE_SCHEMAS["genes"]` declares. Written out rather than imported
# so this script stays runnable without the package installed.
_ORDER = (
    "hgnc_id",
    "symbol",
    "name",
    "aliases",
    "ensembl_gene",
    "ncbi_gene",
    "locus",
    "uniprot",
    "mane_select",
)

_SCOPE_ID = re.compile(r"^\s*- id: (MONDO:\d+)", re.M)


def _first(column: str) -> pl.Expr:
    """The first pipe-separated value of a cell, or the cell itself."""
    return pl.col(column).str.split("|").list.first()


def _refseq(column: str) -> pl.Expr:
    """The RefSeq half of HGNC's "ENST…|NM_…" MANE Select pair.

    Selected by prefix rather than by position: HGNC publishes the Ensembl
    transcript first today, and a value that is only one of the two must not
    silently publish an Ensembl id in a column every consumer reads as RefSeq.
    A pair with no RefSeq half yields null, which `TABLE_SCHEMAS["genes"]`
    permits.
    """
    return (
        pl.col(column)
        .str.split("|")
        .list.eval(pl.element().filter(pl.element().str.starts_with("NM_")))
        .list.first()
    )


def convert(source: Path, out: Path, keep: Collection[str]) -> int:
    """Write the mirror rows for `keep`. Returns the number of rows written."""
    frame = pl.read_csv(source, separator="\t", infer_schema_length=0)
    wanted = frame.filter(pl.col("hgnc_id").is_in(list(keep)))

    rows = (
        wanted.select(
            pl.col("hgnc_id"),
            pl.col("symbol"),
            pl.col("name"),
            pl.col("alias_symbol").alias("aliases"),
            pl.col("ensembl_gene_id").alias("ensembl_gene"),
            pl.col("entrez_id").alias("ncbi_gene"),
            pl.col("location").alias("locus"),
            _first("uniprot_ids").alias("uniprot"),
            _refseq("mane_select").alias("mane_select"),
        )
        .select(_ORDER)
        # Sorted on the published string, matching `sort_order.py`, so the
        # committed file is already in the order validation requires.
        .sort("hgnc_id")
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    rows.write_csv(out, separator="\t", include_header=True)
    return len(rows)


def _in_scope_genes(root: Path) -> set[str]:
    """Every gene either validity mirror names against an in-scope disease."""
    scope = set(_SCOPE_ID.findall((root / "curation" / "chd_scope.yaml").read_text()))
    genes: set[str] = set()
    for name in ("clingen_gene_validity.tsv", "gencc_submissions.tsv"):
        frame = pl.read_csv(root / "mirrors" / name, separator="\t", infer_schema_length=0)
        genes.update(frame.filter(pl.col("disease").is_in(list(scope)))["gene"].to_list())
    return genes


def main(argv: list[str]) -> int:
    if not 2 <= len(argv) <= 3:
        print(__doc__, file=sys.stderr)
        return 2
    source, out = Path(argv[0]), Path(argv[1])
    root = Path(argv[2]) if len(argv) == 3 else Path.cwd()
    written = convert(source, out, keep=_in_scope_genes(root))
    print(f"wrote {written} gene(s) to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
