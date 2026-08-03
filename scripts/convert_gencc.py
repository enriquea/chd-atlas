"""Convert the GenCC submissions-export TSV into the mirror TSV.

GenCC is a pass-through aggregator, not a primary curator. Its own terms:
"The GenCC does not independently verify the submitted information." It
harmonises ~19 member databases' vocabularies (gene, disease, mode of
inheritance, classification) onto one scale and publishes every submitter's
verdict side by side, with none adjudicated. So the mirror's key is
(gene, disease, moi, submitter, sgc_id) -- the submitter is part of the key
on purpose: two submitters disagreeing about one gene-disease pair is the
normal case here, not an error to resolve.

`sgc_id` is in the key too, and every row GenCC publishes is written: a
submitter can also disagree with *itself*. Measured 2026-08-03 against the
real export (30,410 rows, ?format=new), (gene, disease, moi, submitter)
alone still collides on 133 groups, and 70 of those 133 carry more than one
distinct `classification_title` -- e.g. Ambry Genetics submitting both
`Limited` (SGC-102815) and `Strong` (SGC-104042) for HGNC:20226 /
MONDO:0859332 / Autosomal recessive. `version_number` is identical within
every one of the 133 groups, so GenCC marks neither row as superseding the
other, and nothing else in the export says which is current. There is
therefore no principled way to keep one and drop the other without inventing
a verdict GenCC itself does not assert, so this converter keeps every row
GenCC publishes and lets `sgc_id` -- GenCC's own row id, unique across all
30,410 rows -- make the key unique instead.

Run locally by a curator, never by CI or the build -- design decision D6: CI
validates, builds and deploys; it never authors data. The output is committed,
so a refresh is a reviewable diff rather than a silent change under a build.

    curl -sSL -o /tmp/gencc.tsv \\
        'https://search.thegencc.org/download/action/submissions-export-tsv?format=new'
    .venv/bin/python scripts/convert_gencc.py /tmp/gencc.tsv \\
        mirrors/gencc_submissions.tsv

``?format=new`` is mandatory: GenCC's download page states the legacy TSV
format (keyed on a UUID rather than an ``sgc_id``/``version_number`` pair)
"will be removed and unavailable after September 30, 2026." The download
quota is 20 requests per IP per day, a refresh-time concern only, since the
snapshot committed to this repository is what CI and the build read.
"""

from __future__ import annotations

import sys
from pathlib import Path

import polars as pl

# Harmonised (post-GenCC) columns, not the raw `submitted_as_*` ones: the
# latter retain the submitter's own free-text spelling ("Modeare", "STRONG",
# trailing spaces) and would need re-deriving the same vocabulary GenCC has
# already normalised. `gene_curie` and `disease_curie` are GenCC's harmonised
# identifiers -- unlike `disease_original_curie`, which keeps the submitter's
# own OMIM/Orphanet/MONDO id -- so they are the columns that actually satisfy
# `HGNC_PATTERN` / `MONDO_PATTERN`.
_COLUMNS: dict[str, str] = {
    "sgc_id": "sgc_id",
    "gene_curie": "gene",
    "gene_symbol": "gene_symbol",
    "disease_curie": "disease",
    "disease_title": "disease_label",
    "moi_title": "moi",
    "classification_title": "classification",
    "submitter_title": "submitter",
    "submitted_as_date": "submitted_on",
    "submitted_as_public_report_url": "report_url",
}

_SORT_KEY = ("gene", "disease", "moi", "submitter", "sgc_id")

_HGNC_PATTERN = r"^HGNC:\d+$"
_MONDO_PATTERN = r"^MONDO:\d{7}$"


def convert(source: Path, out: Path) -> int:
    """Write the mirror TSV. Returns the number of data rows written."""
    frame = pl.read_csv(source, separator="\t", infer_schema_length=0)
    missing = sorted(set(_COLUMNS) - set(frame.columns))
    if missing:
        raise ValueError(f"{source}: expected columns absent: {missing}")

    total = frame.height

    # Documented no-op today, not dead code: measured 2026-08-03 against the
    # real export, `gene_curie` is HGNC and `disease_curie` is MONDO on all
    # 30,410/30,410 rows, so this drops 0 rows against the file this mirror
    # was built from. It stays because GenCC's harmonisation guarantee is
    # observed, not documented -- unlike `disease_original_curie`, which is a
    # mix of MONDO/OMIM/Orphanet -- and a future refresh that adds a row keyed
    # on another identifier system should fail as one reported count rather
    # than as thousands of TBL004 rows burying everything else.
    curie_ok = frame["gene_curie"].str.contains(_HGNC_PATTERN) & frame[
        "disease_curie"
    ].str.contains(_MONDO_PATTERN)
    dropped_curie = total - int(curie_ok.sum())
    if dropped_curie:
        print(
            f"dropping {dropped_curie} row(s) outside the HGNC/MONDO CURIE forms",
            file=sys.stderr,
        )
    frame = frame.filter(curie_ok)

    # No further filtering: GenCC's row is the unit of publication here, and
    # `sgc_id` in the key (see the module docstring) means every row GenCC
    # publishes -- including a submitter disagreeing with itself -- reaches
    # this file.
    frame = frame.select(
        [pl.col(source_name).alias(target) for source_name, target in _COLUMNS.items()]
    ).sort(_SORT_KEY)

    out.parent.mkdir(parents=True, exist_ok=True)
    frame.write_csv(out, separator="\t")
    return frame.height


def main(argv: list[str]) -> int:
    if len(argv) != 3:
        print(f"usage: {argv[0]} <gencc.tsv> <out.tsv>", file=sys.stderr)
        return 2
    rows = convert(Path(argv[1]), Path(argv[2]))
    print(f"wrote {rows} rows to {argv[2]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
