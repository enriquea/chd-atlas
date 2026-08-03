"""Convert the GenCC submissions-export TSV into the mirror TSV.

GenCC is a pass-through aggregator, not a primary curator. Its own terms:
"The GenCC does not independently verify the submitted information." It
harmonises ~19 member databases' vocabularies (gene, disease, mode of
inheritance, classification) onto one scale and publishes every submitter's
verdict side by side, with none adjudicated. So the mirror's key is
(gene, disease, moi, submitter) -- the submitter is part of the key on
purpose: two submitters disagreeing about one gene-disease pair is the normal
case here, not an error to resolve.

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

The export carries one row per *submission*, not one row per
(gene, disease, moi, submitter): the four harmonised columns used as this
mirror's key still collide on 133 groups (134 extra rows) out of 30,410,
almost all the same submitter resubmitting the same gene/disease/moi under a
new ``sgc_id`` without GenCC superseding the old one -- a data-quality
artifact in the upstream export, not something curated here. `convert`
resolves each collision to the row with the latest ``submitted_as_date``
(ties broken by ``sgc_id``, which is assigned in submission order), so the
committed mirror is unique on the declared key.
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

# The published mirror's key. `sgc_id` is not part of it -- it is used only to
# break a same-day resubmission tie deterministically -- and is dropped before
# writing, since `GENCC_SUBMISSIONS` does not declare it as a column.
_KEY = ("gene", "disease", "moi", "submitter")

_HGNC_PATTERN = r"^HGNC:\d+$"
_MONDO_PATTERN = r"^MONDO:\d{7}$"


def convert(source: Path, out: Path) -> int:
    """Write the mirror TSV. Returns the number of data rows written."""
    frame = pl.read_csv(source, separator="\t", infer_schema_length=0)
    missing = sorted(set(_COLUMNS) - set(frame.columns))
    if missing:
        raise ValueError(f"{source}: expected columns absent: {missing}")

    total = frame.height

    # Defensive, not currently load-bearing: measured 2026-08-03 against the
    # real export, `gene_curie` is HGNC and `disease_curie` is MONDO on all
    # 30,410/30,410 rows -- GenCC's harmonisation step apparently guarantees
    # this for the columns used here (unlike `disease_original_curie`, which
    # is a mix of MONDO/OMIM/Orphanet). Kept so a future refresh that adds a
    # row keyed on another identifier system fails as one reported count
    # rather than as thousands of TBL004 rows burying everything else.
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

    frame = frame.select(
        [pl.col(source_name).alias(target) for source_name, target in _COLUMNS.items()]
        + [pl.col("sgc_id")]
    )

    # Resolve same-key collisions to the most recently submitted row. Sorting
    # ascending by (key..., submitted_on, sgc_id) and keeping the last member
    # of each key group picks the latest `submitted_on`; an exact tie on that
    # (measured: same-day, same-second resubmissions do occur) falls back to
    # the lexicographically last `sgc_id`, which GenCC assigns in submission
    # order, so the tiebreak still favours the newer entry.
    before_dedupe = frame.height
    frame = frame.sort([*_KEY, "submitted_on", "sgc_id"]).unique(
        subset=list(_KEY), keep="last", maintain_order=True
    )
    dropped_collisions = before_dedupe - frame.height
    if dropped_collisions:
        print(
            f"dropping {dropped_collisions} row(s) that collided on {_KEY} "
            "with an older submission from the same submitter",
            file=sys.stderr,
        )

    frame = frame.drop("sgc_id").sort(list(_KEY))

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
