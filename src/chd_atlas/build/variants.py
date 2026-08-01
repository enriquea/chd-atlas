# src/chd_atlas/build/variants.py
"""Variant shards, one gzipped JSON per chromosome.

Gzipped rather than plain because spec 5.1 reserves `.json.gz` for these shards
and the search index and for nothing else: only there does controlling the
compression level pay for what it costs the client. GitHub Pages serves a
`.json.gz` as `application/gzip` and sets no `Content-Encoding` on it, so a
browser receives raw compressed bytes rather than a decoded body and has to
decompress them itself — in a page, by piping the response through
`DecompressionStream("gzip")`. Plain `.json` is gzipped in transit by Pages
anyway, costs the consumer nothing to read and stays `curl`-friendly, which is
why every other artifact is one.

The chromosome is this table's entire addressing scheme: it names the mirror
file, it names the published shard, and a consumer turns a chromosome into a URL
by that rule alone. Nothing upstream ties the two together — `validate_table`
checks the `chrom` column and `unexpected_mirror_entries` checks only that a
shard is a `.tsv` — so making the filename and the column agree is this module's
job, and it refuses to publish rather than publish a shard filed under the wrong
name.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Final

from chd_atlas.build.emit import Emitter
from chd_atlas.tables import TABLE_SCHEMAS, mirror_paths, read_table

# The chromosomes, in karyotype order. The same 25 values as the `chrom` column's
# `allowed` set — a frozenset, which carries no order, so the order lives here
# and `test_the_chromosome_vocabulary_matches_the_schema` fails if the two drift.
#
# Ordered, because the index is what a consumer builds a chromosome picker from,
# and sorting shard names as strings publishes 1, 10, 11, ..., 19, 2, 20 with MT
# ahead of X. A consumer's only fix for that is to reimplement this rule in
# JavaScript against a vocabulary it has to hardcode; publishing the order costs
# nothing and is the same order every genome browser shows.
CHROMOSOMES: Final[tuple[str, ...]] = (*(str(n) for n in range(1, 23)), "X", "Y", "MT")
_ORDER: Final[dict[str, int]] = {name: index for index, name in enumerate(CHROMOSOMES)}


def _chromosome(path: Path) -> str:
    """Which chromosome one mirror shard holds, from its filename.

    `mirrors/variants/<chr>.tsv` is the documented layout and the filename is
    what becomes the URL, so a stem outside the vocabulary is refused rather
    than published. `omics.py` puts its stem through `paths.slug` instead
    because a dataset accession is an open vocabulary and any string has to be
    made path-safe; here the vocabulary is closed and every member of it is
    already path-safe, so membership is the stronger check and the one that
    catches `chr12.tsv`, `12_part2.tsv` and `12 (copy).tsv` — each of which
    slugs into a perfectly servable URL that names no chromosome.
    """
    if path.stem not in _ORDER:
        raise ValueError(
            f"{path}: a variant shard is named for the chromosome it holds, and "
            f"{path.stem!r} is not one of {list(CHROMOSOMES)}"
        )
    return path.stem


def _refuse_foreign_rows(path: Path, chrom: str, rows: list[dict[str, Any]]) -> None:
    """Fail unless every row in the shard is on the chromosome the shard is named for.

    A chr7 row filed in `12.tsv` is otherwise published inside
    `variants/12.json.gz` under a `chrom` of "12", where a client that resolves
    a chromosome to its shard can never reach it and nothing reports it missing:
    the row is in `dist/`, its checksum verifies, and the build is green. The
    same check catches one shard holding two chromosomes, which the
    `<chr>.tsv` layout has no way to publish: one of the two would have to go
    out under the other's name.

    The offenders are collected through a set and sorted before they are named,
    so the message does not vary between runs on data that does not vary.
    """
    foreign = sorted({str(row.get("chrom")) for row in rows} - {chrom})
    if foreign:
        raise ValueError(
            f"{path}: filed under chromosome {chrom!r} but holds rows for {foreign}; "
            f"a consumer resolving a chromosome to its shard would never reach them"
        )


def build_variants(root: Path, emitter: Emitter) -> dict[str, list[dict[str, Any]]]:
    """Emit one shard per chromosome and return each gene's variants.

    The mapping is returned rather than linked because `bundles.py` embeds a
    gene's variants in its bundle, where it links omics rows to a shard. That is
    affordable because every row here is human-curated evidence about congenital
    heart disease, not because variants are bounded per gene: TTN alone carries
    tens of thousands of ClinVar submissions, and a bulk import at that scale
    would make a gene bundle grow without limit and make this the place to link
    rather than embed.
    """
    shards: list[str] = []
    by_gene: dict[str, list[dict[str, Any]]] = {}

    # Every filename is resolved before any file is read, so a misnamed shard
    # fails before this module has written anything. Sorted explicitly rather
    # than leaning on `mirror_paths`, which does sort — but on the filename, in
    # the lexical order this exists to replace.
    shard_paths = [path for path, schema_name in mirror_paths(root) if schema_name == "variants"]
    located = [(_chromosome(path), path) for path in shard_paths]

    for chrom, path in sorted(located, key=lambda pair: _ORDER[pair[0]]):
        frame, _ = read_table(path, TABLE_SCHEMAS["variants"])
        if frame is None:
            # An unreadable mirror is reported by `validate_table` against the
            # same file; failing the build here would report it a second time in
            # a place with less to say about it.
            continue

        # Materialised once: `to_dicts()` builds a fresh list of fresh dicts on
        # every call, so asking again for the per-gene mapping would copy every
        # variant in the atlas a second time. Nothing here mutates a row, so the
        # shard and the bundles can share them.
        rows = frame.to_dicts()
        _refuse_foreign_rows(path, chrom, rows)

        relative = f"variants/{chrom}.json.gz"
        emitter.write_json_gz(relative, {"chrom": chrom, "rows": rows})
        shards.append(relative)

        for row in rows:
            # `gene` is not nullable, so a null here is a curation error that
            # `validate_table` reports as TBL003 and `build_site` refuses to
            # build on. Skipped rather than coerced because `str(None)` is the
            # key "None", which would collect every such row under one fictional
            # gene; the variant is still published in its shard either way.
            gene = row.get("gene")
            if gene is not None:
                by_gene.setdefault(str(gene), []).append(row)

    # Always emitted, even empty: a client that has to probe for 404s to discover
    # which chromosomes exist is a client that breaks the first time one is added.
    emitter.write_json("variants/index.json", {"shards": shards})
    return by_gene
