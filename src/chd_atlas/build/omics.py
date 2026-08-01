# src/chd_atlas/build/omics.py
"""Per-dataset omics shards, plus the bounded summary each gene bundle embeds.

A gene bundle carries counts and a top slice; the full rows live in a shard the
detail page fetches only if someone opens the omics tab. The alternative — one
self-contained bundle per gene — grows without bound as datasets accumulate, and
a visitor reading only the assertion would still pay for every proteomics row.

Three of the four dataset-linked tables name their gene in a column. `phospho`
does not: its rows are keyed on a UniProt accession, and `mirrors/genes.tsv` is
the only thing that maps one back to a gene. `proteomics.gene` is nullable and
so falls to the same join whenever a row is identified only by accession. That
join is the whole risk in this module: every way it can silently return nothing
ends with a gene's evidence absent from its bundle while sitting in a shard that
nothing links to.
"""

from __future__ import annotations

from collections.abc import Mapping
from functools import partial
from pathlib import Path
from typing import Any, Final

from chd_atlas.build.emit import Emitter
from chd_atlas.build.paths import slug
from chd_atlas.tables import TABLE_SCHEMAS, mirror_paths, read_table

# How many rows of each modality a gene bundle embeds. Enough to render a
# meaningful preview without the bundle tracking dataset volume.
TOP_N: Final = 25

# The dataset-linked tables, and the column that names the gene in each — `None`
# where the table has no gene column at all, which is true only of phospho.
_GENE_COLUMN: Final[dict[str, str | None]] = {
    "expression": "gene",
    "profiles": "gene",
    "proteomics": "gene",
    "phospho": None,
}

# Ranking for the embedded slice: FDR ascending, so the most significant row
# leads. `profiles` reports no FDR — it is an abundance table, not a contrast —
# so its slice is the head of the table's own sort order instead.
_RANK_BY: Final = "fdr"


def _genes_by_accession(root: Path) -> dict[str, tuple[str, ...]]:
    """Map each UniProt accession in the gene registry to the genes declaring it.

    Several genes rather than one, because neither the registry nor biology
    constrains it to one: `uniprot` is not part of the table's sort key, so
    nothing checks it for duplicates, and a single UniProt entry covers every
    gene encoding an identical protein — the histone clusters (H4C1..H4C16, all
    P62805) are the standard case. Keeping only one of them would drop the
    others' phospho evidence with nothing reporting it.

    Sorted, because this tuple decides the order genes are summarised in and the
    registry's row order should not. Deduplicated through a list rather than a
    `set` for the same reason: a cluster holds at most a handful of genes, and a
    set would make the pre-sort order depend on `PYTHONHASHSEED` — leaving a
    guard that a regression passes on some fraction of runs instead of one that
    fails on every run. No other order in this module comes out of a hash
    container.
    """
    path = root / "mirrors" / "genes.tsv"
    if not path.is_file():
        return {}
    frame, _ = read_table(path, TABLE_SCHEMAS["genes"])
    if frame is None or not {"uniprot", "hgnc_id"} <= set(frame.columns):
        return {}
    index: dict[str, list[str]] = {}
    for accession, gene in frame.select(["uniprot", "hgnc_id"]).rows():
        if accession is None or gene is None:
            continue
        genes = index.setdefault(str(accession), [])
        if str(gene) not in genes:
            genes.append(str(gene))
    return {accession: tuple(sorted(genes)) for accession, genes in index.items()}


def _genes_for_row(
    row: Mapping[str, Any], column: str | None, by_accession: Mapping[str, tuple[str, ...]]
) -> tuple[str, ...]:
    """Which genes one omics row is evidence about.

    The gene column when the table has one and the cell is filled; the protein
    join otherwise. Reached by more than phospho: `proteomics.gene` is nullable,
    and a proteomics row carrying only an accession is exactly as invisible to
    its gene as a phospho row would be.

    The accession is required to be a string rather than coerced with `str()`:
    coercing turns a null cell into the literal "None", which then matches any
    registry row whose `uniprot` reads "None" and attributes a site to a gene
    that has nothing to do with it. `validate_table` rejects such a registry
    value against `UNIPROT_PATTERN`, but the build reads the mirror directly and
    does not re-run validation.
    """
    if column is not None:
        gene = row.get(column)
        if gene is not None:
            return (str(gene),)
    accession = row.get("protein")
    if not isinstance(accession, str):
        return ()
    exact = by_accession.get(accession)
    if exact is not None:
        return exact
    # A quantified site belongs to an isoform ("Q99593-2") more often than to the
    # canonical accession, while the registry records the gene's canonical one.
    # `UNIPROT_PATTERN` admits both forms, so without this the whole of a gene's
    # phospho evidence can be dropped for a suffix. Tried second so that a
    # registry entry that *is* isoform-specific still wins where it exists.
    return by_accession.get(accession.split("-", 1)[0], ())


def _comparable(value: object) -> tuple[int, float, str]:
    """One totally ordered form of a cell, so a mixed-type key cannot raise.

    Numbers compare as numbers and strings as strings, with nulls first — the
    same order `validate/sort_order.py` enforces on the table itself. Passing
    `str(value)` instead is the shortcut that does not work: it sorts a
    `position` of 100 before one of 9, so the tie-break would not reproduce the
    canonical order it exists to borrow, and *which* rows survive the truncation
    would follow the decimal spelling of an integer.

    Columns are single-typed, so the three tiers never interleave in practice;
    they are what makes the key total rather than a case anything relies on.
    """
    if value is None:
        return (0, 0.0, "")
    if isinstance(value, int | float):
        return (1, float(value), "")
    return (2, 0.0, str(value))


def _rank(row: Mapping[str, Any], sort_key: tuple[str, ...]) -> tuple[Any, ...]:
    """Sort key for the embedded slice: significance first, canonical order after.

    A missing or null FDR ranks last, since a row reporting no FDR is not
    evidence of significance. For `profiles`, which has no such column, that
    applies to every row and the canonical key alone decides the slice.
    """
    fdr = row.get(_RANK_BY)
    significance = float("inf") if fdr is None else float(fdr)
    return (significance, *(_comparable(row.get(field)) for field in sort_key))


def build_omics(root: Path, emitter: Emitter) -> dict[str, dict[str, Any]]:
    """Emit one shard per omics table and return a per-gene, per-modality summary.

    The returned mapping is `{hgnc_id: {modality: {count, shards, top}}}`.
    `count` is every row about the gene, `shards` the files holding them, and
    `top` at most `TOP_N` of the rows themselves.
    """
    by_accession = _genes_by_accession(root)
    summaries: dict[str, dict[str, Any]] = {}

    for path, schema_name in mirror_paths(root):
        if schema_name not in _GENE_COLUMN:
            continue
        frame, _ = read_table(path, TABLE_SCHEMAS[schema_name])
        if frame is None:
            # An unreadable mirror is reported by `validate_table` against the
            # same file; failing the build here would report it a second time in
            # a place with less to say about it.
            continue

        # Materialised once: `to_dicts()` builds a fresh list of fresh dicts on
        # every call, so asking again for the summary would copy every row in the
        # atlas a second time. Nothing here mutates a row, so the shard and the
        # summary can share them.
        rows = frame.to_dicts()

        # Through `slug` for the same reason a gene bundle path is: the stem is a
        # filename, and a space or a colon in one would be published as a URL
        # that has to be escaped before it can be fetched. Every accession
        # grammar in `identifiers.py` passes through unchanged, so this only
        # fires on a filename that should not have been curated. Two stems that
        # slug alike — or that differ only in case — collide, and `Emitter`
        # raises rather than let one shard overwrite the other.
        relative = f"omics/{schema_name}/{slug(path.stem)}.json"
        emitter.write_json(relative, {"table": schema_name, "rows": rows})

        for row in rows:
            for gene in _genes_for_row(row, _GENE_COLUMN[schema_name], by_accession):
                modality = summaries.setdefault(gene, {}).setdefault(
                    schema_name, {"count": 0, "shards": [], "top": []}
                )
                modality["count"] += 1
                if relative not in modality["shards"]:
                    modality["shards"].append(relative)
                modality["top"].append(row)

    # Rank and truncate once, after every shard has contributed, so the slice is
    # the best rows across all datasets rather than the best of whichever file
    # happened to be read last.
    for modalities in summaries.values():
        for schema_name, modality in modalities.items():
            modality["top"].sort(key=partial(_rank, sort_key=TABLE_SCHEMAS[schema_name].sort_key))
            modality["top"] = modality["top"][:TOP_N]
            # `shards` is published as a JSON array, whose order `encode_json`
            # does not touch, and it is accumulated from a directory listing —
            # the case emit.py names as the likeliest way this build loses
            # determinism. `mirror_paths` sorts, but on the filename rather than
            # on the slug of it, and the two need not agree.
            modality["shards"].sort()
    return summaries
