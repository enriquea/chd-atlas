# src/chd_atlas/build/search.py
"""The search index: one flat record per gene, publication and phenotype.

Deliberately a list of records rather than an inverted index. At this corpus size
a client filters the whole array in a fraction of a frame, while an inverted
index is a structure to maintain — kept in step with every payload it summarises,
and debugged when a query misses — before there is the query volume that would
need one.

Variants and datasets are excluded. The variant space is the one that grows
without bound, so indexing it would let the least likely query dominate the size
of the file every visitor downloads before they have typed anything.

Genes are drawn from the assertion set, so a gene carrying omics or variant
evidence but no curated assertion is not searchable. That is `gene_facts`' rule —
the atlas browses curated claims — and matching it here is more than consistency:
`path` advertises `genes/<id>.json`, `build_genes` writes one bundle per gene
`gene_facts` returns, and that is exactly the set below. Widening this to the
gene registry would publish search hits whose only action is a fetch that 404s.
Amendment A29 records the sharper case, which this module does not change: a
`FunctionalEvidence` record about a gene with no assertion reaches no published
file at all, so nothing here can make it findable.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence

from chd_atlas.build.emit import Emitter, Json
from chd_atlas.build.paths import gene_bundle_path
from chd_atlas.corpus import Corpus


def _terms(values: Iterable[str]) -> list[str]:
    """The strings a visitor might type, each kept once, in the order given.

    `dict.fromkeys` rather than a `set`: this list is published, so its order is
    part of the file that every visitor downloads, and set order is neither
    stable across interpreters nor meaningful to a reader. Insertion order is a
    language guarantee, so keeping the first occurrence of each term costs no
    determinism.

    Duplicates arrive on a green gate rather than only through a bypassed one.
    `mirrors/genes.tsv` is free to list a symbol among that gene's own aliases,
    an ontology export routinely repeats a term's label as one of its synonyms,
    and nothing forbids a publication naming one author twice. None of those
    changes a query's outcome — a record matching a term twice is one hit — so
    this buys bytes and a file that reads as intended, not correctness.
    """
    return list(dict.fromkeys(values))


def build_search(
    corpus: Corpus,
    emitter: Emitter,
    symbols: Mapping[str, str],
    aliases: Mapping[str, Sequence[str]],
) -> None:
    """Emit `search/index.json.gz`, the whole searchable surface in one file.

    Each record carries `terms`, the strings a visitor might type, and `path`,
    the payload that answers the query once they pick the result. `label` is what
    the result row displays; `id` is what identifies the thing it names.

    `symbols` and `aliases` come from `mirrors/genes.tsv`, keyed on HGNC id, and
    are taken as `Mapping` because nothing here mutates them — the concrete types
    are invariant, so spelling them `dict` would reject legitimate arguments for
    no benefit. `build_genes` takes `symbols` the same way for the same reason.

    A gene absent from `symbols` keeps its HGNC id as its label, so a result row
    reads "HGNC:4173" rather than blank. Behind `build_site` that fallback is
    unreachable: REF001 makes an asserted gene missing from the registry an
    error, TBL008 makes a missing registry one, and the build refuses on any
    error — so reaching it means the gate was bypassed, and an ugly row a reader
    can still search for beats an unlabelled one.

    Aliases are sorted; authors and synonyms are not. An alias cell is a
    pipe-separated dump whose order is whatever the upstream release emitted, so
    sorting it makes two mirrors listing the same aliases publish identical
    bytes. An author list is data — first author, last author — and a synonym
    list is the curator's, so reordering either would publish something the
    record does not say.

    The two literal filenames below are written by `literature.py`, and this is
    the only place that dependency exists; renaming either payload there means
    renaming it here, and nothing checks the pair. Gene paths carry no such risk
    because `gene_bundle_path` is what `build_genes` writes with.
    """
    records: list[Json] = []

    # Sorted so that the published array does not follow set iteration order,
    # which varies with `PYTHONHASHSEED` between two builds of one commit and
    # would move every checksum downstream of it. The set itself is what makes a
    # gene with several assertions one result rather than one per assertion.
    for gene in sorted({assertion.gene for assertion in corpus.assertions}):
        symbol = symbols.get(gene, gene)
        records.append(
            {
                "kind": "gene",
                "id": gene,
                "label": symbol,
                "terms": _terms([symbol, gene, *sorted(aliases.get(gene, ()))]),
                # `HgncId` reaches this untouched, because `assertion.gene` is
                # annotated as one — unlike `gene_facts`' keys, which are plain
                # `str` and make `bundles.py` recover the NewType before it can
                # name a file. Nothing to assert here, so nothing to get wrong.
                "path": gene_bundle_path(gene),
            }
        )

    # These two arrays are sorted for determinism only. They are already ordered
    # by identifier in `publications.json` and `phenotypes.json`, and a search
    # result is ranked by whatever the client matched on rather than by the order
    # it read them in.
    for publication in sorted(corpus.publications, key=lambda item: item.id):
        records.append(
            {
                "kind": "publication",
                "id": publication.id,
                "label": publication.title,
                "terms": _terms([publication.title, publication.journal, *publication.authors]),
                "path": "publications.json",
            }
        )

    for term in sorted(corpus.phenotypes, key=lambda item: item.id):
        records.append(
            {
                "kind": "phenotype",
                "id": term.id,
                "label": term.label,
                "terms": _terms([term.label, *term.synonyms]),
                "path": "phenotypes.json",
            }
        )

    # Gzipped because this is the one payload downloaded in full by a visitor who
    # has not asked for anything yet, and it is almost entirely repeated English.
    # Emitted even when empty, so a search box fetching it before the first
    # assertion is curated reads an empty array rather than having to tell "no
    # records" from "wrong URL" by interpreting a 404.
    emitter.write_json_gz("search/index.json.gz", {"records": records})
