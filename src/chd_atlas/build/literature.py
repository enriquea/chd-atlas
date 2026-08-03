# src/chd_atlas/build/literature.py
"""The four whole-corpus payloads: publications, featured, phenotypes, datasets.

None of the four is sharded. What makes serving each entire affordable is that
every one of them is bounded by hand curation rather than by an upstream
release: publications, featured entries and the phenotype vocabulary are one
line-reviewed YAML file apiece, and a dataset record is one file per accession
carrying nothing but its own contrast definitions. The measurements those
records interpret are not here — they live in the mirror tables `omics.py` and
`variants.py` shard. Of the four *curated* record types only a dataset nests
another entity, carrying every contrast it declares, so `datasets.json` is the
first that would stop being small if curation scaled. What `_featured` nests is
a second copy of a record `publications.json` already carries in full, bounded
by how many manuscripts a landing page shows.

All four are emitted even when empty: a consumer that has to tell "not curated
yet" from "wrong URL" by reading a 404 will get it wrong, and phenotypes.json is
fetched by any page that renders a lesion facet.
"""

from __future__ import annotations

from pydantic import BaseModel

from chd_atlas.build.emit import Emitter, Json
from chd_atlas.build.paths import DATASETS, FEATURED, PHENOTYPES, PUBLICATIONS, SOURCES
from chd_atlas.corpus import Corpus
from chd_atlas.validate.sources import SourceRegistry


def _dump(model: BaseModel) -> dict[str, Json]:
    """One pydantic model as a JSON-ready dict.

    `mode="json"` so every value is a JSON primitive by construction. Nothing in
    these four models needs the conversion today — their vocabularies are
    `StrEnum`s, which `json.dumps` already writes as their value because a
    `StrEnum` member *is* a `str`, and none of them carries a `date` — so this is
    what keeps that true once one of them gains a field `json.dumps` refuses.

    Typed against `BaseModel` and `Json` rather than `Any`. `emit.Json` exists so
    an unserialisable payload is caught at the call site, and an `Any` in the one
    function every payload here passes through would switch that check off for
    the whole module.
    """
    return model.model_dump(mode="json")


def _featured(corpus: Corpus) -> list[Json]:
    """The featured manuscripts, in display order, each carrying its publication.

    Inline rather than by PMID because the landing page is meant to be one fetch:
    resolving the reference client-side would make first paint wait on a second
    request for a list already known at build time. The key it replaces held that
    PMID and `Publication.id` carries it, so overwriting loses nothing.

    Sorted on `(order, publication)`. `FeaturedFile.display_order_is_unique` makes
    the first element decide, and `curation/featured.yaml` is the only file
    `corpus.featured` is ever loaded from, so that per-file guarantee is the
    corpus-wide one. The PMID is there so a tie — reachable only by building a
    `Corpus` directly — is broken by content rather than by the order two records
    happened to be listed in.
    """
    publications = {publication.id: publication for publication in corpus.publications}

    entries: list[Json] = []
    for entry in sorted(corpus.featured, key=lambda item: (item.order, item.publication)):
        cited = publications.get(entry.publication)
        # Refused rather than published as the bare PMID. That fallback would make
        # `publication` an object for every resolvable entry and a string for this
        # one, so a page reading `entry.publication.title` gets `undefined` with
        # nothing raised — in the payload the landing page renders without a
        # second fetch. `validate_references` reports the same entry as REF002,
        # and in the one case it does not — a corpus that failed to load, where
        # referential checks are skipped — the load failure is itself an error.
        # `build_site` refuses to build on any error, so reaching this means
        # validation was bypassed, not that a curator hit an unusual case.
        if cited is None:
            raise ValueError(
                f"featured entry {entry.order} cites {entry.publication}, which is not in "
                f"curation/publications.yaml; the landing page embeds the publication "
                f"itself and has nothing to render for it"
            )
        payload = _dump(entry)
        payload["publication"] = _dump(cited)
        entries.append(payload)
    return entries


def build_sources(registry: SourceRegistry, emitter: Emitter) -> None:
    """Emit `sources.json`: what the atlas mirrors, and on whose terms.

    A licensing obligation rather than a curated payload. `mirrors/sources.yaml`
    records HPO as `redistribution: permitted_with_attribution`, and the site
    republishes its term labels and synonyms verbatim in `phenotypes.json` and in
    the search index. Attribution therefore has to travel *with the data*: a
    consumer fetching JSON never opens the README, and the repository's only
    `LICENSE` is Apache-2.0, which a reader will reasonably take to cover the
    published files as well.

    The whole registry is published rather than only the sources whose terms
    demand it. Which resources a claim rests on is provenance a scientific reader
    is entitled to, and filtering to the ones with obligations would make the
    file a compliance artifact instead of a statement of where the data came
    from — and would silently drop a source the day someone relaxes its terms.

    Sorted by id so two builds of one commit agree, for the same reason every
    other array here is sorted.
    """
    emitter.write_json(
        SOURCES,
        {"sources": [_dump(source) for source in sorted(registry.sources, key=lambda s: s.id)]},
    )


def build_literature(corpus: Corpus, emitter: Emitter) -> None:
    """Emit publications.json, featured.json, phenotypes.json and datasets.json."""
    # Resolved before anything is written, so an unresolvable citation fails with
    # the output directory untouched rather than after two of the four payloads
    # are already on disk. Same reason `variants.py` resolves every shard filename
    # before it reads a file.
    featured = _featured(corpus)

    # Every array is ordered by identifier, lexically — which puts PMID:10 and
    # GSE100 ahead of PMID:9 and GSE9. Deliberate: a PMID is issued when PubMed
    # indexes a record and an accession when data is deposited, so even in
    # numeric order neither ranks anything a reader asked for — not publication
    # year, not relevance — and nothing displays these arrays in identifier order
    # at all. A page ranks publications by year, or in the curator's featured
    # order. What the sort must do is agree between two builds of one commit —
    # `encode_json` orders dict keys and never list elements — and agree with
    # `derive.py`, which orders a gene's publication list the same way. Making
    # only these numeric would leave the two disagreeing about the same PMIDs.
    #
    # Emitted from `corpus.publications` rather than from the lookup `_featured`
    # keys on: a dict keyed by PMID silently drops the second of two records
    # sharing one. `PublicationFile.ids_are_unique` forbids that pair, so this is
    # about what happens when it is reached anyway — a record published twice is
    # visible to whoever reads the file, one dropped is not.
    emitter.write_json(
        PUBLICATIONS,
        {
            "publications": [
                _dump(publication)
                for publication in sorted(corpus.publications, key=lambda item: item.id)
            ]
        },
    )

    emitter.write_json(FEATURED, {"featured": featured})

    emitter.write_json(
        PHENOTYPES,
        {"phenotypes": [_dump(term) for term in sorted(corpus.phenotypes, key=lambda t: t.id)]},
    )

    emitter.write_json(
        DATASETS,
        {"datasets": [_dump(dataset) for dataset in sorted(corpus.datasets, key=lambda d: d.id)]},
    )
