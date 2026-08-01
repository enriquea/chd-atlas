# src/chd_atlas/build/bundles.py
"""The gene index and the per-gene bundles.

A bundle is what replaces a backend: a gene detail page is one fetch of one
small file, joined here so that no client ever has to assemble a gene's
assertions, functional records, variants and omics summaries itself.

`genes/index.json` is the other half, and the only payload every visitor
downloads before they have picked a gene. It carries what the browser ranks and
filters on plus the path to fetch the rest, and none of the evidence itself —
that is what keeps it small as curation grows, and why the bundle is a separate
file rather than a section of this one.

Two rules decide what goes where.

A contested gene must never be displayed as settled. `strongest()` ranks
classifications on one linear scale where `definitive` outranks `refuted`, while
ClinGen treats disputed and refuted as a separate axis rather than weaker rungs
of the same ladder, so `headline_confidence` on its own hides a refutation
behind a stronger assertion about the same gene. `has_conflicting_evidence` is
the other half of that pair; `_headline` below is the single place either is
written, so a payload cannot come to carry the confidence without the flag.

Every path the index advertises must be a file that was written. `bundle` is the
URL a detail page is fetched from and nothing downstream checks it — the
manifest records what was written, not what was promised, so a dead link would
pass every existing gate. One `gene_bundle_path` call per gene supplies both the
index entry and the write, which is what makes the two agree by construction.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import Any, cast

from chd_atlas.build.derive import GeneFacts, gene_facts
from chd_atlas.build.emit import Emitter, Json
from chd_atlas.build.omics import ModalitySummary
from chd_atlas.build.paths import gene_bundle_path
from chd_atlas.corpus import Corpus
from chd_atlas.identifiers import HgncId
from chd_atlas.models.assertion import GeneDiseaseAssertion
from chd_atlas.models.functional import FunctionalEvidence


def _records_by_gene[RecordT: (GeneDiseaseAssertion, FunctionalEvidence)](
    records: Iterable[RecordT],
) -> dict[str, list[dict[str, Json]]]:
    """Group one kind of curated record by gene, ordered by id, JSON-ready.

    Ordered by id rather than left in corpus order, because these are published
    as JSON arrays and `encode_json`'s `sort_keys` orders dict keys and nothing
    else: whatever order they reach the emitter in is the order a consumer
    downloads. Corpus order is the order `load_curation` globbed the files in,
    so without this a record moved from `TBX5.yaml` into a second file would
    reorder a bundle that is otherwise unchanged.

    `mode="json"` so every value is a JSON primitive by construction.
    `GeneDiseaseAssertion` carries two `date` fields, which `json.dumps` refuses
    outright — without this the build dies part-way through, with `dist/` already
    half written. `FunctionalEvidence` carries none today and would serialise
    either way, which is exactly why the conversion is applied to both rather
    than to the one that needs it.

    The type parameter is constrained to the two record types rather than bound
    to `BaseModel`, because the body reads `.gene` and `.id`, which `BaseModel`
    does not declare. Constrained means mypy checks this body once per record
    type, so a third record type has to be added here deliberately.
    """
    grouped: dict[str, list[RecordT]] = {}
    for record in records:
        grouped.setdefault(record.gene, []).append(record)
    return {
        gene: [record.model_dump(mode="json") for record in sorted(items, key=lambda r: r.id)]
        for gene, items in grouped.items()
    }


def _summaries(modalities: Mapping[str, ModalitySummary]) -> Mapping[str, Json]:
    """One gene's omics summaries, retyped rather than rebuilt.

    `ModalitySummary` is a `TypedDict`, and a `TypedDict` is assignable only to
    `Mapping[str, object]` — never to `Mapping[str, Json]`, whatever its fields
    hold. That is a structural fact about `TypedDict`, not a judgement about
    these values, so `emit.Json` cannot accept one and the cast is the whole fix.

    The alternative is to rebuild the summary from the three fields it declares
    today. That type-checks with nothing cast, and it is worse: a field added to
    `ModalitySummary` would be dropped here silently, and the failure would be a
    gene page missing data that the shard it links to still contains. Passing the
    value through instead makes the same change either work or fail loudly, since
    `encode_json` raises on anything `json.dumps` cannot serialise.
    """
    return cast(Mapping[str, Json], modalities)


def _headline(gene: str, symbol: str, fact: GeneFacts) -> dict[str, Json]:
    """What the index row and the bundle must say identically about one gene.

    One function rather than two literals so the browse row and the page it
    opens cannot disagree, and — the reason it exists — so `headline_confidence`
    and `has_conflicting_evidence` are written together. Dropping the flag from
    one payload alone is not an edit that can be made here by accident.

    `.value` on each enum is redundant at runtime, since a `StrEnum` member *is*
    a `str` and `json.dumps` writes it as its value either way. It is here so
    that a vocabulary changed to a plain `Enum` fails type checking rather than
    silently publishing "Classification.DEFINITIVE".
    """
    return {
        "gene": gene,
        "symbol": symbol,
        "headline_confidence": fact.headline_confidence.value,
        "has_conflicting_evidence": fact.has_conflicting_evidence,
        "lesion_groups": [group.value for group in fact.lesion_groups],
    }


def build_genes(
    corpus: Corpus,
    emitter: Emitter,
    symbols: Mapping[str, str],
    omics: Mapping[str, Mapping[str, ModalitySummary]],
    variants: Mapping[str, list[dict[str, Any]]],
) -> None:
    """Emit `genes/index.json` and one bundle per gene carrying an assertion.

    `symbols` comes from `mirrors/genes.tsv`, keyed on HGNC id. A gene absent
    from it keeps its HGNC id as its display label: the browse row and the page
    heading both render this value, so blank reads as a broken site while
    "HGNC:4173" reads as a symbol that has not been mirrored yet and is still
    something a reader can search for.

    `omics` and `variants` are what `build_omics` and `build_variants` returned,
    taken as `Mapping` because nothing here mutates them. Genes they carry that
    hold no assertion are ignored, which is the same rule `gene_facts` applies:
    the atlas browses curated claims, and a gene with nothing asserted has no
    confidence to display. Their rows are still published in the shards those two
    modules wrote — the gene index simply does not link to them.
    """
    facts = gene_facts(corpus)
    assertions = _records_by_gene(corpus.assertions)
    functional = _records_by_gene(corpus.functional)

    index: list[Json] = []
    # Sorted here as well as in `gene_facts`, which returns its genes in this
    # same order. Redundant today, and kept because this array's order is what
    # consumers download while that one is an internal detail of another module
    # that no test pins.
    for gene in sorted(facts):
        fact = facts[gene]
        # `HgncId` rather than `str` at this boundary: `gene_facts` keys on a
        # plain `str` by choice, but every key it returns came from
        # `assertion.gene`, which pydantic validated against `HGNC_PATTERN`. The
        # NewType is being recovered, not asserted. `gene_bundle_path` requires
        # it so that only HGNC ids can name a file in `genes/`, where a
        # `ContrastId` would slug to a name that collides with one on a
        # case-insensitive filesystem.
        bundle = gene_bundle_path(HgncId(gene))
        headline = _headline(gene, symbols.get(gene, gene), fact)
        rows = variants.get(gene, [])

        index.append(
            {
                **headline,
                "confidence_by_lesion_group": {
                    group.value: classification.value
                    for group, classification in fact.confidence_by_lesion_group.items()
                },
                "evidence_counts": {
                    evidence_class.value: count
                    for evidence_class, count in fact.evidence_counts.items()
                },
                "assertion_count": fact.assertion_count,
                "functional_count": fact.functional_count,
                # Counted from the same list the bundle embeds, so the browse
                # promise and the page cannot differ.
                "variant_count": len(rows),
                "bundle": bundle,
            }
        )

        emitter.write_json(
            bundle,
            {
                **headline,
                "publications": list(fact.publications),
                "assertions": assertions.get(gene, []),
                "functional": functional.get(gene, []),
                # Embedded rather than linked. Not because variants are bounded
                # per gene — nothing bounds them, and TTN alone carries tens of
                # thousands of ClinVar records — but because this atlas curates
                # variants by hand as CHD evidence, so the count is a curation
                # policy rather than a property of the data. If that policy
                # changes, this has to become a link the way omics rows are.
                "variants": rows,
                # Always present, empty when the gene has none, so a consumer
                # reads one shape rather than guarding for a missing key.
                "omics": _summaries(omics.get(gene, {})),
            },
        )

    # Written last, so a gene whose bundle could not be written — two genes
    # normalising to one filename is what `Emitter` refuses — fails the build
    # with no index at all, rather than leaving one that advertises a file
    # nothing wrote.
    #
    # Emitted even when there is nothing to list: a consumer that has to tell
    # "no genes curated yet" from "wrong URL" by reading a 404 will get it wrong,
    # and this is the first thing the browse page fetches.
    emitter.write_json("genes/index.json", {"genes": index})
