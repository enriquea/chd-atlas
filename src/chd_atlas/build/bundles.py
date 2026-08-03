# src/chd_atlas/build/bundles.py
"""The gene index and the per-gene bundles.

A bundle is what replaces a backend: a gene detail page is one fetch of one
small file, joined here so that no client ever has to assemble a gene's
assertions, functional records, variants and omics summaries itself.

`genes/index.json` is the other half: the browse payload, downloaded by every
visitor before they have picked a gene. It carries what the browser ranks and
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
URL a detail page is fetched from, and nothing downstream can catch a wrong one:
`Emitter.checksums` records what the build wrote, never what a payload promised,
so a dead link is published with every checksum verifying and every gate green.
One `gene_bundle_path` call per gene supplies both the index entry and the
write, which is what makes the two agree by construction rather than by review.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import Any, cast

from chd_atlas.build.derive import GeneFacts, gene_facts
from chd_atlas.build.emit import Emitter, Json
from chd_atlas.build.omics import ModalitySummary
from chd_atlas.build.paths import gene_bundle_path
from chd_atlas.build.validity import GeneValidity, ValidityRecord, uncurated
from chd_atlas.corpus import Corpus
from chd_atlas.identifiers import HgncId
from chd_atlas.models.assertion import LesionAssertion
from chd_atlas.models.functional import FunctionalEvidence


def _records_by_gene[RecordT: (LesionAssertion, FunctionalEvidence)](
    records: Iterable[RecordT],
) -> dict[str, list[dict[str, Json]]]:
    """Group one kind of curated record by gene, ordered by id, JSON-ready.

    Ordered by id rather than left in corpus order, because these are published
    as JSON arrays and `encode_json`'s `sort_keys` orders dict keys and nothing
    else: whatever order they reach the emitter in is the order a consumer
    downloads. Corpus order is `load_curation`'s: record files sorted by
    filename — `_record_files` sorts the glob — and then each file's own order
    within it, which makes it an artefact of how the curation is filed rather
    than of the records. Without this, moving a record from `TBX5.yaml` into a
    second file would reorder a bundle that is otherwise unchanged.

    `mode="json"` so every value is a JSON primitive by construction.
    `LesionAssertion` carries two `date` fields, which `json.dumps` refuses
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


def _validity_record(record: ValidityRecord) -> dict[str, Json]:
    """One mirrored classification, JSON-ready.

    Every record in the published array carries the same eleven keys, whichever
    source it came from. `sop`, `classification_date` and `gcep` are ClinGen-only
    and `submitter` is GenCC-only on `ValidityRecord` itself (see its
    docstring), and each is `null` on a record from the other source rather than
    left out of the object. The alternative — omitting a source's absent fields
    — would publish an array whose objects differ in shape, which is a trap for
    any consumer that reads one field off `records[0]` and expects the same key
    to exist on `records[1]`.

    `classification` publishes `null` for GenCC's `Supportive`, which
    `ValidityRecord.classification` already represents as `None` — the mapping
    exception documented on that field — while `classification_term` always
    carries the authority's own word, mapped or not.
    """
    return {
        "source": record.source.value,
        "classification_term": record.classification_term,
        "classification": record.classification.value if record.classification else None,
        "disease": record.disease,
        "disease_label": record.disease_label,
        "moi": record.moi,
        "sop": record.sop,
        "classification_date": record.classification_date,
        "gcep": record.gcep,
        "report_url": record.report_url,
        "submitter": record.submitter,
    }


def _validity(gene_validity: GeneValidity) -> dict[str, Json]:
    """The bundle's `validity` block: how curated the gene is, and by whom.

    `records` is serialised in the order `build.validity.gene_validity` already
    sorted it into — content order, not encounter order, so two builds of one
    commit agree on it. Never re-sorted or re-derived here: that module already
    did the work of deciding the order, and doing it again risks disagreeing
    with it instead of reusing it.

    `state` and `has_source_discordance` duplicate the top-level
    `validity_state` and `has_source_discordance` that `_headline` already
    writes into this same bundle. That is deliberate rather than an oversight:
    the top-level fields are what the browse row and the bundle publish
    identically (`_headline`'s whole reason to exist), while this block is the
    self-contained provenance object — a consumer who only wants "who curated
    this gene, and what did each of them say" reads `validity` alone and never
    has to cross-reference the flat fields beside it.
    """
    return {
        "state": gene_validity.state.value,
        "has_source_discordance": gene_validity.has_source_discordance,
        "records": [_validity_record(record) for record in gene_validity.records],
    }


def _headline(gene: str, symbol: str, fact: GeneFacts) -> dict[str, Json]:
    """What the index row and the bundle must say identically about one gene.

    One function rather than two literals so the browse row and the page it
    opens cannot disagree, and — the reason it exists — so `headline_confidence`,
    `validity_state`, `has_conflicting_evidence` and `has_source_discordance` are
    written together. Dropping any one of them from one payload alone is not an
    edit that can be made here by accident.

    `.value` publishes the vocabulary's string rather than leaning on a `StrEnum`
    member *being* a `str`. Measured, not assumed: mypy accepts either spelling
    inside a payload this size, and `json.dumps` writes the same bytes for both
    today. What differs is a vocabulary that stops being a `StrEnum` — the bare
    member would then fail type checking here and be refused by `json.dumps`,
    while `.value` keeps publishing the same JSON. Neither is a trap; `.value` is
    the one that states at the call site what the file will contain.

    `headline_confidence` publishes `null` for a gene no authority has assessed
    — `fact.headline_confidence` is `None` in that case, per `derive.gene_facts`
    — rather than coercing it to a string, which would invent a confidence
    nobody stated.
    """
    return {
        "gene": gene,
        "symbol": symbol,
        "headline_confidence": (
            fact.headline_confidence.value if fact.headline_confidence else None
        ),
        "validity_state": fact.validity_state.value,
        "has_conflicting_evidence": fact.has_conflicting_evidence,
        "has_source_discordance": fact.has_source_discordance,
        "lesion_groups": [group.value for group in fact.lesion_groups],
    }


def build_genes(
    corpus: Corpus,
    emitter: Emitter,
    symbols: Mapping[str, str],
    omics: Mapping[str, Mapping[str, ModalitySummary]],
    variants: Mapping[str, list[dict[str, Any]]],
    validity: dict[str, GeneValidity],
) -> None:
    """Emit `genes/index.json` and one bundle per gene carrying an assertion.

    `symbols` comes from `mirrors/genes.tsv`, keyed on HGNC id. A gene absent
    from it keeps its HGNC id as its display label: the browse row and the page
    heading both render this value, so blank reads as a broken site while
    "HGNC:4173" reads as a symbol that has not been mirrored yet and is still
    something a reader can search for.

    `validity` is passed straight to `gene_facts` — see its docstring for what
    a gene absent from it publishes. The bundle additionally embeds the whole
    `GeneValidity` as a `validity` object (`_validity`), which `gene_facts`
    does not carry forward because the browse layer never needs the individual
    mirrored records, only the fields it already reduces `GeneValidity` to.

    `omics` and `variants` are what `build_omics` and `build_variants` returned,
    taken as `Mapping` because nothing here mutates them. Genes they carry that
    hold no assertion are ignored, which is the same rule `gene_facts` applies:
    the atlas browses curated claims, and a gene with nothing asserted has no
    confidence to display. Their rows are still published in the shards those two
    modules wrote — the gene index simply does not link to them.

    That last sentence covers `omics` and `variants` and nothing else. A
    `FunctionalEvidence` record about a gene with no assertion is worse off: no
    other build module reads `corpus.functional` — `derive.py` reads it only to
    count records per gene, and this is its only writer — so it reaches no
    published file at all rather than an unlinked one. That is evidence loss,
    not a missing link. It still does not raise here, because such a record is
    legal under every validator in the project and `gene_facts` is what decides
    which genes exist; the cost is pinned instead, by
    `test_a_gene_with_evidence_but_no_assertion_is_not_published`.
    """
    facts = gene_facts(corpus, validity)
    assertions = _records_by_gene(corpus.assertions)
    functional = _records_by_gene(corpus.functional)

    index: list[Json] = []
    # Sorted here as well as in `gene_facts`, which returns its genes in this
    # same order — pinned by `test_derive.py`'s
    # `test_assertions_are_counted_per_gene_and_the_genes_come_back_sorted`.
    # Redundant today, and kept anyway: this array is what consumers download,
    # so its order is the published contract and should not rest on a guarantee
    # another module is free to withdraw.
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
        # Looked up again rather than threaded through `GeneFacts`: `gene_facts`
        # already reduced `GeneValidity` to the scalar fields the browse layer
        # needs (`headline_confidence`, `validity_state`, ...) and does not carry
        # the record list forward, so the bundle -- the one payload that needs
        # the records themselves -- reads the same `validity` mapping this
        # function was handed, with the same `uncurated()` fallback `gene_facts`
        # documents for a gene neither mirror curates.
        gene_validity = validity.get(gene, uncurated())

        index.append(
            {
                **headline,
                "confidence_by_lesion_group": {
                    group.value: classification.value
                    for group, classification in fact.confidence_by_lesion_group.items()
                },
                # Written immediately beside the map it qualifies, for the reason
                # `_headline` exists: a confidence a reader can act on must never
                # be published without the flag that says it is contested. At
                # gene level one function guarantees that; here it is adjacency
                # plus the test that reads both, because the two are different
                # shapes and cannot share a writer.
                "conflicting_lesion_groups": [
                    group.value for group in fact.conflicting_lesion_groups
                ],
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
                # The full mirrored picture, alongside the flat fields
                # `headline` already carries -- see `_validity`'s docstring for
                # why both are published rather than one replacing the other.
                "validity": _validity(gene_validity),
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
