# src/chd_atlas/build/derive.py
"""Per-gene facts derived from the corpus, never written back to it.

Spec 5.2: the source of truth stays normalized and the API is denormalized for
consumption. Nothing here touches the filesystem, so it is cheap to test and
cannot disagree with what validation saw — it reads the same in-memory `Corpus`
rather than re-reading the files behind it.

Every collection returned is ordered by content rather than by encounter, which
is what lets the build be byte-identical. The two tuples are the ones that
matter: whatever emits them writes them as JSON *arrays*, and `encode_json`'s
`sort_keys` orders dict keys and nothing else, so an array's order here is the
order a consumer downloads. `publications` in particular is accumulated through
a `set`, whose iteration order for strings varies with `PYTHONHASHSEED` —
unsorted, it checksums differently from one process to the next, which no local
re-run reproduces. The two dicts are ordered for readability only; `sort_keys`
re-orders those on the way out regardless.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Collection
from dataclasses import dataclass

from chd_atlas.build.validity import GeneValidity, uncurated
from chd_atlas.corpus import Corpus
from chd_atlas.models.assertion import LesionAssertion
from chd_atlas.vocab import (
    AtlasCuration,
    Classification,
    EvidenceClass,
    LesionGroup,
    ValiditySource,
    ValidityState,
    has_conflicting_evidence,
    strongest,
)


@dataclass(frozen=True)
class GeneFacts:
    """What the browse layer needs about one gene without opening its bundle."""

    gene: str
    # `None` for a gene no authority has assessed -- deliberately not
    # `Classification.NO_KNOWN_ASSOCIATION`, which is itself an assessed verdict
    # ("a panel looked and found nothing") and would assert a conclusion nobody
    # reached. Sourced from `build/validity.py`'s mirrored ClinGen/GenCC records,
    # never from a curated assertion: the atlas mirrors validity, it does not
    # author it.
    headline_confidence: Classification | None
    # How well curated the gene is: an expert panel, only submitters, or nobody.
    # Published as a fact rather than left inferable from `headline_confidence`
    # being `None`, because `None` alone cannot distinguish "no authority has
    # assessed this gene" from a future state where the atlas has some other
    # reason to withhold a headline.
    validity_state: ValidityState
    # Paired with `headline_confidence` at every display site. `strongest` ranks
    # on one linear scale, so a gene carrying both a definitive and a refuted
    # record resolves to definitive and buries the refutation; this is what
    # stops the browse layer showing a contested gene as settled.
    has_conflicting_evidence: bool
    # True iff one mirrored source contests the gene while a *different* source
    # supports it -- `build/validity.py::_has_source_discordance`. A narrower
    # signal than `has_conflicting_evidence`, which also catches one source
    # split against itself; this one says the two authorities disagree with each
    # other, not merely that the evidence is mixed.
    has_source_discordance: bool
    lesion_groups: tuple[LesionGroup, ...]
    # The same pairing as `headline_confidence`/`has_conflicting_evidence`, one
    # level down -- but grouped, not independently derived. ClinGen and GenCC
    # classify a gene against a disease, never against a specific lesion, so the
    # mirrored classifications carry no per-group information of their own: this
    # collapses the *same* `strongest()` of the gene's mirrored classifications
    # under every lesion group the curated assertions name for the gene, which is
    # why every group a contested gene names ends up in `conflicting_lesion_groups`
    # together -- there is no finer-grained signal to divide them with. Empty
    # exactly when `headline_confidence` is `None`, for the same reason.
    confidence_by_lesion_group: dict[LesionGroup, Classification]
    # Names every group `confidence_by_lesion_group` carries when the gene-level
    # `has_conflicting_evidence` is true, and none when it is false. A list
    # beside the map rather than a flag inside it: additive, so a consumer
    # already reading `confidence_by_lesion_group` is unaffected, and it mirrors
    # how `lesion_groups` already sits beside the same map. Issue #4.
    conflicting_lesion_groups: tuple[LesionGroup, ...]
    evidence_counts: dict[EvidenceClass, int]
    assertion_count: int
    # Whether the atlas has curated evidence of its own for this gene, as
    # opposed to republishing an expert panel's classification. True of 1 of the
    # 92 genes published today. See `vocab.AtlasCuration`.
    atlas_curation: AtlasCuration
    # Every functional record about the gene, not only those an assertion cites.
    # `referential.py` requires a *cited* record to be about the asserting gene
    # (REF008) but never requires a curated record to be cited at all, so the two
    # populations genuinely differ. This counts the same set the gene bundle is
    # specified to list — every record about the gene — because a browse-level
    # count smaller than the page it links to reads as missing data.
    functional_count: int
    publications: tuple[str, ...]


def gene_facts(
    corpus: Corpus,
    validity: dict[str, GeneValidity],
    published: Collection[str],
) -> dict[str, GeneFacts]:
    """Derive one `GeneFacts` per published gene.

    Keyed on HGNC id. `published` is `build.validity.published_genes()`'s return
    -- design decision D21: the atlas publishes a gene when a ClinGen expert
    panel classifies it `Limited` or better for an in-scope disease, or when two
    or more GenCC submitters independently assert it and no ClinGen panel
    contests it.

    That population is neither the asserted genes nor the registry. This
    function used to key on `{assertion.gene for assertion in corpus.assertions}`
    and justify it with "a gene with nothing asserted has no confidence to
    display". Since the validity backbone landed that sentence has been false:
    every field below that carries a confidence -- `headline_confidence`,
    `validity_state`, `has_conflicting_evidence`, `has_source_discordance` and
    `confidence_by_lesion_group` -- comes from `validity`, and none of them from
    a curated assertion. A gene an expert panel has classified has exactly a
    confidence to display, and 91 of the 92 genes published today have no
    assertion at all.

    An asserted gene outside `published` gets no facts, so no bundle and no
    page. That is D37: a curator's assertion does not admit a gene to the
    published set, and candidate genes belong in their own labelled section. The
    2026-08-06 widening does not soften this -- it moves the authority from one
    external source to two, never to the atlas.

    `validity` is `build.validity.gene_validity()`'s return. A gene absent from
    it gets `uncurated()`, which cannot happen for a member of `published` --
    `published_genes` derives from the same mapping -- and is kept because a
    caller may pass a wider population.

    The three collection fields that come from assertions (`lesion_groups`,
    `confidence_by_lesion_group`, `evidence_counts`) are empty for an uncurated
    gene, which is the honest answer: the atlas has recorded no lesion for it.
    """
    by_gene: dict[str, list[LesionAssertion]] = {}
    for assertion in corpus.assertions:
        by_gene.setdefault(assertion.gene, []).append(assertion)

    # Annotated `str` rather than left to inference: `record.gene` is `HgncId`, a
    # NewType over `str`, so an inferred `Counter[HgncId]` cannot be indexed by
    # the plain `str` keys of `by_gene` under --strict. The two are the same
    # value at runtime; this keeps the key type consistent with the returned
    # `dict[str, GeneFacts]` rather than pushing the NewType through the API.
    functional_counts: Counter[str] = Counter(record.gene for record in corpus.functional)

    facts: dict[str, GeneFacts] = {}
    # `sorted` rather than `published`'s own iteration order, which makes the
    # returned mapping's key order a contract this function keeps rather than an
    # accident of the `set[str]` `published_genes` hands it. That matters because
    # the return is public: `build_genes` hands the same mapping to
    # `build_gene_pages` and `build_gene_index_page`, and a caller is entitled to
    # iterate it without sorting first.
    #
    # It is **not** what keeps the build byte-identical, and this comment claimed
    # it was until the claim was measured. Every consumer that turns this mapping
    # into published bytes re-sorts: `bundles.py::build_genes`, `build_gene_pages`
    # and `build_gene_index_page` all iterate `sorted(facts)`, and the browse
    # page's facets go through sets that `pages.py` sorts itself. Measured
    # 2026-08-04, one process per seed with this `sorted` dropped: full builds
    # under `PYTHONHASHSEED` 0, 2, 12345 and 99 were byte-identical to each other
    # and to a sorted build, across all 57 emitted files. Keep the sort for the
    # contract; do not cite it as the determinism guard, or removing a caller's
    # `sorted` will look safe.
    #
    # Pinned by `test_assertions_are_counted_per_gene_and_the_genes_come_back_sorted`,
    # which hands in a reverse-ordered `list` rather than a `set` so the mutant
    # fails on every run instead of on a fraction of the seeds.
    for gene in sorted(published):
        assertions = by_gene.get(gene, [])
        gene_validity = validity.get(gene, uncurated())
        # Only a record with a mapped `classification` takes a side --
        # `Classification | None` on `ValidityRecord`, `None` where GenCC's
        # `Supportive` declines to grade the association at all.
        mirrored = [
            record.classification
            for record in gene_validity.records
            if record.classification is not None
        ]

        groups: set[LesionGroup] = set()
        counts: Counter[EvidenceClass] = Counter()
        publications: set[str] = set()
        for assertion in assertions:
            groups.update(assertion.lesion_groups)
            for evidence in assertion.evidence:
                counts[evidence.evidence_class] += 1
                publications.add(evidence.publication)

        # Sorted once and reused, so the array and the breakdown cannot come to
        # disagree about the order of the same set of groups.
        ordered_groups = sorted(groups, key=lambda group: group.value)

        # **The headline is the ADMITTING EXPERT PANEL's grade, or nothing.**
        # It was `strongest()` over every mirrored classification -- ClinGen's
        # and GenCC's together -- until 2026-08-06. That was harmless while the
        # gate required a ClinGen `Definitive`, because a chip could not appear
        # without a chartered panel behind it. Widening the gate to admit genes
        # on GenCC agreement broke that: measured on the widened corpus, **24 of
        # 93 genes** would carry a chip stronger than any ClinGen grade, and
        # five -- ELN, GDF1, MMP21, PKD1L1, TBX1 -- would show a green
        # `definitive` chip with no ClinGen record at all.
        #
        # So a gene no expert panel graded publishes `None`, exactly as an
        # uncurated gene does and for the same stated reason: coercing it to a
        # string would invent a confidence nobody stated. What the submitters
        # said is not lost -- it is in `validity.records` and in
        # `admitted_by.submitters` -- it is simply not the atlas's headline,
        # because there is no single authority to headline. GDF1 is why: its
        # in-scope submissions run from G2P `Definitive` to Illumina `No Known
        # Disease Relationship`, and a max over submitters publishes that as
        # settled.
        #
        # `validity_state` sits beside this and says `submitter_curated`, which
        # is what distinguishes "no panel has graded this" from "nobody has
        # assessed this gene at all".
        panel_graded = [
            record.classification
            for record in gene_validity.records
            if record.source is ValiditySource.CLINGEN and record.classification is not None
        ]
        headline = strongest(panel_graded) if panel_graded else None
        contested = has_conflicting_evidence(mirrored)

        facts[gene] = GeneFacts(
            gene=gene,
            headline_confidence=headline,
            validity_state=gene_validity.state,
            has_conflicting_evidence=contested,
            has_source_discordance=gene_validity.has_source_discordance,
            lesion_groups=tuple(ordered_groups),
            confidence_by_lesion_group=(
                {group: headline for group in ordered_groups} if headline is not None else {}
            ),
            # The mirrors classify per gene-disease, not per lesion group, so
            # every group the curated assertions name shares the one gene-wide
            # answer: all of `ordered_groups` when the gene is contested, none
            # when it is not.
            conflicting_lesion_groups=tuple(ordered_groups) if contested else (),
            evidence_counts={
                item: counts[item] for item in sorted(counts, key=lambda item: item.value)
            },
            assertion_count=len(assertions),
            atlas_curation=(AtlasCuration.CURATED if assertions else AtlasCuration.NOT_YET_CURATED),
            functional_count=functional_counts[gene],
            publications=tuple(sorted(publications)),
        )
    return facts
