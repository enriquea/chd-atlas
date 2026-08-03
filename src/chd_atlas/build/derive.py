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
from dataclasses import dataclass

from chd_atlas.build.validity import GeneValidity, uncurated
from chd_atlas.corpus import Corpus
from chd_atlas.models.assertion import LesionAssertion
from chd_atlas.vocab import (
    Classification,
    EvidenceClass,
    LesionGroup,
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
    # Every functional record about the gene, not only those an assertion cites.
    # `referential.py` requires a *cited* record to be about the asserting gene
    # (REF008) but never requires a curated record to be cited at all, so the two
    # populations genuinely differ. This counts the same set the gene bundle is
    # specified to list — every record about the gene — because a browse-level
    # count smaller than the page it links to reads as missing data.
    functional_count: int
    publications: tuple[str, ...]


def gene_facts(corpus: Corpus, validity: dict[str, GeneValidity]) -> dict[str, GeneFacts]:
    """Derive one `GeneFacts` per gene that carries at least one assertion.

    Keyed on HGNC id. Genes present in `mirrors/genes.tsv` but carrying no
    assertion are deliberately absent: the atlas browses curated claims, and a
    gene with nothing asserted has no confidence to display.

    `validity` is `build.validity.gene_validity()`'s return: one `GeneValidity`
    per gene either mirror curates, within CHD scope. `headline_confidence`,
    `has_conflicting_evidence`, `validity_state`, `has_source_discordance` and
    the two lesion-group breakdowns all come from it rather than from any
    curated assertion -- the atlas mirrors gene-disease validity, it does not
    author it. A gene absent from `validity` gets `uncurated()`: no authority has
    assessed it, which is a fact worth publishing rather than an absent key.
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
    # `sorted` rather than insertion order, which would follow `corpus.assertions`
    # and so the filenames it was loaded from. That is stable, so this is not what
    # keeps the build byte-identical today; it is what makes the ordering rule at
    # the top of this module true of the returned dict as well, and what stops a
    # consumer that iterates `facts.items()` into a JSON array having its gene
    # index reordered by an unrelated file rename.
    for gene in sorted(by_gene):
        assertions = by_gene[gene]
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

        # `strongest` raises on an empty sequence, so it is only ever called
        # behind `if mirrored`. Empty `mirrored` publishes `None`, matching
        # `uncurated()` and a mirror that curated the gene under a term that
        # maps to no rung at all (GenCC's `Supportive`, mapped to `None`).
        headline = strongest(mirrored) if mirrored else None
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
            functional_count=functional_counts[gene],
            publications=tuple(sorted(publications)),
        )
    return facts
