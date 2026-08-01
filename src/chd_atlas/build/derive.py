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

from chd_atlas.corpus import Corpus
from chd_atlas.models.assertion import GeneDiseaseAssertion
from chd_atlas.vocab import (
    Classification,
    EvidenceClass,
    LesionGroup,
    has_conflicting_evidence,
    strongest,
)


@dataclass(frozen=True)
class GeneFacts:
    """What the browse layer needs about one gene without opening its bundle."""

    gene: str
    headline_confidence: Classification
    # Paired with `headline_confidence` at every display site. `strongest` ranks
    # on one linear scale, so a gene carrying both a definitive and a refuted
    # assertion resolves to definitive and buries the refutation; this is what
    # stops the browse layer showing a contested gene as settled.
    has_conflicting_evidence: bool
    lesion_groups: tuple[LesionGroup, ...]
    confidence_by_lesion_group: dict[LesionGroup, Classification]
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


def gene_facts(corpus: Corpus) -> dict[str, GeneFacts]:
    """Derive one `GeneFacts` per gene that carries at least one assertion.

    Keyed on HGNC id. Genes present in `mirrors/genes.tsv` but carrying no
    assertion are deliberately absent: the atlas browses curated claims, and a
    gene with nothing asserted has no confidence to display.
    """
    by_gene: dict[str, list[GeneDiseaseAssertion]] = {}
    for assertion in corpus.assertions:
        by_gene.setdefault(assertion.gene, []).append(assertion)

    # Annotated `str` rather than left to inference: `record.gene` is `HgncId`, a
    # NewType over `str`, so an inferred `Counter[HgncId]` cannot be indexed by
    # the plain `str` keys of `by_gene` under --strict. The two are the same
    # value at runtime; this keeps the key type consistent with the returned
    # `dict[str, GeneFacts]` rather than pushing the NewType through the API.
    functional_counts: Counter[str] = Counter(record.gene for record in corpus.functional)

    facts: dict[str, GeneFacts] = {}
    for gene, assertions in by_gene.items():
        classifications = [assertion.classification for assertion in assertions]

        per_group: dict[LesionGroup, list[Classification]] = {}
        counts: Counter[EvidenceClass] = Counter()
        publications: set[str] = set()
        for assertion in assertions:
            for group in assertion.lesion_groups:
                per_group.setdefault(group, []).append(assertion.classification)
            for evidence in assertion.evidence:
                counts[evidence.evidence_class] += 1
                publications.add(evidence.publication)

        # Sorted once and reused, so the array and the breakdown cannot come to
        # disagree about the order of the same set of groups.
        ordered_groups = sorted(per_group, key=lambda group: group.value)

        # `strongest` raises on an empty sequence. It cannot be reached with one
        # here: `by_gene` gains a key only when an assertion is appended, and
        # `per_group` only when a classification is, while the schema requires at
        # least one `lesion_groups` entry per assertion.
        facts[gene] = GeneFacts(
            gene=gene,
            headline_confidence=strongest(classifications),
            has_conflicting_evidence=has_conflicting_evidence(classifications),
            lesion_groups=tuple(ordered_groups),
            confidence_by_lesion_group={
                group: strongest(per_group[group]) for group in ordered_groups
            },
            evidence_counts={
                item: counts[item] for item in sorted(counts, key=lambda item: item.value)
            },
            assertion_count=len(assertions),
            functional_count=functional_counts[gene],
            publications=tuple(sorted(publications)),
        )
    return facts
