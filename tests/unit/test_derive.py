# tests/unit/test_derive.py
from datetime import date
from pathlib import Path

from chd_atlas.build.derive import gene_facts
from chd_atlas.corpus import Corpus
from chd_atlas.models.assertion import Evidence, GeneDiseaseAssertion, SupplementaryLocator
from chd_atlas.models.functional import FunctionalEvidence
from chd_atlas.vocab import Classification, EvidenceClass, LesionGroup


def _evidence(**overrides: object) -> Evidence:
    payload: dict[str, object] = {
        "publication": "PMID:8988165",
        "evidence_class": "genetic_case",
        "locator": SupplementaryLocator(file="T.xlsx"),
        "summary": "s",
        "strength": "strong",
    }
    payload.update(overrides)
    return Evidence.model_validate(payload)


def _assertion(**overrides: object) -> GeneDiseaseAssertion:
    payload: dict[str, object] = {
        "id": "CHDA:AST:0000001",
        "gene": "HGNC:11604",
        "phenotypes": ["HP:0001631"],
        "lesion_groups": ["septal"],
        "classification": "definitive",
        "inheritance": ["AD"],
        "mechanism": "haploinsufficiency",
        "syndromic": "both",
        "evidence": [_evidence()],
        "source_tier": "own_curation",
        "curator": "c",
        "curated_on": date(2026, 7, 1),
        "last_reviewed": date(2026, 7, 1),
    }
    payload.update(overrides)
    return GeneDiseaseAssertion.model_validate(payload)


def _functional(**overrides: object) -> FunctionalEvidence:
    payload: dict[str, object] = {
        "id": "CHDA:FUN:0000001",
        "gene": "HGNC:11604",
        "organism": "NCBITaxon:10090",
        "perturbation": "knockout",
        "zygosity": "heterozygous",
        "cardiac_phenotype": ["MP:0000284"],
        "phenocopies_human": "partial",
        "rescue_outcome": "not_attempted",
        "publication": "PMID:8988165",
    }
    payload.update(overrides)
    return FunctionalEvidence.model_validate(payload)


def _corpus(**overrides: object) -> Corpus:
    payload: dict[str, object] = {"root": Path("."), "assertions": (_assertion(),)}
    payload.update(overrides)
    return Corpus(**payload)  # type: ignore[arg-type]


def test_headline_confidence_is_the_strongest_classification() -> None:
    corpus = _corpus(
        assertions=(
            _assertion(id="CHDA:AST:0000001", classification="moderate"),
            _assertion(id="CHDA:AST:0000002", classification="definitive"),
        )
    )

    facts = gene_facts(corpus)

    assert facts["HGNC:11604"].headline_confidence == Classification.DEFINITIVE


def test_the_conflict_flag_is_set_by_contradiction_and_by_nothing_else() -> None:
    """Spec 5.2: a contested gene must never be displayed as settled.

    `strongest` ranks refuted below definitive on one linear scale, so rank alone
    resolves the first gene here to `definitive` and buries the refutation. The
    flag is what stops the browse layer presenting it as settled.

    Both directions in one test, because a flag pinned to a constant is the
    obvious way this breaks and only one of the two assertions notices each
    constant: pinned true, the uncontested corpus is the sole witness; pinned
    false, the contested one is.
    """
    contested = _corpus(
        assertions=(
            _assertion(id="CHDA:AST:0000001", classification="definitive"),
            _assertion(id="CHDA:AST:0000002", classification="refuted"),
        )
    )

    facts = gene_facts(contested)

    assert facts["HGNC:11604"].headline_confidence == Classification.DEFINITIVE
    assert facts["HGNC:11604"].has_conflicting_evidence is True
    assert gene_facts(_corpus())["HGNC:11604"].has_conflicting_evidence is False


def test_each_gene_is_derived_only_from_its_own_assertions() -> None:
    """No accumulator may leak from one gene into the next.

    The governing requirement is that a gene is never displayed as better
    evidenced than its evidence, and the way this module would break it is not a
    ranking error but a scoping one: an accumulator hoisted out of the per-gene
    loop, or an inner loop reading `corpus.assertions` instead of this gene's
    slice. Either publishes `HGNC:11604` — curated `limited` — as `definitive`,
    because some *other* gene in the corpus is definitive.

    A single-gene fixture cannot see any of that, and every other test in this
    file uses one. The two genes here are given disjoint values for each derived
    field so a leak cannot coincidentally agree with the truth, and the
    assertions are on the gene derived *second*, which is the one a leak reaches.
    """
    corpus = _corpus(
        assertions=(
            _assertion(
                id="CHDA:AST:0000001",
                gene="HGNC:11599",
                classification="definitive",
                lesion_groups=["avsd"],
                evidence=[
                    _evidence(publication="PMID:1001", evidence_class="genetic_segregation")
                ],
            ),
            _assertion(
                id="CHDA:AST:0000002",
                gene="HGNC:11604",
                classification="limited",
                lesion_groups=["septal"],
                evidence=[_evidence(publication="PMID:2002", evidence_class="genetic_case")],
            ),
        )
    )

    facts = gene_facts(corpus)

    second = facts["HGNC:11604"]
    assert second.gene == "HGNC:11604"
    assert second.headline_confidence == Classification.LIMITED
    assert second.lesion_groups == (LesionGroup.SEPTAL,)
    assert second.publications == ("PMID:2002",)
    assert second.evidence_counts == {EvidenceClass.GENETIC_CASE: 1}

    # The gene derived first keeps its own values too, so nothing is fixed by
    # leaking in the other direction.
    assert facts["HGNC:11599"].headline_confidence == Classification.DEFINITIVE
    assert facts["HGNC:11599"].publications == ("PMID:1001",)


def test_confidence_is_broken_down_per_lesion_group() -> None:
    """The two halves of the fan-out, neither of which a disjoint fixture reaches.

    An assertion listing several lesion groups scores *every* one of them — a
    curator writing `[septal, conotruncal]` is claiming the gene for both, and
    dropping the tail would silently narrow the claim. And a group reaching more
    than one classification takes the strongest, not the first seen, which is why
    the weaker of septal's two is ordered first here: taking first-seen would
    publish septal as `limited` and this fixture is the only thing that notices.
    """
    corpus = _corpus(
        assertions=(
            _assertion(
                id="CHDA:AST:0000001",
                lesion_groups=["septal", "conotruncal"],
                classification="limited",
            ),
            _assertion(
                id="CHDA:AST:0000002", lesion_groups=["septal"], classification="definitive"
            ),
            _assertion(id="CHDA:AST:0000003", lesion_groups=["avsd"], classification="moderate"),
        )
    )

    facts = gene_facts(corpus)

    assert facts["HGNC:11604"].confidence_by_lesion_group == {
        LesionGroup.AVSD: Classification.MODERATE,
        LesionGroup.CONOTRUNCAL: Classification.LIMITED,
        LesionGroup.SEPTAL: Classification.DEFINITIVE,
    }
    # The array and the breakdown describe the same set of groups in the same
    # order. Pinned together so neither can start reporting a group the other
    # does not.
    assert facts["HGNC:11604"].lesion_groups == tuple(
        facts["HGNC:11604"].confidence_by_lesion_group
    )


def test_evidence_is_counted_per_class() -> None:
    corpus = _corpus(
        assertions=(
            _assertion(
                evidence=[
                    _evidence(evidence_class="genetic_case"),
                    _evidence(evidence_class="genetic_case"),
                    _evidence(evidence_class="genetic_segregation"),
                ]
            ),
        )
    )

    facts = gene_facts(corpus)

    assert facts["HGNC:11604"].evidence_counts == {
        EvidenceClass.GENETIC_CASE: 2,
        EvidenceClass.GENETIC_SEGREGATION: 1,
    }


def test_cited_publications_are_collected_and_deduplicated() -> None:
    corpus = _corpus(
        assertions=(
            _assertion(
                evidence=[
                    _evidence(publication="PMID:8988165"),
                    _evidence(publication="PMID:8988165"),
                ]
            ),
        )
    )

    assert gene_facts(corpus)["HGNC:11604"].publications == ("PMID:8988165",)


def test_a_corpus_with_no_assertions_derives_nothing() -> None:
    assert gene_facts(Corpus(root=Path("."))) == {}


def test_cited_publications_come_back_sorted() -> None:
    """They are collected through a `set`, and published as a JSON array.

    Array order survives `encode_json`, whose `sort_keys` orders dict keys and
    nothing else, so an unsorted tuple here reaches `genes/<HGNC>.json` verbatim.
    Set iteration order for strings depends on `PYTHONHASHSEED`, so dropping the
    sort yields a build that checksums differently on the CI runner than on the
    machine trying to reproduce it. No other assertion in this file notices.

    The fixture is six publications because this guard is probabilistic, not
    deterministic: *n* strings can land in an order that happens to equal sorted
    order, in which case the unsorted code passes. That chance falls sharply with
    *n*. Measured here, one process per seed, counting seeds on which removing
    the `sorted()` still passes: 48/200 at three publications, 1/200 at five,
    0/1000 at the six below. Nothing pins `PYTHONHASHSEED` in CI or in a
    `conftest.py`, so at three this was a guard that would have waved a
    regression through roughly one run in four.
    """
    corpus = _corpus(
        assertions=(
            _assertion(
                evidence=[
                    _evidence(publication="PMID:9"),
                    _evidence(publication="PMID:8988165"),
                    _evidence(publication="PMID:11729"),
                    _evidence(publication="PMID:3"),
                    _evidence(publication="PMID:40404"),
                    _evidence(publication="PMID:777"),
                ]
            ),
        )
    )

    # Lexicographic, not numeric: these are strings, and that is what the site
    # serves.
    assert gene_facts(corpus)["HGNC:11604"].publications == (
        "PMID:11729",
        "PMID:3",
        "PMID:40404",
        "PMID:777",
        "PMID:8988165",
        "PMID:9",
    )


def test_lesion_groups_are_sorted_rather_than_first_seen() -> None:
    """Also published as a JSON array, so its order is part of the output.

    Sorted by value rather than by declaration order in `LesionGroup`, so two
    genes whose assertions were curated in a different order still list the same
    groups the same way.
    """
    corpus = _corpus(
        assertions=(
            _assertion(id="CHDA:AST:0000001", lesion_groups=["septal", "avsd"]),
            _assertion(id="CHDA:AST:0000002", lesion_groups=["conotruncal"]),
        )
    )

    assert gene_facts(corpus)["HGNC:11604"].lesion_groups == (
        LesionGroup.AVSD,
        LesionGroup.CONOTRUNCAL,
        LesionGroup.SEPTAL,
    )


def test_assertions_are_counted_per_gene_and_the_genes_come_back_sorted() -> None:
    """The returned mapping is ordered by gene id, not by the order they arrived.

    Nothing else in this file looked at `facts`' own key order, so dropping the
    `sorted(by_gene)` in `gene_facts` passed the whole suite — which is what the
    second assertion is here for. The assertions below are deliberately curated
    HGNC:11604 first and HGNC:11599 last, so insertion order and sorted order
    disagree and the check cannot pass by coincidence. A consumer that iterates
    `facts.items()` into a JSON array would otherwise have its gene index
    reordered by an unrelated file rename.
    """
    corpus = _corpus(
        assertions=(
            _assertion(id="CHDA:AST:0000001"),
            _assertion(id="CHDA:AST:0000002"),
            _assertion(id="CHDA:AST:0000003", gene="HGNC:11599"),
        )
    )

    facts = gene_facts(corpus)

    assert facts["HGNC:11604"].assertion_count == 2
    assert facts["HGNC:11599"].assertion_count == 1
    assert list(facts) == ["HGNC:11599", "HGNC:11604"]


def test_every_functional_record_about_a_gene_is_counted() -> None:
    """Not only the ones an assertion cites.

    `referential.py` requires a cited record to be about the gene asserting it,
    but nothing requires a curated record to be cited at all. The gene bundle
    embeds every record about the gene, so counting only the cited ones would
    print a number in the browse index smaller than the list on the page it
    links to.
    """
    corpus = _corpus(
        functional=(
            _functional(id="CHDA:FUN:0000001"),
            _functional(id="CHDA:FUN:0000002"),
            _functional(id="CHDA:FUN:0000003", gene="HGNC:11599"),
        )
    )

    facts = gene_facts(corpus)

    assert facts["HGNC:11604"].functional_count == 2
    # HGNC:11599 has a functional record and no assertion, which is the state the
    # docstring calls deliberately absent: keying off the union of assertions and
    # functional records is a coherent alternative design, and this is what says
    # it is not the one chosen. A gene with nothing asserted has no confidence to
    # display, so there is no headline it could be given.
    assert "HGNC:11599" not in facts
    # And the other end of the count: a gene with no functional record at all
    # reports zero rather than raising, which is what the `Counter` buys over a
    # plain dict.
    assert gene_facts(_corpus())["HGNC:11604"].functional_count == 0
