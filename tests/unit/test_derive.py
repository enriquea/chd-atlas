# tests/unit/test_derive.py
from datetime import date
from pathlib import Path

from chd_atlas.build.derive import gene_facts
from chd_atlas.build.validity import GeneValidity, ValidityRecord
from chd_atlas.corpus import Corpus
from chd_atlas.models.assertion import Evidence, GeneDiseaseAssertion, SupplementaryLocator
from chd_atlas.models.functional import FunctionalEvidence
from chd_atlas.vocab import (
    Classification,
    EvidenceClass,
    LesionGroup,
    ValiditySource,
    ValidityState,
)


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


def _validity_record(**overrides: object) -> ValidityRecord:
    payload: dict[str, object] = {
        "source": ValiditySource.CLINGEN,
        "classification_term": "Definitive",
        "classification": Classification.DEFINITIVE,
        "disease": "MONDO:0007732",
        "disease_label": "Holt-Oram syndrome",
        "moi": "AD",
        "report_url": None,
    }
    payload.update(overrides)
    return ValidityRecord(**payload)  # type: ignore[arg-type]


def _gene_validity(**overrides: object) -> GeneValidity:
    payload: dict[str, object] = {
        "state": ValidityState.EXPERT_CURATED,
        "records": (_validity_record(),),
        "has_source_discordance": False,
    }
    payload.update(overrides)
    return GeneValidity(**payload)  # type: ignore[arg-type]


def test_headline_confidence_is_the_strongest_mirrored_classification() -> None:
    """`strongest()` over the gene's *mirrored* records, never its assertions.

    Two `ValidityRecord`s standing in for two ClinGen curations of the same
    gene-disease pair over time, moderate then definitive — the shape
    `gene_validity()` actually produces when one mirror row supersedes another.
    """
    validity = {
        "HGNC:11604": _gene_validity(
            records=(
                _validity_record(classification=Classification.MODERATE),
                _validity_record(classification=Classification.DEFINITIVE),
            )
        )
    }

    facts = gene_facts(_corpus(), validity)

    assert facts["HGNC:11604"].headline_confidence == Classification.DEFINITIVE


def test_headline_confidence_comes_from_the_mirror_not_the_curated_record() -> None:
    """The curated assertion's own `classification` must never leak into it.

    The assertion below is curated `refuted`; the mirror says `definitive`. A
    derivation that still reads `assertion.classification` anywhere publishes
    `refuted` here — this only passes if the mirror is the sole source.
    """
    corpus = _corpus(assertions=(_assertion(classification="refuted"),))
    validity = {
        "HGNC:11604": _gene_validity(
            records=(_validity_record(classification=Classification.DEFINITIVE),)
        )
    }

    facts = gene_facts(corpus, validity)

    assert facts["HGNC:11604"].headline_confidence == Classification.DEFINITIVE


def test_an_uncurated_gene_publishes_no_headline_and_no_breakdown() -> None:
    """A gene absent from `validity` gets no headline, not a fabricated one.

    `None` here is not `Classification.NO_KNOWN_ASSOCIATION` — that
    classification is itself an assessed verdict ("a panel looked and found
    nothing"), and no authority has assessed this gene at all. The per-group
    breakdown is empty for the identical reason: there is no mirrored
    classification to collapse under any group.
    """
    fact = gene_facts(_corpus(), validity={})["HGNC:11604"]

    assert fact.headline_confidence is None
    assert fact.validity_state is ValidityState.UNCURATED
    assert fact.has_conflicting_evidence is False
    assert fact.has_source_discordance is False
    assert fact.confidence_by_lesion_group == {}
    assert fact.conflicting_lesion_groups == ()


def test_the_conflict_flag_is_set_by_contradiction_and_by_nothing_else() -> None:
    """Spec 5.2: a contested gene must never be displayed as settled.

    `strongest` ranks refuted below definitive on one linear scale, so rank alone
    resolves the mirrored records here to `definitive` and buries the
    refutation. The flag is what stops the browse layer presenting it as
    settled.

    Both directions in one test, because a flag pinned to a constant is the
    obvious way this breaks and only one of the two cases notices each
    constant: pinned true, the contested gene is the sole witness; pinned
    false, the uncontested one is.
    """
    contested = {
        "HGNC:11604": _gene_validity(
            records=(
                _validity_record(
                    source=ValiditySource.CLINGEN, classification=Classification.DEFINITIVE
                ),
                _validity_record(
                    source=ValiditySource.GENCC,
                    classification=Classification.REFUTED,
                    submitter="s",
                ),
            )
        )
    }

    facts = gene_facts(_corpus(), contested)

    assert facts["HGNC:11604"].headline_confidence == Classification.DEFINITIVE
    assert facts["HGNC:11604"].has_conflicting_evidence is True
    assert gene_facts(_corpus(), {})["HGNC:11604"].has_conflicting_evidence is False


def test_each_gene_is_derived_only_from_its_own_assertions_and_validity() -> None:
    """No accumulator may leak from one gene into the next.

    The governing requirement is that a gene is never displayed as better
    evidenced than its evidence, and the way this module would break it is not a
    ranking error but a scoping one: an accumulator hoisted out of the per-gene
    loop, or a `validity` lookup that used the wrong key, or an inner loop
    reading `corpus.assertions` instead of this gene's slice. Any of these
    publishes `HGNC:11604` — mirrored `limited` — as `definitive`, because some
    *other* gene in the corpus or in `validity` is definitive.

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
                lesion_groups=["avsd"],
                evidence=[_evidence(publication="PMID:1001", evidence_class="genetic_segregation")],
            ),
            _assertion(
                id="CHDA:AST:0000002",
                gene="HGNC:11604",
                lesion_groups=["septal"],
                evidence=[_evidence(publication="PMID:2002", evidence_class="genetic_case")],
            ),
        )
    )
    validity = {
        "HGNC:11599": _gene_validity(
            records=(_validity_record(classification=Classification.DEFINITIVE),)
        ),
        "HGNC:11604": _gene_validity(
            records=(_validity_record(classification=Classification.LIMITED),)
        ),
    }

    facts = gene_facts(corpus, validity)

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


def test_confidence_by_lesion_group_applies_the_mirrored_headline_to_every_declared_group() -> None:
    """The mirrors classify a gene against a disease, never against a lesion.

    So the per-group breakdown cannot show finer-grained confidence than the
    gene-wide mirrored answer: every lesion group a curated assertion names for
    the gene — three of them here, spread across two assertions — publishes the
    *same* `strongest()` of the gene's mirrored classifications. Nothing about
    any one assertion's own `classification` field enters this at all.
    """
    corpus = _corpus(
        assertions=(
            _assertion(id="CHDA:AST:0000001", lesion_groups=["septal", "conotruncal"]),
            _assertion(id="CHDA:AST:0000002", lesion_groups=["avsd"]),
        )
    )
    validity = {
        "HGNC:11604": _gene_validity(
            records=(
                _validity_record(classification=Classification.MODERATE),
                _validity_record(
                    source=ValiditySource.GENCC,
                    classification=Classification.DEFINITIVE,
                    submitter="s",
                ),
            )
        )
    }

    facts = gene_facts(corpus, validity)

    assert facts["HGNC:11604"].confidence_by_lesion_group == {
        LesionGroup.AVSD: Classification.DEFINITIVE,
        LesionGroup.CONOTRUNCAL: Classification.DEFINITIVE,
        LesionGroup.SEPTAL: Classification.DEFINITIVE,
    }


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

    facts = gene_facts(corpus, {})

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

    assert gene_facts(corpus, {})["HGNC:11604"].publications == ("PMID:8988165",)


def test_a_corpus_with_no_assertions_derives_nothing() -> None:
    assert gene_facts(Corpus(root=Path(".")), {}) == {}


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
    assert gene_facts(corpus, {})["HGNC:11604"].publications == (
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

    assert gene_facts(corpus, {})["HGNC:11604"].lesion_groups == (
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

    facts = gene_facts(corpus, {})

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

    facts = gene_facts(corpus, {})

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
    assert gene_facts(_corpus(), {})["HGNC:11604"].functional_count == 0
