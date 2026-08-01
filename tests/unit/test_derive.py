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


def test_a_contested_gene_carries_the_conflict_flag() -> None:
    """Spec 5.2: a contested gene must never be displayed as settled.

    `strongest` ranks refuted below definitive on one linear scale, so rank alone
    resolves this gene to `definitive` and buries the refutation. The flag is what
    stops the browse layer presenting it as settled.
    """
    corpus = _corpus(
        assertions=(
            _assertion(id="CHDA:AST:0000001", classification="definitive"),
            _assertion(id="CHDA:AST:0000002", classification="refuted"),
        )
    )

    facts = gene_facts(corpus)

    assert facts["HGNC:11604"].headline_confidence == Classification.DEFINITIVE
    assert facts["HGNC:11604"].has_conflicting_evidence is True


def test_an_uncontested_gene_does_not_carry_the_flag() -> None:
    facts = gene_facts(_corpus())
    assert facts["HGNC:11604"].has_conflicting_evidence is False


def test_confidence_is_broken_down_per_lesion_group() -> None:
    corpus = _corpus(
        assertions=(
            _assertion(
                id="CHDA:AST:0000001", lesion_groups=["septal"], classification="definitive"
            ),
            _assertion(
                id="CHDA:AST:0000002", lesion_groups=["conotruncal"], classification="limited"
            ),
        )
    )

    facts = gene_facts(corpus)

    assert facts["HGNC:11604"].confidence_by_lesion_group == {
        LesionGroup.CONOTRUNCAL: Classification.LIMITED,
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
    Set iteration order for strings depends on `PYTHONHASHSEED`, which makes that
    a build that checksums differently on the CI runner than on the machine that
    reproduces it. Dropping the sort leaves every other assertion in this file
    green, so this is the one that pins it.
    """
    corpus = _corpus(
        assertions=(
            _assertion(
                evidence=[
                    _evidence(publication="PMID:9"),
                    _evidence(publication="PMID:8988165"),
                    _evidence(publication="PMID:11729"),
                ]
            ),
        )
    )

    assert gene_facts(corpus)["HGNC:11604"].publications == (
        "PMID:11729",
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


def test_assertions_are_counted_per_gene() -> None:
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

    assert gene_facts(corpus)["HGNC:11604"].functional_count == 2


def test_a_gene_with_no_functional_records_counts_zero() -> None:
    assert gene_facts(_corpus())["HGNC:11604"].functional_count == 0
