# tests/unit/test_derive.py
from datetime import date
from pathlib import Path

from chd_atlas.build.derive import gene_facts
from chd_atlas.build.validity import GeneValidity, ValidityRecord
from chd_atlas.corpus import Corpus
from chd_atlas.models.assertion import Evidence, LesionAssertion, SupplementaryLocator
from chd_atlas.models.functional import FunctionalEvidence
from chd_atlas.vocab import (
    AtlasCuration,
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


def _assertion(**overrides: object) -> LesionAssertion:
    payload: dict[str, object] = {
        "id": "CHDA:AST:0000001",
        "gene": "HGNC:11604",
        "phenotypes": ["HP:0001631"],
        "lesion_groups": ["septal"],
        "inheritance": ["AD"],
        "mechanism": "haploinsufficiency",
        "syndromic": "both",
        "evidence": [_evidence()],
        "curator": "c",
        "curated_on": date(2026, 7, 1),
        "last_reviewed": date(2026, 7, 1),
    }
    payload.update(overrides)
    return LesionAssertion.model_validate(payload)


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

    facts = gene_facts(_corpus(), validity, published={"HGNC:11604"})

    assert facts["HGNC:11604"].headline_confidence == Classification.DEFINITIVE


def test_an_uncurated_gene_publishes_no_headline_and_no_breakdown() -> None:
    """A gene absent from `validity` gets no headline, not a fabricated one.

    `None` here is not `Classification.NO_KNOWN_ASSOCIATION` — that
    classification is itself an assessed verdict ("a panel looked and found
    nothing"), and no authority has assessed this gene at all. The per-group
    breakdown is empty for the identical reason: there is no mirrored
    classification to collapse under any group.
    """
    fact = gene_facts(_corpus(), validity={}, published={"HGNC:11604"})["HGNC:11604"]

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

    facts = gene_facts(_corpus(), contested, published={"HGNC:11604"})

    assert facts["HGNC:11604"].headline_confidence == Classification.DEFINITIVE
    assert facts["HGNC:11604"].has_conflicting_evidence is True
    assert (
        gene_facts(_corpus(), {}, published={"HGNC:11604"})["HGNC:11604"].has_conflicting_evidence
        is False
    )


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

    facts = gene_facts(corpus, validity, published={"HGNC:11599", "HGNC:11604"})

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
    *same* headline. Nothing about any one assertion's own `classification` field
    enters this at all.

    **The fixture also pins which records the headline is taken over, and that
    is what it measures since 2026-08-06.** ClinGen grades this gene `moderate`
    and a GenCC submitter grades it `definitive`, so a `strongest()` over every
    mirrored record publishes `definitive` and the panel-only rule publishes
    `moderate`. The two answers differ here by construction — before the gate
    widened they could not, because every published gene was ClinGen
    `definitive` and a max over a superset returned the same value either way.
    That is the degeneracy CLAUDE.md §4.30 records: two figures equal in the
    corpus are one figure to every test.
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

    facts = gene_facts(corpus, validity, published={"HGNC:11604"})

    assert facts["HGNC:11604"].confidence_by_lesion_group == {
        LesionGroup.AVSD: Classification.MODERATE,
        LesionGroup.CONOTRUNCAL: Classification.MODERATE,
        LesionGroup.SEPTAL: Classification.MODERATE,
    }
    assert facts["HGNC:11604"].headline_confidence is Classification.MODERATE


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

    facts = gene_facts(corpus, {}, published={"HGNC:11604"})

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

    assert gene_facts(corpus, {}, published={"HGNC:11604"})["HGNC:11604"].publications == (
        "PMID:8988165",
    )


def test_an_empty_published_set_derives_nothing() -> None:
    """The degenerate input: nothing published, so nothing derived.

    This used to read "a corpus with no assertions derives nothing", which D21
    made false in the direction that matters --
    `test_a_published_gene_with_no_assertion_still_gets_facts` is the same
    empty corpus with one published gene, and it derives one fact. What is left
    to pin is that an empty population returns an empty mapping rather than
    raising, since `published` is now the only thing that decides.
    """
    assert gene_facts(Corpus(root=Path(".")), {}, published=set()) == {}


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
    assert gene_facts(corpus, {}, published={"HGNC:11604"})["HGNC:11604"].publications == (
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

    assert gene_facts(corpus, {}, published={"HGNC:11604"})["HGNC:11604"].lesion_groups == (
        LesionGroup.AVSD,
        LesionGroup.CONOTRUNCAL,
        LesionGroup.SEPTAL,
    )


def test_assertions_are_counted_per_gene_and_the_genes_come_back_sorted() -> None:
    """The returned mapping is ordered by gene id, not by the order they arrived.

    Nothing else in this file looks at `facts`' own key order, so dropping the
    `sorted(published)` in `gene_facts` passes the whole suite bar the last
    assertion below. What the build actually passes is a `set` from
    `published_genes`, whose iteration order for strings varies with
    `PYTHONHASHSEED`; a set here would therefore catch that mutation on some
    seeds and wave it through on others, which §4.12 of `CLAUDE.md` records as
    not being a guard at all. `published` is declared `Collection[str]`, so this
    hands it an ordered `list` in the reverse of sorted order instead: the
    mutant then fails on every run rather than on a fraction of them.

    A consumer that iterates `facts.items()` into a JSON array would otherwise
    have its gene index reordered between two builds of one commit.
    """
    corpus = _corpus(
        assertions=(
            _assertion(id="CHDA:AST:0000001"),
            _assertion(id="CHDA:AST:0000002"),
            _assertion(id="CHDA:AST:0000003", gene="HGNC:11599"),
        )
    )

    facts = gene_facts(corpus, {}, published=["HGNC:11604", "HGNC:11599"])

    assert facts["HGNC:11604"].assertion_count == 2
    assert facts["HGNC:11599"].assertion_count == 1
    assert list(facts) == ["HGNC:11599", "HGNC:11604"]


def test_a_published_gene_with_no_assertion_still_gets_facts() -> None:
    """The 22-of-23 case, and the whole point of D21.

    Before this, `gene_facts` keyed on the asserted genes, so a gene an expert
    panel calls definitive and the atlas has not curated produced no facts, no
    bundle and no page. Its mirrored confidence still comes from the mirror --
    the assertion was never where it came from.
    """
    corpus = _corpus(assertions=())
    validity = {"HGNC:1": _gene_validity()}

    facts = gene_facts(corpus, validity, published={"HGNC:1"})

    assert set(facts) == {"HGNC:1"}
    fact = facts["HGNC:1"]
    assert fact.headline_confidence is Classification.DEFINITIVE
    assert fact.assertion_count == 0
    assert fact.atlas_curation is AtlasCuration.NOT_YET_CURATED
    assert fact.lesion_groups == ()
    assert fact.confidence_by_lesion_group == {}


def test_atlas_curation_is_curated_exactly_when_an_assertion_exists() -> None:
    """The invariant a consumer filters on, asserted in both directions.

    One-directional would pass on a hardcoded "curated", which is the mutant
    that matters: it would advertise atlas evidence for 22 genes that have none.
    """
    corpus = _corpus(assertions=(_assertion(gene="HGNC:1"),))
    validity = {"HGNC:1": _gene_validity(), "HGNC:2": _gene_validity()}

    facts = gene_facts(corpus, validity, published={"HGNC:1", "HGNC:2"})

    assert facts["HGNC:1"].atlas_curation is AtlasCuration.CURATED
    assert facts["HGNC:1"].assertion_count == 1
    assert facts["HGNC:2"].atlas_curation is AtlasCuration.NOT_YET_CURATED
    assert facts["HGNC:2"].assertion_count == 0


def test_an_asserted_gene_outside_the_published_set_gets_no_facts() -> None:
    """`published` is the population, full stop.

    A curated assertion does not admit a gene: D21 is the gate and D37 keeps
    candidate genes out of the definitive set. Without this, a curator adding an
    assertion for a gene no panel has graded would silently publish it beside
    the panel-adjudicated ones.
    """
    corpus = _corpus(assertions=(_assertion(gene="HGNC:9"),))

    facts = gene_facts(corpus, {"HGNC:1": _gene_validity()}, published={"HGNC:1"})

    assert set(facts) == {"HGNC:1"}


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

    facts = gene_facts(corpus, {}, published={"HGNC:11604"})

    assert facts["HGNC:11604"].functional_count == 2
    # HGNC:11599 carries a functional record and is outside `published`, so it
    # gets no facts: keying off the union of `published` and the functionally
    # evidenced genes is a coherent alternative design, and this is what says it
    # is not the one chosen. Its record is counted for nobody and published
    # nowhere -- the cost `bundles.py` pins with
    # `test_a_gene_with_evidence_but_outside_the_published_set_is_not_published`.
    assert "HGNC:11599" not in facts
    # And the other end of the count: a gene with no functional record at all
    # reports zero rather than raising, which is what the `Counter` buys over a
    # plain dict.
    assert gene_facts(_corpus(), {}, published={"HGNC:11604"})["HGNC:11604"].functional_count == 0
