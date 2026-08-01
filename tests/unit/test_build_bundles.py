# tests/unit/test_build_bundles.py
import json
from datetime import date
from pathlib import Path
from typing import Any

from chd_atlas.build.bundles import build_genes
from chd_atlas.build.emit import Emitter
from chd_atlas.build.omics import ModalitySummary
from chd_atlas.corpus import Corpus
from chd_atlas.models.assertion import Evidence, GeneDiseaseAssertion, SupplementaryLocator
from chd_atlas.models.functional import FunctionalEvidence

# Two real HGNC ids whose lexical order is the reverse of their symbols':
# "HGNC:11604" sorts before "HGNC:4173", while "GATA4" sorts before "TBX5".
# Several tests below rely on that disagreement to tell one ordering rule from
# the other.
TBX5 = "HGNC:11604"
GATA4 = "HGNC:4173"
SYMBOLS = {TBX5: "TBX5", GATA4: "GATA4"}


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
        "gene": TBX5,
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
        "gene": TBX5,
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


def _summary(**overrides: object) -> ModalitySummary:
    payload: dict[str, Any] = {
        "count": 3,
        "shards": ["omics/expression/GSE1000.json"],
        "top": [{"gene": TBX5, "fdr": 0.01}],
    }
    payload.update(overrides)
    return ModalitySummary(**payload)


def _variant(**overrides: object) -> dict[str, Any]:
    payload: dict[str, Any] = {"vrs_id": "ga4gh:VA.x", "gene": TBX5, "chrom": "12", "pos": 100}
    payload.update(overrides)
    return payload


def _read(root: Path, relative: str) -> dict[str, Any]:
    """One emitted artifact, read back from the file that was written."""
    payload: dict[str, Any] = json.loads((root / relative).read_text())
    return payload


def _entries(root: Path) -> list[dict[str, Any]]:
    entries: list[dict[str, Any]] = _read(root, "genes/index.json")["genes"]
    return entries


def _by_gene(root: Path) -> dict[str, dict[str, Any]]:
    return {str(entry["gene"]): entry for entry in _entries(root)}


def test_the_index_carries_what_the_browser_filters_on(tmp_path: Path) -> None:
    emitter = Emitter(root=tmp_path)

    build_genes(_corpus(), emitter, symbols=SYMBOLS, omics={}, variants={})

    entry = _entries(tmp_path)[0]
    assert entry["gene"] == TBX5
    assert entry["symbol"] == "TBX5"
    assert entry["headline_confidence"] == "definitive"
    assert entry["has_conflicting_evidence"] is False
    assert entry["lesion_groups"] == ["septal"]
    assert entry["confidence_by_lesion_group"] == {"septal": "definitive"}
    assert entry["evidence_counts"] == {"genetic_case": 1}
    assert entry["assertion_count"] == 1
    assert entry["functional_count"] == 0
    assert entry["variant_count"] == 0
    assert entry["bundle"] == "genes/HGNC_11604.json"


def test_the_index_carries_nothing_the_browser_does_not_filter_on(tmp_path: Path) -> None:
    """The index is downloaded by every visitor before they pick a gene.

    Its size is the cost of the browse page, so what belongs here is what ranks
    or filters a row and the path to fetch the rest. An assertion, its evidence,
    a variant list or a publication list added here would be downloaded by
    everyone to be read by nobody — the bundle already carries all four.
    """
    emitter = Emitter(root=tmp_path)

    build_genes(_corpus(), emitter, symbols=SYMBOLS, omics={}, variants={})

    assert set(_entries(tmp_path)[0]) == {
        "gene",
        "symbol",
        "headline_confidence",
        "has_conflicting_evidence",
        "lesion_groups",
        "confidence_by_lesion_group",
        "evidence_counts",
        "assertion_count",
        "functional_count",
        "variant_count",
        "bundle",
    }


def test_a_bundle_carries_the_whole_gene_page_and_nothing_more(tmp_path: Path) -> None:
    """A gene detail page is one fetch, so every section it renders is here."""
    emitter = Emitter(root=tmp_path)

    build_genes(_corpus(), emitter, symbols=SYMBOLS, omics={}, variants={})

    assert set(_read(tmp_path, "genes/HGNC_11604.json")) == {
        "gene",
        "symbol",
        "headline_confidence",
        "has_conflicting_evidence",
        "lesion_groups",
        "publications",
        "assertions",
        "functional",
        "variants",
        "omics",
    }


def test_the_index_links_each_gene_to_its_bundle(tmp_path: Path) -> None:
    """The client must never build a bundle path itself."""
    emitter = Emitter(root=tmp_path)

    build_genes(_corpus(), emitter, symbols={}, omics={}, variants={})

    assert _entries(tmp_path)[0]["bundle"] == "genes/HGNC_11604.json"
    assert (tmp_path / "genes" / "HGNC_11604.json").is_file()


def test_every_bundle_path_the_index_advertises_was_written(tmp_path: Path) -> None:
    """Checked over the whole index rather than on one hardcoded name.

    `bundle` is the URL a detail page is fetched from, so an index entry whose
    path was never written is a dead link that nothing else in the build
    notices: the manifest lists what was written, not what was promised. The
    two are produced by one `gene_bundle_path` call per gene, and this is what
    holds them together if that stops being true.
    """
    corpus = _corpus(
        assertions=(
            _assertion(id="CHDA:AST:0000001", gene=TBX5),
            _assertion(id="CHDA:AST:0000002", gene=GATA4),
        )
    )
    emitter = Emitter(root=tmp_path)

    build_genes(corpus, emitter, symbols=SYMBOLS, omics={}, variants={})

    advertised = [str(entry["bundle"]) for entry in _entries(tmp_path)]
    assert len(advertised) == 2
    for relative in advertised:
        assert relative in emitter.checksums, f"{relative} is in the index but was never written"
        assert (tmp_path / relative).is_file(), f"{relative} is in the manifest but not on disk"


def test_a_contested_gene_is_flagged_in_both_the_index_and_the_bundle(
    tmp_path: Path,
) -> None:
    """Spec 5.2: a contested gene must never be displayed as settled.

    `strongest()` ranks on one linear scale where definitive outranks refuted,
    so a gene carrying both resolves to definitive and the refutation is
    invisible in `headline_confidence` alone. The flag is the other half of that
    pair, and this module is its first consumer.

    Both genes are asserted here, contested and uncontested, and both the index
    entry and the bundle are read. A writer that hardcoded either value, or that
    dropped the key and let a reader's `.get(..., False)` fill it in, fails on
    one of the assertions below rather than passing on a fixture that only ever
    shows one answer. The `headline_confidence` assertions are what make the
    point: TBX5 still reads `definitive`, so the flag is carrying information
    that the confidence cannot — and GATA4 reads `moderate`, so a headline
    hardcoded to the value a contested gene happens to resolve to fails here
    too.
    """
    corpus = _corpus(
        assertions=(
            _assertion(id="CHDA:AST:0000001", gene=TBX5, classification="definitive"),
            _assertion(id="CHDA:AST:0000002", gene=TBX5, classification="refuted"),
            _assertion(id="CHDA:AST:0000003", gene=GATA4, classification="moderate"),
        )
    )
    emitter = Emitter(root=tmp_path)

    build_genes(corpus, emitter, symbols=SYMBOLS, omics={}, variants={})

    entries = _by_gene(tmp_path)
    assert entries[TBX5]["headline_confidence"] == "definitive"
    assert entries[TBX5]["has_conflicting_evidence"] is True
    assert entries[GATA4]["headline_confidence"] == "moderate"
    assert entries[GATA4]["has_conflicting_evidence"] is False

    contested = _read(tmp_path, "genes/HGNC_11604.json")
    settled = _read(tmp_path, "genes/HGNC_4173.json")
    assert contested["headline_confidence"] == "definitive"
    assert contested["has_conflicting_evidence"] is True
    assert settled["headline_confidence"] == "moderate"
    assert settled["has_conflicting_evidence"] is False


def test_confidence_is_broken_down_by_lesion_group(tmp_path: Path) -> None:
    """The breakdown is per group, so it disagrees with the headline by design.

    TBX5 here is definitive for septal disease and refuted for conotruncal. A
    browse row that showed the headline against every group would claim the
    refuted association is definitive.
    """
    corpus = _corpus(
        assertions=(
            _assertion(
                id="CHDA:AST:0000001", lesion_groups=["septal"], classification="definitive"
            ),
            _assertion(
                id="CHDA:AST:0000002", lesion_groups=["conotruncal"], classification="refuted"
            ),
        )
    )
    emitter = Emitter(root=tmp_path)

    build_genes(corpus, emitter, symbols=SYMBOLS, omics={}, variants={})

    entry = _entries(tmp_path)[0]
    assert entry["confidence_by_lesion_group"] == {
        "conotruncal": "refuted",
        "septal": "definitive",
    }
    # Sorted by value, and the same set of groups the breakdown keys on.
    assert entry["lesion_groups"] == ["conotruncal", "septal"]
    assert _read(tmp_path, "genes/HGNC_11604.json")["lesion_groups"] == ["conotruncal", "septal"]


def test_evidence_counts_are_carried_per_evidence_class(tmp_path: Path) -> None:
    """Two classes with different counts, so a total or a length cannot pass."""
    corpus = _corpus(
        assertions=(
            _assertion(
                evidence=[
                    _evidence(evidence_class="genetic_case"),
                    _evidence(evidence_class="genetic_case", publication="PMID:9"),
                    _evidence(
                        evidence_class="functional_model",
                        functional_evidence="CHDA:FUN:0000001",
                    ),
                ]
            ),
        ),
        functional=(_functional(),),
    )
    emitter = Emitter(root=tmp_path)

    build_genes(corpus, emitter, symbols=SYMBOLS, omics={}, variants={})

    assert _entries(tmp_path)[0]["evidence_counts"] == {"functional_model": 1, "genetic_case": 2}


def test_the_browse_counts_match_the_bundle_they_link_to(tmp_path: Path) -> None:
    """A count smaller than the page it links to reads as missing data.

    The three counts are the browse layer's promise about what opening a gene
    will show, and each is derived separately from the list the bundle carries —
    `assertion_count` and `functional_count` in `derive.py`, `variant_count`
    here. The two genes differ on all three, so a count wired to the wrong
    collection, or to a constant, fails rather than coinciding.
    """
    corpus = _corpus(
        assertions=(
            _assertion(id="CHDA:AST:0000001", gene=TBX5),
            _assertion(id="CHDA:AST:0000002", gene=TBX5),
            _assertion(id="CHDA:AST:0000003", gene=GATA4),
        ),
        functional=(
            _functional(id="CHDA:FUN:0000001", gene=TBX5),
            _functional(id="CHDA:FUN:0000002", gene=TBX5),
            _functional(id="CHDA:FUN:0000003", gene=TBX5),
        ),
    )
    emitter = Emitter(root=tmp_path)

    build_genes(
        corpus,
        emitter,
        symbols=SYMBOLS,
        omics={},
        variants={
            TBX5: [_variant(vrs_id="ga4gh:VA.1"), _variant(vrs_id="ga4gh:VA.2")],
            GATA4: [_variant(vrs_id="ga4gh:VA.3", gene=GATA4)],
        },
    )

    counted = {
        str(entry["gene"]): (
            entry["assertion_count"],
            entry["functional_count"],
            entry["variant_count"],
        )
        for entry in _entries(tmp_path)
    }
    assert counted == {TBX5: (2, 3, 2), GATA4: (1, 0, 1)}
    for entry in _entries(tmp_path):
        bundle = _read(tmp_path, str(entry["bundle"]))
        assert entry["assertion_count"] == len(bundle["assertions"])
        assert entry["functional_count"] == len(bundle["functional"])
        assert entry["variant_count"] == len(bundle["variants"])


def test_a_gene_absent_from_the_registry_keeps_its_hgnc_id_as_its_label(
    tmp_path: Path,
) -> None:
    """`mirrors/genes.tsv` is a mirror, so a gene can be asserted before it lands.

    The fallback is the id rather than the empty string because the browse row
    and the page heading are rendered from this value: blank reads as a bug in
    the site, while `HGNC:4173` reads as a gene whose symbol is not in yet and
    is still enough to search for. TBX5 is in the registry in the same call, so
    a writer that always falls back fails here too.
    """
    corpus = _corpus(
        assertions=(
            _assertion(id="CHDA:AST:0000001", gene=TBX5),
            _assertion(id="CHDA:AST:0000002", gene=GATA4),
        )
    )
    emitter = Emitter(root=tmp_path)

    build_genes(corpus, emitter, symbols={TBX5: "TBX5"}, omics={}, variants={})

    entries = _by_gene(tmp_path)
    assert entries[TBX5]["symbol"] == "TBX5"
    assert entries[GATA4]["symbol"] == GATA4
    assert _read(tmp_path, "genes/HGNC_4173.json")["symbol"] == GATA4


def test_a_bundle_carries_its_assertions_in_full(tmp_path: Path) -> None:
    """Including the evidence, which is the record a curator is judged on."""
    emitter = Emitter(root=tmp_path)

    build_genes(_corpus(), emitter, symbols=SYMBOLS, omics={}, variants={})

    bundle = _read(tmp_path, "genes/HGNC_11604.json")
    assertion = bundle["assertions"][0]
    assert assertion["id"] == "CHDA:AST:0000001"
    assert assertion["evidence"][0]["summary"] == "s"
    assert assertion["evidence"][0]["locator"] == {
        "kind": "supplementary",
        "file": "T.xlsx",
        "sheet": None,
        "row": None,
    }
    # Dumped in JSON mode, so a date is a string rather than a `date` that
    # `json.dumps` refuses part-way through a build.
    assert assertion["curated_on"] == "2026-07-01"
    assert assertion["classification"] == "definitive"


def test_a_bundle_carries_every_functional_record_about_the_gene(tmp_path: Path) -> None:
    """Every record about the gene, not only those an assertion cites.

    `referential.py` requires a *cited* functional record to be about the
    asserting gene, but never requires a curated record to be cited at all, so
    the two populations differ. `derive.functional_count` counts this same
    population, and the browse count would otherwise exceed the page.
    """
    corpus = _corpus(
        functional=(
            _functional(id="CHDA:FUN:0000001", gene=TBX5),
            _functional(id="CHDA:FUN:0000002", gene=TBX5),
        )
    )
    emitter = Emitter(root=tmp_path)

    build_genes(corpus, emitter, symbols=SYMBOLS, omics={}, variants={})

    bundle = _read(tmp_path, "genes/HGNC_11604.json")
    assert [record["id"] for record in bundle["functional"]] == [
        "CHDA:FUN:0000001",
        "CHDA:FUN:0000002",
    ]
    assert bundle["functional"][0]["organism"] == "NCBITaxon:10090"


def test_a_bundle_carries_the_omics_summaries_verbatim(tmp_path: Path) -> None:
    """The whole `ModalitySummary`, not a re-listing of the fields known today.

    A field added to `ModalitySummary` should reach the bundle without this
    module being edited; whole-value equality is what fails if the summary is
    rebuilt key by key instead of passed through.
    """
    emitter = Emitter(root=tmp_path)

    build_genes(
        _corpus(),
        emitter,
        symbols=SYMBOLS,
        omics={TBX5: {"expression": _summary()}},
        variants={},
    )

    assert _read(tmp_path, "genes/HGNC_11604.json")["omics"] == {
        "expression": {
            "count": 3,
            "shards": ["omics/expression/GSE1000.json"],
            "top": [{"gene": TBX5, "fdr": 0.01}],
        }
    }


def test_a_bundle_embeds_its_variants_rather_than_linking_them(tmp_path: Path) -> None:
    """Affordable because curation bounds the count, not because the data does.

    The order is the one `build_variants` produced — karyotype order across
    shards, then each shard's own validated sort order — and nothing here
    reorders it.
    """
    emitter = Emitter(root=tmp_path)

    build_genes(
        _corpus(),
        emitter,
        symbols=SYMBOLS,
        omics={},
        variants={TBX5: [_variant(vrs_id="ga4gh:VA.2"), _variant(vrs_id="ga4gh:VA.1")]},
    )

    bundle = _read(tmp_path, "genes/HGNC_11604.json")
    assert [row["vrs_id"] for row in bundle["variants"]] == ["ga4gh:VA.2", "ga4gh:VA.1"]
    assert bundle["variants"][0]["pos"] == 100


def test_a_gene_with_no_omics_or_variants_gets_empty_containers(tmp_path: Path) -> None:
    """A consumer should read `bundle["omics"]`, never guard for its absence."""
    emitter = Emitter(root=tmp_path)

    build_genes(_corpus(), emitter, symbols=SYMBOLS, omics={}, variants={})

    bundle = _read(tmp_path, "genes/HGNC_11604.json")
    assert bundle["omics"] == {}
    assert bundle["variants"] == []


def test_a_bundle_holds_only_its_own_genes_evidence(tmp_path: Path) -> None:
    """Every collection in the bundle is per-gene, and each is filtered separately.

    Both genes carry all four kinds of evidence here, so a collection joined on
    the wrong key — or not joined at all — shows up as another gene's records in
    this bundle rather than as an empty one.
    """
    corpus = _corpus(
        assertions=(
            _assertion(id="CHDA:AST:0000001", gene=TBX5),
            _assertion(
                id="CHDA:AST:0000002", gene=GATA4, evidence=[_evidence(publication="PMID:9")]
            ),
        ),
        functional=(
            _functional(id="CHDA:FUN:0000001", gene=TBX5),
            _functional(id="CHDA:FUN:0000002", gene=GATA4),
        ),
    )
    emitter = Emitter(root=tmp_path)

    build_genes(
        corpus,
        emitter,
        symbols=SYMBOLS,
        omics={TBX5: {"expression": _summary()}, GATA4: {"profiles": _summary(count=9)}},
        variants={TBX5: [_variant()], GATA4: [_variant(gene=GATA4, vrs_id="ga4gh:VA.y")]},
    )

    bundle = _read(tmp_path, "genes/HGNC_4173.json")
    assert bundle["gene"] == GATA4
    assert [item["id"] for item in bundle["assertions"]] == ["CHDA:AST:0000002"]
    assert [item["id"] for item in bundle["functional"]] == ["CHDA:FUN:0000002"]
    assert [row["vrs_id"] for row in bundle["variants"]] == ["ga4gh:VA.y"]
    assert list(bundle["omics"]) == ["profiles"]
    assert bundle["publications"] == ["PMID:9"]


def test_the_index_is_ordered_by_hgnc_id_rather_than_by_symbol(tmp_path: Path) -> None:
    """The index is a JSON array, so its order is published and must be fixed.

    `encode_json` sorts dict keys and never array elements, so whatever order
    this list reaches the emitter in is the order every consumer downloads. The
    two genes are given in the reverse of both candidate orders, and their
    symbols rank them the other way round from their ids, so a build that
    emitted encounter order or sorted on the symbol fails here.

    This pins the published contract end to end. It cannot attribute the sort:
    `gene_facts` already returns its genes sorted by id, so `bundles.py`
    sorting them again is defence against that changing rather than what
    produces the order today.
    """
    corpus = _corpus(
        assertions=(
            _assertion(id="CHDA:AST:0000002", gene=GATA4),
            _assertion(id="CHDA:AST:0000001", gene=TBX5),
        )
    )
    emitter = Emitter(root=tmp_path)

    build_genes(corpus, emitter, symbols=SYMBOLS, omics={}, variants={})

    assert [entry["gene"] for entry in _entries(tmp_path)] == [TBX5, GATA4]
    assert [entry["symbol"] for entry in _entries(tmp_path)] == ["TBX5", "GATA4"]


def test_bundle_assertions_and_functional_records_are_ordered_by_id(tmp_path: Path) -> None:
    """Both are JSON arrays whose order survives `encode_json` untouched.

    Given in reverse id order and grouped by encounter, so a build that kept
    corpus order fails on every run rather than on the fraction of runs a
    hash-ordered container would.
    """
    corpus = _corpus(
        assertions=(
            _assertion(id="CHDA:AST:0000003"),
            _assertion(id="CHDA:AST:0000001"),
            _assertion(id="CHDA:AST:0000002"),
        ),
        functional=(
            _functional(id="CHDA:FUN:0000003"),
            _functional(id="CHDA:FUN:0000001"),
            _functional(id="CHDA:FUN:0000002"),
        ),
    )
    emitter = Emitter(root=tmp_path)

    build_genes(corpus, emitter, symbols=SYMBOLS, omics={}, variants={})

    bundle = _read(tmp_path, "genes/HGNC_11604.json")
    assert [item["id"] for item in bundle["assertions"]] == [
        "CHDA:AST:0000001",
        "CHDA:AST:0000002",
        "CHDA:AST:0000003",
    ]
    assert [item["id"] for item in bundle["functional"]] == [
        "CHDA:FUN:0000001",
        "CHDA:FUN:0000002",
        "CHDA:FUN:0000003",
    ]


def test_a_bundle_lists_the_publications_its_evidence_cites(tmp_path: Path) -> None:
    """Sorted lexically, matching `publications.json` rather than contradicting it.

    Six of them because `derive.py` accumulates them through a `set`, whose
    iteration order for strings varies with `PYTHONHASHSEED`: at three PMIDs an
    unsorted tuple still comes out sorted on roughly one seed in four. The
    ordering rule is `derive.py`'s; what this pins is that the bundle publishes
    it rather than re-deriving one of its own.
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
    emitter = Emitter(root=tmp_path)

    build_genes(corpus, emitter, symbols=SYMBOLS, omics={}, variants={})

    assert _read(tmp_path, "genes/HGNC_11604.json")["publications"] == [
        "PMID:11729",
        "PMID:3",
        "PMID:40404",
        "PMID:777",
        "PMID:8988165",
        "PMID:9",
    ]


def test_a_gene_with_omics_and_variants_but_no_assertion_is_not_published(
    tmp_path: Path,
) -> None:
    """It has no confidence to display, so it gets no bundle and no index row.

    `gene_facts` omits such a gene deliberately — the atlas browses curated
    claims — and this module publishes what that returns. The evidence itself is
    still served, in `variants/<chrom>.json.gz` and in the omics shards, but
    nothing in the gene index links to it: reaching it means already knowing the
    HGNC id. That is the accepted cost of the rule, not an oversight, and this
    test is where it is pinned so that changing it is a decision rather than a
    side effect.
    """
    emitter = Emitter(root=tmp_path)

    build_genes(
        _corpus(),
        emitter,
        symbols=SYMBOLS,
        omics={GATA4: {"expression": _summary()}},
        variants={GATA4: [_variant(gene=GATA4)]},
    )

    assert [entry["gene"] for entry in _entries(tmp_path)] == [TBX5]
    assert set(emitter.checksums) == {"genes/index.json", "genes/HGNC_11604.json"}


def test_the_index_is_emitted_for_an_empty_corpus(tmp_path: Path) -> None:
    """A consumer telling "no genes yet" from "wrong URL" by reading a 404 gets
    it wrong, and the browse page fetches this before it knows either."""
    emitter = Emitter(root=tmp_path)

    build_genes(Corpus(root=Path(".")), emitter, symbols={}, omics={}, variants={})

    assert _read(tmp_path, "genes/index.json") == {"genes": []}
    assert set(emitter.checksums) == {"genes/index.json"}
