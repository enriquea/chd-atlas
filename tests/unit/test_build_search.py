# tests/unit/test_build_search.py
import gzip
import json
from datetime import date
from pathlib import Path
from typing import Any

from chd_atlas.build.emit import Emitter
from chd_atlas.build.search import build_search
from chd_atlas.corpus import Corpus
from chd_atlas.models.assertion import GeneDiseaseAssertion
from chd_atlas.models.functional import FunctionalEvidence
from chd_atlas.models.literature import PhenotypeTerm, Publication

TBX5 = "HGNC:11604"
GATA4 = "HGNC:4173"
NKX2_5 = "HGNC:2488"

# Six genes, listed here in the reverse of their lexical order. Six is a measured
# choice, not a round number — see `test_each_kind_is_ordered_by_identifier`.
SIX_GENES = ("HGNC:7881", "HGNC:6188", GATA4, NKX2_5, "HGNC:12873", TBX5)


def _assertion(gene: str = TBX5, identifier: str = "CHDA:AST:0000001") -> GeneDiseaseAssertion:
    return GeneDiseaseAssertion.model_validate(
        {
            "id": identifier,
            "gene": gene,
            "phenotypes": ["HP:0001631"],
            "lesion_groups": ["septal"],
            "classification": "definitive",
            "inheritance": ["AD"],
            "mechanism": "haploinsufficiency",
            "syndromic": "isolated",
            "evidence": [
                {
                    "publication": "PMID:8988165",
                    "evidence_class": "genetic_case",
                    "locator": {"kind": "in_text", "element": "Figure 1"},
                    "summary": "s",
                    "strength": "strong",
                }
            ],
            "source_tier": "own_curation",
            "curator": "c",
            "curated_on": date(2026, 7, 1),
            "last_reviewed": date(2026, 7, 1),
        }
    )


def _publication(pmid: str = "PMID:8988165", **overrides: object) -> Publication:
    payload: dict[str, object] = {
        "id": pmid,
        "title": "Mutations in human TBX5 cause Holt-Oram syndrome",
        "journal": "Nature Genetics",
        # Not in sorted order, and not sortable into one: author order is data,
        # so a build that sorted this list would publish an authorship the
        # record does not claim.
        "authors": ["Li QY", "Basson CT"],
        "year": 1997,
        "study_type": "family_linkage",
    }
    payload.update(overrides)
    return Publication.model_validate(payload)


def _phenotype(
    identifier: str = "HP:0001629",
    label: str = "Ventricular septal defect",
    **overrides: object,
) -> PhenotypeTerm:
    payload: dict[str, object] = {
        "id": identifier,
        "label": label,
        "lesion_group": "septal",
        # Reverse of sorted order, for the same reason as the authors above.
        "synonyms": ["VSD", "Interventricular septal defect"],
    }
    payload.update(overrides)
    return PhenotypeTerm.model_validate(payload)


def _functional(gene: str) -> FunctionalEvidence:
    return FunctionalEvidence.model_validate(
        {
            "id": "CHDA:FUN:0000001",
            "gene": gene,
            "organism": "NCBITaxon:10090",
            "perturbation": "knockout",
            "zygosity": "heterozygous",
            "cardiac_phenotype": ["MP:0000284"],
            "phenocopies_human": "partial",
            "rescue_outcome": "not_attempted",
            "publication": "PMID:8988165",
        }
    )


def _corpus(**overrides: object) -> Corpus:
    payload: dict[str, object] = {
        "root": Path("."),
        # One of each kind. A32: the assertion is what makes a gene record exist
        # at all — `build_search` draws genes from `corpus.assertions`, so a
        # fixture without one leaves every gene assertion in this file vacuous.
        "assertions": (_assertion(),),
        "publications": (_publication(),),
        "phenotypes": (_phenotype(),),
    }
    payload.update(overrides)
    return Corpus(**payload)  # type: ignore[arg-type]


def _records(root: Path) -> list[dict[str, Any]]:
    """The published array, read back from the gzipped file that was written."""
    raw = (root / "search" / "index.json.gz").read_bytes()
    records: list[dict[str, Any]] = json.loads(gzip.decompress(raw))["records"]
    return records


def test_each_kind_is_published_as_one_record_of_what_a_visitor_types(tmp_path: Path) -> None:
    """The whole published contract for all three kinds, as one equality.

    `terms` is the haystack a client filters and `path` is what it fetches when a
    result is picked, so both are exact rather than approximate: `path` is pinned
    as a literal because the client must never build one itself, and the gene's
    is the slugged form, since `genes/HGNC:11604.json` is a URL carrying a
    character no consumer would leave unescaped.

    Term order is part of the file. The gene's aliases arrive unsorted and are
    published sorted, because an alias cell is an upstream dump whose order means
    nothing; the authors and the synonyms arrive unsorted and stay that way,
    because author order is data and synonym order is the curator's. All three
    lists would be identical under a mutation that sorted everything or sorted
    nothing, which is why each is deliberately out of order here.

    A30 asked whether a merged equality can hide a type error. Checked and it
    cannot here: every value in this payload is a string or a list of strings, so
    there is no boolean or number for `{"f": 0} == {"f": False}` to smuggle
    through. Should a count or a flag ever join a record, it needs an `is` beside
    this assertion.
    """
    emitter = Emitter(root=tmp_path)

    build_search(
        _corpus(),
        emitter,
        symbols={TBX5: "TBX5"},
        aliases={TBX5: ["T-box 5", "Holt-Oram syndrome 1"]},
    )

    assert _records(tmp_path) == [
        {
            "kind": "gene",
            "id": TBX5,
            "label": "TBX5",
            "terms": ["TBX5", TBX5, "Holt-Oram syndrome 1", "T-box 5"],
            "path": "genes/HGNC_11604.json",
        },
        {
            "kind": "publication",
            "id": "PMID:8988165",
            "label": "Mutations in human TBX5 cause Holt-Oram syndrome",
            "terms": [
                "Mutations in human TBX5 cause Holt-Oram syndrome",
                "Nature Genetics",
                "Li QY",
                "Basson CT",
            ],
            "path": "publications.json",
        },
        {
            "kind": "phenotype",
            "id": "HP:0001629",
            "label": "Ventricular septal defect",
            "terms": [
                "Ventricular septal defect",
                "VSD",
                "Interventricular septal defect",
            ],
            "path": "phenotypes.json",
        },
    ]


def test_each_kind_is_ordered_by_identifier(tmp_path: Path) -> None:
    """The published order of the one array every visitor downloads.

    Sized so that removing any `sorted()` fails on every run rather than on a
    fraction of them, which is what A33 found missing. The three sorts fail for
    two different reasons, and only one of them is about hash order:

    Publications and phenotypes come from tuples, so a build with the sort
    removed publishes the corpus order. That is deterministic — the corpus order
    here is deliberately neither lexical nor numeric, so the mutant fails on
    every run under every `PYTHONHASHSEED`, and three records of each is already
    conclusive. The PMIDs also pin *lexical* rather than numeric order, the
    decision A25 recorded: PMID:10 and PMID:100 lead PMID:9, and `derive.py`
    orders the same PMIDs the same way.

    Genes come from a set, so the mutant publishes hash order and the guard is
    only as strong as the fixture is wide. Measured on these exact ids in this
    exact corpus order, counting the seeds on which set order already equals
    sorted order and the mutant would therefore survive: 2 genes 96/200, 3 genes
    45/200, 4 genes 1/200, 5 genes 0/200, 6 genes 0/1000. Six is the size A10
    measured to 0/200 on the equivalent guard in `derive.py`, and it holds here.

    What A33 could not have fixed by widening the fixture is the shape of the
    test it describes. Comparing two builds inside one process can never fail:
    `PYTHONHASHSEED` is fixed for the life of the interpreter, so an unsorted set
    iterates identically in both builds however many genes it holds. Order has to
    be asserted against a literal, as below, not against a second run.
    """
    corpus = _corpus(
        assertions=tuple(
            _assertion(gene=gene, identifier=f"CHDA:AST:{index:07d}")
            for index, gene in enumerate(SIX_GENES, start=1)
        ),
        publications=(_publication("PMID:9"), _publication("PMID:10"), _publication("PMID:100")),
        phenotypes=(
            _phenotype("MONDO:0005453", "Atrioventricular septal defect"),
            _phenotype("HP:0001631", "Atrial septal defect"),
            _phenotype("HP:0001629", "Ventricular septal defect"),
        ),
    )
    emitter = Emitter(root=tmp_path)

    build_search(corpus, emitter, symbols={}, aliases={})

    assert [record["id"] for record in _records(tmp_path)] == [
        # Genes, by HGNC id. The corpus lists these in the reverse of this order.
        TBX5,
        "HGNC:12873",
        NKX2_5,
        GATA4,
        "HGNC:6188",
        "HGNC:7881",
        # Publications, lexically by PMID — neither numeric order nor the order
        # the corpus lists them in.
        "PMID:10",
        "PMID:100",
        "PMID:9",
        # Phenotypes, by identifier.
        "HP:0001629",
        "HP:0001631",
        "MONDO:0005453",
    ]


def test_a_gene_is_one_record_however_many_assertions_and_never_repeats_a_term(
    tmp_path: Path,
) -> None:
    """Two ways a term arrives twice, and the one that survives a green gate.

    TBX5 here is mirrored with its own symbol among its aliases, which nothing
    validates and no gate refuses — `mirrors/genes.tsv` is an upstream dump and
    the aliases cell is free text. GATA4 is the other way in: absent from the
    registry, its label falls back to the HGNC id, which `terms` already carries,
    so the id would appear twice. That second case is unreachable behind
    `build_site` — REF001 makes an asserted gene missing from the registry an
    error — so it is the first that justifies deduplicating and the second that
    comes along for free.

    Neither duplicate would change a query's outcome; both would publish bytes
    that say the atlas thinks "TBX5" is two different things you can search for.

    Also pinned here: two assertions about one gene are one result. `id` is the
    gene, not the assertion, so without the set a curator splitting one gene's
    evidence across two records would double it in every search.
    """
    corpus = _corpus(
        assertions=(
            _assertion(gene=TBX5, identifier="CHDA:AST:0000001"),
            _assertion(gene=TBX5, identifier="CHDA:AST:0000002"),
            _assertion(gene=GATA4, identifier="CHDA:AST:0000003"),
        ),
        publications=(),
        phenotypes=(),
    )
    emitter = Emitter(root=tmp_path)

    build_search(corpus, emitter, symbols={TBX5: "TBX5"}, aliases={TBX5: ["TBX5"]})

    records = _records(tmp_path)
    assert [record["id"] for record in records] == [TBX5, GATA4]
    assert records[0]["terms"] == ["TBX5", TBX5]
    assert records[1]["label"] == GATA4
    assert records[1]["terms"] == [GATA4]


def test_only_a_gene_carrying_an_assertion_is_searchable(tmp_path: Path) -> None:
    """Searchable means curated, and it has to, because `path` is a promise.

    A gene record advertises `genes/<id>.json`, and `build_genes` writes one
    bundle per gene `gene_facts` returns — the genes carrying an assertion.
    Drawing this list from the gene registry instead would publish search hits
    whose only action is a fetch that 404s, and nothing downstream could catch
    it: `Emitter.checksums` records what the build wrote, never what a payload
    promised.

    NKX2-5 below is mirrored but uncurated, the ordinary case — the registry
    holds every gene, the atlas curates a few. GATA4 carries functional evidence
    and no assertion, which is A29's case: that record reaches no published file
    at all, so its gene cannot be made findable here either. Both are excluded by
    one rule, and this is the only test in the file that would notice if either
    stopped being.
    """
    corpus = _corpus(
        assertions=(_assertion(gene=TBX5),),
        functional=(_functional(GATA4),),
        publications=(),
        phenotypes=(),
    )
    emitter = Emitter(root=tmp_path)

    build_search(
        corpus,
        emitter,
        symbols={TBX5: "TBX5", GATA4: "GATA4", NKX2_5: "NKX2-5"},
        aliases={NKX2_5: ["CSX"]},
    )

    assert [record["id"] for record in _records(tmp_path)] == [TBX5]


def test_an_empty_corpus_still_publishes_one_gzipped_index_at_a_fixed_path(
    tmp_path: Path,
) -> None:
    """A search box fetching this before the first assertion is curated must read
    an empty array rather than tell "nothing curated" from "wrong URL" by
    interpreting a 404.

    Pins the published contract as well as the degenerate input: one file, at
    that exact path, gzipped — `read_bytes` here is fed straight to
    `gzip.decompress`, so a build that wrote plain JSON fails rather than passing
    on a file that merely parses — under the key `records`, which is what every
    client indexes into.
    """
    emitter = Emitter(root=tmp_path)

    build_search(Corpus(root=Path(".")), emitter, symbols={}, aliases={})

    assert set(emitter.checksums) == {"search/index.json.gz"}
    written = (tmp_path / "search" / "index.json.gz").read_bytes()
    assert json.loads(gzip.decompress(written)) == {"records": []}
