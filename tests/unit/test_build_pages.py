# tests/unit/test_build_pages.py
"""What a gene page must say, especially when the atlas has curated nothing."""

from __future__ import annotations

import re
from dataclasses import replace
from datetime import date
from enum import StrEnum
from pathlib import Path

import pytest

from chd_atlas.build.derive import GeneFacts
from chd_atlas.build.emit import Emitter
from chd_atlas.build.pages import _SCOPE_RULE, build_gene_index_page, build_gene_pages
from chd_atlas.build.validity import GeneValidity, ValidityRecord
from chd_atlas.models.assertion import Evidence, InTextLocator, LesionAssertion
from chd_atlas.models.literature import Publication
from chd_atlas.vocab import (
    AtlasCuration,
    Classification,
    EvidenceClass,
    EvidenceStrength,
    LesionGroup,
    ValiditySource,
    ValidityState,
)

TBX5 = "HGNC:11604"
GATA4 = "HGNC:4173"


def _validity(
    disease_label: str = "Holt-Oram syndrome",
    report_url: str | None = "https://search.clinicalgenome.org/kb/gene-validity/x",
) -> GeneValidity:
    return GeneValidity(
        records=(
            ValidityRecord(
                source=ValiditySource.CLINGEN,
                classification=Classification.DEFINITIVE,
                classification_term="Definitive",
                disease="MONDO:0007732",
                disease_label=disease_label,
                moi="AD",
                sop="SOP11",
                gcep="Syndromic Disorders Gene Curation Expert Panel",
                report_url=report_url,
            ),
        ),
        state=ValidityState.EXPERT_CURATED,
        has_source_discordance=False,
    )


def _facts(
    gene: str,
    curation: AtlasCuration,
    groups: tuple[LesionGroup, ...] = (),
    functional_count: int = 0,
) -> GeneFacts:
    return GeneFacts(
        gene=gene,
        headline_confidence=Classification.DEFINITIVE,
        validity_state=ValidityState.EXPERT_CURATED,
        has_conflicting_evidence=False,
        has_source_discordance=False,
        lesion_groups=groups,
        confidence_by_lesion_group={group: Classification.DEFINITIVE for group in groups},
        conflicting_lesion_groups=(),
        evidence_counts={EvidenceClass.GENETIC_CASE: 1} if groups else {},
        assertion_count=1 if curation is AtlasCuration.CURATED else 0,
        functional_count=functional_count,
        publications=("PMID:8988165",) if curation is AtlasCuration.CURATED else (),
        atlas_curation=curation,
    )


@pytest.fixture
def facts_uncurated() -> dict[str, GeneFacts]:
    return {GATA4: _facts(GATA4, AtlasCuration.NOT_YET_CURATED)}


@pytest.fixture
def facts_curated() -> dict[str, GeneFacts]:
    return {TBX5: _facts(TBX5, AtlasCuration.CURATED, groups=(LesionGroup.SEPTAL,))}


@pytest.fixture
def facts_two() -> dict[str, GeneFacts]:
    return {
        TBX5: _facts(TBX5, AtlasCuration.CURATED, groups=(LesionGroup.SEPTAL,)),
        GATA4: _facts(GATA4, AtlasCuration.NOT_YET_CURATED),
    }


@pytest.fixture
def validity_two() -> dict[str, GeneValidity]:
    """Both genes of `facts_two`, each Definitive for a *different* in-scope disease.

    Different labels on purpose. The `definitive for` column and the rail row it
    mirrors are both per-gene lookups, and a fixture where every gene carried
    the same disease label could not tell a correct lookup from one wired to the
    first gene, to `next(iter(validity))`, or to a constant.
    """
    return {TBX5: _validity(), GATA4: _validity("structural congenital heart disease")}


@pytest.fixture
def assertions() -> dict[str, list[LesionAssertion]]:
    """One syndromic assertion, which is what `LesionAssertion` will accept.

    `syndromic_lists_extracardiac_features` refuses a `syndromic` record that
    names no extracardiac feature, so `extracardiac_features` is not optional
    here -- HP:0001156 is brachydactyly, the limb half of Holt-Oram.
    """
    return {
        TBX5: [
            LesionAssertion(
                id="CHDA:AST:0000001",
                gene=TBX5,
                phenotypes=["HP:0001631"],
                lesion_groups=[LesionGroup.SEPTAL],
                inheritance=["AD"],
                mechanism="haploinsufficiency",
                syndromic="syndromic",
                extracardiac_features=["HP:0001156"],
                curator="Enrique Audain",
                curated_on=date(2026, 7, 31),
                last_reviewed=date(2026, 7, 31),
                evidence=[
                    Evidence(
                        publication="PMID:8988165",
                        evidence_class=EvidenceClass.GENETIC_CASE,
                        locator=InTextLocator(element="Abstract"),
                        summary="A nonsense TBX5 mutation was found in affected members.",
                        strength=EvidenceStrength.STRONG,
                    )
                ],
            )
        ]
    }


@pytest.fixture
def publications() -> dict[str, Publication]:
    return {
        "PMID:8988165": Publication(
            id="PMID:8988165",
            title="Mutations in human TBX5 cause limb and cardiac malformation.",
            journal="Nature genetics",
            year=1997,
            authors=["Basson CT"],
            study_type="family_linkage",
            own_lab=False,
        )
    }


def _page(root: Path, name: str) -> str:
    return (root / "genes" / name).read_text(encoding="utf-8")


def _validity_table(page: str) -> str:
    """The mirrored-validity `<table>` alone, up to its closing tag.

    Scoped rather than checked page-wide because the rail's "this gene as JSON"
    link and the shell's `<nav>` are both `<a href=`, so "no anchor anywhere"
    could never pass on any page. What the report column must not contain is one.
    """
    start = page.index("<h2>Mirrored gene&ndash;disease validity</h2>")
    return page[start : page.index("</table>", start)]


def test_an_uncurated_gene_page_says_the_atlas_has_not_curated_it(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """22 of the 23 genes published today. The page must state it.

    An absent section is indistinguishable from "the atlas looked and found
    nothing", and a reader deciding what a gene means clinically must not have
    to infer which. The mirrored validity table is still present -- that is the
    content the page exists to show.

    Measured: deleting the `_NOT_CURATED` paragraph from `build_gene_pages`
    fails this test and no other in the suite. The rail's `atlas curation` row
    does not rescue it -- that renders `AtlasCuration.NOT_YET_CURATED.value`,
    which is `not_yet_curated` and matches no spelling asserted here.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts_uncurated,
        emitter,
        symbols={GATA4: "GATA4"},
        validity={GATA4: _validity()},
        assertions={},
        publications={},
    )

    page = _page(tmp_path, "HGNC_4173.html")
    assert "not yet curated" in page.lower()
    assert "Holt-Oram syndrome" in page
    assert "GATA4" in page


def test_a_curated_gene_page_carries_its_evidence_quote_and_pmid(
    tmp_path: Path,
    facts_curated: dict[str, GeneFacts],
    assertions: dict[str, list[LesionAssertion]],
    publications: dict[str, Publication],
) -> None:
    """The evidence a curator wrote, on a page, plus the way back to the JSON.

    The bundle link is asserted here rather than in its own test because it is
    the same interpolation the escaping tests already cover; what is pinned is
    the shape -- one `../` out of `genes/` and back in, since the page and its
    bundle are siblings and `paths.py` names both from the site root.

    The PubMed href is pinned for a stronger reason: nothing observed the URL
    `_pubmed` constructs. Asserting `"PMID:8988165" in page` passes on the link
    *text*, which is the identifier verbatim, so the `href` beside it was
    unmeasured -- dropping `.removeprefix("PMID:")` publishes
    `https://pubmed.ncbi.nlm.nih.gov/PMID:8988165/`, a 404 on every evidence row
    of every curated gene page, and the whole suite stayed green (measured
    2026-08-04: 616 passed). With the href asserted below, that mutant fails
    here and nowhere else: 1 failed, 621 passed.

    Asserted as the full attribute rather than as the bare number, because
    "8988165" is a substring of "PMID:8988165" and a substring check on the
    digits alone cannot fail.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts_curated,
        emitter,
        symbols={TBX5: "TBX5"},
        validity={TBX5: _validity()},
        assertions=assertions,
        publications=publications,
    )

    page = _page(tmp_path, "HGNC_11604.html")
    assert "PMID:8988165" in page
    assert 'href="https://pubmed.ncbi.nlm.nih.gov/8988165/"' in page
    assert "A nonsense TBX5 mutation was found in affected members." in page
    assert "Mutations in human TBX5 cause limb and cardiac malformation." in page
    assert "not yet curated" not in page.lower()
    assert 'href="../genes/HGNC_11604.json"' in page


def test_a_symbol_carrying_markup_is_escaped(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """Nothing curated or mirrored is reviewed the way page copy is.

    The symbol reaches the page twice -- the `<h1>` in the rail and the
    `<title>` through `document` -- and only the first is escaped by this
    module. Dropping that one `html.escape` fails the first assertion here
    while the second still passes on the title's own escaping, which is why
    both are asserted.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts_uncurated,
        emitter,
        symbols={GATA4: "<script>x</script>"},
        validity={GATA4: _validity()},
        assertions={},
        publications={},
    )

    page = _page(tmp_path, "HGNC_4173.html")
    assert "<script>x</script>" not in page
    assert "&lt;script&gt;" in page


@pytest.mark.parametrize(
    "report_url",
    ["Pseudoautosomal region, recessive", "javascript:alert(1)", "  https://evil.example"],
)
def test_a_non_http_report_url_never_becomes_a_link(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts], report_url: str
) -> None:
    """Escaping is not sanitising, and neither mirror schema constrains this column.

    The first case is real, not invented: `mirrors/gencc_submissions.tsv` carries
    one row -- SHOX, `MONDO:0009588`, Ambry Genetics -- whose `report_url` is the
    literal string `Pseudoautosomal region, recessive` (measured 2026-08-04, the
    only such row of the file). Langer mesomelic dysplasia is out of CHD scope
    today, so that row reaches no page and the defect is latent; the second case
    is the same column carrying a scheme that executes. `render.py` would escape
    either into a syntactically valid `href` and publish it.

    Both assertions matter. The first is the guard: no anchor at all for this
    row. The second is what stops a "fix" that merely escapes harder -- the raw
    string must not appear inside any attribute, quoted or not, because
    `href="javascript:alert(1)"` needs no unescaped character to run.

    The third case pins that leading whitespace does not smuggle a scheme past
    `str.startswith`; browsers strip it, `startswith` does not, so the em dash is
    the right answer there too.

    Measured 2026-08-04 with the guard removed (`Link(text="open", href=url) if
    url else _EM_DASH`): 3 failed, 619 passed -- one failure per case here, and
    no other test in the suite noticed.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts_uncurated,
        emitter,
        symbols={GATA4: "GATA4"},
        validity={GATA4: _validity(report_url=report_url)},
        assertions={},
        publications={},
    )

    page = _page(tmp_path, "HGNC_4173.html")
    assert "<a href=" not in _validity_table(page)
    assert report_url.strip() not in page


def test_an_uncurated_notice_never_denies_functional_records_the_rail_counts(
    tmp_path: Path,
) -> None:
    """The page must not contradict the column beside it.

    `atlas_curation` is `curated` iff `assertion_count > 0` (`derive.gene_facts`),
    but `_rail` counts `functional_count` from a population that needs no
    assertion at all, so a gene can carry curated functional records and still be
    `NOT_YET_CURATED`. The paragraph used to answer that with "not yet curated
    evidence for this gene ... nothing on this page is the atlas's own
    assessment", one column away from `functional records | 3`.

    `atlas_curation`'s meaning is deliberately unchanged -- it is tested,
    documented and published -- so what is pinned here is the wording: the notice
    names the *lesion assertion* as what is absent, narrows its denial to
    classifications, and acknowledges the records rather than denying them.

    Latent: the committed corpus has no `curation/functional/` directory, so
    every published gene counts zero functional records (measured 2026-08-04).
    Reachable only from this fixture.
    """
    facts = {GATA4: _facts(GATA4, AtlasCuration.NOT_YET_CURATED, functional_count=3)}
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts,
        emitter,
        symbols={GATA4: "GATA4"},
        validity={GATA4: _validity()},
        assertions={},
        publications={},
    )

    page = _page(tmp_path, "HGNC_4173.html")
    # The rail is counting them; the prose must agree.
    assert "<dt>functional records</dt><dd>3</dd>" in page
    assert "not yet curated" in page.lower()
    assert "not yet curated</strong> a lesion assertion" in page
    assert "3 functional records" in page
    # The two sentences that would deny it, in the forms they could come back in.
    assert "nothing on this page is the atlas's own assessment" not in page
    assert "not yet curated</strong> evidence for this gene" not in page


def test_a_mirrored_disease_label_carrying_markup_is_escaped(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """The validity table is filled from a mirror, not from curated text."""
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts_uncurated,
        emitter,
        symbols={GATA4: "GATA4"},
        validity={GATA4: _validity(disease_label="<img src=x onerror=alert(1)>")},
        assertions={},
        publications={},
    )

    assert "<img src=x" not in _page(tmp_path, "HGNC_4173.html")


@pytest.mark.parametrize(
    ("conflicting", "discordant", "shown", "hidden"),
    [
        (True, False, "conflicting evidence", "sources disagree"),
        (False, True, "sources disagree", "conflicting evidence"),
    ],
)
def test_a_contested_gene_is_never_chipped_as_settled(
    tmp_path: Path,
    facts_uncurated: dict[str, GeneFacts],
    conflicting: bool,
    discordant: bool,
    shown: str,
    hidden: str,
) -> None:
    """The green pill tracks the classification, not merely the fact of one.

    `chip-definitive` is the only chip class the stylesheet fills with the
    success colour, so keying it on `headline_confidence is not None` would
    paint a refuted gene green while its own label read "refuted" -- the
    display failure `vocab.strongest` is documented against.

    Parametrised over the two warning flags in opposite states because they are
    separate axes: `has_conflicting_evidence` is any mixed evidence,
    `has_source_discordance` is specifically ClinGen and GenCC disagreeing with
    each other. Each case is the unique killer of dropping one of the two chips
    -- measured 2026-08-04 by replacing each `if` with `if False:` in `_rail`:
    dropping "sources disagree" failed only the second case, dropping
    "conflicting evidence" only the first, and nothing else in this file
    noticed either. The `hidden` assertion is what catches the two flags being
    read into the wrong chip.

    Neither flag is set on any of the 23 genes published today, and no gene
    published today carries a headline other than `definitive` (measured
    2026-08-04 on a real build's `genes/index.json`), so every branch this test
    exercises is latent -- reachable only from a fixture until the mirrors
    change, and invisible to any test that renders the committed corpus.
    """
    facts = {
        GATA4: replace(
            facts_uncurated[GATA4],
            headline_confidence=Classification.REFUTED,
            has_conflicting_evidence=conflicting,
            has_source_discordance=discordant,
        )
    }
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts,
        emitter,
        symbols={GATA4: "GATA4"},
        validity={GATA4: _validity()},
        assertions={},
        publications={},
    )

    page = _page(tmp_path, "HGNC_4173.html")
    # The rendered chip, not the bare class name: `.chip-definitive` is also a
    # rule in the stylesheet every page inlines, so a substring check on the
    # class alone can never fail.
    assert 'class="chip chip-definitive"' not in page
    assert "refuted" in page
    assert shown in page
    assert hidden not in page


def test_every_gene_in_the_facts_gets_exactly_one_page(
    tmp_path: Path,
    facts_two: dict[str, GeneFacts],
    validity_two: dict[str, GeneValidity],
) -> None:
    """Including a gene with no symbol in the registry -- `symbols` is empty here."""
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts_two,
        emitter,
        symbols={},
        validity=validity_two,
        assertions={},
        publications={},
    )

    assert sorted(path.name for path in (tmp_path / "genes").iterdir()) == [
        "HGNC_11604.html",
        "HGNC_4173.html",
    ]


def test_the_browse_table_is_complete_before_any_script_runs(
    tmp_path: Path,
    facts_two: dict[str, GeneFacts],
    validity_two: dict[str, GeneValidity],
) -> None:
    """D29. The rows are rendered by the build, not fetched.

    The mutant this kills is a page that ships an empty tbody and populates it
    from genes/index.json at runtime: it looks identical in a browser and serves
    nothing to curl, to a crawler, or to a reader with scripts disabled.

    Measured 2026-08-04 by rendering `<tbody>` empty and appending the rows from
    an inline `fetch('index.json')` instead: 3 failed, 609 passed. This test
    failed first, on `"TBX5" in page` -- the symbol never reaches the markup at
    all, so the later `fetch(` and `<script>`-ordering assertions are not what
    catches this particular mutant; they catch the narrower one that renders the
    rows *and* then re-fetches. The other two failures were the two browse tests
    below, which is expected: an empty table has no order and no links.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_index_page(
        facts_two, emitter, symbols={TBX5: "TBX5", GATA4: "GATA4"}, validity=validity_two
    )

    page = _page(tmp_path, "index.html")
    assert "TBX5" in page and "GATA4" in page
    assert "fetch(" not in page
    assert page.index("HGNC_11604.html") < page.index("<script>")


def test_browse_rows_are_ordered_by_hgnc_id_against_a_literal(
    tmp_path: Path,
    facts_two: dict[str, GeneFacts],
    validity_two: dict[str, GeneValidity],
) -> None:
    """Asserted against a literal, not by building twice and comparing.

    PYTHONHASHSEED is fixed for the life of an interpreter, so a same-process
    build-twice comparison cannot catch a dropped sort at any fixture size.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_index_page(
        facts_two, emitter, symbols={TBX5: "TBX5", GATA4: "GATA4"}, validity=validity_two
    )

    page = _page(tmp_path, "index.html")
    assert page.index("HGNC_11604.html") < page.index("HGNC_4173.html")


def test_browse_facet_options_are_ordered_against_a_literal(
    tmp_path: Path, validity_two: dict[str, GeneValidity]
) -> None:
    """The `<option>` values are built from sets, so their order is seed-dependent.

    `build_gene_index_page` derives each facet's values with a set comprehension
    over `facts.values()`. Set iteration order for strings varies with
    `PYTHONHASHSEED`, and these options are part of the page's bytes and
    therefore of its checksum and its manifest entry -- so dropping the
    `sorted(values)` makes two builds of one commit differ. Measured 2026-08-04
    on the committed corpus, one process per seed: with the sort dropped,
    `genes/index.html` and `manifest.json` differ between `PYTHONHASHSEED=0` and
    `PYTHONHASHSEED=2`, and no other file does. The whole suite still passed --
    616 passed -- so nothing guarded this. With this test in place the mutant
    fails here and nowhere else: 1 failed, 621 passed.

    Asserted against a literal rather than by building twice: `PYTHONHASHSEED` is
    fixed for the life of an interpreter, so a same-process comparison cannot see
    this at any fixture size (CLAUDE.md §4.13).

    The fixture names **all nine** `LesionGroup` members on one gene, because
    this guard is probabilistic in the same way `test_cited_publications_come_back_sorted`
    is: *n* strings can land in an order that happens to equal sorted order, and
    the unsorted code then passes. Measured with the sort dropped, one process
    per seed, counting seeds on which the mutant survives:

    * the two-value `curation` facet the committed corpus would give
      (`curated`, `not_yet_curated`) -- **105/200**. Sizing this test to the real
      corpus would have produced a guard that waved the regression through on
      better than half of all runs, which §4.12 records as not being a guard.
    * the nine-member `lesion` facet below -- **0/200**, and **0/1000**.

    Nothing pins `PYTHONHASHSEED` in CI or in a `conftest.py`, so that difference
    is the whole value of the fixture. The `curation` assertion is kept beneath
    the `lesion` one anyway: it is the facet a reader actually uses to tell a
    curated gene from a mirrored one, and pinning its two values costs nothing.
    """
    facts = {
        TBX5: _facts(TBX5, AtlasCuration.CURATED, groups=tuple(LesionGroup)),
        GATA4: _facts(GATA4, AtlasCuration.NOT_YET_CURATED),
    }
    emitter = Emitter(root=tmp_path)
    build_gene_index_page(
        facts, emitter, symbols={TBX5: "TBX5", GATA4: "GATA4"}, validity=validity_two
    )

    facets = dict(
        re.findall(r'<select name="([^"]+)"[^>]*>(.*?)</select>', _page(tmp_path, "index.html"))
    )
    assert re.findall(r'<option value="([^"]*)">', facets["lesion"]) == [
        "",
        "anomalous_venous_return",
        "avsd",
        "complex",
        "conotruncal",
        "heterotaxy",
        "lvoto",
        "other",
        "rvoto",
        "septal",
    ]
    assert re.findall(r'<option value="([^"]*)">', facets["curation"]) == [
        "",
        "curated",
        "not_yet_curated",
    ]


def test_the_browse_page_says_whose_classification_the_confidence_column_carries(
    tmp_path: Path,
    facts_two: dict[str, GeneFacts],
    validity_two: dict[str, GeneValidity],
) -> None:
    """23 rows reading `definitive` under a column headed `confidence`, unattributed.

    Every gene page says "Every classification below is an upstream panel's or
    submitter's" and the landing page says "The atlas authors no validity
    classification of its own". This was the only page on the site where a
    panel's call could be read as the atlas's own, and the `atlas curation`
    column that answers the question was unexplained beside it.

    Both halves are asserted because they are separate claims: the attribution,
    and the pointer to the column that says whether the atlas has curated the
    gene at all. `_BROWSE_HEADERS` is checked in the same breath so the sentence
    cannot come to name a column the table stopped rendering.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_index_page(
        facts_two, emitter, symbols={TBX5: "TBX5", GATA4: "GATA4"}, validity=validity_two
    )

    page = _page(tmp_path, "index.html")
    assert "upstream panel's or submitter's" in page
    assert "the atlas authors no validity classification of its own" in page
    assert "<strong>atlas curation</strong> column" in page
    assert "<th>atlas curation</th>" in page
    assert "<th>confidence</th>" in page


def test_no_browse_row_states_a_bare_definitive_without_the_disease_it_is_for(
    tmp_path: Path,
    facts_two: dict[str, GeneFacts],
    validity_two: dict[str, GeneValidity],
) -> None:
    """`KMT2D | definitive` on a site called "CHD Atlas" is a wrong claim.

    ClinGen's assertion is `KMT2D -- Definitive for Kabuki syndrome 1`, made by
    the SCID-CID GCEP. Measured 2026-08-04 against the committed mirrors, the
    browse page rendered 23 rows all reading `definitive`, with no disease column
    at all: 10 of the 23 are qualified by a panel that is not the Congenital
    Heart Disease GCEP, and 10 are definitive for a disease whose label names no
    cardiac feature (Kabuki 1 and 2, Mowat-Wilson, CHARGE, Alagille, Ellis-van
    Creveld twice, Holt-Oram, TARP, NR2F2-related multiple congenital anomalies).

    Asserted per row rather than as "both labels appear somewhere on the page",
    because the failure this guards is a lookup wired to the wrong gene: the
    fixture gives the two genes two different diseases, and a page naming
    Holt-Oram on both rows would satisfy a page-wide substring check while
    telling a reader GATA4 is definitive for Holt-Oram syndrome.

    The disease cell is required to sit in the same `<tr>` as the `definitive`
    cell, so moving the column to the far right -- or into a second table -- fails
    here even though every string is still on the page.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_index_page(
        facts_two, emitter, symbols={TBX5: "TBX5", GATA4: "GATA4"}, validity=validity_two
    )

    page = _page(tmp_path, "index.html")
    assert "<th>definitive for</th>" in page
    rows = re.findall(r"<tr(?: data-[^>]*)?><td>(.*?)</tr>", page)
    assert len(rows) == 2
    expected = {TBX5: "Holt-Oram syndrome", GATA4: "structural congenital heart disease"}
    for row in rows:
        gene = next(g for g in expected if g in row)
        assert f"<td>definitive</td><td>{expected[gene]}</td>" in row, row
        # The other gene's disease must not be on this row at all.
        other = next(expected[g] for g in expected if g != gene)
        assert other not in row


def test_a_gene_definitive_for_two_in_scope_diseases_names_both_in_a_fixed_order(
    tmp_path: Path,
) -> None:
    """Latent today, and sorted so it stays deterministic when it is not.

    Measured 2026-08-04 against the committed mirrors: all 23 published genes
    carry exactly one in-scope ClinGen `Definitive` record naming exactly one
    distinct disease label, so nothing in the real corpus exercises this. It is a
    fixture, not a hypothetical: `curation/chd_scope.yaml` is an editorial file a
    curator widens, and one added MONDO term can give a gene a second qualifying
    record without any code changing.

    Two properties, and they are separate. That **both** labels are named --
    dropping one would publish half the reason the gene is on the site. And that
    they come back in sorted order -- `_definitive_diseases` de-duplicates through
    a `set`, whose iteration order for strings varies with `PYTHONHASHSEED`, so
    an unsorted return makes `genes/index.html` and its manifest checksum differ
    between two builds of one commit. Asserted against a literal in reverse
    alphabetical fixture order rather than by building twice: `PYTHONHASHSEED` is
    fixed for the life of an interpreter, so a same-process comparison cannot see
    a dropped sort at any fixture size (CLAUDE.md §4.13).

    Only a GenCC record calling the gene definitive is also present, and it is
    asserted absent: GenCC admits no gene to this population (D21), so naming its
    disease here would qualify the confidence with a disease that is not the one
    the gate turned on.
    """
    two = GeneValidity(
        records=(
            _validity("Zellweger syndrome").records[0],
            _validity("Alagille syndrome").records[0],
            ValidityRecord(
                source=ValiditySource.GENCC,
                classification=Classification.DEFINITIVE,
                classification_term="Definitive",
                disease="MONDO:0000001",
                disease_label="a submitter's disease",
                moi="AD",
                report_url=None,
                submitter="Ambry Genetics",
            ),
        ),
        state=ValidityState.EXPERT_CURATED,
        has_source_discordance=False,
    )
    facts = {TBX5: _facts(TBX5, AtlasCuration.CURATED, groups=(LesionGroup.SEPTAL,))}
    emitter = Emitter(root=tmp_path)
    build_gene_index_page(facts, emitter, symbols={TBX5: "TBX5"}, validity={TBX5: two})
    build_gene_pages(
        facts,
        emitter,
        symbols={TBX5: "TBX5"},
        validity={TBX5: two},
        assertions={},
        publications={},
    )

    browse = _page(tmp_path, "index.html")
    assert "<td>definitive</td><td>Alagille syndrome; Zellweger syndrome</td>" in browse
    assert "a submitter's disease" not in browse
    gene_page = _page(tmp_path, "HGNC_11604.html")
    assert "<dt>definitive for</dt><dd>Alagille syndrome; Zellweger syndrome</dd>" in gene_page


def test_both_page_kinds_state_the_rule_that_admits_a_gene_to_this_atlas(
    tmp_path: Path,
    facts_two: dict[str, GeneFacts],
    validity_two: dict[str, GeneValidity],
) -> None:
    """No page on the site said what the 23-gene set is. `docs/data-api.md` did.

    A reader of the HTML never sees that document, so the browse page and every
    gene page now carry the rule themselves: a gene is published when a ClinGen
    expert panel classifies it Definitive **for a disease in this atlas's CHD
    scope**, which is not the same as definitive for congenital heart disease.

    Both page kinds are checked from one constant, `pages._SCOPE_RULE`, because
    the rule is one editorial claim and two copies of it are two things that
    drift. What is asserted here is that the constant *reaches both kinds of
    page* -- the drift this cannot catch is the constant itself being reworded,
    which is what the phrase assertions below are for. They are deliberately not
    the whole string: this must fail when the sentence stops distinguishing the
    two claims, not merely when a comma moves.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_index_page(
        facts_two, emitter, symbols={TBX5: "TBX5", GATA4: "GATA4"}, validity=validity_two
    )
    build_gene_pages(
        facts_two,
        emitter,
        symbols={TBX5: "TBX5", GATA4: "GATA4"},
        validity=validity_two,
        assertions={},
        publications={},
    )

    for name in ("index.html", "HGNC_11604.html", "HGNC_4173.html"):
        page = _page(tmp_path, name)
        assert _SCOPE_RULE in page, f"{name} does not state the inclusion rule"
        assert "for a disease in this atlas's CHD scope" in page
        assert "not the same as definitive for congenital heart disease" in page


def test_a_gene_page_names_the_disease_beside_the_chip_that_says_definitive(
    tmp_path: Path,
    facts_uncurated: dict[str, GeneFacts],
) -> None:
    """The rail is the first screen, and the green chip on it said only `definitive`.

    `_rail` renders the symbol, the HGNC id and a `chip-definitive` pill, and on
    a viewport under 46rem `.layout` collapses to one column so those three fill
    the screen on their own. On a site titled "CHD Atlas" a green `definitive`
    with no disease reads as definitive for congenital heart disease; for CHD7 it
    means CHARGE syndrome, graded by the Hearing Loss GCEP.

    The order is asserted, not only the presence: the disease must come *after*
    the chip and before the rest of the rail's rows, so it is read as qualifying
    the chip rather than as one more fact three rows down.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts_uncurated,
        emitter,
        symbols={GATA4: "GATA4"},
        validity={GATA4: _validity("CHARGE syndrome")},
        assertions={},
        publications={},
    )

    page = _page(tmp_path, "HGNC_4173.html")
    assert "<dt>definitive for</dt><dd>CHARGE syndrome</dd>" in page
    assert page.index('chip chip-definitive">definitive') < page.index("CHARGE syndrome")
    assert page.index("CHARGE syndrome") < page.index("<dt>validity</dt>")


def test_every_browse_row_links_to_a_page_that_was_written(
    tmp_path: Path,
    facts_two: dict[str, GeneFacts],
    validity_two: dict[str, GeneValidity],
) -> None:
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts_two,
        emitter,
        symbols={},
        validity={TBX5: _validity(), GATA4: _validity()},
        assertions={},
        publications={},
    )
    build_gene_index_page(facts_two, emitter, symbols={}, validity=validity_two)

    page = _page(tmp_path, "index.html")
    for name in ("HGNC_11604.html", "HGNC_4173.html"):
        assert name in page
        assert f"genes/{name}" in emitter.checksums


def test_a_facet_option_value_is_escaped_like_every_other_published_string(
    tmp_path: Path, facts_two: dict[str, GeneFacts]
) -> None:
    """The one interpolation on the site where escaping rested on provenance.

    Review on #15 measured this: every other value on every page reaches markup
    through a `render.py` primitive, and `<option value="{value}">` did not. It
    was safe only because all 21 members of `LesionGroup`, `Classification`,
    `ValidityState` and `AtlasCuration` happen to be `[a-z_]`. Driving the loop
    with `"><script>alert(1)</script>` emitted a live `<script>` element into
    the page while the same build's table cells stayed correctly escaped.

    The hostile value arrives as a `StrEnum` member rather than a bare `str`,
    because the facet loop reads `.value` off each field -- so a plain string
    never reaches the interpolation and a test using one would fail for the
    wrong reason, which is what happened when this test was first written.
    """

    class Hostile(StrEnum):
        X = '"><script>alert(1)</script>'

    facts = dict(facts_two)
    facts[GATA4] = replace(facts[GATA4], validity_state=Hostile.X)  # type: ignore[arg-type]
    emitter = Emitter(root=tmp_path)

    build_gene_index_page(facts, emitter, symbols={}, validity={})

    page = _page(tmp_path, "index.html")
    assert "<script>alert(1)</script>" not in page
    assert "&lt;script&gt;alert(1)&lt;/script&gt;" in page


def test_every_facet_names_a_data_attribute_the_filter_script_reads(
    tmp_path: Path,
    facts_two: dict[str, GeneFacts],
    validity_two: dict[str, GeneValidity],
) -> None:
    """A facet whose name matches no `data-*` attribute filters nothing, silently.

    `FILTER_SCRIPT` resolves each facet by string concatenation --
    `row.getAttribute('data-' + select.name)` -- so a `<select name="lesions">`
    over rows carrying `data-lesion` reads `null`, matches nothing, and hides
    every row the moment a reader touches that facet. Nothing in the type system
    couples the two spellings and no other test in this file reads either, so
    this is the coupling. Not in the plan; added because the mismatch is exactly
    the silent failure this project is written against.

    The empty `<option value="">` is asserted for the same reason: the script
    reads an empty select value as "no filter", so a facet lacking one has no
    way back to showing every gene. `name="q"` likewise -- the script does
    `form.querySelector('[name=q]').value`, which throws on a missing input and
    leaves every filter dead.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_index_page(
        facts_two, emitter, symbols={TBX5: "TBX5", GATA4: "GATA4"}, validity=validity_two
    )

    page = _page(tmp_path, "index.html")
    assert 'name="q"' in page
    facets = re.findall(r'<select name="([^"]+)"[^>]*>(.*?)</select>', page)
    assert [name for name, _ in facets] == ["lesion", "confidence", "validity", "curation"]
    for _, options in facets:
        assert options.startswith('<option value="">')
    for row in re.findall(r"<tr((?: data-[^>]*)?)>", page):
        if row:
            assert set(re.findall(r"data-([a-z]+)=", row)) == {"search", *(n for n, _ in facets)}


def test_every_browse_control_is_named_for_a_screen_reader(
    tmp_path: Path,
    facts_two: dict[str, GeneFacts],
    validity_two: dict[str, GeneValidity],
) -> None:
    """Five form controls, five `aria-label`s. Raised by review on #14.

    Nothing else names any of them. The search box carries only a placeholder,
    which is announced inconsistently and disappears the moment a character is
    typed; each `<select>` carries only its own first option, and "any lesion"
    describes the control's *current value* rather than what it filters. A
    reader on a screen reader would meet five controls that announce nothing.

    Asserted against literals rather than against whatever the code emits,
    because a test that reads the label out of the page and compares it to
    itself would pass on an empty string. The labels differ deliberately from
    the neutral option text -- "Filter by lesion group" against "any lesion" --
    so an assertion that merely required *some* `aria-label` would not catch a
    label wired to the wrong facet.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_index_page(
        facts_two, emitter, symbols={TBX5: "TBX5", GATA4: "GATA4"}, validity=validity_two
    )

    page = _page(tmp_path, "index.html")
    assert 'aria-label="Search by gene symbol or HGNC id"' in page
    assert re.findall(r'<select name="([^"]+)" aria-label="([^"]+)"', page) == [
        ("lesion", "Filter by lesion group"),
        ("confidence", "Filter by mirrored confidence"),
        ("validity", "Filter by validity state"),
        ("curation", "Filter by atlas curation"),
    ]
