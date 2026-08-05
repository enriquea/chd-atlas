# tests/unit/test_build_pages.py
"""What a gene page must say, especially when the atlas has curated nothing."""

from __future__ import annotations

import re
from dataclasses import replace
from datetime import date
from enum import StrEnum
from pathlib import Path

import pytest

from chd_atlas.build.burden import BurdenRow
from chd_atlas.build.derive import GeneFacts
from chd_atlas.build.emit import Emitter
from chd_atlas.build.pages import (
    _EM_DASH,
    _SCOPE_RULE,
    build_gene_index_page,
    build_gene_pages,
)
from chd_atlas.build.validity import GeneValidity, ValidityRecord
from chd_atlas.models.assertion import Evidence, InTextLocator, LesionAssertion
from chd_atlas.models.cohort import Cohort
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
        burden={},
        cohorts={},
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
        burden={},
        cohorts={},
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
        burden={},
        cohorts={},
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
        burden={},
        cohorts={},
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
        burden={},
        cohorts={},
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
        burden={},
        cohorts={},
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
        burden={},
        cohorts={},
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
        burden={},
        cohorts={},
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
        facts_two,
        emitter,
        symbols={TBX5: "TBX5", GATA4: "GATA4"},
        validity=validity_two,
        burden_counts={},
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
        facts_two,
        emitter,
        symbols={TBX5: "TBX5", GATA4: "GATA4"},
        validity=validity_two,
        burden_counts={},
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
        facts,
        emitter,
        symbols={TBX5: "TBX5", GATA4: "GATA4"},
        validity=validity_two,
        burden_counts={},
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
        facts_two,
        emitter,
        symbols={TBX5: "TBX5", GATA4: "GATA4"},
        validity=validity_two,
        burden_counts={},
    )

    page = _page(tmp_path, "index.html")
    assert "upstream panel's or submitter's" in page
    assert "the atlas authors no validity classification of its own" in page
    assert "<strong>atlas curation</strong> column" in page
    assert '<th scope="col">atlas curation</th>' in page
    assert '<th scope="col">confidence</th>' in page


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
        facts_two,
        emitter,
        symbols={TBX5: "TBX5", GATA4: "GATA4"},
        validity=validity_two,
        burden_counts={},
    )

    page = _page(tmp_path, "index.html")
    assert '<th scope="col">definitive for</th>' in page
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
    build_gene_index_page(
        facts, emitter, symbols={TBX5: "TBX5"}, validity={TBX5: two}, burden_counts={}
    )
    build_gene_pages(
        facts,
        emitter,
        symbols={TBX5: "TBX5"},
        validity={TBX5: two},
        assertions={},
        publications={},
        burden={},
        cohorts={},
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
        facts_two,
        emitter,
        symbols={TBX5: "TBX5", GATA4: "GATA4"},
        validity=validity_two,
        burden_counts={},
    )
    build_gene_pages(
        facts_two,
        emitter,
        symbols={TBX5: "TBX5", GATA4: "GATA4"},
        validity=validity_two,
        assertions={},
        publications={},
        burden={},
        cohorts={},
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
        burden={},
        cohorts={},
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
        burden={},
        cohorts={},
    )
    build_gene_index_page(facts_two, emitter, symbols={}, validity=validity_two, burden_counts={})

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

    build_gene_index_page(facts, emitter, symbols={}, validity={}, burden_counts={})

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
        facts_two,
        emitter,
        symbols={TBX5: "TBX5", GATA4: "GATA4"},
        validity=validity_two,
        burden_counts={},
    )

    page = _page(tmp_path, "index.html")
    assert 'name="q"' in page
    facets = re.findall(r'<select name="([^"]+)"[^>]*>(.*?)</select>', page)
    assert [name for name, _ in facets] == [
        "lesion",
        "confidence",
        "validity",
        "curation",
        "burden",
    ]
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
        facts_two,
        emitter,
        symbols={TBX5: "TBX5", GATA4: "GATA4"},
        validity=validity_two,
        burden_counts={},
    )

    page = _page(tmp_path, "index.html")
    assert 'aria-label="Search by gene symbol or HGNC id"' in page
    assert re.findall(r'<select name="([^"]+)" aria-label="([^"]+)"', page) == [
        ("lesion", "Filter by lesion group"),
        ("confidence", "Filter by mirrored confidence"),
        ("validity", "Filter by validity state"),
        ("curation", "Filter by atlas curation"),
        ("burden", "Filter by whether burden evidence exists"),
    ]


def _burden_row(**overrides: object) -> BurdenRow:
    """TAB2's syndromic loss-of-function row: the unbounded case, from real data."""
    payload: dict[str, object] = {
        "study": "PMID:42230622",
        "gene": GATA4,
        "cohort_stratum": "syndromic",
        "lesion_group": None,
        "variant_class": "snv_indel",
        "consequence_class": "lof",
        "origin": "any",
        "maf_max": 0.001,
        "n_case_carriers": 5,
        "n_cases": 1471,
        "comparator": "control_cohort",
        "n_control_carriers": 0,
        "n_controls": 45082,
        "expected_count": None,
        "effect": None,
        "effect_measure": "odds_ratio",
        "effect_bound": "unbounded_above",
        "ci_low": 28.1,
        "ci_high": None,
        "pvalue": 3.13e-08,
        "pvalue_test": "fisher_exact",
        "pvalue_adjusted": None,
        "pvalue_adjustment": None,
        "case_cohorts": ("cnchd", "ddd"),
        "control_cohorts": ("ukbb",),
        "method_note": None,
        "source": "audain2026_sd3",
    }
    payload.update(overrides)
    return BurdenRow(**payload)  # type: ignore[arg-type]


_COHORTS = {
    "cnchd": Cohort(id="cnchd", name="German Competence Network", description="German registry."),
    "ddd": Cohort(
        id="ddd",
        name="Deciphering Developmental Disorders (DDD)",
        description="Ascertained on developmental disorder, so enriched for syndromic CHD.",
    ),
    "ukbb": Cohort(
        id="ukbb",
        name="UK Biobank",
        description="Adults recruited at 40-69; survivorship bias against severe CHD.",
    ),
    "pcgc": Cohort(id="pcgc", name="PCGC", description="Pediatric Cardiac Genomics Consortium."),
}

_PUBLICATION = Publication(
    id="PMID:42230622",
    title="Assessing the contribution of rare variants to congenital heart disease",
    journal="NPJ genomic medicine",
    year=2026,
    authors=["Audain E", "Hitz MP"],
    study_type="case_control",  # type: ignore[arg-type]
    own_lab=True,
    tests_reported=138609,
)


def _burden_page(
    tmp_path: Path,
    facts: dict[str, GeneFacts],
    rows: list[BurdenRow],
    name: str = "HGNC_4173.html",
    publications: dict[str, Publication] | None = None,
) -> str:
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts,
        emitter,
        symbols={GATA4: "GATA4", TBX5: "TBX5"},
        validity={GATA4: _validity(), TBX5: _validity()},
        assertions={},
        publications=publications or {_PUBLICATION.id: _PUBLICATION},
        burden={GATA4: rows},
        cohorts=_COHORTS,
    )
    return _page(tmp_path, name)


def test_an_effect_size_is_never_rendered_without_the_measure_that_names_it(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """The guard the single `effect` column was chosen against.

    One column carries odds ratios and de novo enrichments alike, which is what
    lets this schema absorb a fifth study without a migration -- and it is the
    one place two incomparable quantities could silently merge. An odds ratio of
    2.45 and a de novo enrichment of 2.45 are different claims, and a cell
    reading `2.45` under a header reading `effect` equates them.

    Both measures are rendered here from otherwise identical rows, so a
    `_effect` that dropped the label would produce two cells a reader could not
    tell apart -- and this test would see one string where it expects two.
    """
    page = _burden_page(
        tmp_path,
        facts_uncurated,
        [
            _burden_row(effect=2.45, effect_bound=None, ci_low=1.2, ci_high=8.1),
            _burden_row(
                cohort_stratum="all",
                comparator="mutation_model",
                n_control_carriers=None,
                n_controls=None,
                control_cohorts=(),
                expected_count=0.42,
                effect=2.45,
                effect_measure="enrichment_ratio",
                effect_bound=None,
                ci_low=1.2,
                ci_high=8.1,
                pvalue_test="poisson",
            ),
        ],
    )

    assert "OR 2.45 (95% CI 1.2–8.1)" in page
    assert "enrichment 2.45 (95% CI 1.2–8.1)" in page
    assert "<td>2.45</td>" not in page


def test_an_unbounded_odds_ratio_renders_its_lower_bound_rather_than_a_blank(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """34 rows of the committed mirror have no finite odds ratio, and they are
    the strongest results in the study.

    `Infinity` cannot be published -- `encode_json` sets `allow_nan=False`, and
    `JSON.parse` rejects it -- so the row carries a null effect and
    `effect_bound`. A page rendering that null as an em dash would throw away the
    finding: "at least 28.1" is what the study showed. The carrier counts beside
    it are what make it readable: 5 of 1,471 against 0 of 45,082.
    """
    page = _burden_page(tmp_path, facts_uncurated, [_burden_row()])

    assert "OR ∞ (95% CI 28.1–∞)" in page
    assert "<td>5 / 1,471</td>" in page
    assert "<td>0 / 45,082</td>" in page


def test_a_gene_with_no_burden_rows_gets_no_section_at_all(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """Unlike the validity table, whose header names the columns and whose
    emptiness is itself an answer, an empty burden table would read as "this gene
    was studied and nothing was found" -- a claim no study made.
    """
    assert "Rare variant burden" not in _burden_page(tmp_path, facts_uncurated, [])


def test_two_studies_sharing_a_cohort_are_declared_not_independent(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """Latent on the committed corpus, which carries one study.

    The atlas computes no pooled statistic because these cohorts overlap, and
    that decision is invisible to a reader comparing two tables by eye. The
    sentence is what makes it visible, and it names the collection rather than
    its id.
    """
    page = _burden_page(
        tmp_path,
        facts_uncurated,
        [_burden_row(), _burden_row(study="PMID:8988165", case_cohorts=("ddd", "pcgc"))],
    )

    assert "These two studies are not independent." in page
    assert "Deciphering Developmental Disorders (DDD)" in page


def test_the_section_says_what_a_missing_row_means_and_what_synonymous_is_for(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """Two sentences that are the difference between a table and a trap.

    The burden matrix is sparse -- 42 of the 145 genes in the committed mirror
    are missing at least one (stratum, consequence) cell -- and a reader meeting
    a gap will read it as "not tested" unless told. Zero of the 1,192 rows have
    no carrier in both groups, so an absent cell means exactly "nobody carried
    one", which is a much stronger statement.

    The synonymous sentence is the other half: a reader who does not know it is
    the negative control reads a null result as a null finding rather than as
    the calibration check it is.
    """
    page = _burden_page(tmp_path, facts_uncurated, [_burden_row()])

    assert "not a null result" in page
    assert "negative control" in page
    assert "uncorrected" in page


def test_the_page_names_how_many_tests_the_study_ran(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """The most serious finding of the 2026-08-05 review.

    Measured on the built site that day: of the 187 burden rows the 23 published
    gene pages render, **32 have p < 0.05 and 3 survive Bonferroni** over the
    study's own 138,609 comparisons. So 29 rows read as significant and are not,
    under a green `definitive` chip, on pages read by clinical geneticists.

    The atlas publishes no corrected p -- the supplement carries none, and
    computing one would be authoring a statistic (D12/D33). Naming the
    denominator is the honest alternative: it is a count the study made, and it
    is what lets a reader apply their own threshold. `0.05 / 138609` is stated
    as the arithmetic it is, not as a threshold the study endorsed.
    """
    page = _burden_page(tmp_path, facts_uncurated, [_burden_row()])

    assert "138,609" in page
    assert "uncorrected" in page
    assert "3.6e-07" in page


def test_the_cohort_caveats_a_curator_wrote_reach_the_reader(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """They reached no published byte until 2026-08-05.

    `curation/cohorts.yaml` records that UK Biobank's participants are adults
    while the CHD cases were largely enrolled in childhood -- a survivorship
    bias that inflates every odds ratio on every page -- and `models/cohort.py`
    says in as many words that `description` is where such caveats belong.
    `cohort_registry` returned `{id: name}` and dropped them.

    This is the project's characteristic failure in its purest form: the work
    was done, it was correct, it was validated, and no reader could reach it.
    """
    page = _burden_page(tmp_path, facts_uncurated, [_burden_row()])

    assert "survivorship bias against severe CHD" in page
    assert "enriched for syndromic CHD" in page
    # In a `<details>`, so the caveats are one click from every table rather than
    # repeated at full length above 23 of them -- and still in the document for
    # a crawler, a `curl`, and a reader with no JavaScript.
    assert "<details" in page


def test_the_page_says_what_was_counted_not_only_how_many(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """Four partition columns reached no page while populated on all 187 rows.

    A section headed "Rare variant burden" never said what *rare* meant, and
    `origin: any` -- a case-control count of variants regardless of inheritance
    -- was indistinguishable from a trio's de novo count. `vocab.VariantOrigin`
    warns about exactly that, and `tables.py` says the partition "exists to stop
    two incomparable rows from *looking* comparable".

    The origin label is asserted in full because the parenthetical is the whole
    point: "any inheritance" alone would still let a reader assume de novo.
    """
    page = _burden_page(tmp_path, facts_uncurated, [_burden_row()])

    assert "MAF below 0.001" in page
    assert "any inheritance (not a de novo test)" in page
    assert "SNVs and indels" in page
    assert "Fisher&#x27;s exact test" in page


def test_a_study_by_an_author_of_this_atlas_says_so(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """`own_lab` reached `publications.json` and zero HTML files.

    Every burden table on the site today comes from one study whose first author
    is this repository's author. A reader being asked to weigh those numbers
    should not have to cross-reference a JSON payload to learn that.
    """
    page = _burden_page(tmp_path, facts_uncurated, [_burden_row()])
    assert "this study is by an author of this atlas" in page

    outside = _PUBLICATION.model_copy(update={"own_lab": False})
    other = _burden_page(
        tmp_path, facts_uncurated, [_burden_row()], publications={outside.id: outside}
    )
    assert "by an author of this atlas" not in other


def test_a_de_novo_row_renders_the_expectation_it_was_compared_against(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """`_count(None, None)` rendered an em dash where the comparator belongs.

    On a `mutation_model` row the modelled expected count is the only thing the
    enrichment was computed against, and an em dash is indistinguishable from a
    control count nobody recorded. Latent until a trio study lands, which is
    precisely why it is pinned now.
    """
    page = _burden_page(
        tmp_path,
        facts_uncurated,
        [
            _burden_row(
                comparator="mutation_model",
                n_control_carriers=None,
                n_controls=None,
                control_cohorts=(),
                expected_count=0.42,
                effect=11.9,
                effect_measure="enrichment_ratio",
                effect_bound=None,
                ci_high=30.0,
                pvalue_test="poisson",
            )
        ],
    )

    assert "0.42 expected" in page
    assert "enrichment 11.9" in page


def test_a_method_note_is_rendered_as_the_footnote_the_schema_promises(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """`tables.py` called `method_note` "rendered verbatim as a row footnote"
    from the day the column existed, and nothing read it.

    The column exists for the CNV case, where "carrier" means different things
    in different papers, so the first curator to need it would have got a green
    build, a correct bundle and no page text.
    """
    page = _burden_page(
        tmp_path,
        facts_uncurated,
        [_burden_row(method_note="Carrier means any exonic overlap of the deletion.")],
    )

    assert "Carrier means any exonic overlap of the deletion." in page


def test_the_provenance_line_unions_every_row_rather_than_reading_the_first(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """`case_cohorts` is a per-row column and was reported per study block.

    A study whose strata drew on different collections published one row's
    cohorts as if they were all of them -- and `shared_cohorts` reads the same
    column with a *union* rule, so the "not independent" notice could name a
    cohort the provenance line did not. Both union now, so they cannot disagree.
    """
    page = _burden_page(
        tmp_path,
        facts_uncurated,
        [
            _burden_row(case_cohorts=("cnchd",)),
            _burden_row(cohort_stratum="all", case_cohorts=("ddd", "pcgc")),
        ],
    )

    # Scoped to the provenance paragraph. Asserted page-wide, this passed with
    # `_names` mutated to print bare ids, because `_cohort_notes` renders the
    # same names in the `<details>` block below the table.
    start = page.index('<p class="provenance">')
    line = page[start : page.index("</p>", start)]
    assert "German Competence Network" in line
    assert "Deciphering Developmental Disorders (DDD)" in line
    assert "PCGC" in line
    assert "cnchd" not in line


def test_a_study_heading_never_renders_a_nameless_author(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """`Publication.authors` constrains the list, not its elements.

    A whitespace-only first author passes every validator and rendered
    "    et al. 2026" -- a study heading naming nobody. `authors[0]` was read
    only by `search.py` before this section existed, so the page is where the
    pre-existing model weakness first became visible. Same shape as the blank
    `symbol` defect `runner.py::_cell` exists to catch.
    """
    nameless = _PUBLICATION.model_copy(update={"authors": ["   "]})
    page = _burden_page(
        tmp_path, facts_uncurated, [_burden_row()], publications={nameless.id: nameless}
    )

    assert "et al." not in page
    assert "PMID:42230622</a>" in page


def test_the_pooling_warning_is_absent_when_there_is_only_one_study(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """It asserted "these cohorts overlap" in the present tense on a page showing
    one table, sending a reader to look for a second study that is not there and
    undermining the caveats that *are* live. The corpus carries one study today,
    so this is the branch every real page takes.
    """
    one = _burden_page(tmp_path, facts_uncurated, [_burden_row()])
    assert "no pooled statistic across studies" not in one

    two = _burden_page(
        tmp_path,
        facts_uncurated,
        [_burden_row(), _burden_row(study="PMID:8988165", case_cohorts=("ddd", "pcgc"))],
    )
    assert "no pooled statistic across studies" in two


def test_the_consequence_column_is_headed_for_the_column_it_renders(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """The header and the cell must name the same field.

    Raised by review on #17. `_BURDEN_HEADERS[1]` read "variant class" while the
    cell rendered `consequence_class`, and the method line above the table
    rendered the real `variant_class` -- so one phrase meant two things on one
    page and was wrong in the more prominent of them. The rename that introduced
    it was itself a fix for an ambiguous header.

    Asserted by pairing the header with a cell whose value could only have come
    from one field: `loss-of-function` is a `ConsequenceClass` member and
    `SNVs and indels` is a `VariantClass` label, so a header claiming the latter
    over a cell holding the former fails here.
    """
    page = _burden_page(tmp_path, facts_uncurated, [_burden_row()])
    # Scoped to the burden section: `page.index("<thead>")` finds the mirrored
    # validity table, which is rendered first. Same trap as the provenance test.
    section = page[page.index("Rare variant burden") :]

    header = section[section.index("<thead>") : section.index("</thead>")]
    assert '<th scope="col">consequence</th>' in header
    assert "variant class" not in header

    body = section[section.index("<tbody>") : section.index("</tbody>")]
    assert "loss-of-function" in body
    # The real `variant_class` is on the page exactly once, in the method line,
    # and never in this table.
    assert "SNVs and indels" not in body
    assert page.count("SNVs and indels") == 1


def test_the_browse_page_shows_how_much_burden_evidence_each_gene_has(
    tmp_path: Path,
    facts_two: dict[str, GeneFacts],
    validity_two: dict[str, GeneValidity],
) -> None:
    """`genes/index.json` published `burden_row_count`; the page showed nothing.

    Measured on the built site 2026-08-05: `genes/index.html` contained the word
    "burden" zero times -- no column, no filter, no count -- while every row of
    the payload behind it carried the number. A reader scanning 23 rows that all
    read `definitive` could not tell that four of them (ISL1, NR2F2, RBM10,
    SMAD2) carry no loss-of-function burden evidence at all.

    The em dash for a gene with none is deliberate and matches the `lesion
    groups` column beside it: `0` reads as a measured zero, while an em dash
    reads as "nothing here", which is what an absent row means.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_index_page(
        facts_two,
        emitter,
        symbols={TBX5: "TBX5", GATA4: "GATA4"},
        validity=validity_two,
        burden_counts={TBX5: 9},
    )

    page = _page(tmp_path, "index.html")

    assert '<th scope="col">burden rows</th>' in page
    rows = {
        match.group(1): match.group(0)
        for match in re.finditer(r'<tr data-search="(hgnc:\d+)[^"]*".*?</tr>', page)
    }
    assert "<td>9</td>" in rows["hgnc:11604"]
    assert f"<td>{_EM_DASH}</td>" in rows["hgnc:4173"]

    # The facet and the attribute must spell the same word; see
    # `test_every_facet_names_a_data_attribute_the_filter_script_reads`.
    assert 'data-burden="yes"' in rows["hgnc:11604"]
    assert 'data-burden="no"' in rows["hgnc:4173"]


def test_the_rail_counts_burden_rows_and_names_the_publications_it_counts(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """Two rail defects, both raised by review 2026-08-05.

    The rail read `publications | 0` on 22 of the 23 published genes while the
    burden section one column away cited a linked, PubMed-referenced study. The
    count was right -- `fact.publications` is what this atlas's *own* assertion
    evidence cites -- and the label was wrong, so the label changed rather than
    the number.

    And the rail carried no burden count at all, though `burden_row_count`
    exists and is published in every browse row.
    """
    page = _burden_page(
        tmp_path, facts_uncurated, [_burden_row(), _burden_row(cohort_stratum="all")]
    )

    assert "<dt>curated publications</dt><dd>0</dd>" in page
    assert "<dt>burden rows</dt><dd>2</dd>" in page
    assert "<dt>publications</dt>" not in page


def test_a_published_correction_is_rendered_and_names_its_method(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """The column that answers the sharpest finding of the 2026-08-05 review.

    That review measured 187 uncorrected p-values on gene pages with nothing to
    judge them against. The atlas still computes no correction -- that would be
    authoring a statistic (D12/D33) -- but PMID:34324492 publishes one, and it
    **changes the conclusion**: CHD7 has a raw permutation p of 0.0068 and a
    family-wise corrected p of 0.991. Without this column the page shows a
    number that reads as significant for a gene the study found nothing for.

    The method is named in the cell because "0.991" alone does not say what it
    was corrected against, and a family-wise permutation correction and a
    Bonferroni factor are different claims.
    """
    page = _burden_page(
        tmp_path,
        facts_uncurated,
        [
            _burden_row(
                variant_class="cnv_deletion",
                effect=None,
                effect_measure=None,
                effect_bound=None,
                ci_low=None,
                ci_high=None,
                pvalue=0.0068,
                pvalue_test="permutation",
                pvalue_adjusted=0.991,
                pvalue_adjustment="familywise_permutation",
            )
        ],
    )

    assert '<th scope="col">corrected p</th>' in page
    assert "<td>0.991 (family-wise)</td>" in page
    assert "<td>0.0068</td>" in page


def test_the_method_line_tells_a_corrected_study_apart_from_an_uncorrected_one(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """A sentence true of one study is false of another.

    PMID:42230622 publishes no correction, so the reader needs the denominator
    and the advice to judge against the whole scan. PMID:34324492 corrects every
    row, so that same advice would point them at the column they should *not* be
    reading. The block says whichever is true of its own rows.
    """
    corrected = _burden_page(
        tmp_path,
        facts_uncurated,
        [_burden_row(pvalue_adjusted=0.991, pvalue_adjustment="bonferroni")],
    )
    assert "read it rather than the raw p" in corrected
    assert "judge the raw one against the whole scan" not in corrected

    plain = _burden_page(tmp_path, facts_uncurated, [_burden_row()])
    assert "judge the raw one against the whole scan" in plain
    assert "read it rather than the raw p" not in plain
