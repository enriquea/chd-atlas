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
from chd_atlas.build.concordance import cohort_families, evidence_axes
from chd_atlas.build.derive import GeneFacts
from chd_atlas.build.emit import Emitter
from chd_atlas.build.pages import (
    _EM_DASH,
    _SCOPE_RULE,
    build_gene_index_page,
    build_gene_pages,
)
from chd_atlas.build.paths import gene_page_path
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
        "count_unit": "individuals",
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


_OTHER_PUBLICATIONS = {
    "PMID:34324492": Publication(
        id="PMID:34324492",
        title="Integrative analysis of genomic variants",
        journal="PLoS genetics",
        year=2021,
        authors=["Audain E", "Wilsdon A"],
        study_type="meta_analysis",  # type: ignore[arg-type]
        own_lab=True,
        tests_reported=11515,
    ),
}


def _burden_page(
    tmp_path: Path,
    facts: dict[str, GeneFacts],
    rows: list[BurdenRow],
    name: str = "HGNC_4173.html",
    publications: dict[str, Publication] | None = None,
    families: tuple[frozenset[str], ...] | None = None,
    axes: tuple[tuple[str, str], ...] | None = None,
) -> str:
    """Render one gene page.

    `families` and `axes` default to the derivation the runner performs, so a
    test that does not care about the matrix still renders a real one rather
    than an empty section -- and a test that does care can pin them.
    """
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
        families=families if families is not None else cohort_families(rows),
        axes=axes if axes is not None else evidence_axes(rows),
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
    assert "<td>5 carriers / 1,471</td>" in page
    assert "<td>0 carriers / 45,082</td>" in page


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


def test_the_section_says_what_a_missing_row_means_without_asserting_why(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """The sentence that is the difference between a table and a trap -- and the
    mechanism it must NOT assert.

    The burden matrix is sparse, and a reader meeting a gap will read it as "not
    tested" unless told. So the page says an absent cell is one the study did not
    report.

    **It must not say why.** Until 2026-08-05 it asserted "no variant of that
    class was seen in either group, so there was nothing to compare" -- measured
    over PMID:42230622 alone, where 0 of 1,192 rows have no carrier on either
    side. That is false for the other two studies now rendering tables beside it:
    PMID:34324492 tests one consequence class by construction, and PMID:40127276
    observed 14,364 synonymous variants (its own Dataset S4) and published no
    synonymous row because its gene-level table covers only damaging classes. A
    reader told "no variant of that class was seen" would conclude that study
    found none.
    """
    page = _burden_page(tmp_path, facts_uncurated, [_burden_row()])
    section = page[page.index("Rare variant burden") :]

    assert "did not report" in section
    assert "does not always mean" in section
    # The refuted mechanism, in the exact words that were wrong.
    assert "no variant of that class was seen in either group" not in section


def test_the_synonymous_negative_control_is_explained_only_where_one_exists(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """Conditional, for the reason `_POOLING_NOTICE` and `_composite_note` are.

    A reader who does not know the synonymous row is a negative control reads a
    null result as a null finding rather than as the calibration check it is --
    so where one exists, the page says so.

    Where none exists, the page must not tell a reader to weigh one. Measured
    2026-08-05: 69 of the 290 published burden rows are synonymous and **every
    one comes from PMID:42230622**; the other two studies publish no negative
    control at all. Rendered unconditionally, this paragraph promised a check
    that two of the three tables on the page cannot offer.
    """
    without = _burden_page(tmp_path, facts_uncurated, [_burden_row(consequence_class="lof")])
    assert "negative control" not in without

    with_syn = _burden_page(
        tmp_path / "syn",
        facts_uncurated,
        [_burden_row(consequence_class="lof"), _burden_row(consequence_class="synonymous")],
    )
    assert "negative control" in with_syn
    assert "uncorrected" in with_syn


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
    # Scoped to the method line. This asserted `"uncorrected" in page` until
    # 2026-08-05 and passed on the word from the *synonymous* paragraph, which is
    # a different sentence entirely -- the trap
    # `test_the_consequence_column_is_headed_for_the_column_it_renders` records.
    # When that paragraph became conditional the assertion failed, which is the
    # only reason the mis-scoping was found.
    method = page[page.index('<p class="method">') :]
    method = method[: method.index("</p>")]

    assert "138,609" in method
    assert "3.6e-07" in method
    assert "judge the raw one against the whole scan" in method


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

    Two of the three studies published today are by this repository's author, and
    a reader being asked to weigh their numbers should not have to
    cross-reference a JSON payload to learn that.

    This said "every burden table on the site today comes from one study whose
    first author is this repository's author" until 2026-08-05. PMID:40127276
    (Sierant MC et al., `own_lab: false`) made that false for 15 of the 23
    published gene pages, 8 of which now render three studies with the
    declaration on two of them. The sentence was the written record of when the
    line applies, so a stale version of it is worse than none -- it describes a
    corpus that no longer exists and invites the conclusion that the disclosure
    is unconditional.
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


def test_the_pooling_notice_is_the_matrix_caption_and_no_longer_conditional(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """**Relocated, deliberately, and the condition went with it.**

    It used to render only where a gene carried two studies, because asserting
    "these cohorts overlap" above a single table sent a reader hunting for a
    second study that was not there.

    The evidence matrix changed what it qualifies. A summary that *counts
    datasets* is adjacent to pooling and invites exactly the arithmetic this
    sentence forbids -- so it is now the caption directly beneath the matrix,
    and it belongs there whenever the matrix does, one study or three. Folding
    it while introducing the thing it guards against would have been the worst
    available pairing.
    """
    one = _burden_page(tmp_path, facts_uncurated, [_burden_row()])
    assert "no pooled statistic across studies" in one
    # Immediately after the matrix, not four paragraphs above it.
    assert one.index("</table>") < one.index("no pooled statistic across studies")
    # **And it must be true of a one-study page.** It used to assert "these
    # cohorts overlap" in the present tense, which on this page sent a reader
    # hunting for a second study that is not there -- the defect the original
    # conditional existed to fix, reintroduced by relocating it. It is a
    # statement of policy now.
    assert "these cohorts overlap" not in one

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

    # Past the evidence matrix: it is now the first table in this section and
    # has its own `<thead>`. The comment above already records that
    # `page.index("<thead>")` finds the validity table; the matrix is the same
    # trap one layer in.
    section = section[section.index("</table>") :]
    header = section[section.index("<thead>") : section.index("</thead>")]
    assert '<th scope="col">consequence</th>' in header
    assert "variant class" not in header

    body = section[section.index("<tbody>") : section.index("</tbody>")]
    assert "loss-of-function" in body
    # The real `variant_class` is on the page exactly once, in the method line,
    # and never in this table.
    assert "SNVs and indels" not in body
    # Twice on the page, and neither is in this table: once as an evidence-matrix
    # column head, and once in the folded method line. It read `== 1` until the
    # matrix landed, when the column head made it 2 -- the number is asserted
    # rather than loosened so a third occurrence has to be a deliberate change.
    assert page.count("SNVs and indels") == 2


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


def test_every_count_cell_names_what_it_counted(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """`_count`'s rule, and `_effect`'s rule applied one column over.

    Until 2026-08-05 the header read "cases (carriers / n)" and every cell was a
    bare "5 / 1,471", which was true only because every study curated so far
    counted people. PMID:40127276's case-control arm counts *alleles* -- measured
    against its own Dataset S4 in `CountUnit`, where its D-Mis `Obs` tracks the
    variant count in 245 of 248 genes and the distinct-proband count in only
    235 -- so under the old header its rows would have claimed 21,768 people
    sequenced where 11,555 were.

    The three units render in one table here because that is the real case: one
    study contributes a de novo arm counting mutations against trios and a
    case-control arm counting alleles, into the same table, and a per-column
    header cannot label either.

    Scoped to the burden `<tbody>`. A page-wide assertion passes with the unit
    dropped from the cells, because `_BURDEN_PREAMBLE` names all three words
    higher up the page -- the trap `test_the_consequence_column_is_headed_for_
    the_column_it_renders` records one section above.
    """
    rows = [
        _burden_row(consequence_class="lof", count_unit="individuals"),
        _burden_row(consequence_class="missense_damaging", count_unit="alleles"),
        _burden_row(
            consequence_class="missense_all",
            count_unit="de_novo_mutations",
            origin="de_novo",
            comparator="mutation_model",
            n_control_carriers=None,
            n_controls=None,
            control_cohorts=(),
            expected_count=0.166,
            effect=None,
            effect_measure=None,
            effect_bound=None,
            ci_low=None,
        ),
    ]
    section = _burden_page(tmp_path, facts_uncurated, rows)
    section = section[section.index("Rare variant burden") :]
    # Past the evidence matrix, whose cells also carry counts.
    section = section[section.index("</table>") :]
    body = section[section.index("<tbody>") : section.index("</tbody>")]

    assert "<td>5 carriers / 1,471</td>" in body
    assert "<td>5 alleles / 1,471</td>" in body
    # The denominator of a de novo row is trios, and the word travels with it:
    # 5 mutations in 1,471 *families* is not 5 in 1,471 alleles.
    assert "<td>5 de novo / 1,471 trios</td>" in body

    # The control column takes the same unit, so an allele row does not describe
    # its cases in alleles and its controls in people.
    assert "<td>0 alleles / 45,082</td>" in body
    assert "<td>0 carriers / 45,082</td>" in body

    # No cell anywhere in the table is a bare "n / n": that spelling is what the
    # column meant before the unit existed, and it is the mutant this kills.
    assert "<td>5 / 1,471</td>" not in body
    assert "<td>0 / 45,082</td>" not in body


def test_a_count_unit_nothing_has_taught_the_renderer_still_names_itself(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """`_COUNT_UNIT_LABEL` is hand-maintained, so it can fall behind `CountUnit`.

    A member added to the enum and forgotten here must degrade to its raw token
    rather than to silence: "5 genomes / 1,471" is ugly and true, while "5 /
    1,471" is the exact false claim this column was added to prevent. Same shape
    as `_study_label` falling back to the bare PMID rather than an em dash.

    Reached only behind a bypassed gate -- `validate_table` refuses a value
    outside `CountUnit`, and `build_site` refuses the corpus -- which is why it
    is pinned here rather than left to be discovered.
    """
    page = _burden_page(tmp_path, facts_uncurated, [_burden_row(count_unit="genomes")])
    section = page[page.index("Rare variant burden") :]
    # Past the evidence matrix, which is now the first table in this section and
    # has a `<tbody>` of its own. Slicing on the first one silently measured the
    # matrix instead -- the same trap this file records twice already.
    section = section[section.index("</table>") :]
    body = section[section.index("<tbody>") : section.index("</tbody>")]

    assert "<td>5 genomes / 1,471</td>" in body
    assert "<td>5 / 1,471</td>" not in body


def test_the_composite_row_is_named_as_a_union_of_the_two_below_it(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """Three consequence rows must not read as three independent findings.

    PMID:40127276 reports `damaging (LOF + missense)` alongside its two
    components, and the composite is the analysis its 60 genes are defined by.
    Measured on CHD7: the damaging de novo row is 20 mutations and the two rows
    below it are 16 and 4 of *those same* 20, so a reader adding them up gets 40
    from 20 variants.

    The note is conditional, and both halves are asserted. A study reporting only
    components must not be told its rows decompose something that is not there --
    the defect `_POOLING_NOTICE` was made conditional for, where an
    unconditional sentence sent readers hunting for a second study that did not
    exist.
    """
    components = [
        _burden_row(consequence_class="lof"),
        _burden_row(consequence_class="missense_damaging"),
    ]
    composite = [_burden_row(consequence_class="damaging"), *components]

    with_composite = _burden_page(tmp_path, facts_uncurated, composite)
    assert "union</strong> of the loss-of-function" in with_composite

    without = _burden_page(tmp_path / "b", facts_uncurated, components)
    assert "union</strong> of the loss-of-function" not in without


def test_every_particular_caveat_survives_outside_the_fold(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """**The fold rule needs a test, not good intentions.**

    A caveat may fold if it is *general* -- true on every gene page, teaching
    how to read the table. It must stay visible if it is *particular* -- it
    fired because of this gene's data and changes what a reader concludes about
    this gene.

    This project has twice measured that a caveat nobody reads is nearly a
    caveat that does not exist, so folding a particular one later has to turn
    this red. Asserted against the page with every `<details>` block removed,
    which is what "unfolded" means to the reader who never clicks.

    The multiple-testing warning is in this list deliberately. It was the
    sharpest finding of the 2026-08-05 review -- 32 of 187 published rows clear
    0.05 and 3 survive Bonferroni over the study's own 138,609 comparisons --
    so `_method_line` returns it separately from the design clause, which does
    fold. Reuniting them would bury it one release after it was added.
    """
    rows = [
        _burden_row(consequence_class="lof"),
        _burden_row(
            study="PMID:34324492",
            consequence_class="lof",
            case_cohorts=("taa_cases",),
            method_note="Cases include 777 sporadic thoracic aortic aneurysm probands.",
        ),
    ]
    page = _burden_page(
        tmp_path,
        facts_uncurated,
        rows,
        publications={_PUBLICATION.id: _PUBLICATION, **_OTHER_PUBLICATIONS},
    )
    unfolded = re.sub(r"<details.*?</details>", "", page, flags=re.S)

    # The atlas's claim about itself.
    assert "computes none of them" in unfolded
    # The relocated pooling caption -- it qualifies the matrix directly above it.
    assert "count those children twice" in unfolded
    # A row's own method note: this is where the TAA contamination, which is not
    # congenital heart disease, reaches a reader at all.
    assert "thoracic aortic aneurysm" in unfolded
    # The multiple-testing denominator.
    assert "138,609" in unfolded


def test_the_general_reading_notes_are_folded_out_of_the_way(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """The other half of the rule, so "fold everything" fails as loudly as
    "fold nothing".

    These three paragraphs are identical on all 23 pages and teach how to read
    the table rather than saying anything about this gene. Nothing is lost by
    folding them: the unit word is still in every count cell, and the matrix
    renders `not tested` as its own state rather than as a gap.
    """
    page = _burden_page(tmp_path, facts_uncurated, [_burden_row(consequence_class="synonymous")])
    unfolded = re.sub(r"<details.*?</details>", "", page, flags=re.S)

    for general in ("every cell names what it counted", "negative control"):
        assert general in page, f"{general!r} must still be on the page"
        assert general not in unfolded, f"{general!r} must be behind the fold"


def test_the_matrix_names_an_untested_design_rather_than_leaving_a_blank(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """The holes are the point, and a blank cell is indistinguishable from a bug.

    Columns come from the whole corpus, so a design nobody ran for this gene
    renders as absent instead of vanishing. Nobody has published a CNV de novo
    analysis for any gene, and that quadrant is empty on every page: the shape
    of the literature is itself information.
    """
    page = _burden_page(
        tmp_path,
        facts_uncurated,
        [_burden_row(consequence_class="lof")],
        families=(frozenset({"PMID:42230622"}), frozenset({"PMID:34324492"})),
        axes=(("snv_indel", "control_cohort"), ("cnv_deletion", "control_cohort")),
    )
    section = page[page.index("Rare variant burden") :]
    matrix = section[section.index('<table class="matrix">') : section.index("</table>")]

    # **The cell reads short and explains itself on hover.** The full sentence
    # set the column width for the whole table and pushed every cell wide enough
    # to crowd its neighbours, so it moved to the `title` -- where it is still
    # reachable, which is the half of this that must not regress.
    assert ">not tested</span>" in matrix
    assert 'title="this dataset did not test this gene"' in matrix

    # Every *data* cell says what it is in words and is never simply blank. The
    # one bare `<td>` is the corner of a table headed on both axes, which
    # carries no value by construction.
    assert matrix.count("<td></td>") == 1
    assert "<td></td>" not in matrix[matrix.index("<tbody>") :]


def test_the_gene_and_its_symbol_both_link_to_the_page_and_the_strip_does_not(
    tmp_path: Path,
    facts_two: dict[str, GeneFacts],
    validity_two: dict[str, GeneValidity],
) -> None:
    """Where the affordance belongs, and where it does not.

    Both name cells link: `HGNC:11604` is the identifier a reader cites and
    `TBX5` is the name they scan for, and linking only the id put the affordance
    on the string nobody reaches for.

    **The strip is deliberately not a link.** It was one briefly, on the
    reasoning that it summarises the gene page's matrix -- but the row already
    links to that page twice, and a third link wrapped around a tally puts the
    affordance on a number rather than on a name. Each dot keeps its `title`,
    which is the detail the cell owes a reader.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_index_page(
        facts_two,
        emitter,
        symbols={GATA4: "GATA4", TBX5: "TBX5"},
        validity=validity_two,
        burden_counts={},
        concordance={
            gene: {"tested": 0, "enriched": 0, "corrected": 0, "families": []} for gene in facts_two
        },
    )
    page = _page(tmp_path, "index.html")
    href = f'href="../{gene_page_path(TBX5)}"'

    # Twice: once from the id, once from the symbol.
    assert page.count(href) == 2
    assert f"<a {href}>{TBX5}</a>" in page
    assert f"<a {href}>TBX5</a>" in page

    # And the strip carries no anchor at all.
    strip = page[page.index('<span class="strip">') :]
    strip = strip[: strip.index("</span></span>")]
    assert "<a " not in strip
    assert "strip-tally" in strip


def test_the_matrix_cell_drops_the_interval_but_never_the_measure(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """The one property `_effect_compact` trades width for, and the one it does not.

    The matrix cell is a summary, so it drops the confidence interval: that
    qualifies the estimate rather than naming it, it is in the row's own table a
    screen below, and carrying it ran a cell to 61 characters and set the column
    width for the whole grid.

    **The measure is not negotiable.** `_effect`'s rule is that no branch omits
    it, because an odds ratio of 3.1 and a de novo enrichment of 3.1 are
    different claims and a bare `3.1` equates them. A mutation matrix on
    2026-08-05 found that dropping it from the compact form survived the entire
    suite, which is the only reason this test exists.

    An unbounded effect still renders `∞` rather than a blank, for the reason
    `_effect` does: it is the strongest result in the data, and the one
    `allow_nan=False` refuses to publish as a number.
    """
    rows = [
        _burden_row(consequence_class="lof", effect=6.53, effect_bound=None, ci_low=2.1),
        _burden_row(consequence_class="missense_damaging"),  # the unbounded fixture
    ]
    page = _burden_page(tmp_path, facts_uncurated, rows)
    section = page[page.index("Rare variant burden") :]
    matrix = section[section.index('<table class="matrix">') : section.index("</table>")]

    # The measure travels with every number.
    assert "OR 6.53" in matrix or "OR ∞" in matrix
    assert (
        re.search(r'class="cell (corrected|nominal)"[^>]*>(?!OR|enrichment|rate)', matrix) is None
    )

    # The interval does not -- it is in the table below and on the cell's title.
    assert "95% CI" not in matrix
    assert "95% CI" in section[section.index("</table>") :]


def test_the_browse_headers_and_cells_line_up(
    tmp_path: Path,
    facts_two: dict[str, GeneFacts],
    validity_two: dict[str, GeneValidity],
) -> None:
    """Each column holds what its header says it holds.

    **`data_table` zips headers to cells by position and cannot detect a
    mismatch**, so inserting a header without moving its cell silently heads one
    column for a field two places away. That happened: "burden across studies"
    went in after `definitive for` while the strip stayed after `atlas
    curation`, and the whole suite passed -- 774 tests, none of which read a
    cell against its own header.

    This project already has `test_the_consequence_column_is_headed_for_the_
    column_it_renders` for exactly this failure on the burden table. The browse
    table had no equivalent, which is why the defect survived a build, a review
    of the rendered page, and a push.

    Asserted by pairing each header with a value only that column can produce,
    so a swap of any two columns fails rather than merely a shift of one.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_index_page(
        facts_two,
        emitter,
        symbols={GATA4: "GATA4", TBX5: "TBX5"},
        validity=validity_two,
        burden_counts={TBX5: 4},
        concordance={
            gene: {"tested": 2, "enriched": 1, "corrected": 1, "families": []} for gene in facts_two
        },
    )
    page = _page(tmp_path, "index.html")

    headers = [
        re.sub(r"<[^>]+>", "", cell)
        for cell in re.findall(r"<th[^>]*>.*?</th>", page[: page.index("</thead>")])
    ]
    row = next(r for r in re.findall(r"<tr[^>]*>(.*?)</tr>", page, re.S) if f">{TBX5}<" in r)
    cells = [re.sub(r"<[^>]+>", "", cell) for cell in re.findall(r"<td>(.*?)</td>", row, re.S)]

    assert len(headers) == len(cells), "a header was added without a cell, or the reverse"
    by_header = dict(zip(headers, cells, strict=True))

    assert by_header["gene"] == TBX5
    assert by_header["symbol"] == "TBX5"
    assert by_header["confidence"] == "definitive"
    # Only the strip renders a tally; only `validity` renders a validity state.
    assert "tested" in by_header["burden across studies"]
    assert by_header["validity"] == "expert_curated"
    assert by_header["atlas curation"] in {"curated", "not_yet_curated"}
    assert by_header["burden rows"] == "4"


def test_the_gene_page_matrix_carries_its_own_key_and_caption(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """**The browse page had these and the gene page did not.**

    That is the wrong way round. A reader arriving from a search lands on the
    gene page, and KDM6A's matrix there is entirely hollow -- four `not tested`
    and `no enrichment` cells -- beside a green `definitive` chip. Without the
    caption, that reads as the data contradicting the classification. It does
    not: KDM6A causes Kabuki syndrome, and burden tests at these cohort sizes
    routinely detect nothing for genes with overwhelming family evidence.

    Found by an adversarial review of the promotion diff, on the page where the
    omission would have done the damage rather than on the one where it would
    not have.
    """
    page = _burden_page(tmp_path, facts_uncurated, [_burden_row(consequence_class="lof")])
    section = page[page.index("Rare variant burden") :]

    # The caption that stops an empty matrix reading as a verdict.
    assert "not evidence against a gene" in section
    # And a key, so the four states are nameable rather than only coloured.
    assert "tested, no enrichment detected" in section
    assert "not tested by that dataset" in section
    # It sits with the matrix, not at the far end of the page.
    assert section.index("</table>") < section.index("not evidence against a gene")


def test_a_family_wise_corrected_p_is_not_labelled_q(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """`q` denotes a false-discovery quantity, and one of the two methods is not.

    The mirror carries `benjamini_hochberg`, which is an FDR, and
    `familywise_permutation`, which is not. Labelling both `q` made the matrix
    contradict the same page's own table, which names the second "family-wise" a
    screen below -- so the page said two different things about one number.

    Measured 2026-08-05: 29 cells rendered `q`, some over a family-wise p.
    """
    # Different *studies*, so the two land in different cells: one cell shows one
    # row, and putting both in the same one would render only the smaller p.
    fdr = _burden_row(
        consequence_class="lof",
        pvalue=0.001,
        pvalue_adjusted=0.01,
        pvalue_adjustment="benjamini_hochberg",
    )
    familywise = _burden_row(
        study="PMID:34324492",
        consequence_class="lof",
        pvalue=0.001,
        pvalue_adjusted=0.02,
        pvalue_adjustment="familywise_permutation",
        case_cohorts=("decipher",),
        control_cohorts=("gain_controls",),
    )
    page = _burden_page(
        tmp_path,
        facts_uncurated,
        [fdr, familywise],
        publications={_PUBLICATION.id: _PUBLICATION, **_OTHER_PUBLICATIONS},
    )
    section = page[page.index("Rare variant burden") :]
    matrix = section[section.index('<table class="matrix">') : section.index("</table>")]

    assert "q 0.01" in matrix
    assert "q 0.02" not in matrix
    assert "corrected p 0.02" in matrix


def test_the_browse_page_refuses_a_gene_it_has_no_concordance_for(
    tmp_path: Path,
    facts_two: dict[str, GeneFacts],
    validity_two: dict[str, GeneValidity],
) -> None:
    """The guard `build_genes` has, applied to the layer that lacked it.

    `bundles._concordance_for` raises on a published gene the mapping does not
    cover, because defaulting it would publish `"tested": 0` -- byte-identical to
    a measured "no study reported this gene". The browse page did
    `(concordance or {}).get(gene, {})` and rendered exactly that zero, silently.
    `runner.py` says the two layers cannot disagree; the guard existed on one of
    them.

    `None` remains allowed and means "this build has no burden data at all",
    which is what most tests of this page pass. A *mapping* that omits a
    published gene is the error: it means the caller derived it over a different
    population.

    Found by an adversarial review of the promotion diff, and it survived the
    mutation matrix until this test existed.
    """
    emitter = Emitter(root=tmp_path)

    with pytest.raises(KeyError, match="no concordance derived for published gene"):
        build_gene_index_page(
            facts_two,
            emitter,
            symbols={GATA4: "GATA4", TBX5: "TBX5"},
            validity=validity_two,
            burden_counts={},
            concordance={GATA4: {"tested": 0, "enriched": 0, "corrected": 0, "families": []}},
        )
