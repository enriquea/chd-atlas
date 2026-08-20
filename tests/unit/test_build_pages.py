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
    _phase_sentence,
    build_gene_index_page,
    build_gene_pages,
)
from chd_atlas.build.paths import gene_page_path
from chd_atlas.build.profiles import (
    EMPTY_EXPRESSION_PROFILE,
    DatasetProfileEntry,
    ExpressionProfile,
    PhaseInfo,
    Placement,
    Specificity,
    StageProfileEntry,
    TissueProfileEntry,
)
from chd_atlas.build.validity import GeneValidity, ValidityRecord
from chd_atlas.models.assertion import Evidence, InTextLocator, LesionAssertion
from chd_atlas.models.cohort import Cohort
from chd_atlas.models.dataset import Dataset, Stage
from chd_atlas.models.literature import Publication
from chd_atlas.models.phases import CardiacPhase, CardiacPhaseFile, EndBasis
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
    no_association_reported_by: tuple[str, ...] = (),
) -> GeneFacts:
    """A published gene's facts. The third axis defaults off and is opt-in.

    `no_association_reported_by` carries both halves of the pair: the flag is
    derived from it here rather than passed separately, so a fixture cannot
    construct the impossible state -- an authority named while the flag reads
    false -- that `derive.gene_facts` refuses to produce.
    """
    return GeneFacts(
        gene=gene,
        headline_confidence=Classification.DEFINITIVE,
        validity_state=ValidityState.EXPERT_CURATED,
        has_conflicting_evidence=False,
        has_source_discordance=False,
        has_no_association_report=bool(no_association_reported_by),
        no_association_reported_by=no_association_reported_by,
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
    # Scoped to the lesion-assertion axis, like line 399's own check, rather
    # than the bare phrase: this gene carries no profiles data either (the
    # call above passes no `profiles=`), so its page also carries the
    # unrelated "not yet curated a developmental expression profile" notice
    # -- a second, independent axis legitimately reusing the same template
    # `_NOT_CURATED` established. A bare page-wide check cannot tell the two
    # apart; CLAUDE.md section 4.19 is this exact failure shape.
    assert "not yet curated</strong> a lesion assertion" not in page
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

        `chip-definitive` is the top of the stylesheet's grade ramp and the only
        class filled with the full success colour, so keying it on
        `headline_confidence is not None` would paint a refuted gene green while its
        own label read "refuted" -- the display failure `vocab.strongest` is
        documented against. Since 2026-08-06 the class comes from
        `render.GRADE_CHIP_KIND`, keyed on the grade's own value, and `refuted` maps
        to `warn`: the assertion below is that the rail's chip is the caution colour,
        not merely that it is not green.

        **The slice is the rail, not the page.** `render.grade_legend` renders one
        chip of every class -- including `chip-definitive` -- inside a `<details>` on
        every page, so a page-wide `not in` can never fail again and passed only
        because the legend did not exist when it was written. CLAUDE.md section 4.19
        is this exact shape: a page-wide assertion is not a section assertion.

        Parametrised over the two warning flags in opposite states because they are
        separate axes: `has_conflicting_evidence` is any mixed evidence,
        `has_source_discordance` is specifically ClinGen and GenCC disagreeing with
        each other. Each case is the unique killer of dropping one of the two chips
        -- measured 2026-08-04 by replacing each `if` with `if False:` in `_rail`:
        dropping "sources disagree" failed only the second case, dropping
        "conflicting evidence" only the first, and nothing else in this file
        noticed either. The `hidden` assertion is what catches the two flags being
        read into the wrong chip.

    Neither flag is set on any of the 92 genes published today (measured
        2026-08-06 on a real build's `genes/index.json`), and no published gene
        carries a contested headline -- the floor excludes those rungs and the
        ClinGen veto excludes the rest -- so every branch this test exercises is
        latent, reachable only from a fixture until the mirrors change and invisible
        to any test that renders the committed corpus.
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
    rail = re.search(r'<aside class="rail">.*?</aside>', page, re.S)
    assert rail is not None, "the gene page no longer renders a rail"
    # The rendered chip, not the bare class name: `.chip-definitive` is also a
    # rule in the stylesheet every page inlines, so a substring check on the
    # class alone can never fail.
    assert 'class="chip chip-definitive"' not in rail.group()
    assert 'class="chip chip-warn">refuted</span>' in rail.group()
    assert shown in rail.group()
    assert hidden not in rail.group()


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

    The column is headed `graded for` since 2026-08-06, not `definitive for`: the
    widened gate publishes 43 `limited`, 9 `moderate` and 1 `strong` gene beside
    the 23 `definitive` ones, and a column head naming one rung above a cell
    naming another is a header that contradicts its own data.

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
    assert '<th scope="col">graded for</th>' in page
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
    they come back in sorted order -- `_graded_diseases` de-duplicates through
    a `set`, whose iteration order for strings varies with `PYTHONHASHSEED`, so
    an unsorted return makes `genes/index.html` and its manifest checksum differ
    between two builds of one commit. Asserted against a literal in reverse
    alphabetical fixture order rather than by building twice: `PYTHONHASHSEED` is
    fixed for the life of an interpreter, so a same-process comparison cannot see
    a dropped sort at any fixture size (CLAUDE.md §4.13).

    Only a GenCC record calling the gene definitive is also present, and it is
    asserted absent: `_graded_diseases` names the diseases the *admitting panel*
    graded, and this gene was admitted on a ClinGen record, so naming a
    submitter's disease here would qualify the confidence with a disease that is
    not the one the gate turned on.

    That rationale read "GenCC admits no gene to this population (D21)" until
    2026-08-06, which the widening made false -- GenCC agreement now admits 16 of
    the 92. The assertion it justifies is unchanged and still correct; only the
    reason was stale. A GenCC-admitted gene has no headline grade at all, so
    `_graded_diseases` returns `()` for it and never reaches this branch.
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
    """No page on the site said what the published set is. `docs/data-api.md` did.

    A reader of the HTML never sees that document, so the browse page and every
    gene page now carry the rule themselves: a gene is published when a ClinGen
    expert panel classifies it Limited or better for a disease an external
    authority treats as CHD, or when two GenCC submitters agree and no panel
    disputes it. Neither is the same as definitive for congenital heart disease.

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
        # **The published provenance claim, corrected 2026-08-06.** Every page
        # used to say a gene is definitive "for a disease in this atlas's CHD
        # scope", which told 24 published pages that this atlas decides what
        # counts as congenital heart disease. It does not, and its owner -- not
        # a clinician -- said so: scope is now a selection among terms external
        # authorities already use, and `curation/chd_scope.yaml` names the
        # authority for every term.
        #
        # The negative assertion is the load-bearing half. Without it a future
        # edit could restore the old sentence beside the new one and satisfy
        # every positive check here.
        assert "this atlas's CHD scope" not in page
        assert "an external authority treats as congenital heart disease" in page
        assert "No disease is in scope on this atlas's own judgement" in page
        assert "not the same as definitive for congenital heart disease" in page
        # **Both warrants, or the sentence describes a gate the site does not
        # have.** Stating only the ClinGen half would tell a reader that every
        # gene here was graded by an expert panel; 16 of the 92 were not, and
        # they are the ones whose provenance most needs stating.
        #
        # **Sliced to the scope-rule element, and that is what makes these
        # assertions able to fail.** Page-wide they could not: the grade key's
        # `ungraded` gloss carries the phrase "two or more Gene Curation
        # Coalition submitters independently assert it" further down every page,
        # so rewriting `_SCOPE_RULE` to state a **one**-submitter gate -- the
        # site describing a gate that admits on a single laboratory's word --
        # **survived all 796 tests** when measured 2026-08-06. The assertion
        # above it, `_SCOPE_RULE in page`, cannot catch a rewording either: it
        # compares the module's constant to itself and passes on any wording at
        # all. CLAUDE.md section 4.19.
        rule = re.search(r'<p class="scope-rule">.*?</p>', page, re.S)
        assert rule, f"{name} carries no scope-rule element"
        sentence = rule.group(0)
        assert "Limited or better" in sentence
        assert "two or more Gene Curation Coalition submitters" in sentence
        assert "one or more Gene Curation Coalition submitters" not in sentence
        assert "no ClinGen panel disputes it" in sentence


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


def test_the_browse_strips_four_glyphs_are_all_named_in_its_legend(
    tmp_path: Path,
    facts_two: dict[str, GeneFacts],
    validity_two: dict[str, GeneValidity],
) -> None:
    """A key naming three of four states leaves the fourth as an unexplained mark.

    Measured by mutation 2026-08-05: replacing every label in `_STRIP_LEGEND`
    with a single letter survived this whole file. The gene page's matrix legend
    was guarded and the browse page's strip legend was not, which is the same
    asymmetry — and the same direction — as the missing matrix caption an
    adversarial review caught one release earlier.

    **Literals, and each caption asserted adjacent to its own glyph.** The
    original version of this test read the captions out of
    `EVIDENCE_STATE_LABELS` and checked membership, which passes under any
    permutation of that constant because the needles come from the thing being
    permuted. Measured 2026-08-06: swapping two entries survived all 787 tests
    and published a browse page whose key captions the untested glyph "tested,
    no enrichment detected" -- the `not_tested`-reads-as-`no_enrichment`
    collapse this project forbids, on a green build.
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
    legend = page[page.index("</table>") :]

    for css, caption in (
        ("dot full", "enriched, and survives that study&#x27;s own correction"),
        ("dot half", "enriched nominally, or no correction published"),
        ("dot none", "tested, no enrichment detected"),
        ("dot untested", "not tested by that dataset"),
    ):
        pair = re.search(rf'<span class="{re.escape(css)}">[^<]*</span>\s*([^<]+)<', legend)
        assert pair, f"the strip legend no longer draws {css!r} with a caption"
        assert pair.group(1).strip() == caption, (
            f"{css!r} is captioned {pair.group(1).strip()!r}, expected {caption!r}"
        )


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
    # Each swatch beside its own caption, not merely both strings present: a
    # membership check passes with the two swapped, and a swapped key on this
    # page tells a reader an untested cell was tested and found nothing. No
    # assertion named `cell-key` existed at all until review 2026-08-06.
    for css, caption in (
        ("cell-key no-enrichment", "tested, no enrichment detected"),
        ("cell-key not-tested", "not tested by that dataset"),
    ):
        pair = re.search(rf'<span class="{re.escape(css)}">[^<]*</span>\s*([^<]+)<', section)
        assert pair, f"the matrix key no longer draws {css!r} with a caption"
        assert pair.group(1).strip() == caption, (
            f"{css!r} is captioned {pair.group(1).strip()!r}, expected {caption!r}"
        )
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


def _limited(gene: str) -> tuple[GeneFacts, GeneValidity]:
    """A gene whose only ClinGen record grades it `Limited`.

    **Not `Definitive`, and that is the whole point of the fixture.** Every other
    gene fixture in this file is graded `Definitive`, so "the headline's grade"
    and "`DEFINITIVE`" are the same value and a lookup keyed on the wrong one is
    invisible -- CLAUDE.md section 4.30, third recurrence. Measured 2026-08-06:
    three mutants survived the whole suite before this existed, including
    `_graded_diseases` keyed back on `DEFINITIVE` and the rail label pinned to
    the literal `definitive for`.

    `confidence_by_lesion_group` is emptied because `_facts` derives it from
    `groups`, which is empty here; leaving `_facts`' default would publish a map
    keyed on lesion groups the gene does not declare.
    """
    limited = ValidityRecord(
        source=ValiditySource.CLINGEN,
        classification=Classification.LIMITED,
        classification_term="Limited",
        disease="MONDO:0007732",
        disease_label="atrial septal defect",
        moi="AD",
        sop="SOP11",
        gcep="Congenital Heart Disease Gene Curation Expert Panel",
        report_url=None,
    )
    return (
        replace(
            _facts(gene, AtlasCuration.NOT_YET_CURATED),
            headline_confidence=Classification.LIMITED,
            confidence_by_lesion_group={},
        ),
        GeneValidity(
            records=(limited,),
            state=ValidityState.EXPERT_CURATED,
            has_source_discordance=False,
        ),
    )


def _submitter_admitted(gene: str) -> tuple[GeneFacts, GeneValidity]:
    """A gene no panel graded, admitted by two GenCC submitters agreeing."""
    submissions = tuple(
        ValidityRecord(
            source=ValiditySource.GENCC,
            classification=Classification.LIMITED,
            classification_term="Limited",
            disease="MONDO:0007732",
            disease_label="atrial septal defect",
            moi="Autosomal dominant",
            report_url=None,
            submitter=submitter,
        )
        for submitter in ("Ambry Genetics", "G2P")
    )
    return (
        replace(
            _facts(gene, AtlasCuration.NOT_YET_CURATED),
            headline_confidence=None,
            validity_state=ValidityState.SUBMITTER_CURATED,
            confidence_by_lesion_group={},
        ),
        GeneValidity(
            records=submissions,
            state=ValidityState.SUBMITTER_CURATED,
            has_source_discordance=False,
        ),
    )


def test_the_rail_names_the_grade_the_panel_gave_and_the_disease_it_gave_it_for(
    tmp_path: Path,
) -> None:
    """One claim in two halves, and before 2026-08-06 they could disagree.

    The rail's first row read `definitive for` unconditionally while
    `_graded_diseases` returned only `DEFINITIVE` records, so a `limited` gene
    rendered `definitive for | —`: a label naming a grade the panel did not give,
    over an em dash where the disease belonged. 69 of the 92 genes published
    today are in that position.

    Both halves now come from `fact.headline_confidence`, so this asserts the
    label, the disease and the chip class **as one string in one element**. A
    page-wide check on each separately would pass with the lookup keyed back on
    `DEFINITIVE`, because `render.grade_legend` puts the word "definitive" on
    every page (CLAUDE.md section 4.19).
    """
    fact, validity = _limited(TBX5)
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        {TBX5: fact},
        emitter,
        symbols={TBX5: "TBX5"},
        validity={TBX5: validity},
        assertions={},
        publications={},
        burden={},
        cohorts={},
    )

    page = _page(tmp_path, "HGNC_11604.html")
    rail = re.search(r'<aside class="rail">.*?</aside>', page, re.S)
    assert rail is not None
    rail_html = rail.group()

    assert 'class="chip chip-limited">limited</span>' in rail_html
    assert "<dt>limited for</dt><dd>atrial septal defect</dd>" in rail_html
    assert "definitive for" not in rail_html
    assert _EM_DASH not in rail_html.split("<dl>")[1].split("</dd>")[0]


def test_a_gene_no_panel_graded_names_the_submitters_that_admitted_it(
    tmp_path: Path,
) -> None:
    """The other warrant, on the page. 16 of the 92 genes published are here.

    `headline_confidence` is `null`, so there is no grade to head the row with
    and no disease a panel named. The row becomes `admitted by` and carries the
    submitters, which is the only place on the page a reader can see *why* the
    gene is here at all -- and the chip must read `not classified` in the
    outline style rather than borrowing any rung of the grade ramp, because "no
    panel looked" is not "graded weakest".
    """
    fact, validity = _submitter_admitted(GATA4)
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        {GATA4: fact},
        emitter,
        symbols={GATA4: "GATA4"},
        validity={GATA4: validity},
        assertions={},
        publications={},
        burden={},
        cohorts={},
    )

    page = _page(tmp_path, "HGNC_4173.html")
    rail = re.search(r'<aside class="rail">.*?</aside>', page, re.S)
    assert rail is not None
    rail_html = rail.group()

    assert 'class="chip chip-ungraded">not classified</span>' in rail_html
    assert "<dt>admitted by</dt><dd>Ambry Genetics; G2P</dd>" in rail_html
    for rung in ("definitive", "strong", "moderate", "limited"):
        assert f"chip-{rung}" not in rail_html


def test_both_page_kinds_gloss_every_grade_they_can_render(tmp_path: Path) -> None:
    """A colour ramp with no key is a ranking the reader has to guess at.

    Until 2026-08-06 the site needed no key: every published gene was
    `definitive` and the only chip was green. It now renders four rungs and an
    ungraded state, and `limited` -- a panel saying the case is not yet made --
    is the plurality at 43 of 92.

    Asserted on **both** page kinds, because a reader arriving from a search
    engine lands on the gene page and never sees the browse page. That is
    CLAUDE.md section 4.31 exactly: the browse page carried the evidence legend
    and the gene page did not, and the gene page is where the omission did the
    damage.

    Each rung's gloss is checked by a distinctive phrase rather than in full, so
    this fails when a rung stops being explained and not when a comma moves.

    **The label, the chip class and the gloss are asserted as one triple, and
    that is the whole point of this test.** It checked only that seven phrases
    appeared *somewhere on the page* until 2026-08-06, which says nothing about
    which rung each is attached to. Measured that day: painting `moderate` with
    `strong`'s chip -- inverting the colour ramp mid-ladder for the 10 genes on
    those two rungs -- **survived all 796 tests**, as did swapping two glosses.
    A ramp whose colours are checked only for presence is a ramp with no guard
    at all, and the chips are the part a reader reads first.

    Scoped to the `<dl>` inside the grade key rather than to the page, because
    every one of these words also occurs in the rail, the browse table or the
    scope rule; a page-wide substring check passes against a legend rendering
    nothing (CLAUDE.md section 4.19).
    """
    fact, validity = _limited(TBX5)
    emitter = Emitter(root=tmp_path)
    build_gene_index_page(
        {TBX5: fact}, emitter, symbols={TBX5: "TBX5"}, validity={TBX5: validity}, burden_counts={}
    )
    build_gene_pages(
        {TBX5: fact},
        emitter,
        symbols={TBX5: "TBX5"},
        validity={TBX5: validity},
        assertions={},
        publications={},
        burden={},
        cohorts={},
    )

    # (label a reader sees, chip class that colours it, phrase unique to its gloss)
    rungs = (
        ("definitive", "chip-definitive", "upheld over time"),
        ("strong", "chip-strong", "without the years of replication"),
        ("moderate", "chip-moderate", "several unrelated probands"),
        ("limited", "chip-limited", "not yet made"),
        ("no known association", "chip-warn", "found no evidence of a relationship"),
        (
            "not classified",
            "chip-ungraded",
            "two or more Gene Curation Coalition submitters independently assert it",
        ),
    )
    for name in ("index.html", "HGNC_11604.html"):
        page = _page(tmp_path, name)
        assert "What these classifications mean" in page, f"{name} carries no grade key"
        block = re.search(r'<details class="grade-key">.*?</details>', page, re.S)
        assert block, f"{name} carries no grade key"
        key = block.group(0)
        legend = re.search(r"<dl>(.*?)</dl>", key, re.S)
        assert legend, f"{name} carries no grade key rows"
        pairs = re.findall(r"<dt>(.*?)</dt><dd>(.*?)</dd>", legend.group(1), re.S)
        assert len(pairs) == len(rungs), f"{name} glosses {len(pairs)} rungs, expected {len(rungs)}"
        for (label, chip, gloss), (dt, dd) in zip(rungs, pairs, strict=True):
            assert chip in dt, f"{name}: {label} is not chipped {chip} -- got {dt}"
            assert f">{label}<" in dt, f"{name}: expected the label {label}, got {dt}"
            assert gloss in dd, f"{name}: {label} is not glossed by {gloss!r}"

        # The scope-qualified claim, on both page kinds. Unqualified it was false
        # for 5 published genes carrying an out-of-scope ClinGen `Disputed`.
        assert "disputed or refuted by an expert panel for a disease" in key
        assert "not published here at all" not in page


# --- Developmental expression -----------------------------------------------
#
# Every fixture below builds `ExpressionProfile` TypedDicts directly, the same
# choice `_burden_row` makes for `BurdenRow`: this is a rendering test, and
# going through `build/profiles.py`'s real pipeline (a TSV on disk, a quantile
# grid, `gene_expression_profiles` itself) would test that module a second
# time rather than the English this one chooses for its output.


def _tissue_entry(
    tissue: str = "Heart",
    median: float = 100.0,
    unit: str = "rpkm",
    n_samples: int = 3,
    placement: Placement | None = None,
    not_placed_reason: str | None = None,
) -> TissueProfileEntry:
    return TissueProfileEntry(
        tissue=tissue,
        median_abundance=median,
        unit=unit,
        n_samples=n_samples,
        placement=placement,
        not_placed_reason=not_placed_reason,
    )


def _placement(
    median_percentile: int = 97,
    q25_percentile: int | None = 92,
    q75_percentile: int | None = 99,
    n_genes: int = 19842,
) -> Placement:
    return Placement(
        q25_percentile=q25_percentile,
        median_percentile=median_percentile,
        q75_percentile=q75_percentile,
        median_abundance=100.0,
        unit="rpkm",
        n_samples=3,
        n_genes=n_genes,
        method="lowest percentile of a tied breakpoint (bisect_left)",
    )


def _specificity(
    tau: float = 0.483,
    tissues: tuple[str, ...] = ("Brain", "Heart", "Liver"),
    highest_in: str | None = "Heart",
) -> Specificity:
    return Specificity(
        tau=tau,
        scale="log2(x+1)",
        method="tau (Yanai et al. 2005)",
        tissues=tissues,
        n_tissues=len(tissues),
        highest_in=highest_in,
        medians=dict.fromkeys(tissues, 10.0),
    )


def _phase_info(
    outcome: str | None = "matched",
    phase_ids: tuple[str, ...] = ("septation",),
    reason: str | None = None,
) -> PhaseInfo:
    return PhaseInfo(outcome=outcome, phase_ids=phase_ids, reason=reason)


def _stage_entry(
    stage: str | None = "7wpc",
    phase: PhaseInfo | None = None,
    specificity: Specificity | None = None,
    specificity_unavailable_reason: str | None = None,
    tissues: tuple[TissueProfileEntry, ...] = (),
) -> StageProfileEntry:
    return StageProfileEntry(
        stage=stage,
        phase=phase if phase is not None else _phase_info(),
        specificity=specificity,
        specificity_unavailable_reason=specificity_unavailable_reason,
        tissues=list(tissues),
    )


def _dataset_profile_entry(
    dataset: str = "E-MTAB-6814",
    quantile_shard: str | None = "omics/profile_quantiles/E-MTAB-6814.json",
    stages: tuple[StageProfileEntry, ...] = (),
) -> DatasetProfileEntry:
    return DatasetProfileEntry(dataset=dataset, quantile_shard=quantile_shard, stages=list(stages))


def _expression_profile(datasets: tuple[DatasetProfileEntry, ...] = ()) -> ExpressionProfile:
    return ExpressionProfile(datasets=list(datasets))


def _profile_dataset(
    accession: str = "E-MTAB-6814",
    *,
    cardiac_tissues: tuple[str, ...] = ("Heart",),
    detection_floor: float = 1.0,
) -> Dataset:
    """A minimal, validly-constructed `Dataset`, mirroring `test_build_profiles.

    py`'s own `_dataset()` factory field for field, so this fixture is proven
    to satisfy `Dataset`'s `a_profile_dataset_is_fully_declared` validator
    rather than guessed at independently.
    """
    return Dataset(
        id=accession,
        archive="arrayexpress",
        technology="bulk_rnaseq",
        design="profile",
        tissue="whole embryo",
        developmental_stage="embryonic",
        organism="NCBITaxon:9606",
        n_samples=3,
        licence="CC BY 4.0",
        contrasts=[],
        cardiac_tissues=cardiac_tissues,
        detection_floor=detection_floor,
        floor_source="source methods, section 4",
        quantile_estimator="linear",
        stages=(Stage(token="7wpc", wpc=7.0, order=1),),
    )


def _cardiac_phase(
    phase_id: str = "atrial_septum_morphogenesis",
    *,
    label: str = "Atrial septum morphogenesis",
    end_wpc: float | None = 6.29,
) -> CardiacPhase:
    """One curated phase, `end_basis` derived from whether an end is given.

    Derived rather than passed, so a fixture cannot construct the state
    `CardiacPhase.end_fields_match_end_basis` forbids and then be read as
    evidence about a state the real vocabulary can hold. `end_wpc=None` is the
    `NOT_STATED` case -- `heart_looping` in `curation/cardiac_phases.yaml`,
    the one phase whose end the source never states.
    """
    stated = end_wpc is not None
    return CardiacPhase(
        id=phase_id,
        go_id="GO:0060413",
        label=label,
        start_wpc=3.71,
        end_wpc=end_wpc,
        start_carnegie_stage="CS12",
        end_carnegie_stage="CS17" if stated else None,
        start_hsapdv_id="HsapDv:0000019",
        end_hsapdv_id="HsapDv:0000024" if stated else None,
        end_basis=EndBasis.STATED if stated else EndBasis.NOT_STATED,
    )


def _cardiac_phases(*phases: CardiacPhase) -> CardiacPhaseFile:
    """The curated phase vocabulary a page bands its trajectory with."""
    return CardiacPhaseFile(
        attributed_to="Buijtendijk MFJ, Barnett P, van den Hoff MJB (2020)",
        citation="PMID:32048790",
        phases=list(phases or (_cardiac_phase(),)),
    )


def _expression_page(
    tmp_path: Path,
    profiles: dict[str, ExpressionProfile],
    datasets: dict[str, Dataset] | None = None,
    facts: dict[str, GeneFacts] | None = None,
    name: str = "HGNC_4173.html",
    phases: CardiacPhaseFile | None = None,
) -> str:
    """One gene page, built only to exercise the expression section.

    Every other input is the smallest fixture `build_gene_pages` accepts: no
    burden, no curated assertions, one mirrored validity record. `facts`
    defaults to a single uncurated GATA4 so a caller testing one gene need not
    build a `GeneFacts` by hand.

    `phases` defaults to `None` -- the corpus with no curated phase
    vocabulary at all -- rather than to `_cardiac_phases()`, so a test that
    wants a banded trajectory has to ask for one. Every expression test
    written before the trajectory landed therefore still exercises exactly
    what it exercised then.
    """
    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        facts or {GATA4: _facts(GATA4, AtlasCuration.NOT_YET_CURATED)},
        emitter,
        symbols={GATA4: "GATA4", TBX5: "TBX5"},
        validity={GATA4: _validity(), TBX5: _validity()},
        assertions={},
        publications={},
        burden={},
        cohorts={},
        profiles=profiles,
        datasets=datasets or {},
        phases=phases,
    )
    return _page(tmp_path, name)


def _slice_between(text: str, start_marker: str, end_markers: tuple[str, ...]) -> str:
    """`text` after `start_marker`, up to the first of `end_markers` (or the end)."""
    start = text.index(start_marker) + len(start_marker)
    rest = text[start:]
    positions = [rest.index(marker) for marker in end_markers if marker in rest]
    return rest[: min(positions)] if positions else rest


def _expression_section_text(page: str) -> str:
    """The `Developmental expression` section alone, up to the next `<h2>`.

    Sliced for the reason `_validity_table` is sliced: a page-wide assertion
    cannot tell this section's own claim from an unrelated occurrence of the
    same words elsewhere on the page (CLAUDE.md section 4.19) -- the research
    notice, the scope rule and the burden section all sit on the same page.
    """
    return _slice_between(page, "<h2>Developmental expression</h2>", ("<h2>",))


def test_a_gene_with_no_expression_data_says_the_atlas_has_not_curated_it(
    tmp_path: Path,
) -> None:
    """`{"datasets": []}` -- `EMPTY_EXPRESSION_PROFILE` itself, not a look-alike
    this test built by hand, so the page and `bundles.py`'s fallback are proven
    to agree on what "empty" means.

    **Chosen wording, and the rejected one.** "The atlas has not yet curated a
    profile for this gene" is the true claim: no dataset mentions this gene at
    all, so there is no specific dataset in view to say "does not cover it"
    about, and that alternative would invent a check this atlas never ran. The
    section still renders -- never an absent one, the `_not_curated` reason:
    an absent section cannot be told apart from "the atlas looked and found
    nothing".
    """
    page = _expression_page(tmp_path, {GATA4: EMPTY_EXPRESSION_PROFILE})
    section = _expression_section_text(page)

    assert "not yet curated" in section
    assert "does not cover" not in section
    assert "<h3>" not in section, "an empty profile must name no dataset"
    assert "<table" not in section


def test_the_gloss_never_calls_a_gene_heart_preferential_on_another_organs_peak(
    tmp_path: Path,
) -> None:
    """Fixture: tau peaking in liver. A page-wide assertion cannot catch this.

    Same shape as `_validity_table()` -- a page-wide check once passed with the
    renderer mutated to print bare ids, because the same strings rendered in a
    <details> block further down. Slice to the section before asserting.

    Both genes share one dataset whose only declared `cardiac_tissues` entry is
    Heart, so the two pages differ only in `highest_in` -- never in whether
    Heart itself is cardiac -- which is what makes this a test of the gate and
    not of the fixture.
    """
    dataset = _profile_dataset(cardiac_tissues=("Heart",))
    heart_peak = _expression_profile(
        (
            _dataset_profile_entry(
                stages=(_stage_entry(specificity=_specificity(highest_in="Heart")),)
            ),
        )
    )
    liver_peak = _expression_profile(
        (
            _dataset_profile_entry(
                stages=(_stage_entry(specificity=_specificity(highest_in="Liver")),)
            ),
        )
    )

    emitter = Emitter(root=tmp_path)
    build_gene_pages(
        {
            TBX5: _facts(TBX5, AtlasCuration.CURATED),
            GATA4: _facts(GATA4, AtlasCuration.NOT_YET_CURATED),
        },
        emitter,
        symbols={TBX5: "TBX5", GATA4: "GATA4"},
        validity={TBX5: _validity(), GATA4: _validity()},
        assertions={},
        publications={},
        burden={},
        cohorts={},
        profiles={TBX5: heart_peak, GATA4: liver_peak},
        datasets={"E-MTAB-6814": dataset},
    )

    heart_section = _expression_section_text(_page(tmp_path, "HGNC_11604.html"))
    liver_section = _expression_section_text(_page(tmp_path, "HGNC_4173.html"))

    assert "heart-preferential" in heart_section
    assert "peaks in Heart" in heart_section

    assert "heart-preferential" not in liver_section
    assert "peaks in Liver" in liver_section
    assert "does not treat as a cardiac tissue" in liver_section


def test_the_gloss_says_neither_organ_when_the_peak_is_tied(tmp_path: Path) -> None:
    """`highest_in: null` -- organs tied at the peak. Neither wording applies.

    Distinct from the cardiac/non-cardiac test above: this fixture's dataset
    still declares Heart as cardiac, so a gate that fired on "any tissues at
    all" rather than specifically on `highest_in is None` would pass the other
    test and fail only this one.
    """
    profile = _expression_profile(
        (_dataset_profile_entry(stages=(_stage_entry(specificity=_specificity(highest_in=None)),)),)
    )
    page = _expression_page(
        tmp_path, {GATA4: profile}, {"E-MTAB-6814": _profile_dataset(cardiac_tissues=("Heart",))}
    )
    section = _expression_section_text(page)

    assert "heart-preferential" not in section
    assert "peaks in" not in section
    assert "tied" in section
    # Rule 3 still holds in the tied case: the number, its scale and its organ
    # list are not conditioned on the peak being resolvable.
    assert "τ = 0.483 (log2(x+1) scale) across Brain, Heart, Liver (3 organs sampled)." in section


def test_tau_never_renders_without_its_scale_and_its_organ_list(tmp_path: Path) -> None:
    """No branch omits the qualifier -- the `_effect`/`count_unit` precedent.

    Pinned as one literal composed sentence rather than as three separate
    membership checks, so dropping the scale, dropping a tissue, or reordering
    the tissue list all fail this one assertion instead of needing three.
    `test_both_page_kinds_state_the_rule_that_admits_a_gene_to_this_atlas`'s
    own lesson: pin published wording against a literal, never against the
    constant that produced it.
    """
    profile = _expression_profile(
        (
            _dataset_profile_entry(
                stages=(
                    _stage_entry(
                        specificity=_specificity(
                            tau=0.483, tissues=("Brain", "Heart", "Liver"), highest_in="Heart"
                        )
                    ),
                )
            ),
        )
    )
    section = _expression_section_text(_expression_page(tmp_path, {GATA4: profile}))

    assert "τ = 0.483 (log2(x+1) scale) across Brain, Heart, Liver (3 organs sampled)." in section


def test_a_below_floor_tissue_names_the_organ_and_the_floor_value_and_never_says_not_detected(
    tmp_path: Path,
) -> None:
    """A negative assertion on published bytes, sharing a fixture with the
    positive wording it replaced.

    Grep the BUILT PAGE, never the source: `pages.py` builds copy from adjacent
    string literals, so a banned phrase split across two source lines returns
    nothing from a source grep while rendering verbatim on the page. That has
    happened here, twice, and the negative assertion is what caught it.
    """
    profile = _expression_profile(
        (
            _dataset_profile_entry(
                stages=(
                    _stage_entry(
                        tissues=(
                            _tissue_entry(
                                tissue="Liver",
                                median=0.3,
                                unit="rpkm",
                                placement=None,
                                not_placed_reason="below_detection_floor",
                            ),
                        )
                    ),
                )
            ),
        )
    )
    page = _expression_page(
        tmp_path, {GATA4: profile}, {"E-MTAB-6814": _profile_dataset(detection_floor=1.0)}
    )
    section = _expression_section_text(page)

    assert "below the detection floor in whole Liver at this stage" in section
    assert "detection floor 1 rpkm" in section
    assert "not detected" not in page.lower()


def test_n_samples_travels_with_every_figure_and_distinguishes_a_single_observation(
    tmp_path: Path,
) -> None:
    """The schema permits `n_samples` = 1, so "median" can be a single
    observation, and the page must say so -- and say something different for
    three.
    """
    profile = _expression_profile(
        (
            _dataset_profile_entry(
                stages=(
                    _stage_entry(
                        tissues=(
                            _tissue_entry(tissue="Heart", n_samples=3, placement=_placement()),
                            _tissue_entry(
                                tissue="Liver",
                                n_samples=1,
                                median=0.3,
                                placement=None,
                                not_placed_reason="below_detection_floor",
                            ),
                        )
                    ),
                )
            ),
        )
    )
    section = _expression_section_text(
        _expression_page(
            tmp_path, {GATA4: profile}, {"E-MTAB-6814": _profile_dataset(detection_floor=1.0)}
        )
    )

    assert "n=3 samples" in section
    assert "n=1 sample)" in section
    assert "n=1 samples" not in section


def test_tau_undefined_at_one_organ_says_so_rather_than_a_number(tmp_path: Path) -> None:
    """`specificity: null` (one organ sampled) beside a stage with a real tau,
    on the same gene page -- the fixture pair CLAUDE.md section 4.14 asks for.
    """
    profile = _expression_profile(
        (
            _dataset_profile_entry(
                stages=(
                    _stage_entry(
                        stage="one_organ",
                        specificity=None,
                        specificity_unavailable_reason="one_organ_sampled",
                    ),
                    _stage_entry(
                        stage="two_organs",
                        specificity=_specificity(tissues=("Heart", "Liver"), highest_in="Heart"),
                    ),
                )
            ),
        )
    )
    section = _expression_section_text(_expression_page(tmp_path, {GATA4: profile}))

    one_organ = _slice_between(section, "<h4>one_organ</h4>", ("<h4>",))
    two_organs = _slice_between(section, "<h4>two_organs</h4>", ("<h4>",))

    assert "τ is not available: only one organ was sampled at this stage" in one_organ
    assert "τ =" not in one_organ
    assert "τ = 0.483" in two_organs


def test_every_stage_outcome_reads_as_a_state_never_a_bug_or_blank(tmp_path: Path) -> None:
    """Pins all four `PhaseOutcome` members, so a gene page is proven correct
    for a stage in each state -- including MATCHED with only one phase, which
    `test_a_matched_stage_with_several_phases_names_every_one` deliberately
    does NOT also cover, per the fixture-diversity rule: a fixture where every
    MATCHED stage carries the same number of phases could not distinguish
    "renders the one phase" from "renders `phase_ids[0]` and drops the rest".
    """
    profile = _expression_profile(
        (
            _dataset_profile_entry(
                stages=(
                    _stage_entry(
                        stage="s_matched",
                        phase=_phase_info(
                            outcome="matched", phase_ids=("outflow_tract_septation",)
                        ),
                    ),
                    _stage_entry(
                        stage="s_outside",
                        phase=_phase_info(
                            outcome="outside_window",
                            phase_ids=(),
                            reason="outside the curated window",
                        ),
                    ),
                    _stage_entry(
                        stage="s_postnatal",
                        phase=_phase_info(outcome="post_natal", phase_ids=(), reason="post-natal"),
                    ),
                    _stage_entry(
                        stage="s_undeclared",
                        phase=_phase_info(
                            outcome="undeclared",
                            phase_ids=(),
                            reason="stage not declared by this dataset",
                        ),
                    ),
                )
            ),
        )
    )
    section = _expression_section_text(_expression_page(tmp_path, {GATA4: profile}))

    matched = _slice_between(section, "<h4>s_matched</h4>", ("<h4>",))
    outside = _slice_between(section, "<h4>s_outside</h4>", ("<h4>",))
    postnatal = _slice_between(section, "<h4>s_postnatal</h4>", ("<h4>",))
    undeclared = _slice_between(section, "<h4>s_undeclared</h4>", ("<h4>",))

    assert "Developmental phase: <strong>outflow tract septation</strong>." in matched
    assert "Developmental phase: outside the curated window." in outside
    assert "Developmental phase: post-natal." in postnatal
    assert "Developmental phase: stage not declared by this dataset." in undeclared
    for block in (matched, outside, postnatal, undeclared):
        assert "Developmental phase: .</p>" not in block, "a blank reason is a rendering bug"


def test_a_matched_stage_with_several_phases_names_every_one() -> None:
    """A stage inside an overlap must name ALL of its phases, not just the
    first -- the unique killer of a renderer that reads only `phase_ids[0]`.

    Deliberately three phases, not two: with two, a truncating renderer that
    always shows exactly one and a correct renderer that drops the last of
    two could both pass a weaker assertion. Three phases makes "shows one",
    "shows two of three" and "shows all three" three distinct, checkable
    outcomes, and the singular/plural noun is checked on both sides so a
    mutant hard-coding "phases" (plural) cannot pass the single-phase test
    above either.
    """
    phase = _phase_info(
        outcome="matched",
        phase_ids=("atrial_septation", "outflow_tract_septation", "ventricular_septation"),
    )
    sentence = _phase_sentence(phase)
    assert sentence == (
        "Developmental phases: <strong>atrial septation</strong>, "
        "<strong>outflow tract septation</strong>, "
        "<strong>ventricular septation</strong>."
    )


def test_the_bulk_dilution_caveat_is_unconditional_not_only_beside_bad_news(
    tmp_path: Path,
) -> None:
    """`_POOLING_NOTICE` is the recorded precedent: made conditional once, it
    read false in the position that motivated the fix. A caveat present only
    beside a below-floor tissue would read as an excuse for that gene's own
    numbers rather than as the general fact it is.
    """
    all_detected = _expression_profile(
        (
            _dataset_profile_entry(
                stages=(_stage_entry(tissues=(_tissue_entry(placement=_placement()),)),)
            ),
        )
    )
    caveat = "Every measurement below is from whole, bulk tissue, never a single cell type."

    good_news_section = _expression_section_text(_expression_page(tmp_path, {GATA4: all_detected}))
    assert caveat in good_news_section

    below_floor = _expression_profile(
        (
            _dataset_profile_entry(
                stages=(
                    _stage_entry(
                        tissues=(
                            _tissue_entry(
                                median=0.3,
                                placement=None,
                                not_placed_reason="below_detection_floor",
                            ),
                        )
                    ),
                )
            ),
        )
    )
    bad_news_section = _expression_section_text(
        _expression_page(
            tmp_path,
            {GATA4: below_floor},
            {"E-MTAB-6814": _profile_dataset(detection_floor=1.0)},
        )
    )
    assert caveat in bad_news_section


def test_percentiles_are_captioned_as_not_comparable_only_when_more_than_one_is_shown(
    tmp_path: Path,
) -> None:
    """Rule 8: not comparable across organs or stages -- and only said where
    the page actually lays out more than one for a reader to compare.

    Two placements, from two different organs, so the count this gate reads is
    proven to sum across tissues rather than merely check "is there a
    placement at all".
    """
    caveat = "not comparable across organs or across developmental stages"

    one_placement = _expression_profile(
        (
            _dataset_profile_entry(
                stages=(_stage_entry(tissues=(_tissue_entry(placement=_placement()),)),)
            ),
        )
    )
    lone_section = _expression_section_text(_expression_page(tmp_path, {GATA4: one_placement}))
    assert caveat not in lone_section
    # The denominator travels with the percentile even when the caveat does not.
    assert "97 of 19,842 genes" in lone_section

    two_placements = _expression_profile(
        (
            _dataset_profile_entry(
                stages=(
                    _stage_entry(
                        tissues=(
                            _tissue_entry(tissue="Heart", placement=_placement()),
                            _tissue_entry(
                                tissue="Liver", placement=_placement(median_percentile=40)
                            ),
                        )
                    ),
                )
            ),
        )
    )
    busy_section = _expression_section_text(_expression_page(tmp_path, {GATA4: two_placements}))
    assert caveat in busy_section


def test_a_tissue_stage_and_dataset_token_carrying_markup_is_escaped(tmp_path: Path) -> None:
    """Tissue and stage tokens are mirrored third-party strings, and this is the
    one artifact kind where a `<` in an organ label is the difference between a
    string and a script tag.
    """
    profile = _expression_profile(
        (
            _dataset_profile_entry(
                dataset="<i>E-MTAB-x</i>",
                quantile_shard=None,
                stages=(
                    _stage_entry(
                        stage="<b>7wpc</b>",
                        tissues=(
                            _tissue_entry(tissue="<script>xss</script>", placement=_placement()),
                        ),
                    ),
                ),
            ),
        )
    )
    page = _expression_page(tmp_path, {GATA4: profile})

    assert "<script>xss</script>" not in page
    assert "&lt;script&gt;xss&lt;/script&gt;" in page
    assert "<b>7wpc</b>" not in page
    assert "&lt;b&gt;7wpc&lt;/b&gt;" in page
    assert "<i>E-MTAB-x</i>" not in page
    assert "&lt;i&gt;E-MTAB-x&lt;/i&gt;" in page


def test_the_dataset_heading_links_to_its_own_percentile_grid(tmp_path: Path) -> None:
    """D39(b)'s auditability reaching a reader, not only a program: the same
    relative-path convention `_rail`'s "this gene as JSON" link uses.
    """
    profile = _expression_profile((_dataset_profile_entry(),))
    page = _expression_page(tmp_path, {GATA4: profile})

    assert '<a href="../omics/profile_quantiles/E-MTAB-6814.json">' in page


# --- The phase-banded trajectory (D43) --------------------------------------
#
# Every fixture below differs from every other in the one value the code under
# test reads: how many heart stages carry a percentile (three, one, none) and
# whether the phase a stage names has a curated end. CLAUDE.md section 4.36 is
# the reason that is spelled out here -- a page fixture whose genes all share
# the value under test measures nothing, and it has cost this repository four
# separate defects.

_HEART_DATASET = {"E-MTAB-6814": _profile_dataset(cardiac_tissues=("heart",))}


def _heart_series(
    *stages: tuple[str, float | None],
    phase_ids: tuple[str, ...] = ("atrial_septum_morphogenesis",),
    unplaced_reason: str = "below_detection_floor",
) -> ExpressionProfile:
    """A heart series over `stages`, in the order given.

    A `float` is a median this dataset placed against its percentile grid; a
    `None` is a stage sampled and not placed, for `unplaced_reason`. The order
    given is the published order -- `gene_expression_profiles` emits stages
    by `Stage.order`, so a fixture's own sequence is the trajectory's.

    `unplaced_reason` is a parameter rather than a constant because the reason
    changes what a page is entitled to say: below a floor is a low reading
    about the gene, and a missing percentile grid is a hole in the reference
    that says nothing about the gene at all.
    """
    entries = [
        _stage_entry(
            stage=token,
            phase=_phase_info(phase_ids=phase_ids),
            tissues=(
                _tissue_entry(
                    tissue="heart",
                    median=median if median is not None else 0.3,
                    unit="tpm",
                    placement=_placement() if median is not None else None,
                    not_placed_reason=None if median is not None else unplaced_reason,
                ),
            ),
        )
        for token, median in stages
    ]
    return _expression_profile((_dataset_profile_entry(stages=tuple(entries)),))


def _dataset_lede(section: str) -> str:
    """One dataset block's opening: its chart, or the sentence instead of one.

    Sliced off the stage tables the way `_validity_table` is sliced, and for
    the same reason (CLAUDE.md section 4.19). `_BULK_DILUTION_NOTICE` sits
    above every dataset block and itself contains "not evidence that", so a
    section-wide assertion on that phrase passes whether or not the sentence
    under test was rendered at all.
    """
    return _slice_between(section, "</h3>", ("<details",))


def test_every_chart_carries_an_atlas_specific_axis(tmp_path: Path) -> None:
    """D43: a chart earns its place by an axis the source's browser lacks.

    Replaces `test_the_section_never_renders_a_chart_only_text_and_tables`,
    which asserted the proxy "never an `<svg>`". The proxy forbade this
    trajectory -- which no external browser has -- while permitting the
    rendering that hid the alphabetical-ordering defect for three releases.

    A bare re-plot of the source's own curve, with no phase band, must fail
    here.
    """
    profile = _heart_series(("4wpc", 5.0), ("5wpc", 40.0), ("6wpc", 275.0))
    section = _expression_section_text(
        _expression_page(tmp_path, {GATA4: profile}, _HEART_DATASET, phases=_cardiac_phases())
    )

    assert "<svg" in section
    assert "chart-band" in section, "a trajectory with no phase band is the re-plot D43 forbids"

    # The same three measurements with no curated phase vocabulary behind them
    # are exactly the source's own curve, and are refused.
    bare = _expression_section_text(
        _expression_page(tmp_path, {GATA4: profile}, _HEART_DATASET, phases=None)
    )
    assert "<svg" not in bare


def test_a_phase_the_source_never_ended_bands_nothing(tmp_path: Path) -> None:
    """`end_basis: not_stated` has no curated width, so it has no band.

    `heart_looping` is the live case: the review states when it starts and
    never when it stops, so `CardiacPhaseFile.phases_for` excludes it at every
    wpc. A band drawn to an invented edge would assert exactly the boundary
    that field exists to refuse -- and would do it in the one encoding a
    reader takes at a glance.

    The guard is here rather than trusted to `phases_for` upstream: a stage
    carrying an unended phase's id is a state the real pipeline never
    produces, which is why this fixture builds one by hand. CLAUDE.md section
    4.28 -- a guard added to one layer is not a guard.
    """
    vocabulary = _cardiac_phases(
        _cardiac_phase(),
        _cardiac_phase("heart_looping", label="Heart looping", end_wpc=None),
    )
    both = _heart_series(
        ("4wpc", 5.0),
        ("5wpc", 40.0),
        ("6wpc", 275.0),
        phase_ids=("atrial_septum_morphogenesis", "heart_looping"),
    )
    section = _expression_section_text(
        _expression_page(tmp_path, {GATA4: both}, _HEART_DATASET, phases=vocabulary)
    )
    assert section.count('class="chart-band"') == 1
    assert "Atrial septum morphogenesis" in section
    assert "Heart looping" not in section

    # And with nothing else to band, there is no admissible chart at all.
    alone = _heart_series(
        ("4wpc", 5.0), ("5wpc", 40.0), ("6wpc", 275.0), phase_ids=("heart_looping",)
    )
    unbanded = _expression_section_text(
        _expression_page(tmp_path, {GATA4: alone}, _HEART_DATASET, phases=vocabulary)
    )
    assert "<svg" not in unbanded


def test_a_stage_below_the_floor_is_an_axis_tick_and_never_a_zero(tmp_path: Path) -> None:
    """Below the floor, not sampled, and zero are three different claims.

    Plotting a below-floor stage at zero would put it on the axis as the
    lowest *measurement*; leaving a gap would make it indistinguishable from a
    stage nobody sampled.

    The counts are asserted, not only the class name: a renderer that ticked
    every stage, or that ticked the below-floor stage *and* plotted it, would
    satisfy a bare `in` check.
    """
    profile = _heart_series(("4wpc", 5.0), ("5wpc", None), ("6wpc", 40.0), ("7wpc", 275.0))
    section = _expression_section_text(
        _expression_page(tmp_path, {GATA4: profile}, _HEART_DATASET, phases=_cardiac_phases())
    )

    assert "chart-absent" in section
    assert section.count('class="chart-absent"') == 1
    assert section.count('class="chart-point"') == 3


def test_the_line_never_spans_a_stage_the_atlas_did_not_place(tmp_path: Path) -> None:
    """A below-floor stage breaks the line; it is not drawn straight over.

    Measured 2026-08-20 on the committed corpus: 20 of the 85 charted genes
    carry a below-floor stage strictly between two plotted points, and TBX1
    has 11 of them against 5 placed stages. A single line through the placed
    points runs smooth and high across exactly the stages where this dataset
    measured below its own floor -- the picture asserting a continuity the
    figures deny.

    The fixture has two runs and a lone leading point, so "one line through
    everything", "splits but drops a run" and "splits correctly" are three
    distinguishable outcomes rather than two.
    """
    profile = _heart_series(
        ("4wpc", 5.0),
        ("5wpc", None),
        ("6wpc", 40.0),
        ("7wpc", 275.0),
        ("8wpc", None),
        ("9wpc", 60.0),
        ("10wpc", 90.0),
    )
    section = _expression_section_text(
        _expression_page(tmp_path, {GATA4: profile}, _HEART_DATASET, phases=_cardiac_phases())
    )

    drawn = re.findall(r'<polyline class="chart-line" points="([^"]+)"', section)
    assert len(drawn) == 2, "one line per run of adjacent placed stages"
    assert [len(run.split(" ")) for run in drawn] == [2, 2]
    # The stage that broke the line is still a point on the page, and the two
    # stages that broke it are still ticked.
    assert section.count('class="chart-point"') == 5
    assert section.count('class="chart-absent"') == 2


def test_a_gene_placed_at_one_stage_gets_markers_and_no_line(tmp_path: Path) -> None:
    """Tier 3. TFAP2B is placed at 1 of 19 stages in the committed corpus.

    A polyline through one point draws a trend the data does not contain.
    """
    profile = _heart_series(("4wpc", 5.0), ("5wpc", None), ("6wpc", None))
    section = _expression_section_text(
        _expression_page(tmp_path, {GATA4: profile}, _HEART_DATASET, phases=_cardiac_phases())
    )

    assert "<polyline" not in section
    assert "chart-point" in section

    # Two adjacent placed stages is the case the constant actually gates, and
    # the one above cannot reach it: a run of a single point is dropped
    # whatever the threshold says, so `_STAGES_FOR_A_TRAJECTORY = 1` survives
    # a one-stage fixture. Two measurements draw a segment that claims
    # something about the interval between them.
    pair = _heart_series(("4wpc", 5.0), ("5wpc", 40.0), ("6wpc", None))
    boundary = _expression_section_text(
        _expression_page(tmp_path, {GATA4: pair}, _HEART_DATASET, phases=_cardiac_phases())
    )

    assert "<polyline" not in boundary
    assert boundary.count('class="chart-point"') == 2


def test_a_gene_never_above_the_floor_gets_a_sentence_and_no_chart(tmp_path: Path) -> None:
    """D42, tier 4, and the reason it is not cosmetic.

    Measured 2026-08-20: seven published genes are below the floor in heart at
    every one of 19 stages -- SEMA3E, ZIC3, USP44, DAW1, FGF8, GDF1, NODAL.
    ZIC3 and NODAL are ClinGen definitive. An empty chart beside a definitive
    chip asserts the thing the bulk-dilution caveat denies.

    Asserted against literals, not against the constant that produced them: a
    test that imports the string it asserts on compares the module to itself
    and passes for any rewording.

    Sliced to the dataset block's own lede before asserting, because
    `_BULK_DILUTION_NOTICE` carries "not evidence that" over every gene page
    on this site, charted or not.
    """
    profile = _heart_series(("4wpc", None), ("5wpc", None), ("6wpc", None))
    section = _dataset_lede(
        _expression_section_text(
            _expression_page(tmp_path, {GATA4: profile}, _HEART_DATASET, phases=_cardiac_phases())
        )
    )

    assert "<svg" not in section
    assert "detection floor at every" in section
    assert "not evidence that" in section
    assert "heart" in section
    assert "3 stages sampled" in section


def test_a_gap_in_the_reference_is_never_reported_as_a_low_reading(tmp_path: Path) -> None:
    """ "Nothing placed" is not one fact, and only one of them is about the gene.

    `below_detection_floor` says this dataset measured the gene under its own
    floor, which is what earns the dilution argument beside it. Every other
    gap -- no percentile grid published for this organ and stage, no floor
    declared, the dataset unregistered -- is a hole in the *reference*, and
    reporting it as a low reading asserts a measurement this dataset never
    made. Shares its shape with the below-floor test above and differs only in
    the recorded reason, which is the value the branch reads.
    """
    profile = _heart_series(
        ("4wpc", None),
        ("5wpc", None),
        ("6wpc", None),
        unplaced_reason="no_quantile_grid",
    )
    section = _dataset_lede(
        _expression_section_text(
            _expression_page(tmp_path, {GATA4: profile}, _HEART_DATASET, phases=_cardiac_phases())
        )
    )

    assert "<svg" not in section
    assert "detection floor at every" not in section
    assert "no complete percentile grid is published for this organ at this stage" in section


def test_a_chart_never_replaces_the_figures_it_summarises(tmp_path: Path) -> None:
    """The evidence-loss mutant.

    A chart that removed the table would take every exact figure out of the
    HTML -- this repository's characteristic defect, curated work reaching no
    page. The table moves into `<details>`; it does not leave.
    """
    profile = _heart_series(("4wpc", 5.0), ("5wpc", 40.0), ("6wpc", 275.0))
    section = _expression_section_text(
        _expression_page(tmp_path, {GATA4: profile}, _HEART_DATASET, phases=_cardiac_phases())
    )

    assert "<svg" in section
    folded = re.search(r'<details class="stage-figures">.*?</details>', section, re.S)
    assert folded is not None, "the stage figures must be folded, not deleted"
    # Not merely present, and not merely somewhere on the page: every figure
    # the 21 blocks carried is inside the fold the chart summarises.
    assert "<table" in folded.group(0)
    assert "275 tpm (n=3 samples)" in folded.group(0)
    assert "97 of 19,842 genes" in folded.group(0)
    # And the section-wide `<details` assertion this replaces is not enough on
    # its own: `_EXPRESSION_READING_NOTES` renders one a few lines above,
    # so a renderer that dropped the wrapper entirely would still satisfy it
    # (CLAUDE.md section 4.19).
    unfolded = re.sub(r"<details.*?</details>", "", section, flags=re.S)
    assert "<table" not in unfolded
