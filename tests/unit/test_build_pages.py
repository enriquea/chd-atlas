# tests/unit/test_build_pages.py
"""What a gene page must say, especially when the atlas has curated nothing."""

from __future__ import annotations

import re
from dataclasses import replace
from datetime import date
from enum import StrEnum
from pathlib import Path
from types import EllipsisType
from typing import Final

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


_SIERANT_PUBLICATION = {
    "PMID:40127276": Publication(
        id="PMID:40127276",
        title="Landscape of rare variants in congenital heart disease",
        journal="Nature genetics",
        year=2025,
        authors=["Sierant MC"],
        study_type="case_control",  # type: ignore[arg-type]
        own_lab=False,
        tests_reported=4128,
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


def _burden_section_text(page: str) -> str:
    """The `Rare variant burden` section alone, up to the next `<h2>`.

    Sliced for the reason `_validity_table` and `_expression_section_text` are
    sliced (CLAUDE.md section 4.19). The developmental-expression section that
    follows this one draws a chart of its own on nearly every page, so a
    page-wide `"<svg" not in page` can be satisfied -- or defeated -- by a
    picture the burden section did not draw.
    """
    return _slice_between(page, "<h2>Rare variant burden</h2>", ("<h2>",))


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


# --- The effect-and-uncertainty forest ------------------------------------
#
# Every fixture below carries **both** values of whatever it is testing.
# CLAUDE.md section 4.36 is why that is spelled out rather than assumed: a
# fixture whose rows all share the value under test measures nothing, and it
# has cost this repository four separate defects -- the last of them in a page
# fixture exactly like these.


def _forest_figures(section: str) -> list[str]:
    """Each `<figure class="forest">` on the page, whole, in document order."""
    return re.findall(r'<figure class="forest".*?</figure>', section, flags=re.S)


def test_two_effect_measures_never_share_one_axis(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """The merge `effect_measure` exists to prevent, and the order it is drawn in.

    779 corpus rows carry an odds ratio and 75 an enrichment ratio. An odds
    ratio of 3.1 and a de novo enrichment of 3.1 are different claims.

    **The panels' order is a published fact and was pinned by nothing.**
    `pages.py` builds the facet list as `sorted({row.effect_measure ...})`
    over a `set`, and removing that `sorted` survived the whole suite: the
    assertion below reads `sorted(set(figures))`, which sorts its own actual
    value and therefore cannot see order at all.
    `test_build_is_reproducible.py` cannot see it either -- it builds twice in
    one process, where `PYTHONHASHSEED` is fixed for the interpreter's life --
    so two builds of one commit on two machines would publish 25 gene pages
    with their forest panels in opposite orders while every checksum verified.

    **The real vocabulary is far too small to fixture safely, and this is
    the case CLAUDE.md section 4.42 says a synthetic set is for.** A probe
    constructing the exact `set` under 300 explicit `PYTHONHASHSEED` values
    measured the published pair `{enrichment_ratio, odds_ratio}` iterating
    pre-sorted **149/300 (49.7%)** -- a coin flip, so a literal-order
    assertion on the real pair is not a guard. All three members of
    `EffectMeasure` measured 73/300 (24.3%), and there is no fourth to widen
    into. Eight synthetic tokens measured **0/300**, matching
    `test_prf008_message_lists_every_offending_unit_in_sorted_order`'s fix for
    the identical constraint. `_effect` resolves an unrecognised measure to
    its own raw token by design, so the tokens render rather than being
    silently dropped, and `_burden_row` builds a `BurdenRow` directly -- which
    vocabulary a token belongs to is `tables.py`'s check, not this one's.
    """
    section = _burden_section_text(
        _burden_page(
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
    )

    figures = re.findall(r'data-effect-measure="([^"]+)"', section)
    assert sorted(set(figures)) == ["enrichment_ratio", "odds_ratio"]
    assert len(figures) == len(set(figures)), "one panel carried two measures"
    # And the two panels are two pictures, not one picture drawn twice: the
    # rows are otherwise identical, so a facet keyed on `study` alone would
    # emit a single figure holding both numbers on one axis.
    assert len(_forest_figures(section)) == 2

    # The order, on a set wide enough that a coincidentally-sorted iteration
    # is not what is being measured. Asserted as the published sequence, not
    # as `sorted(...)` of it.
    synthetic = ["golf", "alpha", "hotel", "charlie", "echo", "bravo", "foxtrot", "delta"]
    ordered = _burden_section_text(
        _burden_page(
            tmp_path,
            facts_uncurated,
            [
                _burden_row(
                    consequence_class=token,
                    effect=2.45,
                    effect_measure=token,
                    effect_bound=None,
                    ci_low=1.2,
                    ci_high=8.1,
                )
                for token in synthetic
            ],
        )
    )
    assert re.findall(r'data-effect-measure="([^"]+)"', ordered) == sorted(synthetic)


def test_a_study_publishing_no_interval_still_gets_a_panel(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """Otherwise the only picture on TBX5's page is of the null study.

    PMID:40127276 reports de novo LOF enriched 297x, q 5.6e-08, surviving its
    own correction, and publishes no intervals. Drawing only where intervals
    exist demotes the surviving finding to text.
    """
    section = _burden_section_text(
        _burden_page(
            tmp_path,
            facts_uncurated,
            [
                _burden_row(
                    study="PMID:40127276",
                    cohort_stratum="all",
                    consequence_class="lof",
                    origin="de_novo",
                    comparator="mutation_model",
                    n_control_carriers=None,
                    n_controls=None,
                    control_cohorts=(),
                    expected_count=0.0128,
                    effect=297.419440226872,
                    effect_measure="enrichment_ratio",
                    effect_bound=None,
                    ci_low=None,
                    ci_high=None,
                    pvalue=1.35e-09,
                    pvalue_test="poisson",
                    pvalue_adjusted=5.58e-08,
                    pvalue_adjustment="benjamini_hochberg",
                ),
            ],
            publications={_PUBLICATION.id: _PUBLICATION, **_SIERANT_PUBLICATION},
        )
    )

    assert "chart-nointerval" in section
    assert "no interval published" in section
    # The panel is drawn, not skipped: the surviving finding is a picture.
    figures = _forest_figures(section)
    assert len(figures) == 1
    # **And 1 is on the axis.** This panel's only row is enriched 297x, so an
    # axis fitted to the values present would start above 1 and leave nothing
    # to say which side of "no enrichment" the row falls on. The null line is
    # the one drawn from the top of the plot rather than along its foot.
    null = re.search(r'<line class="chart-axis" x1="([\d.]+)" y1="16\.0"', figures[0])
    assert null is not None, "the panel drew no null line"
    assert 244.0 <= float(null.group(1)) <= 470.0


def test_an_unbounded_effect_draws_an_arrow_and_never_a_ceiling(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """`effect_bound: unbounded_above` publishes a null effect and null
    ci_high. Those rows carry the strongest signals, and `ci_low` is the
    finding -- never a blank, an em dash, or an invented ceiling."""
    section = _burden_section_text(_burden_page(tmp_path, facts_uncurated, [_burden_row()]))

    # **And the bar reaches the axis edge rather than a number nobody
    # published.** Asserted against literals in the figure's own coordinate
    # system, not against the constants that produced them (CLAUDE.md section
    # 4.38): 470.0 is the right-hand edge of the plot area and 26.0 is the
    # first row's centre line. A renderer that invented a ceiling -- say
    # `ci_high = ci_low * 10` -- would end the bar short of the edge.
    #
    # The class is `chart-arrow-open`, not `chart-arrow`: this study published
    # no correction for this row, and an arrow carries the same fill rule as a
    # circle. `"chart-arrow" in section` was the first spelling of the check
    # above it and is now a substring of the hollow class, so it is asserted
    # here as a whole attribute instead.
    figure = _forest_figures(section)[0]
    assert '<polygon class="chart-arrow-open" points="470.0,26.0 463.0,22.5 463.0,29.5"/>' in figure
    assert 'class="chart-arrow"' not in figure
    # **Scoped to this row's own bar, not to the figure.** `'x2="470.0"' in
    # figure` was the first spelling of this and a mutant giving the row a
    # ceiling of `ci_low * 10` survived it: the horizontal axis line runs to
    # the same edge and satisfied the check on the bar's behalf (CLAUDE.md
    # section 4.19). Matched at `y1="26.0"`, the first row's centre.
    bar = re.search(
        r'<line class="chart-line" x1="[\d.]+" y1="26\.0" x2="([\d.]+)" y2="26\.0"/>', figure
    )
    assert bar is not None, "the unbounded row drew no interval bar"
    assert bar.group(1) == "470.0", "the bar stopped short of the axis edge"
    # The number the study did publish is still the whole finding.
    assert "OR ∞ (95% CI 28.1–∞)" in section


def test_a_union_row_is_marked_as_one(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """`damaging` is LOF + damaging missense; `all` is syndromic +
    nonsyndromic. Stacked as sibling rows they would read as independent
    findings, which the section's prose already denies."""
    section = _burden_section_text(
        _burden_page(
            tmp_path,
            facts_uncurated,
            [
                # A union on the consequence axis and on the stratum axis at
                # once; a union on the stratum axis alone; and a row that is
                # neither, so "tag every row" fails here too.
                _burden_row(
                    cohort_stratum="all",
                    consequence_class="damaging",
                    effect=2.45,
                    effect_bound=None,
                    ci_low=1.2,
                    ci_high=8.1,
                ),
                _burden_row(
                    cohort_stratum="all",
                    consequence_class="lof",
                    effect=3.1,
                    effect_bound=None,
                    ci_low=1.4,
                    ci_high=9.0,
                ),
                _burden_row(
                    cohort_stratum="syndromic",
                    consequence_class="lof",
                    effect=4.2,
                    effect_bound=None,
                    ci_low=1.9,
                    ci_high=11.0,
                ),
                # A second study whose every row is `all cases` and which
                # publishes neither stratum below it, and no composite.
                _burden_row(
                    study="PMID:40127276",
                    cohort_stratum="all",
                    consequence_class="lof",
                    effect=14.96,
                    effect_bound=None,
                    ci_low=None,
                    ci_high=None,
                ),
                _burden_row(
                    study="PMID:40127276",
                    cohort_stratum="all",
                    consequence_class="missense_damaging",
                    effect=2.24,
                    effect_bound=None,
                    ci_low=None,
                    ci_high=None,
                ),
            ],
        )
    )

    assert "chart-union" in section
    # Two panels, in study order: the second study's `all cases` rows are the
    # whole cohort rather than the sum of two rows on the page, because it
    # publishes no syndromic or non-syndromic row at all -- which is exactly
    # PMID:40127276's shape on all 25 genes it covers. Tagging them `union`
    # would point a reader at rows that are not there, so the tag is
    # conditional on the components being in the same panel, as
    # `_composite_note` is conditional on the study reporting both.
    figures = _forest_figures(section)
    assert [figure.count('class="chart-union"') for figure in figures] == [0, 2]
    # The two components are named, so the tag is a claim a reader can check.
    assert "loss-of-function and damaging-missense rows together" in figures[1]
    assert "syndromic and non-syndromic rows together" in figures[1]


def test_the_union_key_never_promises_a_component_row_the_panel_does_not_draw(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """`_union_kinds` tests the components panel-wide; a tag is worn per row.

    Measured 2026-08-20 on the built corpus: 288 rows carry the tag, **0 of
    them with neither component on the panel** -- so the tag itself is never
    false -- but on 10 panels a tagged row has only one of its two components
    drawn. TBX5's Audain odds-ratio panel is one: `all cases · loss-of-
    function` is tagged and the key says "`all cases` is the syndromic and
    non-syndromic rows together", while the panel carries `syndromic · loss-
    of-function` and no non-syndromic one. All 10 are a cell the study
    published no row for at all -- a missing (stratum, consequence) cell is
    not a null result -- so the row below the union is part of it, not all of
    it, and the key said otherwise.

    The tag stays: `all cases · loss-of-function` really does contain the
    syndromic row beneath it, and dropping the tag would publish the two as
    independent findings, which is the defect the tag exists to prevent. It
    is the key's reconciliation instruction that has to become conditional.

    Three panels, because the corpus can only show one of the two axes: all
    10 measured cases split on the **stratum** axis, so a fixture built from
    the corpus alone would leave the consequence branch unmeasured.
    """
    enrichment: dict[str, object] = {
        "study": "PMID:40127276",
        "comparator": "mutation_model",
        "n_control_carriers": None,
        "n_controls": None,
        "control_cohorts": (),
        "expected_count": 0.42,
        "effect_measure": "enrichment_ratio",
        "pvalue_test": "poisson",
    }
    placed: dict[str, object] = {"effect_bound": None, "ci_low": 1.2, "ci_high": 8.1}
    section = _burden_section_text(
        _burden_page(
            tmp_path,
            facts_uncurated,
            [
                # Panel one, Audain odds ratios: `all cases · lof` is tagged
                # and both of its components are drawn. Nothing to warn about.
                _burden_row(cohort_stratum="all", effect=2.45, **placed),
                _burden_row(cohort_stratum="syndromic", effect=3.1, **placed),
                _burden_row(cohort_stratum="nonsyndromic", effect=4.2, **placed),
                # Panel two, Sierant enrichment ratios: `syndromic ·
                # damaging` is tagged on the consequence axis with only one of
                # its two component consequences drawn. No published panel is
                # in this state; every measured case splits on the other axis.
                #
                # **The missing component is on the panel under a different
                # stratum**, and no row here is `all cases`, so the stratum
                # union is off and this measures the consequence branch alone.
                # Without it, dropping `_fully_split`'s `other.cohort_stratum
                # == row.cohort_stratum` -- checking the union across strata
                # rather than within one -- survived all 1,050 tests: this
                # panel's rows all shared a stratum, so holding it fixed or
                # not made no difference to any of them (CLAUDE.md 4.36).
                _burden_row(
                    cohort_stratum="syndromic",
                    consequence_class="damaging",
                    effect=5.0,
                    **placed,
                    **enrichment,
                ),
                _burden_row(
                    cohort_stratum="syndromic",
                    consequence_class="lof",
                    effect=6.0,
                    **placed,
                    **enrichment,
                ),
                _burden_row(
                    cohort_stratum="nonsyndromic",
                    consequence_class="missense_damaging",
                    effect=4.0,
                    **placed,
                    **enrichment,
                ),
                # Panel three, Sierant odds ratios: TBX5's real shape -- a
                # stratum union with the non-syndromic side missing.
                _burden_row(study="PMID:40127276", cohort_stratum="all", effect=7.0, **placed),
                _burden_row(
                    study="PMID:40127276", cohort_stratum="syndromic", effect=8.0, **placed
                ),
            ],
            publications={_PUBLICATION.id: _PUBLICATION, **_SIERANT_PUBLICATION},
        )
    )

    # PMID:40127276 sorts first and `enrichment_ratio` before `odds_ratio`.
    consequence_short, stratum_short, complete = _forest_figures(section)
    warning = "a tagged row is <strong>not fully split</strong>"
    assert warning in consequence_short
    assert warning in stratum_short
    assert warning not in complete
    # The tag itself is untouched on every one of them: the row really does
    # contain the row below it, and saying otherwise would publish two
    # overlapping rows as independent findings.
    assert all(
        figure.count('class="chart-union"') >= 1
        for figure in (consequence_short, stratum_short, complete)
    )


def test_an_arrow_carries_the_same_correction_fill_as_a_circle(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """An arrow is a marker, so its fill has to mean what a circle's means.

    `_forest_estimate` bypassed `_mark_class` for both arrow branches and
    drew every one of them in `chart-arrow`, whose declaration was
    byte-identical to `chart-point`'s -- the fill that says "survived this
    study's own correction". Measured 2026-08-20 on the built corpus: **205
    of the 854 plotted marks were arrows and not one had survived any
    correction.** 30 were corrected and failed, with q as high as 1.0, and
    175 came from a study that published no correction at all; 13 of the
    corrected-and-failed sat on a ClinGen `definitive` gene.

    GATA4 was the worst of them and is the shape of the fixture below: three
    marks on one panel, all solid, no hollow mark anywhere to read them
    against, and the third was `p 1 · q 1`.

    The fixture carries all three fills in one panel, because a fixture whose
    rows share the value under test measures nothing (CLAUDE.md section 4.36)
    -- and the filled branch has **no row in the committed corpus at all**, so
    the surviving unbounded row here is the only thing that proves it is
    reachable rather than dead.
    """
    page = _burden_page(
        tmp_path,
        facts_uncurated,
        [
            # Unbounded above, no correction published: hollow.
            _burden_row(cohort_stratum="all", consequence_class="lof"),
            # Unbounded above and the study's own correction survived:
            # filled. No published row is in this state.
            _burden_row(
                cohort_stratum="syndromic",
                consequence_class="missense_damaging",
                pvalue_adjusted=0.0037,
                pvalue_adjustment="benjamini_hochberg",
            ),
            # An effect of exactly zero -- the other arrow branch, at the
            # opposite edge -- corrected and failed: hollow.
            _burden_row(
                cohort_stratum="nonsyndromic",
                consequence_class="lof",
                effect=0.0,
                effect_bound=None,
                ci_low=0.0,
                ci_high=3.2,
                pvalue_adjusted=1.0,
                pvalue_adjustment="benjamini_hochberg",
            ),
            # The study's own negative control, surviving and unbounded:
            # muted, never in the result colour.
            _burden_row(
                cohort_stratum="all",
                consequence_class="synonymous",
                pvalue_adjusted=0.0041,
                pvalue_adjustment="benjamini_hochberg",
            ),
        ],
    )

    figure = _forest_figures(_burden_section_text(page))[0]
    assert figure.count('class="chart-arrow-open"') == 2
    assert figure.count('class="chart-arrow"') == 1
    assert figure.count('class="chart-control-arrow"') == 1
    # And the hollow form is visibly hollow rather than a second name for the
    # same paint. Asserted on the published stylesheet -- the class is only
    # worth anything if it resolves to a different fill, and `chart-arrow`'s
    # declaration was byte-identical to `chart-point`'s for three releases.
    assert ".chart-arrow-open { fill: var(--bg); stroke: var(--link); stroke-width: 1.6; }" in page
    assert ".chart-control-arrow { fill: var(--muted); }" in page
    # A key that names only the hollow mark leaves the reader to infer the
    # complement, and for every arrow on the site that inference was wrong.
    assert "A <strong>filled</strong> marker is a row that <strong>survived</strong>" in figure


def test_each_fill_in_the_key_is_a_fill_the_reader_can_see(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """`_forest_caption` is conditional clause by clause; the fills were not.

    A legend entry for a glyph that is not on the panel sends a reader
    hunting for it, which is `_POOLING_NOTICE`'s recorded lesson (CLAUDE.md
    section 4.27). Measured 2026-08-20 on the built corpus: **11 of the 137
    panels printed the hollow clause with no hollow marker drawn**, across 10
    genes -- CHD7, GATA4, KMT2D, MYH7, NKX2-5, NODAL, RBFOX2, RNF40, WDR5,
    ZIC3, six of them ClinGen `definitive`. On 9 of the 11 the panel did draw
    arrows, so the key affirmatively told the reader that everything visible
    had failed its correction when the arrows carried no such claim.

    Both directions are asserted, from two panels that differ in exactly this
    -- an all-surviving panel and an all-failing one.
    """
    survived: dict[str, object] = {
        "pvalue_adjusted": 0.0037,
        "pvalue_adjustment": "benjamini_hochberg",
    }
    section = _burden_section_text(
        _burden_page(
            tmp_path,
            facts_uncurated,
            [
                _burden_row(
                    cohort_stratum="all",
                    effect=2.45,
                    effect_bound=None,
                    ci_low=1.2,
                    ci_high=8.1,
                    **survived,
                ),
                _burden_row(
                    study="PMID:40127276",
                    cohort_stratum="all",
                    effect=3.1,
                    effect_bound=None,
                    ci_low=1.4,
                    ci_high=9.0,
                ),
            ],
            publications={_PUBLICATION.id: _PUBLICATION, **_SIERANT_PUBLICATION},
        )
    )

    # PMID:40127276 sorts before PMID:42230622, so the failing panel is first.
    failing, surviving = _forest_figures(section)
    hollow = "A <strong>hollow</strong> marker is a row that <strong>did not survive</strong>"
    filled = "A <strong>filled</strong> marker is a row that <strong>survived</strong>"
    assert filled in surviving and hollow not in surviving
    assert hollow in failing and filled not in failing
    # The glyph each clause names really is the one drawn beside it.
    assert 'class="chart-point-open"' not in surviving
    assert 'class="chart-point"' not in failing


def test_marker_fill_encodes_survival_of_the_studys_own_correction(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """Filled where the study's own correction was survived, hollow where not.

    Never a correction this atlas computed. A row whose study published no
    correction at all is hollow for the same honest reason -- nothing says it
    survived one.

    The fixture carries one surviving and one non-surviving row, because a
    fixture whose rows all share the value under test measures nothing.
    """
    page = _burden_page(
        tmp_path,
        facts_uncurated,
        [
            _burden_row(
                cohort_stratum="all",
                effect=2.45,
                effect_bound=None,
                ci_low=1.2,
                ci_high=8.1,
                pvalue_adjusted=0.0037,
                pvalue_adjustment="benjamini_hochberg",
            ),
            _burden_row(
                cohort_stratum="syndromic",
                effect=3.1,
                effect_bound=None,
                ci_low=1.4,
                ci_high=9.0,
                pvalue_adjusted=0.42,
                pvalue_adjustment="benjamini_hochberg",
            ),
            _burden_row(
                cohort_stratum="nonsyndromic",
                effect=4.2,
                effect_bound=None,
                ci_low=1.9,
                ci_high=11.0,
            ),
        ],
    )
    section = _burden_section_text(page)

    assert "chart-point-open" in section
    assert 'class="chart-point"' in section
    figure = _forest_figures(section)[0]
    # One survived; the other two did not, and one of those two is hollow for
    # the second honest reason -- its study published no correction at all.
    assert figure.count('class="chart-point"') == 1
    assert figure.count('class="chart-point-open"') == 2
    assert "did not survive" in figure
    # And the hollow form resolves to a visibly different paint, asserted on
    # the published stylesheet exactly as the arrow half asserts its own
    # (CLAUDE.md section 4.31: a guard on one of a pair is evidence about
    # one). `.chart-arrow-open` was pinned here and `.chart-point-open` was
    # not, in the release that cites this rule -- and pointing it at
    # `var(--link)` leaves all 604 hollow circles on the site indistinguishable
    # from the 250 filled ones.
    assert ".chart-point-open { fill: var(--bg); stroke: var(--link); stroke-width: 1.6; }" in page


def test_a_gene_with_no_plottable_row_gets_a_sentence_and_no_chart(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """D42. Measured: KLF13, CFC1 and MYH11 carry only CNV rows with no effect
    measure at all; CRIPTO has no burden row. Four of 92."""
    section = _burden_section_text(
        _burden_page(
            tmp_path,
            facts_uncurated,
            [
                _burden_row(
                    study="PMID:34324492",
                    variant_class="cnv_deletion",
                    cohort_stratum="all",
                    consequence_class="all_coding",
                    effect=None,
                    effect_measure=None,
                    effect_bound=None,
                    ci_low=None,
                    ci_high=None,
                    pvalue=0.0068,
                    pvalue_test="binomial",
                    pvalue_adjusted=0.99,
                    pvalue_adjustment="familywise_permutation",
                ),
            ],
            publications={**_OTHER_PUBLICATIONS},
        )
    )

    assert "<svg" not in section
    assert "no effect measure" in section
    # And the figures it does have stay in plain sight rather than folding
    # behind a summary that promises a picture there is none of.
    assert "<table" in re.sub(r"<details.*?</details>", "", section, flags=re.S)


def test_the_forest_carries_no_summary_diamond(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """D33: no pooled statistic, ever. These cohorts overlap -- DDD
    contributes cases to more than one cited paper -- so a pooled estimate
    would count the same children twice. The forest idiom invites a diamond
    and a reader expects one; its absence is deliberate and captioned.

    **Captioned in two places, and only one of them was guarded.** The
    figure's `<title>` is what a screen reader announces before any of the
    marks, and it makes the same claim in its own words. Measured 2026-08-20,
    rewriting it to announce "and a pooled summary" survived the whole suite:
    the assertion below is scoped to `<figure class="forest">`, which contains
    the caption as well, so the caption answered for the title. A reader using
    assistive technology would have been told this atlas publishes the one
    statistic D33 forbids -- and told it first.
    """
    section = _burden_section_text(
        _burden_page(
            tmp_path,
            facts_uncurated,
            [_burden_row(effect=2.45, effect_bound=None, ci_low=1.2, ci_high=8.1)],
        )
    )

    assert "chart-diamond" not in section
    assert "no pooled" in section
    # Scoped to the panel that invites the diamond. `_POOLING_NOTICE` already
    # carries "no pooled statistic across studies" as the matrix caption a
    # screen above, so the section-wide assertion above passes with or without
    # a word from this figure (CLAUDE.md section 4.19).
    figure = _forest_figures(section)[0]
    assert "no pooled" in figure
    # And the same scoping again, one level in: the caption and the title are
    # two sentences and a check over the figure is answered by either.
    announced = re.search(r"<title>(.*?)</title>", figure, re.S)
    assert announced is not None, "the figure announces nothing to a screen reader"
    assert "no pooled summary" in announced.group(1)
    assert "a pooled summary" not in announced.group(1).replace("no pooled summary", "")


def test_a_forest_never_replaces_the_figures_it_summarises(
    tmp_path: Path, facts_uncurated: dict[str, GeneFacts]
) -> None:
    """The evidence-loss mutant.

    A picture that took the table with it would remove every exact count,
    interval and p-value from the HTML -- this repository's characteristic
    defect. The table folds; it does not leave.
    """
    section = _burden_section_text(
        _burden_page(
            tmp_path,
            facts_uncurated,
            [_burden_row(effect=2.45, effect_bound=None, ci_low=1.2, ci_high=8.1)],
        )
    )

    assert "<svg" in section
    folded = re.search(r'<details class="study-figures">.*?</details>', section, re.S)
    assert folded is not None, "the study's table must be folded, not deleted"
    assert "<table" in folded.group(0)
    assert "5 carriers / 1,471" in folded.group(0)
    assert "OR 2.45 (95% CI 1.2–8.1)" in folded.group(0)
    # A section-wide `<details` check would not be enough: the reading notes
    # and the how-this-study-counted block are both `<details>` a few lines
    # above (CLAUDE.md section 4.19). Nor is `"<table" not in unfolded`: the
    # evidence matrix is a table too, and it is outside the fold on purpose.
    # Asserted on a count only this study's own table renders -- the matrix
    # cell carries the effect and one statistic and no count at all.
    unfolded = re.sub(r"<details.*?</details>", "", section, flags=re.S)
    assert "5 carriers / 1,471" not in unfolded


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
    *stages: tuple[str, float | None | EllipsisType],
    phase_ids: tuple[str, ...] = ("atrial_septum_morphogenesis",),
    unplaced_reason: str = "below_detection_floor",
) -> ExpressionProfile:
    """A heart series over `stages`, in the order given.

    Three states, because the page owes a reader three different sentences
    and `_tissue_medians`' docstring insists the last two "must not render the
    same way":

    * a `float` -- a median this dataset placed against its percentile grid;
    * `None` -- a stage sampled and **not placed**, for `unplaced_reason`,
      which `_trajectory` ticks on the axis;
    * `...` -- a stage this dataset has **no heart row for at all**, which
      gets no tick and nothing else, and is why the stage is still in
      `entry["stages"]` and still occupies an x position.

    The order given is the published order -- `gene_expression_profiles` emits
    stages by `Stage.order`, so a fixture's own sequence is the trajectory's.

    `unplaced_reason` is a parameter rather than a constant because the reason
    changes what a page is entitled to say: below a floor is a low reading
    about the gene, and a missing percentile grid is a hole in the reference
    that says nothing about the gene at all.
    """
    entries = [
        _stage_entry(
            stage=token,
            phase=_phase_info(phase_ids=phase_ids),
            tissues=()
            if median is ...
            else (
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

    **And the sentence that replaces the refused chart has to be true**, which
    is a second assertion and not the same one. `_trajectory` returns `""` for
    three different reasons and this half exercises the third; a caller that
    routes all three into one sentence publishes, of a gene placed at every
    stage on the page below, that no measurement there is placed at all. That
    is what shipped until 2026-08-20: the refusal was guarded, the replacement
    was not.
    """
    profile = _heart_series(("4wpc", 5.0), ("5wpc", 40.0), ("6wpc", 275.0))
    section = _expression_section_text(
        _expression_page(tmp_path, {GATA4: profile}, _HEART_DATASET, phases=_cardiac_phases())
    )

    assert "<svg" in section
    assert "chart-band" in section, "a trajectory with no phase band is the re-plot D43 forbids"

    # The same three measurements with no curated phase vocabulary behind them
    # are exactly the source's own curve, and are refused.
    bare = _dataset_lede(
        _expression_section_text(
            _expression_page(tmp_path, {GATA4: profile}, _HEART_DATASET, phases=None)
        )
    )
    assert "<svg" not in bare
    # Every one of these three stages is placed, and the table below prints
    # its percentile. Neither of the other two refusals' sentences may be
    # published here.
    assert "no measurement there is placed" not in bare
    assert "detection floor at every" not in bare
    assert "curated cardiac phase" in bare, "the refusal has to name the reason it refused"
    # And it is never left dangling. The reason list is empty in this branch,
    # so a sentence built to end in one ends "&mdash; ." -- an em dash, a
    # space and a full stop, published on 55 pages by the mutant above.
    assert "&mdash; ." not in bare

    # The third fixture, and the one whose sentence was worst. A series with
    # a below-floor stage *between* two placed ones has exactly one recorded
    # gap reason, so a sentence keyed on the reason set alone reads it as
    # "every stage is below the floor" -- in bold, on a gene whose page
    # prints two percentiles immediately below it. "Nothing placed" and "one
    # gap, below the floor" are different facts and only the first earns that
    # sentence.
    mixed = _dataset_lede(
        _expression_section_text(
            _expression_page(
                tmp_path,
                {GATA4: _heart_series(("4wpc", 5.0), ("5wpc", None), ("6wpc", 275.0))},
                _HEART_DATASET,
                phases=None,
            )
        )
    )
    assert "<svg" not in mixed
    assert "detection floor at every" not in mixed
    assert "no measurement there is placed" not in mixed


def _trajectory_figure(section: str) -> str:
    """The phase-banded trajectory alone, out of the section it sits in.

    Sliced rather than asserted page-wide for the reason `_validity_table` is
    (CLAUDE.md section 4.19): the organ small multiples are `<svg>`s on the
    same page, carry their own `chart-label` text and print their own axis
    range in words, so a section-wide check on either cannot tell which
    picture answered it. The trajectory is the 560-wide one.
    """
    found = re.findall(r'<svg class="chart" viewBox="0 0 560 150".*?</svg>', section, flags=re.S)
    assert len(found) == 1, f"expected one trajectory, found {len(found)}"
    return found[0]


def test_the_trajectory_says_what_its_vertical_axis_spans(tmp_path: Path) -> None:
    """A curve with no numbers on it is a shape, not a measurement.

    The trajectory contained **zero `<text>` elements** -- no ticks, no
    labels, no numbers -- and its vertical axis is fitted to this gene's own
    placed heart values, so every gene's curve filled the same 104 px whatever
    it spanned. Measured 2026-08-20 on the built corpus: PKD1L1 runs 1 to 2
    tpm and TBX20 runs 1 to 526 tpm, and the two were drawn as the same
    full-height excursion; 8 of the 85 charted genes span under 5-fold and 39
    span 20-fold or more. The small-multiples caption one screen below
    condemns exactly this -- "a per-organ axis would rescale every line to its
    own range and hide exactly that difference" -- while printing its own
    range.

    Two genes with different spans, because a fixture whose genes share the
    value under test measures nothing (CLAUDE.md section 4.36): with one
    gene, a renderer that labelled the axis with any constant would pass.

    The numbers are `_fmt`'s, not `charts.coordinate`'s. `coordinate` is the
    *pixel* formatter -- it fixes to a tenth, so it would publish a median of
    275 tpm as `275.0` and one of 4,761 as `4761.0` beside `_fmt`'s `4,761`
    in the table below. The literals here are what `_fmt` produces and what
    `coordinate` does not.
    """
    wide = _expression_section_text(
        _expression_page(
            tmp_path,
            {GATA4: _heart_series(("4wpc", 5.0), ("5wpc", 40.0), ("6wpc", 275.0))},
            _HEART_DATASET,
            phases=_cardiac_phases(),
        )
    )
    narrow = _expression_section_text(
        _expression_page(
            tmp_path,
            {GATA4: _heart_series(("4wpc", 1.0), ("5wpc", 1.5), ("6wpc", 2.0))},
            _HEART_DATASET,
            phases=_cardiac_phases(),
        )
    )

    # Pinned as whole elements in the figure's own coordinate system: 42.0 is
    # the right edge of the margin `_PLOT_LEFT` reserves, 15.4 is the top of
    # the value scale and 119.4 its foot, and `text-anchor="end"` is what
    # keeps a four-digit label off the plot.
    top = '<text class="chart-label" x="42.0" y="15.4" text-anchor="end">'
    foot = '<text class="chart-label" x="42.0" y="119.4" text-anchor="end">'
    assert f"{top}275</text>" in _trajectory_figure(wide)
    assert f"{foot}5</text>" in _trajectory_figure(wide)
    assert f"{top}2</text>" in _trajectory_figure(narrow)
    assert f"{foot}1</text>" in _trajectory_figure(narrow)

    # And the caption says the same span in words, as the sibling figure's
    # already does -- a reader who cannot see the axis labels still gets it.
    #
    # Stated as the **axis's** range and not as the medians'. `_axis_bounds`
    # widens a flat series by half a decade either side, because `LogScale`
    # refuses a zero-span axis, and 3 of the 85 charted genes are flat: for
    # TFAP2B, placed at 1 tpm and nowhere else, "median abundance in whole
    # heart, 0.316 to 3.16 tpm" attributes to its medians a spread neither
    # end of which anyone measured. The third fixture is that gene's shape.
    assert "the axis runs 5 to 275 and is fitted to" in wide
    assert "the axis runs 1 to 2 and is fitted to" in narrow

    flat = _expression_section_text(
        _expression_page(
            tmp_path,
            {GATA4: _heart_series(("4wpc", 1.0), ("5wpc", 1.0), ("6wpc", 1.0))},
            _HEART_DATASET,
            phases=_cardiac_phases(),
        )
    )
    assert "the axis runs 0.316 to 3.16 and is fitted to" in flat
    assert "median abundance in whole heart, tpm, on a log scale" in flat


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


def test_the_band_caption_names_only_the_breaks_this_figure_actually_draws(
    tmp_path: Path,
) -> None:
    """Two ways a trajectory's line can break, and the caption named one of
    them unconditionally -- including on the 55 pages that draw no tick at all.

    Measured 2026-08-20 on the built corpus: 85 pages carry a trajectory and
    every one of them said "the stage is ticked on the axis instead"; **55
    drew zero ticks**, and 55 of the 59 with a visibly broken line were among
    them. TBX5 is one -- two polylines, a visible gap, and the only sentence
    about gaps described a mark that is not on the page.

    The second cause was never named at all: this dataset has **no row** for
    an organ at some stages, identically for every gene (heart has none at
    `school age child` or `elderly`, forebrain none at `6 week post
    conception`, testis none at `neonate` or `school age child`), and that
    breaks the line with no tick under it. `_tissue_medians`' own docstring
    insists the two "must not render the same way"; in the picture they did.

    Four fixtures, one per corner: a fixture where every case shares the
    value under test measures nothing (CLAUDE.md section 4.36).

    **And a fifth, on the axis the four could not see.** A tick is drawn for
    every stage this dataset sampled and did not place, whatever the recorded
    reason; the clause explaining it was conditioned on `value is None` and
    so named the detection floor for all of them. A stage unplaced because no
    percentile grid was published is ticked and captioned as a low reading,
    three lines above a table cell reading "no complete percentile grid is
    published for this organ at this stage" -- the page contradicting itself,
    and in the direction that invents a measurement. `_heart_series` has
    carried `unplaced_reason` since the sentence tier was written and no test
    had pointed it at a caption.
    """

    def caption(
        *stages: tuple[str, float | None | EllipsisType],
        reason: str = "below_detection_floor",
    ) -> str:
        return _expression_section_text(
            _expression_page(
                tmp_path,
                {GATA4: _heart_series(*stages, unplaced_reason=reason)},
                _HEART_DATASET,
                phases=_cardiac_phases(),
            )
        )

    tick = "the stage is ticked on the axis instead"
    no_row = "no row for this organ at that stage"

    # Ticked and unbroken: the floor clause, and nothing about a missing row.
    ticked = caption(("4wpc", None), ("5wpc", 5.0), ("6wpc", 40.0), ("7wpc", 275.0))
    assert tick in ticked
    assert no_row not in ticked

    # TBX5's shape: a line broken by a stage this dataset has no row for, and
    # not one tick anywhere. The floor clause described a mark that is absent.
    unsampled = caption(("4wpc", 5.0), ("5wpc", 40.0), ("6wpc", ...), ("7wpc", 275.0))
    assert tick not in unsampled
    assert no_row in unsampled
    # And it does not point back at a clause this figure does not carry.
    # "A break with no tick under it is different again" beside no tick at
    # all is the same defect one sentence further on.
    assert "different again" not in unsampled

    # Both causes at once -- the caption owes a reader both sentences.
    both = caption(
        ("4wpc", 5.0),
        ("5wpc", 40.0),
        ("6wpc", None),
        ("7wpc", 275.0),
        ("8wpc", ...),
        ("9wpc", 12.0),
    )
    assert tick in both
    assert no_row in both
    assert "different again" in both

    # Neither: an unbroken, fully placed series says nothing about either.
    clean = caption(("4wpc", 5.0), ("5wpc", 40.0), ("6wpc", 275.0))
    assert tick not in clean
    assert no_row not in clean

    # The same ticked shape, unplaced for a reason that is not the floor. The
    # mark is drawn either way, so the caption still owes the reader an
    # explanation of it -- it may not be the floor's.
    grid = caption(
        ("4wpc", None), ("5wpc", 5.0), ("6wpc", 40.0), ("7wpc", 275.0), reason="no_quantile_grid"
    )
    assert grid.count('class="chart-absent"') == 1, "the tick is drawn whatever the reason"
    # Scoped to the caption: `_percentile_cell` prints the same clause in the
    # stage table below, from the same map and deliberately so, and answers a
    # section-wide check on it whatever the caption said.
    caption_only = _slice_between(grid, "Left to right", ("</p>",))
    assert "Below this dataset's detection floor" not in caption_only
    assert "no complete percentile grid is published for this organ at this stage" in caption_only
    # And the mark is still accounted for rather than left on the axis with
    # nothing said about it, which is the defect the floor clause was written
    # for in the first place.
    assert "ticked" in caption_only

    # And the tick is placed off the bottom of the scale rather than at a
    # value on it. Measured: KIF20A's axis starts at 42 tpm against a 1 tpm
    # floor, so its 6 ticks sat where the axis reads 38.6 tpm; the words are
    # what stop a reader interpolating them.
    assert "below the foot of the scale, not at a value on it" in ticked


def test_the_line_never_spans_a_stage_the_atlas_did_not_place(tmp_path: Path) -> None:
    """A below-floor stage breaks the line; it is not drawn straight over, is
    ticked rather than plotted at zero, and its tick never lands on a real
    measurement.

    Measured 2026-08-20 on the committed corpus: 20 of the 85 charted genes
    carry a below-floor stage strictly between two plotted points, and TBX1
    has 11 of them against 5 placed stages. A single line through the placed
    points runs smooth and high across exactly the stages where this dataset
    measured below its own floor -- the picture asserting a continuity the
    figures deny.

    The fixture has two runs and a lone leading point, so "one line through
    everything", "splits but drops a run" and "splits correctly" are three
    distinguishable outcomes rather than two.

    **Absorbs `test_a_stage_below_the_floor_is_an_axis_tick_and_never_a_zero`,
    whose fixture this one is a strict superset of.** Measured 2026-08-20:
    all three mutants that test could kill are killed here too, so it was
    paying for a build without being the unique killer of anything. Its
    rationale is the record and moves with it -- *below the floor*, *not
    sampled* and *zero* are three different claims: plotting a below-floor
    stage at zero puts it on the axis as the lowest measurement this atlas
    made, and leaving a gap makes it indistinguishable from a stage nobody
    sampled. The counts are asserted rather than the class name alone,
    because a renderer that ticked every stage, or that ticked the
    below-floor stage *and* plotted it, satisfies a bare `in` check.

    **And `_FLOOR_GAP` is what keeps the tick off the lowest real point.**
    That constant's whole stated purpose is to stop the atlas's weakest
    statement about a figure -- "we do not vouch for this one" -- landing on
    the same pixel row as its lowest real one. Nothing guarded it: the lowest
    placed point always sits at exactly `_PLOT_BOTTOM - _FLOOR_GAP`, so
    ticking a below-floor stage there instead of on the axis puts the two
    marks on one row on every charted page -- cy 116.0 on HGNC:10249 --
    which is the conflation `_percentile_cell` refuses in words.
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

    # No tick shares a row with a measurement. Read off the published `cy`
    # rather than recomputed, so the assertion cannot agree with the renderer
    # by sharing its arithmetic.
    figure = _trajectory_figure(section)
    rows = {
        css: {cy for cy in re.findall(rf'<circle class="{css}" cx="[^"]*" cy="([^"]*)"', figure)}
        for css in ("chart-point", "chart-absent")
    }
    assert rows["chart-absent"] == {"124.0"}, "a tick belongs on the axis"
    assert "116.0" in rows["chart-point"], "the lowest placed point sits one _FLOOR_GAP up"
    assert not rows["chart-absent"] & rows["chart-point"], (
        "a below-floor tick shares a row with a real measurement, which is the "
        "conflation _FLOOR_GAP exists to prevent"
    )


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


@pytest.mark.parametrize(
    ("reason", "present", "absent"),
    [
        pytest.param(
            "below_detection_floor",
            ("detection floor at every", "not evidence that", "heart", "3 stages sampled"),
            (),
            id="below this dataset's own floor",
        ),
        pytest.param(
            "no_quantile_grid",
            ("no complete percentile grid is published for this organ at this stage",),
            ("detection floor at every",),
            id="a gap in the reference",
        ),
    ],
)
def test_nothing_placed_says_which_of_the_two_facts_it_is(
    tmp_path: Path, reason: str, present: tuple[str, ...], absent: tuple[str, ...]
) -> None:
    """D42, tier 4, and the reason it is not cosmetic.

    Measured 2026-08-20: seven published genes are below the floor in heart at
    every one of 19 stages -- SEMA3E, ZIC3, USP44, DAW1, FGF8, GDF1, NODAL.
    ZIC3 and NODAL are ClinGen definitive. An empty chart beside a definitive
    chip asserts the thing the bulk-dilution caveat denies.

    **"Nothing placed" is not one fact, and only one of them is about the
    gene.** `below_detection_floor` says this dataset measured the gene under
    its own floor, which is what earns the dilution argument beside it. Every
    other gap -- no percentile grid for this organ and stage, no floor
    declared, the dataset unregistered -- is a hole in the *reference*, and
    reporting it as a low reading asserts a measurement this dataset never
    made. The two cases were two tests with identical bodies differing only
    in this parameter, which the second's own docstring said; they are one
    parametrised test.

    Asserted against literals, not against the constant that produced them: a
    test that imports the string it asserts on compares the module to itself
    and passes for any rewording.

    Sliced to the dataset block's own lede before asserting, because
    `_BULK_DILUTION_NOTICE` carries "not evidence that" over every gene page
    on this site, charted or not.
    """
    profile = _heart_series(("4wpc", None), ("5wpc", None), ("6wpc", None), unplaced_reason=reason)
    section = _dataset_lede(
        _expression_section_text(
            _expression_page(tmp_path, {GATA4: profile}, _HEART_DATASET, phases=_cardiac_phases())
        )
    )

    assert "<svg" not in section
    for phrase in present:
        assert phrase in section
    for phrase in absent:
        assert phrase not in section


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


# --- Organ small multiples (D40/D43) -----------------------------------------
#
# The fixtures below all sample at least one cardiac and one non-cardiac organ,
# at genuinely different magnitudes. CLAUDE.md section 4.36 is why that is
# spelled out: a fixture whose organs all share the value under test measures
# nothing, and it has cost this repository four separate defects. A panel where
# every organ rendered `chart-cardiac` and one where none did would be
# indistinguishable on a heart-only fixture.


def _organ_panel(
    *organs: tuple[str, tuple[float | None | EllipsisType, ...]],
    stages: tuple[str, ...] = ("4wpc", "5wpc", "6wpc"),
    highest_in: str | None = "heart",
    phase: PhaseInfo | None = None,
    unplaced_reason: str = "below_detection_floor",
) -> ExpressionProfile:
    """Several organs measured across the same stages, each at its own magnitude.

    The same three states `_heart_series` carries, for the same reason. A
    `float` is a median this dataset placed against its percentile grid; a
    `None` is a stage this organ was sampled at and not placed, below the
    detection floor; `...` is a stage this dataset has **no row for that
    organ** at, which is a third fact and breaks a panel's line with nothing
    to show for it -- E-MTAB-6814's real shape, identically for every gene.

    Every stage carries a real `specificity`, unlike `_heart_series`: tau is
    what this panel exists to make visible, and `_small_multiples` draws
    nothing for a dataset that reports none.

    `phase` defaults to `_phase_info()`, whose id is not in `_cardiac_phases()`
    -- so a panel built here draws no trajectory unless a caller asks for a
    phase the vocabulary actually names. That default is what every test
    written before the trajectory landed relies on; pass one to get bands.

    `unplaced_reason` is a parameter for the reason it is one on
    `_heart_series`: a `None` here is a stage this dataset sampled and did
    not place, and *why* decides what the caption beside the panel is
    entitled to say about it. Below a floor is a low reading about the gene;
    a missing percentile grid is a hole in the reference and says nothing
    about the gene at all.
    """
    names = tuple(name for name, _ in organs)
    entries = tuple(
        _stage_entry(
            stage=token,
            phase=phase,
            specificity=_specificity(tissues=names, highest_in=highest_in),
            tissues=tuple(
                _tissue_entry(
                    tissue=name,
                    median=values[index] if isinstance(values[index], float) else 0.3,
                    unit="tpm",
                    placement=_placement() if isinstance(values[index], float) else None,
                    not_placed_reason=(
                        None if isinstance(values[index], float) else unplaced_reason
                    ),
                )
                for name, values in organs
                if values[index] is not ...
            ),
        )
        for index, token in enumerate(stages)
    )
    return _expression_profile((_dataset_profile_entry(stages=entries),))


def _spark_caption_text(section: str) -> str:
    """The small-multiples caption alone, out of the section it sits in.

    Scoped for the reason `_validity_table` is (CLAUDE.md section 4.19). Every
    gap clause this caption can print is also printed by `_percentile_cell` in
    the stage table a few lines below -- they share `_PLACEMENT_GAP_CLAUSE` on
    purpose, so that a caption and a cell cannot spell one gap two ways -- and
    a section-wide check on any of them is answered by the table whether or
    not the caption said anything at all. Measured 2026-08-20: deleting the
    caption's whole non-floor clause left a section-wide assertion green.
    """
    return _slice_between(section, "One panel per organ with a placed measurement", ("</p>",))


def _spark_panels(section: str) -> dict[str, str]:
    """Each organ's own `<figure class="spark">`, keyed by the caption naming it.

    Scoped rather than page-wide, for the reason `_validity_table` is scoped
    (CLAUDE.md section 4.19): a section-wide `chart-cardiac` check passes on
    the strength of a *different* organ's panel, which is the shape of defect
    this repository has shipped twice.
    """
    keyed: dict[str, str] = {}
    for panel in re.findall(r'<figure class="spark".*?</figure>', section, re.S):
        caption = re.search(r"<figcaption>(.*?)</figcaption>", panel, re.S)
        assert caption is not None, "every panel names the organ it draws"
        keyed[caption.group(1)] = panel
    return keyed


_SPARK_PANEL = _organ_panel(
    ("heart", (50.0, 150.0, 275.0)),
    ("kidney", (0.8, 1.0, 1.5)),
    ("liver", (0.5, 1.2, 2.0)),
)


def test_small_multiples_share_one_axis_across_organs(tmp_path: Path) -> None:
    """The shared axis is the whole point, and it has to be the *right* axis.

    Per-organ axes would rescale each line to its own range and show seven
    similar-looking traces; one axis is what makes heart-preference visible
    instead of asserted once per stage.

    **"Every panel declares the same maximum" is satisfied by every wrong
    maximum too**, which is what `len(set(scales)) == 1` was measuring.
    Measured 2026-08-20: emitting `coordinate(scale.low)` in place of
    `coordinate(scale.high)` published `data-scale-high="0.5"` on all five
    TBX5 panels and passed -- an axis attribute naming the smallest value on
    the page as its ceiling, which is the one number a consumer reading it
    would use to rescale. So the value is asserted, not only its uniformity.

    **And the fixture's tallest organ is not its first**, which is the axis
    the assertion could otherwise not see. `drawn` is sorted, so with heart
    both first and tallest, fitting the axis to `drawn[:1]` -- one organ's
    range, the very thing the caption says a per-organ axis would do --
    leaves `high` unchanged and survives. Here liver is tallest and heart is
    first, so the shared maximum can only come from pooling every organ.

    The printed range is asserted beside it because `high` alone does not pin
    `low`: the caption is where a reader learns what the axis spans, and it
    is the only place `low` is published at all.
    """
    # heart first alphabetically, liver tallest, kidney lowest: no single
    # organ carries both ends, so neither bound can come from one panel.
    profile = _organ_panel(
        ("heart", (50.0, 150.0, 180.0)),
        ("kidney", (0.5, 1.0, 1.5)),
        ("liver", (2.0, 40.0, 275.0)),
    )
    section = _expression_section_text(_expression_page(tmp_path, {GATA4: profile}, _HEART_DATASET))

    scales = re.findall(r'data-scale-high="([^"]+)"', section)
    assert len(scales) > 1
    assert len(set(scales)) == 1, "each organ drew its own axis; heart-preference vanishes"
    # One panel per organ, so a renderer that drew only the tallest -- and
    # therefore trivially shares one axis with itself -- fails here.
    assert len(scales) == 3
    # The tallest placed median anywhere in the dataset, in `coordinate`'s
    # fixed-precision form, and nothing else.
    assert scales[0] == "275.0"
    # Both ends, in `_fmt`'s form, from the caption that publishes them.
    assert "(0.5 to 275 tpm, log scale)" in _spark_caption_text(section)


def test_only_the_declared_cardiac_organ_is_marked_as_one(tmp_path: Path) -> None:
    """Which organ is cardiac is this atlas's declaration, not the source's,
    and it is what makes the panel non-redundant under D43.

    Both directions, on one fixture and one page build. A fixture whose
    organs all render identically measures nothing -- this repository has
    shipped four defects of exactly that shape -- and the positive half was
    its own test until 2026-08-20, when it was measured not to be the unique
    killer of anything: `test_a_spark_line_never_spans_a_stage_the_atlas_did_
    not_place` kills the all-control mutant too, matching on the same class
    name. Only the negative half is unique, and it is worth nothing without
    the positive beside it.
    """
    section = _expression_section_text(
        _expression_page(tmp_path, {GATA4: _SPARK_PANEL}, _HEART_DATASET)
    )
    panels = _spark_panels(section)

    assert set(panels) == {"heart", "kidney", "liver"}
    assert "chart-cardiac" in panels["heart"]
    assert "chart-cardiac" not in panels["liver"]
    assert "chart-control" in panels["liver"]
    assert "chart-cardiac" not in panels["kidney"]


def test_small_multiples_are_absent_when_tau_is_undefined(tmp_path: Path) -> None:
    """One organ sampled means no specificity to show, so no panel is drawn.

    D42: where there is nothing to plot the page says so; it never renders an
    empty frame.

    The one organ this fixture samples is the **cardiac** one, and it is
    placed. So a build that dropped this branch would draw a real panel
    carrying `chart-cardiac`, not an empty frame that a bare "is there a
    figure" check could miss.

    **The sentence the dataset block opens with is asserted too**, because
    this fixture is also the smallest case of D42's other refusal and the two
    share a page. One heart stage, placed, and no curated phase vocabulary:
    no trajectory can be drawn, and the reason is the missing band, not a
    missing measurement. Until 2026-08-20 the page said "no measurement there
    is placed against this dataset's percentile grid" of the one stage whose
    percentile it printed in the table below, and said it of "1 stages".
    """
    profile = _expression_profile(
        (
            _dataset_profile_entry(
                stages=(
                    _stage_entry(
                        stage="4wpc",
                        specificity=None,
                        specificity_unavailable_reason="one_organ_sampled",
                        tissues=(
                            _tissue_entry(tissue="heart", unit="tpm", placement=_placement()),
                        ),
                    ),
                )
            ),
        )
    )
    section = _expression_section_text(_expression_page(tmp_path, {GATA4: profile}, _HEART_DATASET))

    assert "chart-cardiac" not in section
    assert 'class="sparks"' not in section, "an empty frame is what D42 forbids"
    assert "only one organ was sampled" in section

    lede = _dataset_lede(section)
    assert "no measurement there is placed" not in lede
    assert "detection floor at every" not in lede
    assert "1 stages sampled" not in lede, "a count with no singular form"
    assert "1 stage sampled" in lede


def test_the_spark_caption_names_both_ways_a_panel_line_breaks(tmp_path: Path) -> None:
    """One of the two causes was named; the other is the commoner one.

    The caption said a line breaks because "a stage below the floor" is not
    drawn. A line also breaks where this dataset has **no row for that organ
    at that stage**, and in a spark there is no tick, so the two are
    pixel-for-pixel identical. Measured 2026-08-20 on the built corpus: **87
    of the 92 genes have at least one panel split by a no-row stage** --
    testis on 82 pages, heart on 72, forebrain on 69 -- and the caption
    attributed every one of those gaps to the detection floor.

    Both clauses conditional on the cause being in *this* gene's panels, and
    four fixtures so no clause is measured only in the state where it fires.

    **A fifth for the reason axis**, which `_line_breaks` could not see: it
    reports *that* a placed run was interrupted, and the caption read the
    first half of that pair as "below the floor". A panel split by a stage
    with no percentile grid behind it was captioned as a low reading, the
    same conflation `_percentile_cell` refuses in words and `_band_caption`
    made on the picture above.
    """

    def caption(
        *organs: tuple[str, tuple[float | None | EllipsisType, ...]],
        reason: str = "below_detection_floor",
    ) -> str:
        return _expression_section_text(
            _expression_page(
                tmp_path,
                {
                    GATA4: _organ_panel(
                        *organs,
                        stages=("4wpc", "5wpc", "6wpc", "7wpc"),
                        unplaced_reason=reason,
                    )
                },
                _HEART_DATASET,
            )
        )

    floor = "below this dataset's detection floor</strong>, which this atlas does not place"
    no_row = "no row for that organ</strong> and nothing was measured to place"

    below = caption(("heart", (50.0, None, 150.0, 275.0)), ("liver", (0.5, 1.2, 2.0, 2.4)))
    assert floor in below
    assert no_row not in below

    missing = caption(("heart", (50.0, ..., 150.0, 275.0)), ("liver", (0.5, 1.2, 2.0, 2.4)))
    assert no_row in missing
    assert floor not in missing

    both = caption(("heart", (50.0, None, 150.0, 275.0)), ("liver", (0.5, ..., 2.0, 2.4)))
    assert floor in both
    assert no_row in both

    unbroken = caption(("heart", (50.0, 90.0, 150.0, 275.0)), ("liver", (0.5, 1.2, 2.0, 2.4)))
    assert floor not in unbroken
    assert no_row not in unbroken

    # The same broken shape, broken for a reason that is not the floor.
    grid = _spark_caption_text(
        caption(
            ("heart", (50.0, None, 150.0, 275.0)),
            ("liver", (0.5, 1.2, 2.0, 2.4)),
            reason="no_quantile_grid",
        )
    )
    assert floor not in grid
    assert "no complete percentile grid is published for this organ at this stage" in grid
    # The sentence the two clauses hang off stays regardless: what is drawn
    # is what this atlas placed, whether or not anything is missing.
    assert "Only medians this dataset placed against its percentile grid are drawn" in unbroken


def test_the_spark_caption_opens_with_a_claim_true_of_pages_that_drop_a_panel(
    tmp_path: Path,
) -> None:
    """ "One panel per organ this dataset sampled for this gene" was false on
    exactly the pages that then said so three sentences later.

    Measured 2026-08-20 on the built corpus: the opening clause is on 91
    pages and **25 of them also say "No panel is drawn for ..."**. The page
    corrected itself, which is not the same as being right: the first
    sentence is what a reader takes the grid to mean, and an organ dropped
    for having nothing placed is exactly the organ they would otherwise
    conclude was never looked at.

    Asserted as a **negative on the built page** as well as a positive, which
    is the assertion that fails when a deleted sentence comes back (CLAUDE.md
    section 4.35), and against literals rather than the constant that
    produced them (section 4.38).

    **Absorbs `test_an_organ_with_nothing_placed_is_named_rather_than_
    silently_dropped`**, measured 2026-08-20 not to be the unique killer of
    anything: dropping the `undrawn` list entirely is caught here and by
    `test_both_sentences_that_list_every_cardiac_organ_list_them_in_one_order`
    as well. Its record is the evidence-loss argument and moves with it --
    an organ sampled at every stage and placed at none has no line to draw
    and D42 forbids an empty frame for it, so dropping it without a word
    would leave the grid short of panels beside a tau sentence counting
    every organ sampled: the page contradicting itself, and a reader
    concluding this dataset never looked there.
    """
    opening = "One panel per organ with a placed measurement for this gene"
    false_opening = "One panel per organ this dataset sampled for this gene"

    complete = _expression_section_text(
        _expression_page(tmp_path, {GATA4: _SPARK_PANEL}, _HEART_DATASET)
    )
    dropped = _expression_section_text(
        _expression_page(
            tmp_path,
            {
                GATA4: _organ_panel(
                    ("heart", (50.0, 150.0, 275.0)),
                    ("kidney", (0.8, 1.0, 1.5)),
                    # Sampled at every stage and placed at none: no panel.
                    ("liver", (None, None, None)),
                )
            },
            _HEART_DATASET,
        )
    )

    assert "No panel is drawn for" not in complete
    assert "No panel is drawn for liver" in dropped
    assert opening in complete and opening in dropped
    assert false_opening not in complete and false_opening not in dropped
    # And the grid really did drop the panel, so the opening is being read
    # against a page that has one fewer panel than the dataset has organs.
    assert set(_spark_panels(complete)) == {"heart", "kidney", "liver"}
    assert set(_spark_panels(dropped)) == {"heart", "kidney"}


def test_a_spark_line_never_spans_a_stage_the_atlas_did_not_place(tmp_path: Path) -> None:
    """The same rule the trajectory keeps, in the smaller picture beside it.

    Measured 2026-08-20 on the committed corpus: 20 of the 85 charted genes
    carry a below-floor stage strictly between two plotted points. A single
    line through the placed points runs smooth across exactly the stages this
    dataset measured below its own floor -- the picture asserting a continuity
    the figures deny. A guard added to one layer is not a guard.

    Liver is placed at every stage in the same fixture, so "split correctly"
    and "never split anything" are distinguishable.
    """
    profile = _organ_panel(
        ("heart", (50.0, None, 150.0, 275.0)),
        ("liver", (0.5, 0.8, 1.2, 2.0)),
        stages=("4wpc", "5wpc", "6wpc", "7wpc"),
    )
    section = _expression_section_text(_expression_page(tmp_path, {GATA4: profile}, _HEART_DATASET))
    panels = _spark_panels(section)

    assert re.findall(r'<polyline class="chart-cardiac"', panels["heart"]) == [
        '<polyline class="chart-cardiac"'
    ], "the heart run of two adjacent stages is one line, and the lone point is not on it"
    assert '<circle class="chart-cardiac"' in panels["heart"]
    assert len(re.findall(r"<polyline", panels["liver"])) == 1
    assert "<circle" not in panels["liver"]


# --- The order every cardiac organ is named in (CLAUDE.md section 4.36) ------
#
# `_dataset_block` builds `cardiac` as a `frozenset`, and three places walk it:
# the trajectory charts, `_sampled_cardiac`'s sentence and `_spark_caption`'s
# "this atlas declares" clause. All three sort it, and until 2026-08-20 nothing
# could tell. `Dataset.cardiac_tissues` has exactly one member on the committed
# corpus and every fixture above declares one, so a one-element set iterates in
# one order under every seed. Two measurements: replacing all three
# `sorted(cardiac)` calls with a bare `cardiac` leaves the build byte-identical
# under PYTHONHASHSEED 0 and 12345, and mutating them one at a time was killed
# by 0 of the 1,050 tests that existed at 38fb290, the commit before this
# section.
#
# **Seven organs, and the number is measured rather than chosen.** A `frozenset`
# of *two* names iterates in sorted order about half the time, so a two-organ
# fixture would catch a dropped sort on roughly half of CI's random seeds --
# which is not a guard (section 4.12). Measured 2026-08-20 by constructing
# `frozenset(names)` under 400 explicit `PYTHONHASHSEED` values and counting how
# often it iterates pre-sorted:
#
#     2 names (atrium, ventricle)                 197/400   49.2%
#     3 names (aorta, heart, ventricle)            53/400   13.2%
#     4 names (+ atrium)                           17/400    4.2%
#     5 names (+ myocardium, endocardium)           6/400    1.5%
#     6 names (+ septum)                            2/400    0.5%
#     7 names (+ epicardium)   <- these seven       0/400    0.0%
#
# Widening the fixture is the cheap half of section 4.42's remedy and the whole
# of it here, because the rate reaches zero at a size a real dataset could
# plausibly declare. Re-measure rather than cite this if the names change: the
# rate is a property of the exact strings, not of the count.
_CARDIAC_ORGANS: tuple[str, ...] = (
    "ventricle",
    "aorta",
    "septum",
    "myocardium",
    "endocardium",
    "atrium",
    "epicardium",
)

# Declared above in a deliberately unsorted order, so a renderer walking
# `Dataset.cardiac_tissues` itself rather than the sorted frozenset fails too.
_CARDIAC_ORGANS_SORTED: tuple[str, ...] = (
    "aorta",
    "atrium",
    "endocardium",
    "epicardium",
    "myocardium",
    "septum",
    "ventricle",
)

_MANY_CARDIAC_DATASET = {"E-MTAB-6814": _profile_dataset(cardiac_tissues=_CARDIAC_ORGANS)}


def _trajectory_organs(section: str) -> list[str]:
    """Every trajectory's organ, in the order the section draws them.

    Matched on the trajectory's own 560x150 `viewBox` rather than on `<svg`,
    for the reason `_trajectory_figure` is (CLAUDE.md section 4.19): the organ
    small multiples are 96x40 `<svg>`s in the same section, and a *cardiac*
    organ's panel opens `Median abundance in whole <organ>,` exactly as its
    trajectory does -- identical through the comma this pattern stops at. An
    unscoped search would return every cardiac organ twice, interleaved in an
    order the trajectories do not decide.
    """
    return re.findall(
        r'<svg class="chart" viewBox="0 0 560 150"[^>]*><title>'
        r"Median abundance in whole ([^,]+),",
        section,
    )


# Four phases whose id order, label order and walk order are three different
# sequences, so a band list can be attributed to exactly one of them. `covered`
# in `_banded_phases` is a **dict**, filled by walking the stage list, so its own
# key order is insertion order and carries no `PYTHONHASHSEED` risk at all --
# `_banded_phases`' docstring records that measurement and this fixture does not
# repeat it. What the sort buys is the property that docstring claims and
# nothing tested: the band order is a function of the covered ids alone, so it
# does not move when a dataset's stage list does.
_BAND_LABELS: Final = {"alpha": "Zulu", "bravo": "Yankee", "charlie": "Xray", "delta": "Whisky"}

# Deliberately neither sorted nor reverse-sorted: a walk in this order gives a
# fourth sequence again, so "sorted by id" is the only rule that produces the
# literal asserted below.
_BAND_WALK: Final = ("delta", "alpha", "charlie", "bravo")


def _phase_walk(*phase_ids: str) -> ExpressionProfile:
    """One placed heart stage per phase, in the order given."""
    return _expression_profile(
        (
            _dataset_profile_entry(
                stages=tuple(
                    _stage_entry(
                        stage=f"{index + 4}wpc",
                        phase=_phase_info(phase_ids=(phase_id,)),
                        tissues=(
                            _tissue_entry(
                                tissue="heart",
                                median=5.0 * (index + 1),
                                unit="tpm",
                                placement=_placement(),
                            ),
                        ),
                    )
                    for index, phase_id in enumerate(phase_ids)
                )
            ),
        )
    )


def test_the_phase_bands_are_ordered_by_id_not_by_the_stage_walk_that_found_them(
    tmp_path: Path,
) -> None:
    """`_banded_phases` sorts `covered`, and removing that sort changed 85 gene
    pages and the manifest while every test stayed green.

    Not a determinism guard, and the docstring beside the sort says so after
    being measured: `covered` is a dict filled by a deterministic walk, so its
    key order is reproducible with or without the sort. What the sort buys is
    that the band order is a function of the covered phase ids alone rather
    than of the stage list that found them -- a claim that was written down,
    measured once by hand, and pinned by nothing.

    Both halves are asserted, because either alone is weak. The literal alone
    cannot tell "sorted by id" from "this fixture's walk happens to agree";
    the two stage orders alone cannot tell "sorted" from "sorted by label", or
    from any other rule stable under reversal. Together only sorting by id
    produces both results.
    """
    vocabulary = _cardiac_phases(
        *(
            _cardiac_phase(phase_id=phase_id, label=label)
            for phase_id, label in _BAND_LABELS.items()
        )
    )

    def bands(*phase_ids: str) -> list[str]:
        section = _expression_section_text(
            _expression_page(
                tmp_path,
                {GATA4: _phase_walk(*phase_ids)},
                _HEART_DATASET,
                phases=vocabulary,
            )
        )
        return re.findall(r'<rect class="chart-band"[^>]*><title>([^<]*)</title>', section)

    # Sorted by id -- which is neither the walk order (Whisky, Zulu, Xray,
    # Yankee) nor the label order (Whisky, Xray, Yankee, Zulu).
    assert bands(*_BAND_WALK) == ["Zulu", "Yankee", "Xray", "Whisky"]
    # And it does not move when the dataset's stage list does.
    assert bands(*reversed(_BAND_WALK)) == ["Zulu", "Yankee", "Xray", "Whisky"]


def test_the_trajectories_are_drawn_in_one_order_whatever_the_frozenset_does(
    tmp_path: Path,
) -> None:
    """Seven cardiac organs, every one of them charted, pinned to a literal.

    Against a literal rather than `sorted(_CARDIAC_ORGANS)`, for the reason
    `test_build_landing.py` pins published wording against literals (section
    4.38): a test that computes its own expectation from the same expression
    the code uses compares the module to itself and passes for any ordering
    both agree on.
    """
    profile = _organ_panel(
        *((organ, (5.0, 40.0, 275.0)) for organ in _CARDIAC_ORGANS),
        highest_in="ventricle",
        phase=_phase_info(phase_ids=("atrial_septum_morphogenesis",)),
    )
    section = _expression_section_text(
        _expression_page(
            tmp_path, {GATA4: profile}, _MANY_CARDIAC_DATASET, phases=_cardiac_phases()
        )
    )

    assert _trajectory_organs(section) == [
        "aorta",
        "atrium",
        "endocardium",
        "epicardium",
        "myocardium",
        "septum",
        "ventricle",
    ]


def test_both_sentences_that_list_every_cardiac_organ_list_them_in_one_order(
    tmp_path: Path,
) -> None:
    """The two orders a page states in words rather than in geometry.

    One fixture, because the two sentences are published together and only
    together: `_no_trajectory_sentence` renders exactly when no cardiac organ
    has a trajectory, and `_spark_caption`'s "this atlas declares" clause
    exactly when no cardiac organ has a panel either. Every cardiac organ here
    is sampled and below the floor; liver is placed, so the small multiples
    still draw and the caption still has something to caption.

    Both are asserted as one substring apiece rather than as a set of names,
    because the defect is the *order* and a membership check cannot see it.
    """
    profile = _organ_panel(
        *((organ, (None, None, None)) for organ in _CARDIAC_ORGANS),
        ("liver", (0.5, 1.2, 2.0)),
        highest_in="liver",
    )
    section = _expression_section_text(
        _expression_page(tmp_path, {GATA4: profile}, _MANY_CARDIAC_DATASET)
    )

    assert (
        "whole aorta (3 stages sampled), whole atrium (3 stages sampled), "
        "whole endocardium (3 stages sampled), whole epicardium (3 stages sampled), "
        "whole myocardium (3 stages sampled), whole septum (3 stages sampled), "
        "whole ventricle (3 stages sampled)"
    ) in section
    assert (
        "this atlas declares aorta, atrium, endocardium, epicardium, myocardium, "
        "septum, ventricle this dataset's cardiac tissues"
    ) in section
    assert sorted(_spark_panels(section)) == ["liver"], (
        "a cardiac organ drew a panel, so the caption's other branch was taken"
    )

    # The third list in the same caption, from `_sampled_tissues`' own sort
    # over a `set` of the organs the payload names -- a different frozen
    # collection with the same failure mode, and the same eight-name fixture
    # measures 0/400 pre-sorted.
    assert (
        "No panel is drawn for aorta, atrium, endocardium, epicardium, myocardium, "
        "septum, ventricle:"
    ) in section
