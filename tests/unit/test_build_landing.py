# tests/unit/test_build_landing.py
"""`index.html`, the site's front door: the research-use statement and its counts."""

from __future__ import annotations

import json
import re
import shutil
from datetime import date
from pathlib import Path

import pytest

from chd_atlas.build.emit import Emitter
from chd_atlas.build.landing import build_landing
from chd_atlas.build.paths import LANDING
from chd_atlas.build.render import RESEARCH_USE_NOTICE, STYLESHEET
from chd_atlas.build.runner import build_site
from chd_atlas.build.validity import GeneValidity, uncurated
from chd_atlas.corpus import Corpus
from chd_atlas.models.assertion import Evidence, LesionAssertion, SupplementaryLocator

TBX5 = "HGNC:11604"
GATA4 = "HGNC:4173"

# tests/unit/test_build_landing.py -> tests/unit -> tests -> repo root, the same
# climb `tests/unit/test_build_runner.py` uses to find the committed corpus.
REPO = Path(__file__).parent.parent.parent

# The parenthetical said "browsable once curated" until D21, which made it
# false: 23 of these 154 genes are browsable today and 22 of those carry no
# curation from this atlas at all. What decides is ClinGen's grade, and that is
# what it now says. The row's placement and label are pinned below for the
# reason they always were -- this number must never read as published coverage.
#
# It then said "(browsable once ClinGen grades it definitive)", which was false
# in the other direction: measured 2026-08-04 against the committed mirrors, 20
# of these 154 genes already carry a ClinGen `Definitive` grade and are still
# not browsable, because the grade names a disease outside CHD scope -- ELN was
# graded Definitive in 2024, for cutis laxa. The caption told a reader to wait
# for something that had already happened. The gate is definitive **for a
# disease in that scope**, which is what `published_genes` actually tests.
#
# A literal here rather than an import of `landing._MIRRORED_ROW_LABEL`: a test
# that reads the label out of the module and compares it to itself passes on any
# wording at all, including the one this replaces.
_MIRRORED_ROW_LABEL = (
    "Genes with mirrored validity in CHD scope "
    "(browsable once ClinGen grades it definitive for a disease in that scope)"
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
        "gene": TBX5,
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


def _build(
    corpus: Corpus,
    symbols: dict[str, str],
    validity: dict[str, GeneValidity],
    tmp_path: Path,
    published: set[str] | None = None,
) -> str:
    """Render the page.

    `published` defaults to the genes the fixture asserts, which is what the
    page said before D21 and is *not* what it says now -- the default only
    keeps the tests that are about something else (escaping, pluralisation, the
    research-use notice) from having to name a population they do not care
    about. Every test that is about the published-gene figure passes it
    explicitly, and one of them passes a set that disagrees with the assertions
    on purpose.
    """
    emitter = Emitter(root=tmp_path)
    build_landing(
        corpus,
        symbols=symbols,
        validity=validity,
        published=(
            {assertion.gene for assertion in corpus.assertions} if published is None else published
        ),
        emitter=emitter,
    )
    return (tmp_path / LANDING).read_text(encoding="utf-8")


def _prose(text: str) -> str:
    """Whitespace-collapsed text, for asserting a sentence rather than a source line.

    The template wraps long sentences across several source lines the way every
    other line in this file does, and a browser collapses that whitespace when
    it renders — but `"in text"` does not, so a phrase split across a wrap would
    fail a substring check that has nothing to do with what a reader sees. Used
    for prose; the escaping and no-external-request tests below check the raw
    text instead, because a literal `<script` or an unescaped `<` is exactly
    what must not be normalised away.
    """
    return re.sub(r"\s+", " ", text)


def _section(text: str, heading: str) -> str:
    """The slice of the page from one `<h2>` heading up to (not including) the next.

    Lets a test pin *where* a row lives, not only whether the text appears
    somewhere on the page — the defect this file is guarding against was a true
    number in the wrong section, which a bare substring check cannot catch.
    """
    start = text.index(heading)
    nxt = text.find("<h2>", start + len(heading))
    return text[start : nxt if nxt != -1 else len(text)]


@pytest.fixture(scope="module")
def real_build(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """A real build of the committed corpus, for checking the page against `genes/index.json`.

    `curation` and `mirrors` are copied (cheap, and nothing here mutates them);
    `ontologies` is symlinked, the same shortcut `test_build_runner.py` and
    `test_docs_match_the_build.py` use for the same 11 MB of pinned OBO releases.
    Module-scoped because `build_site` costs real time (validation plus every
    builder) and every test in this file that needs a real build asks the same
    question of the same output.
    """
    base = tmp_path_factory.mktemp("landing-real-build")
    root = base / "repo"
    root.mkdir()
    for name in ("curation", "mirrors"):
        shutil.copytree(REPO / name, root / name)
    (root / "ontologies").symlink_to(REPO / "ontologies")
    out = base / "dist"
    build_site(root, out)
    return out


def test_the_page_is_published_through_write_text_and_reaches_the_checksums(
    tmp_path: Path,
) -> None:
    """`build_landing` writes exactly one file, at the path `LANDING` names.

    Through `emitter.checksums` rather than only the file on disk, because that
    mapping is what `write_manifest` publishes as `files` — a page written by a
    method that bypassed `_write` would still be on disk and invisible to the
    manifest, which is the exact failure `Emitter.write_text` exists to close.
    """
    emitter = Emitter(root=tmp_path)
    corpus = Corpus(root=Path("."), assertions=(_assertion(),))

    build_landing(
        corpus,
        symbols={TBX5: "TBX5"},
        validity={TBX5: uncurated()},
        published={TBX5},
        emitter=emitter,
    )

    assert list(emitter.checksums) == ["index.html"]
    assert (tmp_path / "index.html").is_file()


def test_the_published_and_mirrored_counts_are_derived_not_hardcoded(tmp_path: Path) -> None:
    """Three different inputs must publish three different numbers.

    A landing page whose numbers were hardcoded strings would pass a test built
    against one fixture and only that one; varying the fixture and checking the
    rendered counts move with it is what tells a literal from a derivation. Two
    assertions on two genes, three published genes and five mirrored ones --
    none of which match the committed corpus's one assertion, its 23 published
    genes or its 154 mirrored ones -- so a hardcoded "1"/"23"/"154" fails here.

    The three fixture values are deliberately all different, and `published` is
    deliberately not the asserted genes. Since D21 those are separate
    populations -- 22 of the 23 genes the site publishes carry no assertion --
    and the figure wired to the wrong one of the three is the mistake this
    catches: the page used to count `{assertion.gene for ...}` under this very
    label, which would render 2 here instead of 3.
    """
    corpus = Corpus(
        root=Path("."),
        assertions=(_assertion(), _assertion(id="CHDA:AST:0000002", gene=GATA4)),
    )
    validity = {gene: uncurated() for gene in ("HGNC:1", "HGNC:2", "HGNC:3", "HGNC:4", "HGNC:5")}

    text = _build(
        corpus,
        {TBX5: "TBX5", GATA4: "GATA4"},
        validity,
        tmp_path,
        published={TBX5, GATA4, "HGNC:1"},
    )

    assert re.search(r"<dt>Curated gene-disease assertions</dt>\s*<dd>2</dd>", text)
    assert re.search(r"<dt>Genes published</dt>\s*<dd>3</dd>", text)
    assert re.search(
        rf"<dt>{re.escape(_MIRRORED_ROW_LABEL)}</dt>\s*<dd>5</dd>",
        text,
    )
    assert "154" not in text


def test_the_published_gene_count_agrees_with_a_real_build_of_genes_index_json(
    real_build: Path,
) -> None:
    """The page's "Genes published" figure must equal `len(genes/index.json["genes"])`.

    This is the regression the false "154 genes" claim was: `mirrored_gene_count
    = len(validity)` counted every gene either mirror curates, not the genes the
    build actually publishes. `genes/index.json` is keyed on `published`
    (`derive.gene_facts`), which for the committed corpus is 23 genes — not 154,
    and no longer the 1 gene that carries an assertion either. Checked here
    against the real built file rather than against a literal, so it fails the
    same way the original bug would if it recurred, and so that widening the
    gate again moves both numbers together or fails.
    """
    text = (real_build / "index.html").read_text(encoding="utf-8")
    genes_index = json.loads((real_build / "genes" / "index.json").read_text(encoding="utf-8"))
    published = len(genes_index["genes"])

    match = re.search(r"<dt>Genes published</dt>\s*<dd>(\d+)</dd>", text)
    assert match, "the page no longer publishes a 'Genes published' row"
    assert int(match.group(1)) == published


def test_the_landing_page_distinguishes_published_from_curated(tmp_path: Path) -> None:
    """The promotion review forced this distinction into prose once already.

    23 genes are published on an expert panel's classification and 1 carries the
    atlas's own evidence. A single "genes" count would read as coverage the site
    does not have.

    The fixture separates all three numbers that could be wired to this row --
    two assertions, on one gene, against three published genes -- because the
    two plausible miswirings both render a number that is *true of something*:
    `len(corpus.assertions)` gives 2 and `len(published)` gives 3, and either
    would pass a fixture where they happened to coincide. Only 1 is the genes
    the atlas has curated.

    The sentence naming those genes is checked here too: it and the row are one
    derivation in `_render`, and a fixture with more assertions than genes is
    where a second derivation would show up as "2 ... so far: TBX5" beside a
    row saying 2.
    """
    corpus = Corpus(
        root=Path("."),
        assertions=(_assertion(), _assertion(id="CHDA:AST:0000002")),
    )

    text = _build(
        corpus,
        {TBX5: "TBX5"},
        {TBX5: uncurated()},
        tmp_path,
        published={TBX5, GATA4, "HGNC:1"},
    )
    published_section = _section(text, "<h2>What's published</h2>")

    assert re.search(r"<dt>Curated gene-disease assertions</dt>\s*<dd>2</dd>", published_section)
    assert re.search(r"<dt>Genes published</dt>\s*<dd>3</dd>", published_section)
    assert re.search(r"<dt>Genes the atlas has curated</dt>\s*<dd>1</dd>", published_section)
    assert "2 curated gene-disease assertions so far: TBX5 (HGNC:11604)." in _prose(text)


def test_the_page_shares_the_shell_and_the_stylesheet_every_other_page_uses(
    tmp_path: Path,
) -> None:
    """The landing page carries no `<style>` of its own.

    It had one until the gene pages arrived, and three page kinds each carrying
    their own copy is three stylesheets that drift. `render.STYLESHEET` is now
    the only one, so this asserts the rendered page contains it verbatim rather
    than merely that some `<style>` block exists -- a landing page that grew a
    second block back would still satisfy the weaker check.

    The `<nav>` is the other half of the shell, and `root=""` is what makes its
    links right from the site root: a landing page rendered with `"../"` would
    point every visitor's "Genes" link above the deploy root.
    """
    corpus = Corpus(root=Path("."), assertions=(_assertion(),))

    text = _build(corpus, {TBX5: "TBX5"}, {TBX5: uncurated()}, tmp_path)

    assert f"<style>{STYLESHEET}</style>" in text
    assert text.count("<style>") == 1
    assert '<a href="index.html">CHD Atlas</a>' in text
    assert '<a href="genes/index.html">Genes</a>' in text


def test_the_front_page_links_to_the_browsable_gene_index(tmp_path: Path) -> None:
    """ "Browse the data" must offer the page, not only the payload behind it.

    A reader who opens `index.html` is not going to read `genes/index.json`.
    Pinned inside the section rather than anywhere on the page, because the
    `<nav>` links to the same file and would satisfy a bare substring check
    while the browse list stayed JSON-only.
    """
    corpus = Corpus(root=Path("."), assertions=(_assertion(),))

    text = _build(corpus, {TBX5: "TBX5"}, {TBX5: uncurated()}, tmp_path)
    browse = _prose(_section(text, "<h2>Browse the data</h2>"))

    assert 'href="genes/index.html"' in browse
    assert "every gene the atlas publishes, with filters." in browse


def test_the_relabelled_mirrored_validity_row_sits_outside_whats_published(
    tmp_path: Path,
) -> None:
    """The mirrored-validity count must not read as published coverage.

    Pinned by section as well as by label: a future edit that moves the row
    back under "What's published" without touching its label would otherwise
    pass every other test in this file, and that is exactly the shape the
    original defect had -- a true number, correctly derived, in a section that
    implied something false about it.
    """
    corpus = Corpus(root=Path("."), assertions=(_assertion(),))
    validity = {gene: uncurated() for gene in ("HGNC:1", "HGNC:2", "HGNC:3", "HGNC:4", "HGNC:5")}

    text = _build(corpus, {TBX5: "TBX5"}, validity, tmp_path)

    published_section = _section(text, "<h2>What's published</h2>")
    sources_section = _section(text, "<h2>Where this data comes from</h2>")

    assert _MIRRORED_ROW_LABEL not in published_section
    assert "mirrored" not in published_section.lower()
    assert re.search(
        rf"<dt>{re.escape(_MIRRORED_ROW_LABEL)}</dt>\s*<dd>5</dd>",
        sources_section,
    )


def test_a_single_assertion_is_singular_not_plural(tmp_path: Path) -> None:
    """The wording matters as much as the number: "1 assertions" reads as unreviewed."""
    corpus = Corpus(root=Path("."), assertions=(_assertion(),))

    text = _build(corpus, {TBX5: "TBX5"}, {TBX5: uncurated()}, tmp_path)

    assert "1 curated gene-disease assertion so far" in text
    assert "1 curated gene-disease assertions" not in text


def test_no_assertions_yet_reads_as_none_rather_than_as_an_empty_list(tmp_path: Path) -> None:
    """A corpus with nothing curated still has to publish a page, and an honest one."""
    corpus = Corpus(root=Path("."), assertions=())

    text = _build(corpus, {}, {}, tmp_path)

    assert "0 curated gene-disease assertions so far" in text
    assert "none yet" in text
    assert re.search(r"<dt>Genes published</dt>\s*<dd>0</dd>", text)


def test_the_research_use_statement_reaches_the_front_page_exactly_once(
    tmp_path: Path,
) -> None:
    """The development-status notice: exactly two sentences, once, and no more.

    The paragraph this notice used to carry beneath them -- "Today it holds {N}
    curated gene-disease assertion..." -- editorialised with a claim false for
    100% of described genes ("almost every gene this site describes carries an
    upstream expert panel's or submitter's classification and nothing this
    atlas has independently curated"): the committed corpus describes one gene,
    and it is the curated one. That paragraph is gone, and this page has no
    editorial prose left between the notice and "What this is".

    The notice itself now belongs to `render.document`, which puts it on all 25
    pages instead of this one (see `test_build_render.py`). What is still this
    file's business is that the landing page did not keep a second copy: the
    section that used to hold it was removed here rather than left to render
    alongside the shared one, and two copies of the site's only research-use
    statement is how one of them comes to be edited and the other not.

    Counted on `RESEARCH_USE_NOTICE` rather than on the substring "not a
    clinical decision-support tool", because the emphasis markup is part of what
    the owner specified and a page carrying the sentence with the `<strong>`
    dropped would satisfy the looser check.
    """
    corpus = Corpus(root=Path("."), assertions=(_assertion(),))

    text = _build(corpus, {TBX5: "TBX5"}, {TBX5: uncurated()}, tmp_path)

    assert text.count(RESEARCH_USE_NOTICE) == 1
    assert "<h2>Development status &amp; research use</h2>" not in text
    assert "Today it holds" not in text
    assert "independently curated" not in text


def test_a_gene_symbol_is_escaped_before_it_reaches_the_page(tmp_path: Path) -> None:
    """A mirror cell is curated text, not markup — the same discipline `bundles.py`
    and `search.py` already apply to a symbol before it reaches a payload.

    `mirrors/genes.tsv` is reviewed for scientific accuracy, not for HTML safety,
    so nothing upstream of this page guarantees a symbol contains no `<`. Without
    `html.escape` this fixture's symbol would open a real `<script>` tag in the
    page every visitor's browser parses.
    """
    corpus = Corpus(root=Path("."), assertions=(_assertion(),))

    text = _build(corpus, {TBX5: "<script>alert(1)</script>"}, {TBX5: uncurated()}, tmp_path)

    assert "<script>alert(1)</script>" not in text
    assert "&lt;script&gt;" in text


def test_a_gene_absent_from_the_registry_falls_back_to_its_hgnc_id(tmp_path: Path) -> None:
    """The same registry fallback `build_genes` and `search.py` already apply."""
    corpus = Corpus(root=Path("."), assertions=(_assertion(),))

    text = _build(corpus, {}, {TBX5: uncurated()}, tmp_path)

    assert TBX5 in text


def test_the_page_makes_no_external_request(tmp_path: Path) -> None:
    """Self-contained: no script, no stylesheet or font pulled from a CDN.

    A `<a href="https://…">` is fine — that is a link a reader clicks, not a
    request the page issues on load — so this checks for the tags and
    attributes that *do* fetch automatically rather than for the substring
    "http" anywhere in the document.
    """
    corpus = Corpus(root=Path("."), assertions=(_assertion(),))

    text = _build(corpus, {TBX5: "TBX5"}, {TBX5: uncurated()}, tmp_path)

    assert "<script" not in text
    assert "<link" not in text
    assert "@import" not in text
