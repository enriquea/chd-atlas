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
from chd_atlas.build.runner import build_site
from chd_atlas.build.validity import GeneValidity, uncurated
from chd_atlas.corpus import Corpus
from chd_atlas.models.assertion import Evidence, LesionAssertion, SupplementaryLocator

TBX5 = "HGNC:11604"
GATA4 = "HGNC:4173"

# tests/unit/test_build_landing.py -> tests/unit -> tests -> repo root, the same
# climb `tests/unit/test_build_runner.py` uses to find the committed corpus.
REPO = Path(__file__).parent.parent.parent

_MIRRORED_ROW_LABEL = "Genes in scope with mirrored ClinGen/GenCC validity (browsable once curated)"


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
    corpus: Corpus, symbols: dict[str, str], validity: dict[str, GeneValidity], tmp_path: Path
) -> str:
    emitter = Emitter(root=tmp_path)
    build_landing(corpus, symbols=symbols, validity=validity, emitter=emitter)
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

    build_landing(corpus, symbols={TBX5: "TBX5"}, validity={TBX5: uncurated()}, emitter=emitter)

    assert list(emitter.checksums) == ["index.html"]
    assert (tmp_path / "index.html").is_file()


def test_the_published_and_mirrored_counts_are_derived_not_hardcoded(tmp_path: Path) -> None:
    """Two different inputs must publish two different pages, for both gene counts.

    A landing page whose numbers were hardcoded strings would pass a test built
    against one fixture and only that one; varying the fixture and checking the
    rendered counts move with it is what tells a literal from a derivation. Two
    assertions on two different genes, and five mirrored genes -- none of which
    match the committed corpus's one assertion on one gene, or its 154 mirrored
    genes -- so a hardcoded "1"/"1"/"154" would fail here even though nothing
    here matches those numbers.
    """
    corpus = Corpus(
        root=Path("."),
        assertions=(_assertion(), _assertion(id="CHDA:AST:0000002", gene=GATA4)),
    )
    validity = {gene: uncurated() for gene in ("HGNC:1", "HGNC:2", "HGNC:3", "HGNC:4", "HGNC:5")}

    text = _build(corpus, {TBX5: "TBX5", GATA4: "GATA4"}, validity, tmp_path)

    assert re.search(r"<dt>Curated gene-disease assertions</dt>\s*<dd>2</dd>", text)
    assert re.search(r"<dt>Genes published</dt>\s*<dd>2</dd>", text)
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
    build actually publishes. `genes/index.json` is keyed on genes carrying at
    least one curated assertion (`derive.gene_facts`), which for the committed
    corpus is one gene, not 154 — checked here against the real built file
    rather than against a literal, so this fails the same way the original bug
    would have if it recurred.
    """
    text = (real_build / "index.html").read_text(encoding="utf-8")
    genes_index = json.loads((real_build / "genes" / "index.json").read_text(encoding="utf-8"))
    published = len(genes_index["genes"])

    match = re.search(r"<dt>Genes published</dt>\s*<dd>(\d+)</dd>", text)
    assert match, "the page no longer publishes a 'Genes published' row"
    assert int(match.group(1)) == published


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


def test_the_research_use_statement_is_exactly_the_two_sentences_specified(
    tmp_path: Path,
) -> None:
    """The development-status notice: exactly two sentences, nothing else.

    The paragraph this notice used to carry beneath them -- "Today it holds {N}
    curated gene-disease assertion..." -- editorialised with a claim false for
    100% of described genes ("almost every gene this site describes carries an
    upstream expert panel's or submitter's classification and nothing this
    atlas has independently curated"): the committed corpus describes one gene,
    and it is the curated one. That paragraph is gone; this test pins the
    notice section down to the two sentences that replace it, so it cannot grow
    a new editorial claim unnoticed.
    """
    corpus = Corpus(root=Path("."), assertions=(_assertion(),))

    text = _build(corpus, {TBX5: "TBX5"}, {TBX5: uncurated()}, tmp_path)
    notice = _prose(_section(text, "<h2>Development status &amp; research use</h2>"))

    expected = (
        "This atlas is under active development. It is "
        "<strong>not a clinical decision-support tool</strong> and must not be used to make "
        "or guide a diagnostic, treatment or any other clinical decision."
    )
    # `_section` stops at the next `<h2>`, which for the notice is `</section>` first
    # -- strip everything from the closing tag on, so this compares the paragraph
    # only, not the section wrapper around it.
    notice = notice[: notice.index("</section>")]
    assert expected in notice
    assert "Today it holds" not in notice
    assert "independently curated" not in notice


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
