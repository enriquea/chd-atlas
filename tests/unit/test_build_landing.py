# tests/unit/test_build_landing.py
"""`index.html`, the site's front door: the research-use statement and its counts."""

from __future__ import annotations

import re
from datetime import date
from pathlib import Path

from chd_atlas.build.emit import Emitter
from chd_atlas.build.landing import build_landing
from chd_atlas.build.paths import LANDING
from chd_atlas.build.validity import GeneValidity, uncurated
from chd_atlas.corpus import Corpus
from chd_atlas.models.assertion import Evidence, LesionAssertion, SupplementaryLocator

TBX5 = "HGNC:11604"
GATA4 = "HGNC:4173"


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


def test_the_counts_on_the_page_are_read_from_the_corpus_and_validity_not_hardcoded(
    tmp_path: Path,
) -> None:
    """Two different inputs must publish two different pages.

    A landing page whose numbers were hardcoded strings would pass a test built
    against one fixture and only that one; varying the fixture and checking the
    rendered count moves with it is what tells a literal from a derivation. Two
    assertions on two different genes, and five mirrored genes -- none of which
    match the single-assertion, one-gene fixture used elsewhere in this file --
    so a hardcoded "1 curated gene-disease assertion" or "154 genes" would fail
    here even though nothing here matches those numbers.
    """
    corpus = Corpus(
        root=Path("."),
        assertions=(_assertion(), _assertion(id="CHDA:AST:0000002", gene=GATA4)),
    )
    validity = {gene: uncurated() for gene in ("HGNC:1", "HGNC:2", "HGNC:3", "HGNC:4", "HGNC:5")}

    text = _build(corpus, {TBX5: "TBX5", GATA4: "GATA4"}, validity, tmp_path)

    assert "2 curated gene-disease assertions" in text
    assert "5 genes" in text
    assert "1 curated gene-disease assertion " not in text
    assert "154 genes" not in text


def test_a_single_assertion_is_singular_not_plural(tmp_path: Path) -> None:
    """The wording matters as much as the number: "1 assertions" reads as unreviewed."""
    corpus = Corpus(root=Path("."), assertions=(_assertion(),))

    text = _build(corpus, {TBX5: "TBX5"}, {TBX5: uncurated()}, tmp_path)

    assert "1 curated gene-disease assertion " in text
    assert "1 curated gene-disease assertions" not in text


def test_no_assertions_yet_reads_as_none_rather_than_as_an_empty_list(tmp_path: Path) -> None:
    """A corpus with nothing curated still has to publish a page, and an honest one."""
    corpus = Corpus(root=Path("."), assertions=())

    text = _build(corpus, {}, {}, tmp_path)

    assert "0 curated gene-disease assertions" in text
    assert "none yet" in text


def test_the_research_use_statement_is_unambiguous_and_on_the_page(tmp_path: Path) -> None:
    """The one sentence this page exists to make impossible to miss."""
    corpus = Corpus(root=Path("."), assertions=(_assertion(),))

    prose = _prose(_build(corpus, {TBX5: "TBX5"}, {TBX5: uncurated()}, tmp_path))

    assert "not a clinical decision-support tool" in prose
    assert "must not be used to make or guide" in prose


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
