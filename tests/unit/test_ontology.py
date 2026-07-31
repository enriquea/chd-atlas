# tests/unit/test_ontology.py
import shutil
from pathlib import Path

import pytest

from chd_atlas.validate.ontology import OntologyRegistry, validate_terms

FIXTURES = Path(__file__).parent.parent / "fixtures"
FIXTURE = FIXTURES / "mini-hp.obo"
TYPEDEF_FIXTURE = FIXTURES / "mini-so.obo"


@pytest.fixture
def registry(tmp_path: Path) -> OntologyRegistry:
    obo = tmp_path / "ontologies" / "hp-2026-06-01.obo"
    obo.parent.mkdir(parents=True)
    shutil.copy(FIXTURE, obo)
    return OntologyRegistry.from_files(tmp_path, {"HP": "ontologies/hp-2026-06-01.obo"})


@pytest.fixture
def typedef_registry(tmp_path: Path) -> OntologyRegistry:
    obo = tmp_path / "ontologies" / "so-2026-06-01.obo"
    obo.parent.mkdir(parents=True)
    shutil.copy(TYPEDEF_FIXTURE, obo)
    return OntologyRegistry.from_files(tmp_path, {"SO": "ontologies/so-2026-06-01.obo"})


def test_resolves_a_term_present_in_the_pinned_release(registry: OntologyRegistry) -> None:
    assert validate_terms(["HP:0001631"], registry, "curation/x.yaml") == []


def test_reports_a_term_absent_from_the_pinned_release(registry: OntologyRegistry) -> None:
    issues = validate_terms(["HP:9999999"], registry, "curation/x.yaml")
    assert [i.code for i in issues] == ["ONT001"]
    assert "HP:9999999" in issues[0].message


def test_reports_an_obsolete_term(registry: OntologyRegistry) -> None:
    issues = validate_terms(["HP:0001636"], registry, "curation/x.yaml")
    assert [i.code for i in issues] == ["ONT002"]


def test_reports_a_prefix_with_no_pinned_ontology(registry: OntologyRegistry) -> None:
    issues = validate_terms(["MONDO:0007259"], registry, "curation/x.yaml")
    assert [i.code for i in issues] == ["ONT003"]
    assert "MONDO" in issues[0].message


def test_taxon_terms_resolve_against_the_allowlist(registry: OntologyRegistry) -> None:
    assert validate_terms(["NCBITaxon:10090"], registry, "curation/x.yaml") == []
    issues = validate_terms(["NCBITaxon:6239"], registry, "curation/x.yaml")
    assert [i.code for i in issues] == ["ONT001"]


def test_reports_every_bad_term_not_just_the_first(registry: OntologyRegistry) -> None:
    """Five of the other tests pass single-element lists; a first-error-only
    regression would keep them all green."""
    issues = validate_terms(
        ["HP:0001631", "HP:9999999", "HP:0001636", "MONDO:0007259", "HP:0001629"],
        registry,
        "curation/x.yaml",
    )

    assert [i.code for i in issues] == ["ONT001", "ONT002", "ONT003"]


def test_missing_pinned_file_is_reported(tmp_path: Path) -> None:
    registry = OntologyRegistry.from_files(tmp_path, {"HP": "ontologies/absent.obo"})
    assert [i.code for i in registry.load_issues] == ["ONT004"]


def test_unparseable_ontology_is_reported_not_raised(tmp_path: Path) -> None:
    """pronto raises ValueError, SyntaxError or OSError depending on the damage."""
    obo = tmp_path / "ontologies" / "broken.obo"
    obo.parent.mkdir(parents=True)
    obo.write_text("this is not an ontology at all\n")

    registry = OntologyRegistry.from_files(tmp_path, {"HP": "ontologies/broken.obo"})

    assert [i.code for i in registry.load_issues] == ["ONT004"]


def test_corrupt_stanza_is_reported_not_raised(tmp_path: Path) -> None:
    """A corrupt stanza raises SyntaxError, which is not a ValueError subclass."""
    obo = tmp_path / "ontologies" / "corrupt.obo"
    obo.parent.mkdir(parents=True)
    obo.write_text("format-version: 1.2\nontology: hp\n\n[Term]\nid HP:0001631\n")

    registry = OntologyRegistry.from_files(tmp_path, {"HP": "ontologies/corrupt.obo"})

    assert [i.code for i in registry.load_issues] == ["ONT004"]


def test_a_typedef_curie_is_reported_not_raised(typedef_registry: OntologyRegistry) -> None:
    """`curie in ontology` is true for a Typedef, and `get_term` then raises
    KeyError. The atlas only references terms, so a relationship id must fall
    through to ONT001 rather than escaping validate_terms."""
    issues = validate_terms(["SO:has_part"], typedef_registry, "curation/x.yaml")

    assert [i.code for i in issues] == ["ONT001"]
    assert "SO:has_part" in issues[0].message


def test_a_real_term_in_the_typedef_fixture_still_resolves(
    typedef_registry: OntologyRegistry,
) -> None:
    """Guards the Typedef fix against passing because nothing resolves at all."""
    assert validate_terms(["SO:0000704"], typedef_registry, "curation/x.yaml") == []
