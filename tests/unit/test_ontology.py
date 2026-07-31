# tests/unit/test_ontology.py
import shutil
from pathlib import Path

import pytest

from chd_atlas.validate.ontology import OntologyRegistry, validate_terms

FIXTURE = Path(__file__).parent.parent / "fixtures" / "mini-hp.obo"


@pytest.fixture
def registry(tmp_path: Path) -> OntologyRegistry:
    obo = tmp_path / "ontologies" / "hp-2026-06-01.obo"
    obo.parent.mkdir(parents=True)
    shutil.copy(FIXTURE, obo)
    return OntologyRegistry.from_files(tmp_path, {"HP": "ontologies/hp-2026-06-01.obo"})


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


def test_missing_pinned_file_is_reported(tmp_path: Path) -> None:
    registry = OntologyRegistry.from_files(tmp_path, {"HP": "ontologies/absent.obo"})
    assert [i.code for i in registry.load_issues] == ["ONT004"]
