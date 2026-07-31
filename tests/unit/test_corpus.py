# tests/unit/test_corpus.py
from pathlib import Path

from chd_atlas.corpus import load_curation

ASSERTION_YAML = """\
assertions:
  - id: CHDA:AST:0000001
    gene: HGNC:11604
    phenotypes: [HP:0001631]
    lesion_groups: [septal]
    classification: definitive
    inheritance: [AD]
    mechanism: haploinsufficiency
    syndromic: both
    source_tier: own_curation
    curator: 0000-0002-1825-0097
    curated_on: 2026-07-01
    last_reviewed: 2026-07-15
    evidence:
      - publication: PMID:8988165
        evidence_class: genetic_case
        locator: {kind: supplementary, file: TableS3.xlsx, sheet: S3, row: 42}
        summary: De novo nonsense variant in a proband with ASD.
        strength: strong
"""

PUBLICATIONS_YAML = """\
publications:
  - id: PMID:8988165
    title: Mutations in human TBX5 cause limb and cardiac malformation
    journal: Nature Genetics
    year: 1997
    authors: [Basson CT]
    study_type: case_control
"""


def _write_minimal_corpus(root: Path) -> None:
    (root / "curation" / "assertions").mkdir(parents=True)
    (root / "curation" / "assertions" / "TBX5.yaml").write_text(ASSERTION_YAML)
    (root / "curation" / "publications.yaml").write_text(PUBLICATIONS_YAML)


def test_loads_assertions_and_publications(tmp_path: Path) -> None:
    _write_minimal_corpus(tmp_path)

    corpus, issues = load_curation(tmp_path)

    assert issues == []
    assert len(corpus.assertions) == 1
    assert corpus.assertions[0].gene == "HGNC:11604"
    assert len(corpus.publications) == 1


def test_schema_violation_becomes_an_issue_rather_than_an_exception(tmp_path: Path) -> None:
    _write_minimal_corpus(tmp_path)
    broken = ASSERTION_YAML.replace("classification: definitive", "classification: very_sure")
    (tmp_path / "curation" / "assertions" / "TBX5.yaml").write_text(broken)

    corpus, issues = load_curation(tmp_path)

    assert corpus.assertions == ()
    assert len(issues) == 1
    assert issues[0].code == "SCHEMA001"
    assert "classification" in issues[0].message
    assert issues[0].location.endswith("TBX5.yaml")


def test_reports_every_schema_violation_not_just_the_first(tmp_path: Path) -> None:
    _write_minimal_corpus(tmp_path)
    broken = ASSERTION_YAML.replace("classification: definitive", "classification: very_sure")
    broken = broken.replace("mechanism: haploinsufficiency", "mechanism: magic")
    (tmp_path / "curation" / "assertions" / "TBX5.yaml").write_text(broken)

    _, issues = load_curation(tmp_path)

    assert len(issues) == 2


def test_malformed_yaml_becomes_a_parse_issue(tmp_path: Path) -> None:
    _write_minimal_corpus(tmp_path)
    (tmp_path / "curation" / "assertions" / "TBX5.yaml").write_text("assertions: [oops\n")

    _, issues = load_curation(tmp_path)

    assert [issue.code for issue in issues] == ["YAML001"]


def test_missing_curation_directory_is_reported(tmp_path: Path) -> None:
    _, issues = load_curation(tmp_path)

    assert [issue.code for issue in issues] == ["CORPUS001"]
