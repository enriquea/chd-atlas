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
    assert "classification" in issues[0].location


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


def test_non_utf8_bytes_are_reported_rather_than_raised(tmp_path: Path) -> None:
    """A curator name saved as Latin-1 must not crash the whole validation run."""
    _write_minimal_corpus(tmp_path)
    (tmp_path / "curation" / "assertions" / "TBX5.yaml").write_bytes(
        ASSERTION_YAML.replace("0000-0002-1825-0097", "François").encode("latin-1")
    )

    _, issues = load_curation(tmp_path)

    assert [issue.code for issue in issues] == ["YAML001"]


def test_empty_yaml_file_is_reported(tmp_path: Path) -> None:
    """yaml.load("") returns None; that must not read as 'no records, no problem'."""
    _write_minimal_corpus(tmp_path)
    (tmp_path / "curation" / "publications.yaml").write_text("")

    corpus, issues = load_curation(tmp_path)

    assert corpus.publications == ()
    assert [issue.code for issue in issues] == ["SCHEMA001"]


def test_assertions_from_multiple_files_are_concatenated(tmp_path: Path) -> None:
    """A regression to assignment instead of extend would pass every other test."""
    _write_minimal_corpus(tmp_path)
    second = ASSERTION_YAML.replace("CHDA:AST:0000001", "CHDA:AST:0000002").replace(
        "HGNC:11604", "HGNC:4173"
    )
    (tmp_path / "curation" / "assertions" / "GATA4.yaml").write_text(second)

    corpus, issues = load_curation(tmp_path)

    assert issues == []
    assert {a.id for a in corpus.assertions} == {"CHDA:AST:0000001", "CHDA:AST:0000002"}


def test_loads_phenotypes_and_featured(tmp_path: Path) -> None:
    _write_minimal_corpus(tmp_path)
    (tmp_path / "curation" / "phenotypes.yaml").write_text(
        "phenotypes:\n"
        "  - id: HP:0001631\n"
        "    label: Atrial septal defect\n"
        "    lesion_group: septal\n"
    )
    (tmp_path / "curation" / "featured.yaml").write_text(
        "featured:\n"
        "  - publication: PMID:8988165\n"
        "    order: 1\n"
        "    blurb: Founding TBX5 paper.\n"
        "    topic: genomics\n"
    )

    corpus, issues = load_curation(tmp_path)

    assert issues == []
    assert len(corpus.phenotypes) == 1
    assert len(corpus.featured) == 1
