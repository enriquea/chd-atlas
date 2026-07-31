from collections.abc import Callable
from pathlib import Path

import pytest

from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.validate.runner import ValidationReport, validate_repository

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

VALID_SOURCES_YAML = """sources:
  - id: clinvar
    name: ClinVar
    version: "2026-07"
    retrieved_on: 2026-07-15
    url: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/
    licence: public-domain
    redistribution: permitted
"""


def _corpus_load_fails(root: Path) -> None:
    """No curation directory, but a source registry that loads cleanly."""
    (root / "mirrors").mkdir()
    (root / "mirrors" / "sources.yaml").write_text(VALID_SOURCES_YAML)


def _source_registry_load_fails(root: Path) -> None:
    """A curation directory that loads cleanly, but no source registry."""
    (root / "curation").mkdir()


def _both_loads_fail(root: Path) -> None:
    """Neither the corpus nor the source registry is present."""


def test_report_counts_errors_and_warnings_separately() -> None:
    report = ValidationReport(
        issues=[
            ValidationIssue("A001", Severity.ERROR, "x", "boom"),
            ValidationIssue("B001", Severity.WARNING, "y", "hmm"),
        ]
    )
    assert report.error_count == 1
    assert report.warning_count == 1
    assert report.ok is False


def test_empty_report_is_ok() -> None:
    assert ValidationReport(issues=[]).ok is True


def test_issues_are_sorted_and_deduplicated() -> None:
    duplicate = ValidationIssue("A001", Severity.ERROR, "x", "boom")
    report = ValidationReport(
        issues=[
            ValidationIssue("B001", Severity.ERROR, "y", "later"),
            duplicate,
            duplicate,
        ]
    )
    assert [i.code for i in report.issues] == ["A001", "B001"]


def test_validate_repository_reports_a_missing_corpus(tmp_path: Path) -> None:
    report = validate_repository(tmp_path)
    codes = {issue.code for issue in report.issues}
    assert "CORPUS001" in codes
    assert "SRC001" in codes
    assert report.ok is False


def test_a_misnamed_mirror_directory_is_reported(tmp_path: Path) -> None:
    """A typo'd rename must not make the repository look clean."""
    (tmp_path / "mirrors" / "varaints").mkdir(parents=True)

    report = validate_repository(tmp_path)

    assert any(issue.code == "TBL009" for issue in report.issues)


def test_a_shard_directory_replaced_by_a_file_is_reported(tmp_path: Path) -> None:
    (tmp_path / "mirrors").mkdir(parents=True)
    (tmp_path / "mirrors" / "variants").write_text("oops\n")

    report = validate_repository(tmp_path)

    assert any(issue.code == "TBL009" for issue in report.issues)


def test_expected_mirror_entries_are_not_reported(tmp_path: Path) -> None:
    (tmp_path / "mirrors" / "variants").mkdir(parents=True)
    (tmp_path / "mirrors" / "genes.tsv").write_text("hgnc_id\n")
    (tmp_path / "mirrors" / "sources.yaml").write_text("sources: []\n")

    report = validate_repository(tmp_path)

    assert not any(issue.code == "TBL009" for issue in report.issues)


def test_a_missing_shard_directory_is_not_reported(tmp_path: Path) -> None:
    """No variants curated yet is a legitimate state, not an error."""
    (tmp_path / "mirrors").mkdir(parents=True)

    report = validate_repository(tmp_path)

    assert not any(issue.code == "TBL009" for issue in report.issues)


def test_ontology_checks_are_skipped_when_the_registry_fails(tmp_path: Path) -> None:
    """41 ONT003s would otherwise bury the one SRC001 that caused them."""
    report = validate_repository(tmp_path)

    codes = [issue.code for issue in report.issues]
    assert "ONT000" in codes
    assert "ONT003" not in codes


def test_a_failed_registry_does_not_cascade_one_error_per_term(tmp_path: Path) -> None:
    """With terms to check, the unguarded run reports ONT003 once per term."""
    (tmp_path / "curation").mkdir()
    (tmp_path / "curation" / "phenotypes.yaml").write_text(
        "phenotypes:\n"
        "  - id: HP:0001631\n"
        "    label: Atrial septal defect\n"
        "    lesion_group: septal\n"
        "  - id: HP:0001629\n"
        "    label: Ventricular septal defect\n"
        "    lesion_group: septal\n"
    )

    report = validate_repository(tmp_path)

    codes = [issue.code for issue in report.issues]
    assert codes.count("ONT000") == 1
    assert "ONT003" not in codes


@pytest.mark.parametrize(
    ("setup", "expected_warnings", "causing_errors"),
    [
        pytest.param(_corpus_load_fails, {"REF000"}, {"CORPUS001"}, id="corpus-failed"),
        pytest.param(
            _source_registry_load_fails,
            {"SRC000", "ONT000"},
            {"SRC001"},
            id="source-registry-failed",
        ),
        pytest.param(
            _both_loads_fail,
            {"REF000", "SRC000", "ONT000"},
            {"CORPUS001", "SRC001"},
            id="both-failed",
        ),
    ],
)
def test_every_skip_warning_arrives_with_the_error_that_caused_it(
    tmp_path: Path,
    setup: Callable[[Path], None],
    expected_warnings: set[str],
    causing_errors: set[str],
) -> None:
    """`ValidationReport.ok` ignores warnings, so a warnings-only report exits 0.

    That is only safe while every warning means "a check was skipped" and each
    arrives with the error that caused the skip. This pins both halves for all
    three skip paths: the exact set of warnings each state produces, and the
    error that explains them.

    SRC000 and ONT000 share one guard — a source registry that did not load —
    so no repository state triggers either alone; they are pinned as the pair
    they are. REF000 is isolable, and the third state covers all three at once.
    """
    setup(tmp_path)

    report = validate_repository(tmp_path)

    warnings = {issue.code for issue in report.issues if issue.severity is Severity.WARNING}
    errors = {issue.code for issue in report.issues if issue.severity is Severity.ERROR}
    assert warnings == expected_warnings
    assert causing_errors <= errors
    assert report.error_count > 0
    assert report.ok is False


def test_render_separates_the_summary_from_the_issues(tmp_path: Path) -> None:
    """Without a blank line the count runs straight on from the last issue."""
    report = validate_repository(tmp_path)

    lines = report.render().splitlines()
    assert lines[-1] == f"{report.error_count} error(s), {report.warning_count} warning(s)"
    assert lines[-2] == ""


def test_locations_are_repo_relative(tmp_path: Path) -> None:
    """Absolute paths would make CI output differ from a curator's machine."""
    (tmp_path / "mirrors" / "varaints").mkdir(parents=True)

    report = validate_repository(tmp_path)

    assert report.issues
    for issue in report.issues:
        assert not issue.location.startswith(str(tmp_path)), issue.location


def _seeded(tmp_path: Path) -> Path:
    """A corpus with one assertion, enough to make the guards below reachable."""
    (tmp_path / "curation" / "assertions").mkdir(parents=True)
    (tmp_path / "curation" / "assertions" / "TBX5.yaml").write_text(ASSERTION_YAML)
    return tmp_path


def test_a_missing_phenotype_vocabulary_is_an_error(tmp_path: Path) -> None:
    """Without it REF007/REF009/REF010 silently no-op. An error, not a warning,
    because ValidationReport.ok ignores warnings."""
    report = validate_repository(_seeded(tmp_path))

    assert "CUR002" in [issue.code for issue in report.issues]
    assert report.ok is False


def test_a_missing_publication_registry_is_an_error(tmp_path: Path) -> None:
    """Otherwise every citation reports REF002 and nothing names the cause."""
    report = validate_repository(_seeded(tmp_path))

    assert "CUR003" in [issue.code for issue in report.issues]


def test_ontology_checks_are_skipped_when_a_pinned_release_fails(tmp_path: Path) -> None:
    """ONT003 says the prefix is unpinned, which is false when the file is
    merely unreadable - it would send a curator to add a pin that exists."""
    root = _seeded(tmp_path)
    (root / "mirrors").mkdir(parents=True, exist_ok=True)
    (root / "mirrors" / "sources.yaml").write_text(
        "sources:\n"
        "  - id: hpo\n"
        "    name: Human Phenotype Ontology\n"
        '    version: "2026-06-23"\n'
        "    retrieved_on: 2026-07-31\n"
        "    url: https://hpo.jax.org/\n"
        "    licence: custom\n"
        "    redistribution: permitted_with_attribution\n"
        "    ontology_prefix: HP\n"
        "    ontology_file: ontologies/broken.obo\n"
    )
    (root / "ontologies").mkdir()
    (root / "ontologies" / "broken.obo").write_text("not an ontology\n")

    report = validate_repository(root)

    codes = [issue.code for issue in report.issues]
    assert "ONT004" in codes
    assert "ONT000" in codes
    assert "ONT003" not in codes
