from pathlib import Path

from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.validate.runner import ValidationReport, validate_repository


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


def test_locations_are_repo_relative(tmp_path: Path) -> None:
    """Absolute paths would make CI output differ from a curator's machine."""
    (tmp_path / "mirrors" / "varaints").mkdir(parents=True)

    report = validate_repository(tmp_path)

    assert report.issues
    for issue in report.issues:
        assert not issue.location.startswith(str(tmp_path)), issue.location
