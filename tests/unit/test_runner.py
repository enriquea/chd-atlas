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
