from chd_atlas.issues import Severity, ValidationIssue


def test_format_renders_code_location_and_message() -> None:
    issue = ValidationIssue(
        code="REF001",
        severity=Severity.ERROR,
        location="curation/assertions/TBX5.yaml",
        message="unknown gene HGNC:99999",
    )
    assert issue.format() == (
        "ERROR [REF001] curation/assertions/TBX5.yaml: unknown gene HGNC:99999"
    )


def test_issues_sort_deterministically() -> None:
    a = ValidationIssue("REF001", Severity.ERROR, "b.yaml", "second")
    b = ValidationIssue("REF001", Severity.ERROR, "a.yaml", "first")
    assert sorted([a, b]) == [b, a]


def test_issues_are_hashable_so_duplicates_collapse() -> None:
    a = ValidationIssue("SCH001", Severity.WARNING, "x.tsv", "dup")
    b = ValidationIssue("SCH001", Severity.WARNING, "x.tsv", "dup")
    assert len({a, b}) == 1


def test_errors_sort_before_warnings() -> None:
    """Pins intent: adding a severity level whose value breaks this order must fail here."""
    warning = ValidationIssue("A001", Severity.WARNING, "x.yaml", "m")
    error = ValidationIssue("A001", Severity.ERROR, "x.yaml", "m")
    assert sorted([warning, error]) == [error, warning]


def test_severity_renders_as_its_value_not_its_repr() -> None:
    assert str(Severity.ERROR) == "error"
    assert f"{Severity.WARNING}" == "warning"


def test_str_matches_format() -> None:
    issue = ValidationIssue("A001", Severity.ERROR, "x.yaml", "m")
    assert str(issue) == issue.format()
