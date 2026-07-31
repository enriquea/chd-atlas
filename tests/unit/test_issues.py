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
