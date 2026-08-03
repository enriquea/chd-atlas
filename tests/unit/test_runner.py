from collections.abc import Callable
from pathlib import Path

import pytest

from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.validate.runner import ValidationReport, validate_repository

REPO_ROOT = Path(__file__).parent.parent.parent

ASSERTION_YAML = """\
assertions:
  - id: CHDA:AST:0000001
    gene: HGNC:11604
    phenotypes: [HP:0001631]
    lesion_groups: [septal]
    inheritance: [AD]
    mechanism: haploinsufficiency
    syndromic: both
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


def _validity_mirrors_missing(root: Path) -> None:
    """A corpus, gene registry, ID registry and source registry that all load
    cleanly, but neither validity mirror is present.

    Everything else is made to load cleanly on purpose, so that TBL012 is the
    *only* error in the report — unlike the "source-registry-failed" case
    below, where SCP000 rides alongside SRC001 by fixture coincidence, not
    because TBL012 caused it. This isolates SCP000's own causing error.
    """
    (root / "curation").mkdir()
    (root / "curation" / ".id_registry.yaml").write_text("{}\n")
    (root / "mirrors").mkdir()
    (root / "mirrors" / "genes.tsv").write_text("hgnc_id\n")
    (root / "mirrors" / "sources.yaml").write_text(VALID_SOURCES_YAML)


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


def test_a_gap_warning_is_reported_without_blocking_the_build() -> None:
    """The second kind of warning, and why `ok` may ignore one that arrives alone.

    A *skip* warning means a check did not run, and is ignored safely only
    because an error always accompanies it — the three parametrised cases below
    pin that. REF013 is a *gap* warning: the check ran and found curated evidence
    that will reach no reader, and it deliberately arrives alone.

    That is a real weakening of the gate, not an oversight, and it is asserted
    here so that it is a decision on record. The assertion is sound and the site
    is publishable; what is missing is one cell in a mirror, and refusing a whole
    deploy over it is the wrong trade for a corpus curated incrementally. What
    makes it acceptable is the alternative: before REF013 the same gap was
    reported nowhere at all, on a build that said 0 errors and 0 warnings.

    The report is constructed rather than provoked because the claim is about
    `ok`'s semantics, not about when the validator fires — that is
    `test_referential.py`'s subject, where the rule lives.
    """
    gap = ValidationReport(
        issues=[
            ValidationIssue("REF013", Severity.WARNING, "mirrors/genes.tsv", "no uniprot accession")
        ]
    )

    assert gap.warning_count == 1
    assert gap.error_count == 0
    assert gap.ok is True, "a gap must not block a deploy; it must be visible in the report"


@pytest.mark.parametrize(
    ("setup", "expected_warnings", "causing_errors"),
    [
        pytest.param(_corpus_load_fails, {"REF000"}, {"CORPUS001"}, id="corpus-failed"),
        pytest.param(
            _source_registry_load_fails,
            {"SRC000", "ONT000", "SCP000"},
            {"SRC001"},
            id="source-registry-failed",
        ),
        pytest.param(
            _both_loads_fail,
            {"REF000", "SRC000", "ONT000"},
            {"CORPUS001", "SRC001"},
            id="both-failed",
        ),
        pytest.param(
            _validity_mirrors_missing,
            {"SCP000"},
            {"TBL012"},
            id="validity-mirrors-missing",
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

    SCP000 also appears in the "source-registry-failed" case, but not because
    it shares SRC000/ONT000's guard: `_source_registry_load_fails` builds a
    `curation/` directory with no `mirrors/` at all, so the validity-mirror
    read (`_mirrored_validity`) independently returns None. It is a fixture
    coincidence, not a shared cause — `causing_errors` still names only SRC001
    there, and that case's `report.ok is False` holds because of SRC001
    regardless of SCP000.

    "validity-mirrors-missing" is what actually isolates SCP000's own cause:
    every other check is made to load cleanly (gene registry, ID registry,
    source registry all present and valid) and only the two validity mirrors
    are absent, so TBL012 is the *only* error in that report. Before TBL012
    existed, this exact state reported SCP000 as its only warning and
    `report.ok` was True — a corpus whose scope could not be checked at all,
    validating clean. Measured directly: with the TBL012 emission removed from
    `validate_repository`, this case's `report.ok is False` assertion fails
    (`assert True is False`).
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


def test_scope_checks_run_against_the_real_repository() -> None:
    """The committed repository is in scope and reports no scope error.

    A smoke test in the sense CLAUDE.md permits: it is the only thing that
    exercises the wiring from `validate_repository` down through the mirror
    read into `validate_scope_terms`, `scope_candidates` and
    `validate_curation_is_in_scope`.

    Asserting only `report.error_count == 0` is not enough here -- that
    passes identically whether the scope checks ran or were never called.
    Measured directly: with the scope-check call site in `validate_repository`
    reverted entirely, `.venv/bin/python -m pytest` on a naive version of this
    test asserting only `report.error_count == 0` still passed (1 passed) --
    it cannot tell a live wire from a dead one. So this test also asserts that
    a scope check actually produced output.

    SCP003 (a mirrored cardiac-looking disease term that is not on the scope
    list) is emitted whenever such a term exists, and 268 (net 1, label
    keyword) plus 223 (net 2, gene bridge) were measured on 2026-08-03 against
    the committed 68-term `curation/chd_scope.yaml` and the committed
    ClinGen/GenCC mirrors, so it fires today. If the scope list ever grows to
    cover every candidate, replace the `"SCP003" in codes` assertion with one
    that the mirror-reading helper (`validate/runner.py::_mirrored_validity`)
    returned a non-empty result -- do not delete the guard, or this test
    degrades back into the worthless version described above.
    """
    report = validate_repository(REPO_ROOT)

    codes = {issue.code for issue in report.issues}
    assert "SCP003" in codes, "scope checks did not run"
    assert report.error_count == 0, report.render()
