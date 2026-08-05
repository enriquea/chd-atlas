import shutil
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


def test_a_validity_mirror_missing_a_selected_column_is_reported_not_raised(
    tmp_path: Path,
) -> None:
    """Reported by Copilot on the promotion PR, reproduced here.

    `_mirrored_validity` guarded only `"disease" not in frame.columns`, then
    unconditionally did `frame.select(["disease", "disease_label", "gene"])`.
    An upstream rename dropping `disease_label` or `gene` -- not `disease`
    itself -- left `"disease" not in frame.columns` false, so the guard passed
    and `frame.select` raised `polars.exceptions.ColumnNotFoundError`,
    unhandled, straight out of `validate_repository`. That breaks the one
    guarantee this module exists to keep (`CLAUDE.md`: "validators report; the
    build refuses") -- one malformed mirror must cost one issue, not the
    whole run, and `validate_table` already reports exactly this as TBL001
    over the same file, so nothing is lost by skipping the mirror here rather
    than crashing on it.

    Watched failing first: with the fix in `_mirrored_validity` reverted to
    checking only `"disease" not in frame.columns`, this test raised
    `polars.exceptions.ColumnNotFoundError: unable to find column
    "disease_label"` instead of returning a report -- confirmed both via this
    test and via `chd-atlas validate` against a scratch copy of the real
    repository with the same column dropped.

    A copy of the real, valid committed corpus rather than a minimal fixture,
    so `disease_label` is the *only* thing wrong with it and every other
    check -- including the second mirror, GenCC, read successfully -- still
    runs, isolating this guard's own effect from `TBL012`'s ("both mirrors
    unreadable") different failure mode.
    """
    source = tmp_path / "repo"
    for name in ("curation", "mirrors"):
        shutil.copytree(REPO_ROOT / name, source / name)
    (source / "ontologies").symlink_to(REPO_ROOT / "ontologies")

    clingen = source / "mirrors" / "clingen_gene_validity.tsv"
    header, *rows = clingen.read_text(encoding="utf-8").splitlines()
    columns = header.split("\t")
    label_index = columns.index("disease_label")
    del columns[label_index]
    lines = ["\t".join(columns)]
    for row in rows:
        cells = row.split("\t")
        del cells[label_index]
        lines.append("\t".join(cells))
    clingen.write_text("\n".join(lines) + "\n", encoding="utf-8")

    report = validate_repository(source)

    codes = [issue.code for issue in report.issues]
    assert "TBL001" in codes
    assert "TBL012" not in codes, "GenCC alone should still populate the mirrored-validity map"
    assert report.ok is False


def _real_repo(tmp_path: Path) -> Path:
    """A copy of the committed repository, which validates at the baseline.

    Copied rather than fixtured because the defect below is about *scale* -- one
    error naming a cause versus one error per row naming symptoms -- and a
    two-row fixture cannot show the difference.
    """
    for name in ("curation", "mirrors", "ontologies"):
        shutil.copytree(REPO_ROOT / name, tmp_path / name)
    return tmp_path


def test_an_absent_cohort_registry_reports_its_cause_once_not_every_row(
    tmp_path: Path,
) -> None:
    """CUR004, and the guarantee in CLAUDE.md §1 it exists to keep.

    "A missing registry must not report as hundreds of dangling references."
    The first version of `validate_burden_references` passed
    `{str(c.id) for c in corpus.cohorts}` unconditionally, so an absent
    `cohorts.yaml` -- which `load_curation` treats as optional and reports
    nothing for -- produced four BUR009 errors naming four cohorts and nothing
    naming the file. Measured 2026-08-05 before the fix: 4 errors, all BUR009.
    """
    root = _real_repo(tmp_path)
    (root / "curation" / "cohorts.yaml").unlink()

    report = validate_repository(root)
    codes = [issue.code for issue in report.issues if issue.severity is Severity.ERROR]

    assert codes == ["CUR004"]
    assert report.ok is False


def test_an_unreadable_gene_registry_does_not_multiply_by_the_burden_mirror(
    tmp_path: Path,
) -> None:
    """The same cascade from the other side, and the one that grew 146x.

    A header-only `mirrors/genes.tsv` reads fine and carries the `hgnc_id`
    column, so `_known_genes` returns an empty set rather than None and TBL008
    does not fire. Before `mirrors/burden.tsv` existed that input produced
    **1** error (REF001, the one curated assertion's gene). Measured 2026-08-05
    with the burden table added and the registries passed unconditionally:
    **146** -- 1 REF001 plus one BUR011 per distinct gene in the mirror -- none
    of which named the registry.

    Passing `known_genes or None` restores it to 1. The count is asserted, not
    just the code, because the defect was never a wrong code: BUR011 was right
    about all 145 genes. It was the burying that was wrong.
    """
    root = _real_repo(tmp_path)
    genes = root / "mirrors" / "genes.tsv"
    genes.write_text(genes.read_text(encoding="utf-8").splitlines()[0] + "\n", encoding="utf-8")

    report = validate_repository(root)
    codes = [issue.code for issue in report.issues if issue.severity is Severity.ERROR]

    assert codes == ["REF001"]
