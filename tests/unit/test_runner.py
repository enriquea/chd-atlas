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


def _gene_registry_missing(root: Path) -> None:
    """Validity mirrors present, gene registry absent. Isolates GEN000's cause.

    The inverse of `_validity_mirrors_missing`: here there *is* symbol work to
    do -- both validity mirrors print a `gene_symbol` for every row -- and
    nothing to resolve it against, which is the only state in which GEN000 is
    supposed to fire. Both mirrors carry a header alone, which is enough for
    `_mirror_gene_labels` to read the columns and find no rows... so a single
    real row is written into ClinGen's, because a fixture with nothing to check
    would take the "no work, no skip" branch and pin the opposite of what this
    case is for.
    """
    (root / "curation").mkdir()
    (root / "curation" / ".id_registry.yaml").write_text("{}\n")
    (root / "mirrors").mkdir()
    (root / "mirrors" / "sources.yaml").write_text(VALID_SOURCES_YAML)
    (root / "mirrors" / "clingen_gene_validity.tsv").write_text(
        "gene\tgene_symbol\tdisease\tdisease_label\tmoi\tsop\tclassification\t"
        "classification_date\tgcep\treport_url\n"
        "HGNC:11604\tTBX5\tMONDO:0007732\tHolt-Oram syndrome\tAD\tSOP9\tDefinitive\t"
        "2023-01-01\tCHD\thttps://x\n"
    )
    (root / "mirrors" / "gencc_submissions.tsv").write_text(
        "sgc_id\tgene\tgene_symbol\tdisease\tdisease_label\tmoi\tclassification\t"
        "submitter\tsubmitted_on\treport_url\n"
    )


_GENES_TSV_HEADER = (
    "hgnc_id\tsymbol\tname\taliases\tensembl_gene\tncbi_gene\tlocus\tuniprot\tmane_select\n"
)
_PROFILES_TSV_HEADER = "dataset\tgene\ttissue\tstage\tmedian_abundance\tunit\tq25\tq75\tn_samples\n"
_EMPTY_CLINGEN_TSV = (
    "gene\tgene_symbol\tdisease\tdisease_label\tmoi\tsop\tclassification\t"
    "classification_date\tgcep\treport_url\n"
)
_EMPTY_GENCC_TSV = (
    "sgc_id\tgene\tgene_symbol\tdisease\tdisease_label\tmoi\tclassification\t"
    "submitter\tsubmitted_on\treport_url\n"
)
# A minimal, valid profile-design `Dataset` YAML, transcribed from
# `test_validate_profiles.py::_profile_dataset`'s field values so a profiles
# row naming this accession, tissue and stage does not also trip PRF004/005
# or REF005 -- noise this file's fixtures otherwise isolate away from.
_PROFILE_DATASET_YAML = """\
id: GSE999999
archive: geo
technology: bulk_rnaseq
design: profile
tissue: Heart
developmental_stage: embryonic
organism: NCBITaxon:9606
n_samples: 12
licence: CC BY 4.0
contrasts: []
cardiac_tissues: [Heart]
detection_floor: 1.0
floor_source: test fixture
quantile_estimator: linear
stages:
  - token: 7wpc
    wpc: 7.0
"""


def _repo_with_profiles_but_no_quantiles(root: Path) -> Path:
    """Profiles committed, quantiles not yet -- the ordinary state mid-curation.

    Absent is not unreadable, and only unreadable raises TBL000 (`mirror_paths`
    in `tables.py` yields a shard only when `path.is_file()`). Without a
    partner error, PRF000 would be a lone warning; `ValidationReport.ok`
    ignores warnings; the build would publish every gene with its percentile
    band silently missing and nothing on the site saying a figure is missing.

    Everything else here is made to load cleanly -- the gene registry (which
    also carries the one gene the profiles row cites, so PRF007 does not
    additionally fire), the curated dataset record (so neither PRF004/005 nor
    REF005 fires), the id registry, the source registry and both (empty)
    validity mirrors -- following `_gene_registry_missing`'s and
    `_validity_mirrors_missing`'s precedent, so PRF000/PRF010 are isolated to
    their own causing pair rather than riding along with unrelated noise.
    """
    (root / "curation" / "datasets").mkdir(parents=True)
    (root / "curation" / "datasets" / "GSE999999.yaml").write_text(_PROFILE_DATASET_YAML)
    (root / "curation" / ".id_registry.yaml").write_text("{}\n")
    (root / "mirrors").mkdir()
    (root / "mirrors" / "genes.tsv").write_text(
        _GENES_TSV_HEADER
        + "\t".join(["HGNC:11604", "TBX5", "T-box transcription factor 5", "", "", "", "", "", ""])
        + "\n"
    )
    (root / "mirrors" / "sources.yaml").write_text(VALID_SOURCES_YAML)
    (root / "mirrors" / "clingen_gene_validity.tsv").write_text(_EMPTY_CLINGEN_TSV)
    (root / "mirrors" / "gencc_submissions.tsv").write_text(_EMPTY_GENCC_TSV)
    (root / "mirrors" / "profiles").mkdir()
    (root / "mirrors" / "profiles" / "GSE999999.tsv").write_text(
        _PROFILES_TSV_HEADER
        + "\t".join(["GSE999999", "HGNC:11604", "Heart", "7wpc", "10.0", "tpm", "", "", "5"])
        + "\n"
    )
    return root


def _repo_with_no_profiles(root: Path) -> Path:
    """No `mirrors/profiles/` at all -- nothing here for PRF000 to check.

    Not simply an empty `tmp_path`: PRF000 must be provably absent when a
    `mirrors/` directory exists for other tables too, not merely when the
    whole repository happens to be empty by coincidence.
    """
    (root / "mirrors").mkdir(parents=True)
    return root


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
            {"SRC000", "ONT000", "SCP000", "SCP005"},
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
            {"SCP000", "SCP005"},
            {"TBL012"},
            id="validity-mirrors-missing",
        ),
        # GEN000's own cause, isolated the same way TBL012 isolates SCP000's.
        # A missing gene registry is TBL008, and without that error this skip
        # would arrive alone in a report `ok` calls clean -- a repository whose
        # every symbol was unresolvable, validating green. Issue #33.
        pytest.param(
            _gene_registry_missing,
            {"GEN000"},
            {"TBL008"},
            id="gene-registry-missing",
        ),
        # PRF000's own cause, isolated the same way TBL008 isolates GEN000's.
        # Profiles committed, quantiles not yet -- the ordinary state
        # mid-curation -- and without the partner error PRF010, this skip
        # would arrive alone in a report `ok` calls clean: every gene loses
        # its percentile band, silently, on a green build.
        pytest.param(
            _repo_with_profiles_but_no_quantiles,
            {"PRF000"},
            {"PRF010"},
            id="profiles-quantiles-absent",
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


def test_a_false_scope_attribution_is_caught_on_the_branch_that_actually_runs(
    tmp_path: Path,
) -> None:
    """SCP005 through `validate_repository`, on the mirror-READABLE branch.

    **The branch that runs on every real build was killed by 0 of 789 tests.**
    `validate_repository` calls `validate_scope_attribution` twice -- once when
    the mirrors are unreadable, once when they are -- and only the first was
    defended. Deleting the second call left the whole suite green while
    `chd-atlas validate` reported 0 errors on a corpus with an invented
    attribution, and the atlas would publish 24 pages asserting an external
    warrant nothing had checked. Found by adversarial review of #30; it is the
    §4.28 shape, a guard tested on the layer it was not needed.

    The forgery used here is the one the check exists for: a term whose only
    warrant is a commercial laboratory's GenCC submission, relabelled as
    ClinGen's chartered CHD panel. The unit test proves the function refuses it;
    this proves the function is reached.
    """
    root = _real_repo(tmp_path)
    scope = root / "curation" / "chd_scope.yaml"
    text = scope.read_text(encoding="utf-8")
    forged = text.replace(
        "    admitted_by: gencc_submitter\n    attributed_to: Ambry Genetics\n",
        "    admitted_by: clingen_chd_panel\n    attributed_to: Ambry Genetics\n",
        1,
    )
    assert forged != text, "no commercial-lab attribution left to forge; update this fixture"
    scope.write_text(forged, encoding="utf-8")

    report = validate_repository(root)

    assert "SCP005" in {issue.code for issue in report.issues}
    assert report.ok is False


def test_an_absent_quantile_mirror_is_an_error_not_a_lone_warning(tmp_path: Path) -> None:
    """Absent is not unreadable, and only unreadable raises TBL000.

    Profiles committed, quantiles not yet -- the ordinary state mid-curation.
    Without a partner error PRF000 is a lone warning, `ok` ignores warnings,
    the build publishes, and every gene loses its percentile band with nothing
    on the site saying a figure is missing.
    """
    root = _repo_with_profiles_but_no_quantiles(tmp_path)
    report = validate_repository(root)
    codes = {issue.code for issue in report.issues}
    assert "PRF000" in codes
    assert "PRF010" in codes  # the partner error
    assert not report.ok


def test_no_profiles_means_no_skip_warning_at_all(tmp_path: Path) -> None:
    """Emit a skip only when there was work to skip."""
    report = validate_repository(_repo_with_no_profiles(tmp_path))
    assert "PRF000" not in {issue.code for issue in report.issues}


@pytest.mark.parametrize(
    ("include_validity_mirrors", "expect_tbl012"),
    [
        pytest.param(True, False, id="mirrors-readable"),
        pytest.param(False, True, id="mirrors-missing"),
    ],
)
def test_prf_references_run_on_both_scope_branches(
    tmp_path: Path, include_validity_mirrors: bool, expect_tbl012: bool
) -> None:
    """`validate_profile_references` has one call site, and both branches of
    the `mirrored` if/else below it must still reach it.

    Section 4.34: a validator called on two branches is tested on the branch
    you were thinking about. `validate_scope_attribution` is called once on
    each side of this exact if/else, and only the mirror-unreadable side had a
    test -- deleting the mirror-readable call passed all 789 tests (see
    `test_a_false_scope_attribution_is_caught_on_the_branch_that_actually_runs`
    above). `validate_profile_references`'s inputs (`datasets`, `known_genes`,
    `phases`) do not depend on `mirrored` at all, so this module places one
    call site *before* that split rather than duplicating it inside both
    branches -- but "one call site, unconditional" is a claim about the
    source, and this is the measurement: a profiles row citing an
    unregistered gene must still report PRF007 whichever way the split
    resolves, proving a future refactor that "helpfully" moves the call inside
    just one branch -- the shape of the #30 bug -- would be caught.
    """
    (tmp_path / "curation").mkdir()
    (tmp_path / "mirrors" / "profiles").mkdir(parents=True)
    (tmp_path / "mirrors" / "profiles" / "GSE999999.tsv").write_text(
        _PROFILES_TSV_HEADER
        + "\t".join(["GSE999999", "HGNC:99999", "Heart", "7wpc", "10.0", "tpm", "", "", "5"])
        + "\n"
    )
    # A registry with a real row, NOT a bare header. An empty-but-readable
    # registry reaches `known_genes or None` as None, which skips the gene
    # check entirely -- so a header-only fixture would make PRF007 fire for the
    # wrong reason (the registry "did not load") and would keep passing even if
    # the unregistered-gene rule were deleted. The row below is a gene the
    # profiles mirror does NOT cite, so PRF007 fires on HGNC:99999 specifically.
    (tmp_path / "mirrors" / "genes.tsv").write_text(
        _GENES_TSV_HEADER
        + "\t".join(["HGNC:11604", "TBX5", "T-box 5", "", "ENSG00000089225", "6910", "", "", ""])
        + "\n"
    )
    if include_validity_mirrors:
        (tmp_path / "mirrors" / "clingen_gene_validity.tsv").write_text(_EMPTY_CLINGEN_TSV)
        (tmp_path / "mirrors" / "gencc_submissions.tsv").write_text(_EMPTY_GENCC_TSV)

    report = validate_repository(tmp_path)

    codes = {issue.code for issue in report.issues}
    assert ("TBL012" in codes) is expect_tbl012, "fixture landed on the wrong scope branch"
    assert "PRF007" in codes, "validate_profile_references did not run on this branch"


def test_validate_profiles_internal_consistency_runs_through_the_repository(
    tmp_path: Path,
) -> None:
    """`validate_profiles`'s own checks (PRF001/002/003/008) are unit-tested
    directly in `test_validate_profiles.py`; this proves the *wiring* reaches
    them too, the same distinction `test_scope_checks_run_against_the_real_
    repository` draws for the scope validators.

    Called unconditionally, right beside `validate_burden`'s own unconditional
    call -- deleting either is a defect of the same shape, and PRF008 is the
    cheapest of the four to provoke: one dataset's `profiles` rows reporting
    two different units needs no `profile_quantiles` shard at all.
    """
    (tmp_path / "mirrors" / "profiles").mkdir(parents=True)
    (tmp_path / "mirrors" / "profiles" / "GSE999999.tsv").write_text(
        _PROFILES_TSV_HEADER
        + "\t".join(["GSE999999", "HGNC:11604", "Heart", "7wpc", "10.0", "tpm", "", "", "5"])
        + "\n"
        + "\t".join(["GSE999999", "HGNC:11604", "Liver", "7wpc", "10.0", "rpkm", "", "", "5"])
        + "\n"
    )

    report = validate_repository(tmp_path)

    assert "PRF008" in {issue.code for issue in report.issues}


def test_prf009_fires_against_a_real_published_gene_population(tmp_path: Path) -> None:
    """`_gate_published_genes` computes an actual population, not a permanent no-op.

    PRF009 (a published gene missing a `profiles` row in a cell that was
    assayed) can only ever fire when `published_genes` is non-empty -- an
    implementation that always handed `validate_profile_references` an empty
    set would make this check permanently dead, and every other test in this
    file is deliberately built so PRF009 never fires, which would hide exactly
    that. This builds the smallest repository where `build.validity.
    published_genes` admits a real gene -- one ClinGen `Definitive` record for
    a disease `chd_scope.yaml` places in scope -- and proves that population
    reaches `validate_profile_references` through `validate_repository`, not
    merely through the private helper in isolation.
    """
    root = tmp_path
    (root / "curation").mkdir()
    (root / "curation" / "chd_scope.yaml").write_text(
        "diseases:\n"
        "  - id: MONDO:0007732\n"
        "    label: Holt-Oram syndrome\n"
        "    reason: ClinGen classifies this gene-disease pair as CHD.\n"
        "    admitted_by: clingen_chd_panel\n"
        "    attributed_to: Congenital Heart Disease\n"
        "    admitted_on: 2026-01-01\n"
    )
    (root / "mirrors").mkdir()
    (root / "mirrors" / "genes.tsv").write_text(
        _GENES_TSV_HEADER
        + "\t".join(["HGNC:11604", "TBX5", "T-box transcription factor 5", "", "", "", "", "", ""])
        + "\n"
    )
    (root / "mirrors" / "clingen_gene_validity.tsv").write_text(
        "gene\tgene_symbol\tdisease\tdisease_label\tmoi\tsop\tclassification\t"
        "classification_date\tgcep\treport_url\n"
        "HGNC:11604\tTBX5\tMONDO:0007732\tHolt-Oram syndrome\tAD\tSOP9\tDefinitive\t"
        "2023-01-01\tCongenital Heart Disease\thttps://x\n"
    )
    (root / "mirrors" / "gencc_submissions.tsv").write_text(_EMPTY_GENCC_TSV)
    (root / "mirrors" / "profile_quantiles").mkdir(parents=True)
    (root / "mirrors" / "profile_quantiles" / "GSE999999.tsv").write_text(
        "dataset\ttissue\tstage\tpercentile\tvalue\tunit\tn_genes\n"
        "GSE999999\tHeart\t7wpc\t50\t10.0\ttpm\t20000\n"
    )

    report = validate_repository(root)

    issue = next(i for i in report.issues if i.code == "PRF009")
    assert "HGNC:11604" in issue.message


def test_an_unmapped_classification_does_not_crash_validate_repository(tmp_path: Path) -> None:
    """`_gate_published_genes` must never let `build.validity.gene_validity` raise
    out of `validate_repository`.

    `gene_validity` is written for post-gate data and raises `ValueError` by
    design when a mirror carries a classification term the atlas vocabulary
    has not mapped yet -- reproduced directly against a synthetic frame while
    writing `_gate_published_genes` (see its docstring). That is the right
    behaviour for `build_site`, which must refuse rather than publish a
    confidence nobody vetted, but this call sits inside `validate_repository`
    to compute one WARNING-level check's population, and "validators report,
    they do not raise" (CLAUDE.md section 1) must hold regardless of what a
    third-party mirror contains -- including mid-curation, before
    `validate_table` has even finished reporting the same value as TBL004.
    """
    root = tmp_path
    (root / "curation").mkdir()
    (root / "curation" / "chd_scope.yaml").write_text(
        "diseases:\n"
        "  - id: MONDO:0007732\n"
        "    label: Holt-Oram syndrome\n"
        "    reason: ClinGen classifies this gene-disease pair as CHD.\n"
        "    admitted_by: clingen_chd_panel\n"
        "    attributed_to: Congenital Heart Disease\n"
        "    admitted_on: 2026-01-01\n"
    )
    (root / "mirrors").mkdir()
    (root / "mirrors" / "genes.tsv").write_text(_GENES_TSV_HEADER)
    (root / "mirrors" / "clingen_gene_validity.tsv").write_text(
        "gene\tgene_symbol\tdisease\tdisease_label\tmoi\tsop\tclassification\t"
        "classification_date\tgcep\treport_url\n"
        "HGNC:11604\tTBX5\tMONDO:0007732\tHolt-Oram syndrome\tAD\tSOP9\t"
        "Not A Real Classification\t2023-01-01\tCongenital Heart Disease\thttps://x\n"
    )
    (root / "mirrors" / "gencc_submissions.tsv").write_text(_EMPTY_GENCC_TSV)

    report = validate_repository(root)  # must return, never raise

    assert isinstance(report, ValidationReport)
    # The malformed value is still visible in the report, via the schema
    # check that runs independently of `_gate_published_genes` -- proving the
    # defensive fallback silences the crash, not the evidence of the defect.
    assert any(issue.code == "TBL004" for issue in report.issues)


def test_a_programmer_error_in_the_gate_propagates_rather_than_silencing_prf009(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The sibling of the test above, and the case it does NOT cover.

    `_gate_published_genes` absorbs `KeyError` and `ValueError` because those
    are the *data* failures its docstring names, and each already arrives with
    its own ERROR from `validate_table` -- so degrading to an empty published
    set cannot hide anything.

    A blanket `except Exception` also absorbed programmer errors, and those
    have no accompanying validation error to make them visible. Measured
    2026-08-19: injecting a `TypeError` into `gene_validity` -- the shape of an
    ordinary refactor regression -- produced `0 error(s), 4 warning(s)` with
    `ok=True` while PRF009 checked nothing at all. The report was green and one
    validator was dead.

    So this asserts the opposite of its sibling on purpose: a `ValueError` is
    absorbed and a `TypeError` is not. Raised by review on PR #39; the
    docstring's own guarantee had covered data failures and not code ones.
    """

    def boom(*args: object, **kwargs: object) -> dict[str, object]:
        raise TypeError("simulated refactor regression, not a data problem")

    monkeypatch.setattr("chd_atlas.validate.runner.gene_validity", boom)
    with pytest.raises(TypeError, match="simulated refactor regression"):
        validate_repository(Path("."))
