import shutil
from pathlib import Path

from chd_atlas.validate.runner import validate_repository

REPO_ROOT = Path(__file__).parent.parent


def test_the_committed_repository_validates_cleanly() -> None:
    report = validate_repository(REPO_ROOT)
    assert report.ok, report.render()


def test_the_seed_corpus_is_not_empty() -> None:
    from chd_atlas.corpus import load_curation

    corpus, _ = load_curation(REPO_ROOT)
    assert corpus.assertions
    assert corpus.publications
    assert corpus.featured


def test_the_exported_schemas_are_committed() -> None:
    """test_committed_schemas_match_the_models returns early when schemas/ is
    absent, so without this it passes vacuously and never checks for drift."""
    from chd_atlas.schema_export import EXPORTED_MODELS

    committed = {path.name for path in (REPO_ROOT / "schemas").glob("*.schema.json")}
    assert committed == {f"{name}.schema.json" for name in EXPORTED_MODELS}


def test_a_misnamed_assertion_in_the_real_corpus_fails_validation(tmp_path: Path) -> None:
    """The case that motivated the stray-entry sweep: renaming an assertion to
    .yml made it invisible to every check while validate still exited 0."""
    for item in ("curation", "mirrors", "ontologies"):
        shutil.copytree(REPO_ROOT / item, tmp_path / item)
    assert validate_repository(tmp_path).ok, "the copied corpus should start clean"

    assertions = tmp_path / "curation" / "assertions"
    original = next(assertions.glob("*.yaml"))
    original.rename(original.with_suffix(".yml"))

    report = validate_repository(tmp_path)

    assert report.ok is False
    assert "CUR001" in [issue.code for issue in report.issues]


def test_missing_validity_mirrors_fails_validation(tmp_path: Path) -> None:
    """Reproduces the gap found reviewing Task 10's wiring: with both validity
    mirrors gone and everything else in the corpus valid, `_mirrored_validity`
    returns None and only SCP000 -- a WARNING -- was reported. `ok` ignores
    warnings, so `build_site` (which gates on `report.ok`) would have
    published a site where scope, and so which genes exist at all, could not
    be checked, while `chd-atlas validate` still exited 0. TBL012 makes the
    missing mirrors an ERROR, matching how TBL008 covers a missing gene
    registry the same way.

    `is False`, not `== False`: `{"f": 0} == {"f": False}` is `True` in
    Python, and this project pins boolean contracts with `is`.
    """
    for item in ("curation", "mirrors", "ontologies"):
        shutil.copytree(REPO_ROOT / item, tmp_path / item)
    assert validate_repository(tmp_path).ok, "the copied corpus should start clean"

    (tmp_path / "mirrors" / "clingen_gene_validity.tsv").unlink()
    (tmp_path / "mirrors" / "gencc_submissions.tsv").unlink()

    report = validate_repository(tmp_path)

    codes = [issue.code for issue in report.issues]
    assert report.ok is False
    assert "TBL012" in codes
    assert "SCP000" in codes
