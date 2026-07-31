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
