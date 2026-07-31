import json
from pathlib import Path

from chd_atlas.schema_export import EXPORTED_MODELS, export_schemas


def test_writes_one_schema_file_per_model(tmp_path: Path) -> None:
    written = export_schemas(tmp_path)
    assert {path.name for path in written} == {
        f"{name}.schema.json" for name in EXPORTED_MODELS
    }
    for path in written:
        assert path.is_file()


def test_emitted_schema_is_valid_json_with_a_title(tmp_path: Path) -> None:
    export_schemas(tmp_path)
    payload = json.loads(
        (tmp_path / "assertion_file.schema.json").read_text(encoding="utf-8")
    )
    assert payload["title"] == "AssertionFile"


def test_export_is_byte_identical_across_runs(tmp_path: Path) -> None:
    export_schemas(tmp_path)
    first = (tmp_path / "assertion_file.schema.json").read_bytes()
    export_schemas(tmp_path)
    assert (tmp_path / "assertion_file.schema.json").read_bytes() == first


def test_committed_schemas_match_the_models() -> None:
    """Fails when a model changes without re-running `chd-atlas schemas export`."""
    committed = Path("schemas")
    if not committed.is_dir():
        return
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        fresh = Path(tmp)
        export_schemas(fresh)
        for path in sorted(fresh.glob("*.schema.json")):
            target = committed / path.name
            assert target.is_file(), f"{path.name} is not committed"
            assert target.read_text(encoding="utf-8") == path.read_text(encoding="utf-8"), (
                f"{path.name} is stale; run `chd-atlas schemas export`"
            )
