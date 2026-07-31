from pathlib import Path

from typer.testing import CliRunner

from chd_atlas.cli import app

runner = CliRunner()


def test_validate_exits_nonzero_on_an_invalid_repository(tmp_path: Path) -> None:
    result = runner.invoke(app, ["validate", "--root", str(tmp_path)])
    assert result.exit_code == 1
    assert "CORPUS001" in result.stdout
    assert "error(s)" in result.stdout


def test_schemas_export_writes_files(tmp_path: Path) -> None:
    result = runner.invoke(app, ["schemas", "export", "--target", str(tmp_path)])
    assert result.exit_code == 0
    assert (tmp_path / "assertion_file.schema.json").is_file()
    assert "6 schema(s)" in result.stdout


def test_help_lists_both_commands() -> None:
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "validate" in result.stdout
    assert "schemas" in result.stdout
