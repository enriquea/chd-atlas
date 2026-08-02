# tests/unit/test_cli_build.py
import shutil
from pathlib import Path
from typing import Never

import pytest
from typer.testing import CliRunner

from chd_atlas.build.runner import BuildRefused
from chd_atlas.cli import app

REPO = Path(__file__).parent.parent.parent
runner = CliRunner()


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    """A copy of the real committed corpus, which validates cleanly.

    `ontologies` is symlinked rather than copied: it is 11 MB of pinned OBO
    releases that nothing here changes, and the copy would cost more than the
    build does.
    """
    source = tmp_path / "repo"
    source.mkdir()
    for name in ("curation", "mirrors"):
        shutil.copytree(REPO / name, source / name)
    (source / "ontologies").symlink_to(REPO / "ontologies")
    return source


def _never_called(*args: object, **kwargs: object) -> Never:
    raise AssertionError("build_site ran; the argument guard should have returned first")


def test_build_writes_the_site_and_reports_what_it_wrote(repo: Path, tmp_path: Path) -> None:
    """The one end-to-end pass: a real corpus through the real builder.

    Everything else in this file injects `build_site`'s outcome, because the
    command's own job is arguments, exit codes and messages — what it builds is
    `test_build_runner.py`'s subject. This test is what keeps that substitution
    honest: it is the only thing here that would notice the command calling
    `build_site` with the wrong arguments, or the two drifting apart.

    The count is checked against the manifest rather than asserted as a literal,
    so the number the curator is told is the number of files that exist.
    """
    out = tmp_path / "dist"

    result = runner.invoke(app, ["build", "--root", str(repo), "--out", str(out)])

    assert result.exit_code == 0
    assert (out / "manifest.json").is_file()
    built = [path for path in out.rglob("*") if path.is_file()]
    assert f"wrote {len(built)} file(s) to {out}" in result.stdout


def test_a_refused_build_exits_1_and_says_why(
    repo: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Exit 1 and the report, because a curator ran this instead of `validate`.

    The exit code is what a deploy step branches on, and echoing the message is
    what saves the curator a second command. Both are pinned here: an exit of 0
    on a refusal would upload nothing while reporting success, and a swallowed
    message would send them to `validate` to find out what this already knew.

    The refusal is injected rather than provoked. That `build_site` refuses a
    genuinely broken corpus, and writes nothing when it does, is its own
    guarantee and is tested there; reproducing it here would spend a second full
    validation to re-answer a question this command does not ask.
    """
    rendered = "refusing to build: 1 validation error(s)\nSCHEMA001 curation/assertions/BROKEN.yaml"

    def refuse(*args: object, **kwargs: object) -> Never:
        raise BuildRefused(rendered)

    monkeypatch.setattr("chd_atlas.cli.build_site", refuse)
    out = tmp_path / "dist"

    result = runner.invoke(app, ["build", "--root", str(repo), "--out", str(out)])

    assert result.exit_code == 1
    assert "refusing to build" in result.stdout
    assert "SCHEMA001" in result.stdout
    assert not out.exists()


def test_a_root_that_is_not_a_directory_exits_2_without_building(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Exit 2 means "you pointed me at the wrong place", which is not exit 1.

    CI must be able to tell a mistyped `--root` from a repository that genuinely
    fails, because the first is a pipeline bug and the second is a curation one.
    `validate` makes the same distinction for the same reason.

    `build_site` is replaced with a tripwire rather than left alone: the guard's
    whole job is to return *before* the build, and an implementation that
    checked the path afterwards would still exit 2 while having spent a full
    validation on a directory that does not exist.
    """
    monkeypatch.setattr("chd_atlas.cli.build_site", _never_called)
    missing = tmp_path / "nope"

    result = runner.invoke(app, ["build", "--root", str(missing), "--out", str(tmp_path / "dist")])

    assert result.exit_code == 2
    assert str(missing) in result.stdout
    assert "not a directory" in result.stdout
