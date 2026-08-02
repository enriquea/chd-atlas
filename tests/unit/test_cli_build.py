# tests/unit/test_cli_build.py
import shutil
from pathlib import Path
from typing import Never

import pytest
from typer.testing import CliRunner

from chd_atlas.build.runner import BuildRefused
from chd_atlas.cli import app
from chd_atlas.tables import TABLE_SCHEMAS
from chd_atlas.validate.runner import validate_repository

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


@pytest.mark.parametrize("shape", ["existing file", "dangling symlink"])
def test_an_out_that_cannot_hold_a_site_is_refused_before_the_build(
    repo: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, shape: str
) -> None:
    """A destination that cannot work is a bad argument, not a bad repository.

    Both shapes exited 1 with a traceback before this guard existed, which is
    the code meaning "this repository has errors" — so a mistyped `--out` was
    indistinguishable from a corpus that genuinely fails, and the two are fixed
    by different people.

    The dangling symlink is the case a check written the obvious way misses:
    `Path.exists()` follows the link, so a broken one answers False to both
    `exists()` and `is_dir()`. Testing it is what makes `is_symlink()` in the
    guard load-bearing rather than decorative.

    `build_site` is a tripwire because exit 2 alone does not distinguish this
    guard from the `OSError` handler below it — measured: with only the handler,
    both shapes still exit 2. What the guard uniquely buys is refusing *before*
    a second and a half of validation is spent on an argument that cannot work,
    and saying which argument it was.
    """
    monkeypatch.setattr("chd_atlas.cli.build_site", _never_called)
    out = tmp_path / "out"
    if shape == "existing file":
        out.write_text("not a directory\n")
    else:
        out.symlink_to(tmp_path / "nowhere")

    result = runner.invoke(app, ["build", "--root", str(repo), "--out", str(out)])

    assert result.exit_code == 2, f"{shape} reported as a repository failure"
    assert "exists and is not a directory" in result.stdout
    assert str(out) in result.stdout


def test_a_destination_that_cannot_be_written_exits_2_rather_than_crashing(
    repo: Path, tmp_path: Path
) -> None:
    """The shapes no pre-check can see, caught where `schemas export` catches them.

    A read-only parent is indistinguishable from a writable one until the write
    is attempted, and the same is true of a full disk or a path that becomes
    unwritable while the build runs. Left uncaught these surfaced as a
    `NotADirectoryError` traceback and exit 1 — a stack trace where a sentence
    belongs, under the code that means the curation is wrong.
    """
    parent = tmp_path / "readonly"
    parent.mkdir()
    parent.chmod(0o555)
    out = parent / "site"

    try:
        result = runner.invoke(app, ["build", "--root", str(repo), "--out", str(out)])
    finally:
        parent.chmod(0o755)

    assert result.exit_code == 2
    assert "could not write the site to" in result.stdout
    assert str(out) in result.stdout


def test_a_build_guard_that_validate_does_not_cover_reports_instead_of_crashing(
    repo: Path, tmp_path: Path
) -> None:
    """The guards in `build/` are reachable by ordinary curation, not only by a bypass.

    `validate/` checks a variant shard's `chrom` column but not its filename, so
    `mirrors/variants/chr12.tsv` — schema-correct, wrongly named — passes at 0
    errors and 0 warnings and then raises inside `build_variants`. Measured; it
    is why the guard docstrings no longer claim reaching one means the gate was
    bypassed.

    Before this was handled, typer printed a rich traceback and exited 1: the
    right code, attached to a stack trace instead of a sentence, under a message
    that named no file. A curator's first move would have been to run `validate`,
    which reports nothing.

    The partial-site warning is asserted because it is the operational half. This
    failure arrives after writing has begun, so unlike a `BuildRefused` there may
    be a `dist/` on disk — one with no manifest, which a deploy step must not
    upload.
    """
    columns = [column.name for column in TABLE_SCHEMAS["variants"].columns]
    shards = repo / "mirrors" / "variants"
    shards.mkdir(parents=True, exist_ok=True)
    (shards / "chr12.tsv").write_text("\t".join(columns) + "\n")
    assert validate_repository(repo).ok is True, "the premise is that validation passes"
    out = tmp_path / "dist"

    result = runner.invoke(app, ["build", "--root", str(repo), "--out", str(out)])

    assert result.exit_code == 1
    assert "chr12" in result.stdout, "the message must name the offending shard"
    assert "must not be published" in result.stdout
    assert result.exception is None or isinstance(result.exception, SystemExit), (
        "a curation problem must not surface as a traceback"
    )


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
