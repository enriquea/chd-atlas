# tests/unit/test_fs.py
import os
from pathlib import Path

import pytest

from chd_atlas.fs import list_dir, write_bytes_atomically


def test_lists_entries_in_sorted_order(tmp_path: Path) -> None:
    (tmp_path / "b.tsv").write_text("")
    (tmp_path / "a.tsv").write_text("")

    entries, issues = list_dir(tmp_path, "TBL009")

    assert [entry.name for entry in entries] == ["a.tsv", "b.tsv"]
    assert issues == []


def test_reports_an_unreadable_directory_instead_of_raising(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The guarantee is that one unreadable directory costs one issue, not the run.

    The failure is injected rather than provoked with `chmod 000` because mode
    bits do not restrain root, and CI may run as root in a container — the test
    would then pass by not reproducing the condition at all.
    """

    def deny(self: Path) -> object:
        raise PermissionError(13, "Permission denied")

    monkeypatch.setattr(Path, "iterdir", deny)

    entries, issues = list_dir(tmp_path, "TBL009")

    assert entries == []
    assert [issue.code for issue in issues] == ["TBL009"]
    assert "could not list directory" in issues[0].message


@pytest.mark.skipif(os.geteuid() == 0, reason="mode bits do not restrain root")
def test_reports_a_genuinely_unreadable_directory(tmp_path: Path) -> None:
    """The same guarantee against a real permission error, not an injected one."""
    locked = tmp_path / "locked"
    locked.mkdir()
    locked.chmod(0o000)
    try:
        entries, issues = list_dir(locked, "CUR001")
    finally:
        locked.chmod(0o755)

    assert entries == []
    assert [issue.code for issue in issues] == ["CUR001"]


@pytest.mark.parametrize("name", [".DS_Store", "Thumbs.db"])
def test_file_manager_droppings_are_not_listed(tmp_path: Path, name: str) -> None:
    """Opening the folder in Finder must not fail the build."""
    (tmp_path / name).write_text("")
    (tmp_path / "genes.tsv").write_text("")

    entries, issues = list_dir(tmp_path, "TBL009")

    assert [entry.name for entry in entries] == ["genes.tsv"]
    assert issues == []


def test_a_dotfile_the_sweeps_expect_is_still_listed(tmp_path: Path) -> None:
    """`.id_registry.yaml` is a dotfile that must stay visible."""
    (tmp_path / ".id_registry.yaml").write_text("")

    entries, _ = list_dir(tmp_path, "CUR001")

    assert [entry.name for entry in entries] == [".id_registry.yaml"]


def test_atomic_write_replaces_content_completely(tmp_path: Path) -> None:
    path = tmp_path / "out.json"
    path.write_bytes(b"old content that is longer than the new")

    write_bytes_atomically(path, b"new")

    assert path.read_bytes() == b"new"


def test_atomic_write_leaves_no_temporary_file_behind(tmp_path: Path) -> None:
    path = tmp_path / "out.json"

    write_bytes_atomically(path, b"payload")

    assert [entry.name for entry in tmp_path.iterdir()] == ["out.json"]


def test_a_failed_write_leaves_the_original_intact_and_no_debris(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Truncate-then-write would leave a half-written file the drift test misreads.

    The schema drift check compares committed bytes against freshly generated
    ones, so a truncated file reads as merely stale and sends a curator to re-run
    an export rather than telling them the file on disk is not a schema at all.

    The failure is injected at `os.replace` rather than at the write itself
    because `io.BufferedWriter` is an immutable C type and cannot be patched --
    the test would then fail on the patch and never exercise the guard.
    """
    path = tmp_path / "out.json"
    path.write_bytes(b"original")

    def explode(source: object, target: object) -> None:
        raise OSError(28, "No space left on device")

    monkeypatch.setattr("chd_atlas.fs.os.replace", explode)

    with pytest.raises(OSError):
        write_bytes_atomically(path, b"replacement")

    assert path.read_bytes() == b"original"
    # The temporary must be cleaned up, or a failed build leaves debris that the
    # stray-entry sweeps would later report as an unexpected file.
    assert [entry.name for entry in tmp_path.iterdir()] == ["out.json"]
