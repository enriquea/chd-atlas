# src/chd_atlas/fs.py
"""Filesystem access: reads that report failure instead of raising it, writes that land whole."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path
from typing import Final

from chd_atlas.issues import Severity, ValidationIssue

# Files a file manager writes just by opening a directory. They are never curated
# content and never committed — .gitignore excludes them — so reporting one as an
# unexpected entry would fail validation on a macOS curator's machine for a reason
# that has nothing to do with the data, and would clear itself on the next run.
#
# Named explicitly rather than filtered as "any dotfile", because
# curation/.id_registry.yaml is a dotfile the sweeps must keep seeing.
OS_METADATA: Final[frozenset[str]] = frozenset({".DS_Store", "Thumbs.db"})


def list_dir(directory: Path, code: str) -> tuple[list[Path], list[ValidationIssue]]:
    """Sorted entries of one directory, or one issue naming why they are unavailable.

    ``iterdir`` is the only filesystem call the stray-entry sweeps make, and it
    fails for reasons that have nothing to do with the curated data: a directory
    whose mode was set to 000, a stale network mount, a symlink to nowhere. An
    uncaught raise there aborts validation of the entire repository over one
    unreadable directory — the failure mode every other reader in this package is
    written to avoid, and the reason ``read_table`` and ``_read_yaml`` exist in
    the shape they do.

    Returning an empty listing alongside the issue is safe because the caller is
    looking for entries it does not expect: with nothing listed it finds nothing,
    and the issue returned here is what keeps that silence from reading as
    success.

    `OS_METADATA` entries are filtered out. Every caller loads real content by
    glob and uses this listing only to find entries nothing claims, so dropping
    file-manager droppings here cannot hide curated work.
    """
    try:
        return sorted(
            entry for entry in directory.iterdir() if entry.name not in OS_METADATA
        ), []
    except OSError as exc:
        issue = ValidationIssue(
            code,
            Severity.ERROR,
            str(directory),
            f"could not list directory: {exc}",
        )
        return [], [issue]


def write_bytes_atomically(path: Path, payload: bytes) -> None:
    """Replace ``path`` with ``payload`` in a single step.

    Opening the destination for writing truncates it before the payload lands,
    so a write that failed part-way left a truncated file while its siblings kept
    their previous content. Every consumer of these files — the schema drift
    test, the build manifest — compares bytes, and would report the corrupted
    file as merely stale rather than as not a file of that kind at all.

    Writing a sibling temporary and renaming it leaves the target as either its
    old content or the complete new content, never a prefix of the new one.
    """
    # A sibling so `os.replace` stays within one filesystem, where it is atomic.
    descriptor, name = tempfile.mkstemp(dir=path.parent, prefix=f"{path.name}.", suffix=".tmp")
    temporary = Path(name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
        # `mkstemp` creates at 0600, which would leave a published artifact
        # owner-only on whichever machine generated it.
        os.chmod(temporary, 0o644)
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
