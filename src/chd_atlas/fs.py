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

    Opening the destination for writing truncates it before the payload lands, so
    a write that failed part-way left a file that is neither its old content nor
    its new one. What that costs depends on the caller, and for one of them it is
    permanent:

    - `curation/.id_registry.yaml` is a counter. A truncated write rewinds every
      prefix to zero, and the next `allocate` then reissues identifiers already
      in use — the one thing this package's ID module promises never happens, and
      the one failure here that committed history cannot undo.
    - The committed `schemas/*.schema.json` are compared byte-for-byte by the
      drift test. A truncated one reads as merely stale, sending a curator to
      re-run the export rather than telling them the file is not a schema at all.
    - A `dist/` build artifact costs only a re-run. Its manifest checksum is
      taken from the complete payload, so a short file fails verification at the
      consumer instead of passing as valid content.

    Writing a sibling temporary and renaming it leaves the target as either its
    old content or the complete new content, never a prefix of the new one.

    ``path.parent`` must exist. `mkstemp` would otherwise raise
    `FileNotFoundError` naming the temporary it was about to create — a path the
    caller never chose and cannot correlate with anything they typed.

    Raising, rather than returning `ValidationIssue`s the way `list_dir` above
    does, is deliberate. `list_dir` swallows its `OSError` because an unreadable
    directory is a finding *about the curated data*, and one finding must cost
    one issue rather than the whole run. A failed write is not a finding about
    the data — there is no report for it to appear in, and the caller is not
    validating anything. `cli.py` catches `OSError` from `export_schemas` and
    exits 2, reporting "you pointed me somewhere I cannot write". Folding a
    failed write into a validation report instead would turn a lost file into a
    clean exit code, which is the silent loss this function exists to prevent.
    """
    # A sibling so `os.replace` stays within one filesystem, where it is atomic.
    descriptor, name = tempfile.mkstemp(dir=path.parent, prefix=f"{path.name}.", suffix=".tmp")
    temporary = Path(name)
    try:
        # Binary, and the payload is bytes, so the caller's newlines reach the
        # disk untouched. `open(..., "w", encoding="utf-8")` defaults to
        # `newline=None`, which rewrites every "\n" to `os.linesep`. On Windows
        # that emits CRLF, and since the drift test compares bytes it would then
        # report every committed schema as stale on every run.
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
        # `mkstemp` creates at 0600, which would leave a committed, world-readable
        # artifact owner-only on whichever machine generated it. Set the mode
        # explicitly so it does not depend on how the file happened to be made.
        #
        # Before the replace, not after: afterwards there is a window in which
        # the file is visible at its published name with the wrong mode, and the
        # chmod would apply to whatever now sits at `path` rather than to the
        # bytes just written. Reversing the two still passes every test here.
        os.chmod(temporary, 0o644)
        # No fsync: this guards against a failed or interrupted write, not
        # against power loss. These files are regenerated from committed sources
        # by a single command, so an unclean shutdown costs a re-run rather than
        # data — not worth an fsync on every write.
        os.replace(temporary, path)
    # BaseException, not Exception: a Ctrl-C between mkstemp and replace would
    # otherwise strand a `.tmp` in `curation/`, where the CUR001 stray-entry
    # sweep reports it as an unexpected file. `.gitignore` has no `*.tmp` rule
    # and `list_dir` filters only `OS_METADATA`, so nothing else would catch it.
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
