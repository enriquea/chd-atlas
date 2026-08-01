# src/chd_atlas/build/manifest.py
"""The manifest: what was built, from which commit, and a checksum per file.

This is what a consumer verifies a download against and what a cache keys on,
so `files` maps each published URL to the sha256 of the exact bytes served at
it — the compressed bytes for a `.gz`, not the JSON inside. A checksum that
disagrees with the file is worse than no manifest at all: a consumer verifying
a download rejects a file that is correct.

Nothing here varies between two builds of the same commit. A build timestamp
was considered and rejected: it would make the output irreproducible, so no
test could assert byte-identical rebuilds and the checksums would stop being a
pure function of the inputs. A consumer wanting a publication date reads the
commit's.

The manifest cannot list itself — its own checksum would have to be computed
from bytes containing that checksum. The exclusion is not a filter but an
ordering: `write_json` encodes the payload before it records the new entry, so
`manifest.json` enters `checksums` a moment after the bytes naming every other
file are already fixed.
"""

from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Final

from chd_atlas.build.emit import Emitter
from chd_atlas.corpus import Corpus

SCHEMA_VERSION: Final = "1.0"


def source_commit(root: Path) -> str | None:
    """The commit the build was made from, or None when there is not one.

    None rather than a raise: building from an unpacked tarball is legitimate
    and should produce a complete site, just one that cannot state its
    provenance. The same answer covers a `root` that does not exist and a
    checkout with no commits yet — `git -C` failing to change directory records
    exactly as much provenance as a missing `.git` does, so neither earns a
    branch of its own here.

    The `returncode` check is load-bearing rather than defensive, which is not
    obvious: in a repository with no commits `git rev-parse HEAD` exits 128 but
    prints the literal string `HEAD` on *stdout* (measured, git 2.50). Without
    the check this would publish `"source_commit": "HEAD"` — a provenance claim
    resolving to nothing, which is worse than the null a consumer can test for.

    `check=False` is not a second decision: `check=True` raises
    `CalledProcessError`, a `SubprocessError`, which the handler below already
    catches, so the two spellings behave identically. `timeout` is a real one —
    without it a git that blocks on a dead network mount or a credential prompt
    hangs the build with nothing to read.
    """
    try:
        completed = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=False,
            timeout=10,
        )
    # OSError is git absent from PATH, which a build must survive for the same
    # reason a tarball build must.
    except (OSError, subprocess.SubprocessError):
        return None
    if completed.returncode != 0:
        return None
    # `or None`, because "" is falsy but still reads as a recorded commit to a
    # consumer testing the key against null.
    return completed.stdout.strip() or None


def write_manifest(corpus: Corpus, emitter: Emitter, commit: str | None) -> None:
    """Write `manifest.json` and seal the emitter. The build ends here.

    `emitter.checksums` is complete only once every other builder has run, so
    this must be the last thing a build emits. Sealing is what makes that a rule
    rather than an intention — see `Emitter.seal`.
    """
    emitter.write_json(
        "manifest.json",
        {
            "schema_version": SCHEMA_VERSION,
            "source_commit": commit,
            "counts": {
                "assertions": len(corpus.assertions),
                "datasets": len(corpus.datasets),
                "featured": len(corpus.featured),
                "functional": len(corpus.functional),
                "phenotypes": len(corpus.phenotypes),
                "publications": len(corpus.publications),
            },
            # A snapshot, not an ordering. Sorting here was specified and then
            # dropped: `encode_json` passes `sort_keys=True`, which orders keys
            # at every level, so an unsorted mapping encodes to the same bytes.
            # Measured rather than reasoned — restoring the sort changes no
            # published byte and no test can see it — because a sort that buys
            # nothing has already been documented twice in this build as the
            # thing that keeps it reproducible.
            #
            # The copy is equally undetectable and is kept as hygiene, which is
            # the honest description: the payload is encoded before `write_json`
            # records anything, so aliasing `checksums` here would not put
            # `manifest.json` into its own `files` today. It would if that
            # ordering ever changed, and a mapping handed out is one a caller
            # can mutate.
            "files": dict(emitter.checksums),
        },
    )
    emitter.seal()
