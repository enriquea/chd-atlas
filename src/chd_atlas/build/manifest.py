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

import os
import subprocess
from pathlib import Path
from typing import Final

from chd_atlas.build.emit import Emitter
from chd_atlas.corpus import Corpus

# What a consumer can rely on across two builds, in the usual major.minor sense
# — spelled out here because the field was published for three releases without
# anyone saying what it meant.
#
# MINOR rises when a field is added and nothing is removed or repurposed: a
# consumer written against the previous version keeps working, and one that
# wants the new field can test for the version rather than for the key.
# MAJOR rises when a field changes shape or leaves, which breaks such a reader.
#
# 1.1 added `genes` to every omics shard row (issue #3) and
# `conflicting_lesion_groups` to every gene index row (issue #4). Both are
# additive, so 1.0 readers are unaffected.
#
# 2.0 removed `classification` and `source_tier` from the curated assertion and
# moved gene-disease validity to mirrored, attributed records (design decision
# D12: the atlas no longer authors a validity call of its own). A 1.x reader
# looking for a classification on the assertion finds none, which is a
# breaking change and so MAJOR rather than MINOR, even though the release's
# other change — the gene bundle's new `validity` object, carrying every
# mirrored record and its provenance — is purely additive on its own.
#
# 2.1 added `status`, so a programmatic consumer can read the atlas's own
# readiness without scraping `index.html`'s prose for it. Nothing existing
# changed shape or left, so this is MINOR despite following a MAJOR release —
# each release's letter is decided against what it changes, never against the
# rung before it.
SCHEMA_VERSION: Final = "2.1"

# What `status` publishes today. A literal rather than something derived from
# the corpus, unlike every field in `counts`: there is no measurement of "is
# this atlas ready for clinical use" to compute, only a decision to record, and
# recording it as a constant is what stops a future release forgetting to flip
# it — the two-word status and the paragraph on `index.html` explaining what it
# concretely means are edited in the same commit, or neither is trustworthy.
STATUS: Final = "in-development"


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

    The toplevel comparison is what makes "not a checkout" mean it. `git -C X`
    searches *upward* from X, so a `--root` that is merely unpacked inside some
    unrelated repository resolves that repository's HEAD, and this published a
    real, verifiable sha naming a commit containing none of the built data —
    strictly worse than the `"HEAD"` the returncode check exists to prevent,
    because a false provenance claim that resolves is one nobody thinks to
    check. Two environment variables do the same thing more quietly: `GIT_DIR`
    and `GIT_WORK_TREE` override `-C` outright, so they are cleared rather than
    trusted. Reproduced before fixing, and pinned by
    `test_a_root_inside_an_unrelated_checkout_has_no_provenance`.
    """
    # `-C` still points git at `root`; `--show-toplevel` then reports which
    # repository it actually found, which is the answer that has to be checked.
    environment = {key: value for key, value in os.environ.items() if not key.startswith("GIT_")}
    try:
        located = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "--show-toplevel", "HEAD"],
            capture_output=True,
            text=True,
            check=False,
            timeout=10,
            env=environment,
        )
    # OSError is git absent from PATH, which a build must survive for the same
    # reason a tarball build must.
    except (OSError, subprocess.SubprocessError):
        return None
    if located.returncode != 0:
        return None

    # Line-based, not whitespace-based. `--show-toplevel` prints a path, and a
    # path may contain spaces — `~/My Projects/atlas` is ordinary on macOS.
    # Splitting on whitespace made that three tokens instead of two and returned
    # `None`, so a perfectly valid checkout published no provenance and did it
    # by failing into the same answer that means "not a checkout at all".
    # Reported by review on #7 and reproduced before fixing.
    #
    # The commit is the last line and never contains whitespace; everything
    # before it is the path, rejoined, so even a newline in a directory name
    # compares correctly rather than silently answering `None`.
    lines = located.stdout.splitlines()
    if len(lines) < 2:
        return None
    toplevel, commit = "\n".join(lines[:-1]), lines[-1]
    try:
        if Path(toplevel).resolve() != root.resolve():
            return None
    except OSError:
        return None
    # `or None`, because "" is falsy but still reads as a recorded commit to a
    # consumer testing the key against null.
    return commit.strip() or None


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
            "status": STATUS,
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
