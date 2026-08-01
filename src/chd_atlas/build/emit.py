# src/chd_atlas/build/emit.py
"""The only place build output is serialised.

Two builds of the same commit must produce byte-identical output, because the
manifest records a sha256 per file and a checksum that changes on every run
verifies nothing — every consumer caching on it re-downloads the whole site on
each deploy.

What this module removes are the two sources of drift a caller cannot see: dict
iteration order, which follows insertion and so differs between two code paths
that assembled the same content, and gzip's embedded timestamp, which changes
from one second to the next.

Array order is not one of them. `sort_keys` orders dict keys and does nothing to
list elements, so a list that reaches here in a varying order is published in a
varying order. That is the likeliest way this build loses determinism, and the
caller is the only one who can fix it: anything derived from a `set`, from
`polars.group_by` or `DataFrame.unique()` (neither preserves order without
`maintain_order=True`), or from a directory listing must be sorted before it is
handed to `write_json`.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from chd_atlas.fs import write_bytes_atomically

# What `json.dumps` actually accepts, spelled out so a wrong payload is caught at
# the call site rather than raising part-way through a build that has already
# written half of `dist/`. A `date`, a `Decimal`, a `Path`, a `set`, or a pydantic
# model passed instead of its `model_dump()` are all ordinary mistakes here.
#
# `Mapping`/`Sequence` rather than `dict`/`list`: the concrete types are
# invariant, which would reject an ordinary `dict[str, str]` argument. `dict[str,
# Any]` still satisfies this, so `model_dump()` callers are unaffected.
#
# Restricting keys to `str` also closes a trap that `sort_keys` cannot: dumps
# coerces non-string keys, so {9: "b", 10: "a"} sorts to 9, 10 while the already
# stringified {"9": "b", "10": "a"} sorts to 10, 9 — and {True: "a", 1: "b"}
# collapses to a single entry. Mixed key types raise TypeError outright.
#
# One honest hole: `bytes` satisfies `Sequence[int]` and type-checks here even
# though `json.dumps` rejects it. Worth it for everything above.
type Json = str | int | float | bool | None | Mapping[str, Json] | Sequence[Json]


def encode_json(payload: Json) -> bytes:
    """One canonical byte string per payload.

    `sort_keys` because dict order otherwise follows insertion order, so the same
    content assembled by a different code path would checksum differently.
    `ensure_ascii=False` so a term such as "Folie à deux" is stored as itself
    rather than as an escape sequence — the file is served as UTF-8 either way,
    and the readable form is what a `curl` user gets.

    `allow_nan=False` because Python otherwise emits bare `NaN` and `Infinity`
    tokens, which RFC 8259 does not define and `JSON.parse` rejects: a single
    non-finite float silently takes down every page that fetches the artifact,
    while the output stays perfectly deterministic and its checksum verifies.
    `tables.py` refuses these at the input boundary, but build-time arithmetic —
    a ratio, a mean, an aggregate over an empty group — manufactures them from
    data that never passed through table validation. Failing the build loudly is
    the point; do not drop the keyword to make a red build go green.
    """
    text = json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False)
    return (text + "\n").encode("utf-8")


def compress(raw: bytes) -> bytes:
    """Gzip with no embedded timestamp.

    `gzip.compress` writes the current time into the header, so identical content
    compressed a second later differs byte for byte. `mtime=0` removes the only
    non-deterministic field in the format.

    The deflate stream under that header is stable for a given zlib, but is not
    guaranteed across zlib builds: a toolchain upgrade or a switch to zlib-ng
    moves every `.gz` checksum in the manifest at once. That is one deploy's
    worth of cache churn rather than corruption, and it is the expected reason
    for such a change — not a regression in this module.
    """
    buffer = io.BytesIO()
    with gzip.GzipFile(fileobj=buffer, mode="wb", compresslevel=9, mtime=0) as handle:
        handle.write(raw)
    return buffer.getvalue()


def checksum(raw: bytes) -> str:
    """The `sha256:<hex>` identifier the manifest publishes for one artifact."""
    return f"sha256:{hashlib.sha256(raw).hexdigest()}"


@dataclass
class Emitter:
    """Writes every build artifact and remembers what it wrote.

    `checksums` maps each output path to the checksum of the bytes served at it,
    which is what `manifest.py` publishes. Callers must supply those paths
    relative to `root` and separated by forward slashes on every platform,
    because they are URLs to the consumer; nothing here rewrites a path, and
    `_write` rejects the shapes it cannot serve.
    """

    root: Path
    checksums: dict[str, str] = field(default_factory=dict)

    def write_json(self, relative: str, payload: Json) -> None:
        self._write(relative, encode_json(payload))

    def write_json_gz(self, relative: str, payload: Json) -> None:
        self._write(relative, compress(encode_json(payload)))

    def _write(self, relative: str, raw: bytes) -> None:
        segments = relative.split("/")
        # An absolute or dot-bearing path desynchronises the manifest from the
        # disk rather than failing: "/tmp/a.json" splits to a leading "" that
        # `joinpath` discards, so the bytes land at `root/tmp/a.json` while
        # `checksums` advertises "/tmp/a.json" — a published URL that 404s,
        # carrying a checksum matching no file anywhere. ".." escapes `root`
        # outright, and "" reaches `mkstemp` as a directory three frames down.
        # These paths are built from curated identifiers, and `Doi` alone
        # (`^10\.\d{4,9}/\S+$`) admits both slashes and dot segments.
        if any(segment in ("", ".", "..") for segment in segments):
            raise ValueError(
                f"output path must be relative with no empty or dot segments: {relative!r}"
            )
        # Checked before the write, so that the first artifact is the one that
        # survives. Two records normalising to one filename is a curation bug;
        # left alone, one record's page vanishes while the manifest still lists
        # every file and every checksum still verifies. `checksums` is also the
        # only record of what was built, so by the time the manifest runs the
        # overwrite is no longer detectable.
        if relative in self.checksums:
            raise ValueError(f"{relative} written twice; the first artifact would be lost")
        path = self.root.joinpath(*segments)
        path.parent.mkdir(parents=True, exist_ok=True)
        write_bytes_atomically(path, raw)
        # Checksum the serialised artifact, not the payload: for a .gz the two
        # differ, and only the former is what a consumer can verify against what
        # it downloaded. There is no read-back here — `raw` is those exact bytes.
        self.checksums[relative] = checksum(raw)
