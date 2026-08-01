# src/chd_atlas/build/emit.py
"""The only place build output is serialised.

Two builds of the same commit must produce byte-identical output, because the
manifest records a sha256 per file and a checksum that changes on every run
verifies nothing. Exactly two things break that — dictionary ordering and gzip's
embedded timestamp — and both are handled here so no caller has to remember.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from chd_atlas.fs import write_bytes_atomically


def encode_json(payload: Any) -> bytes:
    """One canonical byte string per payload.

    `sort_keys` because dict order otherwise follows insertion order, so the same
    content assembled by a different code path would checksum differently.
    `ensure_ascii=False` so a term such as "Folie à deux" is stored as itself
    rather than as an escape sequence — the file is served as UTF-8 either way,
    and the readable form is what a `curl` user gets.
    """
    text = json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    return text.encode("utf-8")


def compress(raw: bytes) -> bytes:
    """Gzip with no embedded timestamp.

    `gzip.compress` writes the current time into the header, so identical content
    compressed a second later differs byte for byte. `mtime=0` removes the only
    non-deterministic field in the format.
    """
    buffer = io.BytesIO()
    with gzip.GzipFile(fileobj=buffer, mode="wb", compresslevel=9, mtime=0) as handle:
        handle.write(raw)
    return buffer.getvalue()


def digest(raw: bytes) -> str:
    return f"sha256:{hashlib.sha256(raw).hexdigest()}"


@dataclass
class Emitter:
    """Writes every build artifact and remembers what it wrote.

    `files` maps each repo-relative output path to the checksum of the bytes as
    served, which is what `manifest.py` publishes. Paths are stored with forward
    slashes regardless of platform, because they are URLs to the consumer.
    """

    root: Path
    files: dict[str, str] = field(default_factory=dict)

    def write_json(self, relative: str, payload: Any) -> None:
        self._write(relative, encode_json(payload))

    def write_json_gz(self, relative: str, payload: Any) -> None:
        self._write(relative, compress(encode_json(payload)))

    def _write(self, relative: str, raw: bytes) -> None:
        path = self.root.joinpath(*relative.split("/"))
        path.parent.mkdir(parents=True, exist_ok=True)
        write_bytes_atomically(path, raw)
        # Checksum the bytes on disk, not the payload: for a .gz artifact those
        # differ, and only the former is what a consumer can verify.
        self.files[relative] = digest(raw)
