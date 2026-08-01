# tests/unit/test_emit.py
import gzip
import itertools
import json
import time
from pathlib import Path

import pytest

from chd_atlas.build.emit import Emitter, compress, digest, encode_json


def test_encoding_is_stable_regardless_of_key_insertion_order() -> None:
    """Two dicts with the same content must encode to the same bytes.

    Python preserves insertion order, so a payload assembled by a different code
    path would otherwise checksum differently while meaning exactly the same.
    """
    first = encode_json({"b": 1, "a": 2})
    second = encode_json({"a": 2, "b": 1})

    assert first == second


def test_encoding_keeps_non_ascii_readable() -> None:
    raw = encode_json({"name": "Folie à deux"})
    assert "Folie à deux" in raw.decode("utf-8")


def test_encoding_ends_with_a_newline() -> None:
    assert encode_json({"a": 1}).endswith(b"\n")


def test_compression_embeds_no_timestamp(monkeypatch: pytest.MonkeyPatch) -> None:
    """gzip writes the current time into its header by default.

    Left alone, two builds of identical content differ byte for byte, every
    manifest checksum changes on every run, and the determinism test cannot hold.

    The clock is advanced between the two calls, and that is what gives this test
    teeth. gzip stores its timestamp as whole seconds, so two real calls a
    microsecond apart write the same header: against a compressor that *does*
    embed the clock, the two payloads still come out identical and the assertion
    below passes. Removing the fake clock leaves a test that cannot fail.
    """
    ticks = itertools.count(1_700_000_000, 3600)
    monkeypatch.setattr(time, "time", lambda: float(next(ticks)))

    first = compress(b"payload")
    second = compress(b"payload")

    assert first == second
    assert gzip.decompress(first) == b"payload"


def test_digest_is_prefixed_and_hex() -> None:
    value = digest(b"payload")
    assert value.startswith("sha256:")
    assert len(value) == len("sha256:") + 64


def test_emitter_writes_json_and_records_its_checksum(tmp_path: Path) -> None:
    emitter = Emitter(root=tmp_path)

    emitter.write_json("genes/index.json", {"genes": []})

    written = (tmp_path / "genes" / "index.json").read_bytes()
    assert json.loads(written) == {"genes": []}
    assert emitter.files == {"genes/index.json": digest(written)}


def test_emitter_checksums_the_bytes_actually_served(tmp_path: Path) -> None:
    """For a .gz file that is the compressed bytes, not the JSON inside.

    The checksum exists so a consumer can verify what it downloaded; checksumming
    the pre-compression payload would make it unverifiable against the file.
    """
    emitter = Emitter(root=tmp_path)

    emitter.write_json_gz("variants/1.json.gz", {"variants": []})

    written = (tmp_path / "variants" / "1.json.gz").read_bytes()
    assert gzip.decompress(written) == encode_json({"variants": []})
    assert emitter.files == {"variants/1.json.gz": digest(written)}


def test_emitter_creates_missing_parent_directories(tmp_path: Path) -> None:
    emitter = Emitter(root=tmp_path)

    emitter.write_json("a/b/c/deep.json", {})

    assert (tmp_path / "a" / "b" / "c" / "deep.json").is_file()
