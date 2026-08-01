# tests/unit/test_emit.py
import gzip
import itertools
import json
import time
from pathlib import Path

import pytest

from chd_atlas.build.emit import Emitter, checksum, compress, encode_json


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


def test_encoded_bytes_are_exactly_the_published_form() -> None:
    """A known answer, because the published bytes are a contract.

    Every other assertion here compares `encode_json` against itself and so
    survives a change to `indent`, which would re-checksum every artifact in the
    atlas and invalidate every consumer's cache while the suite stayed green.
    Spelling the bytes out is what makes that change visible.
    """
    assert encode_json({"b": [1, 2], "a": "x"}) == (
        b'{\n  "a": "x",\n  "b": [\n    1,\n    2\n  ]\n}\n'
    )


@pytest.mark.parametrize("value", [float("nan"), float("inf"), float("-inf")])
def test_non_finite_floats_are_refused(value: float) -> None:
    """`NaN` and `Infinity` are not JSON; `JSON.parse` throws on them.

    Python's own loader accepts these tokens, so a round-trip through `json.loads`
    reports success while every browser fetching the artifact fails. The output
    would be perfectly deterministic and its checksum would verify — the page
    would simply be broken. `tables.py` rejects non-finite floats at the input
    boundary, but a build-time ratio or an aggregate over an empty group can
    produce one that never passed through table validation.
    """
    with pytest.raises(ValueError):
        encode_json({"af": value})


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


def test_compression_level_is_pinned() -> None:
    """Changing the level re-checksums every `.gz` in the atlas at once.

    Asserted through the header's XFL byte, which gzip sets to 2 for best
    compression and 0 for anything between fast and best, rather than through the
    compressed bytes: the deflate stream is stable only for a given zlib, so a
    golden copy of it would fail on a toolchain upgrade or a switch to zlib-ng —
    cache churn rather than a defect. The header byte is written by Python and
    does not move.
    """
    assert compress(b"payload")[8] == 2


def test_checksum_is_prefixed_and_hex() -> None:
    value = checksum(b"payload")
    assert value.startswith("sha256:")
    assert len(value) == len("sha256:") + 64


def test_checksum_is_a_known_sha256() -> None:
    """The shape assertion above is also satisfied by `sha512(...)[:64]`.

    Which algorithm is used is published in the manifest and verified against by
    consumers, so it is pinned to a known answer rather than to a length.
    """
    assert checksum(b"payload") == (
        "sha256:239f59ed55e737c77147cf55ad0c1b030b6d7ee748a7426952f9b852d5a935e5"
    )


def test_emitter_writes_json_and_records_its_checksum(tmp_path: Path) -> None:
    emitter = Emitter(root=tmp_path)

    emitter.write_json("genes/index.json", {"genes": []})

    written = (tmp_path / "genes" / "index.json").read_bytes()
    assert json.loads(written) == {"genes": []}
    assert emitter.checksums == {"genes/index.json": checksum(written)}


def test_emitter_checksums_the_bytes_actually_served(tmp_path: Path) -> None:
    """For a .gz file that is the compressed bytes, not the JSON inside.

    The checksum exists so a consumer can verify what it downloaded; checksumming
    the pre-compression payload would make it unverifiable against the file.
    """
    emitter = Emitter(root=tmp_path)

    emitter.write_json_gz("variants/1.json.gz", {"variants": []})

    written = (tmp_path / "variants" / "1.json.gz").read_bytes()
    assert gzip.decompress(written) == encode_json({"variants": []})
    assert emitter.checksums == {"variants/1.json.gz": checksum(written)}


def test_emitter_creates_missing_parent_directories(tmp_path: Path) -> None:
    emitter = Emitter(root=tmp_path)

    emitter.write_json("a/b/c/deep.json", {})

    assert (tmp_path / "a" / "b" / "c" / "deep.json").is_file()


def test_writing_one_path_twice_is_refused(tmp_path: Path) -> None:
    """A collision means two records normalised to one filename — a curation bug.

    Overwriting loses one record's page while the manifest still lists every file
    and every checksum still verifies, and `checksums` is the only record of what
    was built, so nothing downstream can notice.
    """
    emitter = Emitter(root=tmp_path)
    emitter.write_json("genes/index.json", {"genes": ["first"]})

    with pytest.raises(ValueError, match="written twice"):
        emitter.write_json("genes/index.json", {"genes": ["second"]})

    # The first artifact is the one that must survive, so the check has to happen
    # before the write rather than after it.
    written = (tmp_path / "genes" / "index.json").read_bytes()
    assert json.loads(written) == {"genes": ["first"]}


@pytest.mark.parametrize(
    "relative",
    ["/tmp/absolute.json", "../escape.json", "a/../b.json", "", "a//b.json", "./a.json"],
)
def test_paths_that_would_desynchronise_the_manifest_are_refused(
    tmp_path: Path, relative: str
) -> None:
    """An absolute path publishes a URL that does not exist.

    `joinpath` discards the empty leading segment of "/tmp/absolute.json", so the
    file lands under `root` while the manifest advertises the absolute path — a
    404 carrying a checksum that matches nothing. ".." writes outside `root`
    entirely. These paths are built from curated identifiers, so the shapes are
    reachable from data rather than only from a typo.
    """
    emitter = Emitter(root=tmp_path)

    with pytest.raises(ValueError, match="must be relative"):
        emitter.write_json(relative, {})

    assert emitter.checksums == {}
    assert list(tmp_path.iterdir()) == []
