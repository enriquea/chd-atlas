# tests/unit/test_emit.py
import gzip
import itertools
import json
import time
from pathlib import Path

import pytest

from chd_atlas.build.emit import Emitter, checksum, compress, encode_json


def test_encoding_keeps_non_ascii_readable() -> None:
    """The one encoding property the known-answer test below cannot pin.

    An escaped payload is still valid JSON of the same value, so `ensure_ascii`
    needs a fixture that actually carries a non-ASCII character; the ASCII golden
    bytes are identical either way.
    """
    raw = encode_json({"name": "Folie à deux"})
    assert "Folie à deux" in raw.decode("utf-8")


def test_encoded_bytes_are_exactly_the_published_form() -> None:
    """A known answer, because the published bytes are a contract.

    Any assertion that compares `encode_json` against itself survives a change to
    `indent`, which would re-checksum every artifact in the atlas and invalidate
    every consumer's cache while the suite stayed green. Spelling the bytes out is
    what makes that change visible.

    Three properties are pinned by this one literal, which is why each no longer
    has a test of its own:

    * `indent`, above;
    * `sort_keys` — the payload is written with its keys in the opposite order to
      the one asserted, and Python preserves insertion order, so a payload
      assembled by a different code path would otherwise checksum differently
      while meaning exactly the same;
    * the trailing newline, which is the last byte below.
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


def test_checksum_is_a_known_sha256() -> None:
    """Which algorithm is used is published in the manifest and verified against.

    Pinned to a known answer rather than to a shape, because a shape assertion —
    the `sha256:` prefix, then 64 hex characters — is equally satisfied by
    `sha512(...)[:64]` and by an `md5` padded out to the same length. The prefix
    and the length are both implied by the literal below, so the separate shape
    test this replaces added no detection at all.
    """
    assert checksum(b"payload") == (
        "sha256:239f59ed55e737c77147cf55ad0c1b030b6d7ee748a7426952f9b852d5a935e5"
    )


def test_emitter_writes_json_and_records_its_checksum(tmp_path: Path) -> None:
    """The happy path, together with the directories it has to create on the way.

    `a/b/c/deep.json` is nested more than one level under `root`, which is what
    separates a recursive `mkdir` from a single-level one: `genes/index.json`
    alone is still written successfully against `parents=False`.
    """
    emitter = Emitter(root=tmp_path)

    emitter.write_json("genes/index.json", {"genes": []})
    emitter.write_json("a/b/c/deep.json", {})

    written = (tmp_path / "genes" / "index.json").read_bytes()
    assert json.loads(written) == {"genes": []}
    assert (tmp_path / "a" / "b" / "c" / "deep.json").is_file()
    assert list(emitter.checksums) == ["genes/index.json", "a/b/c/deep.json"]
    assert emitter.checksums["genes/index.json"] == checksum(written)


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


# Both orderings, because they fail differently. The stored key is casefolded, so
# a guard that looked up the incoming path *without* casefolding it would still
# match when the second path is the all-lowercase one — and would silently miss
# when the second path carries the capitals. Only the second row catches that.
@pytest.mark.parametrize(
    ("first", "second"),
    [
        ("genes/HGNC_11604.json", "genes/hgnc_11604.json"),
        ("genes/hgnc_11604.json", "genes/HGNC_11604.json"),
    ],
)
def test_two_paths_differing_only_in_case_are_refused(
    tmp_path: Path, first: str, second: str
) -> None:
    """A case-insensitive filesystem keeps one file for the two, silently.

    `paths.slug` is injective over strings, but case is what APFS and NTFS
    discard: `HGNC:11604` is a valid `HgncId` and `hgnc_11604` a valid
    `ContrastId`, and they slug to these two paths. Left unguarded on macOS the
    second write replaces the first while `checksums` advertises both, so a
    published URL serves another record's bytes under a checksum that cannot
    verify. The Linux CI would build two files from the same input and never
    reproduce it.
    """
    emitter = Emitter(root=tmp_path)
    emitter.write_json(first, {"record": "first"})

    with pytest.raises(ValueError, match="differ only in case"):
        emitter.write_json(second, {"record": "second"})

    # As with the exact-duplicate guard, the first artifact is the one that
    # survives, so the check must precede the write.
    assert list(emitter.checksums) == [first]
    assert [p.name for p in (tmp_path / "genes").iterdir()] == [first.split("/")[1]]
    assert json.loads((tmp_path / first).read_bytes()) == {"record": "first"}


def test_the_case_guard_covers_paths_the_emitter_was_constructed_with(tmp_path: Path) -> None:
    """The guard has to hold for every path the emitter knows, not just its own writes.

    `checksums` is a constructor argument, so an emitter can start out already
    knowing about files. Seeding `_casefolded` only from `_write` left those
    entries visible to the exact-duplicate guard and invisible to this one —
    weakening precisely the direction that fails silently, since an exact
    duplicate is refused on every filesystem while a case collision is refused
    on none and then loses a file on macOS and Windows.

    The second half is the state the first cannot reach: a mapping handed in
    that already contains the collision. Accepting it would leave an emitter
    whose own invariant is false from the first line, so it is refused at
    construction rather than at whichever later write happens to notice.
    """
    emitter = Emitter(root=tmp_path, checksums={"genes/HGNC_11604.json": "sha256:0"})

    with pytest.raises(ValueError, match="differ only in case"):
        emitter.write_json("genes/hgnc_11604.json", {"record": "second"})
    assert not (tmp_path / "genes").exists()

    with pytest.raises(ValueError, match="differ only in case"):
        Emitter(root=tmp_path, checksums={"a.json": "sha256:0", "A.json": "sha256:1"})


def test_the_case_guard_reports_a_different_failure_than_the_exact_one(tmp_path: Path) -> None:
    """Two distinct bugs: one path built twice, versus two paths one disk merges.

    The first is a curation duplicate and reproduces everywhere; the second only
    loses data on a case-insensitive filesystem. A single shared message would
    send a curator looking for the wrong thing.
    """
    emitter = Emitter(root=tmp_path)
    emitter.write_json("genes/A.json", {})

    with pytest.raises(ValueError, match="written twice") as exact:
        emitter.write_json("genes/A.json", {})
    with pytest.raises(ValueError, match="differ only in case") as folded:
        emitter.write_json("genes/a.json", {})

    # The case message must name both paths; the curator cannot find the other
    # one by grepping for a filename that is not in the error.
    assert "genes/a.json" in str(folded.value)
    assert "genes/A.json" in str(folded.value)
    assert "differ only in case" not in str(exact.value)


# Four rows, each of which is the sole survivor of some way of weakening the
# guard, so none can be dropped without a mutant going unnoticed:
#
#   "/tmp/absolute.json"  the leading empty segment
#   ""                    the whole path empty — reaches `mkstemp` as a directory
#   "a/../b.json"         a dot-dot *inside* the path, which a `startswith("..")`
#                         guard waves through where a leading "../" is caught
#   "./a.json"            a single dot, the only segment kind the other three miss
#
# Two earlier rows are gone because they were strictly weaker: "../escape.json"
# is caught by every mutation "a/../b.json" catches and one fewer, and
# "a//b.json" adds nothing to "/tmp/absolute.json" and "" together.
@pytest.mark.parametrize("relative", ["/tmp/absolute.json", "", "a/../b.json", "./a.json"])
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
