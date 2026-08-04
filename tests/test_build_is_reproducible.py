# tests/test_build_is_reproducible.py
import shutil
from pathlib import Path

import pytest

from chd_atlas.build.runner import build_site

REPO = Path(__file__).parent.parent


@pytest.fixture(scope="module")
def builds(
    tmp_path_factory: pytest.TempPathFactory,
) -> tuple[Path, Path, Path, dict[str, str], dict[str, str]]:
    """The real corpus, built twice into two fresh directories.

    Module-scoped because a build costs about a second and a half, almost all of
    it validation, and both tests below need the same two builds. `ontologies` is
    symlinked rather than copied: 11 MB of pinned OBO releases that nothing here
    edits.

    Asserts nothing. A fixture that checked `first == second` would turn every
    determinism regression into a fixture ERROR on both tests rather than a
    FAILURE on the one that owns the claim, which is a worse report for the
    same information.
    """
    base = tmp_path_factory.mktemp("reproducible")
    root = base / "repo"
    root.mkdir()
    for name in ("curation", "mirrors"):
        shutil.copytree(REPO / name, root / name)
    (root / "ontologies").symlink_to(REPO / "ontologies")

    first_out, second_out = base / "first", base / "second"
    first = build_site(root, first_out)
    second = build_site(root, second_out)
    return root, first_out, second_out, first, second


def test_two_builds_of_one_commit_are_byte_identical(
    builds: tuple[Path, Path, Path, dict[str, str], dict[str, str]],
) -> None:
    """The manifest publishes a sha256 per file, and that is a promise.

    A checksum that changes between two builds of one commit verifies nothing,
    and every consumer caching on it re-downloads the whole site on each deploy.

    **What this does not guard, so nobody deletes the guards that do.** It cannot
    detect hash-order nondeterminism — an array built from a `set` and published
    unsorted — because `PYTHONHASHSEED` is fixed for the life of an interpreter,
    so two builds in one process iterate an unsorted set identically. That is
    amendment A33/A34, and it holds at any corpus size: this test builds twice in
    one process, and nothing about the corpus can change that.

    A second argument used to stand beside it and no longer does. It read that
    the committed corpus has one gene, so a one-element set has one iteration
    order under every seed, and that dropping the gene sort in `search.py` and
    building under `PYTHONHASHSEED` 0 and 12345 produced byte-identical output.
    That was measured against a build whose search index was keyed on the
    asserted genes. D31 rekeyed it on `published`, which is 23 genes, and the
    same measurement rerun on 2026-08-04 now gives two different files —
    `search/index.json.gz` and, through it, `manifest.json`. The corpus-size
    escape hatch is gone; the in-process one above is the whole reason this test
    cannot see hash order, which is why it is stated first.

    Hash order is guarded where a fixture can be sized to expose it —
    `test_build_search.py`, whose six-gene fixture kills that mutant on 30/30
    seeds, and likewise `test_derive.py` and `test_build_omics.py`. What only a
    whole build can show is a timestamp anywhere in the output, and gzip's
    embedded mtime; both are caught here, because the builds are seconds apart
    and gzip's MTIME has one-second resolution. Verified by mutating `mtime=0`
    to `mtime=None` in a fresh interpreter, which fails this test — the first
    attempt at that measurement edited the file after the module was already
    imported and proved nothing.
    """
    _, first_out, second_out, first, second = builds

    assert first == second, "the two builds disagree about what they wrote"
    for relative in first:
        assert (first_out / relative).read_bytes() == (second_out / relative).read_bytes(), (
            f"{relative} differs between two builds of one commit"
        )


def test_rebuilding_in_place_is_stable_but_leaves_a_retired_artifact(
    builds: tuple[Path, Path, Path, dict[str, str], dict[str, str]],
) -> None:
    """Two claims about building on top of a previous build: one good, one not.

    CI builds into a fresh checkout, but a curator rebuilds in place, and the
    checksums must not move when they do.

    The second claim is a real limitation, pinned because it is a limitation and
    not because it is wanted. `build_site` writes files; it does not own the
    directory and never deletes from it. So a gene retired from curation keeps
    its bundle at the URL it always had — on disk, absent from the new manifest,
    and still served by any host that lists the directory rather than reading
    the manifest.

    Not fixed by clearing `out` first: `--out` is a caller-supplied path, and a
    build command that empties a directory it was merely pointed at is a worse
    failure than a stale page. The deploy workflow is where this is handled, by
    publishing the built directory as a fresh artifact rather than syncing into
    a previous one — and this test is what makes that a requirement on Task 14
    rather than an assumption it happens to satisfy.
    """
    root, first_out, _, first, _second = builds
    retired = first_out / "genes" / "HGNC_99999.json"
    retired.write_text('{"gene": "retired between builds"}\n')

    again = build_site(root, first_out)

    assert again == first
    assert retired.is_file()
    assert "genes/HGNC_99999.json" not in again
