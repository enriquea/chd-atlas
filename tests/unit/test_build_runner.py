# tests/unit/test_build_runner.py
import gzip
import hashlib
import json
import shutil
import subprocess
from pathlib import Path
from typing import Any

import pytest

from chd_atlas.build.runner import BuildRefused, _gene_registry, build_site
from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.validate.runner import ValidationReport

REPO = Path(__file__).parent.parent.parent

# Every artifact the published API promises. Listed rather than derived, because
# this is the contract: a builder dropped from `build_site` would otherwise take
# its artifact out of the site and out of the assertion in the same edit.
PUBLISHED = (
    "manifest.json",
    "index.html",
    "genes/index.json",
    "publications.json",
    "featured.json",
    "phenotypes.json",
    "datasets.json",
    "sources.json",
    "variants/index.json",
    "search/index.json.gz",
)


@pytest.fixture
def repo(tmp_path: Path) -> Path:
    """A copy of the real committed corpus, which validates cleanly.

    `curation` and `mirrors` are copied because tests edit them; `ontologies` is
    symlinked because it is 11 MB of pinned OBO releases that nothing here
    changes, and copying it per test would cost more than the builds do.
    """
    source = tmp_path / "repo"
    source.mkdir()
    for name in ("curation", "mirrors"):
        shutil.copytree(REPO / name, source / name)
    (source / "ontologies").symlink_to(REPO / "ontologies")
    return source


def _break(repo: Path) -> None:
    """Make the corpus fail validation with one SCHEMA001.

    An empty `assertions` list, because `AssertionFile` requires at least one
    record — a file that parses as YAML and fails as a record, so the failure is
    the validator's rather than the loader's.
    """
    (repo / "curation" / "assertions" / "BROKEN.yaml").write_text("assertions: []\n")


def _search_records(out: Path) -> list[dict[str, Any]]:
    raw = (out / "search" / "index.json.gz").read_bytes()
    records: list[dict[str, Any]] = json.loads(gzip.decompress(raw))["records"]
    return records


def test_the_registry_reader_cleans_the_alias_cell_and_skips_an_unlabelled_row(
    tmp_path: Path,
) -> None:
    """The two ways a mirror cell reaches the site as something a curator did not write.

    Read directly rather than through `build_site`, because both inputs are what
    the gate exists to prevent and one of them the gate does not catch.

    Whitespace and empty segments are reachable on a fully green build: nothing
    validates the contents of the `aliases` cell, so `T-box 5|  TBX5B  ||Chr12`
    is a legal mirror row, and unstripped it publishes a search term with
    leading spaces plus an empty term that every "contains" query matches.

    A null `symbol` is not reachable behind `build_site` — it is non-nullable, so
    TBL003 refuses the build first — and the row is skipped rather than
    published because `str(None)` is the string "None", which `build_genes` would
    render as the gene's symbol in the browse row and as the page heading. The
    gene keeps its HGNC id as its label instead, which is what `search.py`'s own
    fallback does with the same reasoning. Same class as A12's null protein
    matching an unrelated gene.
    """
    mirror = tmp_path / "mirrors"
    mirror.mkdir()
    (mirror / "genes.tsv").write_text(
        "hgnc_id\tsymbol\tname\taliases\tensembl_gene\tncbi_gene\tlocus\tuniprot\tmane_select\n"
        "HGNC:11604\tTBX5\tT-box transcription factor 5\tT-box 5|  TBX5B  ||Chr12\t\t\t\t\t\n"
        "HGNC:4173\t\tGATA binding protein 4\t\t\t\t\t\t\n"
        # A single space, which TBL003 does not catch: `read_table` maps only the
        # empty string to null. Published, it is a browse row and a search result
        # that render as nothing.
        "HGNC:9999\t \tSome gene\t\t\t\t\t\t\n"
    )

    registry = _gene_registry(tmp_path)

    assert list(registry) == ["HGNC:11604"]
    assert registry["HGNC:11604"].aliases == ("T-box 5", "TBX5B", "Chr12")
    assert registry["HGNC:11604"].name == "T-box transcription factor 5"


def test_the_build_publishes_every_artifact_and_a_manifest_that_verifies(
    repo: Path, tmp_path: Path
) -> None:
    """The whole contract in one build, checked against the bytes on disk.

    Every checksum is verified rather than merely shaped, because that is the
    manifest's entire purpose and A43 records what the plan's suite missed: with
    every value replaced by 64 zeros, tests asserting only the `sha256:` prefix
    all passed. A manifest whose checksums are wrong is worse than none — a
    consumer verifying a download rejects a file that is correct.

    The file set is checked in both directions. A file on disk that the manifest
    omits is served with no checksum and invisible to any consumer that
    verifies; a manifest entry with no file is a published URL that 404s. Only
    `manifest.json` may be in one and not the other, since it cannot contain its
    own checksum.
    """
    out = tmp_path / "dist"

    written = build_site(repo, out)

    for relative in PUBLISHED:
        assert relative in written, f"{relative} is promised by the API and was not built"

    manifest = json.loads((out / "manifest.json").read_text())
    for relative, advertised in manifest["files"].items():
        actual = hashlib.sha256((out / relative).read_bytes()).hexdigest()
        assert advertised == f"sha256:{actual}", f"{relative} does not match its checksum"

    on_disk = {path.relative_to(out).as_posix() for path in out.rglob("*") if path.is_file()}
    assert on_disk == set(manifest["files"]) | {"manifest.json"}

    # The return value accounts for the manifest as well as everything in it, so
    # a caller checking what was built does not have to know that one artifact
    # is missing from the manifest by construction.
    assert set(written) == on_disk


def test_a_build_from_a_checkout_stamps_the_manifest_with_that_commit(
    repo: Path, tmp_path: Path
) -> None:
    """The provenance branch that only ever runs in the deploy, and nowhere else.

    Every other test here builds from a `tmp_path` copy that is not a git
    repository, so `source_commit` returns `None` in all of them and the manifest
    records no commit. That is the opposite of what happens in CI, where the
    build always runs inside a checkout — so the wiring the published site
    actually depends on was reachable by no test, and a regression in how
    `build_site` obtains the commit would have shipped green.

    Made a checkout here rather than asserting against this repository's own
    HEAD, which would pass just as well against a `source_commit` that ignored
    `root` entirely.
    """
    subprocess.run(["git", "init", "-q"], cwd=repo, check=True)
    subprocess.run(["git", "add", "-A"], cwd=repo, check=True)
    subprocess.run(
        ["git", "-c", "user.email=t@t", "-c", "user.name=t", "commit", "-qm", "corpus"],
        cwd=repo,
        check=True,
    )
    head = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=repo, capture_output=True, text=True, check=True
    ).stdout.strip()
    out = tmp_path / "dist"

    build_site(repo, out)

    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["source_commit"] == head
    assert len(head) >= 40


def test_a_repository_that_fails_validation_is_refused_before_anything_is_written(
    repo: Path, tmp_path: Path
) -> None:
    """Publishing from a corpus that fails its own gate defeats having one.

    Two halves, and the second is the one that matters on a deploy: the failure
    must name the validation errors, so a curator who ran `build` knows what to
    fix without being told to go and run a second command; and nothing may be
    written, not even an empty directory, because a half-built `dist/` is
    something a deploy step will happily upload.
    """
    _break(repo)
    out = tmp_path / "dist"

    with pytest.raises(BuildRefused) as raised:
        build_site(repo, out)

    assert "SCHEMA001" in str(raised.value)
    assert not out.exists()


@pytest.mark.parametrize(
    ("code", "severity", "refuses"),
    [("SCHEMA001", Severity.ERROR, True), ("REF000", Severity.WARNING, False)],
)
def test_the_gate_refuses_on_an_error_and_builds_through_a_warning(
    repo: Path,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    code: str,
    severity: Severity,
    refuses: bool,
) -> None:
    """Which condition the gate reads, decided here and not reimplemented.

    A stricter gate makes a repository CI calls green unpublishable; a looser one
    publishes a corpus the validator rejected. `ValidationReport.ok` ignores
    warnings by design, and that half cannot be provoked from a fixture: every
    warning this project emits is a "skipped checks" notice that fires *because*
    something else errored, so no corpus produces warnings alone. Leaving it
    untested is how a gate tightened to `warning_count == 0` reaches `main`.

    The report is therefore supplied rather than provoked, which makes this test
    about the condition alone. That the real validator's verdict is what reaches
    the gate is a different claim, pinned by the two tests above: a genuinely
    broken corpus refuses and names its SCHEMA001, and a genuinely clean one
    builds. Running the real validator again here cost 3 of the file's 11
    seconds and duplicated both.
    """
    report = ValidationReport(issues=[ValidationIssue(code, severity, str(repo), "constructed")])
    assert report.ok is not refuses
    monkeypatch.setattr("chd_atlas.build.runner.validate_repository", lambda _: report)
    out = tmp_path / "dist"

    if refuses:
        with pytest.raises(BuildRefused, match=code):
            build_site(repo, out)
        assert not out.exists()
    else:
        build_site(repo, out)
        assert (out / "manifest.json").is_file()


def test_every_gene_label_the_registry_holds_reaches_the_site(repo: Path, tmp_path: Path) -> None:
    """Symbol, approved name and aliases, from `mirrors/genes.tsv` to both payloads.

    The committed mirror has a null `aliases` cell, so the split is exercised
    only if this test supplies one. That matters more than it looks: the cell is
    pipe-separated and handing it over unsplit type-checks — `frame.to_dicts()`
    yields `dict[str, Any]`, and an `Any` satisfies `tuple[str, ...]` silently —
    so the mistake this asserts against is one mypy cannot catch. A reader that
    skipped the split publishes one term reading "T-box 5|Chr12q24.1", and a
    reader that passed the raw string publishes a term per character.

    The name is asserted because it is the column every valid mirror is required
    to have (TBL003) and the one A38 found unreachable: a visitor typing
    "T-box transcription factor" must find TBX5.
    """
    mirror = repo / "mirrors" / "genes.tsv"
    header, row = mirror.read_text().splitlines()[:2]
    columns = header.split("\t")
    cells = row.split("\t")
    cells[columns.index("aliases")] = "T-box 5|Chr12q24.1"
    mirror.write_text("\t".join(columns) + "\n" + "\t".join(cells) + "\n")
    out = tmp_path / "dist"

    build_site(repo, out)

    index = json.loads((out / "genes" / "index.json").read_text())
    assert [entry["symbol"] for entry in index["genes"]] == ["TBX5"]

    gene = next(record for record in _search_records(out) if record["kind"] == "gene")
    assert gene["label"] == "TBX5"
    assert gene["terms"] == [
        "TBX5",
        "HGNC:11604",
        "T-box transcription factor 5",
        "Chr12q24.1",
        "T-box 5",
    ]


def test_the_site_carries_the_terms_of_everything_it_republishes(
    repo: Path, tmp_path: Path
) -> None:
    """Attribution has to travel with the data, not sit in a README.

    `mirrors/sources.yaml` records HPO as `permitted_with_attribution`, and the
    build republishes its term labels and synonyms verbatim in `phenotypes.json`
    and in the search index. A consumer fetching JSON never opens the README, and
    the repository's only `LICENSE` is Apache-2.0 — which a reader will
    reasonably take to cover the published files too.

    So the obligation is asserted where it is discharged: every source the
    registry knows appears in the published payload, with its licence and its
    redistribution terms intact.
    """
    out = tmp_path / "dist"

    build_site(repo, out)

    published = json.loads((out / "sources.json").read_text())["sources"]
    by_id = {source["id"]: source for source in published}
    assert "hpo" in by_id, "the ontology whose content is republished must be named"
    assert by_id["hpo"]["redistribution"] == "permitted_with_attribution"
    assert by_id["hpo"]["licence"], "a source with terms must carry them"
    for source in published:
        assert source["name"] and source["url"] and source["licence"]
    # Ordered, because this is a published array like every other.
    assert [source["id"] for source in published] == sorted(by_id)
