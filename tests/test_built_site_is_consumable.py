# tests/test_built_site_is_consumable.py
import gzip
import json
import shutil
from pathlib import Path
from typing import Any

import pytest

from chd_atlas.build.render import RESEARCH_USE_NOTICE
from chd_atlas.build.runner import build_site

REPO = Path(__file__).parent.parent


@pytest.fixture(scope="module")
def site(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """One build of the real committed corpus, shared by every test here.

    The unit suites build from fixtures they control. This file exists to ask a
    different question — whether the data actually committed produces a site a
    frontend can consume — so it must build the real thing, and once.
    """
    base = tmp_path_factory.mktemp("consumable")
    root = base / "repo"
    root.mkdir()
    for name in ("curation", "mirrors"):
        shutil.copytree(REPO / name, root / name)
    (root / "ontologies").symlink_to(REPO / "ontologies")
    out = base / "dist"
    build_site(root, out)
    return out


def _read(site: Path, relative: str) -> dict[str, Any]:
    payload: dict[str, Any] = json.loads(site.joinpath(*relative.split("/")).read_text())
    return payload


def _search_records(site: Path) -> list[dict[str, Any]]:
    raw = (site / "search" / "index.json.gz").read_bytes()
    records: list[dict[str, Any]] = json.loads(gzip.decompress(raw))["records"]
    return records


def _advertised(site: Path) -> list[tuple[str, str]]:
    """Every path any payload points a consumer at, paired with who advertised it."""
    references: list[tuple[str, str]] = []

    for relative in _read(site, "manifest.json")["files"]:
        references.append(("manifest.json", relative))

    index = _read(site, "genes/index.json")
    for entry in index["genes"]:
        references.append(("genes/index.json", entry["bundle"]))
        bundle = _read(site, entry["bundle"])
        for modality in bundle["omics"].values():
            for shard in modality["shards"]:
                references.append((entry["bundle"], shard))

    for shard in _read(site, "variants/index.json")["shards"]:
        references.append(("variants/index.json", shard))

    for record in _search_records(site):
        references.append(("search/index.json.gz", record["path"]))

    return references


def test_every_path_the_site_advertises_resolves(site: Path) -> None:
    """The dead-link sweep, across every payload that points at another.

    A dangling path here is the worst kind of failure this project has: the
    build is green, every checksum verifies, and a reader sees a gene with no
    data rather than a broken site. Nothing downstream can catch it —
    `Emitter.checksums` records what was written, never what a payload promised.

    This is the only test that checks the search index's `path` values against
    what was actually built. That matters because `search.py` names
    `publications.json` and `phenotypes.json` as string literals while
    `literature.py` is what writes them, with nothing tying the two together —
    amendment A36, whose structural fix is shared constants and is queued. Until
    then, this is the detection.

    The floors are asserted so the sweep cannot pass by finding nothing, which
    is how a loop over an empty collection reports success. They are floors
    rather than exact counts so that curating more data strengthens this test
    instead of breaking it. Two of the five sources contribute nothing today —
    no gene carries omics summaries and no variant shard exists in the committed
    corpus — so those loops are currently unexercised, and the floors are
    deliberately not raised to hide that.
    """
    references = _advertised(site)
    by_source: dict[str, int] = {}
    for source, _ in references:
        by_source[source] = by_source.get(source, 0) + 1

    assert by_source.get("manifest.json", 0) >= 8
    assert by_source.get("genes/index.json", 0) >= 1
    assert by_source.get("search/index.json.gz", 0) >= 3

    for source, relative in references:
        assert site.joinpath(*relative.split("/")).is_file(), (
            f"{source} advertises {relative}, which was never written"
        )


def test_every_html_page_the_build_writes_carries_the_research_use_notice(site: Path) -> None:
    """Not one page kind. Every page, counted off the filesystem.

    `test_build_render.py` pins the sentence and pins that `document` emits it.
    This is the different question, and the one the defect actually was: *does
    every page the real build writes go through `document`?* Measured 2026-08-04
    on the build immediately before this change, the answer was no in effect --
    the notice lived in `landing.py`, so `grep -l` matched 1 of 25 files and the
    24 carrying gene-level classifications matched none.

    Enumerated with `rglob("*.html")` rather than against a list of expected
    paths, so a page kind added later is covered without anyone remembering to
    add it here -- which is the failure mode this replaces. The count is asserted
    too: an `rglob` that matched nothing would satisfy an all-pages loop
    vacuously, and 25 is what the committed corpus publishes (1 landing + 1
    browse + 23 gene pages).
    """
    pages = sorted(site.rglob("*.html"))

    assert len(pages) == 25, f"expected 25 HTML pages, found {len(pages)}"
    for page in pages:
        assert RESEARCH_USE_NOTICE in page.read_text(encoding="utf-8"), (
            f"{page.relative_to(site)} carries no research-use notice"
        )


def test_no_advertised_path_would_leave_the_site_root(site: Path) -> None:
    """Every path in the output is a URL relative to the site root.

    An absolute path leaking from a build machine resolves on the machine that
    built it and 404s for every visitor, which is a failure only a reader ever
    sees. A `..` segment escapes the published tree entirely, and a backslash is
    a path separator on the filesystem a curator may be using but a literal
    character in a URL.

    Checked over what the payloads advertise rather than over the emitter's
    guard, because the two are different claims: `Emitter._write` refuses these
    shapes for files it writes, and this asserts that no payload *points* at one.
    """
    for source, relative in _advertised(site):
        assert not relative.startswith("/"), f"{source} advertises an absolute path {relative}"
        assert ".." not in relative.split("/"), f"{source} advertises an escaping path {relative}"
        assert "\\" not in relative, f"{source} advertises a backslash in {relative}"
        assert "://" not in relative, f"{source} advertises an off-site URL {relative}"


def test_the_seed_gene_and_its_manuscript_render_from_the_payloads_alone(site: Path) -> None:
    """What a frontend needs on first paint, from the data actually committed.

    Three separate promises, checked against the real corpus rather than a
    fixture: a browse row carries what it ranks on, a gene page is one fetch,
    and the landing page renders a featured manuscript without a second request.

    `featured[].publication` is asserted to be a resolved object, not a bare
    PMID. Amendment A24 records why: the specified fallback made it an object
    for resolvable entries and a string otherwise, so a consumer reading
    `entry.publication.title` got `undefined` with nothing raised — in the one
    payload the landing page renders without a further fetch.

    `headline_confidence` and `validity_state` read `"definitive"` and
    `"expert_curated"` here, not `null`/`"uncurated"`: the committed ClinGen
    mirror curates TBX5 Definitive for Holt-Oram syndrome
    (`MONDO:0007732`), `build_site` reads that mirror, and a gene the mirror
    plainly classifies must not publish as unassessed.
    """
    entry = next(
        item for item in _read(site, "genes/index.json")["genes"] if item["gene"] == "HGNC:11604"
    )
    assert entry["symbol"] == "TBX5"
    assert entry["headline_confidence"] == "definitive"
    assert entry["validity_state"] == "expert_curated"
    assert entry["has_conflicting_evidence"] is False

    bundle = _read(site, entry["bundle"])
    assert bundle["assertions"], "a gene page with no assertions is an empty page"
    assert bundle["publications"], "the page cites nothing"

    featured = _read(site, "featured.json")["featured"][0]
    assert isinstance(featured["publication"], dict), "a bare PMID needs a second fetch to render"
    assert featured["publication"]["title"]
    assert featured["blurb"]

    assert {record["kind"] for record in _search_records(site)} == {
        "gene",
        "publication",
        "phenotype",
    }


def test_tbx5_publishes_its_mirrored_validity_with_the_provenance_behind_it(
    site: Path,
) -> None:
    """The gene bundle's `validity` block, checked against the real mirrors.

    Not a fixture: `mirrors/clingen_gene_validity.tsv` curates TBX5 Definitive
    for Holt-Oram syndrome under SOP11, so the bundle's first record must carry
    that, and it must carry it *with its SOP version* — an attributed
    classification with no SOP is an unqualified claim, since the mirror spans
    SOP4 through SOP12 with no published crosswalk between them.

    `has_source_discordance` is asserted with `is`, not `==`: `{"f": 0} ==
    {"f": False}` is `True` in Python, so an equality alone would silently
    accept an `int` where the published contract is a `bool`.
    """
    bundle = _read(site, "genes/HGNC_11604.json")
    validity = bundle["validity"]

    assert validity["state"] == "expert_curated"
    assert validity["has_source_discordance"] is False

    clingen = validity["records"][0]
    assert clingen["source"] == "clingen"
    for field in ("source", "classification_term", "disease", "sop", "report_url"):
        assert clingen[field], f"{field} is empty on TBX5's ClinGen record"
