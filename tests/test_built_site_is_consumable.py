# tests/test_built_site_is_consumable.py
import gzip
import json
import math
import posixpath
import re
import shutil
from pathlib import Path
from typing import Any

import pytest

from chd_atlas.build.render import RESEARCH_USE_NOTICE
from chd_atlas.build.runner import build_site
from chd_atlas.corpus import load_curation

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
    add it here -- which is the failure mode this replaces.

    **The count is asserted against the published population, not against a
    literal.** An `rglob` matching nothing would satisfy an all-pages loop
    vacuously, so a count is needed; a literal 25 was what the committed corpus
    published until the gate widened on 2026-08-06, and pinning the total meant
    re-editing this line for a change it has nothing to say about. Derived from
    the bundle count instead: one page per published gene, plus the landing page
    and the browse page. That fails on a gene page that was not written and on a
    stray page that was, which is what this guards, and it does not fail merely
    because the population moved.
    """
    pages = sorted(site.rglob("*.html"))
    bundles = sorted((site / "genes").glob("HGNC_*.json"))
    expected = len(bundles) + 2

    assert len(bundles) > 1, "the build published no gene bundles; the fixture is broken"
    assert len(pages) == expected, (
        f"expected {expected} HTML pages (1 landing + 1 browse + {len(bundles)} gene "
        f"pages), found {len(pages)}"
    )
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


def _page_links(site: Path) -> list[tuple[str, str, str]]:
    """Every internal `href` on every built page, with the page that carries it.

    Returns `(page, href, resolved)` where `resolved` is the site-relative path
    the browser would request. External schemes and pure fragments are skipped:
    this sweep is about links the build is responsible for.
    """
    links: list[tuple[str, str, str]] = []
    for path in sorted(site.rglob("*.html")):
        page = path.relative_to(site).as_posix()
        for href in re.findall(r'href="([^"]+)"', path.read_text(encoding="utf-8")):
            if href.startswith(("http://", "https://", "mailto:", "#")):
                continue
            target = posixpath.normpath(posixpath.join(posixpath.dirname(page), href))
            links.append((page, href, target))
    return links


def test_every_link_on_every_page_resolves_to_a_file_the_build_wrote(site: Path) -> None:
    """The HTML half of the sweep above, which had no HTML half at all.

    `_advertised` walks manifest entries, gene bundles, omics and variant shards
    and search `path`s -- every one of them JSON. It never opens an `.html`
    file. This release added 25 pages and, measured here, 100+ internal links,
    none of which anything checked.

    This project has already shipped exactly this failure once: `runner.py`
    records that the nav and the landing page both linked `genes/index.html`
    while no builder wrote one -- green build, verifying checksums, 404 for
    every visitor. It was fixed by hand and no sweep was added.

    Measured before this test existed: changing the browse row's href from
    `f"../{gene_page_path(...)}"` to `gene_page_path(...)` -- which points all
    23 rows at `genes/genes/HGNC_*.html` -- left the whole 624-test suite green.
    `test_every_browse_row_links_to_a_page_that_was_written` could not catch it
    because it asserts a bare substring (`name in page`) and a key in
    `emitter.checksums`, neither of which sees the `../` prefix.

    Resolution is `posixpath.normpath` against the page's own directory, so a
    link that escapes the site root normalises to a path starting `..` and fails
    the membership test rather than silently resolving somewhere on disk.
    """
    written = {path.relative_to(site).as_posix() for path in site.rglob("*") if path.is_file()}
    links = _page_links(site)

    assert links, "no page links found; the sweep would pass vacuously"

    broken = [
        (page, href)
        for page, href, target in links
        # A directory URL is served by its `index.html`, which is how
        # `genes/index.html` is reached as `genes/`.
        if target not in written and f"{target}/index.html" not in written
    ]
    assert not broken, f"pages link to files the build never wrote: {broken}"


def test_every_cohort_a_published_burden_row_cites_resolves_with_its_caveat(site: Path) -> None:
    """A bare cohort id in a bundle must resolve to a description in `cohorts.json`.

    **This is the join the burden columns exist for.** A row names its sample
    collections by bare id -- `["taa_cases"]`, `["gnomad_controls"]` -- and the
    ids carry no meaning on their own. Until schema 2.9 nothing published mapped
    them, so a consumer had the statistics and none of the caveats: measured
    2026-08-06 over a real build, 13 distinct ids appeared across 915 published
    rows and `curation/cohorts.yaml` reached zero published bytes.

    Asserted over a *real* build rather than a fixture, because what could break
    it is not the emitter -- covered by
    `test_cohorts_json_publishes_the_registry_a_burden_row_resolves_against` --
    but the publication gate moving under it. A widened gate can admit a gene
    whose rows cite a collection nobody curated, and that is a dangling id in a
    published payload rather than a build failure. `BUR009` refuses it at the
    gate; this checks the guarantee survives to the bytes.

    `taa_cases` is named explicitly. It is 777 thoracic aortic aneurysm probands
    who **do not have congenital heart disease**, and it is the id whose
    description most changes what a row means -- a reader who resolves it and a
    reader who does not are reading different data.
    """
    registry = {
        cohort["id"]: cohort
        for cohort in json.loads((site / "cohorts.json").read_text(encoding="utf-8"))["cohorts"]
    }
    assert registry, "the build published no cohort registry"

    cited: set[str] = set()
    for bundle in sorted((site / "genes").glob("HGNC_*.json")):
        for row in json.loads(bundle.read_text(encoding="utf-8"))["burden"]:
            cited.update(row["case_cohorts"])
            cited.update(row["control_cohorts"])

    assert cited, "no published burden row names a cohort, so this proves nothing"
    dangling = sorted(cited - set(registry))
    assert not dangling, f"published rows cite cohorts that resolve to nothing: {dangling}"

    # Every resolved record carries the sentence a reader needs, not just a name.
    for identifier in sorted(cited):
        assert registry[identifier]["description"].strip(), f"{identifier} publishes no caveat"

    assert "not congenital heart disease" in registry["taa_cases"]["description"].lower()


def test_a_gene_an_authority_reported_no_association_for_says_so_in_both_payloads(
    site: Path,
) -> None:
    """Issue #13's third axis, on the one published gene that exercises it.

    GDF1 carries G2P `Definitive`, Labcorp `Strong` and Illumina `No Known
    Disease Relationship` in scope. Until schema 2.10 it published
    `has_conflicting_evidence: false` and nothing else, so a consumer was told
    the evidence did not conflict while two bodies disagreed about whether an
    association exists at all. `no_known_association` takes neither side of
    `CONTESTED` -- correctly, since a null result is not a refutation -- and the
    disagreement therefore reached no published byte.

    Asserted on the bundle **and** the browse row, because a flag a reader meets
    on one and not the other is missing exactly where they are standing, and
    asserted against a real build rather than a fixture because what breaks it is
    the gate moving: GDF1 is admitted on GenCC agreement, and it is the only
    published gene on this axis, so a widening or narrowing that drops it leaves
    the axis published on nothing while every unit test still passes.

    `has_conflicting_evidence` is asserted **false** in the same breath. The
    whole design decision is that these are two axes rather than one, and a test
    that checked only the new flag would pass with the two collapsed together.
    """
    bundle = json.loads((site / "genes" / "HGNC_4214.json").read_text(encoding="utf-8"))
    assert bundle["symbol"] == "GDF1", "this test is about GDF1"

    assert bundle["has_no_association_report"] is True
    assert bundle["no_association_reported_by"] == ["Illumina"]
    # The other axis stays untouched: a null result is not a contest.
    assert bundle["has_conflicting_evidence"] is False

    row = next(
        item
        for item in json.loads((site / "genes" / "index.json").read_text(encoding="utf-8"))["genes"]
        if item["gene"] == "HGNC:4214"
    )
    assert row["has_no_association_report"] is True
    assert row["no_association_reported_by"] == ["Illumina"]
    assert row["has_conflicting_evidence"] is False

    # The pair is consistent everywhere, not only on this gene: an authority
    # named while the flag reads false would say two things at once.
    for bundle_path in sorted((site / "genes").glob("HGNC_*.json")):
        payload = json.loads(bundle_path.read_text(encoding="utf-8"))
        assert payload["has_no_association_report"] == bool(payload["no_association_reported_by"])


# Every attribute the chart renderers position anything with. `points` and
# `viewBox` hold several numbers each and are split on whitespace and commas.
_GEOMETRY_ATTRIBUTES = (
    "data-scale-high",
    "viewBox",
    "height",
    "points",
    "width",
    "x1",
    "x2",
    "y1",
    "y2",
    "cx",
    "cy",
    "x",
    "y",
    "r",
)
# Longest-first. Measured: the order is not load-bearing -- the alternation
# backtracks, so `x` failing on the `1` of `x1="..."` lets `x1` match anyway --
# but written this way so the scan does not depend on that.
_GEOMETRY = re.compile(r"\b(?:" + "|".join(_GEOMETRY_ATTRIBUTES) + r')="([^"]*)"')

# What `charts.coordinate` produces: an integer, or an integer and exactly one
# decimal place. A raw `repr` reaches `0.30000000000000004` and `1e-05`.
_FIXED_PRECISION = re.compile(r"^-?\d+(\.\d)?$")


def test_every_geometry_the_pages_publish_has_coordinate_s_fixed_shape(site: Path) -> None:
    """`coordinate` is pinned; that it is *used* was not, until 2026-08-20.

    `test_build_charts.py::test_coordinates_are_fixed_precision_so_two_builds_
    agree` proves the function rounds. It cannot prove any renderer calls it,
    so a new chart interpolating a raw float -- the ordinary way to write an
    SVG attribute -- passes every test in this repository while publishing
    `y="43.33333333333333"`. That is the whole reason `coordinate` exists:
    binary floating point is not identical across platforms at the last bits,
    so a raw repr makes two builds of one commit differ by architecture while
    every checksum still verifies against itself.

    Scanned off the built pages rather than off the source, for the reason
    CLAUDE.md section 4.35 gives: the question is what a reader's browser
    receives, and `grep` over `src/` answers a different one -- these strings
    are assembled from adjacent literals and f-strings that no source-level
    search reconstructs.

    Measured 2026-08-20 over a real build: 94 pages, 41,833 geometry tokens,
    0 non-conforming.

    **The count is derived, never pinned to a literal**, for the reason
    `test_every_html_page_the_build_writes_carries_the_research_use_notice`
    derives its own: a scan matching nothing satisfies a loop vacuously, and a
    literal total is re-edited every time the population moves. Every gene page
    must carry geometry -- measured, all 92 do, though not all for the same
    reason: GDF1 draws three forest panels and no organ chart at all, and
    CFC1's 138 tokens are the fewest. The landing and browse pages carry none
    today and are scanned anyway, so a chart added to either is covered
    without anyone remembering this test.
    """
    pages = sorted(site.rglob("*.html"))
    gene_pages = sorted((site / "genes").glob("HGNC_*.html"))
    assert len(gene_pages) > 1, "the build published no gene pages; the fixture is broken"

    total = 0
    for page in pages:
        tokens = [
            token
            for value in _GEOMETRY.findall(page.read_text(encoding="utf-8"))
            for token in value.replace(",", " ").split()
        ]
        if page in gene_pages:
            assert tokens, f"{page.relative_to(site)} publishes no geometry at all"
        total += len(tokens)
        for token in tokens:
            assert _FIXED_PRECISION.match(token), (
                f"{page.relative_to(site)} publishes the geometry {token!r}, which did not "
                "come through charts.coordinate -- two builds on two platforms can differ"
            )

    assert total > len(gene_pages), f"only {total} geometry tokens over {len(pages)} pages"


# Each `<rect class="chart-band">` a trajectory draws, paired with the phase
# label its own `<title>` names. Anchored on the `<rect>` rather than on the
# bare class name, because `render.STYLESHEET` is inlined verbatim into every
# page and carries `.chart-band { ... }` — a `"chart-band" in page` check
# passes on all 94 pages with every band deleted (CLAUDE.md section 4.35).
_BAND = re.compile(r'<rect class="chart-band"[^>]*><title>([^<]*)</title>')


def test_the_curated_phase_vocabulary_reaches_the_pages_it_is_drawn_on(site: Path) -> None:
    """The wiring, not the renderer: `build_gene_pages(phases=...)` in `runner.py`.

    Every page test hands `phases=` a fixture explicitly, so the keyword
    `runner.py` actually passes was exercised by nothing. Measured 2026-08-20,
    changing that one argument to `None` passed all 1,053 tests while changing
    86 of 198 published files: **85 pages lost their trajectory entirely** and
    30 of them gained, in bold, "This gene reads *below this dataset's
    detection floor at every stage sampled* in whole heart (19 stages
    sampled)" — including HGNC:10249, whose own bundle places 18 of those 19
    heart stages against the percentile grid with the percentiles printed in
    the table directly below. A false claim, on 30 pages, in bold, under a
    green build. This repository's characteristic defect exactly.

    Three assertions, because the wiring can fail in three ways and only the
    first is the mutant above:

    * bands are drawn at all — `phases=None` publishes zero;
    * every label drawn is one a curator wrote in
      `curation/cardiac_phases.yaml` — a renderer inventing its own phase
      names would satisfy the first;
    * the one phase whose end the source never states is banded nowhere.
      `_banded_phases` refuses it and `CardiacPhaseFile.phases_for` refuses
      it, and this is the only check that either rule survives the real
      pipeline rather than a fixture. A band has to run to an edge, so
      drawing one for `heart_looping` would assert in the encoding a reader
      takes at a glance the very boundary `end_basis: not_stated` exists to
      refuse.

    Measured 2026-08-20 on the committed corpus: 85 gene pages, 4 bands each,
    4 distinct labels of the 5 the vocabulary makes eligible.
    `embryonic_heart_tube_morphogenesis` ends before this dataset's first
    sampled stage, so it is eligible and unbanded — which is why the labels
    are asserted as a non-empty subset rather than pinned to a literal set
    that a second dataset would falsify.
    """
    phases = load_curation(REPO)[0].cardiac_phases
    assert phases is not None, "the committed corpus has no phase vocabulary; fixture is broken"
    curated = {phase.label for phase in phases.phases}
    unended = {phase.label for phase in phases.phases if phase.end_wpc is None}
    assert unended, "no unended phase is curated; the third assertion below is vacuous"

    banded: set[str] = set()
    pages_with_a_band = 0
    for page in sorted((site / "genes").glob("HGNC_*.html")):
        labels = _BAND.findall(page.read_text(encoding="utf-8"))
        if labels:
            pages_with_a_band += 1
        banded.update(labels)

    assert pages_with_a_band, "no page carries a phase band; the vocabulary reached no page"
    assert banded <= curated, f"a band names a phase no curator wrote: {sorted(banded - curated)}"
    assert banded, "bands were drawn with no phase named in any of them"
    assert not banded & unended, (
        f"a phase whose end the source never states is banded: {sorted(banded & unended)}"
    )


def _tau_denominator(method: str, n_organs: int) -> float:
    """The divisor the published `method` string *names*, read out of the words.

    Deliberately a parser and not a constant. D39(a)'s whole promise is that a
    consumer who reads `method` and never opens this repository reproduces the
    number; the only way to test that promise is to follow the sentence rather
    than the source. A constant `n - 1` here would agree with `specificity`
    for the same reason the code agrees with itself, and would have passed
    unchanged while the string said "mean over organs" — which is what
    happened, through four review passes.

    Raises rather than defaulting on an unrecognised wording: a rewrite this
    function cannot follow is a rewrite a consumer cannot follow either, and
    silently picking a denominator would restore exactly the hole this closes.
    """
    if "divided by (n - 1)" in method:
        return float(n_organs - 1)
    if "mean over organs" in method:
        return float(n_organs)
    raise AssertionError(
        f"the published tau method names no denominator this test can follow: {method!r}"
    )


def test_every_published_tau_is_reproducible_from_its_own_method_string(site: Path) -> None:
    """D39(a) made executable: `method` + `medians` must reproduce `tau`.

    `specificity` publishes the formula beside the number precisely so a
    consumer can re-derive it without reading `build/profiles.py`. Nothing
    checked that the sentence and the arithmetic agreed, and they did not:
    the string said "mean over organs of (1 - x_i/x_max)" — a divisor of *n* —
    while the code divides by *n − 1*, the standard Yanai τ. Measured
    2026-08-20 over the committed corpus: all 1817 published specificity
    blocks match `sum/(n − 1)` and 4 also match `sum/n`, those four being
    coincidences of ties rather than agreement.

    So this asserts equality *and* that the two candidate divisors are
    distinguishable on this corpus — per section 4.30, two figures that are
    equal in the fixture are one figure to the test. Without the second
    assertion a corpus of two-organ stages (where n − 1 == 1 and the
    difference is a factor of 2 on every block, but a degenerate one-organ
    corpus would collapse) could let a wrong string pass.

    `scale` is asserted alongside because `method` names the transform and the
    two are separately published; a consumer following one and reading the
    other must not get two different answers.
    """
    blocks = 0
    also_matching_plain_mean = 0
    for bundle in sorted((site / "genes").glob("HGNC_*.json")):
        payload = json.loads(bundle.read_text(encoding="utf-8"))
        for dataset in payload["expression_profile"]["datasets"]:
            for stage in dataset["stages"]:
                specificity = stage["specificity"]
                if specificity is None:
                    continue
                blocks += 1
                medians = specificity["medians"]
                method = specificity["method"]

                assert specificity["scale"] == "log2(x+1)"
                assert "x = log2(median+1)" in method
                assert "clamped to 0" in method

                logged = {
                    organ: math.log2(max(value, 0.0) + 1.0) for organ, value in medians.items()
                }
                x_max = max(logged.values())
                total = sum(1.0 - value / x_max for value in logged.values())

                derived = total / _tau_denominator(method, len(medians))
                assert derived == pytest.approx(specificity["tau"], abs=1e-12), (
                    f"{bundle.name}: tau {specificity['tau']} is not what the published "
                    f"method reproduces ({derived}) from the published medians"
                )
                if total / len(medians) == pytest.approx(specificity["tau"], abs=1e-12):
                    also_matching_plain_mean += 1

    assert blocks, "no specificity block was published; this test measured nothing"
    assert also_matching_plain_mean < blocks, (
        "every published tau matches sum/n as well as sum/(n-1); this corpus "
        "cannot distinguish the two divisors, so the assertion above is vacuous"
    )
