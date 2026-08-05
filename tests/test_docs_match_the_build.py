# tests/test_docs_match_the_build.py
"""Pin the consumer documentation against what the build actually publishes.

`docs/data-api.md` opens by claiming every example is copied from a real build
rather than written by hand. Nothing checked that, and it drifted: the
`genes/index.json` example lost `conflicting_lesion_groups` — the field
`schema_version` 1.1 was minted for, and the one the document's own closing
section calls "the one consumer obligation this API imposes".

The consequence was not a missing line. A consumer typing their browse row off
an exhaustive, alphabetically-ordered example ships without the per-group
contested flag, and for a gene that is definitive-and-refuted on septal disease
while merely moderate on conotruncal, what they render inverts the science: the
definitive call reads as the solid one when it is precisely the disputed one.

So this is not a style check. It is the same guarantee the dead-link sweep gives
the payloads, applied to the document a consumer builds against.
"""

import html
import json
import re
import shutil
from pathlib import Path
from typing import Any

import pytest

from chd_atlas.build.runner import build_site

REPO = Path(__file__).parent.parent
DOC = REPO / "docs" / "data-api.md"

# Each entry maps a documented example to the published object it claims to
# describe. Only *exhaustive* examples are listed — blocks written with an
# ellipsis are illustrations of shape, not of field sets, and are excluded here
# rather than silently half-checked.
#
# `heading` locates the fenced block; `select` pulls the object out of the built
# site that the block is an example of.
EXHAUSTIVE: dict[str, str] = {
    "## `manifest.json`": "manifest",
    "## `genes/index.json`": "index_row",
    "## `publications.json`": "publication",
    "## `featured.json`": "featured",
    "## `phenotypes.json`": "phenotype",
    "## `search/index.json.gz`": "search_record",
    "### The bundle's `validity` object: mirrored, attributed, never authored here": (
        "validity_record"
    ),
}


@pytest.fixture(scope="module")
def site(tmp_path_factory: pytest.TempPathFactory) -> Path:
    base = tmp_path_factory.mktemp("docs")
    root = base / "repo"
    root.mkdir()
    for name in ("curation", "mirrors"):
        shutil.copytree(REPO / name, root / name)
    (root / "ontologies").symlink_to(REPO / "ontologies")
    out = base / "dist"
    build_site(root, out)
    return out


def _text_of(page: str) -> str:
    """One HTML page as its visible text, whitespace collapsed.

    Tags stripped rather than parsed: the assertion is about a sentence a
    reader sees, and `<strong>` inside it is exactly what a substring check on
    the raw markup would trip over.
    """
    return re.sub(r"\s+", " ", html.unescape(re.sub(r"<[^>]+>", " ", page))).strip()


def _example(heading: str) -> dict[str, Any]:
    """The first fenced JSON block under one heading, parsed."""
    text = DOC.read_text()
    start = text.index(heading)
    nxt = text.find("\n## ", start + 1)
    section = text[start : nxt if nxt != -1 else len(text)]
    match = re.search(r"```json\n(.*?)```", section, re.S)
    assert match, f"no JSON example under {heading}"
    payload: dict[str, Any] = json.loads(match.group(1))
    return payload


def _published(site: Path, kind: str) -> dict[str, Any]:
    import gzip

    def load(rel: str) -> Any:
        raw = (site / rel).read_bytes()
        return json.loads(gzip.decompress(raw) if rel.endswith(".gz") else raw)

    if kind == "manifest":
        return dict(load("manifest.json"))
    if kind == "index_row":
        # TBX5 by id, not `genes[0]`. The example under that heading is TBX5's
        # row and the first row is now TBX20 (HGNC:11598) -- D21 took the index
        # from 1 gene to 23, sorted by HGNC id. Every row carries the same keys
        # by construction (`bundles.py::_headline` plus one literal dict), so
        # this test would still pass on the wrong row; selecting the documented
        # gene is what keeps it checking the example against what the example
        # claims to be.
        row = next(
            item for item in load("genes/index.json")["genes"] if item["gene"] == "HGNC:11604"
        )
        return dict(row)
    if kind == "publication":
        return dict(load("publications.json")["publications"][0])
    if kind == "featured":
        return dict(load("featured.json")["featured"][0])
    if kind == "phenotype":
        return dict(load("phenotypes.json")["phenotypes"][0])
    if kind == "search_record":
        return dict(load("search/index.json.gz")["records"][0])
    if kind == "validity_record":
        # TBX5 (HGNC:11604) is the committed corpus's one asserted gene, and its
        # first record is ClinGen's -- `build.validity._sort_key` sorts
        # `(source, disease, moi, submitter)`, and "clingen" precedes "gencc"
        # lexically.
        return dict(load("genes/HGNC_11604.json")["validity"]["records"][0])
    raise AssertionError(f"unknown kind {kind}")


@pytest.mark.parametrize(("heading", "kind"), sorted(EXHAUSTIVE.items()))
def test_a_documented_example_carries_every_field_the_build_publishes(
    site: Path, heading: str, kind: str
) -> None:
    """Key sets, not values: the doc illustrates content, but promises shape.

    Values legitimately differ — the example may show a contested gene where the
    committed corpus has none — so comparing them would make the document a
    second copy of the corpus. What a consumer types off the block is the *field
    set*, and that is what has to be exhaustive.
    """
    documented = _example(heading)
    published = _published(site, kind)

    # The example may be a container (`{"genes": [ … ]}`) or the object itself.
    if len(documented) == 1:
        only = next(iter(documented.values()))
        if isinstance(only, list) and only and isinstance(only[0], dict):
            documented = only[0]

    missing = sorted(set(published) - set(documented))
    invented = sorted(set(documented) - set(published))
    assert not missing, f"{heading}: the build publishes {missing}, the doc does not show it"
    assert not invented, f"{heading}: the doc shows {invented}, the build does not publish it"


def test_the_contested_example_is_labelled_with_the_payload_it_comes_from(site: Path) -> None:
    """`confidence_by_lesion_group` and its flag are browse-row fields only.

    The closing section shows the three contested fields in one object and tells
    a consumer to read them together. That object exists in `genes/index.json`
    and nowhere else — the gene bundle carries `has_conflicting_evidence` but
    neither per-group field — while the section above says the gene-level flag
    "appears in both the browse row and the bundle", which invites exactly the
    wrong inference.

    A detail-page author who takes the obligation to the bundle finds nothing to
    implement it with, and cannot derive it either: reconstructing per-group
    confidence means reimplementing `strongest()`'s rank and
    `has_conflicting_evidence`'s treatment of `no_known_association`, neither of
    which is published. So the document has to say where the fields live.
    """
    bundle = json.loads((site / "genes" / "HGNC_11604.json").read_text())
    # By id, for the reason `_published` gives: `genes[0]` is TBX20 since D21,
    # and the bundle this is compared against is TBX5's.
    index_row = next(
        item
        for item in json.loads((site / "genes" / "index.json").read_text())["genes"]
        if item["gene"] == "HGNC:11604"
    )

    assert "conflicting_lesion_groups" in index_row
    assert "conflicting_lesion_groups" not in bundle, "if this changes, the doc must too"
    assert "confidence_by_lesion_group" not in bundle

    text = DOC.read_text()
    section = text[text.index("## Contested genes") :]
    # Specific enough that the earlier sentence "it appears in both the browse
    # row and the bundle" — which is about the gene-level flag and is exactly
    # what misleads here — cannot satisfy it.
    assert "and **nowhere else**" in section, (
        "the section must say the per-group fields are browse-row only"
    )
    assert "conflicting_lesion_groups" in section


def test_the_documented_uncurated_notice_is_the_sentence_the_page_actually_shows(
    site: Path,
) -> None:
    """The doc quotes that paragraph verbatim, so drift makes the doc a lie.

    Both of its clauses are narrower than the obvious wording, and both were
    made narrower deliberately: "not yet curated **a lesion assertion**" rather
    than "evidence", because `atlas_curation` is derived from `LesionAssertion`
    records alone and a gene may carry curated functional evidence while
    reporting `not_yet_curated`; and "**no classification** on this page is the
    atlas's own" rather than "nothing on this page", which would deny that same
    curated work one column from the rail counting it.

    The document was written against the earlier, wider wording and kept it
    through the fix, reintroducing in prose the ambiguity the code had just
    removed. `EXHAUSTIVE` above pins documented *JSON* examples against the
    build and had nothing to say about a quoted sentence; this closes that.
    Raised by review on #14.

    Compared against a real uncurated gene page rather than against
    `pages._NOT_CURATED`, so the assertion fails if the constant stops reaching
    the page at all -- comparing the doc to the constant would pass on a page
    that renders neither.
    """
    page = (site / "genes" / "HGNC_4173.html").read_text()
    assert "not_yet_curated" in page, "GATA4 must still be an uncurated gene for this to test"

    text = DOC.read_text()
    marker = "or, for the 22 genes published today with no curation here, a paragraph saying"
    quote = text[text.index(marker) :]
    quoted = " ".join(
        line.lstrip("> ").strip() for line in quote.splitlines() if line.startswith(">")
    )
    assert quoted, "the doc must still quote the notice"
    # `**bold**` in the quote is `<strong>` on the page; strip both to compare
    # the sentence rather than its emphasis.
    assert quoted.replace("**", "") in _text_of(page)


def test_the_burden_census_in_the_doc_is_the_census_the_build_publishes(site: Path) -> None:
    """Every count in the `burden` section, asserted against a real build.

    **This is the guard for the failure that actually happened.** Commit 770aa2a
    added a third study to the mirror -- 180 rows, 90 of them published, the
    first `mutation_model` comparator the atlas has ever carried -- and edited
    `docs/data-api.md` only to add the `count_unit` table. Every census figure in
    that section, and the sentence "Only `control_cohort` appears today", went on
    describing the two-study corpus. An adversarial review of the promotion diff
    measured them all false before the deploy.

    The numbers were internally consistent with each other, which is what made it
    survive reading: they were exactly the pre-Sierant subset. Only measurement
    catches that, so the doc's counts are now measured rather than read.

    The mirror figures are counted from the file rather than the build because
    the whole point of the paragraph is the gap between them -- the mirror is
    deliberately wider than the publication gate.
    """
    import csv

    doc = DOC.read_text()
    start = doc.index("### The bundle's `burden` array")
    # Whitespace collapsed: these are claims in prose, and where a sentence
    # happens to wrap is not part of the claim. Asserting against the raw text
    # would make this test fail on a reflow and pass on a wrong number.
    section = re.sub(r"\s+", " ", doc[start : doc.index("\n## ", start)])

    published = [
        row
        for path in sorted((site / "genes").glob("HGNC_*.json"))
        for row in json.loads(path.read_text()).get("burden", [])
    ]
    mirror = list(csv.DictReader((REPO / "mirrors" / "burden.tsv").open(), delimiter="\t"))

    studies = sorted({row["study"] for row in published})
    comparators = {
        c: sum(1 for r in published if r["comparator"] == c)
        for c in {r["comparator"] for r in published}
    }
    synonymous = sum(1 for r in published if r["consequence_class"] == "synonymous")
    lof = sum(1 for r in published if r["consequence_class"] == "lof")
    composite = sum(1 for r in published if r["consequence_class"] == "damaging")
    unbounded = sum(1 for r in published if r["effect_bound"] is not None)

    # Each claim, in the exact spelling the document uses.
    assert f"**{len(published)} rows reach the API**" in section
    assert f"{len(studies)} studies ({', '.join(studies)})" in section
    assert f"holds {len(mirror):,} rows for {len({r['gene'] for r in mirror})} genes" in section
    assert f"the other {len(mirror) - len(published):,}" in section

    # The sentence a consumer writes code against. It read "Only
    # `control_cohort` appears today" while 45 rows were `mutation_model`.
    assert f"`control_cohort` ({comparators['control_cohort']} rows)" in section
    assert f"`mutation_model` ({comparators['mutation_model']} rows)" in section

    assert f"{unbounded} published rows today" in section
    assert f"{synonymous} of the {len(published)} published rows are synonymous" in section
    assert f"loss-of-function is {lof} and synonymous {synonymous}" in section
    assert f"{composite} of the {len(published)} rows" in section

    # And the claim that the composite is a union rather than a partition, which
    # is what stops a consumer double-counting when it aggregates.
    assert "double-counts" in section


def test_the_independent_datasets_section_states_the_shape_the_build_publishes(
    site: Path,
) -> None:
    """The census guard, applied to the new section before it can go stale.

    The burden section drifted by a whole study because its counts were written
    once and never re-measured, and the numbers stayed internally consistent so
    reading them found nothing. This section is pinned from the day it lands
    rather than after the same thing happens to it.

    The three display obligations are asserted as *text*, not as numbers,
    because they are what stops the object being rendered as a verdict — and a
    later edit that softens them should fail here rather than ship.
    """
    doc = DOC.read_text()
    start = doc.index("### `independent_datasets`")
    section = re.sub(r"\s+", " ", doc[start : doc.index("\n## ", start)])

    published = [
        json.loads(path.read_text())["independent_datasets"]
        for path in sorted((site / "genes").glob("HGNC_*.json"))
    ]
    states = {entry["state"] for gene in published for entry in gene["families"]}

    # Every state the build can emit is documented; a fifth would fail here.
    assert states <= {"corrected", "nominal", "no_enrichment", "not_tested"}
    for state in states:
        assert f"| `{state}` |" in section, f"{state!r} is emitted and undocumented"

    # `tested` is the denominator, and it is never `len(families)` where any
    # family did not test the gene. Measured on the real corpus rather than
    # asserted, so the doc's claim is checked against what ships.
    assert any(gene["tested"] < len(gene["families"]) for gene in published), (
        "no published gene has an untested family, so the denominator claim is untested"
    )

    assert "never be rendered as one" in section
    assert "never `len(families)`" in section
    assert "KDM6A causes Kabuki syndrome" in section
