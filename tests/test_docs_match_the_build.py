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

from chd_atlas.build.profiles import LOOKUP_RULE, TAU_METHOD, TAU_SCALE
from chd_atlas.build.runner import build_site
from chd_atlas.models.dataset import Dataset

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


def _first_block_quote(text: str) -> str:
    """The first `>` block in `text`, joined, stopping at the blank line after it."""
    lines: list[str] = []
    for line in text.splitlines():
        if line.startswith(">"):
            lines.append(line.lstrip("> ").strip())
        elif lines:
            break
    return " ".join(lines)


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


def test_the_documented_manifest_census_is_the_census_the_build_publishes(site: Path) -> None:
    """`counts` and `schema_version`, in the doc, measured against a real build.

    Values rather than key sets, unlike the shape test above, and deliberately:
    the `manifest.json` block is the one example in this document whose numbers a
    reader takes as a *census of the atlas* rather than as an illustration. It
    said `"publications": 1` for a corpus holding 4 — true when it was written,
    false for three releases after — and nothing could see it, because a wrong
    census stays internally consistent (CLAUDE.md section 4.24).

    `schema_version` is pinned here as well as in `test_build_manifest.py`. That
    file compares the published value against a literal; this one checks the
    *document* names the same version, which is the claim a consumer branches on
    after reading it. A release that bumps the constant and not the prose leaves
    the two disagreeing with nothing to notice.

    Checksums and the commit stay placeholders and are not asserted: both change
    with every commit, so pinning them would guarantee this document is wrong by
    the next one.

    **Parsed, not substring-matched, and the difference is not cosmetic.** The
    first version of this test asserted `f'"{key}": {value}' in section`, which
    is unanchored: with the build at `burden_rows: 29` and the doc still saying
    `290`, the needle `"burden_rows": 29` is a substring of `"burden_rows": 290`
    and the test passes. Any documented figure having the build's figure as a
    prefix slips through — exactly the stale-census failure this guard exists to
    catch. Found by review 2026-08-06. `_example` already parses the fenced
    block, so comparing dicts is both stronger and less code.

    The whole `counts` object is compared, not key by key, so a key the document
    invents or drops fails here too.
    """
    doc = DOC.read_text()
    start = doc.index("## `manifest.json`")
    section = doc[start : doc.index("\n## ", start)]
    documented = _example("## `manifest.json`")
    manifest = json.loads((site / "manifest.json").read_text())

    assert f'"schema_version": "{manifest["schema_version"]}"' in section
    assert documented["counts"] == manifest["counts"], (
        f"the build counts {manifest['counts']}; the documented example says {documented['counts']}"
    )


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
    # **The count is read from the build, never written into the marker.** It
    # was hard-coded as "the 22 genes" here and in the doc, so the assertion
    # pinned a stale number *in place*: the widening moved it to 91 and this
    # test went on passing, because it was comparing the doc to itself.
    uncurated = sum(
        1
        for bundle in (site / "genes").glob("*.json")
        if bundle.name != "index.json"
        and json.loads(bundle.read_text())["atlas_curation"] == "not_yet_curated"
    )
    marker = (
        f"or, for the {uncurated} genes published today with no curation here, a paragraph saying"
    )
    assert marker in text, f"the doc must state the {uncurated} uncurated genes this build has"
    # Only the *first* quoted block after the marker. The doc quotes two notices
    # since 2026-08-06 -- one for a panel-graded gene, one for a
    # submitter-admitted one -- and joining every `>` line to the end of the file
    # welds them into a sentence no page has ever rendered.
    quoted = _first_block_quote(text[text.index(marker) :])
    assert quoted, "the doc must still quote the notice"
    # `**bold**` in the quote is `<strong>` on the page; strip both to compare
    # the sentence rather than its emphasis.
    assert quoted.replace("**", "") in _text_of(page)

    # **The second notice, on a gene no panel graded.** The doc quotes both, and
    # the ungraded one is the half that was wrong: every one of the 16
    # submitter-admitted pages carried "the classification above is an expert
    # panel's" until 2026-08-06, describing a commercial laboratory's submission
    # as a chartered panel's classification. A gene is chosen by reading
    # `validity_state` off the build rather than by naming an HGNC id, so this
    # keeps testing the right gene as the population moves.
    ungraded = next(
        bundle
        for bundle in sorted((site / "genes").glob("*.json"))
        if bundle.name != "index.json"
        and json.loads(bundle.read_text())["validity_state"] == "submitter_curated"
    )
    ungraded_page = _text_of(ungraded.with_suffix(".html").read_text())
    assert "no ClinGen expert panel has graded it" in ungraded_page
    assert "the classification above is an expert panel's" not in ungraded_page.lower()


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

    # The per-study split, and the gene the widened gate admitted with no burden
    # evidence at all. Both are new in schema 2.8 and both are the kind of claim
    # section 4.24 says goes stale silently: the first is four numbers that stay
    # internally consistent while describing a corpus that no longer exists, and
    # the second is a sentence that was simply not true before the gate widened.
    per_study = ", ".join(
        str(sum(1 for r in published if r["study"] == study)) for study in studies
    )
    assert f"{per_study} rows" in section, f"the per-study split is not {per_study}"

    covered = {
        path.stem
        for path in sorted((site / "genes").glob("HGNC_*.json"))
        if json.loads(path.read_text()).get("burden")
    }
    total = len(sorted((site / "genes").glob("HGNC_*.json")))
    assert f"across {len(covered)} of the {total} published genes" in section
    assert (len(covered) < total) == ("carries no burden rows at all" in section), (
        "the doc and the build disagree about whether a published gene lacks burden rows"
    )

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


def test_the_omics_section_documents_profiles_stratified_top_not_ranked_by_significance() -> None:
    """16a: `top` was documented as "ranked by significance" with no
    qualification, which Task 8 made false for `profiles` -- that table has
    no significance column at all, so its slice is stratified: the cardiac
    series leads on derived percentile, with at least one row of every other
    tissue reserved so a gene's own τ stays auditable from the same payload
    that publishes it. A consumer reading the old, unscoped sentence expects
    the 25 highest-ranked rows for every modality and gets something
    deliberately different for this one.
    """
    doc = DOC.read_text()
    start = doc.index("## `omics/<modality>/<accession>.json`")
    section = re.sub(r"\s+", " ", doc[start : doc.index("\n## ", start)])

    # Scoped to the three modalities it is actually true of, not stated as a
    # blanket rule every modality's `top` follows.
    assert "For `expression`, `proteomics` and `phospho`" in section
    assert "ranked by significance" in section
    assert "profiles` ranks `top` differently" in section
    assert "stratified" in section
    assert "median_percentile" in section
    assert "at least one row of every other tissue" in section
    assert "reserved" in section
    assert "D39(b)" in section
    assert "ranking bug" in section


def test_the_datasets_field_list_names_every_field_the_model_declares() -> None:
    """16b: the field list for `datasets.json` must name every field `Dataset`
    publishes, so a field added to the model and not to this sentence fails a
    test instead of aging into a false sentence the way Task 3's six fields
    already did once, silently, before this guard existed (CLAUDE.md section
    4.24 -- true when written, false by the next release, and internally
    consistent throughout, so re-reading it finds nothing).

    Read from `Dataset.model_fields`, never from a real corpus record: no
    dataset has ever been committed, so a real build's `datasets.json` is
    always `{"datasets": []}` and there is no record on disk whose keys this
    test could read instead. `build/literature.py::_dump` calls
    `model_dump(mode="json")` with no field filtering, so the model's own
    declared fields are exactly what a future record will publish -- reading
    them here is not an approximation of the build, it is what the build
    does, which is what makes this pin survive a field being *added* rather
    than only catching one being *removed from the doc*.
    """
    doc = DOC.read_text()
    start = doc.index("## `datasets.json`")
    section = doc[start : doc.index("\n## ", start)]

    for field in Dataset.model_fields:
        assert f"`{field}`" in section, (
            f"Dataset.{field} is not named in the datasets.json field list"
        )


def test_the_expression_profile_layer_is_documented_end_to_end() -> None:
    """16c: the two mirror tables, the `profile_quantiles` payload and the
    `expression_profile` bundle key were entirely undocumented before this.
    Checks the load-bearing claims a consumer needs, sharing one fixture
    because they are all facts about one new section:

    - the method travels with the number -- `LOOKUP_RULE`, `TAU_SCALE` and
      `TAU_METHOD` quoted from the constants that actually ship, not retyped,
      so a wording change to either constant must be echoed here or this
      fails;
    - why the quantile grid is published at all -- a one-level trade, stated
      as one, never as a proof;
    - percentiles are not comparable across organs, and -- more importantly,
      because the layer's headline question rides on it -- not comparable
      across developmental stages;
    - τ carries no adjective and no band (D39(c));
    - τ is blind to where a gene peaks. Verified directly against
      `profiles.specificity` before writing the doc, not reasoned from a
      number handed down: `specificity({"heart": 20.0, "liver": 200.0, "o1":
      5.0, "o2": 5.0, "o3": 5.0, "o4": 5.0, "o5": 5.0}, floor=0.0)` returns
      `tau == 0.6227723406454512` and `highest_in == "liver"`; the same inputs
      score `0.9625` computed by hand on the raw values with no log transform.
      Both round to the figures asserted below.
    """
    doc = DOC.read_text()

    assert "## `omics/profile_quantiles/<accession>.json`" in doc

    start = doc.index("### The bundle's `expression_profile` object")
    section = re.sub(r"\s+", " ", doc[start : doc.index("\n## ", start)])

    assert LOOKUP_RULE in section
    assert TAU_SCALE in section
    assert TAU_METHOD in section

    assert "one-level guarantee" in section
    assert "not an unbounded one" in section
    assert "re-derivable" in section

    assert "not comparable across organs" in section
    assert "not comparable across developmental stages" in section
    assert "matters more" in section

    assert "no adjective, no band" in section

    assert "0.963" in section
    assert "0.623" in section
    assert "liver" in section
    assert "highest_in" in section


def test_the_api_doc_states_the_expression_census_the_build_produces(site: Path) -> None:
    """Section 4.24: when a release adds a data source, every count in every
    prose artifact is suspect by default. Grep the numbers and re-derive them
    from a real build; do not re-read the sentences. Eight false claims
    survived multiple readings last time because they still agreed with each
    other.

    Both figures are 0 on the committed corpus -- there is no `profiles`
    mirror -- and the doc has to say so in exactly those terms, not merely in
    a rounded-off "none yet".
    """
    manifest = json.loads((site / "manifest.json").read_text())
    doc = DOC.read_text()

    assert f"{manifest['counts']['profile_genes']} genes" in doc
    assert f"{manifest['counts']['profile_datasets']} datasets" in doc
