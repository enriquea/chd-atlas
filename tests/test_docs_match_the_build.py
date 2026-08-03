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
        return dict(load("genes/index.json")["genes"][0])
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
    index_row = json.loads((site / "genes" / "index.json").read_text())["genes"][0]

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
