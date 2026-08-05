# tests/unit/test_models_cohort.py
"""The cohort registry, and the committed file it validates."""

from __future__ import annotations

from pathlib import Path

import pytest
from pydantic import ValidationError
from ruamel.yaml import YAML

from chd_atlas.models.cohort import CohortFile

REPO_ROOT = Path(__file__).resolve().parents[2]

_ENTRY = {
    "id": "ddd",
    "name": "Deciphering Developmental Disorders (DDD)",
    "description": "UK-wide study of children with severe undiagnosed developmental disorders.",
}


@pytest.mark.parametrize(
    ("mutation", "match"),
    [
        # The registry's own rule, not pydantic's: two entries under one id make
        # a `;`-joined burden cell ambiguous about which description it means.
        pytest.param({}, "duplicate cohort ids", id="duplicate-id"),
        # `_NonBlankStr` rather than a bare `min_length=1`. Whitespace renders as
        # nothing on a page while passing every presence check, so the entry
        # looks answered instead of missing -- the defect `models/scope.py` and
        # `validate/burden.py::_blank` both exist to catch on their own side.
        pytest.param({"description": "   "}, "at least 1 character", id="blank-description"),
        pytest.param({"name": "\t\n"}, "at least 1 character", id="blank-name"),
        # The grammar is shared with the `;`-joined mirror columns, so an id
        # carrying a separator would split into two ids nothing defines.
        pytest.param({"id": "ddd;ukbb"}, "should match pattern", id="separator-in-id"),
        pytest.param({"id": "DDD"}, "should match pattern", id="capital-in-id"),
    ],
)
def test_a_registry_entry_that_could_not_be_resolved_is_refused(
    mutation: dict[str, str], match: str
) -> None:
    """Each case makes some burden row's `case_cohorts` cell unresolvable."""
    second = {**_ENTRY, **mutation} if mutation else _ENTRY
    with pytest.raises(ValidationError, match=match):
        CohortFile.model_validate({"cohorts": [_ENTRY, second] if not mutation else [second]})


def test_an_unknown_field_is_refused_rather_than_silently_dropped() -> None:
    """`extra="forbid"`. A curator writing `descripton:` would otherwise get an
    entry with no description and no error, which is this project's
    characteristic failure in miniature.
    """
    with pytest.raises(ValidationError, match="Extra inputs"):
        CohortFile.model_validate({"cohorts": [{**_ENTRY, "role": "case"}]})


def test_the_committed_registry_parses_and_covers_every_cited_cohort() -> None:
    """A contract pin: `mirrors/burden.tsv` cites these ids by name.

    Named rather than merely counted, because the point of the file is *which*
    collections a row draws on -- a registry that still had the right number of
    entries after `ddd` was renamed would satisfy a count and break every row
    citing it.

    `taa_cases` is in this set and is deliberately not a CHD cohort: 777 of the
    7,958 cases in PMID:34324492 have sporadic thoracic aortic aneurysm, which
    `curation/chd_scope.yaml` puts outside this atlas. It is registered so that
    every row drawing on it says so, rather than being quietly folded into a
    case count labelled CHD.
    """
    path = REPO_ROOT / "curation" / "cohorts.yaml"
    parsed = CohortFile.model_validate(YAML(typ="safe").load(path.read_text(encoding="utf-8")))

    assert {cohort.id for cohort in parsed.cohorts} == {
        # PMID:42230622
        "cnchd",
        "ddd",
        "nottingham",
        "ukbb",
        # PMID:34324492
        "signature",
        "decipher",
        "isca",
        "taa_cases",
        "gain_controls",
        "wtccc_controls",
        "ottawa_controls",
    }
