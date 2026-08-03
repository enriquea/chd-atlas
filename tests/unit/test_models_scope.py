# tests/unit/test_models_scope.py
"""The scope list is the atlas's editorial claim, so it is shaped like one.

Every entry records who admitted the term and when. That is not bookkeeping:
inclusion is the one thing about a gene the atlas asserts on its own authority
-- every validity classification is mirrored -- so "why is this gene here?" must
have an answer with a name on it.
"""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from chd_atlas.models.scope import ChdScopeFile

ENTRY = {
    "id": "MONDO:0007732",
    "label": "Holt-Oram syndrome",
    "reason": "cardiac septation defects are a defining feature",
    "admitted_by": "Enrique Audain",
    "admitted_on": "2026-08-03",
}


def test_an_entry_requires_a_reason_an_admitter_and_a_date() -> None:
    assert ChdScopeFile.model_validate({"diseases": [ENTRY]}).diseases[0].id == "MONDO:0007732"
    for omitted in ("reason", "admitted_by", "admitted_on"):
        incomplete = {k: v for k, v in ENTRY.items() if k != omitted}
        with pytest.raises(ValidationError):
            ChdScopeFile.model_validate({"diseases": [incomplete]})


def test_a_blank_reason_is_rejected() -> None:
    """An empty reason is worse than no field: it looks answered.

    `min_length` on a stripped value, because a reason of one space passes a
    bare `min_length` and renders as nothing -- the same defect `_cell` in
    `build/runner.py` exists to catch on the mirror side.
    """
    with pytest.raises(ValidationError):
        ChdScopeFile.model_validate({"diseases": [{**ENTRY, "reason": "   "}]})


def test_a_duplicate_term_is_rejected() -> None:
    """Two entries for one term means two reasons, and no rule picks between them."""
    with pytest.raises(ValidationError):
        ChdScopeFile.model_validate({"diseases": [ENTRY, {**ENTRY, "reason": "again"}]})


def test_a_phenotype_term_is_not_a_disease_term() -> None:
    """`DiseaseId` is MONDO-only; an HP term names a feature, not an entity."""
    with pytest.raises(ValidationError):
        ChdScopeFile.model_validate({"diseases": [{**ENTRY, "id": "HP:0001631"}]})


def test_the_committed_scope_list_loads() -> None:
    """A smoke test over the real file, which no unit fixture can stand in for.

    The committed list is the atlas's scope claim. If it stops parsing, every
    gene falls out of scope at once.
    """
    from pathlib import Path

    from ruamel.yaml import YAML

    path = Path(__file__).parent.parent.parent / "curation" / "chd_scope.yaml"
    parsed = ChdScopeFile.model_validate(YAML(typ="safe").load(path.read_text(encoding="utf-8")))
    assert len(parsed.diseases) == 68
