# tests/unit/test_validate_scope.py
"""MONDO is not pinned as an ontology, so scope terms are checked against the
mirror instead of against a resolver: see `validate/scope.py` for why.
"""

from __future__ import annotations

from datetime import date

from chd_atlas.issues import Severity
from chd_atlas.models.scope import ScopeEntry
from chd_atlas.validate.scope import validate_scope_terms

_LOCATION = "curation/chd_scope.yaml"


def _entry(**overrides: object) -> ScopeEntry:
    payload: dict[str, object] = {
        "id": "MONDO:0007732",
        "label": "Holt-Oram syndrome",
        "reason": "cardiac septation defects are a defining feature",
        "admitted_by": "Enrique Audain",
        "admitted_on": date(2026, 8, 3),
    }
    payload.update(overrides)
    return ScopeEntry.model_validate(payload)


def test_a_term_absent_from_the_mirror_is_reported() -> None:
    issues = validate_scope_terms([_entry()], mirrored={}, location=_LOCATION)

    assert [i.code for i in issues] == ["SCP001"]
    assert issues[0].severity is Severity.ERROR
    assert "MONDO:0007732" in issues[0].message


def test_a_disagreeing_label_is_reported_naming_both_labels() -> None:
    mirrored = {"MONDO:0007732": "Holt-Oram syndrome, ClinGen wording"}

    issues = validate_scope_terms([_entry()], mirrored=mirrored, location=_LOCATION)

    assert [i.code for i in issues] == ["SCP002"]
    assert issues[0].severity is Severity.ERROR
    assert "Holt-Oram syndrome" in issues[0].message
    assert "Holt-Oram syndrome, ClinGen wording" in issues[0].message


def test_no_mirror_reports_exactly_one_skip_warning_and_checks_nothing_else() -> None:
    """`mirrored=None` means the validity mirror could not be read at all.

    Checking terms against `{}` would report SCP001 for every one of them,
    burying the single cause -- the same cascade REF000/SRC000/ONT000 exist to
    avoid. So a bad label on an unreadable-mirror run must not also surface as
    SCP002; only the one skip warning is reported.
    """
    entries = [_entry(), _entry(id="MONDO:0005453", label="congenital heart disease")]

    issues = validate_scope_terms(entries, mirrored=None, location=_LOCATION)

    assert len(issues) == 1
    assert issues[0].code == "SCP000"
    assert issues[0].severity is Severity.WARNING


def test_a_matching_term_is_accepted() -> None:
    mirrored = {"MONDO:0007732": "Holt-Oram syndrome"}

    assert validate_scope_terms([_entry()], mirrored=mirrored, location=_LOCATION) == []
