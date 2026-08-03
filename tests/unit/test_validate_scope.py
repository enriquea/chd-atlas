# tests/unit/test_validate_scope.py
"""MONDO is not pinned as an ontology, so scope terms are checked against the
mirror instead of against a resolver: see `validate/scope.py` for why.
"""

from __future__ import annotations

from datetime import date

from chd_atlas.issues import Severity
from chd_atlas.models.scope import ScopeEntry
from chd_atlas.validate.scope import scope_candidates, validate_scope_terms

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


# --- scope_candidates -----------------------------------------------------


def test_net1_finds_a_cardiac_label_and_excludes_admitted_and_non_cardiac_terms() -> None:
    """Net 1 (label keyword) must neither re-report an already-admitted term
    nor flag a term whose label carries no cardiac word.
    """
    mirrored = {
        "MONDO:0001": "atrial septal defect",  # candidate: cardiac, not in scope
        "MONDO:0005453": "congenital heart disease",  # already admitted: excluded
        "MONDO:0002": "diabetes mellitus",  # non-cardiac: excluded
    }
    scope_terms = {"MONDO:0005453"}

    issues = scope_candidates(
        scope_terms=scope_terms,
        mirrored=mirrored,
        diseases_by_gene={},
        in_scope_genes=set(),
        location=_LOCATION,
    )

    assert [i.code for i in issues] == ["SCP003"]
    assert issues[0].severity is Severity.WARNING
    assert "MONDO:0001" in issues[0].message
    assert "MONDO:0005453" not in issues[0].message
    assert "MONDO:0002" not in issues[0].message


def test_net2_finds_a_gene_bridge_candidate_that_net1_misses() -> None:
    """Holt-Oram syndrome's label carries no cardiac keyword, so net 1 is
    blind to it -- this is the measured defect the two-net design fixes. Net
    2 still finds it because TBX5 is already in scope for a different term.
    """
    mirrored = {
        "MONDO:0005453": "congenital heart disease",  # TBX5's in-scope term
        "MONDO:0007732": "Holt-Oram syndrome",  # TBX5's other curation
    }
    scope_terms = {"MONDO:0005453"}
    diseases_by_gene = {"HGNC:11604": {"MONDO:0005453", "MONDO:0007732"}}
    in_scope_genes = {"HGNC:11604"}

    issues = scope_candidates(
        scope_terms=scope_terms,
        mirrored=mirrored,
        diseases_by_gene=diseases_by_gene,
        in_scope_genes=in_scope_genes,
        location=_LOCATION,
    )

    assert [i.code for i in issues] == ["SCP003"]
    assert "net 2" in issues[0].message
    assert "MONDO:0007732" in issues[0].message


def test_at_most_two_issues_with_candidates_sorted() -> None:
    mirrored = {
        "MONDO:0009": "vascular ring",
        "MONDO:0001": "atrial septal defect",
        "MONDO:0005": "aortic coarctation",
    }
    scope_terms: set[str] = set()
    diseases_by_gene = {
        "HGNC:2": {"MONDO:0099", "MONDO:0011"},
        "HGNC:1": {"MONDO:0011"},
    }
    in_scope_genes = {"HGNC:1", "HGNC:2"}

    issues = scope_candidates(
        scope_terms=scope_terms,
        mirrored=mirrored,
        diseases_by_gene=diseases_by_gene,
        in_scope_genes=in_scope_genes,
        location=_LOCATION,
    )

    assert len(issues) == 2
    net1 = next(i for i in issues if "net 1" in i.message)
    net2 = next(i for i in issues if "net 2" in i.message)
    # MONDO:0001 sorts before MONDO:0005, which sorts before MONDO:0009.
    assert net1.message.index("MONDO:0001") < net1.message.index("MONDO:0005")
    assert net1.message.index("MONDO:0005") < net1.message.index("MONDO:0009")
    assert "MONDO:0011" in net2.message
    assert "MONDO:0099" in net2.message
    assert net2.message.index("MONDO:0011") < net2.message.index("MONDO:0099")


def test_nothing_to_report_is_an_empty_list() -> None:
    mirrored = {"MONDO:0005453": "congenital heart disease"}
    scope_terms = {"MONDO:0005453"}
    diseases_by_gene = {"HGNC:11604": {"MONDO:0005453"}}
    in_scope_genes = {"HGNC:11604"}

    issues = scope_candidates(
        scope_terms=scope_terms,
        mirrored=mirrored,
        diseases_by_gene=diseases_by_gene,
        in_scope_genes=in_scope_genes,
        location=_LOCATION,
    )

    assert issues == []
