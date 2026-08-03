# src/chd_atlas/validate/scope.py
"""Checks the curated CHD scope list against the mirrored validity data.

`curation/chd_scope.yaml` is the atlas's own editorial claim about what counts
as congenital heart disease (see `models/scope.py` for why neither mechanical
rule -- expert panel, MONDO `is_a` closure -- works). MONDO is deliberately
*not* pinned as an ontology the way HPO is: the full release would dwarf the
11 MB HPO pin for a resource this module only ever needs to spot-check a
handful of terms in. So a scope term cannot be resolved the way a phenotype
term is in `validate/ontology.py`; instead it is checked against the mirror
-- `mirrors/clingen_gene_validity.tsv` and `mirrors/gencc_submissions.tsv`,
which already carry a `disease`/`disease_label` pair for every term either
authority curates against.

That substitution has a known, accepted blind spot: it does not catch a term
that ClinGen or GenCC still cites but that MONDO has since obsoleted upstream.
Checking against a resolver would catch that; checking against the mirror
cannot, because the mirror is exactly what an obsoleted-but-still-cited term
would still agree with. This is accepted rather than fixed, because the atlas
publishes each scope term with the label its authority gave it -- not with
MONDO's current name for it -- so the claim stays attributable to that
authority either way, obsolete upstream or not.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping

from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.models.scope import ScopeEntry


def validate_scope_terms(
    scope: Iterable[ScopeEntry],
    mirrored: Mapping[str, str] | None,
    location: str,
) -> list[ValidationIssue]:
    """Check every scope term resolves in the mirror, under the label given.

    `mirrored` maps a MONDO term to the label its authority (ClinGen or GenCC)
    published for it; the caller builds this once from both mirror tables
    rather than this function reading them, so it stays a pure check over
    already-loaded data, matching `validate_source_references` and
    `validate_terms` elsewhere in this package.

    `mirrored is None` means neither mirror table could be read at all, not
    that scope has no terms to check against. Checking every entry against an
    empty mapping would report SCP001 once per scope term -- 68 failures
    naming the symptom, none naming the cause -- exactly the cascade
    REF000/SRC000/ONT000 exist to prevent elsewhere in this package. So a
    missing mirror is reported once, as SCP000, and nothing else runs.

    SCP001 (term absent from the mirror) alone cannot catch a copy-pasted
    wrong MONDO id that happens to exist in the mirror under an unrelated
    label -- it would resolve, and the wrong disease would be silently
    admitted. SCP002 exists for that: the label is the only field a human
    reviewing the scope file actually reads, so it is what gets checked
    against the mirror's own label for the same id.
    """
    if mirrored is None:
        return [
            ValidationIssue(
                "SCP000",
                Severity.WARNING,
                location,
                "skipped scope checks: no validity mirror could be read",
            )
        ]

    issues: list[ValidationIssue] = []
    for entry in scope:
        authority_label = mirrored.get(entry.id)
        if authority_label is None:
            issues.append(
                ValidationIssue(
                    "SCP001",
                    Severity.ERROR,
                    location,
                    f"scope term {entry.id} is not in the mirrored ClinGen/GenCC validity data",
                )
            )
            continue
        if entry.label != authority_label:
            issues.append(
                ValidationIssue(
                    "SCP002",
                    Severity.ERROR,
                    location,
                    f"scope term {entry.id} is labelled '{entry.label}' here, "
                    f"but the mirror calls it '{authority_label}'",
                )
            )
    return issues
