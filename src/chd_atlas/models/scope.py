# src/chd_atlas/models/scope.py
"""What counts as congenital heart disease for this atlas.

This file is the answer to "why is this gene here?", and it exists because
neither obvious mechanical rule works. Measured 2026-08-03 against the committed
ClinGen mirror: filtering on the Congenital Heart Disease expert panel drops
TBX5 (curated by Syndromic Disorders), NOTCH1 (Hereditary Cardiovascular
Disease), CHD7 (Hearing Loss), KMT2D and KDM6A (SCID-CID), FLT4 (Prenatal), ELN
and CTNND1 (Craniofacial) -- every one of them classified Definitive. Filtering
by MONDO `is_a` closure fails from the other side: verified against EBI OLS4,
Holt-Oram syndrome (MONDO:0007732) has 31 hierarchical ancestors and
MONDO:0005453 "congenital heart disease" is not among them. It descends from
MONDO:0005267 "heart disorder", which also subsumes cardiomyopathy, arrhythmia
and acquired disease -- so the closure wide enough to catch Holt-Oram stops
meaning CHD.

Scope is therefore a *selection* among terms external authorities already use --
never a judgement authored here. Until 2026-08-06 every entry in this file
carried `admitted_by: Enrique Audain`, which made the atlas the authority on what
counts as congenital heart disease. It is not one, and the project owner said so
in as many words: *"I am not a clinician... We should take a reputed source that
already do this work for us."*

So `admitted_by` is now a closed vocabulary of **external** authorities and a
person's name is not representable in it. `attributed_to` names the specific
panel or submitter, and `validate/scope.py` checks that authority actually uses
the term in a mirror -- so the attribution is measured rather than asserted.
Measured 2026-08-06 across the 68 committed terms: 13 are used by ClinGen's
Congenital Heart Disease GCEP, 10 by another ClinGen expert panel, 45 by GenCC
submitters, and **0 by no authority at all**. Every term already had an external
warrant; only the field recording it was wrong.
"""

from __future__ import annotations

from datetime import date
from enum import StrEnum
from typing import Annotated

from pydantic import BaseModel, ConfigDict, Field, StringConstraints, model_validator

from chd_atlas.duplicates import duplicates
from chd_atlas.identifiers import DiseaseId

# Non-blank after stripping: a value of all whitespace passes a bare
# `min_length` and renders as nothing, which is worse than an absent field
# because it looks answered. Mirrors the defect `_cell` in `build/runner.py`
# exists to catch on the mirror side.
_NonBlankStr = Annotated[str, StringConstraints(strip_whitespace=True, min_length=1)]


class ScopeAuthority(StrEnum):
    """Who treats a disease term as congenital heart disease. Never this atlas.

    **A closed vocabulary precisely so a curator's name cannot be entered.** The
    previous field was free text and held `Enrique Audain` on all 68 terms; an
    enum makes that unrepresentable rather than discouraged, the same reason
    `PUBLICATION_FLOOR` is a rank floor rather than a set of admissible rungs.

    Ordered strongest first. `CLINGEN_CHD_PANEL` is ClinGen's chartered
    Congenital Heart Disease Gene Curation Expert Panel using the term itself --
    the closest thing to an authority purpose-built for this question.
    `CLINGEN_EXPERT_PANEL` is another chartered ClinGen panel, which matters
    because the syndromic CHD genes are curated there: TBX5 by Syndromic
    Disorders, CHD7 by Hearing Loss, KDM6A and KMT2D by SCID-CID.
    `GENCC_SUBMITTER` is a body submitting through GenCC -- an aggregator, so
    `attributed_to` must name the submitter and never GenCC itself.
    """

    CLINGEN_CHD_PANEL = "clingen_chd_panel"
    CLINGEN_EXPERT_PANEL = "clingen_expert_panel"
    GENCC_SUBMITTER = "gencc_submitter"


class ScopeEntry(BaseModel):
    """One disease term in the atlas's scope, and the authority that scopes it.

    `label` is transcribed by the curator from the mirrored source, not
    derived from `id`. A later task cross-checks it against the ClinGen and
    GenCC mirrors (SCP002): without that check, a copy-pasted wrong MONDO id
    that happens to exist would silently admit an unrelated disease under a
    label that looks right.

    `admitted_by` and `attributed_to` replace a free-text field that named the
    project owner. `reason` survives and still says why the term is congenital
    heart disease -- but it is now a transcription of an external authority's
    position, not an argument this atlas makes.
    """

    model_config = ConfigDict(extra="forbid")

    id: DiseaseId
    label: _NonBlankStr
    reason: _NonBlankStr
    admitted_by: ScopeAuthority
    # The specific panel or submitter, verbatim as the mirror spells it, so
    # SCP005 can check the claim rather than take it. "ClinGen" alone is not an
    # answer for a `GENCC_SUBMITTER`: GenCC aggregates, and its largest in-scope
    # submitter is ClinGen itself (measured 2026-08-06: 111 rows, 109 genes), so
    # a row naming the aggregator hides which body actually made the call.
    attributed_to: _NonBlankStr
    admitted_on: date


class ChdScopeFile(BaseModel):
    """Top level of ``curation/chd_scope.yaml``."""

    model_config = ConfigDict(extra="forbid")

    diseases: list[ScopeEntry] = Field(min_length=1)

    @model_validator(mode="after")
    def ids_are_unique(self) -> ChdScopeFile:
        found = duplicates(entry.id for entry in self.diseases)
        if found:
            raise ValueError(f"duplicate scope term ids: {found}")
        return self
