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

Scope is therefore an editorial claim, recorded as one.
"""

from __future__ import annotations

from datetime import date
from typing import Annotated

from pydantic import BaseModel, ConfigDict, Field, StringConstraints, model_validator

from chd_atlas.duplicates import duplicates
from chd_atlas.identifiers import DiseaseId

# Non-blank after stripping: a value of all whitespace passes a bare
# `min_length` and renders as nothing, which is worse than an absent field
# because it looks answered. Mirrors the defect `_cell` in `build/runner.py`
# exists to catch on the mirror side.
_NonBlankStr = Annotated[str, StringConstraints(strip_whitespace=True, min_length=1)]


class ScopeEntry(BaseModel):
    """One disease term admitted into the atlas's scope, and why.

    `label` is transcribed by the curator from the mirrored source, not
    derived from `id`. A later task cross-checks it against the ClinGen and
    GenCC mirrors (SCP002): without that check, a copy-pasted wrong MONDO id
    that happens to exist would silently admit an unrelated disease under a
    label that looks right.
    """

    model_config = ConfigDict(extra="forbid")

    id: DiseaseId
    label: _NonBlankStr
    reason: _NonBlankStr
    admitted_by: _NonBlankStr
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
