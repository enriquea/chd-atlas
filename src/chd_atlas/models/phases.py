# src/chd_atlas/models/phases.py
"""The cardiac morphogenetic phases a developmental stage can fall in.

This is the first developmental-biology claim in the atlas, and it is
transcribed rather than authored: every phase carries the source that defines
its boundaries. `chd_scope.yaml` is the precedent (`cohorts.yaml` is not --
measured 2026-08-14, it carries no attribution at all).

Intervals are half-open, `[start_wpc, end_wpc)`, so a stage sitting exactly on
a boundary falls in exactly one phase. Under closed intervals "8 wpc" would
belong to both the phase ending at 8 and the one starting there, and an
overlap check could not then distinguish a real overlap from adjacency.
"""

from __future__ import annotations

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from chd_atlas.duplicates import duplicates


class CardiacPhase(BaseModel):
    model_config = ConfigDict(extra="forbid")

    id: str = Field(min_length=1)
    label: str = Field(min_length=1)
    start_wpc: float = Field(gt=0)
    end_wpc: float = Field(gt=0)

    @model_validator(mode="after")
    def interval_is_non_empty(self) -> CardiacPhase:
        if self.end_wpc <= self.start_wpc:
            raise ValueError(f"phase {self.id}: end_wpc must exceed start_wpc")
        return self


class CardiacPhaseFile(BaseModel):
    """Top level of ``curation/cardiac_phases.yaml``."""

    model_config = ConfigDict(extra="forbid")

    attributed_to: str = Field(min_length=1)
    citation: str = Field(min_length=1)
    phases: list[CardiacPhase] = Field(default_factory=list)

    @field_validator("attributed_to")
    @classmethod
    def attributed_to_is_not_blank(cls, value: str) -> str:
        """Same rule and reason as `Dataset.licence_is_not_blank`.

        `min_length=1` accepts "   ", which records no provenance at all --
        and provenance is the entire justification for a developmental-biology
        claim reaching this atlas. The atlas transcribes embryology; it does
        not author it, so a blank here is worse than an absent file: it looks
        attributed and is not.
        """
        if not value.strip():
            raise ValueError("attributed_to must not be blank")
        return value

    @field_validator("citation")
    @classmethod
    def citation_is_not_blank(cls, value: str) -> str:
        """Same rule as `attributed_to_is_not_blank`: a blank citation names no source."""
        if not value.strip():
            raise ValueError("citation must not be blank")
        return value

    @model_validator(mode="after")
    def phases_are_unique_and_disjoint(self) -> CardiacPhaseFile:
        found = duplicates(phase.id for phase in self.phases)
        if found:
            raise ValueError(f"duplicate phase ids: {found}")
        ordered = sorted(self.phases, key=lambda phase: (phase.start_wpc, phase.id))
        for earlier, later in zip(ordered, ordered[1:], strict=False):
            if later.start_wpc < earlier.end_wpc:
                raise ValueError(
                    f"phases {earlier.id} and {later.id} overlap: "
                    f"[{earlier.start_wpc}, {earlier.end_wpc}) and "
                    f"[{later.start_wpc}, {later.end_wpc})"
                )
        return self

    def phase_for(self, wpc: float | None) -> CardiacPhase | None:
        """The phase containing `wpc`, or None -- including for a null wpc.

        A null wpc is post-natal, which is outside every cardiac morphogenetic
        phase by definition. Returning None rather than raising is what lets the
        page say "this stage lies outside the curated window" instead of
        dropping the row.
        """
        if wpc is None:
            return None
        for phase in sorted(self.phases, key=lambda item: (item.start_wpc, item.id)):
            if phase.start_wpc <= wpc < phase.end_wpc:
                return phase
        return None
