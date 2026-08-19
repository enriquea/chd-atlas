# src/chd_atlas/models/phases.py
"""The cardiac morphogenetic phases a developmental stage can fall in.

This is the first developmental-biology claim in the atlas, and it is
transcribed rather than authored: every phase carries the source that defines
its boundaries. `chd_scope.yaml` is the precedent (`cohorts.yaml` is not --
measured 2026-08-14, it carries no attribution at all).

Intervals are half-open, `[start_wpc, end_wpc)`, so a stage sitting exactly on
a boundary falls in exactly one phase's own edge case rather than a coin flip.

**A stage may carry SEVERAL phases at once, and this is a fact about the
biology, not a modelling looseness.** An earlier version of this file enforced
a disjoint partition -- one phase per stage -- and that was wrong. Human
cardiac morphogenetic processes run concurrently; a disjoint partition would
force picking one label per stage and silently discarding the others, and the
discarded one may be exactly the lesion-relevant process a reader came for.
Measured from Buijtendijk et al. 2020 (PMID:32048790), against the boundaries
this file transcribes:

| elapsed wpc | processes simultaneously underway |
|---|---|
| 4 | looping, atrial septation, ventricular septation |
| 6 | atrial septation, ventricular septation, outflow tract septation |
| 7 | ventricular septation, outflow tract septation, valve formation |

At 6 wpc an AVSD gene needs "atrial septation" and a TGA gene needs "outflow
tract septation", and both are true at once -- a disjoint model can represent
only one of them. So `CardiacPhaseFile` enforces unique phase **ids** only; it
does not, and must not, enforce that phases are pairwise disjoint.

**Boundaries are computed from HsapDv day-post-fertilization ranges (day /
7.0), never transcribed from a source's ordinal-week prose.** A source stating
"the sixth week of development" means elapsed wpc **5 to 6** (ordinal week N
spans elapsed (N-1) to N), not 6 to 7 -- HsapDv gives Carnegie stage 16 as day
37-42 (5.29-6.00 elapsed wpc), and Buijtendijk et al. 2020 calls CS16 "the
sixth week", which is consistent with the ordinal reading and NOT with reading
"sixth week" as elapsed wpc 6-7. Transcribing the prose directly would shift
every boundary a full Carnegie stage too late, silently and consistently
enough that nothing downstream would catch it -- exactly the kind of error
this file's own docstring exists to name rather than let recur. See
`curation/cardiac_phases.yaml`'s header comment for the same warning at the
point a curator is most likely to re-introduce it.
"""

from __future__ import annotations

from enum import StrEnum

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from chd_atlas.duplicates import duplicates


class EndBasis(StrEnum):
    """Whether a phase's `end_wpc` is transcribed fact or a curator's editorial cap.

    Several of the transcribed processes have no stated end at all -- the
    source describes when they start and simply moves on to the next process
    -- and leaving `end_wpc` undefined is not an option `CardiacPhase` offers
    (`end_wpc` is required, `gt=0`, and must exceed `start_wpc`). Capping such
    a phase at the last Carnegie stage the source discusses (CS23, "the heart
    already resembles the postnatal configuration") is a defensible curator
    choice, but it is a **different kind of fact** from a boundary the source
    itself states, and collapsing the two into one untyped `end_wpc` would let
    a curated guess masquerade as a transcribed measurement -- silently, and
    permanently, since nothing would ever again distinguish them.

    `STATED` means the source names this boundary directly. `CURATOR_CAPPED`
    means the source names no end for this process, and the curator drew a
    line rather than leaving `end_wpc` undefined.
    """

    STATED = "stated"
    CURATOR_CAPPED = "curator_capped"


class CardiacPhase(BaseModel):
    """One named developmental process and the window it occupies.

    `id` is this atlas's own slug, rendered on a gene page by replacing
    underscores with spaces (`build/pages.py::_phase_sentence`) -- distinct
    from `go_id`, which is the transcribed Gene Ontology accession and is
    never rendered directly. Every boundary field below exists so a reader can
    check this atlas's transcription against the source without taking this
    file's word for it: `start_carnegie_stage`/`end_carnegie_stage` name the
    Carnegie stages the source itself cites, `start_hsapdv_id`/`end_hsapdv_id`
    name the HsapDv ontology terms those stages resolve to, and `end_basis`
    says whether `end_wpc` is one of the source's own statements or a
    curator's cap (see `EndBasis`).
    """

    model_config = ConfigDict(extra="forbid")

    id: str = Field(min_length=1)
    go_id: str = Field(min_length=1)
    label: str = Field(min_length=1)
    start_wpc: float = Field(gt=0)
    end_wpc: float = Field(gt=0)
    start_carnegie_stage: str = Field(min_length=1)
    end_carnegie_stage: str = Field(min_length=1)
    start_hsapdv_id: str = Field(min_length=1)
    end_hsapdv_id: str = Field(min_length=1)
    end_basis: EndBasis

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
    def phase_ids_are_unique(self) -> CardiacPhaseFile:
        """Ids must be unique; **intervals need not be disjoint** (see the module docstring).

        An earlier version of this validator also rejected overlapping
        intervals. That was a modelling error, not a data-quality check: human
        cardiac morphogenesis runs several processes at once, so a stage
        legitimately belongs to more than one phase, and rejecting that at
        load time would make it impossible to curate the real vocabulary at
        all rather than catch a curator's mistake.
        """
        found = duplicates(phase.id for phase in self.phases)
        if found:
            raise ValueError(f"duplicate phase ids: {found}")
        return self

    def phases_for(self, wpc: float | None) -> tuple[CardiacPhase, ...]:
        """Every phase containing `wpc`, sorted by `(start_wpc, id)` -- possibly several.

        Empty, never `None`, both for a null wpc and for a real wpc matched by
        no phase: a post-natal stage (`wpc is None`) is outside every cardiac
        morphogenetic phase by definition, and a real wpc outside the curated
        window is a different, real fact. Returning an empty tuple rather than
        raising in either case is what lets a page say "outside the curated
        window" or "post-natal" instead of dropping the row. Sorted
        deterministically so two callers -- and two builds of one commit --
        see phases in the same order regardless of this list's own YAML
        declaration order or of `PYTHONHASHSEED`.
        """
        if wpc is None:
            return ()
        matches = [phase for phase in self.phases if phase.start_wpc <= wpc < phase.end_wpc]
        return tuple(sorted(matches, key=lambda phase: (phase.start_wpc, phase.id)))
