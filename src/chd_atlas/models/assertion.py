# src/chd_atlas/models/assertion.py
from __future__ import annotations

from datetime import date
from typing import Annotated, Final, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

from chd_atlas.identifiers import (
    AccessionId,
    AssertionId,
    FunctionalId,
    HgncId,
    PhenotypeId,
    Pmid,
)
from chd_atlas.vocab import (
    Classification,
    EvidenceClass,
    EvidenceStrength,
    Inheritance,
    LesionGroup,
    Mechanism,
    SourceTier,
    SyndromicStatus,
)

# Quantitative omics results are always derived from a deposited dataset, so the
# accession is required. `regulatory` is deliberately excluded: a reporter or MPRA
# result can rest on a figure in the paper with nothing deposited.
_DATASET_BACKED_CLASSES: Final[frozenset[EvidenceClass]] = frozenset(
    {EvidenceClass.EXPRESSION, EvidenceClass.PROTEOMIC, EvidenceClass.PTM}
)


class SupplementaryLocator(BaseModel):
    """Points at an exact place inside a supplementary archive."""

    model_config = ConfigDict(extra="forbid")

    kind: Literal["supplementary"] = "supplementary"
    file: str = Field(min_length=1)
    sheet: str | None = None
    row: int | None = Field(default=None, ge=1)


class InTextLocator(BaseModel):
    """Points at a figure, table or section in the article body."""

    model_config = ConfigDict(extra="forbid")

    kind: Literal["in_text"] = "in_text"
    element: str = Field(min_length=1)


Locator = Annotated[SupplementaryLocator | InTextLocator, Field(discriminator="kind")]


class Evidence(BaseModel):
    model_config = ConfigDict(extra="forbid")

    publication: Pmid
    evidence_class: EvidenceClass
    locator: Locator
    summary: str = Field(min_length=1)
    strength: EvidenceStrength
    dataset: AccessionId | None = None
    functional_evidence: FunctionalId | None = None

    @model_validator(mode="after")
    def functional_reference_matches_class(self) -> Evidence:
        is_functional = self.evidence_class is EvidenceClass.FUNCTIONAL_MODEL
        if is_functional and self.functional_evidence is None:
            raise ValueError(
                "evidence_class 'functional_model' requires 'functional_evidence'"
            )
        if not is_functional and self.functional_evidence is not None:
            raise ValueError(
                "'functional_evidence' is only valid for evidence_class 'functional_model'"
            )
        return self

    @model_validator(mode="after")
    def omics_evidence_cites_a_dataset(self) -> Evidence:
        if self.evidence_class in _DATASET_BACKED_CLASSES and self.dataset is None:
            raise ValueError(
                f"evidence_class '{self.evidence_class.value}' requires 'dataset'"
            )
        return self


class GeneDiseaseAssertion(BaseModel):
    model_config = ConfigDict(extra="forbid")

    id: AssertionId
    gene: HgncId
    phenotypes: list[PhenotypeId] = Field(min_length=1)
    lesion_groups: list[LesionGroup] = Field(min_length=1)
    classification: Classification
    inheritance: list[Inheritance] = Field(min_length=1)
    mechanism: Mechanism
    syndromic: SyndromicStatus
    evidence: list[Evidence] = Field(min_length=1)
    source_tier: SourceTier
    curator: str = Field(min_length=1)
    curated_on: date
    last_reviewed: date
    extracardiac_features: list[PhenotypeId] = Field(default_factory=list)
    notes: str | None = None

    @model_validator(mode="after")
    def extracardiac_requires_syndromic(self) -> GeneDiseaseAssertion:
        if self.extracardiac_features and self.syndromic is SyndromicStatus.ISOLATED:
            raise ValueError(
                "'extracardiac_features' cannot be set when syndromic is 'isolated'"
            )
        return self

    @model_validator(mode="after")
    def syndromic_lists_extracardiac_features(self) -> GeneDiseaseAssertion:
        if self.syndromic is SyndromicStatus.SYNDROMIC and not self.extracardiac_features:
            raise ValueError(
                "syndromic assertions must list at least one extracardiac feature"
            )
        return self

    @model_validator(mode="after")
    def review_not_before_curation(self) -> GeneDiseaseAssertion:
        if self.last_reviewed < self.curated_on:
            raise ValueError("'last_reviewed' cannot precede 'curated_on'")
        return self


class AssertionFile(BaseModel):
    """Top level of ``curation/assertions/<SYMBOL>.yaml``."""

    model_config = ConfigDict(extra="forbid")

    assertions: list[GeneDiseaseAssertion] = Field(min_length=1)
