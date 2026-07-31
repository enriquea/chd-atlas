# src/chd_atlas/models/functional.py
from __future__ import annotations

from pydantic import BaseModel, ConfigDict, Field, field_validator

from chd_atlas.identifiers import (
    FunctionalId,
    HgncId,
    ModelPhenotypeId,
    Pmid,
    TaxonId,
)
from chd_atlas.vocab import (
    MODEL_ORGANISMS,
    Perturbation,
    PhenocopyAssessment,
    Zygosity,
)


class FunctionalEvidence(BaseModel):
    """Animal- or cell-model evidence, weighted explicitly in classification."""

    model_config = ConfigDict(extra="forbid")

    id: FunctionalId
    gene: HgncId
    organism: TaxonId
    perturbation: Perturbation
    zygosity: Zygosity
    cardiac_phenotype: list[ModelPhenotypeId] = Field(min_length=1)
    phenocopies_human: PhenocopyAssessment
    rescue_performed: bool
    publication: Pmid

    @field_validator("organism")
    @classmethod
    def organism_is_allowlisted(cls, value: str) -> str:
        if value not in MODEL_ORGANISMS:
            raise ValueError(
                f"{value} is not an allowed model organism; "
                f"extend MODEL_ORGANISMS in vocab.py to add one"
            )
        return value


class FunctionalFile(BaseModel):
    """Top level of ``curation/functional/<SYMBOL>.yaml``."""

    model_config = ConfigDict(extra="forbid")

    functional_evidence: list[FunctionalEvidence] = Field(min_length=1)
