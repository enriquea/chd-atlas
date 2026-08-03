# src/chd_atlas/models/functional.py
from __future__ import annotations

from typing import Final

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

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
    RescueOutcome,
    Zygosity,
)

# Morpholino knockdown needs externally-developing, optically accessible embryos.
# It is not a mammalian technique.
#
# Zebrafish, both Xenopus species, and chick. The chick is here because the rule
# above admits it — the embryo develops in ovo and is reachable through a window
# in the shell — and because in-ovo morpholino electroporation is a standard
# technique in avian cardiac development, particularly for cardiac neural crest
# and outflow-tract studies. Gallus gallus is already an allowed model organism,
# so omitting it here rejected a legitimate record: the costlier failure of the
# two, because it blocks curation rather than admitting bad data.
_MORPHOLINO_ORGANISMS: Final[frozenset[str]] = frozenset(
    {"NCBITaxon:7955", "NCBITaxon:8355", "NCBITaxon:8364", "NCBITaxon:9031"}
)

# These perturbations create a stable, heritable genotype, so zygosity is always
# knowable and ClinGen scoring needs it stated.
_GERMLINE_PERTURBATIONS: Final[frozenset[Perturbation]] = frozenset(
    {
        Perturbation.KNOCKOUT,
        Perturbation.KNOCKIN,
        Perturbation.CONDITIONAL_KO,
        Perturbation.PATIENT_VARIANT_KNOCKIN,
    }
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
    rescue_outcome: RescueOutcome
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

    @model_validator(mode="after")
    def morpholino_requires_a_suitable_organism(self) -> FunctionalEvidence:
        if (
            self.perturbation is Perturbation.MORPHOLINO
            and self.organism not in _MORPHOLINO_ORGANISMS
        ):
            raise ValueError(
                f"morpholino knockdown is not valid in {self.organism}; "
                f"it requires an externally-developing embryo"
            )
        return self

    @model_validator(mode="after")
    def germline_perturbation_states_zygosity(self) -> FunctionalEvidence:
        if (
            self.perturbation in _GERMLINE_PERTURBATIONS
            and self.zygosity is Zygosity.NOT_APPLICABLE
        ):
            raise ValueError(
                f"{self.perturbation} is a germline perturbation and must state zygosity"
            )
        return self


class FunctionalFile(BaseModel):
    """Top level of ``curation/functional/<SYMBOL>.yaml``."""

    model_config = ConfigDict(extra="forbid")

    functional_evidence: list[FunctionalEvidence] = Field(min_length=1)
