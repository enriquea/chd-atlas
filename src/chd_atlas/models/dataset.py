# src/chd_atlas/models/dataset.py
from __future__ import annotations

from typing import Final

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from chd_atlas.duplicates import duplicates
from chd_atlas.identifiers import AccessionId, ContrastId, Pmid, TaxonId
from chd_atlas.vocab import Archive, DatasetDesign, Technology

_ACCESSION_PREFIX: Final[dict[Archive, tuple[str, ...]]] = {
    Archive.PRIDE: ("PXD",),
    Archive.GEO: ("GSE",),
    Archive.ARRAYEXPRESS: ("E-",),
    Archive.EGA: ("EGAS", "EGAD"),
}

_MASS_SPEC: Final[frozenset[Technology]] = frozenset(
    {Technology.TMT_MS, Technology.LFQ_MS, Technology.PHOSPHOPROTEOMICS}
)
_SEQUENCING: Final[frozenset[Technology]] = frozenset(
    {Technology.BULK_RNASEQ, Technology.SCRNASEQ, Technology.MICROARRAY}
)

# PRIDE is the mass-spectrometry repository; GEO, ArrayExpress and EGA hold
# sequencing and array data. A dataset filed under the wrong one is a curator
# error, not an unusual deposit.
_ARCHIVE_TECHNOLOGIES: Final[dict[Archive, frozenset[Technology]]] = {
    Archive.PRIDE: _MASS_SPEC,
    Archive.GEO: _SEQUENCING,
    Archive.ARRAYEXPRESS: _SEQUENCING,
    Archive.EGA: _SEQUENCING,
}


class Contrast(BaseModel):
    """A named comparison within a dataset.

    A fold change is uninterpretable without knowing what was compared to what,
    so contrasts are modelled explicitly rather than left as a free-text column.
    """

    model_config = ConfigDict(extra="forbid")

    id: ContrastId
    description: str = Field(min_length=1)
    case_group: str = Field(min_length=1)
    control_group: str = Field(min_length=1)
    statistical_method: str = Field(min_length=1)
    software: str = Field(min_length=1)
    covariates: list[str] = Field(default_factory=list)

    @model_validator(mode="after")
    def groups_differ(self) -> Contrast:
        if self.case_group.strip().casefold() == self.control_group.strip().casefold():
            raise ValueError(
                f"case_group and control_group must differ, both are '{self.case_group}'"
            )
        return self


class Stage(BaseModel):
    """One of a dataset's own stage tokens, with its position in development.

    `wpc` is null post-natally, which is not a gap: a post-natal stage is
    outside every cardiac morphogenetic phase by definition, and the page says
    so rather than dropping the row.

    `order` is the curated chronological position, 1-based and ascending, and
    is what every renderer sorts on. It is **required and not defaulted**,
    because `wpc` cannot do this job: it is null for every post-natal stage,
    so a dataset's post-natal tokens have no numeric order at all. Sorting
    those tokens as strings is what published `elderly` second of eight and
    `4 week post conception` after `19 week post conception` -- deterministic,
    reproducible and wrong, for three releases.

    Chosen over sorting by position in the YAML file, which would make
    declaration order silently load-bearing: a curator tidying the block would
    change published output with nothing flagging it.
    """

    model_config = ConfigDict(extra="forbid")

    token: str = Field(min_length=1)
    wpc: float | None = Field(default=None, gt=0)
    order: int = Field(ge=1)


class Dataset(BaseModel):
    """Top level of ``curation/datasets/<ACCESSION>.yaml``."""

    model_config = ConfigDict(extra="forbid")

    id: AccessionId
    archive: Archive
    technology: Technology
    tissue: str = Field(min_length=1)
    developmental_stage: str = Field(min_length=1)
    organism: TaxonId
    # A dataset with contrasts declares at least a case and a control sample.
    # The same `ge=2` also binds a profile dataset, where the rationale is
    # unstated (see D40/D41 for what a profile dataset must declare instead).
    n_samples: int = Field(ge=2)
    licence: str = Field(min_length=1)
    contrasts: list[Contrast] = Field(default_factory=list)
    design: DatasetDesign
    # D40. A set, not a scalar: a dataset with two cardiac tokens would file the
    # second as a comparison organ, and tau would then measure specificity
    # against the heart's own other chamber -- driving it down for exactly the
    # most cardiac genes.
    cardiac_tissues: tuple[str, ...] = ()
    stages: tuple[Stage, ...] = ()
    # D41. Per dataset, never a constant: a floor of 1.0 is a defensible RPKM
    # convention, is meaningless in `nx`, and is actively wrong in `lfq`, where
    # intensities are log-scale and 0 is a middling abundance.
    detection_floor: float | None = None
    floor_source: str | None = None
    quantile_estimator: str | None = None
    publication: Pmid | None = None

    @field_validator("licence")
    @classmethod
    def licence_is_not_blank(cls, value: str) -> str:
        """Same rule, same reason, as `Source.licence`.

        `min_length=1` accepts "   ", which records no provenance at all. A
        dataset redistributes third-party measurements exactly as a source does,
        so the two must not disagree about what counts as a recorded licence.
        """
        if not value.strip():
            raise ValueError("licence must not be blank")
        return value

    @field_validator("cardiac_tissues")
    @classmethod
    def cardiac_tissues_are_not_blank(cls, value: tuple[str, ...]) -> tuple[str, ...]:
        """Same rule as `licence_is_not_blank`: a blank entry names no tissue at all."""
        if any(not tissue.strip() for tissue in value):
            raise ValueError("cardiac_tissues must not contain a blank entry")
        return value

    @model_validator(mode="after")
    def accession_matches_archive(self) -> Dataset:
        prefixes = _ACCESSION_PREFIX[self.archive]
        if not self.id.startswith(prefixes):
            raise ValueError(f"accession {self.id} does not match archive '{self.archive.value}'")
        return self

    @model_validator(mode="after")
    def contrast_ids_are_unique(self) -> Dataset:
        found = duplicates(contrast.id for contrast in self.contrasts)
        if found:
            raise ValueError(f"duplicate contrast ids: {found}")
        return self

    @model_validator(mode="after")
    def technology_matches_archive(self) -> Dataset:
        if self.technology not in _ARCHIVE_TECHNOLOGIES[self.archive]:
            raise ValueError(
                f"technology '{self.technology.value}' is not deposited in "
                f"archive '{self.archive.value}'"
            )
        return self

    @model_validator(mode="after")
    def contrasts_match_design(self) -> Dataset:
        if self.design is DatasetDesign.CONTRAST and not self.contrasts:
            raise ValueError("design 'contrast' requires at least one contrast")
        if self.design is DatasetDesign.PROFILE and self.contrasts:
            raise ValueError("design 'profile' declares no contrasts")
        return self

    @model_validator(mode="after")
    def a_profile_dataset_is_fully_declared(self) -> Dataset:
        """Everything the derivation needs, required where it is used.

        Left optional on a contrast dataset because none of it applies there;
        required here because each field is the input to a published figure,
        and a null would publish the figure computed against a default nobody
        curated.
        """
        if self.design is not DatasetDesign.PROFILE:
            return self
        if not self.cardiac_tissues:
            raise ValueError("cardiac_tissues must not be empty for a profile dataset")
        if not self.stages:
            raise ValueError("stages must not be empty for a profile dataset")
        if self.detection_floor is None or not self.floor_source:
            raise ValueError("a profile dataset declares detection_floor and floor_source")
        if not self.quantile_estimator:
            raise ValueError("a profile dataset declares quantile_estimator")
        return self

    @model_validator(mode="after")
    def stage_tokens_are_unique(self) -> Dataset:
        found = duplicates(stage.token for stage in self.stages)
        if found:
            raise ValueError(f"duplicate stage tokens: {found}")
        return self

    @model_validator(mode="after")
    def cardiac_tissues_are_unique(self) -> Dataset:
        found = duplicates(self.cardiac_tissues)
        if found:
            raise ValueError(f"duplicate cardiac_tissues: {found}")
        return self
