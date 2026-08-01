# src/chd_atlas/models/dataset.py
from __future__ import annotations

from typing import Final

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from chd_atlas.duplicates import duplicates
from chd_atlas.identifiers import AccessionId, ContrastId, Pmid, TaxonId
from chd_atlas.vocab import Archive, Technology

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
    n_samples: int = Field(ge=2)
    licence: str = Field(min_length=1)
    contrasts: list[Contrast] = Field(min_length=1)
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

    @model_validator(mode="after")
    def accession_matches_archive(self) -> Dataset:
        prefixes = _ACCESSION_PREFIX[self.archive]
        if not self.id.startswith(prefixes):
            raise ValueError(
                f"accession {self.id} does not match archive '{self.archive.value}'"
            )
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
