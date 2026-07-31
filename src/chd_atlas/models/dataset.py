# src/chd_atlas/models/dataset.py
from __future__ import annotations

from collections import Counter
from typing import Final

from pydantic import BaseModel, ConfigDict, Field, model_validator

from chd_atlas.identifiers import AccessionId, ContrastId, Pmid, TaxonId
from chd_atlas.vocab import Archive, Technology

_ACCESSION_PREFIX: Final[dict[Archive, tuple[str, ...]]] = {
    Archive.PRIDE: ("PXD",),
    Archive.GEO: ("GSE",),
    Archive.ARRAYEXPRESS: ("E-",),
    Archive.EGA: ("EGAS", "EGAD"),
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


class Dataset(BaseModel):
    """Top level of ``curation/datasets/<ACCESSION>.yaml``."""

    model_config = ConfigDict(extra="forbid")

    id: AccessionId
    archive: Archive
    technology: Technology
    tissue: str = Field(min_length=1)
    developmental_stage: str = Field(min_length=1)
    organism: TaxonId
    n_samples: int = Field(ge=1)
    licence: str = Field(min_length=1)
    contrasts: list[Contrast] = Field(min_length=1)
    publication: Pmid | None = None

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
        counts = Counter(contrast.id for contrast in self.contrasts)
        duplicate = next((cid for cid, n in counts.items() if n > 1), None)
        if duplicate is not None:
            raise ValueError(f"duplicate contrast {duplicate}")
        return self
