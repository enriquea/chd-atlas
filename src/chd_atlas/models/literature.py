# src/chd_atlas/models/literature.py
from __future__ import annotations

from collections import Counter
from collections.abc import Hashable, Iterable

from pydantic import BaseModel, ConfigDict, Field, model_validator

from chd_atlas.identifiers import PhenotypeId, Pmid
from chd_atlas.vocab import FeaturedTopic, LesionGroup, StudyType


def _first_duplicate(values: Iterable[Hashable]) -> Hashable | None:
    counts = Counter(values)
    return next((value for value, count in counts.items() if count > 1), None)


class Publication(BaseModel):
    model_config = ConfigDict(extra="forbid")

    id: Pmid
    title: str = Field(min_length=1)
    journal: str = Field(min_length=1)
    year: int = Field(ge=1900, le=2100)
    authors: list[str] = Field(min_length=1)
    study_type: StudyType
    doi: str | None = None
    pmcid: str | None = None
    own_lab: bool = False
    cohort_size: int | None = Field(default=None, ge=1)
    ancestry: list[str] = Field(default_factory=list)


class PublicationFile(BaseModel):
    """Top level of ``curation/publications.yaml``."""

    model_config = ConfigDict(extra="forbid")

    publications: list[Publication] = Field(min_length=1)

    @model_validator(mode="after")
    def ids_are_unique(self) -> PublicationFile:
        duplicate = _first_duplicate(p.id for p in self.publications)
        if duplicate is not None:
            raise ValueError(f"duplicate publication {duplicate}")
        return self


class FeaturedManuscript(BaseModel):
    model_config = ConfigDict(extra="forbid")

    publication: Pmid
    order: int = Field(ge=1)
    blurb: str = Field(min_length=1)
    topic: FeaturedTopic


class FeaturedFile(BaseModel):
    """Top level of ``curation/featured.yaml``."""

    model_config = ConfigDict(extra="forbid")

    featured: list[FeaturedManuscript] = Field(min_length=1)

    @model_validator(mode="after")
    def display_order_is_unique(self) -> FeaturedFile:
        duplicate = _first_duplicate(entry.order for entry in self.featured)
        if duplicate is not None:
            raise ValueError(f"duplicate display order {duplicate}")
        return self


class PhenotypeTerm(BaseModel):
    model_config = ConfigDict(extra="forbid")

    id: PhenotypeId
    label: str = Field(min_length=1)
    lesion_group: LesionGroup
    synonyms: list[str] = Field(default_factory=list)


class PhenotypeFile(BaseModel):
    """Top level of ``curation/phenotypes.yaml``."""

    model_config = ConfigDict(extra="forbid")

    phenotypes: list[PhenotypeTerm] = Field(min_length=1)

    @model_validator(mode="after")
    def ids_are_unique(self) -> PhenotypeFile:
        duplicate = _first_duplicate(term.id for term in self.phenotypes)
        if duplicate is not None:
            raise ValueError(f"duplicate phenotype {duplicate}")
        return self
