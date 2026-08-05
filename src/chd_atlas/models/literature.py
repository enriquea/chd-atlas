# src/chd_atlas/models/literature.py
from __future__ import annotations

from pydantic import BaseModel, ConfigDict, Field, model_validator

from chd_atlas.duplicates import duplicates
from chd_atlas.identifiers import Doi, PhenotypeId, Pmcid, Pmid
from chd_atlas.vocab import FeaturedTopic, LesionGroup, StudyType


class Publication(BaseModel):
    model_config = ConfigDict(extra="forbid")

    id: Pmid
    title: str = Field(min_length=1)
    journal: str = Field(min_length=1)
    year: int = Field(ge=1900, le=2100)
    authors: list[str] = Field(min_length=1)
    study_type: StudyType
    doi: Doi | None = None
    pmcid: Pmcid | None = None
    own_lab: bool = False
    cohort_size: int | None = Field(default=None, ge=1)
    ancestry: list[str] = Field(default_factory=list)
    # How many comparisons the study reported, for a reader judging an
    # uncorrected p-value off a burden table.
    #
    # Every burden study this atlas will cite is an exome- or genome-wide scan,
    # and their supplements publish raw p-values: Supplementary Data 3 of
    # PMID:42230622 carries `fet.p_value` and no FDR or q-value column at all.
    # So the atlas cannot republish a corrected p, and computing one would be
    # authoring a statistic of its own (D12/D33).
    #
    # Naming the denominator is the honest alternative -- it is a count the
    # study made, not a judgement the atlas adds, and it is what lets a reader
    # apply their own threshold. Measured 2026-08-05 for PMID:42230622:
    # 46,342 + 46,081 + 46,186 = 138,609 gene x consequence-class x stratum
    # comparisons, of which 187 reach a published page. Without this figure a
    # page shows `p = 0.0015` beside a green `definitive` chip and says nothing
    # about the 138,608 other tests it was drawn from.
    #
    # Optional: a family-based linkage study reports no such count, and
    # `pages._method_line` simply omits the sentence when it is absent.
    tests_reported: int | None = Field(default=None, ge=1)


class PublicationFile(BaseModel):
    """Top level of ``curation/publications.yaml``."""

    model_config = ConfigDict(extra="forbid")

    publications: list[Publication] = Field(min_length=1)

    @model_validator(mode="after")
    def ids_are_unique(self) -> PublicationFile:
        found = duplicates(p.id for p in self.publications)
        if found:
            raise ValueError(f"duplicate publication ids: {found}")
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
        found = duplicates(entry.order for entry in self.featured)
        if found:
            raise ValueError(f"duplicate display order: {found}")
        return self


class PhenotypeTerm(BaseModel):
    """One HPO term the atlas cites, cardiac or not.

    `lesion_group` is `None` for a term that is not itself a cardiac lesion --
    an extracardiac feature registered here only so `validate_labels` can check
    its transcription, the same guarantee every cardiac term already gets.
    `None` here is why `validate/referential.py`'s REF010 does not flag such a
    term when an assertion cites it under `extracardiac_features`: that check
    reads `is not None`, so only a term that *is* a cardiac lesion group can
    make an assertion look internally inconsistent.
    """

    model_config = ConfigDict(extra="forbid")

    id: PhenotypeId
    label: str = Field(min_length=1)
    lesion_group: LesionGroup | None = None
    synonyms: list[str] = Field(default_factory=list)


class PhenotypeFile(BaseModel):
    """Top level of ``curation/phenotypes.yaml``."""

    model_config = ConfigDict(extra="forbid")

    phenotypes: list[PhenotypeTerm] = Field(min_length=1)

    @model_validator(mode="after")
    def ids_are_unique(self) -> PhenotypeFile:
        found = duplicates(term.id for term in self.phenotypes)
        if found:
            raise ValueError(f"duplicate phenotype ids: {found}")
        return self
