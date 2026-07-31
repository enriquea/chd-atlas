# tests/unit/test_literature_models.py
import pytest
from pydantic import ValidationError

from chd_atlas.models.literature import (
    FeaturedFile,
    FeaturedManuscript,
    PhenotypeFile,
    PhenotypeTerm,
    Publication,
    PublicationFile,
)


def _publication(**overrides: object) -> dict[str, object]:
    base: dict[str, object] = {
        "id": "PMID:8988165",
        "title": "Mutations in human TBX5 cause limb and cardiac malformation",
        "journal": "Nature Genetics",
        "year": 1997,
        "authors": ["Basson CT", "Bachinsky DR"],
        "study_type": "case_control",
        "own_lab": False,
    }
    base.update(overrides)
    return base


def test_parses_a_well_formed_publication() -> None:
    publication = Publication.model_validate(_publication())
    assert publication.own_lab is False
    assert publication.cohort_size is None


def test_own_lab_flag_marks_enrique_audain_papers() -> None:
    publication = Publication.model_validate(_publication(own_lab=True))
    assert publication.own_lab is True


def test_rejects_an_implausible_year() -> None:
    with pytest.raises(ValidationError):
        Publication.model_validate(_publication(year=1799))


def test_publication_file_rejects_duplicate_pmids() -> None:
    with pytest.raises(ValidationError, match="duplicate publication"):
        PublicationFile.model_validate(
            {"publications": [_publication(), _publication()]}
        )


def test_featured_file_rejects_duplicate_display_order() -> None:
    entry = {
        "publication": "PMID:8988165",
        "order": 1,
        "blurb": "Founding TBX5 paper.",
        "topic": "genomics",
    }
    with pytest.raises(ValidationError, match="duplicate display order"):
        FeaturedFile.model_validate({"featured": [entry, dict(entry)]})


def test_featured_manuscript_requires_a_blurb() -> None:
    with pytest.raises(ValidationError):
        FeaturedManuscript.model_validate(
            {"publication": "PMID:8988165", "order": 1, "blurb": "", "topic": "genomics"}
        )


def test_phenotype_term_carries_a_lesion_group() -> None:
    term = PhenotypeTerm.model_validate(
        {"id": "HP:0001631", "label": "Atrial septal defect", "lesion_group": "septal"}
    )
    assert term.lesion_group == "septal"


def test_phenotype_file_rejects_duplicate_terms() -> None:
    term = {"id": "HP:0001631", "label": "Atrial septal defect", "lesion_group": "septal"}
    with pytest.raises(ValidationError, match="duplicate phenotype"):
        PhenotypeFile.model_validate({"phenotypes": [term, dict(term)]})
