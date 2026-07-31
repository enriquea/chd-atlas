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
    """Titles differ, so this fails unless the guard keys on id."""
    with pytest.raises(ValidationError, match="duplicate publication"):
        PublicationFile.model_validate(
            {
                "publications": [
                    _publication(),
                    _publication(title="A different paper reusing the same PMID"),
                ]
            }
        )


def test_publication_file_reports_every_duplicate() -> None:
    with pytest.raises(ValidationError) as exc:
        PublicationFile.model_validate(
            {
                "publications": [
                    _publication(id="PMID:1"),
                    _publication(id="PMID:1"),
                    _publication(id="PMID:2"),
                    _publication(id="PMID:2"),
                ]
            }
        )
    message = str(exc.value)
    assert "PMID:1" in message
    assert "PMID:2" in message


def test_featured_file_rejects_duplicate_display_order() -> None:
    """Publications differ, so this fails unless the guard keys on order."""
    first = {"publication": "PMID:8988165", "order": 1, "blurb": "b", "topic": "genomics"}
    second = {**first, "publication": "PMID:9199562"}
    with pytest.raises(ValidationError, match="duplicate display order"):
        FeaturedFile.model_validate({"featured": [first, second]})


def test_featured_file_accepts_one_publication_in_two_slots() -> None:
    first = {"publication": "PMID:8988165", "order": 1, "blurb": "b", "topic": "genomics"}
    second = {**first, "order": 2}
    parsed = FeaturedFile.model_validate({"featured": [first, second]})
    assert len(parsed.featured) == 2


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
    """Labels differ, so this fails unless the guard keys on id."""
    term = {"id": "HP:0001631", "label": "Atrial septal defect", "lesion_group": "septal"}
    other = {**term, "label": "ASD, alternate label"}
    with pytest.raises(ValidationError, match="duplicate phenotype"):
        PhenotypeFile.model_validate({"phenotypes": [term, other]})


def test_publication_accepts_well_formed_doi_and_pmcid() -> None:
    publication = Publication.model_validate(
        _publication(doi="10.1038/ng0197-30", pmcid="PMC1234567")
    )
    assert publication.doi == "10.1038/ng0197-30"
    assert publication.pmcid == "PMC1234567"


def test_publication_rejects_a_malformed_doi() -> None:
    with pytest.raises(ValidationError):
        Publication.model_validate(_publication(doi="not-a-doi"))


def test_publication_rejects_a_malformed_pmcid() -> None:
    with pytest.raises(ValidationError):
        Publication.model_validate(_publication(pmcid="1234567"))
