# tests/unit/test_dataset_model.py
import pytest
from pydantic import ValidationError

from chd_atlas.models.dataset import _ARCHIVE_TECHNOLOGIES, Contrast, Dataset
from chd_atlas.vocab import Archive, Technology


def _contrast(**overrides: object) -> dict[str, object]:
    base: dict[str, object] = {
        "id": "tof_vs_control",
        "description": "Tetralogy of Fallot right ventricle versus non-CHD control",
        "case_group": "TOF RV",
        "control_group": "Non-CHD RV",
        "statistical_method": "limma-voom",
        "software": "limma 3.58.1",
    }
    base.update(overrides)
    return base


def _dataset(**overrides: object) -> dict[str, object]:
    base: dict[str, object] = {
        "id": "PXD012345",
        "archive": "pride",
        "technology": "phosphoproteomics",
        "tissue": "right ventricle",
        "developmental_stage": "infant",
        "organism": "NCBITaxon:9606",
        "n_samples": 24,
        "licence": "CC-BY-4.0",
        "contrasts": [_contrast()],
    }
    base.update(overrides)
    return base


def test_parses_a_well_formed_dataset() -> None:
    dataset = Dataset.model_validate(_dataset())
    assert dataset.contrasts[0].id == "tof_vs_control"
    assert dataset.publication is None


def test_contrast_id_must_be_snake_case() -> None:
    with pytest.raises(ValidationError):
        Contrast.model_validate(_contrast(id="TOF vs Control"))


def test_dataset_requires_at_least_one_contrast() -> None:
    with pytest.raises(ValidationError):
        Dataset.model_validate(_dataset(contrasts=[]))


def test_dataset_rejects_duplicate_contrast_ids() -> None:
    with pytest.raises(ValidationError, match="duplicate contrast ids:.*tof_vs_control"):
        Dataset.model_validate(_dataset(contrasts=[_contrast(), _contrast()]))


def test_dataset_requires_a_licence() -> None:
    with pytest.raises(ValidationError):
        Dataset.model_validate(_dataset(licence=""))


def test_dataset_rejects_an_accession_that_does_not_match_its_archive() -> None:
    with pytest.raises(
        ValidationError, match="accession GSE123456 does not match archive 'pride'"
    ):
        Dataset.model_validate(_dataset(id="GSE123456"))


def test_archive_technology_map_is_total() -> None:
    """A missing key would raise KeyError instead of a clean ValidationError."""
    assert set(_ARCHIVE_TECHNOLOGIES) == set(Archive)
    covered: set[Technology] = set()
    for technologies in _ARCHIVE_TECHNOLOGIES.values():
        covered |= technologies
    assert covered == set(Technology)


def test_rejects_mass_spec_technology_in_a_sequencing_archive() -> None:
    with pytest.raises(ValidationError, match="phosphoproteomics"):
        Dataset.model_validate(_dataset(id="GSE123456", archive="geo"))


def test_rejects_sequencing_technology_in_pride() -> None:
    with pytest.raises(ValidationError, match="bulk_rnaseq"):
        Dataset.model_validate(_dataset(technology="bulk_rnaseq"))


def test_accepts_sequencing_technology_in_geo() -> None:
    dataset = Dataset.model_validate(
        _dataset(id="GSE123456", archive="geo", technology="bulk_rnaseq")
    )
    assert dataset.technology == "bulk_rnaseq"


def test_dataset_with_contrasts_needs_at_least_two_samples() -> None:
    with pytest.raises(ValidationError):
        Dataset.model_validate(_dataset(n_samples=1))


def test_rejects_a_self_comparing_contrast() -> None:
    with pytest.raises(ValidationError, match="case_group"):
        Contrast.model_validate(_contrast(case_group="RV", control_group="rv "))


def test_reports_every_duplicate_contrast_id() -> None:
    with pytest.raises(ValidationError) as exc:
        Dataset.model_validate(
            _dataset(
                contrasts=[
                    _contrast(id="a"),
                    _contrast(id="a"),
                    _contrast(id="b"),
                    _contrast(id="b"),
                ]
            )
        )
    message = str(exc.value)
    assert "a" in message
    assert "b" in message
