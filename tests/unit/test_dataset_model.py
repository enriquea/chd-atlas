# tests/unit/test_dataset_model.py
import pytest
from pydantic import ValidationError

from chd_atlas.models.dataset import Contrast, Dataset


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
    with pytest.raises(ValidationError, match="duplicate contrast"):
        Dataset.model_validate(_dataset(contrasts=[_contrast(), _contrast()]))


def test_dataset_requires_a_licence() -> None:
    with pytest.raises(ValidationError):
        Dataset.model_validate(_dataset(licence=""))


def test_dataset_rejects_an_accession_that_does_not_match_its_archive() -> None:
    with pytest.raises(ValidationError, match="accession GSE123456 does not match archive 'pride'"):
        Dataset.model_validate(_dataset(id="GSE123456"))
