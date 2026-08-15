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
        "design": "contrast",
        "tissue": "right ventricle",
        "developmental_stage": "infant",
        "organism": "NCBITaxon:9606",
        "n_samples": 24,
        "licence": "CC-BY-4.0",
        "contrasts": [_contrast()],
    }
    base.update(overrides)
    return base


def _a_contrast(**overrides: object) -> Contrast:
    base: dict[str, object] = {
        "id": "tof_vs_control",
        "description": "Tetralogy of Fallot right ventricle versus non-CHD control",
        "case_group": "TOF RV",
        "control_group": "Non-CHD RV",
        "statistical_method": "limma-voom",
        "software": "limma 3.58.1",
    }
    base.update(overrides)
    return Contrast(**base)


def _profile_dataset(**overrides: object) -> Dataset:
    base: dict[str, object] = {
        "id": "E-MTAB-6814",
        "archive": "arrayexpress",
        "technology": "bulk_rnaseq",
        "design": "profile",
        "tissue": "Heart",
        "developmental_stage": "embryonic",
        "organism": "NCBITaxon:9606",
        "n_samples": 297,
        "licence": "CC BY 4.0",
        "contrasts": [],
        "cardiac_tissues": ["Heart"],
        "detection_floor": 1.0,
        "floor_source": "source methods, section 4",
        "quantile_estimator": "linear",
        "stages": [{"token": "7wpc", "wpc": 7.0}, {"token": "senior", "wpc": None}],
    }
    base.update(overrides)
    return Dataset(**base)


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
    with pytest.raises(ValidationError, match="accession GSE123456 does not match archive 'pride'"):
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


def test_rejects_a_whitespace_only_licence() -> None:
    """`min_length=1` accepts "   ", which records no provenance at all.

    `Source.licence` already rejects it. A dataset redistributes third-party
    measurements exactly as a source does, so the two must agree on what counts
    as a recorded licence.
    """
    with pytest.raises(ValidationError, match="licence must not be blank"):
        Dataset.model_validate(_dataset(licence="   "))


def test_a_profile_dataset_declares_no_contrast_and_a_contrast_dataset_must() -> None:
    """`design` is explicit so an empty list cannot mean two different things.

    Relaxing `contrasts` to "may be empty" would make a curator who simply
    forgot the contrasts indistinguishable from an abundance atlas that
    correctly has none.
    """
    assert _profile_dataset().contrasts == []

    with pytest.raises(ValidationError, match="design 'profile' declares no contrasts"):
        _profile_dataset(contrasts=[_a_contrast()])

    with pytest.raises(ValidationError, match="design 'contrast' requires at least one contrast"):
        _profile_dataset(design="contrast", contrasts=[])


def test_a_profile_dataset_must_name_a_cardiac_tissue_it_actually_has() -> None:
    """D40: the cardiac token is declared, never matched on the string "heart".

    A set, not a scalar: a dataset with two cardiac tokens would otherwise file
    the second as a comparison organ, and tau would measure specificity against
    the heart's own other chamber.
    """
    assert _profile_dataset().cardiac_tissues == ("Heart",)
    assert _profile_dataset(cardiac_tissues=["Heart", "Atrium"]).cardiac_tissues == (
        "Heart",
        "Atrium",
    )

    with pytest.raises(ValidationError, match="cardiac_tissues must not be empty"):
        _profile_dataset(cardiac_tissues=[])

    with pytest.raises(ValidationError, match="duplicate cardiac_tissues"):
        _profile_dataset(cardiac_tissues=["Heart", "Heart"])


def test_stage_tokens_are_unique_and_wpc_may_be_null_postnatally() -> None:
    stages = {stage.token: stage.wpc for stage in _profile_dataset().stages}
    assert stages["7wpc"] == 7.0
    assert stages["senior"] is None

    with pytest.raises(ValidationError, match="duplicate stage tokens"):
        _profile_dataset(stages=[{"token": "7wpc", "wpc": 7.0}, {"token": "7wpc", "wpc": 9.0}])
