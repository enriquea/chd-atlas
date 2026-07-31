# tests/unit/test_identifiers.py
import pytest
from pydantic import TypeAdapter, ValidationError

from chd_atlas.identifiers import (
    AccessionId,
    AssertionId,
    FunctionalId,
    HgncId,
    ModelPhenotypeId,
    PhenotypeId,
    Pmid,
    UniprotAccession,
)


@pytest.mark.parametrize(
    ("type_", "value"),
    [
        (HgncId, "HGNC:11604"),
        (Pmid, "PMID:8988165"),
        (PhenotypeId, "HP:0001631"),
        (PhenotypeId, "MONDO:0007259"),
        (ModelPhenotypeId, "MP:0010402"),
        (ModelPhenotypeId, "ZP:0000706"),
        (AssertionId, "CHDA:AST:0000001"),
        (FunctionalId, "CHDA:FUN:0000042"),
        (UniprotAccession, "Q99593"),
        (UniprotAccession, "A0A0B4J2D5"),
        (AccessionId, "PXD012345"),
        (AccessionId, "GSE123456"),
        (AccessionId, "E-MTAB-1234"),
    ],
)
def test_accepts_well_formed_identifiers(type_: type, value: str) -> None:
    assert TypeAdapter(type_).validate_python(value) == value


@pytest.mark.parametrize(
    ("type_", "value"),
    [
        (HgncId, "11604"),
        (HgncId, "hgnc:11604"),
        (Pmid, "8988165"),
        (PhenotypeId, "HP:1631"),
        (PhenotypeId, "OMIM:142900"),
        (AssertionId, "CHDA:AST:1"),
        (AssertionId, "CHDA:FUN:0000001"),
        (UniprotAccession, "NOTANACC"),
        (AccessionId, "PXD12345"),
    ],
)
def test_rejects_malformed_identifiers(type_: type, value: str) -> None:
    with pytest.raises(ValidationError):
        TypeAdapter(type_).validate_python(value)
