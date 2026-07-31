# tests/unit/test_identifiers.py
import pytest
from pydantic import TypeAdapter, ValidationError

from chd_atlas.identifiers import (
    AccessionId,
    AssertionId,
    ContrastId,
    FunctionalId,
    HgncId,
    ModelPhenotypeId,
    ModificationId,
    PhenotypeId,
    Pmid,
    SequenceOntologyId,
    TaxonId,
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
        (TaxonId, "NCBITaxon:9606"),
        (ContrastId, "tof_vs_control"),
        (ContrastId, "a1"),
        (SequenceOntologyId, "SO:0001587"),
        (ModificationId, "MOD:00046"),
        (AccessionId, "EGAS00001000123"),
        (AccessionId, "PXD1234567"),
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
        (TaxonId, "9606"),
        (TaxonId, "NCBITaxon:"),
        (ContrastId, "TOF_vs_control"),
        (ContrastId, "_foo"),
        (ContrastId, "foo_"),
        (SequenceOntologyId, "SO:12345"),
        (ModificationId, "MOD:123456"),
        (FunctionalId, "CHDA:AST:0000001"),
        (AccessionId, "EGAS0001"),
    ],
)
def test_rejects_malformed_identifiers(type_: type, value: str) -> None:
    with pytest.raises(ValidationError):
        TypeAdapter(type_).validate_python(value)
