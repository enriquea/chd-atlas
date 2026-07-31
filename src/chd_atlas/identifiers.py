# src/chd_atlas/identifiers.py
"""Constrained string types for community identifiers.

Free-text keys are prohibited anywhere a community identifier exists. Each type
below is a plain ``str`` at runtime, so records stay JSON- and TSV-friendly,
but Pydantic rejects malformed values at load time.
"""

from __future__ import annotations

from typing import Annotated

from pydantic import StringConstraints

HgncId = Annotated[str, StringConstraints(pattern=r"^HGNC:\d+$")]
Pmid = Annotated[str, StringConstraints(pattern=r"^PMID:\d+$")]
PhenotypeId = Annotated[str, StringConstraints(pattern=r"^(HP|MONDO):\d{7}$")]
SequenceOntologyId = Annotated[str, StringConstraints(pattern=r"^SO:\d{7}$")]
ModificationId = Annotated[str, StringConstraints(pattern=r"^MOD:\d{5}$")]
ModelPhenotypeId = Annotated[str, StringConstraints(pattern=r"^(MP|ZP):\d{7}$")]
TaxonId = Annotated[str, StringConstraints(pattern=r"^NCBITaxon:\d+$")]

AssertionId = Annotated[str, StringConstraints(pattern=r"^CHDA:AST:\d{7}$")]
FunctionalId = Annotated[str, StringConstraints(pattern=r"^CHDA:FUN:\d{7}$")]

# Official UniProt accession grammar, with an optional isoform suffix.
UniprotAccession = Annotated[
    str,
    StringConstraints(
        pattern=r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-\d+)?$"
    ),
]

# PRIDE (PXD######), GEO (GSE#####), ArrayExpress (E-XXXX-###), EGA (EGAS###########).
AccessionId = Annotated[
    str,
    StringConstraints(
        pattern=r"^(PXD\d{6}|GSE\d+|E-[A-Z]{4}-\d+|EGA[SD]\d{11})$"
    ),
]

ContrastId = Annotated[str, StringConstraints(pattern=r"^[a-z0-9]+(_[a-z0-9]+)*$")]
