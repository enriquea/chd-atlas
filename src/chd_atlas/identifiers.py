# src/chd_atlas/identifiers.py
"""Constrained string types for community identifiers.

Free-text keys are prohibited anywhere a community identifier exists. Each type
below is a plain ``str`` at runtime, so records stay JSON- and TSV-friendly,
but Pydantic rejects malformed values at load time. Each is defined via
``NewType`` rather than a plain type alias, so mypy also treats them as
mutually non-interchangeable at type-check time even though they share the
same runtime representation.
"""

from __future__ import annotations

from typing import Annotated, Final, NewType

from pydantic import StringConstraints

# Shared with the TSV table schemas in tables.py so each grammar has a single
# definition rather than two literals that can silently drift apart.
HGNC_PATTERN: Final = r"^HGNC:\d+$"
SEQUENCE_ONTOLOGY_PATTERN: Final = r"^SO:\d{7}$"
MODIFICATION_PATTERN: Final = r"^MOD:\d{5}$"

HgncId = NewType("HgncId", Annotated[str, StringConstraints(pattern=HGNC_PATTERN)])
Pmid = NewType("Pmid", Annotated[str, StringConstraints(pattern=r"^PMID:\d+$")])
Doi = NewType("Doi", Annotated[str, StringConstraints(pattern=r"^10\.\d{4,9}/\S+$")])
Pmcid = NewType("Pmcid", Annotated[str, StringConstraints(pattern=r"^PMC\d+$")])
PhenotypeId = NewType(
    "PhenotypeId", Annotated[str, StringConstraints(pattern=r"^(HP|MONDO):\d{7}$")]
)
SequenceOntologyId = NewType(
    "SequenceOntologyId",
    Annotated[str, StringConstraints(pattern=SEQUENCE_ONTOLOGY_PATTERN)],
)
ModificationId = NewType(
    "ModificationId", Annotated[str, StringConstraints(pattern=MODIFICATION_PATTERN)]
)
# HP covers organisms with no dedicated phenotype ontology: human iPSC models, and
# chicken, pig, and Xenopus, for which no MP/ZP equivalent exists.
ModelPhenotypeId = NewType(
    "ModelPhenotypeId", Annotated[str, StringConstraints(pattern=r"^(MP|ZP|HP):\d{7}$")]
)
TaxonId = NewType("TaxonId", Annotated[str, StringConstraints(pattern=r"^NCBITaxon:\d+$")])

AssertionId = NewType(
    "AssertionId", Annotated[str, StringConstraints(pattern=r"^CHDA:AST:\d{7}$")]
)
FunctionalId = NewType(
    "FunctionalId", Annotated[str, StringConstraints(pattern=r"^CHDA:FUN:\d{7}$")]
)

# Official UniProt accession grammar, with an optional isoform suffix.
UniprotAccession = NewType(
    "UniprotAccession",
    Annotated[
        str,
        StringConstraints(
            pattern=(
                r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|"
                r"[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-\d+)?$"
            )
        ),
    ],
)

# PRIDE (PXD######, unpadded past 999999), GEO (GSE#####), ArrayExpress (E-XXXX-###),
# EGA (EGAS###########).
AccessionId = NewType(
    "AccessionId",
    Annotated[
        str,
        StringConstraints(pattern=r"^(PXD\d{6,}|GSE\d+|E-[A-Z]{4}-\d+|EGA[SD]\d{11})$"),
    ],
)

ContrastId = NewType(
    "ContrastId", Annotated[str, StringConstraints(pattern=r"^[a-z0-9]+(_[a-z0-9]+)*$")]
)
