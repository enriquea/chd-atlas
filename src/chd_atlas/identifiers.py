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
MONDO_PATTERN: Final = r"^MONDO:\d{7}$"
# Named for the same reason as the four above: `mirrors/burden.tsv` keys every
# row on the publication it was transcribed from, so this grammar now has a
# reader on the TSV side as well as in `Pmid` below.
PMID_PATTERN: Final = r"^PMID:\d+$"
# A cohort id in `curation/cohorts.yaml`, and the members of the `;`-joined
# `case_cohorts` / `control_cohorts` columns. Same lowercase-token grammar as
# `ContrastId`; `;` is excluded from the token so the joined form below cannot
# be ambiguous.
COHORT_ID_PATTERN: Final = r"^[a-z0-9]+(_[a-z0-9]+)*$"
COHORT_LIST_PATTERN: Final = r"^[a-z0-9]+(_[a-z0-9]+)*(;[a-z0-9]+(_[a-z0-9]+)*)*$"
# Official UniProt accession grammar, with an optional isoform suffix.
UNIPROT_PATTERN: Final = (
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-\d+)?$"
)

HgncId = NewType("HgncId", Annotated[str, StringConstraints(pattern=HGNC_PATTERN)])
Pmid = NewType("Pmid", Annotated[str, StringConstraints(pattern=PMID_PATTERN)])
CohortId = NewType("CohortId", Annotated[str, StringConstraints(pattern=COHORT_ID_PATTERN)])
Doi = NewType("Doi", Annotated[str, StringConstraints(pattern=r"^10\.\d{4,9}/\S+$")])
Pmcid = NewType("Pmcid", Annotated[str, StringConstraints(pattern=r"^PMC\d+$")])
PhenotypeId = NewType(
    "PhenotypeId", Annotated[str, StringConstraints(pattern=r"^(HP|MONDO):\d{7}$")]
)
# Narrower than `PhenotypeId`, which admits HP as well because a curated
# assertion may cite either. A mirrored gene-disease record names a disease
# entity and nothing else, so accepting an HP term here would let a phenotypic
# feature masquerade as the disease a classification was made against.
# `PhenotypeId` is deliberately left alone: narrowing it changes curated records
# and belongs in its own commit.
DiseaseId = NewType("DiseaseId", Annotated[str, StringConstraints(pattern=MONDO_PATTERN)])
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

AssertionId = NewType("AssertionId", Annotated[str, StringConstraints(pattern=r"^CHDA:AST:\d{7}$")])
FunctionalId = NewType(
    "FunctionalId", Annotated[str, StringConstraints(pattern=r"^CHDA:FUN:\d{7}$")]
)

UniprotAccession = NewType(
    "UniprotAccession", Annotated[str, StringConstraints(pattern=UNIPROT_PATTERN)]
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
