# src/chd_atlas/vocab.py
from __future__ import annotations

from collections.abc import Iterable
from enum import StrEnum
from typing import Final


class Classification(StrEnum):
    """ClinGen-style gene-disease validity classification."""

    DEFINITIVE = "definitive"
    STRONG = "strong"
    MODERATE = "moderate"
    LIMITED = "limited"
    DISPUTED = "disputed"
    REFUTED = "refuted"
    NO_KNOWN_ASSOCIATION = "no_known_association"


CLASSIFICATION_RANK: Final[dict[Classification, int]] = {
    Classification.DEFINITIVE: 6,
    Classification.STRONG: 5,
    Classification.MODERATE: 4,
    Classification.LIMITED: 3,
    Classification.DISPUTED: 2,
    Classification.REFUTED: 1,
    Classification.NO_KNOWN_ASSOCIATION: 0,
}


def strongest(classifications: Iterable[Classification]) -> Classification:
    """Return the highest-ranked classification, used to derive headline confidence.

    This does not flag conflicting evidence: DISPUTED and REFUTED still sit on the
    same linear rank as the rest, so a REFUTED classification will not surface if a
    higher-ranked classification is also present for the same gene. Callers that
    need to know whether a gene's evidence is contested must also call
    `has_conflicting_evidence()`.
    """
    items = list(classifications)
    if not items:
        raise ValueError("strongest() requires at least one classification")
    return max(items, key=lambda c: CLASSIFICATION_RANK[c])


CONTESTED: Final[frozenset[Classification]] = frozenset(
    {Classification.DISPUTED, Classification.REFUTED}
)


def has_conflicting_evidence(classifications: Iterable[Classification]) -> bool:
    """True when a gene carries both a supportive and a contesting classification.

    ClinGen treats disputed and refuted as a separate axis rather than weaker
    rungs of the definitive-to-limited ladder, so they cannot be compared on the
    single rank `strongest()` uses. Callers displaying headline confidence must
    call this alongside `strongest()`, so a contested gene is never presented as
    settled.
    """
    items = set(classifications)
    contesting = items & CONTESTED
    supportive = items - CONTESTED - {Classification.NO_KNOWN_ASSOCIATION}
    return bool(contesting and supportive)


class LesionGroup(StrEnum):
    CONOTRUNCAL = "conotruncal"
    SEPTAL = "septal"
    LVOTO = "lvoto"
    RVOTO = "rvoto"
    AVSD = "avsd"
    HETEROTAXY = "heterotaxy"
    ANOMALOUS_VENOUS_RETURN = "anomalous_venous_return"
    COMPLEX = "complex"
    OTHER = "other"


class Inheritance(StrEnum):
    AD = "AD"
    AR = "AR"
    XL = "XL"
    DE_NOVO = "de_novo"
    SOMATIC = "somatic"
    OLIGOGENIC = "oligogenic"
    UNKNOWN = "unknown"


class Mechanism(StrEnum):
    LOSS_OF_FUNCTION = "loss_of_function"
    GAIN_OF_FUNCTION = "gain_of_function"
    DOMINANT_NEGATIVE = "dominant_negative"
    HAPLOINSUFFICIENCY = "haploinsufficiency"
    UNKNOWN = "unknown"


class SyndromicStatus(StrEnum):
    ISOLATED = "isolated"
    SYNDROMIC = "syndromic"
    BOTH = "both"


class EvidenceClass(StrEnum):
    GENETIC_CASE = "genetic_case"
    GENETIC_SEGREGATION = "genetic_segregation"
    GENETIC_CASE_CONTROL = "genetic_case_control"
    FUNCTIONAL_MODEL = "functional_model"
    EXPRESSION = "expression"
    PROTEOMIC = "proteomic"
    PTM = "ptm"
    REGULATORY = "regulatory"


class EvidenceStrength(StrEnum):
    STRONG = "strong"
    MODERATE = "moderate"
    SUPPORTING = "supporting"
    CONTRADICTORY = "contradictory"


class SourceTier(StrEnum):
    OWN_CURATION = "own_curation"
    THIRD_PARTY_LITERATURE = "third_party_literature"
    DATABASE_MIRROR = "database_mirror"


class StudyType(StrEnum):
    TRIO_WES = "trio_wes"
    TRIO_WGS = "trio_wgs"
    CASE_CONTROL = "case_control"
    FUNCTIONAL = "functional"
    PROTEOMICS = "proteomics"
    TRANSCRIPTOMICS = "transcriptomics"
    REVIEW = "review"
    META_ANALYSIS = "meta_analysis"


class FeaturedTopic(StrEnum):
    GENOMICS = "genomics"
    VARIANTS = "variants"
    EXPRESSION = "expression"
    PROTEOMICS = "proteomics"
    REVIEW = "review"


class Perturbation(StrEnum):
    KNOCKOUT = "knockout"
    KNOCKIN = "knockin"
    CONDITIONAL_KO = "conditional_ko"
    MORPHOLINO = "morpholino"
    CRISPR = "crispr"
    OVEREXPRESSION = "overexpression"
    PATIENT_VARIANT_KNOCKIN = "patient_variant_knockin"


class Zygosity(StrEnum):
    HOMOZYGOUS = "homozygous"
    HETEROZYGOUS = "heterozygous"
    COMPOUND_HET = "compound_het"
    HEMIZYGOUS = "hemizygous"
    NOT_APPLICABLE = "n/a"


class PhenocopyAssessment(StrEnum):
    YES = "yes"
    PARTIAL = "partial"
    NO = "no"
    NOT_ASSESSED = "not_assessed"


class RescueOutcome(StrEnum):
    NOT_ATTEMPTED = "not_attempted"
    NO_RESCUE = "no_rescue"
    PARTIAL_RESCUE = "partial_rescue"
    FULL_RESCUE = "full_rescue"


class Archive(StrEnum):
    PRIDE = "pride"
    GEO = "geo"
    ARRAYEXPRESS = "arrayexpress"
    EGA = "ega"


class Technology(StrEnum):
    BULK_RNASEQ = "bulk_rnaseq"
    SCRNASEQ = "scrnaseq"
    MICROARRAY = "microarray"
    TMT_MS = "tmt_ms"
    LFQ_MS = "lfq_ms"
    PHOSPHOPROTEOMICS = "phosphoproteomics"


MODEL_ORGANISMS: Final[dict[str, str]] = {
    "NCBITaxon:9606": "Homo sapiens",
    "NCBITaxon:10090": "Mus musculus",
    "NCBITaxon:10116": "Rattus norvegicus",
    "NCBITaxon:7955": "Danio rerio",
    "NCBITaxon:8355": "Xenopus laevis",
    "NCBITaxon:8364": "Xenopus tropicalis",
    "NCBITaxon:9031": "Gallus gallus",
    "NCBITaxon:9823": "Sus scrofa",
}
