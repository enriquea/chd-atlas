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


# The single place an authority's vocabulary becomes the atlas's. Mirrors store
# what was published; this maps it. There is deliberately no `.get(term,
# default)` anywhere downstream -- a term absent from these dicts is an error,
# because a classification silently coerced to an adjacent rung is a confidence
# nobody asserted.
CLINGEN_CLASSIFICATIONS: Final[dict[str, Classification]] = {
    "Definitive": Classification.DEFINITIVE,
    "Strong": Classification.STRONG,
    "Moderate": Classification.MODERATE,
    "Limited": Classification.LIMITED,
    "Disputed": Classification.DISPUTED,
    "Refuted": Classification.REFUTED,
    "No Known Disease Relationship": Classification.NO_KNOWN_ASSOCIATION,
}

# `Supportive` maps to None rather than to a member. It is GenCC's mapping
# exception for submitters that assert an association without grading evidence,
# so it is not a weaker rung of the ladder and must not be placed on one.
# Callers publish the verbatim term and omit the mapped value, which is honest:
# the submitter said "associated", not "how strongly".
GENCC_CLASSIFICATIONS: Final[dict[str, Classification | None]] = {
    "Definitive": Classification.DEFINITIVE,
    "Strong": Classification.STRONG,
    "Moderate": Classification.MODERATE,
    "Limited": Classification.LIMITED,
    "Disputed Evidence": Classification.DISPUTED,
    "Refuted Evidence": Classification.REFUTED,
    "No Known Disease Relationship": Classification.NO_KNOWN_ASSOCIATION,
    "Supportive": None,
}


class ValidityState(StrEnum):
    """How well curated a gene's validity is, published as a fact.

    `UNCURATED` is not an absence. A consumer cannot tell "no authority has
    assessed this gene" from "the build dropped it" by looking at a missing
    field, and that distinction is exactly the silent evidence loss this project
    treats as its characteristic failure.
    """

    EXPERT_CURATED = "expert_curated"
    SUBMITTER_CURATED = "submitter_curated"
    UNCURATED = "uncurated"


class AtlasCuration(StrEnum):
    """Whether the atlas itself has curated evidence for a gene it publishes.

    D21 publishes a gene on an expert panel's classification, which is not the
    atlas's own work. 22 of the 23 genes published today carry no
    `LesionAssertion` at all, so a consumer that read `headline_confidence`
    alone would take a mirrored ClinGen call for curated atlas content.

    Published as a field rather than left inferable from `assertion_count == 0`
    for the reason `validity_state` is: a consumer filtering a browse list needs
    one key to read, not an arithmetic rule to reimplement, and the two must not
    be able to drift apart.
    """

    CURATED = "curated"
    NOT_YET_CURATED = "not_yet_curated"


class ValiditySource(StrEnum):
    CLINGEN = "clingen"
    GENCC = "gencc"


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


class Redistribution(StrEnum):
    """What the source's licence permits the atlas to republish.

    Read by `validate_source_references`, so a source marked `prohibited` is
    reported rather than silently mirrored.
    """

    PERMITTED = "permitted"
    PERMITTED_WITH_ATTRIBUTION = "permitted_with_attribution"
    NON_COMMERCIAL_ONLY = "non_commercial_only"
    PROHIBITED = "prohibited"


class StudyType(StrEnum):
    TRIO_WES = "trio_wes"
    TRIO_WGS = "trio_wgs"
    FAMILY_LINKAGE = "family_linkage"
    CASE_CONTROL = "case_control"
    FUNCTIONAL = "functional"
    PROTEOMICS = "proteomics"
    TRANSCRIPTOMICS = "transcriptomics"
    REVIEW = "review"
    META_ANALYSIS = "meta_analysis"


class BurdenComparator(StrEnum):
    """What a burden row's statistic compared the case count against.

    This is the discriminator the whole burden schema turns on. Every published
    burden analysis asks the same question -- is this gene hit more often than
    expected? -- and differs only in where "expected" comes from:

    - `CONTROL_COHORT`: an observed rate in sequenced controls (case-control).
    - `MUTATION_MODEL`: a per-gene expectation from a trinucleotide mutation
      rate model, the standard for de novo enrichment in trios, where there is
      no control group at all.
    - `NONE`: a case series. A numerator and a denominator, and nothing to
      compare them with -- so no effect size and no p-value may be published,
      which `validate_burden` enforces rather than trusts.

    A row's design is a property of the *row*, not of the study: one paper
    routinely contributes de novo rows and case-control rows, so this cannot be
    hoisted onto the publication record.
    """

    CONTROL_COHORT = "control_cohort"
    MUTATION_MODEL = "mutation_model"
    NONE = "none"


class EffectMeasure(StrEnum):
    """Which quantity a burden row's `effect` column holds.

    Mandatory whenever an effect is present, because an odds ratio of 3.1 and a
    de novo enrichment of 3.1 are different claims and a column headed "effect"
    holding a bare number silently equates them. `EFFECT_MEASURES` below pins
    which measures each comparator can produce; the renderer additionally must
    never print the number without this label.
    """

    ODDS_RATIO = "odds_ratio"
    ENRICHMENT_RATIO = "enrichment_ratio"
    RATE_RATIO = "rate_ratio"


# Which effect measures each comparator can legitimately yield. A mutation
# model cannot produce an odds ratio -- it has no control odds -- and a control
# cohort cannot produce an "enrichment over expectation" because its expectation
# is observed, not modelled. Read by `validate_burden`; the empty set for `NONE`
# is what makes a case series unable to carry an effect at all.
EFFECT_MEASURES: Final[dict[BurdenComparator, frozenset[EffectMeasure]]] = {
    BurdenComparator.CONTROL_COHORT: frozenset(
        {EffectMeasure.ODDS_RATIO, EffectMeasure.RATE_RATIO}
    ),
    BurdenComparator.MUTATION_MODEL: frozenset({EffectMeasure.ENRICHMENT_RATIO}),
    BurdenComparator.NONE: frozenset(),
}


class EffectBound(StrEnum):
    """Why a row reports no finite effect size though the test ran.

    Fisher's exact test returns an unbounded odds ratio when no control carries
    the variant class, and the published upper confidence bound is likewise
    infinite. Measured 2026-08-04 against Supplementary Data 3 of PMID 42230622:
    927 of 138,609 rows across its three strata report `Infinity`, always in
    `fet.odds_ratio` and `fet.ci_95_upper` together and nowhere else.

    Those rows carry the *strongest* signal in the study -- TAB2, 6 syndromic
    carriers against 0 of 45,082 controls -- so dropping them would discard the
    result the reader most wants. They cannot be published as a number either:
    `encode_json` sets `allow_nan=False`, and `Infinity` is accepted by Python's
    `json.loads` but rejected by `JSON.parse`, so it would break every page
    while the manifest checksum still verified.

    The row therefore carries a null `effect`, a null `ci_high`, and this flag,
    and a page renders "unbounded" rather than a blank or an invented ceiling.

    There is deliberately no `unbounded_below` member. An odds ratio of zero --
    13,430 rows in the same sheet, where no case carries -- is finite, publishes
    as `0.0`, and comes with a finite upper bound, so it needs no flag.
    """

    UNBOUNDED_ABOVE = "unbounded_above"


class CohortStratum(StrEnum):
    """Which slice of a study's cases a burden row counted.

    Separate from `SyndromicStatus`, which records what a curator asserts about
    a gene. This records how a *published analysis* partitioned its cohort, and
    its `ALL` is not that enum's `BOTH`: `ALL` means the undivided case set,
    while `BOTH` means a gene causes isolated and syndromic disease alike.
    """

    ALL = "all"
    SYNDROMIC = "syndromic"
    NONSYNDROMIC = "nonsyndromic"


class VariantClass(StrEnum):
    """What kind of DNA change a burden row counted.

    Kept separate from `ConsequenceClass` rather than merged into one
    "variant category": deletion-and-loss-of-function and SNV-and-loss-of-function
    are both real cells, and merging them would make "all loss-of-function
    evidence for this gene, whatever the variant type" unaskable.
    """

    SNV_INDEL = "snv_indel"
    CNV_DELETION = "cnv_deletion"
    CNV_DUPLICATION = "cnv_duplication"
    SV_OTHER = "sv_other"


class ConsequenceClass(StrEnum):
    """What the counted variants do to the protein.

    `SYNONYMOUS` is a first-class member, not an artifact to be filtered out on
    the way in. Synonymous variants should show no case-control enrichment, so a
    synonymous row is the study's own negative control -- and several published
    genes carry a nominally significant one. Publishing it beside the
    loss-of-function row for the same gene, study and stratum is what lets a
    reader judge how well that comparison was calibrated, which no other CHD
    resource shows. `build/bundles.py` groups on that key for exactly this
    reason; nothing anywhere may drop a row for being synonymous.
    """

    LOF = "lof"
    MISSENSE_DAMAGING = "missense_damaging"
    MISSENSE_ALL = "missense_all"
    SYNONYMOUS = "synonymous"
    ALL_CODING = "all_coding"


class VariantOrigin(StrEnum):
    """Whether the counted variants were required to be de novo.

    Without it, a trio study's de novo loss-of-function count and a
    case-control study's rare-inherited loss-of-function count differ in no
    published column and read as the same measurement.
    """

    DE_NOVO = "de_novo"
    INHERITED = "inherited"
    ANY = "any"


class StatisticalTest(StrEnum):
    """Which test produced a burden row's p-value.

    Deliberately not constrained against `BurdenComparator` the way
    `EFFECT_MEASURES` constrains the measure. A binomial test appears against
    both a control rate and a modelled rate in the published literature, so a
    compatibility table here would encode a rule this project has not measured.
    The one rule that *is* asserted -- a `NONE` comparator carries no test at
    all -- lives in `validate_burden`.
    """

    FISHER_EXACT = "fisher_exact"
    POISSON = "poisson"
    BINOMIAL = "binomial"


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
