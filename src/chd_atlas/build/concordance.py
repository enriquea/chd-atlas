# src/chd_atlas/build/concordance.py
"""Whether independent datasets agree about a gene.

**This module computes no statistic.** It counts datasets, and reads the
direction of statistics the studies themselves published. That is arithmetic
over study identities rather than over their results, so design decision D33 --
no pooled statistic, ever -- stands: these cohorts overlap, and combining their
p-values would count the same children twice.

The distinction is the entire reason `cohort_families` exists rather than a
count of studies. A reader handed "3 of 3 datasets agree" will treat it as
replication, and it is only replication if the three describe different people.

Split from `burden.py` rather than added to it: that module loads and shapes
rows, this one derives a view over them, and the two change for different
reasons.
"""

from __future__ import annotations

from collections.abc import Iterable
from enum import StrEnum
from typing import Final

from chd_atlas.build.burden import BurdenRow, shared_cohorts

# The threshold both states are read against. 0.05 is the *studies'* convention
# and not the atlas's: every p-value and corrected p-value here was published
# against it, and this module chooses no threshold of its own. Named rather than
# inlined so the two call sites below cannot drift apart.
ALPHA: Final = 0.05


class FamilyState(StrEnum):
    """What one cohort family says about one gene.

    `NOT_TESTED` is a distinct member and must never be folded into
    `NO_ENRICHMENT`; `family_state` records why.
    """

    CORRECTED = "corrected"
    NOMINAL = "nominal"
    NO_ENRICHMENT = "no_enrichment"
    NOT_TESTED = "not_tested"


def cohort_families(rows: Iterable[BurdenRow]) -> tuple[frozenset[str], ...]:
    """Studies grouped into sets that do not describe independent people.

    A family is a connected component of the graph whose nodes are studies and
    whose edges are shared sample collections. `shared_cohorts` computes the
    edges; this walks them.

    **Transitivity is why this is a graph walk and not a pairwise read of
    `shared_cohorts`.** If A shares a collection with B and B with C, all three
    describe overlapping people even where A and C share nothing directly. A
    reader told "two independent results" would be wrong, and the pairwise
    mapping alone cannot say so.

    Measured 2026-08-05 over the committed corpus: three studies, three
    singleton families, because they draw on disjoint collections. The
    definition is therefore currently indistinguishable from "one family per
    study" -- and it is written this way so that the day a fourth study reuses
    DDD or PCGC, the merge happens with no code change, and nobody is ever told
    that a reused cohort is independent evidence.

    Sorted on the sorted tuple of study ids. Components come out of a set walk,
    and `encode_json`'s `sort_keys` orders dict keys only, so an array derived
    from a set must be sorted by its caller or the build stops being
    reproducible.
    """
    materialised = list(rows)
    parent: dict[str, str] = {row.study: row.study for row in materialised}

    def find(study: str) -> str:
        # Path halving: the chains here are two or three long, so this is for
        # the property rather than the speed -- `find` must be idempotent for
        # the grouping pass below to be stable.
        while parent[study] != study:
            parent[study] = parent[parent[study]]
            study = parent[study]
        return study

    for left, right in shared_cohorts(materialised):
        parent[find(left)] = find(right)

    grouped: dict[str, set[str]] = {}
    for study in parent:
        grouped.setdefault(find(study), set()).add(study)

    return tuple(
        frozenset(members) for members in sorted(grouped.values(), key=lambda m: sorted(m))
    )


def is_upward(row: BurdenRow) -> bool:
    """Whether the row's own statistic points to enrichment in cases.

    Read from what the study published, in three steps, because the three
    studies do not publish the same columns:

    1. `effect_bound == "unbounded_above"` -- an infinite odds ratio, which
       `allow_nan=False` refuses to publish as a number, so the flag *is* the
       result. 19 published rows.
    2. An effect against 1. Every member of `EFFECT_MEASURES` is a ratio, so one
       comparison covers odds ratios, rate ratios and enrichment ratios alike.
    3. Only where a study published no effect at all -- PMID:34324492's PLINK
       CNV test reports none -- the rates the study *did* publish.

    Step 3 compares two published numbers and authors nothing, but it is the
    closest this module comes to D12's line, so it is documented rather than
    left for a reviewer to infer. Note what it is not: it does not compute a
    ratio, a confidence interval or a p-value, and nothing it returns is
    published as a quantity -- only as a direction that decides a glyph.

    A row with neither an effect nor a usable comparator returns `False`. That
    is the safe default: it means "this does not count as agreement", and a row
    that cannot state its direction should not be counted as supporting.
    """
    if row.effect_bound == "unbounded_above":
        return True
    if row.effect is not None:
        return row.effect > 1
    if row.comparator == "mutation_model" and row.expected_count is not None:
        return row.n_case_carriers > row.expected_count
    if row.n_control_carriers is not None and row.n_controls:
        return (row.n_case_carriers / row.n_cases) > (row.n_control_carriers / row.n_controls)
    return False


def _supports(row: BurdenRow) -> bool:
    """Whether one row is evidence *for* the gene at `ALPHA`.

    `synonymous` is excluded because it is the study's negative control rather
    than a result. A significant synonymous row says that gene's comparison is
    poorly calibrated -- it is a warning about the two rows above it -- so
    counting it as support would invert the meaning of the only calibration
    check the table offers.
    """
    if row.consequence_class == "synonymous":
        return False
    if row.pvalue is None or row.pvalue >= ALPHA:
        return False
    return is_upward(row)


def family_state(rows: Iterable[BurdenRow]) -> FamilyState:
    """The strongest thing one family says about one gene.

    **`NOT_TESTED` is not `NO_ENRICHMENT`, and the distinction prevents real
    harm.** KDM6A is ClinGen `definitive` -- it causes Kabuki syndrome -- and
    shows nothing in either of the two datasets that tested it. Collapsing the
    two states would render that as "0 of 3" beside a green definitive chip, and
    a clinical geneticist reads that as the data contradicting the
    classification. It is not: burden tests at these cohort sizes routinely
    detect nothing for genes with overwhelming family and functional evidence.

    So the tally counts *tested* families, the renderer gives the two states
    different glyphs, and the page carries a fixed caption saying that no
    enrichment here is not evidence against a gene.

    `CORRECTED` outranks `NOMINAL` because fill encodes the correction: one
    glyph carries both statistics, which is what lets a reader see that a result
    is nominally significant and does not survive its own study's correction.
    """
    materialised = list(rows)
    if not materialised:
        return FamilyState.NOT_TESTED

    supporting = [row for row in materialised if _supports(row)]
    if not supporting:
        return FamilyState.NO_ENRICHMENT
    if any(row.pvalue_adjusted is not None and row.pvalue_adjusted < ALPHA for row in supporting):
        return FamilyState.CORRECTED
    return FamilyState.NOMINAL


def gene_concordance(
    rows: Iterable[BurdenRow], families: tuple[frozenset[str], ...]
) -> dict[str, object]:
    """One gene's `independent_datasets` object, JSON-ready.

    `families` is passed in rather than derived here because it is a property of
    the *corpus*, not of the gene: a gene no study tested must still render the
    same family slots as every other gene, or a hole stops looking like a hole.

    Named `independent_datasets` in the payload rather than `concordance`. D12
    says the atlas authors no validity classification of its own, and a
    verdict-shaped key sitting beside a mirrored ClinGen `definitive` invites
    being read as a competing call. The name states what is counted.

    Every key is present on every gene, `families` empty where no study reported
    one -- the rule `burden_payload` follows, for the same reason: an object
    whose shape varies is a trap for a consumer that reads a field off one gene
    and expects it on the next.
    """
    by_study: dict[str, list[BurdenRow]] = {}
    for row in rows:
        by_study.setdefault(row.study, []).append(row)

    entries: list[dict[str, object]] = []
    for family in families:
        member_rows = [row for study in sorted(family) for row in by_study.get(study, ())]
        entries.append({"studies": sorted(family), "state": family_state(member_rows).value})

    states = [entry["state"] for entry in entries]
    return {
        # The denominator, and the reason the tally cannot read as a verdict.
        "tested": sum(1 for state in states if state != FamilyState.NOT_TESTED),
        "enriched": sum(
            1 for state in states if state in {FamilyState.CORRECTED, FamilyState.NOMINAL}
        ),
        "corrected": sum(1 for state in states if state == FamilyState.CORRECTED),
        "families": entries,
    }


def evidence_axes(rows: Iterable[BurdenRow]) -> tuple[tuple[str, str], ...]:
    """The distinct `(variant_class, comparator)` designs in the corpus.

    Passed the *whole* corpus, never one gene's rows. A per-gene derivation
    would give each page a different matrix, and a design nobody ran for that
    gene would vanish rather than render as the hole it is -- which is the same
    confusion between "not tested" and "tested and found nothing" that
    `FamilyState` has four members for.

    Measured 2026-08-05 over the published corpus: three axes -- SNV/indel
    case-control, SNV/indel de novo, CNV deletion case-control. The fourth
    quadrant, CNV de novo, is empty because nobody has published one, and it
    should stay visibly empty: the shape of the literature is itself information,
    and it is the argument for which study to curate next.

    Sorted: derived from a set, and it orders the columns of published HTML.
    """
    return tuple(sorted({(row.variant_class, row.comparator) for row in rows}))
