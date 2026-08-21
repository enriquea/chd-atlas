# src/chd_atlas/build/pages.py
"""One HTML page per published gene, plus the browse page over all of them.

Everything a gene page shows has been published as JSON since the validity
backbone landed. `genes/HGNC_11604.json` carries seven mirrored validity
records -- four of them with a report URL, and all four render as links --
alongside a curated assertion whose evidence quote, class, strength and
publication this module renders, under a summary of its lesion groups,
inheritance, mechanism and syndromic status (measured 2026-08-04 against a real
build of the committed corpus). All of it sat at a URL, with a verifying
checksum, and with nothing rendering it. This module is what makes that content
reachable by a person rather than only by a program.

It does not make *all* of it reachable, and this paragraph claimed it did until
the claim was measured. On the same real build, `HP:0001631` appears once in
`genes/HGNC_11604.json` and **zero** times in `genes/HGNC_11604.html`; the
evidence locator's `Abstract` likewise, 1 and 0. The assertion's phenotype terms
and the in-text locator on every evidence item reach no page. Rendering them is
a product decision rather than a defect to be patched, so what is recorded here
is only that the JSON, not the page, is the complete artifact.

Written through `Emitter.write_text`, never `write_json`. `Json`'s union
includes `str`, so `write_json(path, page)` type-checks and then runs the whole
document through `json.dumps`, publishing it quoted and escaped as a JSON string
literal. mypy cannot catch that; the mistake is which encoder to call.

**The uncurated notice is load-bearing.** 91 of the 92 genes published under D21
carry no `LesionAssertion` (measured 2026-08-04 on `genes/index.json`). Their
pages state that the atlas has not curated them rather than omitting the
evidence section, because an absent section is indistinguishable from "the atlas
looked and found nothing" and a reader deciding what a gene means clinically
must not have to infer which. Its exact wording is load-bearing for a second
reason -- it sits one column from a rail that counts curated content this gene
may have -- and `_not_curated` is where that wording and the reason for it are
recorded.

**Neither page kind states a bare grade.** `_SCOPE_RULE` says what admits a gene
to this site at all, `render.grade_legend` says what each rung means in the
grading authority's own words, and `_graded_diseases` names the disease each
ClinGen panel actually graded -- in the browse table's `graded for` column and in
the gene page's rail, one line under the chip. All three read the same `validity`
mapping the bundles were written from, so the two pages cannot name different
diseases for one gene.

Every value interpolated into a page goes through `render.py`, which escapes
text, links and attributes -- except the three bare interpolations in `_rail`,
which escape themselves and say so. Nothing curated or mirrored is reviewed the
way page copy is, and this is the one artifact kind where a `<` in a disease
label is the difference between a string and a script tag. A mirrored disease
label now reaches two more places for that reason: the rail row and the browse
cell, both of which go through `definition_list` and `data_table` respectively.

The research-use notice is `render.document`'s, not this module's, so it is on
every gene page and on the browse page without either builder mentioning it.

Nothing here reads a clock and every value comes from the caller, so two builds
of one commit render byte-identical pages. `build_gene_pages` iterates its
mapping in sorted order for the reason every other loop in this build does.
"""

from __future__ import annotations

import html
from collections.abc import Iterable, Mapping, Sequence
from typing import Final, NamedTuple

from chd_atlas.build.burden import BurdenRow, _sort_key, shared_cohorts
from chd_atlas.build.charts import LogScale, coordinate, marker, polyline, svg_figure
from chd_atlas.build.concordance import FamilyState, family_state
from chd_atlas.build.derive import GeneFacts
from chd_atlas.build.emit import Emitter, Json
from chd_atlas.build.paths import GENE_INDEX_PAGE, gene_bundle_path, gene_page_path
from chd_atlas.build.profiles import (
    EMPTY_EXPRESSION_PROFILE,
    DatasetProfileEntry,
    ExpressionProfile,
    PhaseInfo,
    PhaseOutcome,
    ProfileGap,
    Specificity,
    StageProfileEntry,
    TissueProfileEntry,
)
from chd_atlas.build.render import (
    EVIDENCE_POWER_CAVEAT,
    FILTER_SCRIPT,
    GRADE_CHIP_KIND,
    Cell,
    Link,
    Markup,
    Row,
    chip,
    data_table,
    definition_list,
    document,
    evidence_legend,
    grade_legend,
)
from chd_atlas.build.validity import GeneValidity, agreeing_submitters, uncurated
from chd_atlas.identifiers import HgncId
from chd_atlas.models.assertion import LesionAssertion
from chd_atlas.models.cohort import Cohort
from chd_atlas.models.dataset import Dataset
from chd_atlas.models.literature import Publication
from chd_atlas.models.phases import CardiacPhase, CardiacPhaseFile
from chd_atlas.vocab import AtlasCuration, Classification, ValiditySource, ValidityState

# What the published set *is*, in one sentence, on both page kinds that show it.
#
# No page on this site said it. `docs/data-api.md` states the rule and a reader
# of the HTML never sees it, so a browse page of 23 rows all reading
# `definitive`, under a title reading "CHD Atlas" and a tagline about congenital
# heart disease, said "these 23 genes are definitive for CHD". Measured
# 2026-08-04 against the committed mirrors, that was false for most of them: 10
# of the 23 were qualified by a panel that is not the Congenital Heart Disease
# GCEP -- SCID-CID for KMT2D and KDM6A, Hearing Loss for CHD7, Kidney Cystic and
# Ciliopathy Disorders for NOTCH2 and EVC, Intellectual Disability and Autism
# for ZEB2, Skeletal Disorders for EVC2, Syndromic Disorders for TBX5 and
# TFAP2B, Prenatal for FLT4 -- and 10 were definitive for a disease whose label
# names no cardiac feature at all (Kabuki syndrome 1 and 2, Mowat-Wilson,
# CHARGE, Alagille, Ellis-van Creveld twice, Holt-Oram, TARP, and NR2F2-related
# multiple congenital anomalies).
#
# The 2026-08-06 widening makes the sentence carry more, not less. Measured on
# the widened corpus: 92 genes, of which 76 are admitted by a ClinGen record at
# or above `Limited` and 16 by two or more GenCC submitters with no ClinGen
# record at all. A reader meeting `limited` needs to know both that a chartered
# panel wrote it and that this atlas did not.
#
# One constant for both pages, in the idiom `_NOT_CURATED` already uses: the
# rule is one editorial claim, and two copies of it are two things that drift.
_SCOPE_RULE: Final = (
    '<p class="scope-rule">A gene is published here when an external authority already '
    "treats it as a congenital heart disease gene: either a <strong>ClinGen expert "
    "panel</strong> classifies it <strong>Limited or better</strong> for a disease that "
    "an external authority treats as congenital heart disease, or <strong>two or more "
    "Gene Curation Coalition submitters</strong> independently assert it and no ClinGen "
    "panel disputes it. "
    "<strong>No disease is in scope on this atlas's own judgement</strong>, and neither "
    "is any gene: every scope term names the authority it came from "
    "(<code>curation/chd_scope.yaml</code>), and every gene page names the authority that "
    "admitted the gene. A classification is not the same as definitive for congenital "
    "heart disease: for many of these genes the disease graded is a syndrome of which a "
    "cardiac lesion is one feature, and the disease the panel actually named is shown "
    "beside the gene.</p>"
)

# **Two notices, because the sentence was true of 23 genes and is now true of
# 76.** Until 2026-08-06 every published gene carried a ClinGen grade, so "the
# classification above is an expert panel's" was unconditionally true. The
# widened gate admits 16 genes no panel has graded at all -- measured on the
# committed corpus, all 16 shipped that sentence, telling a reader a commercial
# laboratory's submission was a chartered panel's classification, one screen
# below a legend saying no panel had graded it. A sentence's truth conditions
# travel with the population it describes.
#
# `validity_state` is the discriminator rather than `headline_confidence is
# None`, because the two are equal on the committed corpus (16 and 16) and equal
# figures are one figure to every test.
_NOT_CURATED: Final = (
    "<p>The atlas has <strong>not yet curated</strong> a lesion assertion for this gene. "
    "The classification above is an expert panel's, mirrored with its provenance "
    "intact; no classification on this page is the atlas's own assessment.</p>"
)

_NOT_CURATED_UNGRADED: Final = (
    "<p>The atlas has <strong>not yet curated</strong> a lesion assertion for this gene, "
    "and <strong>no ClinGen expert panel has graded it</strong> for a disease an external "
    "authority treats as congenital heart disease. It is published because the Gene "
    "Curation Coalition submitters named above independently assert it; their "
    "classifications are mirrored with their provenance intact, and no classification on "
    "this page is the atlas's own assessment.</p>"
)

_EM_DASH: Final = "—"

_VALIDITY_HEADERS: Final = (
    "source",
    "submitter or panel",
    "disease",
    "MOI",
    "classification",
    "SOP",
    "date",
    "report",
)

_EVIDENCE_HEADERS: Final = ("class", "strength", "summary", "publication")

# Three of these read as something else on the same page, which review measured
# on 2026-08-05:
#
# * `cohort` held `all cases / syndromic / non-syndromic` -- strata of the case
#   set -- six lines below a provenance paragraph using "cohort" for the sample
#   *collections* (`cnchd`, `ddd`, `ukbb`), which is also what
#   `curation/cohorts.yaml` is named for. `case group` is the stratum.
# * `variants` held `loss-of-function / synonymous`; in a row whose other four
#   cells are numbers, that header invites reading it as a count. It became
#   "variant class" on 2026-08-05 and that was worse: this column renders
#   `consequence_class`, while the method line above the table renders the real
#   `variant_class` ("SNVs and indels"). One phrase meant two things on one page
#   and was wrong in the more prominent place. Raised by review on #17.
#   `test_the_variant_column_is_headed_for_the_column_it_renders` pins it.
# * `cases` and `controls` hold "carriers / total", and nothing said so. The
#   denominator must travel with the numerator -- six carriers is a different
#   claim in 3,876 people than in 45,082 -- and the header is where a reader
#   looks for what a cell means. `_BURDEN_PREAMBLE` states the rule as well.
_BURDEN_HEADERS: Final = (
    "case group",
    "consequence",
    # Not "cases (carriers / n)". The unit is a property of the *row*, not of
    # the column -- one study's table can hold a de novo arm counting mutations
    # and a case-control arm counting alleles -- so it is rendered per cell by
    # `_count` and the header names only what the column is.
    "cases",
    "controls",
    "effect",
    "p",
    "corrected p",
)

# The word for each count unit, as (numerator, denominator suffix). The
# denominator of a de novo row is *trios*, not alleles or people, because a trio
# is what it takes to call one -- so this pairs both halves rather than labelling
# the numerator and leaving the reader to assume the other matches.
_COUNT_UNIT_LABEL: Final[dict[str, tuple[str, str]]] = {
    "individuals": ("carriers", ""),
    "alleles": ("alleles", ""),
    "de_novo_mutations": ("de novo", " trios"),
}

# How a published correction is named in the cell, so a reader is never left to
# guess what "0.99" was corrected against.
_ADJUSTMENT_LABEL: Final[dict[str, str]] = {
    "bonferroni": "Bonferroni",
    "familywise_permutation": "family-wise",
    "benjamini_hochberg": "Benjamini-Hochberg",
}

_STRATUM_LABEL: Final[dict[str, str]] = {
    "all": "all cases",
    "syndromic": "syndromic",
    "nonsyndromic": "non-syndromic",
}

# The four states, as a reader meets them. `not_tested` is deliberately not a
# synonym for `no result`: KDM6A is ClinGen definitive -- it causes Kabuki
# syndrome -- and shows nothing in either dataset that tested it, so a page that
# merged the two would report "0 of 3" beside a green definitive chip and read
# as the data contradicting the classification.
# The design each comparator names, for the matrix column heads. Reuses the
# discriminator the whole burden schema turns on rather than inventing a
# parallel vocabulary for the page.
# What the glyphs mean, and the sentence that stops the column reading as a
# verdict. **The caption is not optional.** KDM6A is ClinGen definitive -- it
# causes Kabuki syndrome -- and shows nothing in either dataset that tested it,
# so its row reads "0 of 2 tested" beside a green `definitive` chip. Without the
# caption a clinical geneticist reads that as the data contradicting the
# classification; with it, they read it as the power statement it is.
# The matrix's own key and caption. **The browse page carried these and the gene
# page did not**, which is the wrong way round: a reader arriving from a search
# lands on the gene page, and KDM6A's matrix there is entirely hollow beside a
# green `definitive` chip with nothing to explain it. Found by review before the
# deploy, on the page where the omission would have done the damage.
_MATRIX_LEGEND: Final = (
    evidence_legend(swatches=True)
    + '<p class="strip-legend">Rows are <strong>independent cohort families</strong>, not '
    "studies: two papers sharing a sample collection describe the same people and appear "
    "once. " + EVIDENCE_POWER_CAVEAT + "</p>"
)

_STRIP_LEGEND: Final = (
    evidence_legend(swatches=False)
    + '<p class="strip-legend">One glyph per <strong>independent cohort family</strong>, not '
    "per study: two papers sharing a sample collection describe the same people and count "
    "once. "
    + EVIDENCE_POWER_CAVEAT
    + " The tally reads <em>enriched</em> of <em>tested</em>: the first "
    "number counts datasets that found something, the second counts datasets that looked. "
    "Neither is a verdict on the gene.</p>"
)

_COMPARATOR_LABEL: Final[dict[str, str]] = {
    "control_cohort": "case-control",
    "mutation_model": "de novo, trios",
    "none": "case series",
}

_STATE_GLYPH: Final[dict[str, str]] = {
    "corrected": "full",
    "nominal": "half",
    "no_enrichment": "none",
    "not_tested": "untested",
}

# What a matrix cell *reads*, as against what it says on hover. The full
# sentence -- "this dataset did not test this gene" -- set the column width for
# the whole table and pushed the cells wide enough to crowd each other. The
# short form keeps the grid tight; `_STATE_TITLE` still carries the sentence on
# the cell's `title`, so nothing is lost to a reader who wants it.
_STATE_CELL: Final[dict[str, str]] = {
    "not_tested": "not tested",
    "no_enrichment": "tested",
}

_STATE_TITLE: Final[dict[str, str]] = {
    "corrected": "enriched, and survives this study's own correction",
    "nominal": "enriched nominally; not after correction, or no correction published",
    "no_enrichment": "tested, no enrichment detected",
    "not_tested": "this dataset did not test this gene",
}

_CONSEQUENCE_LABEL: Final[dict[str, str]] = {
    # Names its two components rather than saying "damaging", so a reader who
    # meets this row first is told what it is the union of before they reach the
    # two rows below it. `_composite_note` says the rest.
    "damaging": "damaging (LOF + missense)",
    "lof": "loss-of-function",
    "missense_damaging": "missense (damaging)",
    "missense_all": "missense (all)",
    "all_coding": "all coding",
    "synonymous": "synonymous",
}

# Spelled out beside every number, because the abbreviation is the difference
# between two claims. See `_effect`.
_MEASURE_LABEL: Final[dict[str, str]] = {
    "odds_ratio": "OR",
    "enrichment_ratio": "enrichment",
    "rate_ratio": "rate ratio",
}

# What the value 1 *means*, which is not the same sentence for every measure --
# and `_forest_caption` printed one of them on all of them.
#
# Measured 2026-08-20 on the built site: 112 of the 137 forest panels, across
# 88 gene pages, read "Axis: **OR** ... the vertical line is **1**, no
# enrichment". For an odds ratio 1 is *no association*; "no enrichment" is the
# de-novo-against-expectation vocabulary belonging to the 25
# `enrichment_ratio` panels. That is exactly the merge `effect_measure` exists
# to prevent (see `_effect`), committed in the one sentence whose job is
# telling a reader how to read the axis.
#
# Keyed on the same token as `_MEASURE_LABEL` and consulted in the same place,
# so the label and its gloss cannot drift apart again: a fourth measure that
# lands here with a label and no gloss gets `_NULL_LINE_DEFAULT`, which is
# true of any ratio, rather than one vocabulary's words applied to another's
# number. `rate_ratio` carries 0 rows today and is written anyway, for the
# reason `_effect` has no branch that omits its measure.
_NULL_LINE_GLOSS: Final[dict[str, str]] = {
    "odds_ratio": "no association",
    "enrichment_ratio": "no enrichment",
    "rate_ratio": "no difference in rate",
}

# Deliberately not a blank and not a refusal. `_effect` renders an unknown
# measure's raw token rather than nothing, on the rule that a missing label is
# worse than an unfamiliar one; the same rule applies to the gloss, and every
# member of `EffectMeasure` is a ratio, for which 1 is the arms agreeing.
_NULL_LINE_DEFAULT: Final = "no difference between the arms"

# The partition columns that reached no page until 2026-08-05. Without them a
# section headed "Rare variant burden" never says what *rare* means and a
# case-control count of rare inherited variants reads exactly like a trio's de
# novo count.
#
# This said "the four partition columns... each is populated on all 187 published
# rows" and named `lesion_group` among them. **`lesion_group` is null on all 290
# published rows and always has been** -- no curated study stratifies by lesion
# group -- and `maf_max` is null on 45 of them, the de novo rows that applied no
# frequency filter. So the claim was false for two of its four columns.
# `_LESION_LABEL` stays empty rather than being deleted: the column is live in
# the schema, and the first study that uses it needs that map filled.
_VARIANT_CLASS_LABEL: Final[dict[str, str]] = {
    "snv_indel": "SNVs and indels",
    "cnv_deletion": "CNV deletions",
    "cnv_duplication": "CNV duplications",
    "sv_other": "other structural variants",
}

_ORIGIN_LABEL: Final[dict[str, str]] = {
    "de_novo": "de novo only",
    "inherited": "inherited only",
    "any": "any inheritance (not a de novo test)",
    # Says what was excluded and what is unknown, because that is the whole
    # content of the value: de novo variants were removed, and most of what
    # remains was never phased. "inherited" would claim transmission the study
    # did not observe for the singleton majority.
    # No internal semicolon: `_method_line` joins its parts with "; " and the
    # origin labels within a part with ", ", so a label containing one produced
    # "de novo excluded; transmission otherwise known or unknown, de novo only"
    # -- which reads as three clauses at two nesting levels and parses as none.
    "transmitted_or_unphased": "de novo excluded (transmission otherwise known or unknown)",
}

_LESION_LABEL: Final[dict[str, str]] = {}

_TEST_LABEL: Final[dict[str, str]] = {
    "fisher_exact": "Fisher's exact test",
    "poisson": "Poisson test",
    "binomial": "binomial test",
}

# Attribution, in the voice `_validity_section` and `landing.py` already use for
# republished content. The `<h2>` alone said "Rare variant burden" while its
# sibling says "Mirrored gene-disease validity", and the `_not_curated` notice
# above it scopes itself to *classifications* -- so nothing between the heading
# and the table said whose numbers these are.
_MIRRORED_NOTICE: Final = (
    "<p>Every count and statistic below is an upstream study's, republished exactly as "
    "published. The atlas computes none of them and adds no assessment of them.</p>"
)

# Rendered only when a gene carries rows from more than one study. Rendered
# unconditionally it asserted "these cohorts overlap" in the present tense on a
# page showing a single table, sending a reader to look for a second study that
# is not there.
# Now the matrix's caption, so it is rendered on every page with a matrix rather
# than only where a gene carries two studies. That change made the old wording
# false: it asserted "these cohorts overlap" in the present tense, which on a
# page showing one study sent a reader hunting for a second -- the very defect
# the original conditional was added to fix. It is a statement of policy now,
# true whether the page shows one study or three, and `shared_cohorts` still
# names an actual overlap where one exists.
_POOLING_NOTICE: Final = (
    "<p>The atlas computes <strong>no pooled statistic across studies</strong>. "
    "Where two studies draw on the same sample collection they describe partly the same "
    "children, so combining their results would count those children twice; any overlap "
    "between the studies below is named where it occurs.</p>"
)

# What the numbers do not say. Every clause was measured; see `_burden_section`.
#
# The absent-cell sentence used to assert a *mechanism*: "no variant of that
# class was seen in either group, so there was nothing to compare and the study
# reported no row for it." That was measured over PMID:42230622 alone -- 0 of
# its 1,192 rows have no carrier on either side -- and it is false for the other
# two studies now rendering tables beside it. PMID:34324492 tests one
# consequence class by construction, and PMID:40127276 **observed 14,364
# synonymous variants** (its own Dataset S4) and still published no synonymous
# row, because its gene-level table covers only damaging classes. A reader told
# "no variant of that class was seen" would conclude that study found none.
#
# So the sentence now says only what is true of every study: an absent cell is
# something the study did not report, and why is the study's own rule.
_BURDEN_PREAMBLE: Final = (
    "<p><code>cases</code> and <code>controls</code> each show <em>observed / total</em>, "
    "and <strong>every cell names what it counted</strong> &mdash; carriers, alleles, or "
    "de novo mutations. The three are not interchangeable: a study counting alleles "
    "counts a person with two qualifying variants twice, and its denominator is roughly "
    "twice the number of people sequenced. "
    "A consequence class or case group with no row is one the study "
    "<strong>did not report</strong>, which is not a null result and does not always mean "
    "no such variant was seen &mdash; studies differ in which classes their gene-level "
    "tables cover.</p>"
)

# Rendered only where a synonymous row is actually on the page. Rendered
# unconditionally it told a reader to weigh a negative control that 2 of the 3
# published studies do not provide: 69 of the 290 published rows are synonymous
# and every one comes from PMID:42230622.
_SYNONYMOUS_NOTICE: Final = (
    "<p>The <strong>synonymous</strong> row is that study's own negative control &mdash; "
    "synonymous variants should show no enrichment. Read it on the same uncorrected "
    "scale as the rows above it: a <em>strongly</em> enriched synonymous row is a "
    "warning about that gene's comparison, while a nominally significant one is what a "
    "scan of this size produces by chance. <strong>Only a study that publishes one has "
    "it</strong>; a table with no synonymous row offers no such check.</p>"
)

_CITATION_HEADERS: Final = ("id", "title", "year")

# `graded for` sits immediately after `confidence`, because it is what qualifies
# it: the two columns are one claim split in two, and a reader whose eye stops at
# the grade must meet the disease next rather than three columns later.
_BROWSE_HEADERS: Final = (
    "gene",
    "symbol",
    "confidence",
    "graded for",
    # Immediately after the claim, because it is what qualifies it: a reader
    # whose eye stops at `definitive` should meet the evidence next rather than
    # three columns later.
    #
    # Headed for what it *is*, not for what it might mean. "replicated in" was
    # the shorter candidate and was rejected: for a gene showing 0 of 2 it reads
    # as "not replicated", which is a verdict the data do not support -- KDM6A
    # causes Kabuki syndrome and shows nothing in either dataset that tested it.
    # D12 says the atlas authors no validity classification, and a
    # verdict-shaped header beside a mirrored ClinGen `definitive` invites being
    # read as a competing one. This names the evidence type and the axis and
    # leaves the reading to the reader.
    "burden across studies",
    "validity",
    "atlas curation",
    "burden rows",
    "lesion groups",
)


def _pubmed(pmid: str) -> str:
    """The canonical article URL for a `PMID:` identifier.

    A literal rather than a mirrored field: `Publication` records a DOI and a
    PMCID but no PubMed URL, and this is the one link on the page the atlas
    constructs rather than republishes.
    """
    return f"https://pubmed.ncbi.nlm.nih.gov/{pmid.removeprefix('PMID:')}/"


def _graded_diseases(gene_validity: GeneValidity, grade: Classification | None) -> tuple[str, ...]:
    """The in-scope diseases a ClinGen panel graded this gene at `grade`.

    **Keyed on the gene's own headline, not on `DEFINITIVE`.** It was the latter
    until 2026-08-06, when the gate widened past ClinGen `Definitive` and 69 of
    the 92 published genes stopped having one: a function that only ever looks
    for definitive records returns nothing for them, and the caller renders an em
    dash beside a `limited` chip -- the disease the panel actually graded silently
    absent from the one row that exists to name it. That is this project's
    characteristic failure in page copy rather than in a join.

    `grade` is `GeneFacts.headline_confidence`, which `derive.gene_facts` derives
    from ClinGen records alone, so what comes back is always the disease behind
    the chip the reader is looking at. `None` -- a gene admitted on GenCC
    agreement, which no panel graded -- returns empty, and the caller names the
    submitters instead.

    A GenCC submitter's own grade is deliberately never here. It admits a gene
    only by agreement with a second submitter (D21) and is never the headline, so
    naming its disease would qualify a chip with a record that did not produce
    it.

    Returns labels rather than MONDO ids: this is the string a reader reads
    beside the grade, and the id is already on the gene page in
    `_validity_section`, which renders both.

    **Sorted, and de-duplicated through a `set` first.** A gene may hold more
    than one record at its headline grade -- two panels grading it for two
    in-scope diseases, or one disease reached through two rows -- and a `set` of
    strings iterates in an order that varies with `PYTHONHASHSEED`, which would
    make the browse page's bytes and therefore its checksum differ between two
    builds of one commit. Measured 2026-08-06 against the committed mirrors: all
    76 panel-graded genes have exactly one distinct label at their headline
    grade, so the multi-label case is latent and its ordering is guarded by
    fixture alone
    (`test_a_gene_definitive_for_two_in_scope_diseases_names_both_in_a_fixed_order`).
    That name was wrong here until 2026-08-06 -- it cited a
    `..._graded_for_two_...` test that has never existed, which is the one thing
    a cross-reference must not do: it sends the next reader looking for a guard
    they will not find, and deleting the real one then costs nothing.
    """
    if grade is None:
        return ()
    return tuple(
        sorted(
            {
                record.disease_label
                for record in gene_validity.records
                if record.source is ValiditySource.CLINGEN and record.classification is grade
            }
        )
    )


def _not_curated(fact: GeneFacts) -> str:
    """The notice an uncurated gene page carries, and what it must not deny.

    `atlas_curation` is derived from `LesionAssertion`s alone -- `curated` iff
    `assertion_count > 0`, `derive.gene_facts` -- while `_rail` counts
    `functional_count` beside it from a population `referential.py` never
    requires an assertion to cite. The two genuinely differ, so a gene can carry
    functional records the atlas curated and still resolve to
    `NOT_YET_CURATED`.

    That is why this paragraph names the *lesion assertion* as what is absent
    and narrows its second clause to classifications. The sentence it replaces
    said "not yet curated evidence for this gene ... nothing on this page is the
    atlas's own assessment", which such a gene's own rail contradicts one column
    away by counting `functional records | 3`. Denying curated work on the page
    that is rendering it is this project's characteristic failure wearing prose
    instead of a dropped join.

    `atlas_curation`'s meaning is deliberately untouched: `curated iff
    assertion_count > 0` is tested, documented in `vocab.AtlasCuration` and
    published in every bundle and index row. What was wrong was the sentence,
    not the field.

    Latent on the committed corpus, which has no `curation/functional/`
    directory at all, so `corpus.functional` is empty and every published gene
    counts zero (measured 2026-08-04). The second paragraph is reachable only
    from a fixture until a functional record is curated.

    Measured 2026-08-04 by restoring the old paragraph and returning it
    unconditionally: 1 failed, 621 passed --
    `test_an_uncurated_notice_never_denies_functional_records_the_rail_counts`
    is its unique killer, and
    `test_an_uncurated_gene_page_says_the_atlas_has_not_curated_it` does not
    notice, because "not yet curated" is true of both wordings.

    Two further mutants, measured the same day, because the two branches are
    guarded separately: returning `""` for the `functional_count == 0` case
    fails only `test_an_uncurated_gene_page_says_the_atlas_has_not_curated_it`
    (1 failed, 621 passed), and returning `""` unconditionally fails both tests
    (2 failed, 620 passed). Neither branch rests on the other's guard.
    """
    notice = (
        _NOT_CURATED
        if fact.validity_state is ValidityState.EXPERT_CURATED
        else _NOT_CURATED_UNGRADED
    )
    if fact.functional_count == 0:
        return notice
    noun = "record" if fact.functional_count == 1 else "records"
    return notice + (
        f"<p>The {fact.functional_count} functional {noun} counted beside this notice "
        f"<em>are</em> the atlas's own curation, and are published in this gene's JSON.</p>"
    )


def _report_link(url: str | None) -> Cell:
    """A link only for an `http`/`https` report URL; an em dash for anything else.

    Neither mirror schema constrains this column, and the failure is not
    hypothetical: `mirrors/gencc_submissions.tsv` carries one row -- SHOX,
    `MONDO:0009588` Langer mesomelic dysplasia, submitted by Ambry Genetics --
    whose `report_url` is the literal string `Pseudoautosomal region, recessive`
    (measured 2026-08-04, the only such row of the file). That disease is out of
    CHD scope today, so the row reaches no page; the same column would carry a
    `javascript:` URL into an `href` just as readily, and `render.py` escapes an
    attribute value without sanitising the scheme it names.

    Fixed here rather than in `tables.py` on purpose. The mirror is a faithful
    copy of what upstream published and has to stay one -- rejecting the row
    would drop a real GenCC submission over a field no reader needs -- so the
    scheme is decided at the render boundary, where a URL becomes a link.

    A guard on a bypassed gate, in the same idiom as `encode_json`'s
    `allow_nan=False` and `variants.py`'s shard checks: nothing upstream
    promises this, so the page refuses rather than publishes.
    """
    if url and url.startswith(("http://", "https://")):
        return Link(text="open", href=url)
    return _EM_DASH


def _rail(
    gene: str,
    symbol: str,
    fact: GeneFacts,
    diseases: Sequence[str],
    submitters: Sequence[str],
    burden_rows: int,
) -> str:
    """The summary column, populated for every gene whether curated here or not.

    **The chip class is looked up in `GRADE_CHIP_KIND`, keyed on the grade's own
    value.** Until 2026-08-06 it was `"definitive" if headline is DEFINITIVE else
    ""`, which was right while the gate admitted nothing else and became wrong
    the moment it widened: 53 of the 92 genes published carry `strong`,
    `moderate` or `limited`, and a bare `.chip` renders all three in the same
    grey as `not classified`. A reader cannot tell a panel's weakest supportive
    rung from no panel having looked, which is the distinction the widened gate
    turns on. `render.grade_legend` is the other half -- the colours are ordinal
    and mean nothing without ClinGen's own definitions beside them.

    **The first row's label carries the grade, so the label and the disease
    cannot come apart.** It read `definitive for` unconditionally until the same
    date, and `_graded_diseases` was returning the *definitive* records, so on a
    `limited` gene the row was headed `definitive for` beside a `limited` chip.
    Now the label is the headline's own value -- `limited for: atrial septal
    defect` -- and for a gene no panel graded it becomes `admitted by`, naming
    the submitters whose agreement admitted it. Both come from `fact`, so the
    chip and the row below it are one claim.

    It stays first in the list for the reason the browse column sits beside
    `confidence`: the chip one line above says a bare grade, and on this site's
    title and tagline that reads as "definitive for congenital heart disease".
    For KMT2D it means Kabuki syndrome 1. The row is what stops the chip standing
    alone on the first screen.
    """
    headline = fact.headline_confidence
    chips = [
        chip(
            headline.value if headline is not None else "not classified",
            kind=GRADE_CHIP_KIND[headline.value] if headline is not None else "ungraded",
        )
    ]
    if fact.has_conflicting_evidence:
        chips.append(chip("conflicting evidence", kind="warn"))
    if fact.has_source_discordance:
        chips.append(chip("sources disagree", kind="warn"))

    warrant = (
        (f"{headline.value} for", "; ".join(diseases) or _EM_DASH)
        if headline is not None
        else ("admitted by", "; ".join(submitters) or _EM_DASH)
    )
    facts = [
        warrant,
        ("validity", fact.validity_state.value),
        ("atlas curation", fact.atlas_curation.value),
        ("lesion groups", ", ".join(group.value for group in fact.lesion_groups) or _EM_DASH),
        ("assertions", str(fact.assertion_count)),
        ("functional records", str(fact.functional_count)),
        # "curated publications", not "publications". `fact.publications` is the
        # PMIDs this atlas's own assertion evidence cites, and it reads 0 on 22
        # of the published genes -- which sat directly beside a burden section
        # citing a linked, PubMed-referenced study. The count was right and the
        # label was wrong. Raised by review 2026-08-05.
        ("curated publications", str(len(fact.publications))),
        ("burden rows", str(burden_rows)),
    ]
    bundle = gene_bundle_path(HgncId(gene))
    # `html.escape` here rather than a `render.py` primitive: these three are
    # bare interpolations into markup this function assembles, and every other
    # value on the page reaches it through `chip`, `definition_list` or
    # `data_table`, which escape their own arguments.
    return (
        f'<aside class="rail"><h1>{html.escape(symbol)}</h1>'
        f"<p><code>{html.escape(gene)}</code></p><p>{''.join(chips)}</p>"
        f"{definition_list(facts)}"
        f'<p><a href="../{html.escape(bundle)}">this gene as JSON</a></p></aside>'
    )


def _validity_section(gene_validity: GeneValidity) -> str:
    """The mirrored table, verbatim terms and all.

    `classification_term` is what the authority published, not the rung it maps
    onto: GenCC's `Supportive` maps to no rung at all, and a table showing the
    mapped value would render that submitter's row as blank or, worse, as a
    grade nobody gave.
    """
    rows = [
        Row(
            cells=(
                record.source.value,
                record.submitter or record.gcep or _EM_DASH,
                f"{record.disease_label} ({record.disease})",
                record.moi,
                record.classification_term,
                record.sop or _EM_DASH,
                record.classification_date or _EM_DASH,
                _report_link(record.report_url),
            )
        )
        for record in gene_validity.records
    ]
    return (
        "<h2>Mirrored gene&ndash;disease validity</h2>"
        "<p>Every classification below is an upstream panel's or submitter's, "
        "republished with its provenance intact.</p>"
        f"{data_table(_VALIDITY_HEADERS, rows)}"
    )


def _evidence_section(
    records: Sequence[LesionAssertion], publications: Mapping[str, Publication]
) -> str:
    """Each assertion, its evidence, and the publications that evidence cites.

    A cited publication missing from `publications` renders as an em dash rather
    than being dropped from the list: the citation is the curator's, and a row
    that disappears because a lookup missed is the silent evidence loss this
    project treats as its characteristic failure. `referential.py` already
    refuses a corpus that cites an unknown publication, so this is a guard on a
    bypassed gate.
    """
    blocks: list[str] = []
    for assertion in records:
        summary = definition_list(
            [
                ("lesion groups", ", ".join(group.value for group in assertion.lesion_groups)),
                ("inheritance", ", ".join(assertion.inheritance)),
                ("mechanism", assertion.mechanism),
                ("syndromic", assertion.syndromic),
            ]
        )
        items = [
            Row(
                cells=(
                    item.evidence_class.value,
                    item.strength.value,
                    item.summary,
                    Link(text=item.publication, href=_pubmed(item.publication)),
                )
            )
            for item in assertion.evidence
        ]
        blocks.append(
            f"<h3>{html.escape(assertion.id)}</h3>{summary}{data_table(_EVIDENCE_HEADERS, items)}"
        )
    # Sorted: the set is built by iterating assertions, and a `set` of strings
    # iterates in an order that varies with `PYTHONHASHSEED` -- which would make
    # the page's checksum differ between two builds of one commit.
    cited = sorted({item.publication for record in records for item in record.evidence})
    citations = []
    for reference in cited:
        publication = publications.get(reference)
        citations.append(
            Row(
                cells=(
                    reference,
                    publication.title if publication is not None else _EM_DASH,
                    str(publication.year) if publication is not None else _EM_DASH,
                )
            )
        )
    return (
        "<h2>Curated evidence</h2>"
        + "".join(blocks)
        + "<h2>Publications</h2>"
        + data_table(_CITATION_HEADERS, citations)
    )


def _study_label(pmid: str, publications: Mapping[str, Publication]) -> str:
    """ "Audain et al. 2026", or the bare PMID when the record is missing.

    The fallback is a guard on a bypassed gate -- BUR010 reports a burden row
    citing an unregistered study and `build_site` refuses on it -- and it renders
    the PMID rather than an em dash for the reason `build_genes` falls back to an
    HGNC id: an identifier is something a reader can still look up.
    """
    publication = publications.get(pmid)
    if publication is None:
        return pmid
    # `Publication.authors` constrains the *list* to `min_length=1`, not its
    # elements, so a whitespace-only first author passes every validator and
    # renders "    et al. 2026" -- a study heading naming nobody. Same defect
    # `runner.py::_cell` exists to catch on the mirror side, reaching a page for
    # the first time here: `authors[0]` was previously read only by `search.py`.
    author = publication.authors[0].strip()
    if not author:
        return pmid
    return f"{author} et al. {publication.year}"


def _count(carriers: int | None, total: int | None, unit: str) -> str:
    """ "6 carriers / 3,876" -- a numerator, its denominator, and their unit.

    Never the numerator alone. Six carriers is a different claim in 3,876 cases
    than in 45,082, and the two columns of this table hold exactly that contrast.
    The separator is a literal `,` via `:,`, which is locale-independent, so two
    builds on two machines render the same bytes.

    **And never without the unit.** This is `_effect`'s rule applied to the
    count columns, for the same reason and after the same near miss: until
    2026-08-05 the header read "cases (carriers / n)" and every cell was assumed
    to be people. PMID:40127276's case-control arm counts *alleles* -- measured,
    not assumed, in `CountUnit` -- so that header would have claimed 21,768
    people sequenced where 11,555 were, and called an allele count a carrier
    count, on a page whose checksum verified. An unknown unit therefore renders
    its raw token rather than nothing: there is no branch here that omits it.
    """
    if carriers is None or total is None:
        return _EM_DASH
    numerator, denominator = _COUNT_UNIT_LABEL.get(unit, (unit, ""))
    return f"{carriers:,} {numerator} / {total:,}{denominator}"


def _effect(row: BurdenRow) -> str:
    """The effect size, its measure, and its interval -- never a bare number.

    **This is the guard the single `effect` column was chosen against.** One
    column holding both odds ratios and de novo enrichments is what lets this
    schema absorb a fifth study without a migration, and it is also the one place
    two incomparable quantities could silently merge: an odds ratio of 3.1 and an
    enrichment of 3.1 are different claims, and a cell reading "3.1" under a
    header reading "effect" equates them. `_MEASURE_LABEL` is therefore consulted
    on every row and there is no branch that omits it.

    `unbounded` rather than a number, and rather than a blank. Fisher's exact
    test returns an infinite odds ratio where no control carries, which is the
    strongest result in the data and also the one `encode_json`'s
    `allow_nan=False` refuses to publish. The lower bound survives and is the
    whole finding: TAB2's syndromic row says the true odds ratio is at least
    28.1, which a blank cell would have thrown away.
    """
    if row.effect_measure is None:
        return _EM_DASH
    measure = _MEASURE_LABEL.get(row.effect_measure, row.effect_measure)
    if row.effect_bound == "unbounded_above":
        # `∞ (95% CI 0.298–∞)`, not `unbounded (95% CI ≥0.298)`.
        #
        # Measured 2026-08-05 on the built site: 19 cells render an unbounded
        # odds ratio and **11 of them have a lower bound below 1** -- GATA6 and
        # TBX5 at ≥0.298 with p = 0.079, among others. The old wording led with
        # the word "unbounded", which reads as "infinitely enriched", and put
        # the bound that contradicts it in a parenthesis at lower salience, in a
        # `≥` notation that is not how a Fisher interval is written. Eleven times
        # out of nineteen the parenthetical was the entire message.
        #
        # The conventional form puts the interval in the reader's usual shape,
        # so an interval spanning 1 is visible as one. `∞` is safe here and not
        # in the payload: this is HTML text, while `allow_nan=False` governs the
        # JSON, which carries `null` plus `effect_bound`.
        low = _fmt(row.ci_low) if row.ci_low is not None else None
        return f"{measure} ∞" + (f" (95% CI {low}–∞)" if low else "")
    if row.effect is None:
        return _EM_DASH
    value = f"{measure} {_fmt(row.effect)}"
    if row.ci_low is None or row.ci_high is None:
        return value
    return f"{value} (95% CI {_fmt(row.ci_low)}–{_fmt(row.ci_high)})"


def _fmt(value: float) -> str:
    """A statistic at three significant figures, without scientific notation.

    `f"{1810.0:.3g}"` is `1.81e+03`. Measured on the built site, that spelling
    reached four confidence-interval upper bounds -- `1.81e+03`, `2.35e+03`,
    `2.18e+03`, `1.36e+03` -- sitting mid-column beside plain decimals like
    `1.3`. The upper bound is the number that says "do not act on this", and in
    the least legible glyph in the table a skimming reader can take `1.81e+03`
    for `1.81`, turning an uninformative interval into a narrow-looking one.

    Values at or above 1,000 render as separated integers instead: `1,810`.
    Three significant figures are already far more precision than an interval
    that wide carries. `:,` is locale-independent, like `_count`'s.

    Not used for p-values: there an exponent is correct and expected, and
    `2.45e-07` is how a reader wants to see it.
    """
    if abs(value) >= 1000:
        return f"{value:,.0f}"
    return f"{value:.3g}"


def _corrected(row: BurdenRow) -> str:
    """A published corrected p-value, named for the correction that produced it.

    **This column is the answer to the sharpest finding of the 2026-08-05
    review**, which measured 187 uncorrected p-values on gene pages with nothing
    to judge them against. The atlas still computes no correction -- that would
    be authoring a statistic (D12/D33) -- but where a study publishes one, it
    reaches the reader.

    It changes conclusions, not just presentation. CHD7 in PMID:34324492 has a
    raw permutation p of 0.0068 and a family-wise corrected p of 0.99: without
    this column the page shows a number that reads as significant for a gene the
    study found nothing for.

    An em dash where a study published no correction, never a blank: the method
    line above the table then names how many comparisons the study ran, which is
    what a reader has instead.
    """
    if row.pvalue_adjusted is None:
        return _EM_DASH
    label = _ADJUSTMENT_LABEL.get(row.pvalue_adjustment or "", row.pvalue_adjustment or "")
    return f"{row.pvalue_adjusted:.3g} ({label})" if label else f"{row.pvalue_adjusted:.3g}"


def _burden_section(
    rows: Sequence[BurdenRow],
    publications: Mapping[str, Publication],
    cohorts: Mapping[str, Cohort],
    families: tuple[frozenset[str], ...] = (),
    axes: tuple[tuple[str, str], ...] = (),
) -> str:
    """Published rare-variant burden, one table per study.

    Three things this section has to say that the numbers alone do not, each of
    which a reader would otherwise get wrong:

    * **The rows are per study and are never pooled.** The CHD literature reuses
      cohorts, so a pooled p-value counts the same children twice. Declining to
      compute one does not stop a reader doing it by eye, so where two studies
      share a collection the section names it (`shared_cohorts`).
    * **A consequence class with no row is one the study did not report** -- not
      "was tested and found nothing", and *why* it is absent is the study's own
      rule rather than a fact about the data. This bullet claimed the stronger
      thing until 2026-08-05: "had no carrier in either group... measured over
      the committed mirror: zero of its 1,192 rows have no case carrier *and* no
      control carrier". That measurement was over PMID:42230622 alone, and both
      halves have since failed. The mirror now holds 1,475 rows and **one of
      them does have zero carriers on both sides**; and PMID:40127276 observed
      14,364 synonymous variants (its Dataset S4) while publishing no synonymous
      row at all, because its gene-level table covers only damaging classes. The
      matrix is genuinely sparse, so a reader meets a gap often enough for the
      distinction to matter -- which is why the sentence must be weak enough to
      be true of every study on the page.
    * **The synonymous row is the study's own negative control**, and only a
      study that publishes one has it. Synonymous variants should show no
      enrichment; where one does, that gene's comparison is poorly calibrated.
      It is sorted last within each stratum so it reads as what it is -- the row
      that says whether to believe the two above it. 69 of the 290 published
      rows are synonymous and every one is from PMID:42230622, so
      `_SYNONYMOUS_NOTICE` is conditional: a page with no such row must not
      promise a check it cannot offer.

    Returns `""` for a gene with no burden rows rather than an empty section:
    unlike the validity table, whose header names the columns and whose emptiness
    is itself an answer, an empty burden table would say "this gene was studied
    and nothing was found", which is a claim no study made.
    """
    if not rows:
        return ""

    overlaps = shared_cohorts(rows)
    blocks: list[str] = []
    for study in sorted({row.study for row in rows}):
        study_rows = [row for row in rows if row.study == study]
        table = data_table(
            _BURDEN_HEADERS,
            [
                Row(
                    cells=(
                        _STRATUM_LABEL.get(row.cohort_stratum, row.cohort_stratum),
                        _CONSEQUENCE_LABEL.get(row.consequence_class, row.consequence_class),
                        _count(row.n_case_carriers, row.n_cases, row.count_unit),
                        _controls(row),
                        _effect(row),
                        f"{row.pvalue:.3g}" if row.pvalue is not None else _EM_DASH,
                        _corrected(row),
                    )
                )
                for row in study_rows
            ],
        )
        design, warning = _method_line(study_rows, publications.get(study))
        # `design` and `_provenance` fold together: both are identical for a
        # given study on all 91 gene pages carrying one. `warning`, the disclosure, the
        # composite note and the footnotes do not -- each changes what a reader
        # concludes from the table right below it.
        counted = (
            f'<details class="reading-notes"><summary>How this study counted'
            f"</summary>{design}{_provenance(study_rows, cohorts)}</details>"
            if design or _provenance(study_rows, cohorts)
            else ""
        )
        # One panel per effect measure, never one per study: an odds ratio of
        # 3.1 and a de novo enrichment of 3.1 are different claims, and 75 of
        # the 915 published rows are the second kind. `sorted` over a set, for
        # the reason every other loop here is sorted.
        placed = [row for row in study_rows if _plottable(row)]
        measures = sorted({row.effect_measure for row in placed if row.effect_measure is not None})
        forests = "".join(
            _forest(
                [row for row in placed if row.effect_measure == measure],
                measure,
                _study_label(study, publications),
            )
            for measure in measures
        )
        # **Folded only where a panel above it shows the same rows.** THE FOLD
        # RULE says a caveat may fold when it is general and must stay visible
        # when it is particular; a table of this gene's own counts is as
        # particular as anything on the page, and what earns the fold is not
        # that it repeats but that the picture directly above it is drawn from
        # exactly these rows. Where no row is plottable there is no picture, so
        # the table is the only record and stays in plain sight -- 3 published
        # genes are in that state (KLF13, CFC1, MYH11), and folding away the
        # only figures they have would be evidence loss dressed as tidiness.
        figures = (
            f'<details class="study-figures"><summary>every figure from this study'
            f"</summary>{table}</details>"
            if forests
            else table
        )
        blocks.append(
            f"<h3>"
            f'<a href="{html.escape(_pubmed(study))}">'
            f"{html.escape(_study_label(study, publications))}</a></h3>"
            f"{_disclosure(study, publications)}"
            f"{warning}"
            f"{counted}"
            f"{_composite_note(study_rows)}"
            f"{forests or _NO_EFFECT_MEASURE}"
            f"{figures}"
            f"{_footnotes(study_rows)}"
        )

    shared = "".join(
        f'<p class="notice-inline"><strong>These two studies are not independent.</strong> '
        f"{html.escape(_study_label(left, publications))} and "
        f"{html.escape(_study_label(right, publications))} both draw on "
        f"{html.escape(_names(common, cohorts))}, so their results describe partly "
        f"the same people and must not be combined.</p>"
        for (left, right), common in overlaps.items()
    )

    # The pooling sentence is conditional. Rendered unconditionally it asserted
    # "these cohorts overlap" in the present tense on a page showing one table,
    # sending a reader to look for a second study that is not there and
    # undermining the caveats that *are* live.
    pooling = _POOLING_NOTICE if len({row.study for row in rows}) > 1 else ""

    # Conditional for the same reason as `pooling` and `_composite_note`: a page
    # whose tables carry no synonymous row must not tell a reader to weigh one.
    synonymous = (
        _SYNONYMOUS_NOTICE if any(row.consequence_class == "synonymous" for row in rows) else ""
    )

    matrix = _evidence_matrix(rows, families, axes, publications)

    # THE FOLD RULE. A caveat may fold if it is *general* -- true on every gene
    # page, teaching how to read the table. It must stay visible if it is
    # *particular* -- it fired because of this gene's data and changes what a
    # reader concludes about this gene.
    #
    # Folded here: the units glossary, what an absent cell means, and the
    # synonymous negative control. All three are identical on all 91 pages that
    # carry a burden table, and
    # none is lost -- the unit word is still in every cell and the matrix
    # renders `not tested` as its own state.
    #
    # Outside the fold, always: `_MIRRORED_NOTICE` (one line, and it is the
    # atlas's claim about itself), the `shared_cohorts` non-independence notice,
    # the own-lab disclosure, and every row's `method_note` footnote -- which is
    # where PMID:34324492's thoracic-aortic-aneurysm contamination reaches a
    # reader. `test_every_particular_caveat_survives_outside_the_fold` is what
    # stops one of those being tidied inside later.
    reading_notes = (
        '<details class="reading-notes"><summary>How to read these numbers'
        "</summary>" + _BURDEN_PREAMBLE + synonymous + "</details>"
    )

    return (
        "<h2>Rare variant burden</h2>"
        + _MIRRORED_NOTICE
        + matrix
        # RELOCATED, not folded. A summary that counts datasets is adjacent to
        # pooling and invites exactly the arithmetic this sentence forbids, so it
        # belongs directly beneath the matrix rather than four paragraphs above
        # it. Rendered whenever there is a matrix, because the matrix is what it
        # qualifies -- unlike its previous home, where it was conditional on
        # there being more than one study to pool.
        + (_POOLING_NOTICE if matrix else pooling)
        + shared
        + reading_notes
        + "".join(blocks)
        + _cohort_notes(rows, cohorts)
    )


def _controls(row: BurdenRow) -> str:
    """The control column, which a de novo row fills with an expectation.

    `_count(None, None)` renders an em dash, and on a `mutation_model` row that
    dash is where the *entire comparator* belongs: the modelled expected count
    is the only thing the enrichment was computed against, and an em dash is
    indistinguishable from a control count nobody recorded.
    """
    if row.comparator == "mutation_model" and row.expected_count is not None:
        return f"{_fmt(row.expected_count)} expected"
    return _count(row.n_control_carriers, row.n_controls, row.count_unit)


def _disclosure(study: str, publications: Mapping[str, Publication]) -> str:
    """Say so when the study is by an author of this atlas.

    `Publication.own_lab` exists to record exactly this and reached
    `publications.json` and zero HTML files. Two of the three studies published
    today -- `PMID:42230622` and `PMID:34324492` -- are by this repository's
    author, and a reader weighing those numbers should be told without having to
    cross-reference a JSON payload.

    **The conditional is now load-bearing rather than latent.** This docstring
    said "every burden table on the site today is from PMID:42230622" until
    2026-08-05, when `PMID:40127276` (`own_lab: false`) began rendering on 15 of
    the 23 gene pages published then. Re-measured 2026-08-06 against the widened
    population: 91 of the 92 pages carry a burden table, drawn from three studies
    -- PMID:42230622 on 87, PMID:34324492 on 61, PMID:40127276 on 25 -- and 13
    carry all three, with the own-lab declaration on **two** of the three
    (PMID:42230622 and PMID:34324492; PMID:40127276 is not this atlas's). That is
    the arrangement that makes the line informative: an unconditional disclosure
    says nothing, and a stale claim that it is unconditional invites someone to
    make it one. This sentence read "one of the three" until 2026-08-06, which
    contradicted this docstring's own opening paragraph two screens up.
    """
    publication = publications.get(study)
    if publication is None or not publication.own_lab:
        return ""
    return (
        '<p class="disclosure">Declaration: this study is by an author of this atlas. '
        "Its numbers are republished here exactly as published, and this atlas adds no "
        "assessment of them.</p>"
    )


def _labelled(values: Iterable[str | None], labels: Mapping[str, str]) -> list[str]:
    """Distinct values of one column, mapped to display text, sorted.

    Sorted because these come from a `set`, whose iteration order varies with
    `PYTHONHASHSEED` -- and this string is part of a page's bytes and therefore
    its checksum. Sorted on the label rather than the raw value so the rendered
    line reads in the order it is shown.
    """
    return sorted({labels.get(value, value) for value in values if value is not None})


def _method_line(rows: Sequence[BurdenRow], publication: Publication | None) -> tuple[str, str]:
    """What was counted, and against how many other tests.

    **Three of the eight partition columns reached no page.** `variant_class`,
    `origin` and `maf_max` were populated and appeared nowhere, so a section
    headed "Rare variant burden" never said what *rare* meant, and a
    case-control count of rare-inherited variants was indistinguishable from a
    trio's de novo count. `tables.py` states the case against exactly this: the
    partition "exists to stop two incomparable rows from *looking* comparable".

    This said "four" and named `lesion_group` as the fourth until 2026-08-05.
    **`lesion_group` is null on all 290 published rows and always has been** --
    no curated study stratifies by lesion group yet -- so a quarter of the
    stated rationale was for a column that had nothing to render. `_LESION_LABEL`
    is deliberately left empty rather than removed: the column is live in the
    schema and a study that uses it should render, at which point that map needs
    filling. Measured, not assumed; the claim it replaces was not.

    Rendered per study rather than per row, from the distinct values across the
    block. **A field that is multi-valued within a study renders every value**,
    which says what the block contains without claiming which row is which --
    the JSON is row-precise and the preamble says so. That is no longer
    hypothetical: PMID:40127276 contributes two origins and two tests to one
    block, and its `maf_max` is populated on the case-control rows and null on
    the de novo rows, which is why the frequency clause below distinguishes
    "filtered at this threshold" from "not filtered at all" rather than dropping
    the null.

    The multiple-testing sentence is the other half, and the more serious one.
    Measured 2026-08-05 over the built site: 32 of PMID:42230622's 187 published
    rows have p < 0.05 and **3** survive Bonferroni over that study's own
    138,609 tests, so 29 rows read as significant and are not. That study
    publishes no corrected p, and computing one would be authoring a statistic
    (D12/D33) -- but naming the denominator is a fact the study supplies, and it
    is what lets a reader apply their own threshold. The other two studies do
    publish corrections, so they get the other branch below.
    """
    parts: list[str] = []
    for values, labels in (
        ([row.variant_class for row in rows], _VARIANT_CLASS_LABEL),
        ([row.origin for row in rows], _ORIGIN_LABEL),
        ([row.lesion_group for row in rows], _LESION_LABEL),
        ([row.pvalue_test for row in rows], _TEST_LABEL),
    ):
        rendered = _labelled(values, labels)
        if rendered:
            parts.append(", ".join(rendered))
    # A null `maf_max` means "this analysis applied no frequency filter", which
    # is a fact about the rows and not an absence of one. Dropping it and then
    # stating the surviving threshold over the whole block attributed
    # PMID:40127276's MAF < 1e-05 to its 45 de novo rows, which had no filter at
    # all -- a de novo mutation is defined by absence from both parents, not by
    # a population frequency. Measured 2026-08-05: that sentence was wrong on 15
    # of the 23 published gene pages.
    mafs = sorted({row.maf_max for row in rows if row.maf_max is not None})
    unfiltered = any(row.maf_max is None for row in rows)
    if mafs and unfiltered:
        parts.insert(
            1,
            "MAF below "
            + ", ".join(f"{value:g}" for value in mafs)
            + " for the rows that applied a frequency filter, and no filter on the rest",
        )
    elif mafs:
        parts.insert(1, "MAF below " + ", ".join(f"{value:g}" for value in mafs))
    elif unfiltered:
        parts.insert(1, "no frequency filter")

    sentence = html.escape("; ".join(parts)) + "." if parts else ""

    # Three cases, because a sentence true of one study is false of another.
    # PMID:42230622 publishes no correction at all, so the reader needs the
    # denominator; PMID:34324492 corrects every row, so telling them to judge
    # the raw p against the scan would be pointing at a column they should not
    # be reading. A study that corrected only some rows gets both halves.
    correction = ""
    uncorrected = [row for row in rows if row.pvalue is not None and row.pvalue_adjusted is None]
    if publication is not None and publication.tests_reported:
        total = publication.tests_reported
        threshold = 0.05 / total
        correction = f" This study reported {total:,} such comparisons."
        if uncorrected:
            correction += (
                f" Where it published no corrected p, judge the raw one against the whole "
                f"scan rather than against 0.05 &mdash; a Bonferroni threshold over that "
                f"many is p &lt; {threshold:.2g}."
            )
        else:
            correction += (
                " The <strong>corrected p</strong> column carries the study's own "
                "correction over all of them; read it rather than the raw p."
            )
    elif uncorrected:
        correction = " p-values are <strong>uncorrected</strong> as published."

    if not sentence and not correction:
        return ("", "")

    # **The two halves are returned separately because they fold differently.**
    #
    # The design clause ("SNVs and indels; MAF below 0.001; Fisher's exact
    # test") is *general*: identical for a given study on all 23 gene pages, and
    # it teaches how to read the block rather than saying anything about this
    # gene. It folds, with the provenance line, under "How this study counted".
    #
    # The multiple-testing clause does not. It was the sharpest finding of the
    # 2026-08-05 review -- 32 of PMID:42230622's 187 published rows clear 0.05
    # and 3 survive Bonferroni over its own 138,609 comparisons -- so it changes
    # what a reader concludes from every p-value in the block below it. Folding
    # it would reintroduce, one release later, exactly the defect that review
    # found.
    #
    # Returned as a pair rather than as one string the caller splits: this
    # module builds HTML and must never parse it back, and a split on "</p>"
    # would break silently the first time either half gained a paragraph.
    design = f'<p class="method">{sentence}</p>' if sentence else ""
    warning = f'<p class="method">{correction.strip()}</p>' if correction else ""
    return (design, warning)


def _provenance(rows: Sequence[BurdenRow], cohorts: Mapping[str, Cohort]) -> str:
    """The collections a study's rows drew on, unioned across all of them.

    Read off `study_rows[0]` until 2026-08-05, which stated one row's cohorts
    over the whole block. `case_cohorts` and `control_cohorts` are per-row
    columns, so a study whose strata draw on different collections published a
    provenance line the data contradicted -- and `shared_cohorts` reads the same
    columns with a *union* rule, so the "not independent" notice could name a
    cohort the provenance line did not. Both now union, so they cannot disagree.
    """
    cases = sorted({identifier for row in rows for identifier in row.case_cohorts})
    controls = sorted({identifier for row in rows for identifier in row.control_cohorts})
    return (
        f'<p class="provenance">Cases: {html.escape(_names(cases, cohorts))}. '
        f"Controls: {html.escape(_names(controls, cohorts) or 'none (see the method above)')}.</p>"
    )


def _concordance_for(
    concordance: Mapping[str, Mapping[str, Json]] | None, gene: str
) -> Mapping[str, Json]:
    """One gene's concordance, or a refusal naming the gap.

    **`bundles._concordance_for` raised on a missing gene and this page defaulted
    it**, so a mapping built over the wrong population would have made the bundle
    fail loudly and the browse row render "0 of 0 tested" -- a tally byte-identical
    to a measured "no study reported this gene". The guard existed on one of the
    two layers `runner.py` says cannot disagree.

    `None` is still allowed and means "this build has no burden data at all",
    which is what every unit test of this page passes. A *mapping* that omits a
    gene is the error: it means the caller derived it over a different population.
    """
    if concordance is None:
        return {"tested": 0, "enriched": 0, "corrected": 0, "families": []}
    if gene not in concordance:
        raise KeyError(
            f"no concordance derived for published gene {gene}; the mapping was built over "
            f"a different population, and defaulting it would render 'no dataset tested "
            f"this' as if it had been measured"
        )
    return concordance[gene]


def _dot_strip(concordance: Mapping[str, Json]) -> str:
    """One glyph per cohort family, plus a tally naming both statistics.

    **Fill encodes the correction**, so a single glyph carries two facts a
    reader would otherwise need two columns for: that a dataset showed
    enrichment, and whether that survived the study's own correction. TAB2 is
    the case it was designed against -- three families all pointing the same
    way, of which only one survives correction.

    The tally counts *tested* families and says both numbers, because either
    alone misleads. `3 of 3` hides that two do not survive correction;
    `1 corrected` hides that PMID:42230622 publishes no correction on any row
    and so can never earn one.

    A family that did not test the gene gets its own glyph and is excluded from
    the denominator. That is not presentation: KDM6A is ClinGen definitive and
    shows nothing in either dataset that tested it, and "0 of 3" beside a green
    chip reads as the data contradicting the classification, where "0 of 2
    tested" reads as what it is.
    """
    families = concordance.get("families") or ()
    if not isinstance(families, Sequence):
        return _EM_DASH
    glyphs: list[str] = []
    for entry in families:
        if not isinstance(entry, Mapping):
            continue
        state = str(entry.get("state", ""))
        studies = entry.get("studies") or ()
        names = ", ".join(str(study) for study in studies) if isinstance(studies, Sequence) else ""
        title = f"{names}: {_STATE_TITLE.get(state, state)}"
        glyph = _STATE_GLYPH.get(state, "untested")
        body = "&ndash;" if glyph == "untested" else ""
        glyphs.append(f'<span class="dot {glyph}" title="{html.escape(title)}">{body}</span>')
    tested = concordance.get("tested", 0)
    corrected = concordance.get("corrected", 0)
    enriched = concordance.get("enriched", 0)
    tally = f"{enriched} of {tested} tested &middot; {corrected} corrected"
    # Not a link. The row already links to the gene page twice -- from the id and
    # from the symbol -- and a third link wrapped around a tally puts the
    # affordance on a number rather than on a name. Each dot still carries a
    # `title` naming its studies and state, which is the detail this cell owes a
    # reader; the numbers behind it are one click away by the name they scanned
    # for.
    return f'<span class="strip">{"".join(glyphs)}<span class="strip-tally">{tally}</span></span>'


def _evidence_matrix(
    rows: Sequence[BurdenRow],
    families: tuple[frozenset[str], ...],
    axes: tuple[tuple[str, str], ...],
    publications: Mapping[str, Publication],
) -> str:
    """Cohort family against evidence design, with the holes left visible.

    **The empty cells are the point.** Columns come from the whole corpus rather
    than from this gene, so a dataset that did not test this gene renders as an
    absence rather than vanishing -- the same confusion between "not tested" and
    "tested and found nothing" that `FamilyState` has four members for.

    **The guarantee is narrower than an earlier version of this docstring
    claimed.** It said the CNV de novo quadrant "is empty on every page". It is
    on *no* page: `evidence_axes` derives columns from the rows that exist, so a
    design *nobody* ran produces no column at all. What the matrix shows is a
    hole where one family did not run a design another family did; it is silent
    about a design the whole corpus lacks. Measured 2026-08-05: three columns,
    none of them CNV de novo.

    The cell shows the most significant row in it, ordered by corrected p where
    published and by raw p otherwise, tie-broken on
    `(consequence_class, cohort_stratum)` so the choice is deterministic across
    builds. Where a cell holds more than one row the count is stated, so a
    reader knows the cell summarises and the table below is the record.
    """
    if not rows or not families or not axes:
        return ""

    by_study: dict[str, list[BurdenRow]] = {}
    for row in rows:
        by_study.setdefault(row.study, []).append(row)

    header = "".join(
        f'<th scope="col">{html.escape(_VARIANT_CLASS_LABEL.get(variant, variant))}'
        f'<br><span class="sub">{html.escape(_COMPARATOR_LABEL.get(comparator, comparator))}'
        f"</span></th>"
        for variant, comparator in axes
    )

    body: list[str] = []
    for family in families:
        members = [row for study in sorted(family) for row in by_study.get(study, ())]
        label = " + ".join(_study_label(study, publications) for study in sorted(family))
        cells: list[str] = []
        for variant, comparator in axes:
            in_cell = [
                row
                for row in members
                if row.variant_class == variant and row.comparator == comparator
            ]
            cells.append(_matrix_cell(in_cell))
        body.append(f'<tr><th scope="row">{html.escape(label)}</th>{"".join(cells)}</tr>')

    return (
        '<div class="scroll"><table class="matrix">'
        f"<thead><tr><td></td>{header}</tr></thead>"
        f"<tbody>{''.join(body)}</tbody></table></div>" + _MATRIX_LEGEND
    )


def _matrix_cell(rows: Sequence[BurdenRow]) -> str:
    """One family's result for one evidence design.

    An empty cell says `not tested` in words rather than rendering blank: a
    blank cell is indistinguishable from a rendering bug, and this one carries a
    claim -- that this dataset did not run this design on this gene.
    """
    if not rows:
        return (
            f'<td><span class="cell not-tested" title="{_STATE_TITLE["not_tested"]}">'
            f"{_STATE_CELL['not_tested']}</span></td>"
        )

    state = family_state(rows)
    # Most significant first: corrected p where published, raw p otherwise, then
    # a stable tie-break so two builds pick the same row.
    best = sorted(
        rows,
        key=lambda row: (
            row.pvalue_adjusted if row.pvalue_adjusted is not None else (row.pvalue or 1.0),
            row.consequence_class,
            row.cohort_stratum,
        ),
    )[0]

    if state is FamilyState.NO_ENRICHMENT:
        detail = f"n = {len(rows)}" if len(rows) > 1 else "&nbsp;"
        return (
            f'<td><span class="cell no-enrichment" title="{_STATE_TITLE["no_enrichment"]}">'
            f'{_STATE_CELL["no_enrichment"]}<br><span class="sub">{detail}</span></span></td>'
        )

    # **The cell is a summary; the table below is the record.** So it carries the
    # effect and one statistic and stops there. The full spelling -- the
    # confidence interval, the name of the correction -- goes on the `title` and
    # is in the table in full a screen further down.
    #
    # This is a width decision with a correctness edge. The interval and the
    # correction name ran a cell to 61 characters, which set the column width for
    # the whole grid and crowded the cells against each other. What may *not* be
    # dropped is the measure: `_effect` has no branch that omits it, because an
    # odds ratio of 3.1 and a de novo enrichment of 3.1 are different claims, and
    # `_effect_compact` keeps that property.
    kind = "corrected" if state is FamilyState.CORRECTED else "nominal"
    # `_statistic` is shared with the forest panel above, so the two pictures on
    # one page cannot name the same number two different ways; the rule it
    # keeps -- `q` only for a real false-discovery rate -- is recorded there.
    statistic = _statistic(best)
    more = f'<span class="sub"> &middot; {len(rows)} rows</span>' if len(rows) > 1 else ""
    full = f"{_effect(best)}; {_corrected(best) if best.pvalue_adjusted is not None else ''}".strip(
        "; "
    )
    return (
        f'<td><span class="cell {kind}" title="{html.escape(full)}">'
        f"{html.escape(_effect_compact(best))}<br>"
        f'<span class="sub">{html.escape(statistic)}</span>{more}</span></td>'
    )


def _effect_compact(row: BurdenRow) -> str:
    """The effect and its measure, without the interval. For the matrix cell only.

    **The measure survives; only the interval is dropped.** `_effect`'s rule --
    that no branch omits the measure, because an odds ratio of 3.1 and a de novo
    enrichment of 3.1 are different claims -- is the one property this must not
    trade for width. The interval is a different kind of thing: it qualifies the
    estimate rather than naming it, and it is in the row's own table a screen
    below, plus on this cell's `title`.

    An unbounded effect still renders `OR &infin;` rather than a blank, for the
    reason `_effect` does: it is the strongest result in the data and the one
    `allow_nan=False` refuses to publish as a number.
    """
    if row.effect_measure is None:
        return _EM_DASH
    measure = _MEASURE_LABEL.get(row.effect_measure, row.effect_measure)
    if row.effect_bound == "unbounded_above":
        return f"{measure} ∞"
    if row.effect is None:
        return _EM_DASH
    return f"{measure} {_fmt(row.effect)}"


# --- Effect and uncertainty, drawn ------------------------------------------
#
# The evidence matrix a screen above answers *which datasets tested this gene
# and which found enrichment*. It cannot answer **how large** and **how
# certain**, because a swatch has no room for a number and no room for an
# interval. That is what this panel is for, and it is the half of the burden
# layer a reader was previously left to assemble from a nine-row table by eye.
#
# **Every (study, effect measure) pair gets a panel, including the two studies
# that publish no interval at all.** Measured 2026-08-20 over the 915 published
# rows: PMID:42230622 publishes an interval on 100% of its 704 rows and
# PMID:40127276 on 0% of its 150 -- and PMID:40127276 is the study with the
# findings that survive correction. On TBX5 it reports de novo loss-of-function
# enriched 297x at q 5.6e-08, while none of PMID:42230622's eight TBX5 rows
# carries a corrected p at all and seven of their intervals cross 1. The eighth
# does not, and this comment said all eight did until it was measured: the
# non-syndromic synonymous row sits entirely below 1, at OR 0.251 (95% CI
# 0.030-0.924, p 0.0375). It is that study's own negative control, so it is the
# one interval a sentence about this panel least ought to misdescribe.
#
# Drawing only where an interval exists would put a picture of the uncorrected
# result on the page and leave the surviving finding as text -- curated
# evidence visually demoted, which is this repository's characteristic failure
# in a new medium.

# One panel's geometry, in the SVG's own coordinates. The label column is wide
# because the longest label this vocabulary produces -- "non-syndromic ·
# damaging (LOF + missense) union" -- is 47 characters, and a label running
# under the plot area would overlap the very interval it names.
_FOREST_WIDTH: Final = 600
_FOREST_ROW: Final = 20.0
_FOREST_TOP: Final = 26.0
_FOREST_LABEL_X: Final = 2.0
_FOREST_PLOT_LEFT: Final = 244.0
_FOREST_PLOT_RIGHT: Final = 470.0
_FOREST_STAT_X: Final = 478.0
_FOREST_AXIS_GAP: Final = 10.0
_FOREST_FOOT: Final = 16.0

# Text is dropped by an explicit offset rather than by `dominant-baseline`,
# whose initial value differs between renderers: a chart whose labels sit half
# a row too high in one browser and correctly in another is a chart nobody can
# check against the table below it.
_FOREST_TEXT_DROP: Final = 3.4
_FOREST_ARROW: Final = 7.0
_FOREST_ARROW_HALF: Final = 3.5

# A quarter decade of clearance beyond the extreme values, so that no measured
# number lands on an edge of the axis. The edges are **reserved**: an arrow
# there means "past the end of this axis", and a real value sharing that pixel
# would read as one. Same reason `_FLOOR_GAP` exists on the trajectory chart --
# the atlas's weakest statement about a figure and its most extreme real one
# must not land on the same pixel column.
_EFFECT_AXIS_PAD: Final = 10.0**0.25

# What "survived the correction" means, and the only threshold this module
# applies to any number. It is applied to a p-value the **study** corrected and
# published; the atlas corrects nothing (D12/D33), which is why a row whose
# study published no correction at all is drawn hollow rather than filled --
# nothing says it survived one.
_CORRECTED_ALPHA: Final = 0.05

# The two unions this schema contains, each as (the union's own value, the
# components it is the union of). `_composite_note` states the first in words
# and `_CONSEQUENCE_LABEL["damaging"]` names it in the table; the second is
# stated by the counts themselves (measured on TBX5: 3,876 = 1,471 + 2,405).
# The plot has to mark both, for the reason the note exists: bars stacked as
# siblings read as independent findings, and here two of them are inside a
# third.
_CONSEQUENCE_UNION: Final[tuple[str, tuple[str, ...]]] = (
    "damaging",
    ("lof", "missense_damaging"),
)
_STRATUM_UNION: Final[tuple[str, tuple[str, ...]]] = (
    "all",
    ("syndromic", "nonsyndromic"),
)

# What a page says instead of drawing an empty frame (D42). Measured 2026-08-20
# on a real build: 3 published genes -- KLF13, CFC1, MYH11 -- carry only
# PMID:34324492's CNV rows, which publish no effect measure of any kind, and a
# fourth (CRIPTO) carries no burden row at all and so reaches this section
# never. An empty axis under a heading naming a study asserts that the study
# measured something and found it to be nothing, which is not what happened.
_NO_EFFECT_MEASURE: Final = (
    '<p class="method">This study publishes <strong>no effect measure</strong> for this '
    "gene &mdash; no odds ratio, enrichment or rate ratio &mdash; so there is nothing to "
    "place on an effect axis. Its counts and p-values are in the table below.</p>"
)


def _plottable(row: BurdenRow) -> bool:
    """Whether this row can be placed on an effect axis at all.

    An `unbounded_above` row **is** plottable and renders as an arrow: those
    rows carry the strongest signals in the mirror (TAB2: 5 syndromic
    carriers, 0 of 45,082 controls), and `ci_low` is the finding. Excluding
    them was an error in this plan's own first measurements.

    A row with an effect measure, no effect and no unbounded flag is refused
    rather than drawn as a bar with nothing on it. It does not occur in the
    mirror -- `validate/burden.py` pairs `effect` with `effect_measure` -- and
    a panel row with no estimate and no arrow would read as an estimate of
    zero, which is a different and much stronger claim.
    """
    return row.effect_measure is not None and (
        row.effect is not None or row.effect_bound == "unbounded_above"
    )


def _union_kinds(rows: Sequence[BurdenRow]) -> tuple[bool, bool]:
    """Whether each union is *shown as one* in this panel: (consequence, stratum).

    Conditional on the components being in the same panel, exactly as
    `_composite_note` is conditional on the study reporting both. Measured
    2026-08-20: every PMID:40127276 row is `cohort_stratum: all` and that study
    publishes no syndromic or non-syndromic row at all, so its `all` rows are
    the whole cohort rather than the sum of two rows a reader can see. Tagging
    them `union` would point at rows that are not on the page.

    **This answers for the panel, and a tag is worn by a row.** A panel can
    turn a union on for one consequence and then hang the tag on a row whose
    own components are not all drawn -- 10 panels do. `_fully_split` is what
    asks the per-row question, and `_forest_caption` is where the answer
    changes what the key promises.
    """
    composite, parts = _CONSEQUENCE_UNION
    consequence = any(row.consequence_class == composite for row in rows) and any(
        row.consequence_class in parts for row in rows
    )
    whole, strata = _STRATUM_UNION
    stratum = any(row.cohort_stratum == whole for row in rows) and any(
        row.cohort_stratum in strata for row in rows
    )
    return consequence, stratum


def _is_union(row: BurdenRow, kinds: tuple[bool, bool]) -> bool:
    """Whether this row contains other rows in its own panel."""
    consequence, stratum = kinds
    return (consequence and row.consequence_class == _CONSEQUENCE_UNION[0]) or (
        stratum and row.cohort_stratum == _STRATUM_UNION[0]
    )


def _fully_split(row: BurdenRow, rows: Sequence[BurdenRow], kinds: tuple[bool, bool]) -> bool:
    """Whether **every** row this union contains is itself drawn on this panel.

    `_union_kinds` asks its question of the panel; a tag is worn by a row, and
    the two are not the same question. A panel carrying `all cases · missense`
    beside `syndromic · missense` and `non-syndromic · missense` turns the
    stratum union on for the whole panel -- and then `all cases · loss-of-
    function`, whose non-syndromic half the study never published, wears the
    same tag with only one component under it.

    Measured 2026-08-20 on the built corpus: 288 rows carry the tag, **0 of
    them with neither component present**, so the tag is never false; but on
    10 panels a tagged row has one of its two components missing. All 10 are
    a (stratum, consequence) cell absent from the mirror entirely -- and a
    missing cell is not a null result, so nothing here may read as one.

    Both axes are checked with the *other* axis held fixed, because that is
    what the union is: `damaging` for `all cases` is the loss-of-function and
    damaging-missense rows **for `all cases`**, not for any stratum.
    """
    consequence, stratum = kinds
    if consequence and row.consequence_class == _CONSEQUENCE_UNION[0]:
        drawn = {
            other.consequence_class
            for other in rows
            if other.cohort_stratum == row.cohort_stratum
            and other.consequence_class in _CONSEQUENCE_UNION[1]
        }
        if drawn != set(_CONSEQUENCE_UNION[1]):
            return False
    if stratum and row.cohort_stratum == _STRATUM_UNION[0]:
        drawn = {
            other.cohort_stratum
            for other in rows
            if other.consequence_class == row.consequence_class
            and other.cohort_stratum in _STRATUM_UNION[1]
        }
        if drawn != set(_STRATUM_UNION[1]):
            return False
    return True


def _survived_correction(row: BurdenRow) -> bool:
    """Whether this row cleared **its own study's** correction for multiple testing.

    Never a correction this atlas computed -- there is none (D12/D33). A row
    whose study published no correction is `False` for the same honest
    reason a row that failed one is: nothing says it cleared anything.

    One predicate, because it decides the fill of two different glyphs.
    """
    return row.pvalue_adjusted is not None and row.pvalue_adjusted < _CORRECTED_ALPHA


def _mark_class(row: BurdenRow, glyph: str) -> str:
    """Filled where the study's own correction was survived, hollow where not.

    **`glyph` is the shape, and the fill rule is the same for every shape.**
    A panel draws two of them -- `point` for a placed estimate, `arrow` for a
    row the axis has no position for -- and each has the same three states,
    so the class names are composed here rather than written out per shape.
    `_forest_estimate` bypassed this function entirely for its two arrow
    branches and drew every arrow in `chart-arrow`, whose declaration was
    byte-identical to `chart-point`'s. Measured 2026-08-20 on the built
    corpus: **205 of the 854 plotted marks were arrows, not one of them had
    survived any correction, and every one was drawn in the fill that says it
    had.** 30 were corrected and failed, with q as high as 1.0, and 175 came
    from a study publishing no correction at all; **13 of the 30 sat on a
    ClinGen `definitive` gene**, GATA4 among them, whose Sierant enrichment
    panel drew three solid marks of which the third was `p 1 · q 1`. Two
    families of class name, written twice, is how that happened; one function
    is the fix.

    Measured 2026-08-20 over the 915 published rows: 55 carry an adjusted p
    below 0.05 and 45 of those are plottable, all from PMID:40127276.

    A synonymous row that *does* survive is drawn muted rather than in the
    result colour, because it is the study's own negative control and a
    surviving one is a warning about the comparison rather than a finding
    (`_SYNONYMOUS_NOTICE` says the rest). No published row is in that state
    today -- all 261 published synonymous rows come from the study that
    publishes no correction -- so the branch is exercised by its test rather
    than by the corpus. Hollowness is what may not be traded away here: it
    carries the correction, which the muted colour does not.
    """
    if not _survived_correction(row):
        return f"chart-{glyph}-open"
    return f"chart-control-{glyph}" if row.consequence_class == "synonymous" else f"chart-{glyph}"


def _statistic(row: BurdenRow) -> str:
    """The one statistic that speaks for a row: the published correction where
    there is one, and the raw p otherwise.

    `q` only where the correction really is a false-discovery rate. The mirror
    carries two methods -- `benjamini_hochberg`, which is an FDR, and
    `familywise_permutation`, which is not -- and labelling the second `q`
    contradicted the same page's own table, which names it "family-wise" a
    screen below. Measured 2026-08-05: 29 matrix cells rendered `q`, some of
    them over a family-wise p.

    One function for the matrix cell and the forest panel, so the two pictures
    on one page cannot name the same number two different ways.
    """
    if row.pvalue_adjusted is None:
        return f"p {row.pvalue:.3g}" if row.pvalue is not None else _EM_DASH
    if row.pvalue_adjustment == "benjamini_hochberg":
        return f"q {row.pvalue_adjusted:.3g}"
    return f"corrected p {row.pvalue_adjusted:.3g}"


def _forest_statistic(row: BurdenRow) -> str:
    """Both statistics for a plotted row: the raw p, and the correction if any.

    The matrix cell shows one because it has one line; a panel row has a whole
    column, and the pair is what the fill of its marker is derived from. A
    reader who sees a hollow marker beside `p 0.045 · q 0.29` can check the
    encoding against the numbers, which is the only way a legend is worth
    anything.
    """
    raw = f"p {row.pvalue:.3g}" if row.pvalue is not None else _EM_DASH
    if row.pvalue_adjusted is None:
        return raw
    return f"{raw} · {_statistic(row)}"


def _effect_bounds(rows: Sequence[BurdenRow]) -> tuple[float, float]:
    """The axis one panel's rows are placed on, padded clear of both edges.

    **1.0 is always among the values**, so the null line is always on the
    axis: a reader must always be able to see which side of "no enrichment" a
    bar falls on, and a panel whose rows are all above 1 would otherwise draw
    no null at all and read as if every row were enriched relative to the
    weakest of them.

    Only positive values are collected. A zero effect and a zero lower bound
    have no position on a logarithmic axis -- `LogScale` refuses them rather
    than clamping, because a clamped zero reads as "the smallest number here"
    -- and `_forest_estimate` draws them as an arrow off the left edge
    instead. Measured 2026-08-20: 182 of the 915 published rows have an effect
    of exactly zero, every one of them with no case carrier, and 152 of those
    also publish a zero lower bound.

    `_axis_bounds` is what makes a panel whose values are all equal drawable at
    all; the padding is applied after it, so the flat case is widened and then
    padded rather than either alone.
    """
    values = [1.0]
    for row in rows:
        values.extend(
            value
            for value in (row.effect, row.ci_low, row.ci_high)
            if value is not None and value > 0
        )
    low, high = _axis_bounds(values)
    return low / _EFFECT_AXIS_PAD, high * _EFFECT_AXIS_PAD


def _forest_line(x1: float, x2: float, y: float, css_class: str) -> str:
    """One horizontal rule at `y`, from `x1` to `x2`."""
    return (
        f'<line class="{html.escape(css_class)}" x1="{coordinate(x1)}" y1="{coordinate(y)}" '
        f'x2="{coordinate(x2)}" y2="{coordinate(y)}"/>'
    )


def _forest_arrow(x: float, y: float, *, pointing: float, css_class: str) -> str:
    """A triangle at an axis edge: this row continues past the end of the axis.

    `pointing` is +1 for the right-hand edge and -1 for the left. The apex sits
    exactly on the edge, which is why `_EFFECT_AXIS_PAD` keeps every measured
    value away from it.

    `css_class` is required rather than defaulted, so that adding a third
    arrow branch cannot quietly inherit the fill of the first two -- which is
    exactly how every arrow on the site came to be drawn as a surviving
    result. `_mark_class` is what callers pass.
    """
    base = x - pointing * _FOREST_ARROW
    return (
        f'<polygon class="{html.escape(css_class)}" points="{coordinate(x)},{coordinate(y)} '
        f"{coordinate(base)},{coordinate(y - _FOREST_ARROW_HALF)} "
        f'{coordinate(base)},{coordinate(y + _FOREST_ARROW_HALF)}"/>'
    )


# SVG's own closed vocabulary for `text-anchor`. Named rather than
# interpolated freely because this is the one *attribute value* the chart code
# builds from a parameter, and `render.data_table` already records why that is
# the dangerous shape: `html.escape` rewrites the quote characters and leaves
# the space and the `=` alone, so a value reaching an attribute outside its
# quotes opens a second attribute. That primitive *refuses* an attribute name
# rather than escaping it, for exactly this reason; this is the same rule one
# level down.
_TEXT_ANCHORS: Final = frozenset({"start", "middle", "end"})


def _forest_text(x: float, y: float, content: Markup, *, anchor: str) -> str:
    """One label.

    `content` is `Markup` rather than `str`, so the obligation is visible at
    the call site instead of living in this sentence: a row label carries a
    `<tspan>` for its union tag, and escaping here would publish the tag as
    literal text. `render.py` expresses the identical contract the identical
    way -- taking a bare `str` here was an asymmetry a review named, not a
    considered difference.

    `anchor` is checked against SVG's own closed set. Every call site passes a
    literal, so this is a guard on a bypassed gate in `data_table`'s idiom:
    reaching it means a caller derived an attribute value from data, and that
    must fail rather than publish. `raise`, never `assert` -- `-O` strips
    `assert`.
    """
    if anchor not in _TEXT_ANCHORS:
        raise ValueError(
            f"text-anchor {anchor!r} is not one of {sorted(_TEXT_ANCHORS)}; escaping "
            f"cannot make an attribute value safe, because a space in one opens a "
            f"second attribute"
        )
    return (
        f'<text class="chart-label" x="{coordinate(x)}" '
        f'y="{coordinate(y + _FOREST_TEXT_DROP)}" text-anchor="{anchor}">{content.html}</text>'
    )


def _forest_interval(row: BurdenRow, y: float, scale: LogScale) -> str:
    """The row's published interval, or the statement that there is none.

    A row with neither bound gets a hatched bar spanning the whole axis, not a
    short one and not nothing. Measured 2026-08-20 over the 915 published rows:
    211 publish neither bound, and every one of the 55 rows carrying an
    adjusted p below 0.05 is among them. A short bar would invent a precision
    the study never claimed, and no bar at all would make the strongest rows on
    the site the faintest marks on it.

    Only 150 of those 211 reach a panel at all -- the other 61 are
    PMID:34324492's, which publish no effect measure either, so `_plottable`
    refuses them before this is called. The two figures are not the same
    quantity and must not be written as one: 10 of the 55 surviving rows are in
    the 61, so "the plottable subset holds every surviving finding" is false
    while "the rows with no interval do" is true.

    A bound that is absent or zero resolves to the axis edge, where
    `_forest_estimate` has already drawn the arrow that says the row runs past
    it. Measured 2026-08-20 on the published corpus, the two coincide exactly:
    every row with a zero lower bound has an effect of zero, and every row with
    no upper bound is flagged `unbounded_above`.
    """
    if row.ci_low is None and row.ci_high is None:
        return _forest_line(_FOREST_PLOT_LEFT, _FOREST_PLOT_RIGHT, y, "chart-nointerval")
    low = scale.x(row.ci_low) if row.ci_low is not None and row.ci_low > 0 else _FOREST_PLOT_LEFT
    high = (
        scale.x(row.ci_high) if row.ci_high is not None and row.ci_high > 0 else _FOREST_PLOT_RIGHT
    )
    control = row.consequence_class == "synonymous"
    return _forest_line(low, high, y, "chart-control" if control else "chart-line")


def _forest_estimate(row: BurdenRow, y: float, scale: LogScale) -> str:
    """The point estimate, or an arrow where the axis has no position for it.

    **Every branch takes its fill from `_mark_class`**, arrows included. The
    arrow branches carried a hardcoded `chart-arrow` until 2026-08-20 and it
    resolved to the same paint as a survived-its-correction circle; see
    `_mark_class` for what that published.

    Two rows have no point to draw and both are drawn as arrows rather than as
    blanks:

    * `effect_bound: unbounded_above` -- Fisher returns an infinite odds ratio
      where no control carries, and `allow_nan=False` refuses to publish it.
      All 23 such published rows have zero control carriers. The lower bound is
      the finding, and it is the left end of the bar this arrow terminates.
    * an effect of exactly zero -- no case carried, on all 182 published rows
      in that state. Zero has no position on a logarithmic axis, and clamping
      it onto the axis floor would read as "the smallest effect measured here"
      rather than as "none".
    """
    if row.effect is None:
        return _forest_arrow(
            _FOREST_PLOT_RIGHT, y, pointing=1.0, css_class=_mark_class(row, "arrow")
        )
    if row.effect <= 0:
        return _forest_arrow(
            _FOREST_PLOT_LEFT, y, pointing=-1.0, css_class=_mark_class(row, "arrow")
        )
    return marker(scale.x(row.effect), y, css_class=_mark_class(row, "point"))


def _forest_row(row: BurdenRow, index: int, scale: LogScale, *, union: bool) -> str:
    """One burden row: its name, its interval, its estimate and its statistics."""
    y = _FOREST_TOP + index * _FOREST_ROW
    label = (
        f"{_STRATUM_LABEL.get(row.cohort_stratum, row.cohort_stratum)} · "
        f"{_CONSEQUENCE_LABEL.get(row.consequence_class, row.consequence_class)}"
    )
    # A `<tspan>` rather than a second `<text>` at a computed x: the label's
    # width depends on the reader's font, which this build cannot measure, and
    # a tag positioned by arithmetic would land on top of the longest labels.
    tag = '<tspan class="chart-union"> union</tspan>' if union else ""
    return (
        _forest_text(_FOREST_LABEL_X, y, Markup(html.escape(label) + tag), anchor="start")
        + _forest_interval(row, y, scale)
        + _forest_estimate(row, y, scale)
        + _forest_text(
            _FOREST_STAT_X, y, Markup(html.escape(_forest_statistic(row))), anchor="start"
        )
    )


def _forest_caption(
    rows: Sequence[BurdenRow], measure: str, bounds: tuple[float, float], kinds: tuple[bool, bool]
) -> str:
    """The key. Every clause names a mark that is actually on this panel.

    Conditional clause by clause, for the reason `_composite_note` and
    `_SYNONYMOUS_NOTICE` are conditional and `_POOLING_NOTICE` had to be
    reworded when it moved (CLAUDE.md section 4.27): a sentence's truth
    conditions travel with its position, and a legend entry for a glyph the
    reader cannot see sends them hunting for it.

    **The hollow clause was the one exception, in the function whose docstring
    says there are none.** Measured 2026-08-20 on the built corpus: 11 of the
    137 panels printed it with no hollow marker drawn, across 10 genes --
    CHD7, GATA4, KMT2D, MYH7, NKX2-5, NODAL, RBFOX2, RNF40, WDR5, ZIC3, six
    of them ClinGen `definitive` -- and on 9 of those the panel did draw
    arrows, so the key told the reader that everything visible had failed a
    correction. Every other clause measured 0 captioned-but-not-drawn.

    **And there was no clause for a filled mark at all**, so its meaning had
    to be inferred as the complement of the hollow one. That inference was
    wrong for all 205 arrows on the site until `_mark_class` began deciding
    an arrow's fill; it is stated now rather than left to be worked out.
    The filled clause excludes a surviving *synonymous* row, which is drawn
    muted rather than in the result colour and has its own clause below.

    **The axis clause was wrong in a different way: it named a mark that is
    on the panel and then said the wrong thing about it.** The measure label
    was already dynamic; the gloss beside it was the literal "no enrichment"
    for every measure, so 112 of the 137 panels told the reader that 1 on an
    odds-ratio axis means no enrichment. Both halves now resolve from the same
    token -- see `_NULL_LINE_GLOSS` for the measurement and the third
    vocabulary. The sibling `<title>` on the same figure had always varied
    correctly, which is why reading the panel did not show it.

    **The union clause is conditional twice**: once on a union being shown as
    one at all, and once more on whether every tagged row is fully split.
    Measured 2026-08-20: 10 panels tag a row with only one of its two
    components drawn, TBX5's Audain odds-ratio panel among them, and the key
    told the reader to reconcile it against a row that is not there.

    The pooling clause is **not** conditional, and it is worded as policy
    rather than as an observation about this page. It is here at all because
    the forest idiom ends in a summary diamond and a reader who knows the
    idiom will look for one; its absence has to be an answer rather than a
    gap. It deliberately does not reuse `_POOLING_NOTICE`'s own sentence:
    that constant is the matrix's caption, and a test asserts it survives
    outside the fold, which a second copy of its wording here would satisfy on
    its behalf.
    """
    low, high = bounds
    label = _MEASURE_LABEL.get(measure, measure)
    # Both halves resolve from the same token. The gloss was hardcoded beside
    # this already-dynamic label, which is how 112 odds-ratio panels came to
    # call 1 "no enrichment" -- see `_NULL_LINE_GLOSS`.
    gloss = _NULL_LINE_GLOSS.get(measure, _NULL_LINE_DEFAULT)
    parts = [
        f"Axis: <strong>{html.escape(label)}</strong> on a logarithmic scale, "
        f"{_fmt(low)} to {_fmt(high)}; the vertical line is <strong>1</strong>, "
        f"{html.escape(gloss)}. ",
    ]
    if any(_survived_correction(row) and row.consequence_class != "synonymous" for row in rows):
        parts.append(
            "A <strong>filled</strong> marker is a row that <strong>survived</strong> its own "
            "study's correction for multiple testing. "
        )
    if any(not _survived_correction(row) for row in rows):
        parts.append(
            "A <strong>hollow</strong> marker is a row that <strong>did not survive</strong> its "
            "own study's correction for multiple testing, or whose study published no correction "
            "at all &mdash; the atlas computes none of its own. "
        )
    parts.append(
        "A forest plot usually ends in a summary diamond. There is <strong>no pooled</strong> "
        "estimate here and none is computed anywhere on this site: two studies that share a "
        "sample collection describe partly the same children, and adding their rows would "
        "count them more than once. "
    )
    if any(row.ci_low is None and row.ci_high is None for row in rows):
        parts.append(
            "A <strong>hatched</strong> bar means <strong>no interval published</strong> for "
            "that row, not a wide one. "
        )
    if any(row.effect is None for row in rows):
        parts.append(
            "An arrow to the <strong>right</strong> is an effect unbounded above &mdash; no "
            "control carried the variant &mdash; so the study's lower bound is the whole "
            "finding and there is no point estimate to place. "
        )
    if any(row.effect is not None and row.effect <= 0 for row in rows):
        parts.append(
            "An arrow to the <strong>left</strong> is an effect of zero &mdash; no case "
            "carried the variant &mdash; which a logarithmic axis has no position for. "
        )
    consequence, stratum = kinds
    reasons = []
    if consequence:
        reasons.append(
            "<code>damaging</code> is the loss-of-function and damaging-missense rows together"
        )
    if stratum:
        reasons.append("<code>all cases</code> is the syndromic and non-syndromic rows together")
    if reasons:
        parts.append(
            "A row tagged <strong>union</strong> contains at least one other row in this panel "
            "and is not an independent finding: " + "; ".join(reasons) + ". "
        )
        if any(_is_union(row, kinds) and not _fully_split(row, rows, kinds) for row in rows):
            parts.append(
                "Here a tagged row is <strong>not fully split</strong>: at least one of the "
                "rows it contains is not itself plotted, so what sits under a union on this "
                "panel is part of it and not all of it. The table below carries every row "
                "this study published for this gene. "
            )
    if any(row.consequence_class == "synonymous" for row in rows):
        # Not the words "negative control", deliberately. That teaching is
        # `_SYNONYMOUS_NOTICE`'s and it is behind the fold, where
        # `test_the_general_reading_notes_are_folded_out_of_the_way` asserts it
        # stays; a caption repeating the phrase on all 94 pages would satisfy
        # that test on the notice's behalf and unfold nothing. What this clause
        # owes the reader is what the *muting* means, which nothing else says.
        parts.append(
            "The <strong>muted</strong> row is this study's <strong>synonymous</strong> row, "
            "drawn apart from the results it calibrates; <em>How to read these numbers</em>, "
            "above, says why."
        )
    return f"<figcaption>{''.join(parts).strip()}</figcaption>"


def _forest(rows: Sequence[BurdenRow], measure: str, study_label: str) -> str:
    """One study's rows for one effect measure, as effect and interval.

    **One panel per (study, effect measure), never one per study.** 779 of the
    915 published rows carry an odds ratio and 75 an enrichment ratio, and one
    study publishes both: an odds ratio of 3.1 and a de novo enrichment of 3.1
    are different claims, so placing them on one axis would equate them at the
    only scale a reader actually reads. `_effect` has no branch that omits the
    measure's name for the same reason; this is that rule in geometry.

    Rows are in `_sort_key` order, so the panel reads top to bottom in the
    order the table below it reads -- the composite above its components, the
    synonymous negative control last.

    Raises on an empty panel rather than emitting an axis with nothing on it.
    The caller filters on `_plottable` and never calls this with nothing, so
    this is a guard on a bypassed gate: a bypassed gate must fail rather than
    publish (`raise`, never `assert` -- `-O` strips `assert`).
    """
    if not rows:
        raise ValueError("a forest panel needs at least one row to place")
    ordered = sorted(rows, key=_sort_key)
    bounds = _effect_bounds(ordered)
    scale = LogScale(
        low=bounds[0],
        high=bounds[1],
        left=_FOREST_PLOT_LEFT,
        width=_FOREST_PLOT_RIGHT - _FOREST_PLOT_LEFT,
    )
    kinds = _union_kinds(ordered)
    axis_y = _FOREST_TOP + (len(ordered) - 1) * _FOREST_ROW + _FOREST_AXIS_GAP
    top = _FOREST_TOP - _FOREST_ROW / 2
    null_x = scale.x(1.0)
    body = (
        _forest_line(_FOREST_PLOT_LEFT, _FOREST_PLOT_RIGHT, axis_y, "chart-axis")
        + f'<line class="chart-axis" x1="{coordinate(null_x)}" y1="{coordinate(top)}" '
        f'x2="{coordinate(null_x)}" y2="{coordinate(axis_y)}"/>'
        + _forest_text(null_x, top - 6.0, Markup("1"), anchor="middle")
        + _forest_text(_FOREST_PLOT_LEFT, axis_y + 8.0, Markup(_fmt(bounds[0])), anchor="start")
        + _forest_text(_FOREST_PLOT_RIGHT, axis_y + 8.0, Markup(_fmt(bounds[1])), anchor="end")
        + "".join(
            _forest_row(row, index, scale, union=_is_union(row, kinds))
            for index, row in enumerate(ordered)
        )
    )
    noun = "statistic" if len(ordered) == 1 else "statistics"
    title = (
        f"{study_label}: {len(ordered)} published "
        f"{_MEASURE_LABEL.get(measure, measure)} {noun} for this gene on a logarithmic axis, "
        f"with 95% confidence intervals where this study published them, and no pooled summary."
    )
    figure = svg_figure(
        width=_FOREST_WIDTH,
        height=round(axis_y + _FOREST_FOOT),
        title=title,
        body=body,
    )
    return (
        f'<figure class="forest" data-effect-measure="{html.escape(measure)}">'
        f"{figure}{_forest_caption(ordered, measure, bounds, kinds)}</figure>"
    )


def _composite_note(rows: Sequence[BurdenRow]) -> str:
    """Say that the `damaging` rows are the union of the two below them.

    Without it a reader meets three consequence rows and reads three findings,
    then adds up carrier counts that already include each other. Measured on
    CHD7 in PMID:40127276: the damaging de novo row is 20 mutations, and the
    loss-of-function and damaging-missense rows below it are 16 and 4 of *those
    same* 20.

    The composite is not droppable in favour of its parts -- it carries its own
    p-value, which is not a function of theirs, and it is the analysis the study
    defines its results by -- so the relationship has to be stated rather than
    designed away.

    Conditional on the study actually reporting both a composite and a
    component. Rendered unconditionally it would assert a decomposition that
    does not exist for the two studies whose tables carry no `damaging` row at
    all, which is the defect `_POOLING_NOTICE` was made conditional for.
    """
    consequences = {row.consequence_class for row in rows}
    if "damaging" not in consequences or not consequences & {"lof", "missense_damaging"}:
        return ""
    return (
        '<p class="notice-inline">The <strong>damaging (LOF + missense)</strong> rows are '
        "the <strong>union</strong> of the loss-of-function and damaging-missense rows "
        "below them, not a third independent result &mdash; the same variants are counted "
        "in both. They are shown because this study defines its findings by the composite, "
        "and because its p-value is not derivable from the other two.</p>"
    )


def _footnotes(rows: Sequence[BurdenRow]) -> str:
    """`method_note`, rendered -- which `tables.py` said it was and it was not.

    The column's comment reads "rendered verbatim as a row footnote", and no
    code path read it. It exists for the CNV case, where "carrier" means
    different things in different papers (any overlap, exonic, whole gene), so
    the first curator to need it would have got a green build, a correct bundle
    and no page text.
    """
    notes = sorted({row.method_note.strip() for row in rows if row.method_note})
    if not notes:
        return ""
    items = "".join(f"<li>{html.escape(note)}</li>" for note in notes)
    return f'<ul class="footnotes">{items}</ul>'


def _cohort_notes(rows: Sequence[BurdenRow], cohorts: Mapping[str, Cohort]) -> str:
    """What each collection is, and the caveats that qualify every number here.

    **These reached no published byte.** `curation/cohorts.yaml` records, for
    UK Biobank, that its participants are adults recruited at 40-69 while a
    substantial proportion of the CHD cases were enrolled in childhood -- a
    survivorship bias that inflates every odds ratio on every page carrying a
    UK Biobank comparison -- and that
    the combined cohort is ~92% European. `models/cohort.py` says in as many
    words that `description` is where such caveats belong. `cohort_labels`
    returned `{id: name}` and dropped them, so the curator wrote the caveat, the
    model documented it, and no reader could reach it.

    In a `<details>` element: the caveats are long, they repeat on every page
    carrying a burden table, and
    a reader who has read them once should not have to scroll past them again.
    `<details>` is native HTML, needs no script, and its contents are in the
    document for a crawler, a `curl` and a reader with JavaScript disabled.
    """
    cited = sorted(
        {identifier for row in rows for identifier in (*row.case_cohorts, *row.control_cohorts)}
    )
    items = [
        f"<li><strong>{html.escape(cohorts[identifier].name)}</strong> &mdash; "
        f"{html.escape(cohorts[identifier].description)}</li>"
        for identifier in cited
        if identifier in cohorts
    ]
    if not items:
        return ""
    return (
        '<details class="cohort-notes"><summary>About these cohorts, and what '
        "qualifies every number above</summary><ul>" + "".join(items) + "</ul></details>"
    )


def _names(ids: Sequence[str], cohorts: Mapping[str, Cohort]) -> str:
    """Cohort ids rendered as the collections they name, comma-separated.

    An id with no record renders as the id, which a reader can still look up;
    BUR009 reports the absence and the gate refuses, so this is a guard on a
    bypassed gate.
    """
    record = {identifier: cohorts.get(identifier) for identifier in ids}
    return ", ".join(
        value.name if value is not None else identifier for identifier, value in record.items()
    )


# --- Developmental expression -----------------------------------------------
#
# `ExpressionProfile` (`build/profiles.py`) already carries every figure this
# section renders -- a percentile band, tau, a phase placement -- with every
# qualifier D39 requires attached to it in the JSON. Nothing below computes
# anything; it only chooses the English for a payload whose own module
# docstring says exactly that is a later task's job.
#
# D43 governs this section: **never a panel redundant with the source's own
# browser**. That is D32's actual test -- "if a panel would look the same on
# the Heart Cell Atlas, it does not belong here" -- and it replaces the proxy
# this comment carried until 2026-08-20, which read "sentences and small
# per-organ tables, never a plot".
#
# The proxy was wrong in both directions. It forbade the phase-banded
# trajectory, which no external browser has, while permitting the rendering
# that hid a live defect for the whole life of schema 2.11: stages published
# in alphabetical order, 4 wpc eighth, because 21 individually-correct blocks
# cannot show a trajectory and therefore cannot show a scrambled one. (This
# read "for three releases" until 2026-08-21; the profiles mirror landed on
# 2026-08-19 and nothing published a stage array before it, so 2.11 is the
# only version that ever carried the defect. See `manifest.py`'s 2.12 note.)
#
# A bare heart curve IS what the source's own browser shows and remains
# forbidden. The same curve banded by this atlas's own cardiac-phase
# vocabulary is not. `test_every_chart_carries_an_atlas_specific_axis` is what
# watches for the difference.

# What each `ProfileGap` means, in English. Split into three maps rather than
# one, because `ProfileGap`'s own docstring says only the first two members
# mean exactly the same fact wherever they appear -- `dataset_not_registered`
# and `detection_floor_undeclared` are dataset-level gaps with one meaning in
# both fields that publish them, while `no_quantile_grid` is placement-only and
# `one_organ_sampled`/`undefined` are tau-only. Sharing only the two that
# genuinely share a meaning is `EVIDENCE_STATE_LABELS`'s discipline applied
# here: one vocabulary where the fact is one fact, never a false shared
# spelling for two different ones.
_DATASET_GAP_CLAUSE: Final[dict[str, str]] = {
    ProfileGap.DATASET_NOT_REGISTERED.value: (
        "this dataset is not registered in the atlas's own records"
    ),
    ProfileGap.FLOOR_UNDECLARED.value: "this dataset has not declared a detection floor",
}

_PLACEMENT_GAP_CLAUSE: Final[dict[str, str]] = {
    **_DATASET_GAP_CLAUSE,
    ProfileGap.NO_QUANTILE_GRID.value: (
        "no complete percentile grid is published for this organ at this stage"
    ),
    ProfileGap.NOT_ON_A_LOG_AXIS.value: (
        "this atlas placed the measurement and a logarithmic axis has no "
        "position for it: the median is zero or below"
    ),
}

_SPECIFICITY_GAP_CLAUSE: Final[dict[str, str]] = {
    **_DATASET_GAP_CLAUSE,
    ProfileGap.ONE_ORGAN_SAMPLED.value: (
        "only one organ was sampled at this stage, and τ needs at least two"
    ),
    ProfileGap.UNDEFINED.value: "it could not be computed for this gene at this stage",
}

_EXPRESSION_HEADERS: Final = ("tissue", "median abundance", "percentile")

# The shape `bundles._expression_profile` publishes for a gene no curated bulk
# RNA-seq dataset mentions at all (`EMPTY_EXPRESSION_PROFILE`, imported rather
# than reconstructed here, so the two cannot drift on what "empty" means).
#
# **Chosen over "this dataset does not cover it", and that is a different
# claim.** `datasets == []` means no `profiles.tsv` row anywhere names this
# gene: there is no specific dataset in view to say did not cover it, and
# saying so would invent a check this atlas never ran. This is the `_not_
# curated` idiom, not `_burden_section`'s: an absent section here would read as
# "the atlas looked and found nothing", indistinguishable from what this
# sentence actually means, so the section renders regardless.
_NOT_CURATED_EXPRESSION: Final = (
    "<p>The atlas has <strong>not yet curated</strong> a developmental expression "
    "profile for this gene: no curated bulk RNA-seq dataset mentions it today. That is "
    "a different claim from a specific dataset having looked at this gene and found "
    "nothing.</p>"
)

# Unconditional -- rendered beside a page of nothing but detected, well-placed
# figures exactly as beside a gene with every organ below the floor. The
# precedent for both the wording and the placement decision is
# `_POOLING_NOTICE`: made conditional once, it read false in the position that
# motivated the fix, because a caveat that appears only next to the bad news
# reads as an excuse for it. This one is a general fact about bulk RNA-seq and
# CHD gene biology, true whatever this gene's own numbers say, so it carries no
# clause that could be false in either position.
#
# This is the burden layer's "a missing cell is not a null result" discipline,
# needed more here: that layer leaves a gap as an absence, while this one turns
# a low number into a sentence, which is the more persuasive -- and more
# dangerous -- thing to get wrong.
_BULK_DILUTION_NOTICE: Final = (
    '<p class="notice-inline"><strong>Every measurement below is from whole, bulk '
    "tissue, never a single cell type.</strong> A whole embryonic heart is "
    "cardiomyocytes alongside endocardium, epicardium, cardiac neural crest, cushion "
    "mesenchyme and blood, and many congenital heart disease genes act in one rare "
    "lineage &mdash; the cardiac neural crest and outflow tract, the second heart "
    "field, the conduction system. A gene expressed intensely in one such lineage can "
    "still read low, or fall below the detection floor, once diluted across the whole "
    "organ: a low or absent bulk figure is not evidence that a gene is unimportant to "
    "heart development.</p>"
)

_EXPRESSION_READING_NOTES: Final = (
    "<p>Each organ's median abundance is placed against that organ's own percentile "
    "grid for that dataset and developmental stage, and the number of genes measured "
    "travels with every percentile shown below. No percentile is published below a "
    "dataset's own detection floor, because the rank of an unreliable measurement is "
    "not itself a measurement.</p>"
)

# Conditional on the section actually showing more than one placement.
# Rendered on a page with at most one, this would warn a reader against
# comparing figures the page does not even lay out side by side -- the same
# reason `_composite_note` and `_SYNONYMOUS_NOTICE` are conditional on a row
# that could earn them actually being on the page.
_PERCENTILE_COMPARABILITY_NOTICE: Final = (
    "<p>Percentiles on this page are <strong>not comparable across organs or across "
    "developmental stages</strong>: each organ transcribes a different fraction of the "
    "whole gene set, and the reference distribution itself changes shape as the heart "
    "matures. A high percentile shown for one organ or stage says nothing about a "
    "lower one shown elsewhere on this page.</p>"
)


def _phase_sentence(phase: PhaseInfo) -> str:
    """Where a stage falls in the curated cardiac-phase vocabulary, or why not.

    **Every branch names a state; none renders blank.** A stage the curated
    vocabulary does not cover resolves to `OUTSIDE_WINDOW` with a reason
    naming where it misses -- this function says exactly that rather than
    omitting the line or rendering an empty one, which is the difference
    between a documented gap and what reads as a rendering bug.

    **`MATCHED` may name several phases, never just the first.** Human
    cardiac morphogenesis runs several processes concurrently -- at 6 wpc an
    AVSD gene needs "atrial septation" and a TGA gene needs "outflow tract
    septation", and both are true at once (see `models/phases.py`'s module
    docstring for the measured concurrency table) -- so `profiles.PhaseInfo.
    phase_ids` is a tuple, already sorted deterministically by
    `CardiacPhaseFile.phases_for`, and every id in it is rendered, joined by
    commas, with the noun pluralised for two or more. `phase_ids` is non-empty
    exactly when `reason` is `None` (`profiles.PhaseAssignment` guarantees the
    two are never both set), so an unreachable empty-tuple-under-MATCHED case
    still renders a named state rather than raising. The underscore-to-space
    rewrite is the only transformation applied to each id, because a phase id
    is this atlas's own slug (`curation/cardiac_phases.yaml`), not
    third-party mirrored text.

    Every other branch renders `reason` verbatim, through `html.escape` like
    every other value this function assembles by hand. `profiles.py` already
    wrote those as short, reader-facing clauses -- "outside the curated
    window", "post-natal", "stage not declared by this dataset", "no
    developmental stage recorded for this measurement" -- and a second,
    parallel vocabulary here could only ever repeat them or drift from them.
    """
    if phase["outcome"] == PhaseOutcome.MATCHED.value:
        ids = phase["phase_ids"]
        if not ids:
            return "Developmental phase: <strong>unnamed phase</strong>."
        labels = [html.escape(phase_id.replace("_", " ")) or "unnamed phase" for phase_id in ids]
        noun = "phase" if len(labels) == 1 else "phases"
        joined = ", ".join(f"<strong>{label}</strong>" for label in labels)
        return f"Developmental {noun}: {joined}."
    reason = phase["reason"] or "not available"
    return f"Developmental phase: {html.escape(reason)}."


def _specificity_sentence(spec: Specificity, cardiac_tissues: frozenset[str]) -> str:
    """tau, with its scale and its organ list, and the argmax-gated gloss.

    **No branch omits the scale or the tissue list.** They sit in the one
    f-string that also carries the number, so there is no code path that
    renders tau without them -- the same guarantee `_effect` keeps for its
    measure label, after the same near miss: a bare number under a header that
    does not say what it measures.

    **It said "Expression is heart-preferential at this stage" until
    2026-08-21, and that was an adjective this function's own last paragraph
    already forbade.** The phrase was gated correctly -- on `highest_in` being
    one of this dataset's declared `cardiac_tissues`, never on tau alone,
    because tau measures concentration and not location. What it was not
    gated on was *how much* preference the number supports. Measured over the
    committed corpus: 445 stage blocks carried it, median tau 0.708, but **65
    below 0.30, 30 below 0.20 and 16 below 0.15** -- at the extreme, TAB2 at
    4 wpc with tau 0.046, where heart's 124 leads ovary's 116 and hindbrain's
    115 across seven organs. A 7% lead is an argmax; "preferential" reads as
    a finding.

    So the argmax gate stays and the adjective goes. Three mutually exclusive
    endings -- the argmax is cardiac, the argmax is some other organ, or the
    peak is tied (`highest_in is None`) -- each naming the organ and leaving
    the reader to weigh it against the tau printed in the same sentence.

    No adjective and no band on the number itself: this function states the
    scale, the organs and the argmax and stops there, never "highly specific"
    or "broadly expressed" -- D39(c) reserves that judgement from the atlas.
    """
    tissues = ", ".join(html.escape(tissue) for tissue in spec["tissues"])
    base = (
        f"τ = {_fmt(spec['tau'])} ({html.escape(spec['scale'])} scale) across "
        f"{tissues} ({spec['n_tissues']} organs sampled)."
    )
    highest = spec["highest_in"]
    if highest is None:
        peak = "No single organ has the highest median at this stage; the peak is tied."
    elif highest in cardiac_tissues:
        peak = (
            f"Expression peaks in {html.escape(highest)} at this stage, one of "
            "this dataset's cardiac tissues."
        )
    else:
        peak = (
            f"Expression peaks in {html.escape(highest)} at this stage, which this "
            "dataset does not treat as a cardiac tissue."
        )
    return f'<p class="method">{base} {peak}</p>'


def _specificity_gap_sentence(reason: str | None, floor: float | None) -> str:
    """Why tau is absent for a stage, in the same voice as `_specificity_sentence`.

    `peak_below_detection_floor` is kept out of `_SPECIFICITY_GAP_CLAUSE`
    because it is the one reason with a number to show: the floor value, named
    for the same reason rule 5 requires it beside a tissue's own placement gap
    -- a reader told a figure is missing because of a floor is owed the floor.
    """
    if reason == ProfileGap.PEAK_BELOW_DETECTION_FLOOR.value:
        qualifier = f" (detection floor {_fmt(floor)})" if floor is not None else ""
        return (
            "τ is not available: every organ sampled at this stage is below the "
            f"detection floor{qualifier}."
        )
    clause = _SPECIFICITY_GAP_CLAUSE.get(reason or "", "no reason was recorded for this gap")
    return f"τ is not available: {clause}."


def _abundance(tissue: TissueProfileEntry) -> str:
    """ "100 rpkm (n=3 samples)" -- never the number alone.

    `n_samples` travels with every figure on this page, the `_count`/
    `count_unit` discipline applied here: the schema permits `n_samples` = 1,
    so "median" can be one observation, and a reader must be told rather than
    left to assume a study-sized sample sits behind every row.
    """
    noun = "sample" if tissue["n_samples"] == 1 else "samples"
    unit = tissue["unit"]
    return f"{_fmt(tissue['median_abundance'])} {unit} (n={tissue['n_samples']} {noun})"


def _percentile_cell(tissue: TissueProfileEntry, floor: float | None) -> str:
    """The percentile band, or why there is none -- never the words "not detected".

    **`below_detection_floor` is the one reason with its own required wording
    (rule 5).** "Not detected" is a positive-sounding negative assertion a
    clinical reader takes as evidence against a gene, and this atlas measured
    no such thing: a value below a bulk assay's detection floor is a
    measurement this atlas does not vouch for, not an absence. The floor's own
    value is shown beside it, because a reader cannot judge how far below
    without it.
    """
    placed = tissue["placement"]
    if placed is not None:
        band = ""
        if placed["q25_percentile"] is not None and placed["q75_percentile"] is not None:
            band = f", IQR {placed['q25_percentile']}-{placed['q75_percentile']}"
        return f"{placed['median_percentile']} of {placed['n_genes']:,} genes{band}"
    reason = tissue["not_placed_reason"]
    if reason == ProfileGap.BELOW_DETECTION_FLOOR.value:
        qualifier = (
            f" (detection floor {_fmt(floor)} {tissue['unit']})" if floor is not None else ""
        )
        return f"below the detection floor in whole {tissue['tissue']} at this stage{qualifier}"
    fallback = "no percentile is available for this measurement"
    return _PLACEMENT_GAP_CLAUSE.get(reason or "", fallback)


def _tissue_table(tissues: Sequence[TissueProfileEntry], floor: float | None) -> str:
    """One row per organ measured at one stage. A list, never a cross-organ grid.

    Rows, not columns: laying organs out side by side would itself be the
    layout rule 8 warns against, since a table invites a reader's eye to
    compare adjacent cells. `_PERCENTILE_COMPARABILITY_NOTICE` is what a list
    alone cannot say, and is rendered beside it where it applies.
    """
    rows = [
        Row(cells=(tissue["tissue"], _abundance(tissue), _percentile_cell(tissue, floor)))
        for tissue in tissues
    ]
    return data_table(_EXPRESSION_HEADERS, rows)


def _stage_block(
    stage: StageProfileEntry, cardiac_tissues: frozenset[str], floor: float | None
) -> str:
    """One (dataset, stage) cross-section: its phase, its tau, its organs."""
    label = html.escape(_stage_label(stage))
    spec = stage["specificity"]
    if spec is not None:
        specificity_html = _specificity_sentence(spec, cardiac_tissues)
    else:
        gap = _specificity_gap_sentence(stage["specificity_unavailable_reason"], floor)
        specificity_html = f'<p class="method">{gap}</p>'
    return (
        f"<h4>{label}</h4>"
        f'<p class="method">{_phase_sentence(stage["phase"])}</p>'
        f"{specificity_html}"
        f"{_tissue_table(stage['tissues'], floor)}"
    )


# --- The trajectory chart ----------------------------------------------------

# Below this many placed stages a trajectory is not a trajectory. Two points
# make a line segment, and a line segment is a claim about the interval
# between them that two measurements do not support. Measured over the
# committed corpus 2026-08-20: 82 of 92 published genes clear this in heart,
# 3 fall to markers, and 7 have nothing placed at all.
_STAGES_FOR_A_TRAJECTORY: Final = 3

_CHART_WIDTH: Final = 560
_CHART_HEIGHT: Final = 150
_PLOT_LEFT: Final = 46.0
_PLOT_RIGHT: Final = 550.0
_PLOT_TOP: Final = 12.0
_PLOT_BOTTOM: Final = 124.0

# The lowest placed measurement stops this far above the axis, and the gap is
# load-bearing rather than padding. `_trajectory` ticks a below-floor stage ON
# the axis; without the gap, the atlas's weakest statement about a figure ("we
# do not vouch for this one") and its lowest real one ("this is the smallest
# number we measured") would land on the same pixel row, which is the exact
# conflation `_percentile_cell` refuses to make in words.
_FLOOR_GAP: Final = 8.0

# The value axis's own two labels, in the margin `_PLOT_LEFT` already
# reserves. Without them the 104 px between `_PLOT_TOP` and the foot of the
# scale carried no number at all, and every gene's curve filled all of it:
# measured 2026-08-20, PKD1L1 spans 1 to 2 tpm and TBX20 spans 1 to 526, and
# the two were drawn as the same full-height excursion. Right-anchored, so a
# four-digit median (MYH6 reaches 4,761 tpm in heart) grows leftward into the
# margin rather than across the plot. The drop centres the glyphs on the
# point they label, the way `_FOREST_TEXT_DROP` centres a forest row's.
_AXIS_LABEL_X: Final = _PLOT_LEFT - 4.0
_AXIS_LABEL_DROP: Final = 3.4

# A series whose placed values are all equal has no range, and `LogScale`
# refuses a zero-span axis rather than silently placing every point on one
# pixel. One decade centred on the value is the substitute: half a decade
# either side puts a flat series in the middle of the plot, where it reads as
# "one value, no range", instead of pinned to the floor or the ceiling, where
# it would read as low or high against a scale nothing else occupies.
#
# **The condition is equal values, not few of them**, and the two are easy to
# write as one because they overlap. Measured 2026-08-20 on the committed
# corpus, three genes reach this branch on their heart trajectory: FOXH1 (one
# placed stage, 2 tpm), TFAP2B (one placed stage, 1 tpm) and MMP21 -- which is
# placed at *three* stages, all of them 1.0 tpm. A count of genes with few
# placed stages would be a different set: 9 of the 92 are placed at fewer than
# two heart stages, and 7 of those 9 place nothing at all, so `_trajectory`
# returns before there is an axis to scale. No gene reaches this branch
# through `_small_multiples`, whose axis pools every drawn organ.
_FLAT_SERIES_SPAN: Final = 10.0**0.5


def _even_x(index: int, count: int, left: float, right: float) -> float:
    """The `index`-th of `count` evenly spaced positions between `left` and `right`.

    One implementation for both pictures on this page. The trajectory and the
    organ small multiples plot the *same* stage sequence at different widths,
    and a reader reads down from one to the other; two copies of this
    arithmetic could put stage 5 four-ninths of the way along in one and
    five-ninths in the other, with nothing failing.
    """
    if count < 2:
        return (left + right) / 2
    return left + (right - left) * index / (count - 1)


def _stage_x(index: int, count: int) -> float:
    """Where the `index`-th of `count` stages sits on the horizontal axis.

    The axis is the curated stage **order**, never `wpc`, and the two are
    different claims. Eight of this dataset's 21 stages are post-natal and
    carry no `wpc` at all, so a time-proportional axis could not place them;
    and even among the prenatal ones the sampling is uneven (13, then 16, 18,
    19 wpc). Evenly spaced positions say only "these came in this order",
    which is the whole of what `Stage.order` supports. `_trajectory`'s caption
    says so in words, because a reader cannot see it from the picture.
    """
    return _even_x(index, count, _PLOT_LEFT, _PLOT_RIGHT)


def _half_step(count: int) -> float:
    """Half the gap between two adjacent stages, in pixels.

    A phase band runs half a step outside the first and last stage it covers,
    so that a phase covering exactly one stage draws a visible band rather
    than a zero-width rectangle nobody sees. One published phase does exactly
    that: `outflow_tract_septum_morphogenesis` covers 6 wpc alone.
    """
    if count < 2:
        return (_PLOT_RIGHT - _PLOT_LEFT) / 2
    return (_PLOT_RIGHT - _PLOT_LEFT) / (count - 1) / 2


def _axis_bounds(values: Sequence[float]) -> tuple[float, float]:
    """A log axis wide enough to place every value, and never zero-wide."""
    low, high = min(values), max(values)
    if high > low:
        return low, high
    return low / _FLAT_SERIES_SPAN, high * _FLAT_SERIES_SPAN


def _axis_label(value: float, y: float) -> str:
    """One end of the value axis, named.

    `_fmt`, never `charts.coordinate`. `coordinate` is the *pixel* formatter
    and fixes to a tenth of a pixel, so it would publish a median of 275 tpm
    as `275.0` and one of 4,761 as `4761.0`, beside `_fmt`'s `4,761` in the
    table below the chart -- the same figure spelled two ways on one page.
    """
    return (
        f'<text class="chart-label" x="{coordinate(_AXIS_LABEL_X)}" '
        f'y="{coordinate(y + _AXIS_LABEL_DROP)}" text-anchor="end">'
        f"{html.escape(_fmt(value))}</text>"
    )


def _stage_label(stage: StageProfileEntry) -> str:
    """A stage's own token, or the state a null one is in.

    The same fallback `_stage_block` renders as its heading, taken from one
    place so a chart's axis and the heading under it cannot spell a null stage
    two different ways.
    """
    return stage["stage"] if stage["stage"] is not None else "no stage recorded"


class _Median(NamedTuple):
    """One organ's figure at one stage, and -- if it is not plottable -- why not.

    A `NamedTuple` rather than a bare tuple because the fourth field is the
    one two captions were reading off the wrong axis. `value is None` says a
    stage was sampled and not placed; it does not say the gene reads low
    there, and `_band_caption` and `_spark_caption` both said so. The reason
    now travels *with* the point, so a caption cannot look it up on a
    different pass over the data than the one that decided the point was
    unplaceable -- the discipline `_adjacent_runs` keeps by taking its pixel
    points beside its stage indices.

    `gap` is `None` exactly when `value` is not, and is the raw `ProfileGap`
    token the payload recorded rather than a sentence: what the wording
    should be is `_PLACEMENT_GAP_CLAUSE`'s job, in one place, and this module
    already has three maps saying so.

    `position` rather than `index`, which is `tuple.index` and cannot be a
    field name on a `NamedTuple` at all -- mypy refuses it outright, which is
    the only reason this is not a silently shadowed method.
    """

    position: int
    label: str
    value: float | None
    gap: str | None


def _tissue_medians(entry: DatasetProfileEntry, tissue: str) -> list[_Median]:
    """One organ's series: `(stage index, stage label, plottable median, gap)`.

    `None` is a stage this dataset **sampled and did not place** -- below the
    detection floor, or with no complete percentile grid behind it. Those
    stages stay in the list rather than being filtered out, because a gap in
    the sequence and a stage nobody sampled must not render the same way:
    `_trajectory` ticks the first on the axis and the second is simply not
    here at all.

    A stage where this organ has no row is not in the list, for that reason.
    Measured 2026-08-20 on the committed corpus: heart has a row at 19 of the
    21 stages this dataset publishes, missing at `school age child` and
    `elderly`, for all 92 published genes.

    **What is plotted is exactly what the atlas placed.** A below-floor
    figure has a median too, and plotting it would put a number this atlas
    declines to vouch for on the same line as ones it does -- the discipline
    `_percentile_cell` already keeps in words. A non-positive median folds
    into the same `None` for a second, narrower reason: it cannot be located
    on a log axis at all. Neither loses a figure; every one of them is in the
    table inside the `<details>`.
    """
    series: list[_Median] = []
    for index, stage in enumerate(entry["stages"]):
        for measured in stage["tissues"]:
            if measured["tissue"] != tissue:
                continue
            plottable = measured["median_abundance"]
            placed = measured["placement"] is not None and plottable > 0
            series.append(
                _Median(
                    position=index,
                    label=_stage_label(stage),
                    value=plottable if placed else None,
                    # `or` rather than a bare read: a cell the payload placed
                    # whose median is `<= 0` carries no `not_placed_reason` by
                    # `TissueProfileEntry`'s own contract, and would otherwise
                    # reach `_Median` with neither a value nor a gap -- the one
                    # pairing its docstring rules out.
                    gap=(
                        None
                        if placed
                        else (measured["not_placed_reason"] or ProfileGap.NOT_ON_A_LOG_AXIS.value)
                    ),
                )
            )
    return series


def _tissue_unit(entry: DatasetProfileEntry, tissue: str) -> str:
    """The unit this organ's medians are quoted in, or `""` if it has none.

    Read off the measurements rather than off the `Dataset`, for the same
    reason `_abundance` reads it there: the unit is a property of the mirrored
    row, and a page must not quote a floor in one unit beside a median in
    another.
    """
    for stage in entry["stages"]:
        for measured in stage["tissues"]:
            if measured["tissue"] == tissue:
                return measured["unit"]
    return ""


def _banded_phases(
    entry: DatasetProfileEntry, phases: CardiacPhaseFile | None
) -> list[tuple[CardiacPhase, int, int]]:
    """Each curated phase this dataset's stages fall inside, and the span it covers.

    **A phase whose `end_basis` is `NOT_STATED` is skipped, unconditionally.**
    It has no curated width, so a band drawn for it would have to run to an
    edge nobody stated -- asserting exactly the boundary that field exists to
    refuse, and asserting it in the one encoding a reader takes at a glance
    rather than reads. `heart_looping` is the live case: the source states
    when it begins and never when it ends.

    That duplicates a rule `CardiacPhaseFile.phases_for` already keeps, and
    the duplication is deliberate. `phases_for` is why no stage the real
    pipeline builds ever names an unended phase; this function is why that
    stays true of the *picture* even if some later caller assembles
    `phase_ids` another way. A guard added to one layer is not a guard.

    Sorted by phase id -- and **not as a determinism guard**, which is what
    this sentence claimed until it was measured. `covered` is a dict, filled
    by walking a list of stages and, inside each, `phase_ids` in the
    `(start_wpc, id)` order `phases_for` fixes; no `set` is iterated anywhere,
    so its own key order is already reproducible. Measured 2026-08-20 over six
    explicit `PYTHONHASHSEED` values (0, 1, 7, 42, 12345, 99991): one distinct
    build with the sort and one without.

    What the sort does buy is that the band order is a function of the covered
    ids alone rather than of the stage walk that found them, so it does not
    move when a dataset's stage list does. The two orders are genuinely
    different: dropping the sort changed 85 gene pages and the manifest, and
    reversing this dataset's 21 stages reorders the unsorted one while leaving
    the sorted one alone.
    """
    if phases is None:
        return []
    stated = {phase.id: phase for phase in phases.phases if phase.end_wpc is not None}
    covered: dict[str, list[int]] = {}
    for index, stage in enumerate(entry["stages"]):
        for phase_id in stage["phase"]["phase_ids"]:
            if phase_id in stated:
                covered.setdefault(phase_id, []).append(index)
    return [
        (stated[phase_id], min(covered[phase_id]), max(covered[phase_id]))
        for phase_id in sorted(covered)
    ]


def _phase_bands(entry: DatasetProfileEntry, phases: CardiacPhaseFile | None, count: int) -> str:
    """The shaded stage ranges this atlas's own phase vocabulary names.

    **This is what makes the chart admissible under D43.** The curve itself is
    what the source's own browser already draws; the bands are a claim no
    external browser makes, transcribed from PMID:32048790 into
    `curation/cardiac_phases.yaml` and carried onto every gene page here.
    `_trajectory` returns `""` when this returns `""`, so the rule is enforced
    by construction rather than by a comment asking a later author to keep it.

    Each band carries its own `<title>`, so a phase is nameable where several
    overlap -- and they do overlap, deliberately: three cardiac processes run
    concurrently at 6 wpc. Overlapping bands share one fill and merge into a
    single shaded region marking the morphogenetic window; the caption below
    the figure is what names each of them in text, which is also what a reader
    who cannot hover gets.
    """
    edge = _half_step(count)
    rects = []
    for phase, first, last in _banded_phases(entry, phases):
        left = max(_stage_x(first, count) - edge, 0.0)
        right = min(_stage_x(last, count) + edge, float(_CHART_WIDTH))
        rects.append(
            f'<rect class="chart-band" x="{coordinate(left)}" y="{coordinate(_PLOT_TOP)}" '
            f'width="{coordinate(right - left)}" '
            f'height="{coordinate(_PLOT_BOTTOM - _PLOT_TOP)}">'
            f"<title>{html.escape(phase.label)}</title></rect>"
        )
    return "".join(rects)


def _other_gap_clause(reasons: frozenset[str]) -> str:
    """Every gap reason *except* the detection floor's, spelled out.

    The floor is excluded because it is the one reason with its own required
    wording -- `_percentile_cell`'s rule 5, and the same sentence
    `_band_caption` and `_spark_caption` each carry in full. Everything else
    resolves through `_PLACEMENT_GAP_CLAUSE`, so a caption and the table cell
    three lines below it cannot spell one gap two ways.

    Sorted because `reasons` is a `frozenset` and an unsorted join over one
    would reorder a published sentence between two builds of one commit --
    **and, measured 2026-08-20, that cannot happen today.** `_tissue_entry`
    assigns `dataset_not_registered` and `detection_floor_undeclared` from
    `_dataset_gap`, which reads the `Dataset` alone, so either one applies to
    every cell of that dataset or to none; only `no_quantile_grid` is
    per-cell, and it is reached only where `_dataset_gap` returned nothing.
    Two distinct non-floor reasons therefore cannot co-occur inside one
    `DatasetProfileEntry`, this set holds at most one member, and removing
    the sort changes no published byte. The sort stays as the defensive kind
    section 9 keeps deliberately: this module renders whatever the payload
    says and cannot see an invariant that lives in `profiles.py`. It is not
    guarded by a test, and that is the reason -- `_banded_phases`' own
    docstring records the same correction, made after its sort was claimed as
    a determinism guard and measured not to be one.

    Returns `""` when the floor is the only reason present, which is every
    published page today -- measured 2026-08-20, all 2,817 unplaced cells on
    the committed corpus are below-floor.
    """
    others = sorted(reason for reason in reasons if reason != ProfileGap.BELOW_DETECTION_FLOOR)
    if not others:
        return ""
    return "; ".join(
        _PLACEMENT_GAP_CLAUSE.get(reason, "no reason was recorded for this gap")
        for reason in others
    )


def _band_caption(
    entry: DatasetProfileEntry,
    phases: CardiacPhaseFile | None,
    tissue: str,
    floor: float | None,
    low: float,
    high: float,
    series: Sequence[_Median],
    *,
    lined: bool,
) -> str:
    """What the picture cannot say about itself, in words beside it.

    Four things a reader would otherwise have to guess, and two of them are
    traps. The horizontal axis is stage *order*, not elapsed time, so the gap
    between 13 and 16 wpc is drawn the same width as the gap between 4 and 5.
    Saying so is the difference between an ordinal axis and a false linear
    one.

    **And the vertical axis is fitted to this gene alone**, which is the
    second trap and the one that was unstated. Every gene's curve fills the
    same 104 px whatever it spans: measured 2026-08-20, PKD1L1 runs 1 to 2
    tpm in heart and TBX20 runs 1 to 526, and the two were drawn as the same
    full-height excursion under a caption that named neither number. The
    range is printed here for the reason `_spark_caption` prints its own --
    a shared axis and a per-gene one are different claims, and a reader can
    only tell which they are looking at if the figure says.

    **Both break clauses are conditional on the break being in this figure.**
    The tick clause was not, and it named a mark most pages do not draw:
    measured 2026-08-20, all 85 trajectory pages said "the stage is ticked on
    the axis instead" while **55 drew zero ticks**, and 55 of the 59 with a
    visibly broken line were among them -- TBX5 has two polylines, a visible
    gap, and the only sentence about gaps described a mark that is not there.

    **And the second cause was never named at all.** A line also breaks where
    this dataset has no row for that organ at that stage, which for
    E-MTAB-6814 is the same set of stages for every gene (heart has none at
    `school age child` or `elderly`, forebrain none at 6 wpc, testis none at
    `neonate` or `school age child`). `_tissue_medians`' docstring insists
    that state and a below-floor one "must not render the same way", and in
    the picture they did.

    **It is stated as the axis's range, never as the medians'**, and the two
    differ on exactly the genes where the wording matters most. `_axis_bounds`
    widens a flat series by half a decade either side, because `LogScale`
    refuses a zero-span axis; 3 of the 85 charted genes are flat (TFAP2B and
    MMP21 at 1 tpm, FOXH1 at 2), and "median abundance in whole heart, 0.316
    to 3.16 tpm" would attribute to their medians a spread neither end of
    which anyone measured.
    """
    named = ", ".join(html.escape(phase.label) for phase, _, _ in _banded_phases(entry, phases))
    attribution = html.escape(phases.attributed_to) if phases is not None else ""
    unit = html.escape(_tissue_unit(entry, tissue))
    # Every stage this dataset sampled and did not place is ticked, whatever
    # the recorded reason -- so the clause explaining the mark is keyed on the
    # reasons that actually occur, never on `value is None`. Read off the
    # whole series rather than off `_line_breaks`, which sees only the gaps
    # that interrupt a *line*: a leading or trailing unplaced stage draws a
    # tick and breaks nothing, and it still needs explaining.
    ticked = frozenset(point.gap or "" for point in series if point.value is None)
    floor_clause = (
        f" Below this dataset's detection floor ({_fmt(floor)} {unit}) no percentile is "
        "published and the stage is ticked on the axis instead &mdash; below the foot of "
        "the scale, not at a value on it."
        if floor is not None and ProfileGap.BELOW_DETECTION_FLOOR.value in ticked
        else ""
    )
    other_clause = _other_gap_clause(ticked)
    tick_clause = ""
    if other_clause:
        tick_clause = (
            f" A tick can mean something else too: {other_clause}. The mark says only that "
            "this atlas did not place the stage, never that the gene reads low there."
            if floor_clause
            else f" A ticked stage is one this atlas did not place: {other_clause}. The mark "
            "says nothing about how high the gene reads there."
        )
    _, no_row = _line_breaks(series)
    # Two spellings, because the sentence contrasts itself with the clause
    # above it and that clause is not always there. "A break with no tick
    # under it is different again" on a figure drawing no tick at all sends a
    # reader hunting for one -- the defect this whole caption was rewritten
    # for, reintroduced one sentence later (CLAUDE.md section 4.27).
    gap_clause = ""
    if lined and no_row:
        gap_clause = (
            " A break in the line with no tick under it is different again: this dataset has "
            "<strong>no row for this organ at that stage</strong>, so there is no measurement "
            "to place, low or otherwise."
            if floor_clause
            else " The line breaks where this dataset has <strong>no row for this organ at "
            "that stage</strong>: there is no measurement to place, low or otherwise."
        )
    stages = entry["stages"]
    span = (
        f"{html.escape(_stage_label(stages[0]))} to {html.escape(_stage_label(stages[-1]))}, "
        f"all {len(stages)} of this dataset's stages"
        if stages
        else "this dataset's stages"
    )
    return (
        f'<p class="method">Left to right: {span} in curated developmental '
        "order, evenly spaced &mdash; the axis is <strong>order, not elapsed time</strong>. "
        f"Vertical: median abundance in whole {html.escape(tissue)}, {unit}, on a log scale; "
        f"the axis runs {_fmt(low)} to {_fmt(high)} and is fitted to <strong>this gene</strong>, "
        "so a curve's height compares nothing to another gene's."
        f"{floor_clause}{tick_clause}{gap_clause} Shaded: {named} ({attribution}).</p>"
    )


def _adjacent_runs(
    placed: Sequence[tuple[int, float]], points: Sequence[tuple[float, float]]
) -> list[list[tuple[float, float]]]:
    """`points` split wherever the stages behind them are not consecutive.

    A break is any stage the atlas did not place sitting between two it did --
    below this dataset's detection floor, or with no percentile grid behind
    it -- and also a stage where this organ has no row at all. The two are
    different facts and `_trajectory` renders them differently (the first gets
    an axis tick, the second nothing), but they are the same fact *for the
    line*: there is no measurement here to draw through.

    Takes the pixel points beside the stage indices rather than recomputing
    them, so the run boundaries and the markers cannot come from two different
    passes over the series.
    """
    runs: list[list[tuple[float, float]]] = []
    previous: int | None = None
    for (index, _), point in zip(placed, points, strict=True):
        if previous is None or index != previous + 1:
            runs.append([])
        runs[-1].append(point)
        previous = index
    return runs


def _line_breaks(series: Sequence[_Median]) -> tuple[frozenset[str], bool]:
    """Which causes actually interrupt this organ's drawn line:
    `(the recorded reasons for a stage sampled and not placed, whether a stage
    has no row at all)`.

    **The first half is a set of reasons and was a bare `bool`.** "A stage
    sampled and not placed" is not one fact -- below a detection floor is a
    low reading about the gene, a missing percentile grid is a hole in the
    reference -- and a caller handed only `True` can say nothing but the
    first. `_spark_caption` said exactly that, of every break, whatever the
    payload recorded.

    `_adjacent_runs` splits the line on both and cannot tell them apart,
    because for the *line* they are the same fact. For the *caption* they are
    not: the first is a measurement this atlas declines to vouch for and gets
    a tick on the axis, the second is the dataset never having looked, and
    `_tissue_medians`' docstring is explicit that they "must not render the
    same way". They did -- to a reader, a break with a tick and a break
    without one differed only in a mark the caption promised on every page
    and drew on 30 of 85.

    Only a gap *between two placed stages* counts, for either cause. A stage
    the dataset skipped before the first placed one or after the last breaks
    no line, so nothing on the figure needs explaining.
    """
    recorded = {point.position: point.gap for point in series}
    placed = [point.position for point in series if point.value is not None]
    unplaced: set[str] = set()
    unsampled = False
    for first, second in zip(placed, placed[1:], strict=False):
        for gap in range(first + 1, second):
            if gap in recorded:
                unplaced.add(recorded[gap] or "")
            else:
                unsampled = True
    return frozenset(unplaced), unsampled


def _trajectory(
    entry: DatasetProfileEntry,
    tissue: str,
    dataset: Dataset | None,
    phases: CardiacPhaseFile | None,
) -> str:
    """One organ's whole developmental series as a phase-banded figure.

    Returns `""` -- the caller's cue to say something in words instead -- in
    each of the three cases where a figure would assert more than the data
    does: nothing sampled for this organ, nothing placed at any stage, and no
    curated phase band to draw, which is the bare re-plot D43 forbids.

    **The three are different facts and `_dataset_block` owes a reader a
    different sentence for each**, which it could not do while this returned
    a bare `""` for all three. The band question is answered per *dataset*
    rather than per organ -- `_phase_bands` never reads the tissue -- so the
    caller can ask it directly, and does; what it cannot re-derive cheaply is
    which of the first two applies, so it asks `_anything_placed` instead.

    The line is drawn only from `_STAGES_FOR_A_TRAJECTORY` placed points up.
    Below that the markers stand alone: `charts.polyline` refuses fewer than
    two points outright, and two would draw a segment claiming something about
    the interval between them that two measurements do not support.

    **And it is drawn in runs of adjacent stages, never as one line through
    every placed point.** Measured 2026-08-20: 20 of the 85 charted genes have
    a below-floor stage sitting strictly between two plotted points, and TBX1
    has 11 of them against 5 placed stages. One line through the placed points
    would run smooth and high across exactly the stages where this dataset
    measured below its own floor -- the picture asserting a continuity the
    figures deny, which is the conflation `_percentile_cell` refuses in words
    and `_no_trajectory_sentence` refuses for a whole gene.
    """
    series = _tissue_medians(entry, tissue)
    placed = [(point.position, point.value) for point in series if point.value is not None]
    if not placed:
        return ""
    count = len(entry["stages"])
    bands = _phase_bands(entry, phases, count)
    if not bands:
        return ""

    low, high = _axis_bounds([value for _, value in placed])
    scale = LogScale(low=low, high=high, left=0.0, width=_PLOT_BOTTOM - _FLOOR_GAP - _PLOT_TOP)
    points = [
        (_stage_x(index, count), _PLOT_BOTTOM - _FLOOR_GAP - scale.x(value))
        for index, value in placed
    ]
    axis = (
        f'<line class="chart-axis" x1="{coordinate(_PLOT_LEFT)}" '
        f'y1="{coordinate(_PLOT_BOTTOM)}" x2="{coordinate(_PLOT_RIGHT)}" '
        f'y2="{coordinate(_PLOT_BOTTOM)}"/>'
        + _axis_label(high, _PLOT_BOTTOM - _FLOOR_GAP - scale.x(high))
        + _axis_label(low, _PLOT_BOTTOM - _FLOOR_GAP - scale.x(low))
    )
    ticks = "".join(
        marker(_stage_x(point.position, count), _PLOT_BOTTOM, css_class="chart-absent")
        for point in series
        if point.value is None
    )
    enough = len(points) >= _STAGES_FOR_A_TRAJECTORY
    lines = (
        "".join(
            polyline(run, css_class="chart-line")
            for run in _adjacent_runs(placed, points)
            if len(run) >= 2
        )
        if enough
        else ""
    )
    markers = "".join(marker(x, y, css_class="chart-point") for x, y in points)
    unit = _tissue_unit(entry, tissue)
    floor = dataset.detection_floor if dataset is not None else None
    floor_words = f", above a detection floor of {_fmt(floor)} {unit}" if floor is not None else ""
    title = (
        f"Median abundance in whole {tissue}, {series[0].label} to {series[-1].label}: "
        f"{len(placed)} of {len(series)} sampled developmental stages placed against this "
        f"dataset's percentile grid{floor_words}, banded by this atlas's own curated "
        "cardiac phases."
    )
    figure = svg_figure(
        width=_CHART_WIDTH,
        height=_CHART_HEIGHT,
        title=title,
        body=bands + axis + ticks + lines + markers,
    )
    # `lined` rather than a second look at the data: the caption's clause
    # about a break in the line may only be published if there *is* a line,
    # and only this function knows whether one was drawn. Three of the 85
    # charted genes place too few stages for one.
    return figure + _band_caption(
        entry, phases, tissue, floor, low, high, series, lined=bool(lines)
    )


def _sampled_cardiac(entry: DatasetProfileEntry, cardiac: frozenset[str]) -> list[tuple[str, int]]:
    """Each cardiac organ this dataset sampled for this gene, and how often.

    Sorted, because `cardiac` is a `frozenset` and an unsorted iteration would
    reorder a published sentence between two builds of one commit.
    """
    return [
        (tissue, len(_tissue_medians(entry, tissue)))
        for tissue in sorted(cardiac)
        if _tissue_medians(entry, tissue)
    ]


def _sampled_organs_clause(sampled: Sequence[tuple[str, int]]) -> str:
    """ "whole heart (19 stages sampled)" -- the organ clause both refusals open with.

    One implementation, because the two sentences name the same organs for
    two different reasons and a reader who sees one never sees the other; two
    copies of this could spell the same fact two ways with nothing failing.

    **`stage` has a singular form here and did not before.** The count is a
    real 1 on a live corpus shape -- a dataset sampling one cardiac stage for
    a gene -- and "1 stages sampled" is the `_abundance`/`_count` discipline
    lapsing in the one place nothing was reading.
    """
    return ", ".join(
        f"whole {html.escape(tissue)} ({count} stage{'' if count == 1 else 's'} sampled)"
        for tissue, count in sampled
    )


def _unbanded_sentence(entry: DatasetProfileEntry, cardiac: frozenset[str]) -> str:
    """D43's refusal, in words: placed measurements but no phase to band them with.

    The third of `_trajectory`'s three `""`s, and the only one that is not
    about the gene at all. `_phase_bands` returns nothing when no curated
    cardiac phase with a stated end covers any of this dataset's stages, and
    a curve drawn without one is the re-plot of the source's own browser that
    D43 forbids. That is a fact about this atlas's vocabulary, so this
    sentence claims nothing about the measurements -- every one of which is
    placed, and printed in the table below.

    Reached only where something *is* placed. Where nothing is,
    `_no_trajectory_sentence` has the better sentence whatever the bands say:
    "no trajectory" is then a fact about the gene, and the missing band is a
    second reason a reader does not need.

    Returns `""` when this dataset sampled no cardiac organ for this gene,
    the same degenerate case `_no_trajectory_sentence` returns `""` for and
    for the same reason: there is no organ to make a statement about.
    """
    sampled = _sampled_cardiac(entry, cardiac)
    if not sampled:
        return ""
    return (
        f'<p class="notice-inline">No trajectory is drawn for '
        f"{_sampled_organs_clause(sampled)}: none of this dataset's developmental stages "
        "falls inside a curated cardiac phase this atlas can band, and an unbanded curve "
        "is the source's own plot rather than anything this atlas adds. Every per-stage "
        "figure is in the table below.</p>"
    )


def _anything_placed(entry: DatasetProfileEntry, cardiac: frozenset[str]) -> bool:
    """Does any cardiac organ carry a median this dataset placed?

    The one question that tells `_trajectory`'s second refusal from its third,
    and it is asked of every cardiac organ at once because the answer decides
    one sentence for the whole dataset block.
    """
    return any(
        point.value is not None
        for tissue in sorted(cardiac)
        for point in _tissue_medians(entry, tissue)
    )


def _no_trajectory_sentence(
    entry: DatasetProfileEntry, cardiac: frozenset[str], floor: float | None
) -> str:
    """D42: what a page says instead of drawing an empty chart.

    Measured 2026-08-20 against a real build: seven published genes are below
    the floor in every cardiac tissue at every sampled stage -- SEMA3E, ZIC3,
    USP44, DAW1, FGF8, GDF1, NODAL. Their headline grades are `definitive`
    (ZIC3, NODAL), `moderate` (DAW1, FGF8), `limited` (USP44) and none at all
    (SEMA3E, GDF1). They are laterality and outflow-tract genes acting in rare
    lineages, which is exactly the dilution `_BULK_DILUTION_NOTICE` describes.
    An empty chart beside a definitive chip asserts the thing that caption
    denies, so this tier states the fact and repeats the caveat rather than
    showing a frame with nothing in it.

    **Two branches, because "nothing placed" is not one fact.** All seven of
    those genes are below the floor, and that is what earns the dilution
    argument; a cell left unplaced because no percentile grid was published,
    or because the dataset declares no floor, is a gap in the *reference*
    rather than a low reading, and saying "below the detection floor" of it
    would be false. The second branch names the recorded reason instead and
    claims nothing about the gene.

    Returns `""` when this dataset sampled no cardiac organ for this gene at
    all: there is then no organ to make either statement about, and every
    figure it did publish is in the table below regardless.

    **Reached only where nothing at all is placed, and the first branch is
    false the moment that stops holding.** "At every stage sampled" is read
    off the set of recorded gap *reasons*, which says what the gaps were and
    not how many stages had one; a series placed at 4 and 6 wpc and below the
    floor at 5 has exactly one reason in that set. It published the bold
    claim of a gene whose page printed two percentiles three lines below it
    (2026-08-20). `_dataset_block` is what holds the precondition, by asking
    `_anything_placed` before choosing this sentence over
    `_unbanded_sentence`; the guard is stated here because this is where a
    later author would have to break it.
    """
    sampled = _sampled_cardiac(entry, cardiac)
    if not sampled:
        return ""
    organs = _sampled_organs_clause(sampled)
    named = {tissue for tissue, _ in sampled}
    reasons = {
        measured["not_placed_reason"]
        for stage in entry["stages"]
        for measured in stage["tissues"]
        if measured["tissue"] in named and measured["placement"] is None
    }
    if reasons != {ProfileGap.BELOW_DETECTION_FLOOR.value}:
        clauses = "; ".join(
            _PLACEMENT_GAP_CLAUSE.get(reason or "", "no reason was recorded for this gap")
            for reason in sorted(reason or "" for reason in reasons)
        )
        return (
            f'<p class="notice-inline">No trajectory is drawn for {organs}: no measurement '
            f"there is placed against this dataset's percentile grid &mdash; {clauses}. Every "
            "per-stage figure is in the table below.</p>"
        )
    # The floor's own value is a separate sentence rather than a parenthesis
    # inside the first, because the organ clause already ends in one -- "in
    # whole heart (19 stages sampled) (detection floor 1 tpm)" is what nesting
    # them produced, and a reader stops reading at the second bracket.
    qualifier = (
        f" This dataset's detection floor is {_fmt(floor)} "
        f"{html.escape(_dominant_unit(entry, named))}."
        if floor is not None
        else ""
    )
    return (
        f'<p class="notice-inline">This gene reads <strong>below this dataset\'s detection '
        f"floor at every stage sampled</strong> in {organs}, so there is no trajectory to "
        f"draw.{qualifier} That is "
        "<strong>not evidence that</strong> the gene is unimportant to heart development, and it "
        "is not a measured zero: every figure here is a whole-organ median, and a gene "
        "transcribed intensely in one rare cardiac lineage still reads below a bulk floor once "
        "diluted across the whole organ. Every per-stage figure is in the table below.</p>"
    )


def _dominant_unit(entry: DatasetProfileEntry, tissues: set[str]) -> str:
    """The unit the named organs' medians are quoted in.

    A floor quoted without one is the `_count`/`count_unit` defect in a new
    place; the organs here all come from one dataset, so one unit answers for
    all of them, and the first is read rather than a set assembled.
    """
    for tissue in sorted(tissues):
        unit = _tissue_unit(entry, tissue)
        if unit:
            return unit
    return ""


# --- The organ small multiples ----------------------------------------------

# Small on purpose: this is a grid of up to seven panels side by side
# (`.sparks` in `render.STYLESHEET`), and the *comparison between* them is the
# figure. No panel here is meant to be read on its own -- the trajectory above
# is what reads one organ in detail.
_SPARK_WIDTH: Final = 96
_SPARK_HEIGHT: Final = 40
_SPARK_LEFT: Final = 4.0
_SPARK_RIGHT: Final = 92.0
_SPARK_TOP: Final = 5.0
_SPARK_BOTTOM: Final = 35.0


def _spark_x(index: int, count: int) -> float:
    """Where the `index`-th of `count` stages sits inside one small panel.

    The same ordinal axis `_stage_x` documents, at a twentieth of the width.
    """
    return _even_x(index, count, _SPARK_LEFT, _SPARK_RIGHT)


def _sampled_tissues(entry: DatasetProfileEntry) -> list[str]:
    """Every organ this dataset sampled for this gene, sorted.

    Sorted because it is derived from a `set`: `sort_keys` in the JSON encoder
    has nothing to say about the order panels are laid out in, and two builds
    of one commit must emit the same bytes.
    """
    return sorted({cell["tissue"] for stage in entry["stages"] for cell in stage["tissues"]})


def _small_multiples(entry: DatasetProfileEntry, dataset: Dataset | None) -> str:
    """One sparkline per organ, every panel on a single shared abundance axis.

    **The shared axis is the whole point.** τ answers "is this gene
    concentrated in one organ, and which" -- and this section stated it as a
    number restated once per stage, 19 times on the TBX5 page (its other two
    stages carry no τ) with only the figure changing. Per-organ axes would
    rescale each line to its own range and draw five similar-looking traces,
    which is the same non-answer in a picture.

    One axis is what makes a τ of 0.80 to 1.00 -- TBX5's range across those 19
    stages -- legible as a place rather than a number. Measured 2026-08-20 on
    the committed corpus: TBX5's heart medians run 5 to 275 tpm, and the four
    non-heart organs it places anything in peak at forebrain 2, ovary 3,
    testis 4 and liver 12, so heart's peak sits 22.9x -- 1.36 decades -- above
    the highest of them. **Not "above everything else": heart's own minimum,
    5 tpm at middle adult, is below liver's 12.** The panels are what shows
    that, and seven private axes would show none of it. TBX5 draws five of
    them; hindbrain and kidney place nothing and the caption names them.

    `data-scale-high` is emitted from the axis itself rather than recomputed,
    so a build that gave each organ its own scale would publish a different
    value per panel and could not pass the test asserting they are equal.
    Assert the property, never the appearance.

    **What earns this under D43** -- never a panel redundant with the source's
    own browser -- is `Dataset.cardiac_tissues`. Which organ counts as cardiac
    is this atlas's own declaration in `curation/datasets.yaml`, not the
    source's; the source's browser has the per-organ curves and makes no such
    claim. The caption states that in words, because a stroke weight is not an
    argument.

    Returns `""` in the two cases where a panel would assert more than the
    data does, D42's rule that a page says so rather than drawing an empty
    frame:

    - **No stage reports a τ at all.** τ needs at least two organs, and with
      one sampled there is no specificity to draw a specificity panel about.
      `_SPECIFICITY_GAP_CLAUSE` is what the page says instead, once per stage,
      and it says it already.
    - **No organ has a single placed median.** Nothing to plot.

    **Only medians this dataset placed are drawn**, which is `_tissue_medians`'
    discipline reused rather than restated: a below-floor figure is one this
    atlas declines to vouch for, and putting it on the same line as ones it
    does is the conflation `_percentile_cell` refuses in words. The line
    therefore breaks at an unplaced stage exactly as `_trajectory`'s does --
    measured 2026-08-20, 20 of the 85 charted genes have a below-floor stage
    strictly between two plotted points, and a guard added to one layer is not
    a guard.

    An organ with nothing placed gets no panel and is **named in the caption**
    instead. Silently dropping it would leave five panels beside a τ sentence
    reading "7 organs sampled", the page contradicting itself, and a reader
    concluding this dataset never looked there.
    """
    if all(stage["specificity"] is None for stage in entry["stages"]):
        return ""

    series: dict[str, list[tuple[int, float]]] = {}
    sampled_at: dict[str, int] = {}
    # Kept whole, not just counted: `_spark_caption` has to tell a stage this
    # dataset sampled and did not place from one it has no row for at all,
    # and by the time a series is filtered to its placed points the two are
    # the same absence. Handed down rather than re-derived, the reason
    # `_adjacent_runs` takes its points beside its indices.
    wholes: dict[str, list[_Median]] = {}
    for tissue in _sampled_tissues(entry):
        whole = _tissue_medians(entry, tissue)
        sampled_at[tissue] = len(whole)
        wholes[tissue] = whole
        series[tissue] = [
            (point.position, point.value) for point in whole if point.value is not None
        ]
    drawn = [tissue for tissue in series if series[tissue]]
    if not drawn:
        return ""

    low, high = _axis_bounds([value for tissue in drawn for _, value in series[tissue]])
    scale = LogScale(low=low, high=high, left=0.0, width=_SPARK_BOTTOM - _SPARK_TOP)
    cardiac = frozenset(dataset.cardiac_tissues) if dataset is not None else frozenset()
    count = len(entry["stages"])

    panels = ""
    for tissue in drawn:
        placed = series[tissue]
        points = [
            (_spark_x(index, count), _SPARK_BOTTOM - scale.x(value)) for index, value in placed
        ]
        css = "chart-cardiac" if tissue in cardiac else "chart-control"
        # The same floor `_trajectory` applies, and for the same stated
        # reason: two measurements do not support a claim about the interval
        # between them. This panel drew that segment anyway until 2026-08-21,
        # so the two figures on one page published opposite claims about one
        # pair of numbers -- measured, 13 organ series, and on CFC1 the
        # trajectory drew two bare markers while the heart panel beside it,
        # drawn heavier as the declared cardiac tissue, joined them.
        #
        # Gated on the organ's whole placed series rather than per run, again
        # matching `_trajectory`: a run of two inside a longer series is a
        # segment between two adjacent measurements of a trend the series
        # already evidences, which is a different claim from a two-point
        # series being a trend on its own.
        enough = len(placed) >= _STAGES_FOR_A_TRAJECTORY
        body = "".join(
            polyline(run, css_class=css)
            if enough and len(run) >= 2
            else "".join(marker(x, y, css_class=css) for x, y in run)
            for run in _adjacent_runs(placed, points)
        )
        role = (
            ", this atlas's declared cardiac tissue for this dataset" if tissue in cardiac else ""
        )
        title = (
            f"Median abundance in whole {tissue}{role}: {len(placed)} of "
            f"{sampled_at[tissue]} sampled developmental stages placed, on the abundance "
            "axis every organ here shares."
        )
        panels += (
            f'<figure class="spark" data-scale-high="{coordinate(scale.high)}">'
            + svg_figure(width=_SPARK_WIDTH, height=_SPARK_HEIGHT, title=title, body=body)
            + f"<figcaption>{html.escape(tissue)}</figcaption></figure>"
        )
    undrawn = [tissue for tissue in series if not series[tissue]]
    caption = _spark_caption(entry, drawn, undrawn, cardiac, low, high, wholes)
    return f'<div class="sparks">{panels}</div>{caption}'


def _spark_caption(
    entry: DatasetProfileEntry,
    drawn: Sequence[str],
    undrawn: Sequence[str],
    cardiac: frozenset[str],
    low: float,
    high: float,
    wholes: Mapping[str, Sequence[_Median]],
) -> str:
    """What the grid of panels cannot say about itself, in words beside it.

    Four things, and three of them are traps. The horizontal axis is stage
    *order* rather than elapsed time, exactly as `_band_caption` says of the
    trajectory. The vertical axis is **shared**, which is the entire claim the
    figure makes and the one a reader cannot check by looking. The highlighted
    organ is this atlas's own declaration and not the source's, which is what
    admits the panel under D43. And an organ with no panel was still sampled.

    **The opening clause said "one panel per organ this dataset sampled for
    this gene", and 25 of the 91 pages carrying it went on to say "No panel is
    drawn for ..." three sentences later.** A page that corrects itself is not
    a page that is right: the first sentence is what a reader takes the grid
    to mean, and an organ dropped for having nothing placed is exactly the
    organ they would otherwise conclude was never looked at. It now says what
    the grid actually is -- one panel per organ with a *placed* measurement --
    so the sentence naming the dropped organs adds to it rather than
    contradicting it.

    **A line breaks for two reasons and only one was named.** The caption
    attributed every gap to the detection floor; a panel's line also breaks
    where this dataset has no row for that organ at that stage, and in a
    spark there is no tick, so the two are pixel-for-pixel identical.
    Measured 2026-08-20 on the built corpus: 87 of the 92 genes have at least
    one panel split by a no-row stage -- testis on 82 pages, heart on 72,
    forebrain on 69. Each cause is now its own clause, conditional on
    occurring in this gene's panels, which is `_line_breaks`' whole purpose.

    **The cardiac clause has three positions and no blank one.** A panel where
    the declared cardiac tissue has nothing placed carries no highlight, and
    "Highlighted: heart" beside no highlighted panel is false -- a sentence's
    truth conditions travel with its position (`_POOLING_NOTICE`'s recorded
    lesson). Measured 2026-08-20 on the committed corpus, that is 6 of the 91
    genes with a panel, so it is a live case. It names the declaration anyway,
    because "this atlas calls heart the cardiac tissue and placed nothing
    there" is the fact, and dropping the sentence would leave a reader to
    conclude the atlas never declared one.
    """
    unit = html.escape(_dominant_unit(entry, set(drawn)))
    interrupted = [_line_breaks(wholes[tissue]) for tissue in drawn]
    unplaced: frozenset[str] = frozenset().union(*(reasons for reasons, _ in interrupted))
    no_row = any(unsampled for _, unsampled in interrupted)
    causes = []
    if ProfileGap.BELOW_DETECTION_FLOOR.value in unplaced:
        causes.append(
            "a stage measured <strong>below this dataset's detection floor</strong>, which "
            "this atlas does not place"
        )
    other = _other_gap_clause(unplaced)
    if other:
        causes.append(f"a stage this atlas published no percentile for &mdash; {other}")
    if no_row:
        causes.append(
            "a stage where this dataset has <strong>no row for that organ</strong> and "
            "nothing was measured to place"
        )
    breaks = (
        " Only medians this dataset placed against its percentile grid are drawn, so a line "
        "breaks rather than running through " + ", and through ".join(causes) + "."
        if causes
        else " Only medians this dataset placed against its percentile grid are drawn."
    )
    cardiac_drawn = [tissue for tissue in drawn if tissue in cardiac]
    if cardiac_drawn:
        marked = ", ".join(html.escape(tissue) for tissue in cardiac_drawn)
        noun = "tissue" if len(cardiac_drawn) == 1 else "tissues"
        highlight = (
            f" Drawn heavier: {marked} &mdash; this atlas's own declaration of which organ "
            f"is this dataset's cardiac {noun}, not the source's."
        )
    elif cardiac:
        declared = ", ".join(html.escape(tissue) for tissue in sorted(cardiac))
        noun = "tissue" if len(cardiac) == 1 else "tissues"
        highlight = (
            f" No panel here is drawn heavier: this atlas declares {declared} this "
            f"dataset's cardiac {noun}, and this dataset placed no measurement there."
        )
    else:
        highlight = ""
    if undrawn:
        names = ", ".join(html.escape(tissue) for tissue in undrawn)
        organs = "that organ" if len(undrawn) == 1 else "those organs"
        missing = (
            f" No panel is drawn for {names}: this dataset sampled {organs} and placed no "
            "measurement there against its percentile grid."
        )
    else:
        missing = ""
    return (
        '<p class="method">One panel per organ with a placed measurement for this gene, each '
        "drawn left to right in curated developmental order &mdash; the horizontal axis is "
        "<strong>order, not elapsed time</strong>. <strong>Every panel shares one abundance "
        f"axis</strong> ({_fmt(low)} to {_fmt(high)} {unit}, log scale), which is what makes "
        "a gene concentrated in one organ look different from a broadly expressed one; a "
        "per-organ axis would rescale every line to its own range and hide exactly that "
        f"difference.{highlight}{breaks}{missing} Every figure is in the table below.</p>"
    )


def _dataset_block(
    entry: DatasetProfileEntry, dataset: Dataset | None, phases: CardiacPhaseFile | None
) -> str:
    """One dataset's whole contribution: a chart per cardiac organ, one panel per
    organ sampled, then every figure.

    Two pictures, answering two questions, and neither substitutes for the
    other. `_trajectory` reads one cardiac organ in detail across development;
    `_small_multiples` compares every organ on one shared axis, which is what
    τ asserts once per stage and nothing on the page showed. Each returns `""`
    on its own terms, so a dataset can carry either, both or neither.

    The link is D39(b)'s other half reaching a reader rather than only a
    program: `build_profile_quantiles` publishes the 101-point grid a
    percentile above was read against, and this is the one place a person
    reading the page -- not a script reading the bundle -- can reach it.
    Omitted when `quantile_shard` is `None` (a corpus mid-curation, or the
    gate bypassed) rather than linking a file `Emitter` never wrote.

    `dataset` is looked up by the caller and may be `None` --
    `ProfileGap.DATASET_NOT_REGISTERED` is a real, named state, not a caller
    error, so this degrades to no cardiac tissues and no floor rather than
    raising: the gap is already stated on every cell it affects, and this
    function's job is to render what the payload says, not to police a
    registration mismatch a validator already reports elsewhere.

    **The stage blocks move into a `<details>`; they do not leave.** This
    section rendered 21 near-identical blocks per dataset, one per stage, and
    59% of its text was verbatim repetition -- "below the detection floor in
    whole liver at this stage (detection floor 1 tpm)" appeared 20 times on
    the TBX5 page.

    **A bare percentage is not reproducible, so here is the method that
    produced that one.** Build `dev`, take the TBX5 page from
    `<h2>Developmental expression</h2>` to the end, replace every tag with a
    newline (which puts each table cell on its own line), strip and drop
    blanks, and count the words in lines that duplicate an earlier line:
    1,854 of 3,136 words over 513 lines, 146 of them distinct -- 59.1%. Other
    ways of asking the same page give very different answers, none of them
    wrong: 12% of its characters sit inside a repeated *sentence*, and 89% of
    its tokens sit inside a repeated 4-gram. The figure is about *lines*, and
    the lines are mostly table cells.

    A chart that *replaced* them would take every exact figure out of the
    HTML, which is this repository's characteristic defect (curated work
    reaching no page) wearing a redesign. The summary is what a reader sees
    first; the record is one click away and still complete.
    """
    heading = f"<h3>{html.escape(entry['dataset'])}"
    shard = entry["quantile_shard"]
    if shard:
        heading += f' <a href="../{html.escape(shard)}">percentile grid as JSON</a>'
    heading += "</h3>"
    cardiac = frozenset(dataset.cardiac_tissues) if dataset is not None else frozenset()
    floor = dataset.detection_floor if dataset is not None else None
    # Sorted: `cardiac` is a `frozenset`, and an unsorted iteration would put
    # two organs' charts in a different order between two builds of one commit.
    charts = "".join(_trajectory(entry, tissue, dataset, phases) for tissue in sorted(cardiac))
    # **Three refusals, three sentences.** `_trajectory` returns `""` for
    # three different reasons -- nothing sampled, nothing placed, and no
    # curated phase band -- and a caller that routes them all into one
    # sentence publishes a claim that is false for two of them. Measured
    # 2026-08-20 with `runner.py` passing `phases=None`: 55 pages said "no
    # measurement there is placed against this dataset's percentile grid
    # &mdash; ." with the reason list empty, and 30 said, in bold, that the
    # gene reads below the detection floor at every stage sampled -- one of
    # them with 18 of its 19 heart stages placed and printed in the table
    # below. Nothing was published in either state, because `runner.py`
    # passes the real vocabulary; the sentences were wrong all the same, and
    # a dataset whose stages fall in no *ended* cardiac phase reaches them.
    #
    # Placement first, because it is the fact about the gene. A missing band
    # is a fact about this atlas's own vocabulary, and where nothing is
    # placed a reader is better served by the dilution caveat than by it.
    if charts:
        lede = charts
    elif _anything_placed(entry, cardiac):
        lede = _unbanded_sentence(entry, cardiac)
    else:
        lede = _no_trajectory_sentence(entry, cardiac, floor)
    sparks = _small_multiples(entry, dataset)
    stages = "".join(_stage_block(stage, cardiac, floor) for stage in entry["stages"])
    total = len(entry["stages"])
    noun = "stage" if total == 1 else "stages"
    folded = (
        f'<details class="stage-figures"><summary>every figure, all {total} {noun}'
        f"</summary>{stages}</details>"
    )
    return heading + lede + sparks + folded


def _placement_count(entries: Sequence[DatasetProfileEntry]) -> int:
    """How many organ/stage/dataset cells actually carry a percentile.

    What `_PERCENTILE_COMPARABILITY_NOTICE` is conditioned on: a page with at
    most one placement has nothing side by side for a reader to compare.
    """
    return sum(
        1
        for entry in entries
        for stage in entry["stages"]
        for tissue in stage["tissues"]
        if tissue["placement"] is not None
    )


def _expression_section(
    profile: ExpressionProfile,
    datasets: Mapping[str, Dataset],
    phases: CardiacPhaseFile | None,
) -> str:
    """Developmental expression: percentile, tau and phase, per organ and stage.

    **Renders a section even for a gene with no data, unlike `_burden_
    section`.** A burden table's absence would misstate a null result as
    evidence against a gene it was never tested for, so that section returns
    `""`. There is no equivalent risk here: `datasets == []` states a fact this
    atlas already knows -- no curated dataset mentions this gene -- and
    omitting the section would make that indistinguishable from the atlas
    never having built this feature at all. `_NOT_CURATED_EXPRESSION` follows
    `_not_curated`'s idiom, not `_burden_section`'s.

    `datasets` resolves each entry's own accession to the curated `Dataset` it
    came from, for `cardiac_tissues` (the argmax gate) and `detection_floor`
    (rule 5's required value) -- both dataset-level facts this per-gene
    payload does not repeat. A gene absent from the caller's own `profiles`
    mapping is the caller's job to default to `EMPTY_EXPRESSION_PROFILE`,
    matching `bundles._expression_profile`'s own fallback, so the page and the
    bundle cannot disagree about which genes carry data.
    """
    entries = profile["datasets"]
    if not entries:
        return "<h2>Developmental expression</h2>" + _NOT_CURATED_EXPRESSION

    blocks = "".join(
        _dataset_block(entry, datasets.get(entry["dataset"]), phases) for entry in entries
    )
    notes = _EXPRESSION_READING_NOTES
    if _placement_count(entries) > 1:
        notes += _PERCENTILE_COMPARABILITY_NOTICE
    reading = (
        '<details class="reading-notes"><summary>How to read these figures</summary>'
        f"{notes}</details>"
    )
    return "<h2>Developmental expression</h2>" + _BULK_DILUTION_NOTICE + reading + blocks


def build_gene_pages(
    facts: Mapping[str, GeneFacts],
    emitter: Emitter,
    *,
    symbols: Mapping[str, str],
    validity: Mapping[str, GeneValidity],
    assertions: Mapping[str, list[LesionAssertion]],
    publications: Mapping[str, Publication],
    burden: Mapping[str, Sequence[BurdenRow]],
    cohorts: Mapping[str, Cohort],
    families: tuple[frozenset[str], ...] = (),
    axes: tuple[tuple[str, str], ...] = (),
    profiles: Mapping[str, ExpressionProfile] | None = None,
    datasets: Mapping[str, Dataset] | None = None,
    phases: CardiacPhaseFile | None = None,
) -> None:
    """Emit one HTML page per gene in `facts`.

    `facts` is `build_genes`' return, so the pages render from the same
    derivation the bundles were written from rather than deriving it twice --
    a page and its bundle disagreeing about a gene's confidence is the failure
    this shared input removes.

    `symbols` is the same `{gene: labels.symbol}` projection `build_genes` and
    `build_landing` take; a gene absent from it reads as its own HGNC id, the
    same registry fallback those two apply.

    A gene absent from `validity` gets `uncurated()` -- an empty table under its
    own heading, not a missing section, for the reason `render.data_table`
    renders a header over no rows at all.

    `profiles` and `datasets` both default to `None` -- and to `{}` inside --
    rather than to `{}` in the signature, the same reason `concordance` on
    `build_gene_index_page` does: an empty *literal* dict default is one object
    shared across every call that omits the argument, which is harmless only as
    long as nobody ever mutates it, and this project does not rely on that
    holding forever. A gene absent from `profiles` gets `EMPTY_EXPRESSION_
    PROFILE`, matching `bundles._expression_profile`'s own fallback exactly, so
    the page and the bundle cannot disagree about which genes carry a
    developmental expression profile. `datasets` resolves an entry's own
    accession to the curated record `_expression_section` reads
    `cardiac_tissues` and `detection_floor` from; see that function's docstring
    for why a missing dataset degrades rather than raises.

    `phases` is the same `CardiacPhaseFile` the profiles were derived against,
    threaded down rather than re-read, so a band on a chart and the phase
    sentence under the stage it covers cannot name different vocabularies.
    `None` -- no curated phase file at all -- costs the trajectory, not the
    section: `_trajectory` refuses to draw an unbanded curve, because a curve
    with no band is the panel the source's own browser already publishes, and
    every figure it would have summarised is in the tables regardless.

    Sorted, like every loop in this build that iterates a mapping: `sort_keys`
    orders dict keys in a JSON payload and has nothing to say about the order
    files are written in, and a differing write order would move nothing here
    today but is the habit this project keeps.
    """
    profile_by_gene = profiles or {}
    dataset_registry = datasets or {}
    for gene in sorted(facts):
        fact = facts[gene]
        symbol = symbols.get(gene, gene)
        records = assertions.get(gene, [])
        gene_validity = validity.get(gene, uncurated())
        body = (
            '<div class="layout">'
            + _rail(
                gene,
                symbol,
                fact,
                _graded_diseases(gene_validity, fact.headline_confidence),
                agreeing_submitters(gene_validity),
                len(burden.get(gene, ())),
            )
            + "<div>"
            + _SCOPE_RULE
            + grade_legend()
            + _validity_section(gene_validity)
            + (
                _evidence_section(records, publications)
                if fact.atlas_curation is AtlasCuration.CURATED
                else _not_curated(fact)
            )
            # After the atlas's own section, not between it and the mirrored
            # validity table. `_not_curated` says "the classification above is
            # an expert panel's", and a burden table sitting above that sentence
            # would put nine rows of statistics -- which are not classifications
            # -- inside what it refers to. It also keeps the atlas's own curation
            # adjacent to the notice about whether there is any.
            + _burden_section(burden.get(gene, ()), publications, cohorts, families, axes)
            # Last: burden and expression are both mirrored/derived evidence
            # layers rather than the atlas's own curation, and expression is
            # the newer of the two -- appending keeps every existing section's
            # position, and therefore every existing test slicing this page by
            # position, unchanged.
            + _expression_section(
                profile_by_gene.get(gene, EMPTY_EXPRESSION_PROFILE), dataset_registry, phases
            )
            + "</div></div>"
        )
        emitter.write_text(
            gene_page_path(HgncId(gene)),
            document(title=f"{symbol} — CHD Atlas", root="../", body=body),
        )


def build_gene_index_page(
    facts: Mapping[str, GeneFacts],
    emitter: Emitter,
    *,
    symbols: Mapping[str, str],
    validity: Mapping[str, GeneValidity],
    burden_counts: Mapping[str, int],
    concordance: Mapping[str, Mapping[str, Json]] | None = None,
) -> None:
    """Emit `genes/index.html`: every published gene, filterable in the browser.

    `burden_counts` is `{gene: len(rows)}`, the same number `genes/index.json`
    publishes as `burden_row_count`. The browse page did not mention burden at
    all -- no column, no filter, no count -- while the payload behind it did, so
    a reader scanning 23 identical `definitive` rows could not tell that four of
    them (ISL1, NR2F2, RBM10, SMAD2) carry no loss-of-function evidence at all.
    Raised by review 2026-08-05 and measured on the built site.

    Every row is rendered here, by the build. The inline script only *hides*
    rows, so the page is complete to `curl`, to a crawler, and to a reader with
    JavaScript disabled -- the alternative, an empty `<tbody>` populated from
    `genes/index.json` at runtime, looks identical in a browser and serves
    nothing to any of the three.

    Sorted by HGNC id, the order `genes/index.json` publishes in, so the browse
    page and the payload behind it agree on which gene comes first without
    either reading the other.

    **The attribution line is not decoration.** This page renders 23 rows
    reading `definitive` under a column headed `confidence` (measured
    2026-08-04), and it was the only page on the site where that could be read
    as the atlas's own call: `_validity_section` says "Every classification
    below is an upstream panel's or submitter's" on each gene page, and
    `landing.py` says "The atlas authors no validity classification of its own"
    on the front page. The sentence below says the same thing in the same voice,
    and names `atlas curation` as the column that answers the question the
    confidence column does not -- whether the atlas has curated the gene at all.

    **Whose call it is was only half of it.** All 23 of those rows read
    `definitive` with no disease anywhere on the page, and `_SCOPE_RULE` plus the
    `graded for` column are the other half: a row reading `KMT2D | definitive` on
    a site titled "CHD Atlas" claims KMT2D is definitive for congenital heart
    disease, and ClinGen's actual assertion is `KMT2D -- Definitive for Kabuki
    syndrome 1`, made by the SCID-CID GCEP. `validity` is taken for that column
    alone; `facts` carries no disease, and deriving one from a second read of the
    mirrors would let this page and the gene page it links to name different
    diseases for the same gene.

    Since the gate widened the column answers a second question too. A blank
    `confidence` cell now means something specific -- no ClinGen panel graded the
    gene, which is why `_graded_diseases` returns nothing for it -- and an em
    dash beside a blank is the only pair on this table that is not a gap in the
    data. The paragraph under the heading says so, and `grade_legend` glosses it
    on both pages.

    Two couplings to `render.FILTER_SCRIPT` that nothing in the type system
    holds, both pinned by
    `test_every_facet_names_a_data_attribute_the_filter_script_reads`:

    * each `<select>`'s `name` is resolved into an attribute by concatenation,
      `row.getAttribute('data-' + select.name)`, so the four names below must
      spell the four `attributes` keys exactly. A mismatch reads `null` and
      hides every row the moment that facet is used, while the page looks
      correct until someone touches it.
    * the script matches a facet with `.split(' ').indexOf(want)`, so a value
      carrying a space could never be selected. None can: measured 2026-08-04,
      no member of `LesionGroup`, `Classification`, `ValidityState` or
      `AtlasCuration` contains one. That is what lets `data-lesion` hold a whole
      space-joined list and still match a single group.

    Every `<option>` value comes from a vocabulary member -- an enum in
    `vocab.py`, not a curated or mirrored string -- so nothing reaches the
    `selects` markup below unescaped. Everything else on the page goes through
    `render.py`, which escapes cells and attributes itself.
    """
    rows: list[Row] = []
    for gene in sorted(facts):
        fact = facts[gene]
        symbol = symbols.get(gene, gene)
        confidence = fact.headline_confidence.value if fact.headline_confidence else ""
        groups = " ".join(group.value for group in fact.lesion_groups)
        diseases = _graded_diseases(validity.get(gene, uncurated()), fact.headline_confidence)
        rows.append(
            Row(
                cells=(
                    Link(text=gene, href=f"../{gene_page_path(HgncId(gene))}"),
                    # The symbol links too, and it is the one a reader reaches
                    # for: `TBX5` is the name a geneticist recognises and
                    # `HGNC:11604` is the identifier they cite. Linking only the
                    # id put the affordance on the string nobody scans for.
                    Link(text=symbol, href=f"../{gene_page_path(HgncId(gene))}"),
                    confidence or _EM_DASH,
                    "; ".join(diseases) or _EM_DASH,
                    # Position matters and was wrong for one commit: the header
                    # was inserted after `definitive for` while the cell stayed
                    # after `atlas curation`, so the column was headed for a
                    # field two places away. `data_table` zips headers to cells
                    # by position and cannot detect that -- see
                    # `test_the_browse_headers_and_cells_line_up`.
                    Markup(_dot_strip(_concordance_for(concordance, gene))),
                    fact.validity_state.value,
                    fact.atlas_curation.value,
                    str(burden_counts.get(gene, 0)) if burden_counts.get(gene) else _EM_DASH,
                    groups or _EM_DASH,
                ),
                attributes=(
                    ("search", f"{gene} {symbol}".lower()),
                    ("lesion", groups),
                    ("confidence", confidence),
                    ("validity", fact.validity_state.value),
                    ("curation", fact.atlas_curation.value),
                    ("burden", "yes" if burden_counts.get(gene) else "no"),
                ),
            )
        )

    # Sorted for the reason every other set-derived sequence in this build is:
    # a `set` of strings iterates in an order that varies with PYTHONHASHSEED,
    # and the option order is part of the page's bytes and therefore its
    # checksum.
    # `aria-label` on every control, because the only other thing naming any of
    # them is its own first option ("any lesion") or, for the search box, a
    # placeholder. Neither is a label: a placeholder is announced inconsistently
    # and disappears the moment a character is typed, and a `<select>`'s first
    # option describes the *current value*, not what the control filters. A
    # reader on a screen reader would meet five unlabelled controls. Raised by
    # review on #14; the labels are literals here rather than derived from
    # `label` so that "any lesion" stays the neutral option text while the
    # control is announced as what it does.
    # `html.escape` on the option value, though every value is a vocabulary
    # member. Review on #15 measured this as the one interpolation on the whole
    # site where escaping was decided by provenance rather than by `render.py`:
    # driving this loop with a value of `"><script>alert(1)</script>` emits a
    # live `<script>` element into the page while the same build's table cells
    # stay correctly escaped. The provenance argument holds today -- all 21
    # members of `LesionGroup`, `Classification`, `ValidityState` and
    # `AtlasCuration` are `[a-z_]` -- and an invariant that rests on nobody ever
    # adding a member with a quote in it is not one this project keeps.
    selects = "".join(
        f'<select name="{name}" aria-label="Filter by {aria}">'
        f'<option value="">{label}</option>'
        + "".join(
            f'<option value="{html.escape(value)}">{html.escape(value)}</option>'
            for value in sorted(values)
        )
        + "</select>"
        for name, label, aria, values in (
            (
                "lesion",
                "any lesion",
                "lesion group",
                {g.value for f in facts.values() for g in f.lesion_groups},
            ),
            (
                "confidence",
                "any confidence",
                "mirrored confidence",
                {f.headline_confidence.value for f in facts.values() if f.headline_confidence},
            ),
            (
                "validity",
                "any validity",
                "validity state",
                {f.validity_state.value for f in facts.values()},
            ),
            (
                "curation",
                "any curation",
                "atlas curation",
                {f.atlas_curation.value for f in facts.values()},
            ),
            # A boolean facet rather than a numeric one: `FILTER_SCRIPT` matches
            # with `.split(' ').indexOf(want)`, which compares strings, so a
            # range filter would need the script to change as well.
            (
                "burden",
                "any burden",
                "whether burden evidence exists",
                {"yes" if burden_counts.get(gene) else "no" for gene in facts},
            ),
        )
    )

    body = (
        "<h1>Genes</h1>"
        + _SCOPE_RULE
        + "<p>Every confidence below is an upstream panel's or submitter's, republished "
        "with its provenance intact &mdash; the atlas authors no validity classification "
        "of its own. The <strong>atlas curation</strong> column is what says whether the "
        "atlas has curated a gene itself, and the <strong>graded for</strong> column "
        "names the disease each panel actually graded. A blank confidence means no "
        "ClinGen expert panel has graded the gene for any of the diseases in scope here; "
        "its page names the submitters whose agreement admitted it.</p>"
        + grade_legend()
        + f'<p>Showing <span id="shown">{len(rows)}</span> of {len(rows)} genes.</p>'
        '<form id="filters" class="filters">'
        '<input name="q" type="search" aria-label="Search by gene symbol or HGNC id" '
        'placeholder="symbol or HGNC id">'
        f"{selects}</form>"
        f"{data_table(_BROWSE_HEADERS, rows, table_id='gene-table')}" + _STRIP_LEGEND
    )
    emitter.write_text(
        GENE_INDEX_PAGE,
        document(title="Genes — CHD Atlas", root="../", body=body, script=FILTER_SCRIPT),
    )
