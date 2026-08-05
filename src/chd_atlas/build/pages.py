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

**The uncurated notice is load-bearing.** 22 of the 23 genes published under D21
carry no `LesionAssertion` (measured 2026-08-04 on `genes/index.json`). Their
pages state that the atlas has not curated them rather than omitting the
evidence section, because an absent section is indistinguishable from "the atlas
looked and found nothing" and a reader deciding what a gene means clinically
must not have to infer which. Its exact wording is load-bearing for a second
reason -- it sits one column from a rail that counts curated content this gene
may have -- and `_not_curated` is where that wording and the reason for it are
recorded.

**Neither page kind states a bare `definitive`.** `_SCOPE_RULE` says what admits
a gene to this site at all, and `_definitive_diseases` names the disease each
ClinGen panel actually graded -- in the browse table's `definitive for` column
and in the gene page's rail, one line under the green chip. Both read the same
`validity` mapping the bundles were written from, so the two pages cannot name
different diseases for one gene.

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
from typing import Final

from chd_atlas.build.burden import BurdenRow, shared_cohorts
from chd_atlas.build.concordance import FamilyState, family_state
from chd_atlas.build.derive import GeneFacts
from chd_atlas.build.emit import Emitter, Json
from chd_atlas.build.paths import GENE_INDEX_PAGE, gene_bundle_path, gene_page_path
from chd_atlas.build.render import (
    FILTER_SCRIPT,
    Cell,
    Link,
    Markup,
    Row,
    chip,
    data_table,
    definition_list,
    document,
)
from chd_atlas.build.validity import GeneValidity, uncurated
from chd_atlas.identifiers import HgncId
from chd_atlas.models.assertion import LesionAssertion
from chd_atlas.models.cohort import Cohort
from chd_atlas.models.literature import Publication
from chd_atlas.vocab import AtlasCuration, Classification, ValiditySource

# What the 23-gene set *is*, in one sentence, on both page kinds that show it.
#
# No page on this site said it. `docs/data-api.md` states the rule and a reader
# of the HTML never sees it, so a browse page of 23 rows all reading
# `definitive`, under a title reading "CHD Atlas" and a tagline about congenital
# heart disease, said "these 23 genes are definitive for CHD". Measured
# 2026-08-04 against the committed mirrors, that is false for most of them: 10
# of the 23 are qualified by a panel that is not the Congenital Heart Disease
# GCEP -- SCID-CID for KMT2D and KDM6A, Hearing Loss for CHD7, Kidney Cystic and
# Ciliopathy Disorders for NOTCH2 and EVC, Intellectual Disability and Autism
# for ZEB2, Skeletal Disorders for EVC2, Syndromic Disorders for TBX5 and
# TFAP2B, Prenatal for FLT4 -- and 10 are definitive for a disease whose label
# names no cardiac feature at all (Kabuki syndrome 1 and 2, Mowat-Wilson,
# CHARGE, Alagille, Ellis-van Creveld twice, Holt-Oram, TARP, and NR2F2-related
# multiple congenital anomalies).
#
# One constant for both pages, in the idiom `_NOT_CURATED` already uses: the
# rule is one editorial claim, and two copies of it are two things that drift.
_SCOPE_RULE: Final = (
    '<p class="scope-rule">A gene is published here when a ClinGen expert panel classifies '
    "it <strong>Definitive</strong> for a disease in this atlas's CHD scope "
    "(<code>curation/chd_scope.yaml</code>). That is not the same as definitive for "
    "congenital heart disease: for many of these genes the disease graded is a syndrome "
    "of which a cardiac lesion is one feature, and the disease the panel actually named "
    "is shown beside the gene.</p>"
)

_NOT_CURATED: Final = (
    "<p>The atlas has <strong>not yet curated</strong> a lesion assertion for this gene. "
    "The classification above is an expert panel's, mirrored with its provenance "
    "intact; no classification on this page is the atlas's own assessment.</p>"
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
_STRIP_LEGEND: Final = (
    '<p class="strip-legend">'
    '<span class="dot full"></span> enriched, and survives that study\'s own correction'
    ' &nbsp; <span class="dot half"></span> enriched nominally, or no correction published'
    ' &nbsp; <span class="dot none"></span> tested, no enrichment detected'
    ' &nbsp; <span class="dot untested">&ndash;</span> not tested by that dataset'
    "</p>"
    '<p class="strip-legend">One glyph per <strong>independent cohort family</strong>, not '
    "per study: two papers sharing a sample collection describe the same people and count "
    "once. <strong>No enrichment at these cohort sizes is not evidence against a gene</strong> "
    "&mdash; burden tests routinely detect nothing for genes with overwhelming family and "
    "functional evidence, so read the tally as how much has been looked at, never as a "
    "verdict on the gene.</p>"
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
_POOLING_NOTICE: Final = (
    "<p>The atlas computes <strong>no pooled statistic across studies</strong>: these "
    "cohorts overlap, so combining them would count the same people twice.</p>"
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

# `definitive for` sits immediately after `confidence`, because it is what
# qualifies it: the two columns are one claim split in two, and a reader whose
# eye stops at `definitive` must meet the disease next rather than three columns
# later.
_BROWSE_HEADERS: Final = (
    "gene",
    "symbol",
    "confidence",
    "definitive for",
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


def _definitive_diseases(gene_validity: GeneValidity) -> tuple[str, ...]:
    """The in-scope diseases a ClinGen panel grades this gene `Definitive` for.

    Exactly the records `validity.published_genes` admits the gene on -- source
    `CLINGEN`, classification `DEFINITIVE` -- so what this returns is the reason
    the gene has a page, not a summary of everything either mirror says about it.
    A GenCC submitter calling the gene definitive is deliberately not here: it
    admits no gene to this set (D21), so naming it would qualify the confidence
    with a disease that is not the one the gate turned on.

    Returns labels rather than MONDO ids: this is the string a reader reads
    beside the word `definitive`, and the id is already on the gene page in
    `_validity_section`, which renders both.

    **Sorted, and de-duplicated through a `set` first.** A gene may hold more
    than one qualifying record -- two panels grading the same gene for two
    in-scope diseases, or one disease reached through two rows -- and a `set` of
    strings iterates in an order that varies with `PYTHONHASHSEED`, which would
    make the browse page's bytes and therefore its checksum differ between two
    builds of one commit. Measured 2026-08-04 against the committed mirrors: all
    23 published genes have exactly one qualifying record and one distinct
    label, so the multi-label case is latent and its ordering is guarded by
    fixture alone (`test_a_gene_definitive_for_two_in_scope_diseases_names_both`).

    An empty return is unreachable for a member of the published population,
    which is derived from the same predicate -- a guard on a bypassed gate, in
    the idiom of `_report_link` above: the caller renders an em dash rather than
    an empty cell that reads as "no disease".
    """
    return tuple(
        sorted(
            {
                record.disease_label
                for record in gene_validity.records
                if record.source is ValiditySource.CLINGEN
                and record.classification is Classification.DEFINITIVE
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
    if fact.functional_count == 0:
        return _NOT_CURATED
    noun = "record" if fact.functional_count == 1 else "records"
    return _NOT_CURATED + (
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
    gene: str, symbol: str, fact: GeneFacts, diseases: Sequence[str], burden_rows: int
) -> str:
    """The summary column, populated for every gene whether curated here or not.

    `chip-definitive` is the only chip class the stylesheet fills with the
    success colour, so it is keyed on the classification actually being
    `DEFINITIVE` rather than on there being one at all: `strongest()` resolves a
    gene carrying both a definitive and a refuted record to definitive, and a
    green pill beside the word "refuted" is exactly the "contested gene shown as
    settled" failure the ladder is documented against.

    `diseases` is `_definitive_diseases`' return, and the `definitive for` row is
    the first in the list for the same reason the browse column sits beside
    `confidence`: the green chip one line above it says `definitive` and nothing
    else, and on this site's title and tagline that reads as "definitive for
    congenital heart disease". For KMT2D it means Kabuki syndrome 1. The row is
    what stops the chip standing alone on the first screen.
    """
    headline = fact.headline_confidence
    chips = [
        chip(
            headline.value if headline is not None else "not classified",
            kind="definitive" if headline is Classification.DEFINITIVE else "",
        )
    ]
    if fact.has_conflicting_evidence:
        chips.append(chip("conflicting evidence", kind="warn"))
    if fact.has_source_discordance:
        chips.append(chip("sources disagree", kind="warn"))

    facts = [
        ("definitive for", "; ".join(diseases) or _EM_DASH),
        ("validity", fact.validity_state.value),
        ("atlas curation", fact.atlas_curation.value),
        ("lesion groups", ", ".join(group.value for group in fact.lesion_groups) or _EM_DASH),
        ("assertions", str(fact.assertion_count)),
        ("functional records", str(fact.functional_count)),
        # "curated publications", not "publications". `fact.publications` is the
        # PMIDs this atlas's own assertion evidence cites, and it reads 0 on 22
        # of the 23 published genes -- which sat directly beside a burden section
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
        # given study on all 23 gene pages. `warning`, the disclosure, the
        # composite note and the footnotes do not -- each changes what a reader
        # concludes from the table right below it.
        counted = (
            f'<details class="reading-notes"><summary>How this study counted'
            f"</summary>{design}{_provenance(study_rows, cohorts)}</details>"
            if design or _provenance(study_rows, cohorts)
            else ""
        )
        blocks.append(
            f"<h3>"
            f'<a href="{html.escape(_pubmed(study))}">'
            f"{html.escape(_study_label(study, publications))}</a></h3>"
            f"{_disclosure(study, publications)}"
            f"{warning}"
            f"{counted}"
            f"{_composite_note(study_rows)}"
            f"{table}"
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
    # synonymous negative control. All three are identical on all 23 pages, and
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
    the 23 published gene pages. Eight of those pages now carry three tables with
    the declaration on two, which is the arrangement that makes the line
    informative: an unconditional disclosure says nothing, and a stale claim that
    it is unconditional invites someone to make it one.
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
    than from this gene, so a design nobody ran renders as absent instead of
    vanishing -- the same confusion between "not tested" and "tested and found
    nothing" that `FamilyState` has four members for. Nobody has published a CNV
    de novo analysis for any gene, and that quadrant is empty on every page: the
    shape of the literature is itself information, and it is the argument for
    which study to curate next.

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
        f"<tbody>{''.join(body)}</tbody></table></div>"
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
    statistic = (
        f"q {best.pvalue_adjusted:.3g}"
        if best.pvalue_adjusted is not None
        else (f"p {best.pvalue:.3g}" if best.pvalue is not None else _EM_DASH)
    )
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
    survivorship bias that inflates every odds ratio on all 23 pages -- and that
    the combined cohort is ~92% European. `models/cohort.py` says in as many
    words that `description` is where such caveats belong. `cohort_labels`
    returned `{id: name}` and dropped them, so the curator wrote the caveat, the
    model documented it, and no reader could reach it.

    In a `<details>` element: the caveats are long, they repeat on 23 pages, and
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

    Sorted, like every loop in this build that iterates a mapping: `sort_keys`
    orders dict keys in a JSON payload and has nothing to say about the order
    files are written in, and a differing write order would move nothing here
    today but is the habit this project keeps.
    """
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
                _definitive_diseases(gene_validity),
                len(burden.get(gene, ())),
            )
            + "<div>"
            + _SCOPE_RULE
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
    `definitive for` column are the other half: a row reading `KMT2D |
    definitive` on a site titled "CHD Atlas" claims KMT2D is definitive for
    congenital heart disease, and ClinGen's actual assertion is `KMT2D --
    Definitive for Kabuki syndrome 1`, made by the SCID-CID GCEP. `validity` is
    taken for that column alone; `facts` carries no disease, and deriving one
    from a second read of the mirrors would let this page and the gene page it
    links to name different diseases for the same gene.

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
        diseases = _definitive_diseases(validity.get(gene, uncurated()))
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
                    fact.validity_state.value,
                    fact.atlas_curation.value,
                    Markup(_dot_strip((concordance or {}).get(gene, {}))),
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
        "atlas has curated a gene itself, and the <strong>definitive for</strong> column "
        "names the disease each panel actually graded.</p>"
        f'<p>Showing <span id="shown">{len(rows)}</span> of {len(rows)} genes.</p>'
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
