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

Every value interpolated into a page goes through `render.py`, which escapes
text, links and attributes -- except the three bare interpolations in `_rail`,
which escape themselves and say so. Nothing curated or mirrored is reviewed the
way page copy is, and this is the one artifact kind where a `<` in a disease
label is the difference between a string and a script tag.

Nothing here reads a clock and every value comes from the caller, so two builds
of one commit render byte-identical pages. `build_gene_pages` iterates its
mapping in sorted order for the reason every other loop in this build does.
"""

from __future__ import annotations

import html
from collections.abc import Mapping, Sequence
from typing import Final

from chd_atlas.build.derive import GeneFacts
from chd_atlas.build.emit import Emitter
from chd_atlas.build.paths import GENE_INDEX_PAGE, gene_bundle_path, gene_page_path
from chd_atlas.build.render import (
    FILTER_SCRIPT,
    Cell,
    Link,
    Row,
    chip,
    data_table,
    definition_list,
    document,
)
from chd_atlas.build.validity import GeneValidity, uncurated
from chd_atlas.identifiers import HgncId
from chd_atlas.models.assertion import LesionAssertion
from chd_atlas.models.literature import Publication
from chd_atlas.vocab import AtlasCuration, Classification

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

_CITATION_HEADERS: Final = ("id", "title", "year")

_BROWSE_HEADERS: Final = (
    "gene",
    "symbol",
    "confidence",
    "validity",
    "atlas curation",
    "lesion groups",
)


def _pubmed(pmid: str) -> str:
    """The canonical article URL for a `PMID:` identifier.

    A literal rather than a mirrored field: `Publication` records a DOI and a
    PMCID but no PubMed URL, and this is the one link on the page the atlas
    constructs rather than republishes.
    """
    return f"https://pubmed.ncbi.nlm.nih.gov/{pmid.removeprefix('PMID:')}/"


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


def _rail(gene: str, symbol: str, fact: GeneFacts) -> str:
    """The summary column, populated for every gene whether curated here or not.

    `chip-definitive` is the only chip class the stylesheet fills with the
    success colour, so it is keyed on the classification actually being
    `DEFINITIVE` rather than on there being one at all: `strongest()` resolves a
    gene carrying both a definitive and a refuted record to definitive, and a
    green pill beside the word "refuted" is exactly the "contested gene shown as
    settled" failure the ladder is documented against.
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
        ("validity", fact.validity_state.value),
        ("atlas curation", fact.atlas_curation.value),
        ("lesion groups", ", ".join(group.value for group in fact.lesion_groups) or _EM_DASH),
        ("assertions", str(fact.assertion_count)),
        ("functional records", str(fact.functional_count)),
        ("publications", str(len(fact.publications))),
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


def build_gene_pages(
    facts: Mapping[str, GeneFacts],
    emitter: Emitter,
    *,
    symbols: Mapping[str, str],
    validity: Mapping[str, GeneValidity],
    assertions: Mapping[str, list[LesionAssertion]],
    publications: Mapping[str, Publication],
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
        body = (
            '<div class="layout">'
            + _rail(gene, symbol, fact)
            + "<div>"
            + _validity_section(validity.get(gene, uncurated()))
            + (
                _evidence_section(records, publications)
                if fact.atlas_curation is AtlasCuration.CURATED
                else _not_curated(fact)
            )
            + "</div></div>"
        )
        emitter.write_text(
            gene_page_path(HgncId(gene)),
            document(title=f"{symbol} — CHD Atlas", root="../", body=body),
        )


def build_gene_index_page(
    facts: Mapping[str, GeneFacts], emitter: Emitter, *, symbols: Mapping[str, str]
) -> None:
    """Emit `genes/index.html`: every published gene, filterable in the browser.

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
        rows.append(
            Row(
                cells=(
                    Link(text=gene, href=f"../{gene_page_path(HgncId(gene))}"),
                    symbol,
                    confidence or _EM_DASH,
                    fact.validity_state.value,
                    fact.atlas_curation.value,
                    groups or _EM_DASH,
                ),
                attributes=(
                    ("search", f"{gene} {symbol}".lower()),
                    ("lesion", groups),
                    ("confidence", confidence),
                    ("validity", fact.validity_state.value),
                    ("curation", fact.atlas_curation.value),
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
    selects = "".join(
        f'<select name="{name}" aria-label="Filter by {aria}">'
        f'<option value="">{label}</option>'
        + "".join(f'<option value="{value}">{value}</option>' for value in sorted(values))
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
        )
    )

    body = (
        "<h1>Genes</h1>"
        "<p>Every confidence below is an upstream panel's or submitter's, republished "
        "with its provenance intact &mdash; the atlas authors no validity classification "
        "of its own. The <strong>atlas curation</strong> column is what says whether the "
        "atlas has curated a gene itself.</p>"
        f'<p>Showing <span id="shown">{len(rows)}</span> of {len(rows)} genes.</p>'
        '<form id="filters" class="filters">'
        '<input name="q" type="search" aria-label="Search by gene symbol or HGNC id" '
        'placeholder="symbol or HGNC id">'
        f"{selects}</form>"
        f"{data_table(_BROWSE_HEADERS, rows, table_id='gene-table')}"
    )
    emitter.write_text(
        GENE_INDEX_PAGE,
        document(title="Genes — CHD Atlas", root="../", body=body, script=FILTER_SCRIPT),
    )
