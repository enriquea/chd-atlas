# src/chd_atlas/build/validity.py
"""Gene-disease validity, mirrored and attributed -- never authored here.

Design decision D12 (`docs/superpowers/specs/2026-08-03-validity-backbone-design.md`):
the atlas publishes no validity classification of its own. One curator cannot
out-curate a chartered expert panel, and the resource that tried is the
cautionary case: CHDgene publishes a documented 3-tier ranking under which all
189 genes carry 5/5 stars.

What this module does instead is join the two mirrors -- `mirrors/clingen_gene_
validity.tsv` and `mirrors/gencc_submissions.tsv` -- to `curation/chd_scope.yaml`
and say, per gene, which of three things is true: an expert panel has assessed
it (`EXPERT_CURATED`), only submitters have (`SUBMITTER_CURATED`), or nobody has
(`UNCURATED`). The third is published as a fact rather than left as an absent
key, because a consumer cannot tell "no authority has assessed this gene" from
"the build dropped it" by looking at a missing field -- and that distinction is
exactly the silent evidence loss `CLAUDE.md` calls this project's characteristic
failure.

`has_source_discordance` is decided and returned alongside `records` in the same
`GeneValidity`, for the reason `bundles.py::_headline` writes
`headline_confidence` and `has_conflicting_evidence` together: a payload must
never be able to carry the records without the flag that says whether they
disagree.

Every published array here is sorted by the caller. `encode_json` (`build/
emit.py`) passes `sort_keys=True`, which orders dict *keys* and nothing else --
an array built by iterating a `dict` populated one mirror row at a time checksums
differently between processes unless this module sorts it before returning.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

import polars as pl

from chd_atlas.vocab import (
    CLINGEN_CLASSIFICATIONS,
    CONTESTED,
    GENCC_CLASSIFICATIONS,
    Classification,
    ValiditySource,
    ValidityState,
)


@dataclass(frozen=True)
class ValidityRecord:
    """One authority's verdict on one gene-disease pair.

    `classification_term` is what the authority published, verbatim.
    `classification` is that term mapped onto the atlas's `Classification`
    ladder, or `None` where the term is not a rung on it at all -- GenCC's
    `Supportive`, the mapping exception for a submitter that asserts an
    association without grading its strength. Publishing both lets a consumer
    render the authority's own word and still filter on a common scale, without
    the atlas inventing a confidence nobody stated.

    The last four fields are authority-specific: `sop`, `classification_date`
    and `gcep` come from ClinGen only, `submitter` from GenCC only. Each is
    `None` on a record from the other source rather than the two sources
    getting separate record types, because every consumer of `records` wants
    one homogeneous list to sort, filter and render, not a union type to match
    on first.
    """

    source: ValiditySource
    classification_term: str
    classification: Classification | None
    disease: str
    disease_label: str
    moi: str
    report_url: str | None
    sop: str | None = None
    classification_date: str | None = None
    gcep: str | None = None
    submitter: str | None = None


@dataclass(frozen=True)
class GeneValidity:
    """One gene's whole validity picture.

    How curated it is (`state`), by whom and under what verdict (`records`),
    and whether the mirrors agree it is settled (`has_source_discordance`).
    """

    state: ValidityState
    records: tuple[ValidityRecord, ...]
    has_source_discordance: bool


def uncurated() -> GeneValidity:
    """The published shape for a gene no authority has assessed at all.

    Neither an expert panel nor a submitter.

    A function, not a module-level constant: a shared `GeneValidity` instance
    handed to every uncurated gene would still be one object underneath
    "frozen", and nothing stops a future caller mutating through it or relying
    on identity. Returning a fresh instance per call costs nothing and
    forecloses the question.
    """
    return GeneValidity(
        state=ValidityState.UNCURATED,
        records=(),
        has_source_discordance=False,
    )


def _sort_key(record: ValidityRecord) -> tuple[str, str, str, str]:
    """`(source, disease, moi, submitter or "")`.

    So two builds of one commit agree on array order regardless of which
    mirror row was read first. `submitter` is `None` on every ClinGen record;
    `or ""` sorts those first, which never has to arbitrate a real tie because
    a null-vs-populated `submitter` split never happens within one source --
    only GenCC rows carry it, and every GenCC row has one.
    """
    return (record.source.value, record.disease, record.moi, record.submitter or "")


def _clingen_records(frame: pl.DataFrame, in_scope: set[str]) -> dict[str, list[ValidityRecord]]:
    by_gene: dict[str, list[ValidityRecord]] = {}
    for row in frame.to_dicts():
        if row["disease"] not in in_scope:
            continue
        term = row["classification"]
        if term not in CLINGEN_CLASSIFICATIONS:
            raise ValueError(
                f"unmapped ClinGen classification {term!r} for gene {row['gene']}; "
                f"add it to vocab.CLINGEN_CLASSIFICATIONS rather than coercing it "
                f"to an adjacent rung"
            )
        by_gene.setdefault(row["gene"], []).append(
            ValidityRecord(
                source=ValiditySource.CLINGEN,
                classification_term=term,
                classification=CLINGEN_CLASSIFICATIONS[term],
                disease=row["disease"],
                disease_label=row["disease_label"],
                moi=row["moi"],
                report_url=row["report_url"],
                sop=row["sop"],
                classification_date=row["classification_date"],
                gcep=row["gcep"],
            )
        )
    return by_gene


def _gencc_records(frame: pl.DataFrame, in_scope: set[str]) -> dict[str, list[ValidityRecord]]:
    by_gene: dict[str, list[ValidityRecord]] = {}
    for row in frame.to_dicts():
        if row["disease"] not in in_scope:
            continue
        term = row["classification"]
        if term not in GENCC_CLASSIFICATIONS:
            raise ValueError(
                f"unmapped GenCC classification {term!r} for gene {row['gene']}; "
                f"add it to vocab.GENCC_CLASSIFICATIONS rather than coercing it "
                f"to an adjacent rung"
            )
        by_gene.setdefault(row["gene"], []).append(
            ValidityRecord(
                source=ValiditySource.GENCC,
                classification_term=term,
                classification=GENCC_CLASSIFICATIONS[term],
                disease=row["disease"],
                disease_label=row["disease_label"],
                moi=row["moi"],
                report_url=row["report_url"],
                submitter=row["submitter"],
            )
        )
    return by_gene


def _has_source_discordance(records: tuple[ValidityRecord, ...]) -> bool:
    """True iff some source contests the gene while a *different* source supports it.

    Not "the classifications differ": Definitive and Strong differ in rank but
    both sit outside `CONTESTED`, so two sources reporting one each still agree
    the gene is not contested. Not "one source alone carries both a supportive
    and a contesting call" either -- 90 genes in the committed ClinGen mirror do
    exactly that across two diseases or two GCEPs (measured 2026-08-03), and
    that is ClinGen being internally split, not ClinGen disagreeing with
    anyone. (None of those 90 splits survives the CHD-scope filter this
    module applies before a gene is even keyed -- 13 of the 90 have some
    in-scope ClinGen record, but never both sides of the split within scope,
    so this is a definitional justification, not a case seen in published
    data; measured 2026-08-04 against `gene_validity`'s real output, where
    exactly one gene sets `has_conflicting_evidence` -- LEFTY2, HGNC:3122, via
    a *cross-source* split -- and it also sets `has_source_discordance`.) Only
    a record with a mapped `classification` takes a side: GenCC's
    `Supportive` maps to `None` and is a submitter declining to grade evidence,
    not evidence *for* the gene, and `NO_KNOWN_ASSOCIATION` is a stated absence
    of association, not support for one.

    Implemented as the literal reading of "one authority contesting a gene
    another supports": true iff some `c` in the contesting sources and some `s`
    in the supporting sources are different sources. Both sets are subsets of
    `{ValiditySource.CLINGEN, ValiditySource.GENCC}` here, so on every case this
    module's own tests exercise this agrees with the simpler `contested !=
    supporting`; the pairwise form is kept because it also stays correct in the
    one case the two diverge -- both sources internally split identically
    (`contested == supporting == {CLINGEN, GENCC}`), where each source still
    contests something a *different* source supports, and the simpler
    set-inequality check would miss it.
    """
    contesting = {
        record.source
        for record in records
        if record.classification is not None and record.classification in CONTESTED
    }
    supporting = {
        record.source
        for record in records
        if record.classification is not None
        and record.classification not in CONTESTED
        and record.classification is not Classification.NO_KNOWN_ASSOCIATION
    }
    return any(c != s for c in contesting for s in supporting)


def gene_validity(
    clingen: pl.DataFrame,
    gencc: pl.DataFrame,
    in_scope: set[str],
) -> dict[str, GeneValidity]:
    """One `GeneValidity` per gene either mirror curates within CHD scope.

    `in_scope` is the set of MONDO disease terms from `curation/chd_scope.yaml`
    (`Corpus.chd_scope`) -- the atlas's own editorial claim about what counts as
    congenital heart disease (`models/scope.py`). A row naming any other disease
    is dropped here, before a gene is even keyed, rather than filtered
    downstream: GATA6 is Definitive for GATA6-related congenital heart disease
    (`MONDO:0100540`) and Limited for dilated cardiomyopathy (`MONDO:0005021`,
    not CHD), both current and both curated by different ClinGen expert panels
    (measured 2026-08-03 against the committed mirror). Admitting the second row
    would make the gene read as contested *about CHD* when the two verdicts
    concern different diseases entirely.

    A gene reached through any in-scope ClinGen record is `EXPERT_CURATED`, even
    when GenCC has also curated it -- an expert panel's own verdict does not
    change because submitters have opinions too, and this module publishes
    every in-scope GenCC record regardless, so nothing about the field's wider
    view is hidden by the state a gene is filed under. A gene with only in-scope
    GenCC records is `SUBMITTER_CURATED`. A gene reached through neither mirror
    at all is not a key in the returned dict -- callers wanting the third state
    use `uncurated()` for a gene they otherwise know about (e.g. every gene in
    `mirrors/genes.tsv`), the same way `derive.gene_facts` omits a gene with no
    assertion rather than emitting an empty one.

    An unmapped classification term raises `ValueError` naming the term and the
    gene, rather than falling back to any adjacent rung: a term either mirror
    adds in a later release must reach a human, per `CLAUDE.md`'s standing rule
    that a bypassed gate must fail rather than publish a confidence nobody
    asserted.
    """
    from_clingen = _clingen_records(clingen, in_scope)
    from_gencc = _gencc_records(gencc, in_scope)

    result: dict[str, GeneValidity] = {}
    # Sorted so the dict's own key order -- and so anything a caller builds by
    # iterating it into a JSON array or object -- does not depend on which
    # mirror happened to name a gene first.
    for gene in sorted(set(from_clingen) | set(from_gencc)):
        records = tuple(
            sorted(from_clingen.get(gene, []) + from_gencc.get(gene, []), key=_sort_key)
        )
        result[gene] = GeneValidity(
            state=(
                ValidityState.EXPERT_CURATED
                if gene in from_clingen
                else ValidityState.SUBMITTER_CURATED
            ),
            records=records,
            has_source_discordance=_has_source_discordance(records),
        )
    return result


def published_genes(validity: Mapping[str, GeneValidity]) -> set[str]:
    """The genes the atlas publishes a page for. Design decision D21.

    A gene qualifies when some mirrored record carries `source == CLINGEN` and
    `classification == DEFINITIVE`. Every record reaching here already names a
    disease listed in `curation/chd_scope.yaml` -- `gene_validity` filters on
    that -- so this does not re-check scope.

    **Not** "headline confidence is definitive and validity state is expert
    curated." Measured on the mirrors as committed (2026-08-04), the two rules
    select the identical 23 genes, which is precisely why the difference has to
    be written down. `state` records only that ClinGen has *a* row for the gene,
    never what that row says, so the second rule admits a gene ClinGen graded
    `Limited` on the strength of a GenCC submitter grading it `Definitive`. No
    such gene is in scope today. `test_build_validity.py` constructs one.

    GenCC is excluded because it aggregates rather than adjudicates and says so
    itself. The five in-scope genes it alone calls definitive -- ELN, GDF1,
    MMP21, PKD1L1, TBX1 -- have no ClinGen curation for any disease at all
    except ELN, curated only for cutis laxa, so this is not a gap in the scope
    file that a curator could close. GDF1 is the case that settles it: its
    in-scope submissions run from G2P's `Definitive` to Illumina's `No Known
    Disease Relationship`, and `has_conflicting_evidence` reports `False` for
    that pair because `no_known_association` belongs to neither side of its test
    (spec D34). A 28-gene gate would publish it as settled.
    """
    return {
        gene
        for gene, entry in validity.items()
        if any(
            record.source is ValiditySource.CLINGEN
            and record.classification is Classification.DEFINITIVE
            for record in entry.records
        )
    }
