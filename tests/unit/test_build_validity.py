# tests/unit/test_build_validity.py
"""The three published validity states, and the discordance that travels with them."""

from __future__ import annotations

import polars as pl
import pytest

from chd_atlas.build.validity import GeneValidity, gene_validity, uncurated
from chd_atlas.vocab import ValiditySource, ValidityState

_CLINGEN_COLUMNS = (
    "gene",
    "gene_symbol",
    "disease",
    "disease_label",
    "moi",
    "sop",
    "classification",
    "classification_date",
    "gcep",
    "report_url",
)
_GENCC_COLUMNS = (
    "sgc_id",
    "gene",
    "gene_symbol",
    "disease",
    "disease_label",
    "moi",
    "classification",
    "submitter",
    "submitted_on",
    "report_url",
)

# Real HGNC/MONDO ids and labels, measured 2026-08-03 against the committed
# mirrors -- not synthesised, so a fixture cannot disagree with the mirror it is
# meant to stand in for.
_TBX5 = "HGNC:11604"
_GATA6 = "HGNC:4174"
_GATA6_CHD = "MONDO:0100540"
_GATA6_DCM = "MONDO:0005021"
_KIF26A = "HGNC:20226"


def _clingen_row(**overrides: object) -> dict[str, object]:
    row: dict[str, object] = {
        "gene": _GATA6,
        "gene_symbol": "GATA6",
        "disease": _GATA6_CHD,
        "disease_label": "GATA6-related congenital heart disease",
        "moi": "AD",
        "sop": "SOP9",
        "classification": "Definitive",
        "classification_date": "2023-11-21",
        "gcep": "Congenital Heart Disease",
        "report_url": "https://search.clinicalgenome.org/a",
    }
    row.update(overrides)
    return row


def _gencc_row(**overrides: object) -> dict[str, object]:
    row: dict[str, object] = {
        "sgc_id": "SGC-000001",
        "gene": _GATA6,
        "gene_symbol": "GATA6",
        "disease": _GATA6_CHD,
        "disease_label": "GATA6-related congenital heart disease",
        "moi": "Autosomal dominant",
        "classification": "Definitive",
        "submitter": "Ambry Genetics",
        "submitted_on": "2020-01-01",
        "report_url": None,
    }
    row.update(overrides)
    return row


def _clingen_frame(*rows: dict[str, object]) -> pl.DataFrame:
    return pl.DataFrame(list(rows) or None, schema={c: pl.String for c in _CLINGEN_COLUMNS})


def _gencc_frame(*rows: dict[str, object]) -> pl.DataFrame:
    return pl.DataFrame(list(rows) or None, schema={c: pl.String for c in _GENCC_COLUMNS})


def test_a_clingen_record_marks_expert_curated_and_gencc_records_still_publish() -> None:
    """GenCC rows for a gene ClinGen has also curated are not hidden.

    The point of publishing both mirrors is to show what the whole field says,
    not to hide agreement behind whichever source is "more authoritative".
    """
    clingen = _clingen_frame(_clingen_row(gene=_TBX5, disease="MONDO:0007732", moi="AD"))
    gencc = _gencc_frame(_gencc_row(gene=_TBX5, disease="MONDO:0007732", moi="Autosomal dominant"))

    result = gene_validity(clingen, gencc, in_scope={"MONDO:0007732"})

    tbx5 = result[_TBX5]
    assert tbx5.state is ValidityState.EXPERT_CURATED
    assert {record.source for record in tbx5.records} == {
        ValiditySource.CLINGEN,
        ValiditySource.GENCC,
    }


def test_a_gene_with_only_gencc_records_is_submitter_curated() -> None:
    result = gene_validity(_clingen_frame(), _gencc_frame(_gencc_row()), in_scope={_GATA6_CHD})

    assert result[_GATA6].state is ValidityState.SUBMITTER_CURATED


def test_out_of_scope_disease_rows_are_dropped_entirely() -> None:
    """GATA6 is Definitive for GATA6-related CHD and Limited for dilated cardiomyopathy.

    Both current, curated by different ClinGen expert panels (measured
    2026-08-03 against the committed mirror). Only the first is CHD; admitting
    the second would make the gene read as contested *about CHD* when the two
    verdicts concern different diseases.
    """
    clingen = _clingen_frame(
        _clingen_row(disease=_GATA6_CHD, classification="Definitive", gcep="CHD"),
        _clingen_row(
            disease=_GATA6_DCM,
            disease_label="dilated cardiomyopathy",
            classification="Limited",
            gcep="Dilated Cardiomyopathy",
        ),
    )

    result = gene_validity(clingen, _gencc_frame(), in_scope={_GATA6_CHD})

    gata6 = result[_GATA6]
    assert [record.disease for record in gata6.records] == [_GATA6_CHD]
    assert gata6.has_source_discordance is False


def test_records_and_genes_come_back_sorted_by_content_not_encounter_order() -> None:
    """Sorted against literals, not by building twice and comparing.

    `PYTHONHASHSEED` is fixed for the life of one interpreter, so two calls in
    the same process would iterate any unsorted dict or set identically either
    way -- a same-process comparison can never catch a missing sort (measured
    directly for this exact failure mode elsewhere in this project: dropping a
    sort and building under seeds 0 and 12345 gave byte-identical output over
    the one-gene committed corpus). Pinning the exact published order is the
    only check that can fail when the sort is removed.

    Both mirrors are handed rows out of the order this test expects back:
    ClinGen's two rows are inserted GATA6-CHD-then-dilated-cardiomyopathy, which
    sorts the other way round by `disease`; GenCC's two rows are inserted
    Zeta-then-Alpha, which sorts the other way round by `submitter`. TBX5 is
    included too, so the returned dict's own key order -- HGNC:11604 before
    HGNC:4174 -- is pinned in the same assertion.
    """
    clingen = _clingen_frame(
        _clingen_row(disease=_GATA6_CHD, moi="AD", classification="Definitive"),
        _clingen_row(
            disease=_GATA6_DCM,
            disease_label="dilated cardiomyopathy",
            moi="AD",
            classification="Limited",
        ),
    )
    gencc = _gencc_frame(
        _gencc_row(sgc_id="SGC-1", submitter="Zeta Labs", classification="Strong"),
        _gencc_row(sgc_id="SGC-2", submitter="Alpha Diagnostics", classification="Moderate"),
        _gencc_row(sgc_id="SGC-3", gene=_TBX5, disease="MONDO:0007732", submitter="Ambry"),
    )

    result = gene_validity(clingen, gencc, in_scope={_GATA6_CHD, _GATA6_DCM, "MONDO:0007732"})

    assert list(result) == [_TBX5, _GATA6]
    assert [(r.source.value, r.disease, r.moi, r.submitter) for r in result[_GATA6].records] == [
        ("clingen", _GATA6_DCM, "AD", None),
        ("clingen", _GATA6_CHD, "AD", None),
        ("gencc", _GATA6_CHD, "Autosomal dominant", "Alpha Diagnostics"),
        ("gencc", _GATA6_CHD, "Autosomal dominant", "Zeta Labs"),
    ]


def test_supportive_publishes_the_verbatim_term_with_no_mapped_classification() -> None:
    """GenCC's `Supportive` is a mapping exception, not a rung on the ladder.

    5,274 rows in the committed GenCC mirror carry it, all from Orphanet
    (measured 2026-08-03) -- a submitter that asserts an association without
    grading its strength. Coercing it onto `Classification` would invent a
    confidence nobody stated, so the record must carry the verbatim term with
    `classification` left `None`.
    """
    result = gene_validity(
        _clingen_frame(),
        _gencc_frame(_gencc_row(classification="Supportive", submitter="Orphanet")),
        in_scope={_GATA6_CHD},
    )

    record = result[_GATA6].records[0]
    assert record.classification is None
    assert record.classification_term == "Supportive"


def test_an_unmapped_classification_raises_naming_the_term_and_the_gene() -> None:
    """A term either authority adds in a later release must reach a human.

    `raise`, never `assert`: `-O` strips `assert`, and this project's standing
    rule is that a bypassed gate must fail rather than publish a confidence
    nobody asserted.
    """
    clingen = _clingen_frame(_clingen_row(classification="Probably Fine"))

    with pytest.raises(ValueError, match="Probably Fine") as excinfo:
        gene_validity(clingen, _gencc_frame(), in_scope={_GATA6_CHD})
    assert _GATA6 in str(excinfo.value)


def test_discordance_fires_when_one_source_contests_and_another_supports() -> None:
    clingen = _clingen_frame(_clingen_row(classification="Disputed"))
    gencc = _gencc_frame(_gencc_row(classification="Strong", submitter="Ambry Genetics"))

    result = gene_validity(clingen, gencc, in_scope={_GATA6_CHD})

    assert result[_GATA6].has_source_discordance is True


def test_discordance_does_not_fire_for_differing_supportive_grades_across_sources() -> None:
    """Definitive and Strong differ in rank but agree the gene is supported.

    Worked example, measured 2026-08-03: ClinGen is Definitive for GATA6 /
    GATA6-related CHD, and PanelApp Australia (via GenCC) is Strong for the
    same gene-disease pair. `strongest()` already orders Definitive above
    Strong; this flag is not about rank, only about whether one source contests
    what another supports, and here neither does.
    """
    clingen = _clingen_frame(_clingen_row(classification="Definitive"))
    gencc = _gencc_frame(_gencc_row(classification="Strong", submitter="PanelApp Australia"))

    result = gene_validity(clingen, gencc, in_scope={_GATA6_CHD})

    assert result[_GATA6].has_source_discordance is False


def test_discordance_does_not_fire_when_supportive_is_the_only_record_on_one_side() -> None:
    """A `Supportive` record takes no side, so it cannot supply the "other" side.

    ClinGen contests the gene here; GenCC's only record is `Supportive`, which
    maps to no `Classification` at all. A submitter that declined to grade
    evidence has not thereby voted "supported", so this must not read as
    disagreement between the two sources -- only one of them has taken a
    position.
    """
    clingen = _clingen_frame(_clingen_row(classification="Disputed"))
    gencc = _gencc_frame(_gencc_row(classification="Supportive", submitter="Orphanet"))

    result = gene_validity(clingen, gencc, in_scope={_GATA6_CHD})

    assert result[_GATA6].has_source_discordance is False


def test_kif26a_multiple_gencc_submitters_disagreeing_on_grade_is_not_discordance() -> None:
    """Same-source multi-submitter disagreement on strength is not cross-source discordance.

    Worked example, measured 2026-08-03: KIF26A / MONDO:0859332 carries Ambry
    Genetics at both Limited and Strong, Baylor at Moderate, and Labcorp and
    PanelApp Australia at Strong -- five GenCC submissions, zero ClinGen
    records, and every classification supportive-level. `has_source_discordance`
    requires at least two distinct *sources*; GenCC disagreeing with itself
    five ways over is still one source.
    """
    disease = "MONDO:0859332"
    gencc = _gencc_frame(
        _gencc_row(
            sgc_id="SGC-102815",
            gene=_KIF26A,
            gene_symbol="KIF26A",
            disease=disease,
            disease_label="cortical dysplasia, complex, with other brain malformations 11",
            moi="Autosomal recessive",
            classification="Limited",
            submitter="Ambry Genetics",
        ),
        _gencc_row(
            sgc_id="SGC-104042",
            gene=_KIF26A,
            gene_symbol="KIF26A",
            disease=disease,
            disease_label="cortical dysplasia, complex, with other brain malformations 11",
            moi="Autosomal recessive",
            classification="Strong",
            submitter="Ambry Genetics",
        ),
        _gencc_row(
            sgc_id="SGC-125207",
            gene=_KIF26A,
            gene_symbol="KIF26A",
            disease=disease,
            disease_label="cortical dysplasia, complex, with other brain malformations 11",
            moi="Autosomal recessive",
            classification="Moderate",
            submitter="Baylor College of Medicine Research Center",
        ),
    )

    result = gene_validity(_clingen_frame(), gencc, in_scope={disease})

    kif26a = result[_KIF26A]
    assert kif26a.state is ValidityState.SUBMITTER_CURATED
    assert kif26a.has_source_discordance is False


def test_uncurated_is_the_shape_for_no_authority_and_a_fresh_instance_each_call() -> None:
    """A function, not a shared constant, so no caller can mutate a shared instance."""
    first = uncurated()
    second = uncurated()

    assert first == GeneValidity(
        state=ValidityState.UNCURATED, records=(), has_source_discordance=False
    )
    assert first is not second
