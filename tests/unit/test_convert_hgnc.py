# tests/unit/test_convert_hgnc.py
"""What the HGNC converter must not lose on the way into the mirror.

Every column here is one a downstream join reads. `uniprot` in particular is
what any future PTM or proteomics join keys on, so a gene silently arriving
with a null accession is evidence loss rather than a cosmetic gap.
"""

from pathlib import Path

import polars as pl
from scripts.convert_hgnc import convert

_HEADER = (
    "hgnc_id\tsymbol\tname\talias_symbol\tensembl_gene_id\tentrez_id\t"
    "location\tuniprot_ids\tmane_select\tstatus"
)


def _source(tmp_path: Path, *rows: str) -> Path:
    path = tmp_path / "hgnc.tsv"
    path.write_text("\n".join((_HEADER, *rows)) + "\n", encoding="utf-8")
    return path


def test_the_refseq_half_of_mane_select_is_kept_and_the_ensembl_half_dropped(tmp_path):
    source = _source(
        tmp_path,
        "HGNC:11604\tTBX5\tT-box transcription factor 5\t\tENSG00000089225\t6910\t"
        '12q24.21\tQ99593\t"ENST00000405440.7|NM_181486.4"\tApproved',
    )
    out = tmp_path / "genes.tsv"
    convert(source, out, keep={"HGNC:11604"})

    row = pl.read_csv(out, separator="\t", infer_schema_length=0).row(0, named=True)
    assert row["mane_select"] == "NM_181486.4"
    assert row["uniprot"] == "Q99593"
    assert row["ncbi_gene"] == "6910"
    assert row["locus"] == "12q24.21"
    assert row["aliases"] is None


def test_the_first_uniprot_accession_is_kept_when_hgnc_lists_several(tmp_path):
    """RBM10 is the live instance: HGNC lists `P0DW28|P98175`.

    The mirror's `uniprot` column carries `UNIPROT_PATTERN`, which rejects a
    pipe, so the pair cannot be published as one cell. Publishing null instead
    would drop RBM10 out of every accession-keyed join silently — the failure
    this project is built to avoid — so the first accession is taken and the
    choice is pinned here, where changing it is visible.
    """
    source = _source(
        tmp_path,
        "HGNC:9896\tRBM10\tRNA binding motif protein 10\t\tENSG00000182872\t8241\t"
        'Xp11.3\t"P0DW28|P98175"\t"ENST00000377604.9|NM_005676.5"\tApproved',
    )
    out = tmp_path / "genes.tsv"
    convert(source, out, keep={"HGNC:9896"})

    assert pl.read_csv(out, separator="\t", infer_schema_length=0)["uniprot"][0] == "P0DW28"


def test_genes_outside_the_keep_set_are_not_mirrored(tmp_path):
    source = _source(
        tmp_path,
        "HGNC:11604\tTBX5\tT-box transcription factor 5\t\tENSG00000089225\t6910\t"
        '12q24.21\tQ99593\t"ENST00000405440.7|NM_181486.4"\tApproved',
        "HGNC:5\tA1BG\talpha-1-B glycoprotein\t\tENSG00000121410\t1\t"
        '19q13.43\tP04217\t"ENST00000263100.8|NM_130786.4"\tApproved',
    )
    out = tmp_path / "genes.tsv"
    assert convert(source, out, keep={"HGNC:11604"}) == 1

    frame = pl.read_csv(out, separator="\t", infer_schema_length=0)
    assert frame["hgnc_id"].to_list() == ["HGNC:11604"]


def test_rows_are_sorted_by_hgnc_id_as_a_string(tmp_path):
    """`TABLE_SCHEMAS["genes"].sort_key` is `("hgnc_id",)` and `sort_order.py`
    compares the published strings, so "HGNC:11604" precedes "HGNC:2488".
    Asserted against a literal rather than by sorting the output again, which
    would only prove the test agrees with itself.
    """
    source = _source(
        tmp_path,
        "HGNC:2488\tNKX2-5\tNK2 homeobox 5\t\tENSG00000183072\t1482\t"
        '5q35.1\tP52952\t"ENST00000329198.5|NM_004387.4"\tApproved',
        "HGNC:11604\tTBX5\tT-box transcription factor 5\t\tENSG00000089225\t6910\t"
        '12q24.21\tQ99593\t"ENST00000405440.7|NM_181486.4"\tApproved',
    )
    out = tmp_path / "genes.tsv"
    convert(source, out, keep={"HGNC:2488", "HGNC:11604"})

    frame = pl.read_csv(out, separator="\t", infer_schema_length=0)
    assert frame["hgnc_id"].to_list() == ["HGNC:11604", "HGNC:2488"]
