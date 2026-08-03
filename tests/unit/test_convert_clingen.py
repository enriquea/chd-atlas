# tests/unit/test_convert_clingen.py
"""The ClinGen CSV preamble is a trap, and the trap is silent.

Lines 1-3 are a title block, line 4 and line 6 are `+++` separator rows, and
line 5 is the real header. A `read_csv` with default arguments takes line 1 as
the header and produces ten columns named "CLINGEN GENE DISEASE VALIDITY
CURATIONS", "", "", ... -- no exception, no empty frame, just a table whose every
column name is wrong. Every downstream check then reports a missing column and
nothing names the cause.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from convert_clingen import convert  # noqa: E402

BANNER = (
    '"CLINGEN GENE DISEASE VALIDITY CURATIONS","","","","","","","","",""\n'
    '"FILE CREATED: 2026-08-03","","","","","","","","",""\n'
    '"WEBPAGE: https://search.clinicalgenome.org/kb/gene-validity","","","","","","","","",""\n'
    '"+++++++++++","++++++++++++++","+++++++++++++","++++++++++++++++++","+++++++++",'
    '"+++++++++","++++++++++++++","+++++++++++++","+++++++++++++++++++","+++++++++++++++++++"\n'
    '"GENE SYMBOL","GENE ID (HGNC)","DISEASE LABEL","DISEASE ID (MONDO)","MOI","SOP",'
    '"CLASSIFICATION","ONLINE REPORT","CLASSIFICATION DATE","GCEP"\n'
    '"+++++++++++","++++++++++++++","+++++++++++++","++++++++++++++++++","+++++++++",'
    '"+++++++++","++++++++++++++","+++++++++++++","+++++++++++++++++++","+++++++++++++++++++"\n'
)

ROWS = (
    '"TBX5","HGNC:11604","Holt-Oram syndrome","MONDO:0007732","AD","SOP11","Definitive",'
    '"https://example.org/a","2025-03-25T16:00:00.000Z","Syndromic Disorders GCEP"\n'
    '"GATA4","HGNC:4173","congenital heart disease","MONDO:0005453","AD","SOP9","Limited",'
    '"https://example.org/b","2020-01-01T16:00:00.000Z","Congenital Heart Disease GCEP"\n'
)


def test_the_banner_is_skipped_and_separator_rows_are_dropped(tmp_path: Path) -> None:
    source = tmp_path / "clingen.csv"
    source.write_text(BANNER + ROWS, encoding="utf-8")
    out = tmp_path / "clingen_gene_validity.tsv"

    convert(source, out)

    lines = out.read_text(encoding="utf-8").splitlines()
    assert lines[0].split("\t") == [
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
    ]
    assert len(lines) == 3, "two data rows, no separator rows"
    assert not any("+++" in line for line in lines)


def test_output_is_sorted_by_the_triple(tmp_path: Path) -> None:
    """Sorted on write, so `validate_sort_order` has something to verify.

    The source CSV arrives sorted by gene symbol; the mirror's declared sort key
    is (gene, disease, moi) on the HGNC id. TBX5 is HGNC:11604 and GATA4 is
    HGNC:4173, so symbol order and id order disagree here on purpose -- string
    comparison puts "HGNC:11604" before "HGNC:4173".
    """
    source = tmp_path / "clingen.csv"
    source.write_text(BANNER + ROWS, encoding="utf-8")
    out = tmp_path / "out.tsv"

    convert(source, out)

    genes = [line.split("\t")[0] for line in out.read_text().splitlines()[1:]]
    assert genes == sorted(genes)


def test_a_renamed_upstream_column_stops_the_conversion(tmp_path: Path) -> None:
    """Loud, not partial.

    A renamed upstream column must stop the conversion rather than produce a
    mirror missing a field every downstream check assumes is there.
    """
    import pytest

    source = tmp_path / "clingen.csv"
    source.write_text(BANNER.replace('"MOI"', '"MODE OF INHERITANCE"') + ROWS, encoding="utf-8")

    with pytest.raises(ValueError, match="MOI"):
        convert(source, tmp_path / "out.tsv")
