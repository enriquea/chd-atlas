# tests/unit/test_build_omics.py
import json
from pathlib import Path

import pytest

from chd_atlas.build.emit import Emitter
from chd_atlas.build.omics import TOP_N, _accession_index, build_omics

EXPRESSION_HEADER = (
    "dataset\tcontrast\tgene\tlog2fc\tpvalue\tfdr\tdirection\tn_case\tn_control\ttissue\tstage\n"
)
PROFILES_HEADER = "dataset\tgene\ttissue\tstage\tmedian_abundance\tunit\tq25\tq75\tn_samples\n"
PROTEOMICS_HEADER = (
    "dataset\tcontrast\tprotein\tgene\tlog2fc\tpvalue\tfdr\tdirection\t"
    "n_peptides\tsequence_coverage\tquant_method\n"
)
PHOSPHO_HEADER = (
    "dataset\tcontrast\tsite_id\tprotein\tresidue\tposition\tmod_type\t"
    "flanking_sequence\tlog2fc\tpvalue\tfdr\tprotein_normalized\toccupancy\t"
    "known_kinases\n"
)
GENES_HEADER = (
    "hgnc_id\tsymbol\tname\taliases\tensembl_gene\tncbi_gene\tlocus\tuniprot\tmane_select\n"
)


def _row(
    gene: str,
    fdr: str,
    log2fc: str = "1.0",
    dataset: str = "PXD012345",
    contrast: str = "tof_vs_control",
) -> str:
    return f"{dataset}\t{contrast}\t{gene}\t{log2fc}\t0.001\t{fdr}\tup\t10\t10\tRV\tinfant\n"


def _phospho_row(protein: str, position: int, fdr: str = "0.01") -> str:
    return (
        f"PXD012345\ttof_vs_control\t{protein}_S{position}\t{protein}\tS\t{position}\t"
        f"MOD:00046\t\t1.2\t0.001\t{fdr}\ttrue\t\t\n"
    )


def _repo(tmp_path: Path, rows: str) -> Path:
    (tmp_path / "mirrors" / "expression").mkdir(parents=True)
    (tmp_path / "mirrors" / "expression" / "PXD012345.tsv").write_text(EXPRESSION_HEADER + rows)
    return tmp_path


def _registry(tmp_path: Path, rows: str) -> None:
    """Write `mirrors/genes.tsv`, the only thing that maps an accession to a gene."""
    (tmp_path / "mirrors").mkdir(exist_ok=True)
    (tmp_path / "mirrors" / "genes.tsv").write_text(GENES_HEADER + rows)


def _gene_row(hgnc: str, symbol: str, uniprot: str) -> str:
    return f"{hgnc}\t{symbol}\ta name\t\t\t\t\t{uniprot}\t\n"


def _table(tmp_path: Path, table: str, filename: str, text: str) -> None:
    (tmp_path / "mirrors" / table).mkdir(parents=True, exist_ok=True)
    (tmp_path / "mirrors" / table / filename).write_text(text)


def test_a_shard_holds_the_rows_the_summary_only_counts(tmp_path: Path) -> None:
    """The shard is the half of the design that carries the evidence.

    Most assertions in this file read `emitter.checksums`, which records paths
    and digests — so all of them together verify that shard files exist at the
    right URLs and none of them that the files contain anything. Publishing
    every shard as `{"table": ..., "rows": []}` passes the rest of this suite:
    the bundle still advertises counts and shard links, every checksum still
    verifies, and every gene's omics tab renders blank. That is this module's own
    stated failure — evidence absent while the build reports success — with the
    empty half moved from the summary to the file it links to.

    Reading the file back is also what pins the shard's URL, so the separate test
    that only asserted `"omics/expression/PXD012345.json" in emitter.checksums`
    is gone: every mutation of the path it caught moves this file too, and this
    one additionally notices an empty or misnamed payload at that URL.

    `n_case` is asserted to still be an `int` because JSON has one number type
    and the atlas has coordinates and sample counts in it: a count republished as
    10.0 is not the value that was curated.
    """
    root = _repo(tmp_path, _row("HGNC:11604", "0.01") + _row("HGNC:11604", "0.02"))
    emitter = Emitter(root=tmp_path / "dist")

    build_omics(root, emitter)

    shard = tmp_path / "dist" / "omics" / "expression" / "PXD012345.json"
    payload = json.loads(shard.read_bytes())
    assert payload["table"] == "expression"
    assert [row["gene"] for row in payload["rows"]] == ["HGNC:11604", "HGNC:11604"]
    assert [row["fdr"] for row in payload["rows"]] == [0.01, 0.02]
    assert payload["rows"][0]["n_case"] == 10
    assert isinstance(payload["rows"][0]["n_case"], int)


def test_summarises_each_gene_by_modality(tmp_path: Path) -> None:
    root = _repo(tmp_path, _row("HGNC:11604", "0.01") + _row("HGNC:11604", "0.02"))
    emitter = Emitter(root=tmp_path / "dist")

    summaries = build_omics(root, emitter)

    assert summaries["HGNC:11604"]["expression"]["count"] == 2
    assert summaries["HGNC:11604"]["expression"]["shards"] == ["omics/expression/PXD012345.json"]


def test_the_summary_keeps_only_the_most_significant_rows(tmp_path: Path) -> None:
    """A well-studied gene can carry thousands of rows across datasets.

    Embedding them all is what the linked-shard decision exists to avoid, so the
    summary keeps a bounded, deterministically chosen head and the shard keeps
    everything.

    The rows are written *descending* by FDR, and the whole slice is spelled out
    rather than compared against itself. Written ascending — as the file already
    sits in FDR order — an implementation that never sorts at all still satisfies
    `top[0]["fdr"] < top[-1]["fdr"]`, and so does one that truncates before it
    ranks. Both were measured passing against exactly that assertion.
    """
    rows = "".join(_row("HGNC:11604", f"0.{index:03d}") for index in range(TOP_N + 9, 0, -1))
    root = _repo(tmp_path, rows)
    emitter = Emitter(root=tmp_path / "dist")

    summaries = build_omics(root, emitter)

    top = summaries["HGNC:11604"]["expression"]["top"]
    # Ranked by FDR ascending, so the most significant row leads — and the rows
    # kept are the most significant ones, not the first ones encountered.
    assert [row["fdr"] for row in top] == [float(f"0.{index:03d}") for index in range(1, TOP_N + 1)]
    assert summaries["HGNC:11604"]["expression"]["count"] == TOP_N + 9
    # The other half of the bargain: the bundle is bounded because the shard is
    # not. Truncating the published rows too would make the count a promise the
    # site cannot keep.
    shard = tmp_path / "dist" / "omics" / "expression" / "PXD012345.json"
    assert len(json.loads(shard.read_bytes())["rows"]) == TOP_N + 9


def test_the_slice_is_the_best_rows_across_every_dataset(tmp_path: Path) -> None:
    """Truncating per file would publish the best rows of the last file read.

    A gene's strongest evidence is as likely to sit in the second dataset as the
    first, so the ranking has to happen after every shard has contributed.
    """
    (tmp_path / "mirrors" / "expression").mkdir(parents=True)
    (tmp_path / "mirrors" / "expression" / "PXD000001.tsv").write_text(
        EXPRESSION_HEADER
        + "".join(_row("HGNC:11604", "0.500", dataset="PXD000001") for _ in range(TOP_N))
    )
    (tmp_path / "mirrors" / "expression" / "PXD000002.tsv").write_text(
        EXPRESSION_HEADER
        + "".join(_row("HGNC:11604", f"0.00{index}", dataset="PXD000002") for index in range(1, 4))
    )
    emitter = Emitter(root=tmp_path / "dist")

    summaries = build_omics(tmp_path, emitter)

    top = summaries["HGNC:11604"]["expression"]["top"]
    assert len(top) == TOP_N
    assert [row["dataset"] for row in top[:3]] == ["PXD000002"] * 3


def test_a_repository_with_no_omics_tables_summarises_nothing(tmp_path: Path) -> None:
    (tmp_path / "mirrors").mkdir()
    emitter = Emitter(root=tmp_path / "dist")

    assert build_omics(tmp_path, emitter) == {}
    assert emitter.checksums == {}


def test_tables_that_are_not_dataset_linked_get_no_shard(tmp_path: Path) -> None:
    """`variants` carries a gene column too, and is not omics.

    `mirror_paths` yields every mirror table, so what keeps the variant corpus
    out of `omics/` is the membership test rather than anything about the data.
    """
    _registry(tmp_path, _gene_row("HGNC:11604", "TBX5", "Q99593"))
    _table(
        tmp_path,
        "variants",
        "chr12.tsv",
        "vrs_id\tassembly\tchrom\tpos\tref\talt\thgvs_g\thgvs_c\thgvs_p\t"
        "mane_transcript\tgene\tconsequence\tclinvar_vcv\tclinvar_significance\t"
        "clinvar_review_status\tgnomad_af_popmax\tgnomad_version\tacmg_criteria\t"
        "inheritance\tde_novo\tzygosity\tphenotypes\tpublication\tlocator\t"
        "source_tier\tsource\tregulatory_element\ttarget_gene\tassay_evidence\n"
        "ga4gh:VA.1\tGRCh38\t12\t114400000\tA\tG\t\t\t\t\tHGNC:11604\tSO:0001583\t"
        "\t\t\t\t\t\t\t\t\t\t\t\town_curation\ts\t\t\t\n",
    )
    emitter = Emitter(root=tmp_path / "dist")

    summaries = build_omics(tmp_path, emitter)

    assert emitter.checksums == {}
    assert summaries == {}


def test_an_isoform_specific_registry_entry_wins_over_the_canonical_one(
    tmp_path: Path,
) -> None:
    """Stripping the suffix must be the fallback, not the rule.

    Where the registry names an isoform explicitly, that entry is the curator's
    answer for it; reaching for the canonical accession first would overrule a
    deliberate curation with a guess.
    """
    _registry(
        tmp_path,
        _gene_row("HGNC:11604", "TBX5", "Q99593") + _gene_row("HGNC:99999", "TBX5B", "Q99593-2"),
    )
    _table(tmp_path, "phospho", "PXD012345.tsv", PHOSPHO_HEADER + _phospho_row("Q99593-2", 12))
    emitter = Emitter(root=tmp_path / "dist")

    summaries = build_omics(tmp_path, emitter)

    assert sorted(summaries) == ["HGNC:99999"]


# One row per direction the isoform suffix can disagree in, and one mutant each
# that only that row catches: strip the suffix from the row only and the second
# row reaches no gene; strip it from the registry only and the first does not.
# A third row naming two *different* isoforms ("Q99593-3" against "Q99593-2")
# used to sit here and caught neither of those on its own, so it is gone.
@pytest.mark.parametrize(
    ("registry", "reported"),
    [("Q99593", "Q99593-2"), ("Q99593-2", "Q99593")],
)
def test_a_phospho_row_reaches_its_gene_whichever_side_names_the_isoform(
    tmp_path: Path, registry: str, reported: str
) -> None:
    """Phospho rows are keyed on protein; only genes.tsv can map them to a gene.

    Without the join a gene's phospho evidence would silently be absent from its
    bundle while sitting in a shard nothing links to. The plain case — both sides
    naming the same canonical accession — is exercised by
    `test_the_tie_break_orders_a_position_as_a_number`, so it is not repeated
    here; what needs its own fixture is the suffix disagreeing.

    `UNIPROT_PATTERN` admits the isoform form in `genes.uniprot` exactly as it
    does in `phospho.protein`, and both directions occur: a quantified site is
    reported against the isoform it was mapped to while the registry usually
    records the canonical accession, and a curator recording "Q99593-2" for the
    gene makes every canonical-accession row invisible to it. Matching the two
    literally drops every such row — measured against the unstripped join, both
    fixtures below summarised no gene at all while their shard was still written
    and published.

    Normalising one side only is the half-fix that looks right and is not: it
    answers whichever direction it was written for and silently drops the other.
    The exact-match precedence that must survive all this is covered above.
    """
    _registry(tmp_path, _gene_row("HGNC:11604", "TBX5", registry))
    _table(
        tmp_path,
        "phospho",
        "PXD012345.tsv",
        PHOSPHO_HEADER + _phospho_row(reported, 12),
    )
    emitter = Emitter(root=tmp_path / "dist")

    summaries = build_omics(tmp_path, emitter)

    assert summaries["HGNC:11604"]["phospho"]["count"] == 1


def test_a_proteomics_row_with_no_gene_falls_back_to_the_protein_join(
    tmp_path: Path,
) -> None:
    """`proteomics.gene` is nullable, unlike `expression.gene`.

    A row identified only by its accession is exactly as invisible to its gene as
    a phospho row would be, so dropping it rather than joining on `protein` loses
    real evidence into a shard nothing links to.
    """
    _registry(tmp_path, _gene_row("HGNC:11604", "TBX5", "Q99593"))
    _table(
        tmp_path,
        "proteomics",
        "PXD012345.tsv",
        PROTEOMICS_HEADER
        + "PXD012345\ttof_vs_control\tQ99593\t\t1.0\t0.001\t0.01\tup\t3\t20.0\tdia\n",
    )
    emitter = Emitter(root=tmp_path / "dist")

    summaries = build_omics(tmp_path, emitter)

    assert summaries["HGNC:11604"]["proteomics"]["count"] == 1


def test_a_row_with_no_protein_is_attributed_to_no_gene(tmp_path: Path) -> None:
    """Coercing the cell with `str()` invents a join out of a null.

    A null `protein` becomes the literal "None", which then matches any registry
    row whose `uniprot` reads "None" — and the site is filed under a gene that has
    nothing to do with it. Measured against the coercing join, this fixture
    reported one phospho site for TBX5.
    """
    _registry(tmp_path, _gene_row("HGNC:11604", "TBX5", "None"))
    _table(
        tmp_path,
        "phospho",
        "PXD012345.tsv",
        PHOSPHO_HEADER + "PXD012345\ttof_vs_control\tX_S12\t\tS\t12\tMOD:00046\t\t"
        "1.2\t0.001\t0.01\ttrue\t\t\n",
    )
    emitter = Emitter(root=tmp_path / "dist")

    summaries = build_omics(tmp_path, emitter)

    assert summaries == {}


def test_a_repository_with_no_gene_registry_summarises_no_phospho(
    tmp_path: Path,
) -> None:
    """The join has one input, and it is optional on disk.

    Nothing should raise: the shard is still published, and the rows are simply
    attributed to no gene until the registry names their accessions.
    """
    _table(tmp_path, "phospho", "PXD012345.tsv", PHOSPHO_HEADER + _phospho_row("Q99593", 12))
    emitter = Emitter(root=tmp_path / "dist")

    summaries = build_omics(tmp_path, emitter)

    assert summaries == {}
    assert "omics/phospho/PXD012345.json" in emitter.checksums


def test_an_unreadable_gene_registry_does_not_abort_the_build(tmp_path: Path) -> None:
    """A zero-length `genes.tsv` reads as `None`, not as a table with no rows.

    Distinct from the registry being absent, which is checked above: this path
    reaches `read_table` and comes back empty-handed. Reaching for `.columns` on
    that result raises `AttributeError` out of the build entirely — every shard
    in the atlas lost to one unreadable file that `validate_table` already
    reports against by name.
    """
    (tmp_path / "mirrors").mkdir()
    (tmp_path / "mirrors" / "genes.tsv").write_text("")
    _table(tmp_path, "phospho", "PXD012345.tsv", PHOSPHO_HEADER + _phospho_row("Q99593", 12))
    emitter = Emitter(root=tmp_path / "dist")

    summaries = build_omics(tmp_path, emitter)

    assert summaries == {}
    assert "omics/phospho/PXD012345.json" in emitter.checksums


def test_a_registry_missing_its_uniprot_column_does_not_abort_the_build(
    tmp_path: Path,
) -> None:
    """A column renamed upstream is the failure this guard exists for.

    `read_table` returns the frame whatever its columns are, so selecting a
    column that is not there raises `ColumnNotFoundError` out of `build_omics`
    and takes the whole atlas with it — over a registry that `validate_table`
    reports on precisely, by name, as TBL001. The unreadable-file half of the
    same guard is checked above; this is the half that survives a live upstream
    rename.
    """
    (tmp_path / "mirrors").mkdir()
    (tmp_path / "mirrors" / "genes.tsv").write_text(
        GENES_HEADER.replace("uniprot", "swissprot") + _gene_row("HGNC:11604", "TBX5", "Q99593")
    )
    _table(tmp_path, "phospho", "PXD012345.tsv", PHOSPHO_HEADER + _phospho_row("Q99593", 12))
    emitter = Emitter(root=tmp_path / "dist")

    summaries = build_omics(tmp_path, emitter)

    assert summaries == {}
    assert "omics/phospho/PXD012345.json" in emitter.checksums


def test_a_gene_named_by_a_row_outranks_the_accession_registry(tmp_path: Path) -> None:
    """`proteomics` carries both a gene and a protein, and they can disagree.

    The column is the curated answer for that row; the registry is an inference
    from an accession that several genes may claim. Routing every proteomics row
    through the join regardless would silently reattribute the row below to a
    gene its own dataset does not name — and the only proteomics test above uses
    a null gene, so nothing else exercises the filled cell at all.
    """
    _registry(tmp_path, _gene_row("HGNC:99999", "OTHER", "Q99593"))
    _table(
        tmp_path,
        "proteomics",
        "PXD012345.tsv",
        PROTEOMICS_HEADER
        + "PXD012345\ttof_vs_control\tQ99593\tHGNC:11604\t1.0\t0.001\t0.01\tup\t3\t20.0\tdia\n",
    )
    emitter = Emitter(root=tmp_path / "dist")

    summaries = build_omics(tmp_path, emitter)

    assert sorted(summaries) == ["HGNC:11604"]


def test_a_table_with_no_fdr_is_ranked_by_its_own_sort_key_nulls_first(
    tmp_path: Path,
) -> None:
    """`profiles` is an abundance table: it reports no significance at all.

    Ranking on a column that is not there has to fall through to something
    defined rather than to whichever row polars happened to yield first, so the
    tie-break is not a tie-break here but the whole ranking — which is why this
    is the fixture that pins it. `profiles.stage` is also the only nullable
    column in any omics sort key, and nulls first is what
    `validate/sort_order.py::_precedes` treats as canonical.

    Both halves in one fixture of three rows, written in the reverse of the
    table's sort key so that leaving nulls last, ordering the key the other way
    and not sorting at all each produce a different answer from the one below.
    `tissue` precedes `stage` in `profiles.sort_key`, so LV sorts ahead of both
    RV rows regardless of their stage, and the null stage leads within RV.
    """
    _table(
        tmp_path,
        "profiles",
        "GSE000001.tsv",
        PROFILES_HEADER
        + "GSE000001\tHGNC:11604\tRV\tinfant\t5.0\ttpm\t\t\t3\n"
        + "GSE000001\tHGNC:11604\tRV\t\t9.0\ttpm\t\t\t3\n"
        + "GSE000001\tHGNC:11604\tLV\tinfant\t7.0\ttpm\t\t\t3\n",
    )
    emitter = Emitter(root=tmp_path / "dist")

    summaries = build_omics(tmp_path, emitter)

    top = summaries["HGNC:11604"]["profiles"]["top"]
    assert [(row["tissue"], row["stage"]) for row in top] == [
        ("LV", "infant"),
        ("RV", None),
        ("RV", "infant"),
    ]


def test_the_tie_break_orders_a_position_as_a_number(tmp_path: Path) -> None:
    """Rows tie on FDR constantly, so the tie-break decides which ones survive.

    It borrows the table's own sort key, which `validate/sort_order.py` enforces
    with typed comparisons — `position` 9 before 100. Stringifying the values
    reverses that, and the rows written here are ordered so that both the
    stringified key and no key at all produce [100, 11, 9].
    """
    _registry(tmp_path, _gene_row("HGNC:11604", "TBX5", "Q99593"))
    _table(
        tmp_path,
        "phospho",
        "PXD012345.tsv",
        PHOSPHO_HEADER
        + _phospho_row("Q99593", 100)
        + _phospho_row("Q99593", 11)
        + _phospho_row("Q99593", 9),
    )
    emitter = Emitter(root=tmp_path / "dist")

    summaries = build_omics(tmp_path, emitter)

    top = summaries["HGNC:11604"]["phospho"]["top"]
    assert [row["position"] for row in top] == [9, 11, 100]


def test_a_row_reporting_no_fdr_ranks_below_every_row_that_does(
    tmp_path: Path,
) -> None:
    """An empty `fdr` cell is a row with no claim to significance, not the best one.

    `expression.fdr` is non-nullable, but the build reads the mirror directly and
    never re-runs `validate_table`, so an empty cell reaches the ranking. Ranking
    it first would put the one row carrying no statistic at the head of the
    slice; ranking it as zero does the same thing. Both are killed by spelling
    the whole order out, and nothing else in this file mixes rows with and
    without an FDR inside one modality.

    The three rows also carry three different contrasts, ordered so that the
    table's own sort key disagrees with the FDR order. That is what pins
    significance as the *first* component of `_rank` rather than merely a
    component of it: with the canonical key promoted above significance the rows
    come back a_vs, b_vs, c_vs — which is FDR order 0.02, None, 0.01 — and every
    row sharing one contrast, as they did before, cannot tell the two apart.
    """
    root = _repo(
        tmp_path,
        _row("HGNC:11604", "0.02", contrast="a_vs_control")
        + _row("HGNC:11604", "", contrast="b_vs_control")
        + _row("HGNC:11604", "0.01", contrast="c_vs_control"),
    )
    emitter = Emitter(root=tmp_path / "dist")

    summaries = build_omics(root, emitter)

    top = summaries["HGNC:11604"]["expression"]["top"]
    assert [row["fdr"] for row in top] == [0.01, 0.02, None]


def test_shard_paths_go_through_the_identifier_path_rule(tmp_path: Path) -> None:
    """The shard path embeds a filename, and a filename is not a URL.

    Every accession grammar passes `slug` unchanged, so this only fires on a name
    that should not have been curated — but published raw, a space reaches the
    manifest as a URL that has to be escaped before anything can fetch it.
    """
    (tmp_path / "mirrors" / "expression").mkdir(parents=True)
    (tmp_path / "mirrors" / "expression" / "PXD 012345.tsv").write_text(
        EXPRESSION_HEADER + _row("HGNC:11604", "0.01")
    )
    emitter = Emitter(root=tmp_path / "dist")

    build_omics(tmp_path, emitter)

    assert list(emitter.checksums) == ["omics/expression/PXD_012345.json"]


def test_the_shards_of_one_modality_are_listed_in_a_content_determined_order(
    tmp_path: Path,
) -> None:
    """`shards` is a JSON array, and `encode_json` does not order arrays.

    An array assembled in filesystem order is the failure emit.py names as the
    likeliest way this build stops being byte-identical. `mirror_paths` sorts the
    paths, but on the filename rather than on the slug of it: a space sorts before
    a hyphen and the underscore it becomes sorts after one, so these two shards
    are discovered in the opposite of the order they are published in.
    """
    for filename in ("PXD 000002.tsv", "PXD-000001.tsv"):
        (tmp_path / "mirrors" / "expression").mkdir(parents=True, exist_ok=True)
        (tmp_path / "mirrors" / "expression" / filename).write_text(
            EXPRESSION_HEADER + _row("HGNC:11604", "0.01")
        )
    emitter = Emitter(root=tmp_path / "dist")

    summaries = build_omics(tmp_path, emitter)

    assert summaries["HGNC:11604"]["expression"]["shards"] == [
        "omics/expression/PXD-000001.json",
        "omics/expression/PXD_000002.json",
    ]


def test_genes_sharing_an_accession_are_summarised_in_a_stable_order(
    tmp_path: Path,
) -> None:
    """One UniProt entry covers every gene encoding an identical protein.

    The histone clusters are the standard case, and the six genes below are one:
    H4C1 through H4C16 all resolve to P62805. Nothing in the schema makes
    `uniprot` unique — it is not part of the sort key — so a map of accession to
    *one* gene keeps whichever row was read last and drops the rest of the
    cluster's phospho evidence in silence. Asserting the whole key list rather
    than a membership check is what catches that here, which is why the separate
    two-gene test that only checked both genes were present is gone.

    That the accession maps to several genes is half of it; the order it maps to
    them in is the other half, because that order decides which gene `summaries`
    gains first, and so its order.

    The registry is written in reverse HGNC order, which is what gives this
    assertion teeth: dropping the `sorted()` leaves the file's own row order, and
    this fails on every run rather than on a fraction of them. Accumulating the
    cluster through a `set` instead would make the same regression depend on
    `PYTHONHASHSEED` — measured one process per seed, an unsorted set still came
    out in sorted order on 103/200 seeds at two genes, 48/200 at three and 26/200
    at six, because these ids differ only in their last digit and their hashes
    are not independent. A guard that waves a determinism regression through one
    run in eight is why the map does not use one.
    """
    _registry(
        tmp_path,
        "".join(
            _gene_row(f"HGNC:{4781 + index}", f"H4C{index}", "P62805")
            for index in reversed(range(6))
        ),
    )
    _table(tmp_path, "phospho", "PXD012345.tsv", PHOSPHO_HEADER + _phospho_row("P62805", 2))
    emitter = Emitter(root=tmp_path / "dist")

    summaries = build_omics(tmp_path, emitter)

    assert list(summaries) == [f"HGNC:{4781 + index}" for index in range(6)]


def test_an_unreadable_mirror_does_not_stop_the_others(tmp_path: Path) -> None:
    """A zero-length shard makes polars raise, and one bad file is not the build.

    `validate_table` reports the same file against the same path, so raising here
    would report it twice and abort every shard behind it.
    """
    _table(tmp_path, "expression", "PXD000001.tsv", "")
    _table(tmp_path, "expression", "PXD000002.tsv", EXPRESSION_HEADER + _row("HGNC:11604", "0.01"))
    emitter = Emitter(root=tmp_path / "dist")

    summaries = build_omics(tmp_path, emitter)

    assert list(emitter.checksums) == ["omics/expression/PXD000002.json"]
    assert summaries["HGNC:11604"]["expression"]["count"] == 1


def test_the_accession_index_survives_a_repeated_pair_and_a_blank_cell(tmp_path: Path) -> None:
    """Two registry states nothing forbids, each of which corrupts the join quietly.

    `mirrors/genes.tsv` has no uniqueness rule over `(hgnc_id, uniprot)`, so a
    generated registry may list the same pair twice — a merge of two upstream
    releases, or a gene reported under two transcripts. Without the
    duplicate check in `_record`, that gene is filed twice under one accession
    and every phospho row on it counts double, so a bundle advertises twice the
    evidence the shard holds.

    `uniprot` is nullable, and `str(None)` is the string "None". Without the null
    skip, every gene with a blank cell is filed under one fictional accession
    "None" — so a row whose protein failed to parse into anything else joins to
    all of them at once, which is amendment A12(b) in a second place.

    Read through `build_omics` rather than against the index directly, because
    the count in the summary is what a reader sees and the doubling is only
    visible there.
    """
    _registry(
        tmp_path,
        _gene_row("HGNC:11604", "TBX5", "Q99593")
        + _gene_row("HGNC:11604", "TBX5", "Q99593")
        + _gene_row("HGNC:4173", "GATA4", ""),
    )
    _table(
        tmp_path,
        "phospho",
        "PXD012345.tsv",
        PHOSPHO_HEADER + _phospho_row("Q99593", 100) + _phospho_row("Q99593", 200),
    )
    emitter = Emitter(root=tmp_path / "dist")

    summaries = build_omics(tmp_path, emitter)

    assert summaries["HGNC:11604"]["phospho"]["count"] == 2, "the repeated pair doubled the count"

    # The blank cell is checked against the index rather than the summary: with
    # the null skip removed the gene is filed under the literal key "None", and
    # nothing surfaces until a row's accession canonicalises to that same string
    # — at which point every gene with a blank cell joins it at once.
    index = _accession_index(tmp_path)
    assert "None" not in index.literal
    assert "None" not in index.canonical
    assert index.literal["Q99593"] == ("HGNC:11604",)


def test_a_phospho_shard_can_be_filtered_to_the_gene_its_bundle_counted(tmp_path: Path) -> None:
    """Issue #3: the rows a summary counts are now selectable out of the shard.

    `phospho` has no gene column at all, and the only route from a row's
    accession to a gene is `mirrors/genes.tsv`, which the site does not publish.
    So a consumer following `shards` from a bundle could fetch the file and had
    no way to pick out the rows it had just been told the count of — the `count`
    was a promise the shard could not keep.

    Each row now carries the attribution the summary was built from, and this
    asserts the property that makes it worth publishing: **filtering the shard on
    `genes` reproduces the count exactly**. Not that the field is present — that
    a consumer following the documented route arrives at the same number.

    Two genes on one accession is the case a scalar field would have got wrong.
    P62805 is a histone cluster in reality; A12c records the build dropping all
    but one of them once already.
    """
    _registry(
        tmp_path,
        _gene_row("HGNC:11604", "TBX5", "Q99593") + _gene_row("HGNC:4173", "GATA4", "Q99593"),
    )
    _table(
        tmp_path,
        "phospho",
        "PXD012345.tsv",
        PHOSPHO_HEADER + _phospho_row("Q99593", 100) + _phospho_row("Q99593", 200),
    )
    emitter = Emitter(root=tmp_path / "dist")

    summaries = build_omics(tmp_path, emitter)

    shard = json.loads((tmp_path / "dist" / "omics" / "phospho" / "PXD012345.json").read_text())
    for gene in ("HGNC:11604", "HGNC:4173"):
        selected = [row for row in shard["rows"] if gene in row["genes"]]
        assert len(selected) == summaries[gene]["phospho"]["count"], (
            f"a consumer filtering the shard for {gene} does not get the count it was promised"
        )
    # Both genes on one accession, so neither may displace the other.
    assert shard["rows"][0]["genes"] == ["HGNC:11604", "HGNC:4173"]


def test_an_expression_row_carries_the_same_attribution_as_a_phospho_row(
    tmp_path: Path,
) -> None:
    """Written for every modality, so a consumer filters one way everywhere.

    `expression` rows have their own `gene` column and did not need this. Adding
    it anyway is what makes the guarantee uniform: `count` is derived from the
    published attribution for every modality, so there is no table where a
    consumer must know to filter differently, and no table where the row's own
    column and the attribution could come to disagree.
    """
    _registry(tmp_path, _gene_row("HGNC:11604", "TBX5", "Q99593"))
    repo = _repo(tmp_path, _row(gene="HGNC:11604", fdr="0.01"))
    emitter = Emitter(root=tmp_path / "dist")

    summaries = build_omics(repo, emitter)

    shard = json.loads((tmp_path / "dist" / "omics" / "expression" / "PXD012345.json").read_text())
    assert shard["rows"][0]["genes"] == ["HGNC:11604"]
    selected = [row for row in shard["rows"] if "HGNC:11604" in row["genes"]]
    assert len(selected) == summaries["HGNC:11604"]["expression"]["count"]
