# tests/unit/test_build_omics.py
from pathlib import Path

from chd_atlas.build.emit import Emitter
from chd_atlas.build.omics import TOP_N, build_omics

EXPRESSION_HEADER = (
    "dataset\tcontrast\tgene\tlog2fc\tpvalue\tfdr\tdirection\t"
    "n_case\tn_control\ttissue\tstage\n"
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


def _row(gene: str, fdr: str, log2fc: str = "1.0", dataset: str = "PXD012345") -> str:
    return (
        f"{dataset}\ttof_vs_control\t{gene}\t{log2fc}\t0.001\t{fdr}\tup\t"
        f"10\t10\tRV\tinfant\n"
    )


def _phospho_row(protein: str, position: int, fdr: str = "0.01") -> str:
    return (
        f"PXD012345\ttof_vs_control\t{protein}_S{position}\t{protein}\tS\t{position}\t"
        f"MOD:00046\t\t1.2\t0.001\t{fdr}\ttrue\t\t\n"
    )


def _repo(tmp_path: Path, rows: str) -> Path:
    (tmp_path / "mirrors" / "expression").mkdir(parents=True)
    (tmp_path / "mirrors" / "expression" / "PXD012345.tsv").write_text(
        EXPRESSION_HEADER + rows
    )
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


def test_emits_one_shard_per_dataset_table(tmp_path: Path) -> None:
    root = _repo(tmp_path, _row("HGNC:11604", "0.01"))
    emitter = Emitter(root=tmp_path / "dist")

    build_omics(root, emitter)

    assert "omics/expression/PXD012345.json" in emitter.checksums


def test_summarises_each_gene_by_modality(tmp_path: Path) -> None:
    root = _repo(tmp_path, _row("HGNC:11604", "0.01") + _row("HGNC:11604", "0.02"))
    emitter = Emitter(root=tmp_path / "dist")

    summaries = build_omics(root, emitter)

    assert summaries["HGNC:11604"]["expression"]["count"] == 2
    assert summaries["HGNC:11604"]["expression"]["shards"] == [
        "omics/expression/PXD012345.json"
    ]


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
    rows = "".join(
        _row("HGNC:11604", f"0.{index:03d}") for index in range(TOP_N + 9, 0, -1)
    )
    root = _repo(tmp_path, rows)
    emitter = Emitter(root=tmp_path / "dist")

    summaries = build_omics(root, emitter)

    top = summaries["HGNC:11604"]["expression"]["top"]
    # Ranked by FDR ascending, so the most significant row leads — and the rows
    # kept are the most significant ones, not the first ones encountered.
    assert [row["fdr"] for row in top] == [
        float(f"0.{index:03d}") for index in range(1, TOP_N + 1)
    ]
    assert summaries["HGNC:11604"]["expression"]["count"] == TOP_N + 9


def test_the_slice_is_the_best_rows_across_every_dataset(tmp_path: Path) -> None:
    """Truncating per file would publish the best rows of the last file read.

    A gene's strongest evidence is as likely to sit in the second dataset as the
    first, so the ranking has to happen after every shard has contributed.
    """
    (tmp_path / "mirrors" / "expression").mkdir(parents=True)
    (tmp_path / "mirrors" / "expression" / "PXD000001.tsv").write_text(
        EXPRESSION_HEADER
        + "".join(
            _row("HGNC:11604", "0.500", dataset="PXD000001") for _ in range(TOP_N)
        )
    )
    (tmp_path / "mirrors" / "expression" / "PXD000002.tsv").write_text(
        EXPRESSION_HEADER
        + "".join(
            _row("HGNC:11604", f"0.00{index}", dataset="PXD000002")
            for index in range(1, 4)
        )
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


def test_phospho_rows_reach_a_gene_through_the_uniprot_column(tmp_path: Path) -> None:
    """Phospho rows are keyed on protein; only genes.tsv can map them to a gene.

    Without the join a gene's phospho evidence would silently be absent from its
    bundle while sitting in a shard nothing links to.
    """
    _registry(tmp_path, _gene_row("HGNC:11604", "TBX5", "Q99593"))
    _table(tmp_path, "phospho", "PXD012345.tsv", PHOSPHO_HEADER + _phospho_row("Q99593", 12))
    emitter = Emitter(root=tmp_path / "dist")

    summaries = build_omics(tmp_path, emitter)

    assert summaries["HGNC:11604"]["phospho"]["count"] == 1


def test_a_phospho_row_naming_an_isoform_still_reaches_its_gene(tmp_path: Path) -> None:
    """`UNIPROT_PATTERN` admits "Q99593-2", and phospho data routinely uses it.

    A quantified site is reported against the isoform it was mapped to, while the
    registry records the gene's canonical accession. Matching the two literally
    drops every such row: measured against the unstripped join, this fixture
    summarised no gene at all while the shard was still written and published.
    """
    _registry(tmp_path, _gene_row("HGNC:11604", "TBX5", "Q99593"))
    _table(
        tmp_path, "phospho", "PXD012345.tsv", PHOSPHO_HEADER + _phospho_row("Q99593-2", 12)
    )
    emitter = Emitter(root=tmp_path / "dist")

    summaries = build_omics(tmp_path, emitter)

    assert summaries["HGNC:11604"]["phospho"]["count"] == 1


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
        _gene_row("HGNC:11604", "TBX5", "Q99593")
        + _gene_row("HGNC:99999", "TBX5B", "Q99593-2"),
    )
    _table(
        tmp_path, "phospho", "PXD012345.tsv", PHOSPHO_HEADER + _phospho_row("Q99593-2", 12)
    )
    emitter = Emitter(root=tmp_path / "dist")

    summaries = build_omics(tmp_path, emitter)

    assert sorted(summaries) == ["HGNC:99999"]


def test_an_accession_shared_by_several_genes_reaches_all_of_them(
    tmp_path: Path,
) -> None:
    """One UniProt entry covers every gene encoding an identical protein.

    The histone clusters are the standard case: H4C1 through H4C16 all resolve to
    P62805. Nothing in the schema makes `uniprot` unique — it is not part of the
    sort key — so a map of accession to *one* gene keeps whichever row was read
    last and drops the rest of the cluster's phospho evidence in silence.
    """
    _registry(
        tmp_path,
        _gene_row("HGNC:4781", "H4C1", "P62805") + _gene_row("HGNC:4787", "H4C2", "P62805"),
    )
    _table(tmp_path, "phospho", "PXD012345.tsv", PHOSPHO_HEADER + _phospho_row("P62805", 2))
    emitter = Emitter(root=tmp_path / "dist")

    summaries = build_omics(tmp_path, emitter)

    assert summaries["HGNC:4781"]["phospho"]["count"] == 1
    assert summaries["HGNC:4787"]["phospho"]["count"] == 1


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
        PHOSPHO_HEADER
        + "PXD012345\ttof_vs_control\tX_S12\t\tS\t12\tMOD:00046\t\t"
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


def test_a_table_with_no_fdr_column_is_ranked_by_its_own_sort_key(
    tmp_path: Path,
) -> None:
    """`profiles` is an abundance table: it reports no significance at all.

    Ranking on a column that is not there has to fall through to something
    defined rather than to whichever row polars happened to yield first, so the
    rows below are written in the reverse of the table's sort key.
    """
    _table(
        tmp_path,
        "profiles",
        "GSE000001.tsv",
        PROFILES_HEADER
        + "GSE000001\tHGNC:11604\tRV\tinfant\t5.0\ttpm\t\t\t3\n"
        + "GSE000001\tHGNC:11604\tLV\tinfant\t9.0\ttpm\t\t\t3\n",
    )
    emitter = Emitter(root=tmp_path / "dist")

    summaries = build_omics(tmp_path, emitter)

    top = summaries["HGNC:11604"]["profiles"]["top"]
    assert [row["tissue"] for row in top] == ["LV", "RV"]


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
    """That order decides which gene `summaries` gains first, and so its order.

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
    _table(
        tmp_path, "expression", "PXD000002.tsv", EXPRESSION_HEADER + _row("HGNC:11604", "0.01")
    )
    emitter = Emitter(root=tmp_path / "dist")

    summaries = build_omics(tmp_path, emitter)

    assert list(emitter.checksums) == ["omics/expression/PXD000002.json"]
    assert summaries["HGNC:11604"]["expression"]["count"] == 1
