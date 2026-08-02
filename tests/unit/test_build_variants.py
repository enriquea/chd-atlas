# tests/unit/test_build_variants.py
import gzip
import json
from pathlib import Path

import pytest

from chd_atlas.build.emit import Emitter
from chd_atlas.build.variants import build_variants
from chd_atlas.tables import TABLE_SCHEMAS

HEADER = (
    "vrs_id\tassembly\tchrom\tpos\tref\talt\thgvs_g\thgvs_c\thgvs_p\tmane_transcript\t"
    "gene\tconsequence\tclinvar_vcv\tclinvar_significance\tclinvar_review_status\t"
    "gnomad_af_popmax\tgnomad_version\tacmg_criteria\tinheritance\tde_novo\tzygosity\t"
    "phenotypes\tpublication\tlocator\tsource_tier\tsource\tregulatory_element\t"
    "target_gene\tassay_evidence\n"
)


def _row(
    vrs_id: str = "ga4gh:VA.x",
    chrom: str = "12",
    pos: str = "114400000",
    gene: str = "HGNC:11604",
    frequency: str = "",
) -> str:
    """One variant row, built field by field so it cannot drift out of alignment."""
    fields = [
        vrs_id, "GRCh38", chrom, pos, "A", "G", "", "", "", "",
        gene, "SO:0001583", "", "", "", frequency, "", "", "", "",
        "", "", "", "", "third_party_literature", "clinvar", "", "", "",
    ]  # fmt: skip
    return "\t".join(fields) + "\n"


ROW = _row()


def _shard(root: Path, name: str, rows: str) -> None:
    (root / "mirrors" / "variants").mkdir(parents=True, exist_ok=True)
    (root / "mirrors" / "variants" / name).write_text(HEADER + rows)


def _repo(tmp_path: Path) -> Path:
    (tmp_path / "mirrors" / "variants").mkdir(parents=True)
    (tmp_path / "mirrors" / "variants" / "12.tsv").write_text(HEADER + ROW)
    return tmp_path


def _shards(tmp_path: Path) -> list[str]:
    index = json.loads((tmp_path / "dist" / "variants" / "index.json").read_text())
    shards: list[str] = index["shards"]
    return shards


def test_the_fixture_matches_the_schema_it_stands_in_for() -> None:
    """Every other assertion here reads a column by name, so alignment is load-bearing."""
    assert HEADER.rstrip("\n").split("\t") == list(TABLE_SCHEMAS["variants"].column_names)
    assert len(ROW.rstrip("\n").split("\t")) == len(TABLE_SCHEMAS["variants"].columns)


def test_emits_one_gzipped_shard_per_chromosome(tmp_path: Path) -> None:
    emitter = Emitter(root=tmp_path / "dist")

    build_variants(_repo(tmp_path), emitter)

    assert "variants/12.json.gz" in emitter.checksums
    written = (tmp_path / "dist" / "variants" / "12.json.gz").read_bytes()
    assert json.loads(gzip.decompress(written))["rows"][0]["gene"] == "HGNC:11604"


def test_a_shard_declares_the_chromosome_it_holds(tmp_path: Path) -> None:
    emitter = Emitter(root=tmp_path / "dist")

    build_variants(_repo(tmp_path), emitter)

    written = (tmp_path / "dist" / "variants" / "12.json.gz").read_bytes()
    assert json.loads(gzip.decompress(written))["chrom"] == "12"


def test_returns_the_variants_belonging_to_each_gene(tmp_path: Path) -> None:
    emitter = Emitter(root=tmp_path / "dist")

    by_gene = build_variants(_repo(tmp_path), emitter)

    assert [row["vrs_id"] for row in by_gene["HGNC:11604"]] == ["ga4gh:VA.x"]


def test_an_index_lists_the_shards_that_exist(tmp_path: Path) -> None:
    """The client must not guess which chromosomes were emitted.

    Fetching variants/3.json.gz for a chromosome with no curated variants is a
    404 the page has to handle; listing what exists removes the guess.
    """
    emitter = Emitter(root=tmp_path / "dist")

    build_variants(_repo(tmp_path), emitter)

    assert _shards(tmp_path) == ["variants/12.json.gz"]


def test_no_variants_curated_yet_emits_an_empty_index_not_a_crash(tmp_path: Path) -> None:
    """This is the repository's current state, so it must be the tested one."""
    (tmp_path / "mirrors").mkdir()
    emitter = Emitter(root=tmp_path / "dist")

    assert build_variants(tmp_path, emitter) == {}
    assert _shards(tmp_path) == []


def test_the_index_is_ordered_karyotypically_not_lexically(tmp_path: Path) -> None:
    """Sorting shard names as strings publishes 1, 10, 11, ..., 2, 20 — and MT before X.

    A consumer rendering a chromosome picker straight from this list would show
    that order, and the fix on its side is to reimplement the karyotype rule in
    JavaScript. Five chromosomes chosen so that lexical order differs from this
    expectation in both the numeric and the non-numeric part.
    """
    _shard(tmp_path, "10.tsv", _row(chrom="10"))
    _shard(tmp_path, "2.tsv", _row(chrom="2"))
    _shard(tmp_path, "X.tsv", _row(chrom="X"))
    _shard(tmp_path, "1.tsv", _row(chrom="1"))
    _shard(tmp_path, "MT.tsv", _row(chrom="MT"))
    emitter = Emitter(root=tmp_path / "dist")

    build_variants(tmp_path, emitter)

    assert _shards(tmp_path) == [
        "variants/1.json.gz",
        "variants/2.json.gz",
        "variants/10.json.gz",
        "variants/X.json.gz",
        "variants/MT.json.gz",
    ]


def test_a_genes_variants_follow_the_same_order_the_index_publishes(tmp_path: Path) -> None:
    """Across shards karyotypically, within a shard in the file's own sorted order.

    Both halves matter: the list is embedded verbatim in a gene bundle, whose
    bytes are checksummed, so any order that depends on which file the
    filesystem happened to yield first is a checksum that changes for free.
    """
    _shard(tmp_path, "10.tsv", _row(vrs_id="ga4gh:VA.c", chrom="10", pos="200"))
    _shard(
        tmp_path,
        "2.tsv",
        _row(vrs_id="ga4gh:VA.a", chrom="2", pos="100")
        + _row(vrs_id="ga4gh:VA.b", chrom="2", pos="300"),
    )
    emitter = Emitter(root=tmp_path / "dist")

    by_gene = build_variants(tmp_path, emitter)

    assert [row["vrs_id"] for row in by_gene["HGNC:11604"]] == [
        "ga4gh:VA.a",
        "ga4gh:VA.b",
        "ga4gh:VA.c",
    ]


def test_each_gene_gets_only_its_own_variants(tmp_path: Path) -> None:
    _shard(
        tmp_path,
        "12.tsv",
        _row(vrs_id="ga4gh:VA.a", gene="HGNC:11604")
        + _row(vrs_id="ga4gh:VA.b", pos="114400001", gene="HGNC:2201"),
    )
    emitter = Emitter(root=tmp_path / "dist")

    by_gene = build_variants(tmp_path, emitter)

    assert [row["vrs_id"] for row in by_gene["HGNC:11604"]] == ["ga4gh:VA.a"]
    assert [row["vrs_id"] for row in by_gene["HGNC:2201"]] == ["ga4gh:VA.b"]


def test_a_row_with_no_gene_is_published_but_belongs_to_no_gene(tmp_path: Path) -> None:
    """`gene` is not nullable, so a null one is a curation error that
    `validate_repository` reports and `build_site` refuses to build on. Reaching
    here regardless must not file the row under the gene "None", which is what
    coercing the key with `str()` produces — a fictional gene whose bundle would
    collect every unattributed variant in the atlas."""
    _shard(tmp_path, "12.tsv", _row(gene=""))
    emitter = Emitter(root=tmp_path / "dist")

    by_gene = build_variants(tmp_path, emitter)

    assert by_gene == {}
    written = (tmp_path / "dist" / "variants" / "12.json.gz").read_bytes()
    assert json.loads(gzip.decompress(written))["rows"][0]["gene"] is None


def test_a_shard_holding_another_chromosomes_rows_is_refused(tmp_path: Path) -> None:
    """The filename is what names the shard, so it must be what the rows say.

    Left alone, a chr7 variant filed in 12.tsv is published inside
    variants/12.json.gz under a `chrom` of "12": a client that maps a
    chromosome to its shard can never reach it, and nothing reports it missing.

    Checked before the shard is written rather than after: a refusal that leaves
    that artifact on disk has still published the thing this rule exists to
    prevent, and a later build reading `dist/` cannot tell it apart.
    """
    _shard(tmp_path, "12.tsv", _row(chrom="12") + _row(chrom="7", pos="1000"))
    emitter = Emitter(root=tmp_path / "dist")

    # On the rule's own words: the filename rule's message quotes the whole
    # chromosome vocabulary, so matching on "'7'" alone would be satisfied by
    # either rule firing.
    with pytest.raises(ValueError, match="holds rows for") as raised:
        build_variants(tmp_path, emitter)

    assert "'7'" in str(raised.value)
    # The filename is the one actionable token: this raises as a traceback in CI
    # rather than as a rendered report, so a message without it names no file to
    # go and fix.
    assert "12.tsv" in str(raised.value)
    assert not (tmp_path / "dist" / "variants" / "12.json.gz").exists()


def test_the_refusal_names_every_foreign_chromosome_in_a_fixed_order(tmp_path: Path) -> None:
    """The offending values go through a set, whose iteration order follows the
    interpreter's hash seed and so differs between two runs on identical data.
    Sorting them is what makes the message reproducible; five values are enough
    that an unsorted list matching this one by chance is a 1-in-120 event."""
    _shard(
        tmp_path,
        "12.tsv",
        "".join(
            _row(chrom=chrom, pos=str(1000 + index))
            for index, chrom in enumerate(["7", "X", "MT", "2", "21"])
        ),
    )
    emitter = Emitter(root=tmp_path / "dist")

    with pytest.raises(ValueError) as raised:
        build_variants(tmp_path, emitter)

    assert "holds rows for ['2', '21', '7', 'MT', 'X']" in str(raised.value)


def test_a_row_with_no_chromosome_is_reported_as_that(tmp_path: Path) -> None:
    """`chrom` is not nullable, so a null cell is a curation error `validate_table`
    reports as TBL003. Reaching here regardless, the message must not call the
    missing value a chromosome: `str(None)` names "None", which is no chromosome
    in any genome and sends a curator looking through the file for a value that
    is not in it. It is also the coercion this module refuses to make on `gene`
    twenty lines further down."""
    _shard(tmp_path, "12.tsv", _row(chrom=""))
    emitter = Emitter(root=tmp_path / "dist")

    with pytest.raises(ValueError, match="holds rows with no chromosome") as raised:
        build_variants(tmp_path, emitter)

    assert "None" not in str(raised.value)


def test_a_filename_that_is_not_a_chromosome_is_refused(tmp_path: Path) -> None:
    """`variants/chr12.json.gz` is a URL no consumer of the schema would construct.

    The shard is left empty, and the match is on the filename rule's own words,
    so that only that rule can satisfy this test. Filled with chr12 rows the
    shard would be refused for holding rows the filename does not name, and the
    test would pass just as well with the filename never checked at all — which
    is how it was first written and what mutating the check away exposed.

    A second, valid shard is what makes the last assertion mean anything: with
    one shard in the repository `dist/` is empty whether the name is resolved up
    front or in the middle of the loop, and resolving it in the loop publishes
    variants/1.json.gz before refusing the build.
    """
    _shard(tmp_path, "chr12.tsv", "")
    _shard(tmp_path, "1.tsv", _row(chrom="1"))
    emitter = Emitter(root=tmp_path / "dist")

    with pytest.raises(ValueError, match="'chr12' is not one of") as raised:
        build_variants(tmp_path, emitter)

    assert "chr12.tsv" in str(raised.value)
    # Refused before anything is written, rather than part-way through dist/.
    assert not (tmp_path / "dist").exists()


def test_a_chromosome_in_the_wrong_case_is_refused(tmp_path: Path) -> None:
    """`x` is not `X`, and case is the axis this codebase is most exposed on.

    `emit.py` carries a whole guard for it because macOS keeps one file where
    Linux keeps two, so `x.tsv` is a file a curator plausibly creates and CI
    plausibly treats as a second chromosome. Matching on the rule's own words
    because accepting the stem case-insensitively does not merely publish
    `variants/x.json.gz` — it degrades this clear refusal into a bare
    `KeyError: 'x'` thrown out of a sort key with nothing to say about mirrors.
    """
    _shard(tmp_path, "x.tsv", _row(chrom="X"))
    emitter = Emitter(root=tmp_path / "dist")

    with pytest.raises(ValueError, match="'x' is not one of"):
        build_variants(tmp_path, emitter)


def test_a_shard_with_a_header_and_no_rows_is_still_published(tmp_path: Path) -> None:
    """A curated-but-empty chromosome is a fact about the mirror, not an error."""
    _shard(tmp_path, "21.tsv", "")
    emitter = Emitter(root=tmp_path / "dist")

    assert build_variants(tmp_path, emitter) == {}
    assert _shards(tmp_path) == ["variants/21.json.gz"]
    written = (tmp_path / "dist" / "variants" / "21.json.gz").read_bytes()
    assert json.loads(gzip.decompress(written)) == {"chrom": "21", "rows": []}


def test_only_the_variant_mirrors_are_read(tmp_path: Path) -> None:
    """`mirror_paths` yields every table in the repository, and none of the others
    is named for a chromosome. Without the filter, `mirrors/genes.tsv` — which
    every real repository has, and which is yielded before any shard — is handed
    to the filename rule and refuses the build on a repository that is entirely
    correct. Both layouts are represented here: the flat tables and the
    per-accession shards are separate branches of `mirror_paths`.
    """
    _shard(tmp_path, "12.tsv", ROW)
    (tmp_path / "mirrors" / "genes.tsv").write_text(
        "hgnc_id\tsymbol\tname\taliases\tensembl_gene\tncbi_gene\tlocus\tuniprot\tmane_select\n"
        "HGNC:11604\tTBX5\tT-box transcription factor 5\t\t\t\t\t\t\n"
    )
    (tmp_path / "mirrors" / "expression").mkdir()
    (tmp_path / "mirrors" / "expression" / "PXD012345.tsv").write_text(
        "dataset\tcontrast\tgene\tlog2fc\tpvalue\tfdr\tdirection\tn_case\tn_control\t"
        "tissue\tstage\n"
    )
    emitter = Emitter(root=tmp_path / "dist")

    by_gene = build_variants(tmp_path, emitter)

    assert _shards(tmp_path) == ["variants/12.json.gz"]
    assert list(by_gene) == ["HGNC:11604"]


def test_an_unreadable_shard_leaves_the_rest_of_the_build_standing(tmp_path: Path) -> None:
    """`validate_table` reports the same file with far more to say about it."""
    _shard(tmp_path, "12.tsv", ROW)
    (tmp_path / "mirrors" / "variants" / "13.tsv").write_text("")
    emitter = Emitter(root=tmp_path / "dist")

    by_gene = build_variants(tmp_path, emitter)

    assert _shards(tmp_path) == ["variants/12.json.gz"]
    assert [row["vrs_id"] for row in by_gene["HGNC:11604"]] == ["ga4gh:VA.x"]


def test_a_non_finite_frequency_fails_the_build_instead_of_publishing_it(tmp_path: Path) -> None:
    """`NaN` is what R's `write.table` emits for a missing statistic, and no JSON
    parser accepts the bare token. `validate_repository` reports it as TBL010 and
    `build_site` refuses to run, so reaching `encode_json`'s `allow_nan=False`
    means validation was bypassed — it must still fail rather than publish."""
    _shard(tmp_path, "12.tsv", _row(frequency="NaN"))
    emitter = Emitter(root=tmp_path / "dist")

    with pytest.raises(ValueError):
        build_variants(tmp_path, emitter)

    assert not (tmp_path / "dist" / "variants" / "12.json.gz").exists()
