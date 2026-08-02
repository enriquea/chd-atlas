# tests/unit/test_referential.py
from datetime import date
from pathlib import Path

from chd_atlas.corpus import Corpus
from chd_atlas.issues import Severity
from chd_atlas.models.assertion import (
    Evidence,
    GeneDiseaseAssertion,
    SupplementaryLocator,
)
from chd_atlas.models.dataset import Dataset
from chd_atlas.models.functional import FunctionalEvidence
from chd_atlas.models.literature import FeaturedManuscript, PhenotypeTerm, Publication
from chd_atlas.tables import TABLE_SCHEMAS
from chd_atlas.validate.referential import (
    validate_mirror_references,
    validate_ptm_evidence_is_reachable,
    validate_references,
)


def _publication(pmid: str = "PMID:8988165") -> Publication:
    return Publication(
        id=pmid,
        title="T",
        journal="J",
        year=1997,
        authors=["A"],
        study_type="case_control",
    )


def _evidence(**overrides: object) -> Evidence:
    payload: dict[str, object] = {
        "publication": "PMID:8988165",
        "evidence_class": "genetic_case",
        "locator": SupplementaryLocator(file="T.xlsx"),
        "summary": "s",
        "strength": "strong",
    }
    payload.update(overrides)
    return Evidence.model_validate(payload)


def _assertion(**overrides: object) -> GeneDiseaseAssertion:
    payload: dict[str, object] = {
        "id": "CHDA:AST:0000001",
        "gene": "HGNC:11604",
        "phenotypes": ["HP:0001631"],
        "lesion_groups": ["septal"],
        "classification": "definitive",
        "inheritance": ["AD"],
        "mechanism": "haploinsufficiency",
        "syndromic": "both",
        "evidence": [_evidence()],
        "source_tier": "own_curation",
        "curator": "c",
        "curated_on": date(2026, 7, 1),
        "last_reviewed": date(2026, 7, 1),
    }
    payload.update(overrides)
    return GeneDiseaseAssertion.model_validate(payload)


def _corpus(**overrides: object) -> Corpus:
    payload: dict[str, object] = {
        "root": Path("."),
        "assertions": (_assertion(),),
        "publications": (_publication(),),
    }
    payload.update(overrides)
    return Corpus(**payload)  # type: ignore[arg-type]


def test_accepts_a_referentially_intact_corpus() -> None:
    assert validate_references(_corpus(), known_genes={"HGNC:11604"}) == []


def test_reports_an_assertion_gene_absent_from_the_gene_registry() -> None:
    issues = validate_references(_corpus(), known_genes=set())
    assert [i.code for i in issues] == ["REF001"]
    assert "HGNC:11604" in issues[0].message


def test_gene_checks_are_skipped_when_the_registry_is_unavailable() -> None:
    """None means genes.tsv did not load; hundreds of REF001s would bury the cause."""
    assert validate_references(_corpus(), known_genes=None) == []


def test_reports_evidence_citing_an_unknown_publication() -> None:
    corpus = _corpus(assertions=(_assertion(evidence=[_evidence(publication="PMID:1")]),))
    issues = validate_references(corpus, known_genes={"HGNC:11604"})
    assert [i.code for i in issues] == ["REF002"]


def test_reports_an_unresolvable_functional_evidence_reference() -> None:
    evidence = _evidence(evidence_class="functional_model", functional_evidence="CHDA:FUN:0000009")
    corpus = _corpus(assertions=(_assertion(evidence=[evidence]),))
    issues = validate_references(corpus, known_genes={"HGNC:11604"})
    assert [i.code for i in issues] == ["REF003"]


def test_resolves_a_present_functional_evidence_reference() -> None:
    functional = FunctionalEvidence.model_validate(
        {
            "id": "CHDA:FUN:0000009",
            "gene": "HGNC:11604",
            "organism": "NCBITaxon:10090",
            "perturbation": "knockout",
            "zygosity": "heterozygous",
            "cardiac_phenotype": ["MP:0010402"],
            "phenocopies_human": "partial",
            "rescue_outcome": "not_attempted",
            "publication": "PMID:8988165",
        }
    )
    evidence = _evidence(evidence_class="functional_model", functional_evidence="CHDA:FUN:0000009")
    corpus = _corpus(assertions=(_assertion(evidence=[evidence]),), functional=(functional,))
    assert validate_references(corpus, known_genes={"HGNC:11604"}) == []


def test_reports_a_functional_record_about_a_different_gene() -> None:
    functional = FunctionalEvidence.model_validate(
        {
            "id": "CHDA:FUN:0000009",
            "gene": "HGNC:4173",
            "organism": "NCBITaxon:10090",
            "perturbation": "knockout",
            "zygosity": "heterozygous",
            "cardiac_phenotype": ["MP:0010402"],
            "phenocopies_human": "partial",
            "rescue_outcome": "not_attempted",
            "publication": "PMID:8988165",
        }
    )
    evidence = _evidence(evidence_class="functional_model", functional_evidence="CHDA:FUN:0000009")
    corpus = _corpus(assertions=(_assertion(evidence=[evidence]),), functional=(functional,))

    issues = validate_references(corpus, known_genes={"HGNC:11604", "HGNC:4173"})

    assert [i.code for i in issues] == ["REF008"]
    assert "HGNC:4173" in issues[0].message


def test_reports_a_functional_record_with_an_unknown_gene_and_publication() -> None:
    functional = FunctionalEvidence.model_validate(
        {
            "id": "CHDA:FUN:0000009",
            "gene": "HGNC:99999",
            "organism": "NCBITaxon:10090",
            "perturbation": "knockout",
            "zygosity": "heterozygous",
            "cardiac_phenotype": ["MP:0010402"],
            "phenocopies_human": "partial",
            "rescue_outcome": "not_attempted",
            "publication": "PMID:999",
        }
    )
    corpus = _corpus(functional=(functional,))

    issues = validate_references(corpus, known_genes={"HGNC:11604"})

    assert sorted(i.code for i in issues) == ["REF001", "REF002"]


def test_reports_a_featured_manuscript_citing_an_unknown_publication() -> None:
    featured = FeaturedManuscript(publication="PMID:999", order=1, blurb="b", topic="genomics")
    issues = validate_references(_corpus(featured=(featured,)), known_genes={"HGNC:11604"})
    assert [i.code for i in issues] == ["REF002"]


def test_reports_evidence_citing_an_unknown_dataset() -> None:
    corpus = _corpus(assertions=(_assertion(evidence=[_evidence(dataset="PXD999999")]),))
    issues = validate_references(corpus, known_genes={"HGNC:11604"})
    assert [i.code for i in issues] == ["REF004"]


def _septal_term() -> PhenotypeTerm:
    return PhenotypeTerm.model_validate(
        {"id": "HP:0001631", "label": "Atrial septal defect", "lesion_group": "septal"}
    )


def test_phenotype_consistent_with_declared_lesion_group_passes() -> None:
    corpus = _corpus(phenotypes=(_septal_term(),))
    assert validate_references(corpus, known_genes={"HGNC:11604"}) == []


def test_reports_a_phenotype_contradicting_its_lesion_group() -> None:
    """REF009 joins it: septal is missing *and* nothing declared justifies conotruncal.

    Both statements are independently true of this record, and together they say
    what to do — drop conotruncal, add septal — where REF007 alone leaves the
    unjustified group standing.
    """
    corpus = _corpus(
        assertions=(_assertion(lesion_groups=["conotruncal"]),),
        phenotypes=(_septal_term(),),
    )
    issues = validate_references(corpus, known_genes={"HGNC:11604"})
    assert [i.code for i in issues] == ["REF007", "REF009"]
    assert "septal" in issues[0].message


def test_unmapped_phenotype_is_not_reported_as_inconsistent() -> None:
    """A phenotype absent from curation/phenotypes.yaml has no group to contradict."""
    corpus = _corpus(phenotypes=())
    assert validate_references(corpus, known_genes={"HGNC:11604"}) == []


def test_reports_an_unjustified_lesion_group() -> None:
    corpus = _corpus(
        assertions=(_assertion(lesion_groups=["septal", "conotruncal"]),),
        phenotypes=(_septal_term(),),
    )

    issues = validate_references(corpus, known_genes={"HGNC:11604"})

    assert [i.code for i in issues] == ["REF009"]
    assert "conotruncal" in issues[0].message


def test_unmapped_phenotype_exempts_the_unjustified_group_check_but_says_so() -> None:
    """The exemption stands; what changes is that it stops being silent.

    An unmapped phenotype could legitimately justify the extra group, so REF009
    must not fire. Left at that, one term missing from the register disables all
    lesion-group checking for the assertion with nothing in the report to say a
    check stopped running. REF012 names the term that caused it.
    """
    corpus = _corpus(
        assertions=(
            _assertion(
                phenotypes=["HP:0001631", "HP:0001629"],
                lesion_groups=["septal", "conotruncal"],
            ),
        ),
        phenotypes=(_septal_term(),),
    )

    issues = validate_references(corpus, known_genes={"HGNC:11604"})

    assert [i.code for i in issues] == ["REF012"]
    assert "HP:0001629" in issues[0].message


def test_an_empty_phenotype_register_reports_nothing_per_record() -> None:
    """CUR002 names an absent register once; REF012 would repeat it per phenotype.

    This is the cascade guard the rest of this module is built around: a check
    that cannot run must be reported once, at its cause, not once per record that
    could not be checked.
    """
    corpus = _corpus(
        assertions=(_assertion(phenotypes=["HP:0001631", "HP:0001629"]),),
        phenotypes=(),
    )

    assert validate_references(corpus, known_genes={"HGNC:11604"}) == []


def test_reports_a_cardiac_lesion_listed_as_an_extracardiac_feature() -> None:
    corpus = _corpus(
        assertions=(_assertion(syndromic="syndromic", extracardiac_features=["HP:0001631"]),),
        phenotypes=(_septal_term(),),
    )

    issues = validate_references(corpus, known_genes={"HGNC:11604"})

    assert [i.code for i in issues] == ["REF010"]


def _dataset(accession: str = "PXD012345", **overrides: object) -> Dataset:
    payload: dict[str, object] = {
        "id": accession,
        "archive": "pride",
        "technology": "phosphoproteomics",
        "tissue": "right ventricle",
        "developmental_stage": "infant",
        "organism": "NCBITaxon:9606",
        "n_samples": 24,
        "licence": "CC-BY-4.0",
        "contrasts": [
            {
                "id": "tof_vs_control",
                "description": "TOF right ventricle versus control",
                "case_group": "TOF RV",
                "control_group": "Control RV",
                "statistical_method": "limma-voom",
                "software": "limma 3.58.1",
            }
        ],
    }
    payload.update(overrides)
    return Dataset.model_validate(payload)


def test_reports_the_same_accession_declared_by_two_dataset_files() -> None:
    """Each dataset is its own file, so nothing else can see the collision.

    `validate_ids` checks uniqueness only for atlas-minted IDs; an accession is
    a third-party identifier and never reaches it. Two files declaring PXD012345
    with different sample counts both load, and every mirror row citing that
    accession then resolves to whichever the filesystem happened to yield first.
    """
    corpus = Corpus(
        root=Path("/repo"),
        datasets=(
            _dataset(n_samples=24),
            _dataset(n_samples=48, tissue="left ventricle"),
        ),
    )

    issues = validate_references(corpus, known_genes=None)

    assert [i.code for i in issues] == ["REF011"]
    assert issues[0].location == "dataset PXD012345"


def test_accepts_distinct_dataset_accessions() -> None:
    corpus = Corpus(
        root=Path("/repo"),
        datasets=(_dataset("PXD012345"), _dataset("PXD999999")),
    )

    assert validate_references(corpus, known_genes=None) == []


DATASET_YAML = """\
id: PXD012345
archive: pride
technology: phosphoproteomics
tissue: right ventricle
developmental_stage: infant
organism: NCBITaxon:9606
n_samples: 24
licence: CC-BY-4.0
contrasts:
  - id: tof_vs_control
    description: TOF right ventricle versus control
    case_group: TOF RV
    control_group: Control RV
    statistical_method: limma-voom
    software: limma 3.58.1
"""

EXPRESSION_TSV_HEADER = (
    "dataset\tcontrast\tgene\tlog2fc\tpvalue\tfdr\tdirection\tn_case\tn_control\ttissue\tstage\n"
)

# profiles is the one dataset-linked schema with no contrast column.
PROFILES_TSV_HEADER = "dataset\tgene\ttissue\tstage\tmedian_abundance\tunit\tq25\tq75\tn_samples\n"
VALID_EXPRESSION_ROW = (
    "PXD012345\ttof_vs_control\tHGNC:11604\t1.2\t0.001\t0.01\tup\t10\t10\tRV\tinfant\n"
)


def _mirror_repo(tmp_path: Path, row: str) -> Path:
    from chd_atlas.corpus import load_curation

    (tmp_path / "curation" / "datasets").mkdir(parents=True)
    (tmp_path / "curation" / "datasets" / "PXD012345.yaml").write_text(DATASET_YAML)
    (tmp_path / "curation" / "publications.yaml").write_text(
        "publications:\n"
        "  - id: PMID:8988165\n"
        "    title: T\n"
        "    journal: J\n"
        "    year: 1997\n"
        "    authors: [A]\n"
        "    study_type: case_control\n"
    )
    (tmp_path / "mirrors" / "expression").mkdir(parents=True)
    (tmp_path / "mirrors" / "expression" / "PXD012345.tsv").write_text(EXPRESSION_TSV_HEADER + row)
    load_curation(tmp_path)
    return tmp_path


def _profiles_repo(tmp_path: Path, row: str) -> Path:
    """A repository whose expression shard is clean, so only profiles can fail."""
    root = _mirror_repo(tmp_path, VALID_EXPRESSION_ROW)
    (root / "mirrors" / "profiles").mkdir(parents=True)
    (root / "mirrors" / "profiles" / "PXD012345.tsv").write_text(PROFILES_TSV_HEADER + row)
    return root


def test_mirror_rows_referencing_a_known_contrast_pass(tmp_path: Path) -> None:
    from chd_atlas.corpus import load_curation

    root = _mirror_repo(tmp_path, VALID_EXPRESSION_ROW)
    corpus, _ = load_curation(root)

    assert validate_mirror_references(root, corpus) == []


def test_reports_a_mirror_row_referencing_an_unknown_dataset(tmp_path: Path) -> None:
    from chd_atlas.corpus import load_curation

    root = _mirror_repo(
        tmp_path,
        "PXD999999\ttof_vs_control\tHGNC:11604\t1.2\t0.001\t0.01\tup\t10\t10\tRV\tinfant\n",
    )
    corpus, _ = load_curation(root)

    issues = validate_mirror_references(root, corpus)

    assert [i.code for i in issues] == ["REF005"]


def test_reports_a_mirror_row_referencing_an_unknown_contrast(tmp_path: Path) -> None:
    from chd_atlas.corpus import load_curation

    root = _mirror_repo(
        tmp_path,
        "PXD012345\tmystery_contrast\tHGNC:11604\t1.2\t0.001\t0.01\tup\t10\t10\tRV\tinfant\n",
    )
    corpus, _ = load_curation(root)

    issues = validate_mirror_references(root, corpus)

    assert [i.code for i in issues] == ["REF006"]
    assert "mystery_contrast" in issues[0].message


def test_a_profiles_shard_is_checked_although_it_has_no_contrast_column(
    tmp_path: Path,
) -> None:
    """profiles rows are two columns wide here; indexing row[1] would raise."""
    from chd_atlas.corpus import load_curation

    root = _profiles_repo(tmp_path, "PXD012345\tHGNC:11604\tRV\tinfant\t12.5\ttpm\t8.0\t20.0\t24\n")
    corpus, _ = load_curation(root)

    assert validate_mirror_references(root, corpus) == []


def test_reports_a_profiles_row_referencing_an_unknown_dataset(tmp_path: Path) -> None:
    from chd_atlas.corpus import load_curation

    root = _profiles_repo(tmp_path, "PXD999999\tHGNC:11604\tRV\tinfant\t12.5\ttpm\t8.0\t20.0\t24\n")
    corpus, _ = load_curation(root)

    issues = validate_mirror_references(root, corpus)

    assert [i.code for i in issues] == ["REF005"]
    assert "PXD999999" in issues[0].message


def test_a_null_dataset_cell_is_not_reported_as_a_missing_dataset(tmp_path: Path) -> None:
    """A null dataset is TBL003's job; REF005 would name a dataset called None."""
    from chd_atlas.corpus import load_curation

    root = _mirror_repo(
        tmp_path,
        "\ttof_vs_control\tHGNC:11604\t1.2\t0.001\t0.01\tup\t10\t10\tRV\tinfant\n",
    )
    corpus, _ = load_curation(root)

    assert validate_mirror_references(root, corpus) == []


def test_a_null_contrast_cell_is_not_reported_as_a_missing_contrast(tmp_path: Path) -> None:
    """Likewise TBL003's job; REF006 would quote a contrast called None."""
    from chd_atlas.corpus import load_curation

    root = _mirror_repo(
        tmp_path,
        "PXD012345\t\tHGNC:11604\t1.2\t0.001\t0.01\tup\t10\t10\tRV\tinfant\n",
    )
    corpus, _ = load_curation(root)

    assert validate_mirror_references(root, corpus) == []


def _registry(root: Path, uniprot: str) -> None:
    """A one-row `mirrors/genes.tsv` for HGNC:11604 with the given accession."""
    columns = [column.name for column in TABLE_SCHEMAS["genes"].columns]
    cells = dict.fromkeys(columns, "")
    cells["hgnc_id"] = "HGNC:11604"
    cells["symbol"] = "TBX5"
    cells["name"] = "T-box transcription factor 5"
    cells["uniprot"] = uniprot
    mirrors = root / "mirrors"
    mirrors.mkdir(parents=True, exist_ok=True)
    (mirrors / "genes.tsv").write_text(
        "\t".join(columns) + "\n" + "\t".join(cells[c] for c in columns) + "\n"
    )


def test_ptm_evidence_with_no_accession_to_join_by_is_reported(tmp_path: Path) -> None:
    """The only route from a phospho row to a gene is one nullable mirror cell.

    `mirrors/phospho/` carries no gene column at all, so `genes.tsv:uniprot` is
    the whole join. Blank, and the gene's PTM sites are published in a shard and
    summarised in no bundle: `"omics": {}` on the page while the browse row
    beside it advertises the PTM evidence the assertion claims. Reproduced end to
    end on a build reporting 0 errors and 0 warnings before this rule existed,
    which is what makes it the project's characteristic failure rather than an
    inconvenience.

    Both rows matter. With the accession present there must be no issue at all —
    otherwise the rule fires on every correctly curated gene and a curator learns
    to ignore it, which is worse than not reporting.
    """
    ptm = _evidence(evidence_class="ptm", dataset="PXD012345")
    corpus = _corpus(assertions=(_assertion(evidence=[ptm]),))

    _registry(tmp_path, uniprot="")
    issues = validate_ptm_evidence_is_reachable(tmp_path, corpus)
    assert [issue.code for issue in issues] == ["REF013"]
    assert issues[0].severity is Severity.WARNING
    # The gene and the assertion, because the curator has to find both to fix it.
    assert "HGNC:11604" in issues[0].message
    assert "CHDA:AST:0000001" in issues[0].message

    _registry(tmp_path, uniprot="Q99593")
    assert validate_ptm_evidence_is_reachable(tmp_path, corpus) == []


def test_a_gene_asserting_no_ptm_evidence_is_never_asked_for_an_accession(
    tmp_path: Path,
) -> None:
    """`uniprot` is nullable because most genes never need it.

    Only PTM evidence depends on it: `proteomics` rows carry their own `gene`
    column, so an accession is one of two routes there rather than the only one,
    and `expression`/`profiles` never use it. A rule that asked every gene for an
    accession would report most of the registry, which is the cascade REF000 and
    SRC000 exist to prevent.
    """
    _registry(tmp_path, uniprot="")

    others = [
        _evidence(evidence_class="genetic_case"),
        # Dataset-backed and accession-adjacent, and still not covered: a
        # `proteomics` row carries its own `gene` column, so an accession is one
        # of two routes there rather than the only one.
        _evidence(evidence_class="proteomic", dataset="PXD012345"),
        _evidence(evidence_class="expression", dataset="GSE1000"),
    ]
    for evidence in others:
        corpus = _corpus(assertions=(_assertion(evidence=[evidence]),))
        assert validate_ptm_evidence_is_reachable(tmp_path, corpus) == [], evidence.evidence_class
