from pathlib import Path

from chd_atlas.validate.sources import load_sources

SOURCES_YAML = """\
sources:
  - id: clinvar
    name: ClinVar
    version: "2026-07"
    retrieved_on: 2026-07-15
    url: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/
    licence: public-domain
    redistribution: permitted
  - id: hpo
    name: Human Phenotype Ontology
    version: "2026-06-01"
    retrieved_on: 2026-07-15
    url: https://hpo.jax.org/
    licence: HPO-1.0
    redistribution: permitted
    ontology_prefix: HP
    ontology_file: ontologies/hp-2026-06-01.obo
"""


def _write_sources(root: Path, text: str = SOURCES_YAML) -> None:
    (root / "mirrors").mkdir(parents=True, exist_ok=True)
    (root / "mirrors" / "sources.yaml").write_text(text)


def test_loads_the_source_registry(tmp_path: Path) -> None:
    _write_sources(tmp_path)

    registry, issues = load_sources(tmp_path)

    assert issues == []
    assert registry.by_id["clinvar"].licence == "public-domain"
    assert registry.ontology_files() == {"HP": "ontologies/hp-2026-06-01.obo"}


def test_missing_licence_is_an_error(tmp_path: Path) -> None:
    _write_sources(tmp_path, SOURCES_YAML.replace("    licence: public-domain\n", ""))

    _, issues = load_sources(tmp_path)

    assert [i.code for i in issues] == ["SRC002"]
    assert "licence" in issues[0].message


def test_missing_registry_file_is_an_error(tmp_path: Path) -> None:
    (tmp_path / "mirrors").mkdir()

    _, issues = load_sources(tmp_path)

    assert [i.code for i in issues] == ["SRC001"]


def test_duplicate_source_ids_are_rejected(tmp_path: Path) -> None:
    _write_sources(tmp_path, SOURCES_YAML + SOURCES_YAML.split("sources:\n")[1])

    _, issues = load_sources(tmp_path)

    assert any(i.code == "SRC002" and "duplicate source" in i.message for i in issues)


def test_ontology_prefix_requires_a_pinned_file(tmp_path: Path) -> None:
    _write_sources(
        tmp_path,
        SOURCES_YAML.replace("    ontology_file: ontologies/hp-2026-06-01.obo\n", ""),
    )

    _, issues = load_sources(tmp_path)

    assert [i.code for i in issues] == ["SRC002"]
    assert "ontology_file" in issues[0].message


def test_source_referenced_by_a_record_must_have_a_registry_entry(tmp_path: Path) -> None:
    from chd_atlas.validate.sources import validate_source_references

    _write_sources(tmp_path)
    registry, _ = load_sources(tmp_path)

    issues = validate_source_references({"mirrors/variants/12.tsv": {"gnomad"}}, registry)

    assert [i.code for i in issues] == ["SRC003"]
    assert "licence unknown" in issues[0].message


def test_registered_sources_pass_the_reference_check(tmp_path: Path) -> None:
    from chd_atlas.validate.sources import validate_source_references

    _write_sources(tmp_path)
    registry, _ = load_sources(tmp_path)

    assert validate_source_references({"mirrors/variants/12.tsv": {"clinvar"}}, registry) == []
