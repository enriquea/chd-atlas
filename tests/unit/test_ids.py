from pathlib import Path

from chd_atlas.validate.ids import allocate, load_id_registry, validate_ids

REGISTRY_YAML = "prefixes:\n  AST: 3\n  FUN: 1\n"


def _write_registry(root: Path, text: str = REGISTRY_YAML) -> None:
    (root / "curation").mkdir(parents=True, exist_ok=True)
    (root / "curation" / ".id_registry.yaml").write_text(text)


def test_loads_the_counter_file(tmp_path: Path) -> None:
    _write_registry(tmp_path)
    registry, issues = load_id_registry(tmp_path)
    assert issues == []
    assert registry.prefixes == {"AST": 3, "FUN": 1}


def test_missing_registry_is_an_error(tmp_path: Path) -> None:
    (tmp_path / "curation").mkdir()
    _, issues = load_id_registry(tmp_path)
    assert [i.code for i in issues] == ["ID001"]


def test_allocate_returns_the_next_id_and_advances_the_counter(tmp_path: Path) -> None:
    _write_registry(tmp_path)
    registry, _ = load_id_registry(tmp_path)
    assert allocate(registry, "AST") == "CHDA:AST:0000004"
    assert registry.prefixes["AST"] == 4


def test_allocate_starts_a_new_prefix_at_one(tmp_path: Path) -> None:
    _write_registry(tmp_path, "prefixes: {}\n")
    registry, _ = load_id_registry(tmp_path)
    assert allocate(registry, "AST") == "CHDA:AST:0000001"


def test_accepts_ids_within_the_counter(tmp_path: Path) -> None:
    _write_registry(tmp_path)
    registry, _ = load_id_registry(tmp_path)
    assert validate_ids(["CHDA:AST:0000001", "CHDA:AST:0000003"], registry) == []


def test_reports_an_id_above_the_counter(tmp_path: Path) -> None:
    _write_registry(tmp_path)
    registry, _ = load_id_registry(tmp_path)
    issues = validate_ids(["CHDA:AST:0000009"], registry)
    assert [i.code for i in issues] == ["ID002"]


def test_reports_a_duplicate_id(tmp_path: Path) -> None:
    _write_registry(tmp_path)
    registry, _ = load_id_registry(tmp_path)
    issues = validate_ids(["CHDA:AST:0000001", "CHDA:AST:0000001"], registry)
    assert [i.code for i in issues] == ["ID003"]


def test_gaps_are_allowed_because_deleted_ids_are_retired(tmp_path: Path) -> None:
    _write_registry(tmp_path)
    registry, _ = load_id_registry(tmp_path)
    assert validate_ids(["CHDA:AST:0000003"], registry) == []
