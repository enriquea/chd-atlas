from pathlib import Path

import pytest
from ruamel.yaml.representer import RepresenterError

from chd_atlas.issues import Severity
from chd_atlas.validate.ids import (
    allocate,
    load_id_registry,
    save_id_registry,
    validate_ids,
)

REGISTRY_YAML = "prefixes:\n  AST: 3\n  FUN: 1\n"


def _write_registry(root: Path, text: str = REGISTRY_YAML) -> None:
    (root / "curation").mkdir(parents=True, exist_ok=True)
    (root / "curation" / ".id_registry.yaml").write_text(text)


def test_loads_the_counter_file(tmp_path: Path) -> None:
    _write_registry(tmp_path)
    registry, issues = load_id_registry(tmp_path)
    assert issues == []
    assert registry is not None
    assert registry.prefixes == {"AST": 3, "FUN": 1}


def test_missing_registry_is_an_error(tmp_path: Path) -> None:
    (tmp_path / "curation").mkdir()
    registry, issues = load_id_registry(tmp_path)
    assert registry is None
    assert [i.code for i in issues] == ["ID001"]


def test_allocate_returns_the_next_id_and_advances_the_counter(tmp_path: Path) -> None:
    _write_registry(tmp_path)
    registry, _ = load_id_registry(tmp_path)
    assert registry is not None
    assert allocate(registry, "AST") == "CHDA:AST:0000004"
    assert registry.prefixes["AST"] == 4


def test_allocate_starts_a_new_prefix_at_one(tmp_path: Path) -> None:
    _write_registry(tmp_path, "prefixes: {}\n")
    registry, _ = load_id_registry(tmp_path)
    assert registry is not None
    assert allocate(registry, "AST") == "CHDA:AST:0000001"


def test_accepts_ids_within_the_counter(tmp_path: Path) -> None:
    _write_registry(tmp_path)
    registry, _ = load_id_registry(tmp_path)
    assert registry is not None
    assert validate_ids(["CHDA:AST:0000001", "CHDA:AST:0000003"], registry) == []


def test_reports_an_id_above_the_counter(tmp_path: Path) -> None:
    _write_registry(tmp_path)
    registry, _ = load_id_registry(tmp_path)
    assert registry is not None
    issues = validate_ids(["CHDA:AST:0000009"], registry)
    assert [i.code for i in issues] == ["ID002"]


def test_reports_a_duplicate_id(tmp_path: Path) -> None:
    _write_registry(tmp_path)
    registry, _ = load_id_registry(tmp_path)
    assert registry is not None
    issues = validate_ids(["CHDA:AST:0000001", "CHDA:AST:0000001"], registry)
    assert [i.code for i in issues] == ["ID003"]


def test_gaps_are_allowed_because_deleted_ids_are_retired(tmp_path: Path) -> None:
    _write_registry(tmp_path)
    registry, _ = load_id_registry(tmp_path)
    assert registry is not None
    assert validate_ids(["CHDA:AST:0000003"], registry) == []


@pytest.mark.parametrize(
    "text",
    [
        "just-a-scalar\n",
        "- a\n- b\n",
        "prefixes: not-a-mapping\n",
        "prefixes: [1, 2]\n",
        "prefixes:\n  AST: three\n",
        "prefixes:\n  AST:\n",
        "prefixes:\n  AST: [1]\n",
        "prefixes:\n  AST: 3.9\n",
        "prefixes:\n  AST: true\n",
        "prefixes:\n  AST: -5\n",
        "",
    ],
)
def test_malformed_registry_is_reported_rather_than_raised(tmp_path: Path, text: str) -> None:
    _write_registry(tmp_path, text)

    registry, issues = load_id_registry(tmp_path)

    assert registry is None
    assert issues
    assert all(i.severity is Severity.ERROR for i in issues)


def test_out_of_range_date_counter_is_reported(tmp_path: Path) -> None:
    """ruamel builds a datetime.date itself and raises a bare ValueError."""
    _write_registry(tmp_path, "prefixes:\n  AST: 2026-13-45\n")

    registry, issues = load_id_registry(tmp_path)

    assert registry is None
    assert [i.code for i in issues] == ["YAML001"]


def test_absent_prefixes_key_is_a_valid_empty_registry(tmp_path: Path) -> None:
    _write_registry(tmp_path, "{}\n")

    registry, issues = load_id_registry(tmp_path)

    assert issues == []
    assert registry is not None
    assert registry.prefixes == {}


@pytest.mark.parametrize(
    "identifier",
    [
        "CHDA:AST",
        "garbage",
        "CHDA:AST:0000001:extra",
        ":::",
        "CHDA:AST:abc",
        "CHDA:AST:0000001 ",
        "CHDA:AST:-1",
        "CHDA:AST:0000000",
    ],
)
def test_malformed_identifier_is_reported_rather_than_raised(
    tmp_path: Path, identifier: str
) -> None:
    _write_registry(tmp_path)
    registry, _ = load_id_registry(tmp_path)
    assert registry is not None

    issues = validate_ids([identifier], registry)

    assert [i.code for i in issues] == ["ID005"]


def test_save_round_trips_and_is_byte_stable(tmp_path: Path) -> None:
    _write_registry(tmp_path)
    registry, _ = load_id_registry(tmp_path)
    assert registry is not None
    allocate(registry, "AST")
    allocate(registry, "ZZZ")

    save_id_registry(tmp_path, registry)
    first = (tmp_path / "curation" / ".id_registry.yaml").read_bytes()
    reloaded, issues = load_id_registry(tmp_path)
    assert issues == []
    assert reloaded is not None
    assert reloaded.prefixes == {"AST": 4, "FUN": 1, "ZZZ": 1}

    save_id_registry(tmp_path, reloaded)
    assert (tmp_path / "curation" / ".id_registry.yaml").read_bytes() == first


def test_saved_registry_is_world_readable(tmp_path: Path) -> None:
    """A committed artifact; mkstemp's 0600 would leave it owner-only."""
    _write_registry(tmp_path)
    registry, _ = load_id_registry(tmp_path)
    assert registry is not None

    save_id_registry(tmp_path, registry)

    path = tmp_path / "curation" / ".id_registry.yaml"
    assert path.stat().st_mode & 0o777 == 0o644


def test_save_leaves_the_old_counter_intact_if_the_write_fails(tmp_path: Path) -> None:
    """A truncate-then-fail would leave a zero-byte counter, which is data loss."""
    _write_registry(tmp_path)
    path = tmp_path / "curation" / ".id_registry.yaml"
    before = path.read_bytes()

    class Unserialisable:
        pass

    registry, _ = load_id_registry(tmp_path)
    assert registry is not None
    registry.prefixes["BAD"] = Unserialisable()  # type: ignore[assignment]

    # Narrowed from a blind `Exception` because ruff's B017 rejects that, and
    # naming the error proves the dump is what failed rather than the setup.
    with pytest.raises(RepresenterError):
        save_id_registry(tmp_path, registry)

    assert path.read_bytes() == before
