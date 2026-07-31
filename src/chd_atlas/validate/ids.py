"""Atlas identifier allocation and consistency.

IDs are allocated monotonically from a committed counter. A merge conflict on
the counter file is the intended signal that two branches allocated
concurrently and one must renumber before merge. Deleting a record retires its
ID permanently: the counter never rewinds.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from ruamel.yaml import YAML
from ruamel.yaml.error import YAMLError

from chd_atlas.issues import Severity, ValidationIssue

_ID_WIDTH = 7


@dataclass
class IdRegistry:
    prefixes: dict[str, int] = field(default_factory=dict)


def load_id_registry(root: Path) -> tuple[IdRegistry, list[ValidationIssue]]:
    path = root / "curation" / ".id_registry.yaml"
    if not path.is_file():
        return IdRegistry(), [
            ValidationIssue("ID001", Severity.ERROR, str(path), "ID registry not found")
        ]

    yaml = YAML(typ="safe")
    try:
        raw = yaml.load(path.read_text(encoding="utf-8")) or {}
    except (YAMLError, UnicodeDecodeError) as exc:
        return IdRegistry(), [
            ValidationIssue("YAML001", Severity.ERROR, str(path), f"could not read YAML: {exc}")
        ]

    prefixes = {str(k): int(v) for k, v in (raw.get("prefixes") or {}).items()}
    return IdRegistry(prefixes=prefixes), []


def allocate(registry: IdRegistry, prefix: str) -> str:
    """Return the next unused ID for ``prefix`` and advance the counter."""
    next_ordinal = registry.prefixes.get(prefix, 0) + 1
    registry.prefixes[prefix] = next_ordinal
    return f"CHDA:{prefix}:{next_ordinal:0{_ID_WIDTH}d}"


def save_id_registry(root: Path, registry: IdRegistry) -> None:
    path = root / "curation" / ".id_registry.yaml"
    yaml = YAML()
    yaml.default_flow_style = False
    with path.open("w", encoding="utf-8") as handle:
        yaml.dump({"prefixes": dict(sorted(registry.prefixes.items()))}, handle)


def validate_ids(ids: list[str], registry: IdRegistry) -> list[ValidationIssue]:
    issues: list[ValidationIssue] = []

    for identifier, count in sorted(Counter(ids).items()):
        if count > 1:
            issues.append(
                ValidationIssue(
                    "ID003",
                    Severity.ERROR,
                    identifier,
                    f"identifier used {count} times; atlas IDs must be unique",
                )
            )

    for identifier in sorted(set(ids)):
        _, prefix, ordinal = identifier.split(":")
        ceiling = registry.prefixes.get(prefix, 0)
        if int(ordinal) > ceiling:
            issues.append(
                ValidationIssue(
                    "ID002",
                    Severity.ERROR,
                    identifier,
                    f"ordinal exceeds the allocated counter for '{prefix}' ({ceiling}); "
                    f"bump curation/.id_registry.yaml when allocating",
                )
            )

    return issues
