"""Atlas identifier allocation and consistency.

IDs are allocated monotonically from a committed counter. A merge conflict on
the counter file is the intended signal that two branches allocated
concurrently and one must renumber before merge. Deleting a record retires its
ID permanently: the counter never rewinds.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from io import BytesIO
from pathlib import Path

from ruamel.yaml import YAML
from ruamel.yaml.error import YAMLError

from chd_atlas.fs import write_bytes_atomically
from chd_atlas.issues import Severity, ValidationIssue

_ID_WIDTH = 7


@dataclass
class IdRegistry:
    prefixes: dict[str, int] = field(default_factory=dict)


def load_id_registry(root: Path) -> tuple[IdRegistry | None, list[ValidationIssue]]:
    """Read the committed counter file.

    Returns `None` rather than an empty registry when the file cannot be read,
    because the two are not interchangeable: an empty registry has every
    ceiling at zero, so a caller that treated a failed load as "nothing
    allocated yet" would report every atlas ID in the corpus as exceeding its
    counter and name the real cause nowhere.
    """
    path = root / "curation" / ".id_registry.yaml"
    if not path.is_file():
        return None, [ValidationIssue("ID001", Severity.ERROR, str(path), "ID registry not found")]

    yaml = YAML(typ="safe")
    try:
        raw = yaml.load(path.read_text(encoding="utf-8"))
        return _registry_from(str(path), raw)
    # ruamel builds a `datetime.date` directly and raises a bare ValueError on
    # an out-of-range literal such as `2026-13-45`, which YAMLError does not
    # cover. UnicodeDecodeError is already a ValueError subclass; it is kept
    # listed to document that non-UTF-8 bytes are handled here too. The shape
    # checks are inside the guard as well so that no reader defect can escape
    # this function as an exception.
    # OSError covers what read_text can raise before parsing: an unreadable
    # mode, a dangling symlink, a directory where a file was expected. One
    # bad file must be one issue, not a traceback aborting the whole run.
    except (YAMLError, UnicodeDecodeError, ValueError, OSError) as exc:
        return None, [
            ValidationIssue("YAML001", Severity.ERROR, str(path), f"could not read YAML: {exc}")
        ]


def _registry_from(location: str, raw: object) -> tuple[IdRegistry | None, list[ValidationIssue]]:
    """Validate the parsed document's shape before trusting it as a counter.

    Every rejection returns `None`. Coercing instead would rewind the counter:
    `int(3.9)` is 3 and `int(True)` is 1, either of which reissues an ID that
    was already allocated — the one thing this module promises never happens.
    """
    if raw is None:
        # An empty or whitespace-only file parses cleanly as `None`, so it is a
        # shape problem rather than a parse failure.
        return None, [
            ValidationIssue(
                "ID004",
                Severity.ERROR,
                location,
                "ID registry is empty; expected a mapping",
            )
        ]

    if not isinstance(raw, dict):
        return None, [
            ValidationIssue(
                "ID004",
                Severity.ERROR,
                location,
                f"ID registry must be a mapping, found {type(raw).__name__}",
            )
        ]

    if "prefixes" not in raw:
        # The legitimate "nothing allocated yet" state. A present-but-null
        # `prefixes` is not this state: it is a mapping key with the wrong
        # type, and falls through to the check below.
        return IdRegistry(), []

    counters = raw["prefixes"]
    if not isinstance(counters, dict):
        return None, [
            ValidationIssue(
                "ID004",
                Severity.ERROR,
                location,
                f"'prefixes' must be a mapping, found {type(counters).__name__}",
            )
        ]

    issues: list[ValidationIssue] = []
    prefixes: dict[str, int] = {}
    for key, value in counters.items():
        prefix = str(key)
        # `isinstance(True, int)` is True in Python, so booleans need an
        # explicit exclusion or `AST: true` loads as the counter 1.
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            issues.append(
                ValidationIssue(
                    "ID004",
                    Severity.ERROR,
                    location,
                    f"counter for '{prefix}' must be a non-negative integer, found {value!r}",
                )
            )
            continue
        prefixes[prefix] = value

    if issues:
        return None, issues
    return IdRegistry(prefixes=prefixes), []


def allocate(registry: IdRegistry, prefix: str) -> str:
    """Return the next unused ID for ``prefix`` and advance the counter."""
    next_ordinal = registry.prefixes.get(prefix, 0) + 1
    registry.prefixes[prefix] = next_ordinal
    return f"CHDA:{prefix}:{next_ordinal:0{_ID_WIDTH}d}"


def save_id_registry(root: Path, registry: IdRegistry) -> None:
    """Write the counter file atomically.

    A counter file that loses its contents is a silent rewind of every prefix to
    zero, after which `allocate` reissues identifiers that are already in use.
    `write_bytes_atomically` is what prevents a partial write from doing that;
    the mechanism lives there rather than being restated here.

    The dump goes through an in-memory stream rather than straight at the
    destination on purpose. Serialising first means a value ruamel cannot
    represent raises before any file exists — no temporary to clean up, and the
    committed counter never opened. Handing the destination to `yaml.dump`
    directly would be shorter and would reintroduce exactly the truncation this
    guards against.
    """
    path = root / "curation" / ".id_registry.yaml"
    yaml = YAML()
    yaml.default_flow_style = False

    # BytesIO rather than StringIO so the encoding stays ruamel's decision: the
    # YAML instance already carries `encoding = "utf-8"` and applies it here.
    stream = BytesIO()
    yaml.dump({"prefixes": dict(sorted(registry.prefixes.items()))}, stream)
    write_bytes_atomically(path, stream.getvalue())


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
        # Unpacking blindly raises on anything that is not exactly three
        # colon-separated fields, and `int()` raises on a non-numeric ordinal.
        # The digit test also rejects a trailing space, which would otherwise
        # validate clean and defeat the duplicate check above, since ID003 keys
        # on the raw string. The `< 1` clause rejects ordinal zero, which
        # `allocate` starting at 1 means was never issued.
        #
        # `isdecimal()` rather than `isdigit()`: the latter is True for
        # characters such as the superscript "²" that `int()` then refuses,
        # so it would leave exactly the crash this check exists to prevent.
        parts = identifier.split(":")
        if len(parts) != 3 or not parts[2].isdecimal() or int(parts[2]) < 1:
            issues.append(
                ValidationIssue(
                    "ID005",
                    Severity.ERROR,
                    identifier,
                    "malformed atlas identifier",
                )
            )
            continue

        _, prefix, ordinal = parts
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
