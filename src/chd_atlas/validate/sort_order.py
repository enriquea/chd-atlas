# src/chd_atlas/validate/sort_order.py
"""Enforces each table's declared sort key.

A stable sort order is what keeps agent-inserted rows producing small, local
diffs instead of whole-file reshuffles, which is what makes the review digest
in plan 4 meaningful.
"""

from __future__ import annotations

from pathlib import Path

from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.tables import TableSchema, read_table


def validate_sort_order(path: Path, schema: TableSchema) -> list[ValidationIssue]:
    frame, _ = read_table(path, schema)
    if frame is None:
        # validate_table already reported the read failure; don't duplicate it.
        return []
    if any(column not in frame.columns for column in schema.sort_key):
        # Missing columns are reported by validate_table; nothing to check here.
        return []
    if frame.height < 2:
        return []

    keys = frame.select(schema.sort_key).rows()
    for index in range(1, len(keys)):
        if _precedes(keys[index], keys[index - 1]):
            return [
                ValidationIssue(
                    code="SORT001",
                    severity=Severity.ERROR,
                    location=f"{path}:row {index + 2}",
                    message=(
                        f"table is not sorted by {list(schema.sort_key)}: "
                        f"{keys[index]} at row {index + 2} follows "
                        f"{keys[index - 1]} at row {index + 1}"
                    ),
                )
            ]
    return []


def _precedes(left: tuple[object, ...], right: tuple[object, ...]) -> bool:
    """True when ``left`` sorts before ``right``, treating null as smallest."""
    for a, b in zip(left, right, strict=True):
        if a == b:
            continue
        if a is None:
            return True
        if b is None:
            return False
        return bool(a < b)  # type: ignore[operator]
    return False
