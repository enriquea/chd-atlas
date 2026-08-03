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

    # `iter_rows` rather than `rows`, which builds the whole list of tuples
    # before the first comparison. Measured on a 200,000-row `variants` table:
    # that list alone is 36 MB, on top of the frame `read_table` already holds,
    # and it is discarded after the first out-of-order pair. Variants are the
    # table this project says grows without bound, so the copy is paid for by
    # the largest table and read by a check that usually stops early.
    #
    # Only ever two rows are live here, which is all a sortedness check needs.
    previous: tuple[object, ...] | None = None
    for index, current in enumerate(frame.select(schema.sort_key).iter_rows()):
        if previous is not None and _precedes(current, previous):
            # +2: one for the header line, one because `index` is 0-based, so
            # the number matches what a curator sees in an editor.
            return [
                ValidationIssue(
                    code="SORT001",
                    severity=Severity.ERROR,
                    location=f"{path}:row {index + 2}",
                    message=(
                        f"table is not sorted by {list(schema.sort_key)}: "
                        f"{current} at row {index + 2} follows "
                        f"{previous} at row {index + 1}"
                    ),
                )
            ]
        previous = current
    return []


def _precedes(left: tuple[object, ...], right: tuple[object, ...]) -> bool:
    """True when ``left`` sorts before ``right``, treating null as smallest.

    `build/omics.py`'s `_rank_key` implements the same total order for the same
    reason and is deliberately not shared with this: it must never raise, since
    a ranking that dies mid-build publishes nothing, while this must raise on a
    column pair it cannot compare, because a sort order it cannot verify is one
    it must not call clean. The two must agree on the null rule regardless, so a
    reader changing either should read the other — noted here as well as there,
    since a change is as likely to start from this side.
    """
    for a, b in zip(left, right, strict=True):
        if a == b:
            continue
        if a is None:
            return True
        if b is None:
            return False
        return bool(a < b)  # type: ignore[operator]
    return False
