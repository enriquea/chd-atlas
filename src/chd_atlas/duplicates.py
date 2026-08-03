# src/chd_atlas/duplicates.py
"""Duplicate detection shared by the record models."""

from __future__ import annotations

from collections import Counter
from collections.abc import Hashable, Iterable


def duplicates[T: Hashable](values: Iterable[T]) -> list[T]:
    """Every value appearing more than once, in first-seen order.

    Returns all of them rather than just the first: the corpus loader
    accumulates errors so one run reports every problem, and a per-file check
    that stopped at the first duplicate would make a curator fix them serially.
    """
    counts = Counter(values)
    return [value for value, count in counts.items() if count > 1]
