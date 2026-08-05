# scripts/burden_mirror.py
"""Shared writer for `mirrors/burden.tsv`, which now holds more than one study.

`convert_burden.py` wrote the whole file. With a second study that stops being
safe: whichever converter ran last would silently delete the other's rows, and
the loss would look exactly like a successful run.

`merge_study` instead replaces only the rows carrying its own `study` value and
leaves every other study's rows alone. That makes each converter idempotent and
order-independent -- a curator can re-run either, in either order, and get the
same file -- which is the property that stops the mirror depending on the
sequence somebody happened to type.

Imported by the per-study converters. Kept out of `src/` because it is curator
tooling: design decision D6 says CI validates, builds and deploys, and never
authors data, so nothing the build imports may depend on this.
"""

from __future__ import annotations

from pathlib import Path

# The order `TABLE_SCHEMAS["burden"]` declares. Written out rather than imported
# so the converters stay runnable without the package installed, matching
# `convert_hgnc.py` and `convert_burden.py`.
ORDER: tuple[str, ...] = (
    "study",
    "gene",
    "cohort_stratum",
    "lesion_group",
    "variant_class",
    "consequence_class",
    "origin",
    "maf_max",
    "count_unit",
    "n_case_carriers",
    "n_cases",
    "comparator",
    "n_control_carriers",
    "n_controls",
    "expected_count",
    "effect",
    "effect_measure",
    "effect_bound",
    "ci_low",
    "ci_high",
    "pvalue",
    "pvalue_test",
    "pvalue_adjusted",
    "pvalue_adjustment",
    "case_cohorts",
    "control_cohorts",
    "method_note",
    "source",
)

# `TABLE_SCHEMAS["burden"].sort_key`, which `validate_sort_order` enforces.
_SORT_KEY: tuple[str, ...] = (
    "gene",
    "study",
    "cohort_stratum",
    "lesion_group",
    "variant_class",
    "consequence_class",
    "origin",
    "maf_max",
)


def _sort_value(row: dict[str, str], column: str) -> object:
    """One sort-key cell, compared the way `validate_sort_order` compares it.

    That module reads the values polars parsed, so `maf_max` is a float there
    and an empty cell is null. Sorting `maf_max` as text would disagree the
    moment a second frequency threshold appears -- "0.001" sorts after "0.0001"
    as a string and before it as a number -- and an empty cell must sort first,
    matching `_precedes`' null-is-smallest rule.
    """
    value = row[column]
    if column == "maf_max":
        # (present?, number) so an empty cell sorts before every real one.
        return (1, float(value)) if value else (0, 0.0)
    return (1, value) if value else (0, "")


def merge_study(out: Path, study: str, rows: list[dict[str, str]]) -> tuple[int, int]:
    """Replace `study`'s rows in `out` with `rows`. Returns (kept, written).

    `kept` counts rows belonging to other studies that survived untouched --
    printed by every converter, because "0 kept" on a mirror that should hold a
    second study is the one number that reveals a merge which silently became an
    overwrite.
    """
    kept: list[dict[str, str]] = []
    if out.is_file():
        lines = out.read_text(encoding="utf-8").splitlines()
        header = lines[0].split("\t")
        for line in lines[1:]:
            if not line:
                continue
            existing = dict(zip(header, line.split("\t"), strict=False))
            # A row whose file predates a column simply has no value for it;
            # `.get` keeps the merge working across a schema addition rather
            # than raising and leaving the mirror half-written.
            if existing.get("study") != study:
                kept.append({name: existing.get(name, "") for name in ORDER})

    merged = sorted(
        kept + rows, key=lambda row: tuple(_sort_value(row, name) for name in _SORT_KEY)
    )
    body = "\n".join("\t".join(row[name] for name in ORDER) for row in merged)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\t".join(ORDER) + "\n" + body + "\n", encoding="utf-8")
    return len(kept), len(rows)
