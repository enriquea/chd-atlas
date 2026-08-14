# src/chd_atlas/validate/genes.py
"""Symbol-to-id resolution, checked. Issue #33: report, never drop.

`GEN000` follows `REF000`/`SRC000`/`ONT000`/`SCP000`: when the registry cannot
be read, say "checks skipped" once instead of reporting every symbol in the
corpus as unresolvable. A missing reference table must never present as hundreds
of dangling references -- that buries the real error under cascading false
positives, and distinguishing "the data is wrong" from "the reference data never
loaded" is one of this project's two load-bearing guarantees.
"""

from __future__ import annotations

from collections.abc import Iterable

from chd_atlas.genes import GeneRegistry, Resolution
from chd_atlas.issues import Severity, ValidationIssue


def validate_symbol_resolution(
    symbols: Iterable[tuple[str, str]],
    registry: GeneRegistry | None,
    location: str,
) -> list[ValidationIssue]:
    """Every symbol a symbol-keyed source names must resolve to exactly one gene.

    `symbols` is `(symbol, where)` pairs -- the caller says where each came from,
    because "SYMBOL column, row 412" is what a curator needs and this function
    has no way to know it.

    **GEN001 (unmapped) and GEN002 (ambiguous) are the whole point of the
    issue.** A source that publishes symbols alone -- CHDgene is the likely
    first, see #31 -- gets joined onto HGNC ids somewhere, and the failure mode
    is not a crash: it is a row quietly not appearing. `scripts/convert_burden.py`
    already demonstrates the discipline, naming all 25 Excel date-corrupted
    symbols (`2001-03-01` = MARCH1/MARC1 collapsed) on every run, which is the
    only reason anyone knows about them.

    Warnings rather than errors, deliberately. An unresolvable symbol in an
    upstream file is not a defect in *this* repository's curation, and refusing
    to build over one would make a third party's naming choice able to take the
    site down. The gate is 0 errors; what this must do is be impossible to miss.
    """
    wanted = sorted(set(symbols))
    if registry is None:
        # **No skip warning when there was nothing to skip.** A skip means "a
        # check that would have run did not", and with no symbols to resolve no
        # check would have run. Reporting one anyway would put a lone warning in
        # a report that has no error to explain it -- and `ValidationReport.ok`
        # ignores warnings, so a warnings-only report exits 0. Pinned by
        # `test_every_skip_warning_arrives_with_the_error_that_caused_it`, which
        # caught exactly this.
        return _skipped(location, "symbol resolution") if wanted else []

    issues: list[ValidationIssue] = []
    for symbol, where in wanted:
        found = registry.resolve(symbol)
        if found.resolution is Resolution.UNKNOWN:
            issues.append(
                ValidationIssue(
                    code="GEN001",
                    severity=Severity.WARNING,
                    location=where,
                    message=(
                        f"gene symbol {symbol!r} resolves to no gene in mirrors/genes.tsv; "
                        "the row cannot be keyed on an HGNC id and would be dropped"
                    ),
                )
            )
        elif found.resolution is Resolution.AMBIGUOUS:
            issues.append(
                ValidationIssue(
                    code="GEN002",
                    severity=Severity.WARNING,
                    location=where,
                    message=(
                        f"gene symbol {symbol!r} is claimed by "
                        f"{', '.join(found.candidates)}; refusing to guess which"
                    ),
                )
            )
    return issues


def validate_mirror_symbols(
    labelled: Iterable[tuple[str, str, str]],
    registry: GeneRegistry | None,
    location: str,
) -> list[ValidationIssue]:
    """Every symbol our own mirrors print must resolve back to the id beside it.

    `labelled` is `(hgnc_id, symbol, where)` -- what a validity mirror already
    carries in both columns, so this is a *self-consistency* check on data that
    is already keyed correctly rather than a join.

    **GEN003 exists because the registry goes stale and nothing else notices.**
    Measured 2026-08-14 against the committed mirrors: `mirrors/genes.tsv`
    records HGNC:1152 as `POPDC1` with aliases `POP1|HBVES`, while ClinGen labels
    the same gene `BVES` -- and `BVES` is not among those aliases (`HBVES` is a
    different string). Nothing published is wrong today, because every join in
    this repository already keys on the id. But a symbol-keyed source naming
    `BVES` would resolve to nothing and the gene would vanish silently, which is
    the exact failure #33 exists to prevent, sitting in committed data.

    The remedy is a curator regenerating the mirror with
    `scripts/convert_hgnc.py` -- HGNC renames genes and `aliases` is where the
    old and new names are supposed to meet. Design decision D6 keeps that a
    human action rather than something a build does, so this reports and stops.

    A symbol resolving to a *different* gene is the serious case and is reported
    the same way rather than more loudly, because it is not reachable from a
    correct mirror: both columns come from the same upstream row. It is here so
    that a hand-edited mirror cannot introduce it undetected.
    """
    rows = sorted(set(labelled))
    if registry is None:
        # Same rule as above: a skip is only a skip when work existed.
        return _skipped(location, "mirror symbol") if rows else []

    known = {gene for gene, _, _ in rows} & _registry_ids(registry)
    issues: list[ValidationIssue] = []
    for gene, symbol, where in rows:
        # Only genes the registry actually carries. A mirror row for a gene
        # outside the 154-gene registry is not a naming disagreement -- the
        # registry is deliberately narrower than the mirrors, and reporting
        # those would drown the one real finding in 3,000 non-findings.
        if gene not in known:
            continue
        found = registry.resolve(symbol)
        if found.gene == gene:
            continue
        detail = (
            f"resolves to {found.gene}"
            if found.gene is not None
            else f"resolves to nothing ({found.resolution.value})"
        )
        issues.append(
            ValidationIssue(
                code="GEN003",
                severity=Severity.WARNING,
                location=where,
                message=(
                    f"{gene} is labelled {symbol!r} here, but that symbol {detail} "
                    "in mirrors/genes.tsv; regenerate the registry with "
                    "scripts/convert_hgnc.py so a symbol-keyed source can resolve it"
                ),
            )
        )
    return issues


def _skipped(location: str, what: str) -> list[ValidationIssue]:
    """The one GEN000, in one place, so both entry points word it identically."""
    return [
        ValidationIssue(
            code="GEN000",
            severity=Severity.WARNING,
            location=location,
            message=f"mirrors/genes.tsv could not be read; {what} checks skipped",
        )
    ]


def _registry_ids(registry: GeneRegistry) -> set[str]:
    """Every HGNC id the registry carries, via its own public lookup.

    Reaches through `resolve` on the approved symbols rather than touching a
    private attribute, so this module has no back channel into `GeneRegistry`'s
    internals and the two can be changed apart.
    """
    return {
        found.gene
        for found in (registry.resolve(symbol) for symbol in registry.approved_symbols())
        if found.gene is not None
    }
