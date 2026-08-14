# src/chd_atlas/genes.py
"""The one place a gene symbol becomes an HGNC id. Issue #33.

**Every cross-source join in this atlas keys on HGNC id, never on symbol**, and
this module exists so that a source which publishes symbols alone can be brought
onto that key in exactly one place, under exactly one set of rules.

The owner's instruction, 2026-08-06:

> there will be a lot of gene mapping across sources. you need to ensure this
> occur at the gene id level and not at the gene symbol name, since the later is
> more prone to mismatch due to aliases.

The hazard is measured, not hypothetical. Across the two committed validity
mirrors (2026-08-14) **six HGNC ids carry two different symbols**, because
ClinGen and GenCC disagree about the approved name for one gene: HGNC:1152
(`BVES`/`POPDC1`), HGNC:17272 (`CENPJ`/`CPAP`), HGNC:17861
(`IFT54`/`TRAF3IP1`), HGNC:26742 (`ASPNAT`/`NAT8L`), HGNC:28178
(`CCDC115`/`VMA22`) and HGNC:7687 (`COXFA4`/`NDUFA4`). A symbol-keyed map built
first-seen-wins picks whichever the read order surfaced, which is silent and not
even stable across a change in read order.

**An unresolved symbol is reported, never dropped.** That is the load-bearing
rule and the reason `resolve` returns a `Resolution` rather than `str | None`: a
caller that cannot tell "no such symbol" from "this symbol names two genes"
cannot report either, and a gene that fails to map and simply vanishes is
precisely the silent evidence loss this project treats as its characteristic
failure -- the build stays green and the gene is not there.

**Ambiguity is refused, never guessed.** `scripts/convert_burden.py` declines to
match aliases at all, and its stated reason is exactly right: "an alias that is
another gene's current symbol would silently attribute one gene's burden to
another". Measured against the committed registry, that collision does not occur
today -- 154 approved symbols, 309 alias symbols, **zero** overlap between the
two sets and **zero** aliases claimed by more than one gene -- but the guard is
structural rather than a bet on that staying true, because the registry is
regenerated from upstream and the number that matters is the one after the next
refresh.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from enum import StrEnum


class Resolution(StrEnum):
    """Why `GeneRegistry.resolve` returned what it did."""

    APPROVED = "approved"
    """The symbol is some gene's current HGNC-approved symbol."""

    ALIAS = "alias"
    """The symbol is an alias or previous symbol of exactly one gene."""

    UNKNOWN = "unknown"
    """No gene in the registry claims this symbol. Report it; do not drop the row."""

    AMBIGUOUS = "ambiguous"
    """More than one gene claims it. Report it; never pick one."""


@dataclass(frozen=True)
class Resolved:
    """One symbol lookup: the outcome, and the id if there is exactly one.

    `gene` is `None` for `UNKNOWN` and for `AMBIGUOUS`, and `candidates` carries
    every id that claimed the symbol so a report can name them. The two are kept
    apart deliberately -- a caller that unpacked a bare `str | None` would treat
    "nobody claims this" and "two genes claim this" identically, and those need
    opposite responses from a curator.
    """

    symbol: str
    resolution: Resolution
    gene: str | None = None
    candidates: tuple[str, ...] = ()

    @property
    def ok(self) -> bool:
        return self.gene is not None


class GeneRegistry:
    """`mirrors/genes.tsv` as a symbol -> HGNC id lookup.

    Built from the rows of the registry rather than reading the file, so this
    stays a pure structure over already-loaded data -- the same shape
    `validate_scope_terms` and `validate_source_references` have, and what lets
    it be tested without a filesystem.

    **Approved symbols win over aliases, and the precedence is not a tie-break
    to be tidied away.** HGNC retires a symbol by making it a *previous* symbol
    of the gene that kept the name, so a string that is one gene's approved
    symbol and another's alias names the first gene today. Resolving to the
    alias would attribute a row to a gene the source did not mean.
    """

    def __init__(self, rows: Iterable[Mapping[str, object]]) -> None:
        approved: dict[str, str] = {}
        aliases: dict[str, set[str]] = {}
        for row in rows:
            gene = _text(row.get("hgnc_id"))
            if gene is None:
                # No id is nothing to key on. `_gene_registry` in
                # `build/runner.py` drops these too, and TBL003 reports the
                # blank, so this is a bypassed-gate path rather than a silent
                # policy of its own.
                continue
            symbol = _text(row.get("symbol"))
            if symbol is not None:
                approved[symbol] = gene
            raw = _text(row.get("aliases"))
            for alias in (part.strip() for part in (raw or "").split("|")):
                if alias:
                    aliases.setdefault(alias, set()).add(gene)
        self._approved = approved
        self._aliases = aliases

    def __len__(self) -> int:
        return len(self._approved)

    def approved_symbols(self) -> tuple[str, ...]:
        """Every current HGNC-approved symbol the registry carries, sorted.

        Public so that a caller needing the registry's population does not reach
        into `_approved`. Sorted because a caller iterating it may put the result
        in a payload, and `encode_json`'s `sort_keys` orders dict keys only.
        """
        return tuple(sorted(self._approved))

    def resolve(self, symbol: str) -> Resolved:
        """One symbol, resolved or reported. Never guesses, never raises."""
        name = symbol.strip()
        if not name:
            return Resolved(symbol=symbol, resolution=Resolution.UNKNOWN)
        gene = self._approved.get(name)
        if gene is not None:
            return Resolved(
                symbol=name, resolution=Resolution.APPROVED, gene=gene, candidates=(gene,)
            )
        claimed = self._aliases.get(name)
        if not claimed:
            return Resolved(symbol=name, resolution=Resolution.UNKNOWN)
        if len(claimed) > 1:
            return Resolved(
                symbol=name,
                resolution=Resolution.AMBIGUOUS,
                candidates=tuple(sorted(claimed)),
            )
        only = next(iter(claimed))
        return Resolved(symbol=name, resolution=Resolution.ALIAS, gene=only, candidates=(only,))


def _text(value: object) -> str | None:
    """One registry cell as text, with blank and whitespace-only treated as absent.

    The same rule `build/runner.py::_cell` applies, and for the same reason: a
    single space survives a `min_length` check and renders as nothing, so it
    would enter the index as a symbol matching a query nobody typed.
    """
    if value is None:
        return None
    text = str(value).strip()
    return text or None
