# tests/unit/test_genes.py
"""Symbol -> HGNC id: resolved, or reported. Never guessed, never dropped. Issue #33."""

from __future__ import annotations

from chd_atlas.genes import GeneRegistry, Resolution
from chd_atlas.issues import Severity
from chd_atlas.validate.genes import validate_mirror_symbols, validate_symbol_resolution

# Two genes and a third that shares nothing, so every lookup below has a wrong
# answer available to it. `POPDC1` carries the real registry's own aliases, and
# `BVES` is deliberately absent from them -- that is the committed state this
# whole issue is about.
_ROWS = (
    {"hgnc_id": "HGNC:1152", "symbol": "POPDC1", "aliases": "POP1|HBVES"},
    {"hgnc_id": "HGNC:11604", "symbol": "TBX5", "aliases": "HOS"},
    {"hgnc_id": "HGNC:4173", "symbol": "GATA4", "aliases": ""},
)


def test_a_symbol_resolves_by_approved_name_then_by_alias() -> None:
    """Approved wins, an alias resolves, and each says which it was.

    The `Resolution` is not decoration. A caller that only got the id back could
    not tell a current symbol from a retired one, and those want different
    treatment when a source disagrees with the registry about a gene's name.
    """
    registry = GeneRegistry(_ROWS)
    assert len(registry) == 3

    approved = registry.resolve("TBX5")
    assert approved.resolution is Resolution.APPROVED
    assert approved.gene == "HGNC:11604"
    assert approved.ok is True

    alias = registry.resolve("HOS")
    assert alias.resolution is Resolution.ALIAS
    assert alias.gene == "HGNC:11604"

    # Surrounding whitespace is a formatting artefact of a source file, not a
    # different gene; a bare dict lookup would miss it and drop the row.
    assert registry.resolve("  TBX5  ").gene == "HGNC:11604"


def test_an_unknown_symbol_is_reported_and_never_resolved_to_a_guess() -> None:
    """The load-bearing rule: a symbol that maps to nothing must not vanish.

    `BVES` is the live case and the reason this issue exists. ClinGen labels
    HGNC:1152 `BVES`; `mirrors/genes.tsv` calls it `POPDC1` with aliases
    `POP1|HBVES` -- and `HBVES` is a *different string*. A symbol-keyed join on
    `BVES` therefore resolves to nothing, and a build that dropped the row would
    stay green with the gene simply absent.
    """
    registry = GeneRegistry(_ROWS)

    unknown = registry.resolve("BVES")
    assert unknown.resolution is Resolution.UNKNOWN
    assert unknown.gene is None
    assert unknown.ok is False
    assert unknown.candidates == ()

    assert registry.resolve("").resolution is Resolution.UNKNOWN
    assert registry.resolve("   ").resolution is Resolution.UNKNOWN


def test_a_symbol_two_genes_claim_is_refused_rather_than_picked() -> None:
    """Ambiguity names every candidate and resolves to none of them.

    `scripts/convert_burden.py` declines to match aliases at all, and its stated
    reason is this exact hazard: "an alias that is another gene's current symbol
    would silently attribute one gene's burden to another". Measured 2026-08-14
    the collision does not occur in the committed registry -- 154 approved
    symbols, 309 aliases, zero overlap -- so only a fixture can reach it, and it
    is guarded structurally rather than on the bet that the next upstream
    refresh keeps that true.

    **Approved beats alias, and that is not the same case.** HGNC retires a
    symbol by making it a *previous* symbol of the gene that kept the name, so a
    string that is one gene's approved symbol and another's alias names the
    first gene today -- resolving to the alias would attribute the row to a gene
    the source did not mean. That must resolve, not report.
    """
    shared = GeneRegistry(
        (
            {"hgnc_id": "HGNC:1", "symbol": "AAA", "aliases": "SHARED"},
            {"hgnc_id": "HGNC:2", "symbol": "BBB", "aliases": "SHARED"},
        )
    )
    ambiguous = shared.resolve("SHARED")
    assert ambiguous.resolution is Resolution.AMBIGUOUS
    assert ambiguous.gene is None
    assert ambiguous.candidates == ("HGNC:1", "HGNC:2")

    # One gene's approved symbol is another's alias: approved wins outright.
    contested = GeneRegistry(
        (
            {"hgnc_id": "HGNC:1", "symbol": "NAME", "aliases": ""},
            {"hgnc_id": "HGNC:2", "symbol": "OTHER", "aliases": "NAME"},
        )
    )
    resolved = contested.resolve("NAME")
    assert resolved.resolution is Resolution.APPROVED
    assert resolved.gene == "HGNC:1"


def test_a_row_with_no_id_is_skipped_without_taking_the_registry_with_it() -> None:
    """A blank `hgnc_id` has nothing to key on; every other row still loads."""
    registry = GeneRegistry(({"hgnc_id": "  ", "symbol": "GHOST", "aliases": ""}, *_ROWS))
    assert len(registry) == 3
    assert registry.resolve("GHOST").resolution is Resolution.UNKNOWN
    assert registry.resolve("TBX5").gene == "HGNC:11604"


def test_gen003_reports_a_mirror_symbol_the_registry_cannot_resolve() -> None:
    """The check that fires on the committed data, and what it must not report.

    Three rows, one finding. `POPDC1`/`BVES` is the disagreement; `TBX5` agrees
    and must stay silent; and HGNC:9999 is a gene the registry does not carry at
    all, which is **not** a naming disagreement -- `mirrors/genes.tsv` is
    deliberately narrower than the validity mirrors, and reporting those would
    bury the one real finding under thousands.
    """
    issues = validate_mirror_symbols(
        (
            ("HGNC:1152", "BVES", "mirrors/clingen_gene_validity.tsv"),
            ("HGNC:11604", "TBX5", "mirrors/clingen_gene_validity.tsv"),
            ("HGNC:9999", "NOTINREGISTRY", "mirrors/gencc_submissions.tsv"),
        ),
        GeneRegistry(_ROWS),
        "mirrors/genes.tsv",
    )

    assert [issue.code for issue in issues] == ["GEN003"]
    assert issues[0].severity is Severity.WARNING
    assert "HGNC:1152" in issues[0].message
    assert "'BVES'" in issues[0].message
    assert "convert_hgnc.py" in issues[0].message


def test_a_mirror_symbol_naming_a_different_gene_is_reported_as_such() -> None:
    """The serious case: the symbol resolves, but to somebody else.

    Unreachable from a correctly generated mirror -- both columns come from one
    upstream row -- so this exists to stop a hand-edited mirror introducing it
    undetected, and the message must name the gene it actually resolved to
    rather than saying only that something is wrong.
    """
    issues = validate_mirror_symbols(
        (("HGNC:1152", "TBX5", "mirrors/clingen_gene_validity.tsv"),),
        GeneRegistry(_ROWS),
        "mirrors/genes.tsv",
    )
    assert [issue.code for issue in issues] == ["GEN003"]
    assert "resolves to HGNC:11604" in issues[0].message


def test_an_unreadable_registry_reports_one_skip_and_not_one_per_symbol() -> None:
    """GEN000, the REF000/SRC000/ONT000/SCP000 precedent.

    A missing reference table must say "checks skipped" once. Reporting every
    symbol as unresolvable would bury the real error -- the registry itself --
    under one warning per row, with nothing naming the cause. Both entry points
    take the same branch, because a caller that got a clean report from one and
    a cascade from the other would trust the wrong one.
    """
    skipped = validate_mirror_symbols(
        (("HGNC:1152", "BVES", "x"), ("HGNC:11604", "TBX5", "y")), None, "mirrors/genes.tsv"
    )
    assert [issue.code for issue in skipped] == ["GEN000"]
    assert skipped[0].severity is Severity.WARNING
    assert "skipped" in skipped[0].message

    resolution = validate_symbol_resolution(
        (("BVES", "x"), ("HOS", "y")), None, "mirrors/genes.tsv"
    )
    assert [issue.code for issue in resolution] == ["GEN000"]


def test_gen001_and_gen002_report_an_unmapped_and_an_ambiguous_symbol() -> None:
    """The entry point a symbol-keyed source will use. No caller in the build yet.

    **Stated plainly because it matters:** nothing in `build/` calls
    `validate_symbol_resolution` today, since every source currently mirrored --
    ClinGen and GenCC -- is already keyed on HGNC id. It exists so the machinery
    is in place before the first symbol-only source lands (CHDgene is the likely
    one, #31), and it is exercised here rather than left as a guard nobody has
    watched fail.

    A resolvable symbol produces no issue at all: a validator that reported
    every row would be noise, and the report has to be readable to be read.
    """
    issues = validate_symbol_resolution(
        (("BVES", "chdgene.csv:row 4"), ("SHARED", "chdgene.csv:row 9"), ("TBX5", "ok")),
        GeneRegistry(
            (
                *_ROWS,
                {"hgnc_id": "HGNC:1", "symbol": "AAA", "aliases": "SHARED"},
                {"hgnc_id": "HGNC:2", "symbol": "BBB", "aliases": "SHARED"},
            )
        ),
        "mirrors/genes.tsv",
    )

    assert sorted(issue.code for issue in issues) == ["GEN001", "GEN002"]
    unmapped = next(issue for issue in issues if issue.code == "GEN001")
    ambiguous = next(issue for issue in issues if issue.code == "GEN002")
    assert unmapped.location == "chdgene.csv:row 4"
    assert "would be dropped" in unmapped.message
    assert "HGNC:1, HGNC:2" in ambiguous.message
    assert "refusing to guess" in ambiguous.message
