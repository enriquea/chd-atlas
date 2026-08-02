# src/chd_atlas/build/runner.py
"""Orchestration: gate on validation, then emit every artifact in order.

This is the only module that decides *whether* to build. Everything under
`build/` publishes what it is given; the judgement that the corpus is fit to
publish is made once, here, and made before the first byte is written.

The gate is `validate_repository` itself rather than a rule of this module's
own. A stricter rule here would make a repository that CI calls green
unpublishable; a looser one would publish a corpus the validator rejected. Both
are drift, and the only way to not drift is to ask the same question the same
way — hence `report.ok`, which ignores warnings by design, rather than a
condition assembled here.

Order matters in exactly one place: the manifest is last, because
`emitter.checksums` is complete only once every other builder has run. That is
no longer a convention to remember — `write_manifest` seals the emitter, so a
builder added below it fails the build rather than publishing a file the
manifest does not list.
"""

from __future__ import annotations

from pathlib import Path

from chd_atlas.build.bundles import build_genes
from chd_atlas.build.emit import Emitter
from chd_atlas.build.literature import build_literature
from chd_atlas.build.manifest import source_commit, write_manifest
from chd_atlas.build.omics import build_omics
from chd_atlas.build.search import GeneLabels, build_search
from chd_atlas.build.variants import build_variants
from chd_atlas.corpus import load_curation
from chd_atlas.tables import TABLE_SCHEMAS, read_table
from chd_atlas.validate.runner import validate_repository


class BuildRefused(Exception):
    """The repository did not pass validation, so nothing was written.

    Publishing from a corpus that fails its own gate would defeat the point of
    having one: the site would carry claims the validator rejected, and nobody
    would find out until a reader did. The message carries the rendered report,
    because a curator who ran `build` needs to know what to fix without being
    told to go and run a second command.
    """


def _cell(value: object) -> str | None:
    """One mirror cell as text, with blank and whitespace-only treated as absent.

    `read_table` maps the empty string to null, so TBL003 catches `symbol` left
    empty — and catches nothing when the cell holds a single space. A generated
    mirror produces exactly that, and the space survived every check to be
    published as the gene's browse label and its search result: a row rendering
    as nothing, on a build reporting zero errors, against a `docs/data-api.md`
    that promises `symbol` is "always a non-empty string you can render and
    search". Found by adversarial review of #2.

    Stripped rather than reported, because the fallback the caller already has —
    keep the HGNC id — is better than either publishing the blank or refusing to
    build over one stray character.
    """
    if value is None:
        return None
    text = str(value).strip()
    return text or None


def _gene_registry(root: Path) -> dict[str, GeneLabels]:
    """What `mirrors/genes.tsv` says about each gene, keyed on HGNC id.

    One `GeneLabels` per row rather than parallel mappings, so a gene cannot end
    up with a symbol and no name. `build_genes` needs only the symbol and is
    handed a projection of this rather than a second read of the mirror.

    The aliases cell is split here, and this is the only place that knows the
    separator. The mirror stores them pipe-separated because the schema keeps
    every TSV one value per cell rather than inventing a nested format, so
    unpacking that is the reader's job. Handing the raw cell to `GeneLabels`
    would type-check — `frame.to_dicts()` yields `dict[str, Any]`, and an `Any`
    satisfies `tuple[str, ...]` silently — and would publish a search index of
    single characters; the guard in `GeneLabels.__post_init__` is what catches
    that, and it exists because this line is where the mistake is made.

    Returning empty for a missing or unreadable mirror is unreachable behind
    `build_site`, which refuses on any error and so never reaches here without
    TBL008 having fired. It is not a fallback anyone should rely on: it is what
    keeps this function total for a direct caller, in the same way
    `search.py`'s symbol fallback is.
    """
    path = root / "mirrors" / "genes.tsv"
    if not path.is_file():
        return {}
    frame, _ = read_table(path, TABLE_SCHEMAS["genes"])
    if frame is None or "hgnc_id" not in frame.columns:
        return {}

    registry: dict[str, GeneLabels] = {}
    for row in frame.to_dicts():
        gene = _cell(row.get("hgnc_id"))
        symbol = _cell(row.get("symbol"))
        # Both are non-nullable in the schema, so a blank here is a bypassed
        # TBL003. A row with no id has nothing to key on and no row to correct it
        # against; a row with no symbol would label the gene with the string
        # "None" in the browse row and the page heading alike. Dropping the row
        # instead lets `build_genes` and `search.py` fall back to the HGNC id,
        # which renders and searches.
        if gene is None or symbol is None:
            continue
        raw = row.get("aliases")
        registry[gene] = GeneLabels(
            symbol=symbol,
            name=_cell(row.get("name")),
            aliases=tuple(part.strip() for part in str(raw).split("|") if part.strip())
            if raw
            else (),
        )
    return registry


def build_site(root: Path, out: Path) -> dict[str, str]:
    """Build the published API into `out`. Returns path-to-checksum per file.

    Raises `BuildRefused` before writing anything if validation reports errors.
    Nothing is written on *that* refusal — not a partial `dist/`, not an empty
    directory — because a half-built site a deploy step then uploads is worse
    than no build at all.

    That is not a guarantee about every failure, and saying so once read as one.
    The builders below carry their own guards, and `validate/` does not yet check
    everything they do — a variant shard named `chr12.tsv` validates clean and
    raises here — so a `ValueError` from one of them can arrive after earlier
    builders have written. `out` may then hold a partial site with no manifest.
    The manifest is the signal: it is written last and seals the emitter, so a
    `dist/` without one is a build that did not finish. `cli.build` says so in
    as many words rather than leaving a caller to infer it.

    The returned mapping covers every artifact including `manifest.json`, which
    the manifest itself cannot list. `dict()` is hygiene and nothing more: the
    emitter is local and sealed, so it is unreachable once this returns and no
    test can tell a copy from the original. Measured, not assumed — the mutant
    returning `emitter.checksums` directly survives the whole suite, which is
    what an equivalent mutant looks like rather than a gap in the tests.
    """
    report = validate_repository(root)
    if not report.ok:
        raise BuildRefused(
            f"refusing to build: {report.error_count} validation error(s)\n{report.render()}"
        )

    # Issues are discarded rather than re-reported: `validate_repository` above
    # loaded the same corpus and its report is what the gate just passed, so a
    # second copy here could only ever say the same thing or disagree with the
    # decision already made.
    corpus, _ = load_curation(root)
    genes = _gene_registry(root)

    emitter = Emitter(root=out)
    omics = build_omics(root, emitter)
    variants = build_variants(root, emitter)
    build_genes(
        corpus,
        emitter,
        symbols={gene: labels.symbol for gene, labels in genes.items()},
        omics=omics,
        variants=variants,
    )
    build_literature(corpus, emitter)
    build_search(corpus, emitter, genes=genes)
    # Last, and enforced as last: this seals the emitter.
    write_manifest(corpus, emitter, commit=source_commit(root))
    return dict(emitter.checksums)
