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
from chd_atlas.build.landing import build_landing
from chd_atlas.build.literature import build_literature, build_sources
from chd_atlas.build.manifest import source_commit, write_manifest
from chd_atlas.build.omics import build_omics
from chd_atlas.build.pages import build_gene_index_page, build_gene_pages
from chd_atlas.build.search import GeneLabels, build_search
from chd_atlas.build.validity import gene_validity, published_genes
from chd_atlas.build.variants import build_variants
from chd_atlas.corpus import Corpus, load_curation
from chd_atlas.models.assertion import LesionAssertion
from chd_atlas.tables import TABLE_SCHEMAS, read_table
from chd_atlas.validate.runner import validate_repository
from chd_atlas.validate.sources import load_sources


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
        # Both `hgnc_id` and `symbol` are non-nullable in the schema, so a blank
        # here is a bypassed TBL003. The two blanks are not the same failure,
        # though. A row with no id has nothing to key the registry on and no
        # row to correct it against, so it is dropped — there is no gene here to
        # attach `name` or `aliases` to. A row with no symbol is still a real,
        # validated gene: `name` and `aliases` are perfectly good and reached
        # this reader having passed every check that applies to them, so
        # dropping the row along with the blank symbol would discard both for
        # no reason a curator chose. The symbol falls back to the HGNC id
        # instead — the same fallback `build_genes` and `search.py` already
        # apply to a gene missing from this registry entirely — so the row, its
        # name and its aliases all still reach `genes/index.json` and the
        # search index, rather than the gene silently losing everything but its
        # id. Found by adversarial review after release: the previous version
        # of this guard dropped the row and, with it, a name and aliases that
        # had already been validated and were fit to publish.
        if gene is None:
            continue
        raw = row.get("aliases")
        registry[gene] = GeneLabels(
            symbol=_cell(row.get("symbol")) or gene,
            name=_cell(row.get("name")),
            aliases=tuple(part.strip() for part in str(raw).split("|") if part.strip())
            if raw
            else (),
        )
    return registry


def _assertions_by_gene(corpus: Corpus) -> dict[str, list[LesionAssertion]]:
    """Every curated assertion, grouped by the gene it is about, ordered by id.

    `pages.py` renders one `<h3>` block per assertion in the order it is handed
    them, so this order is part of a page's bytes and therefore its checksum.
    Ordered by id rather than left in corpus order for the reason
    `bundles._records_by_gene` gives: corpus order is `load_curation`'s — record
    files sorted by filename, then each file's own order within it — which makes
    it an artefact of how the curation is filed. Without this, moving an
    assertion from `TBX5.yaml` into a second file would reorder a page that is
    otherwise unchanged.

    Not shared with `bundles._records_by_gene`, which groups the same records
    the same way and then dumps each one to JSON: these builders need the models
    themselves, and a helper returning both shapes would be doing two jobs.

    Genes outside the published population are grouped here too and never looked
    up — `build_gene_pages` iterates `facts`, so what it renders is decided by
    the population, not by this mapping.
    """
    grouped: dict[str, list[LesionAssertion]] = {}
    for assertion in corpus.assertions:
        grouped.setdefault(assertion.gene, []).append(assertion)
    return {gene: sorted(items, key=lambda record: record.id) for gene, items in grouped.items()}


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

    # Unreachable behind the gate above in the same way `_gene_registry`'s
    # empty-mirror fallback is: TBL012 is an error, so `validate_repository`
    # already refused a repository whose validity mirrors do not read, and this
    # line never runs on one that reaches here. Read again anyway and `raise`
    # rather than assume it: `-O` strips `assert`, and this project keeps a
    # guard on every bypassed gate (`Emitter.seal`, `encode_json`'s
    # `allow_nan=False`) rather than trusting "the gate already checked this" to
    # stay true as the two modules evolve apart.
    scope_terms = {str(entry.id) for entry in corpus.chd_scope}
    clingen, _ = read_table(
        root / "mirrors" / "clingen_gene_validity.tsv", TABLE_SCHEMAS["clingen_validity"]
    )
    gencc, _ = read_table(
        root / "mirrors" / "gencc_submissions.tsv", TABLE_SCHEMAS["gencc_submissions"]
    )
    if clingen is None or gencc is None:
        raise ValueError("a validity mirror could not be read; the gate should have refused first")
    validity = gene_validity(clingen, gencc, in_scope=scope_terms)

    # Computed once and reused by both callers below: `build_genes` needs the
    # symbol alone, and so does `build_landing`, which is not a second read of
    # the mirror either — both project the same `genes` registry the same way.
    symbols = {gene: labels.symbol for gene, labels in genes.items()}

    # D21's population, computed once for the same reason `symbols` is, and with
    # more riding on it: `build_genes` publishes one index row and one bundle per
    # member, `build_search` indexes one record per member so that every one of
    # those bundles is reachable from the search box, and `build_landing` prints
    # how many there are. Deriving it three times would let the front page state a
    # figure the browse payload contradicts, or let the search index advertise a
    # `genes/<id>.json` no builder wrote — the shape the "154 genes published"
    # defect already had once.
    published = published_genes(validity)

    emitter = Emitter(root=out)
    omics = build_omics(root, emitter)
    variants = build_variants(root, emitter)
    # `facts` rather than a second `gene_facts` call below: the pages and the
    # bundles render from one derivation, so a page cannot state a confidence the
    # bundle it links to contradicts. See `build_genes`' docstring.
    facts = build_genes(
        corpus,
        emitter,
        symbols=symbols,
        omics=omics,
        variants=variants,
        validity=validity,
        published=published,
    )
    build_literature(corpus, emitter)
    # Read again rather than threaded down from the gate: `validate_repository`
    # discards its registry, and re-reading one small YAML costs less than
    # widening the gate's return type to carry a value only the build wants.
    # It cannot fail here — SRC001 is an error, so the gate refused already.
    registry, _ = load_sources(root)
    build_sources(registry, emitter)
    build_search(corpus, emitter, genes=genes, published=published)
    build_landing(corpus, symbols=symbols, validity=validity, published=published, emitter=emitter)
    # The HTML over everything above. Wired here and nowhere else: until this
    # call existed, `pages.py` was imported by its unit test alone, `build_landing`
    # and the shared `<nav>` both linked to `genes/index.html`, and no build wrote
    # one — a green build, every checksum verifying, and every visitor clicking
    # "Genes" served a 404.
    #
    # Keyed by id rather than by position, because both builders look records up
    # by identifier: `publications` resolves the PMID an assertion's evidence
    # cites, the same lookup `literature._featured` performs over the same tuple.
    build_gene_pages(
        facts,
        emitter,
        symbols=symbols,
        validity=validity,
        assertions=_assertions_by_gene(corpus),
        publications={publication.id: publication for publication in corpus.publications},
    )
    # `validity` again, and the same object `build_gene_pages` was handed: the
    # browse row's `definitive for` cell and the gene page's `definitive for`
    # rail row must name one disease, and they do because both read this mapping
    # rather than either re-deriving it from the mirrors.
    build_gene_index_page(facts, emitter, symbols=symbols, validity=validity)
    # Last, and enforced as last: this seals the emitter.
    write_manifest(corpus, emitter, commit=source_commit(root))
    return dict(emitter.checksums)
