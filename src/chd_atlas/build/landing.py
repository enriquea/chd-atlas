# src/chd_atlas/build/landing.py
"""`index.html`: the site's front door, and the one page a person opens directly.

Everything else under `build/` is JSON for a program to fetch, or a gene page
`pages.py` renders. Before this module existed the live site measured `200` on
`manifest.json` and `genes/index.json` and `404` on `/` — the data API worked
and nothing told a visitor what they were looking at, so this is not a cosmetic
gap: a reader who follows the repository link lands nowhere.

The document scaffolding is not this module's. `render.document` supplies the
`<!doctype>`, the `<head>`, the navigation and the one `<style>` block, whose
contents are `render.STYLESHEET` — so this page, `genes/index.html` and every
gene page share a single stylesheet rather than carrying three that drift. What
stays here is the page's content. Self-contained by construction, like every
other artifact this build writes: `document` inlines that style and emits no
`<link>`, `<script src>` or `@import`, so the page renders with no network
request beyond the one that fetched it.

Deterministic for the reason every other builder here is. No build timestamp;
every number on the page is read from what `build_site` already assembled —
the corpus, and the same `published` and `validity` objects `build_genes` was
handed — never hardcoded.

Three gene populations reach this page, and they are three different claims:

- **Genes published** is `len(published)`, D21's population: one gene per member
  of the set `build_genes` keys `genes/index.json` on, each of them a gene a
  ClinGen expert panel calls definitive for an in-scope disease. 23 in the
  committed corpus. Counted from the *same object* `build_genes` was handed,
  not re-derived here, so the front page and the browse payload cannot state
  different numbers.
- **Genes the atlas has curated** is the distinct genes named by
  `corpus.assertions` — the genes carrying evidence authored here, 1 in the
  committed corpus. It is derived once and shared with the sentence that names
  them, so the figure and the list cannot disagree. Both rows are shown because
  22 of the 23 published genes carry no assertion at all: a single "genes"
  figure would read as coverage this site does not have.
- `len(validity)` is neither of those — every gene either mirror curates, 154 in
  the committed corpus. It is kept, but under "Where this data comes from" and
  labelled so it cannot be read as what this site browses.

`tests/unit/test_build_landing.py` guards those three against each other by
fixture rather than by recomputation, and the description here said otherwise
until it was checked. Exactly one of its tests recomputes a number and compares
it against the rendered text —
`test_the_published_gene_count_agrees_with_a_real_build_of_genes_index_json`,
which reads `len(genes/index.json["genes"])` off a real build. Every other
count assertion in the file is a literal.

They are literals on purpose, and the fixture is what makes them strong: it is
hand-sized so that all three populations differ — two assertions, on one gene,
against three published genes and five mirrored ones — so a figure wired to the
wrong population renders a number that is true of *something* and still fails.
None of the three matches the committed corpus's 1 / 23 / 154 either, and
`assert "154" not in text` is what catches the count being hardcoded back.

Every value that reaches the page from curated or mirrored text — a gene
symbol, an HGNC id — goes through `html.escape`. Nothing curated here is
reviewed the way handwritten page copy is, and `bundles.py` and `search.py`
already treat curator- and mirror-supplied strings as data rather than markup
for the same reason; a `<` in a gene symbol is unlikely, not impossible, and
this is the one artifact in the build where that distinction is the difference
between a string and a script tag.
"""

from __future__ import annotations

import html
from collections.abc import Collection, Mapping, Sequence

from chd_atlas.build.emit import Emitter
from chd_atlas.build.paths import LANDING
from chd_atlas.build.render import document
from chd_atlas.build.validity import GeneValidity
from chd_atlas.corpus import Corpus

# Read-only, and not a claim this build can verify — unlike every count on the
# page, which is derived from the corpus this module was handed. Kept as a
# literal because there is no data source for it to drift out of sync with; a
# wrong repository link is a code review problem, not a silent-evidence-loss
# one.
_REPOSITORY_URL = "https://github.com/enriquea/chd-atlas"

# A caption that told a reader to wait for something that had already happened.
#
# It read "(browsable once ClinGen grades it definitive)", and the gate is not
# "definitive" -- it is definitive *for one of the diseases this atlas calls
# CHD*. `validity.published_genes` requires a ClinGen `Definitive` record, and
# `gene_validity` has already dropped every record naming an out-of-scope
# disease before that test runs. Measured 2026-08-04 against the committed
# mirrors: of the 154 genes counted by this row, **20** already carry a ClinGen
# `Definitive` grade and are still not browsable, because the grade is for a
# disease outside CHD scope -- ACTC1, ACVR1, ANKRD11, BMPR1A, BMPR2, COL1A2,
# ELN, FLNA, FOXP1, KDR, MED13L, MYBPC3, MYH11, MYH7, NOTCH1, PDGFRA, POPDC1,
# RECQL4, SCN5A, SMAD6. ELN was graded Definitive in 2024, for cutis laxa.
#
# A constant rather than a literal in the template because the label must render
# on one line for the row to read as one label, and the sentence no longer fits
# in 100 characters of source.
_MIRRORED_ROW_LABEL = (
    "Genes with mirrored validity in CHD scope "
    "(browsable once ClinGen grades it definitive for a disease in that scope)"
)


def _plural(count: int, noun: str) -> str:
    """The only pluralisation this page needs: `1 assertion` vs `2 assertions`."""
    return f"{count} {noun}" if count == 1 else f"{count} {noun}s"


def _asserted_genes(genes: Sequence[str], symbols: Mapping[str, str]) -> str:
    """Every gene carrying a curated assertion, escaped and comma-joined.

    Takes the genes rather than deriving them, so the "Genes the atlas has
    curated" figure and this sentence come from one set: a row saying 2 beside a
    sentence naming 3 genes is the shape of drift this removes.

    `genes` is sorted by HGNC id upstream, the same order `genes/index.json`
    publishes in, so this sentence and that file agree on which gene is named
    first without either reading the other. `symbols.get(gene, gene)` is the
    same registry fallback `build_genes` and `search.py` already apply to a gene
    absent from `mirrors/genes.tsv` — a symbol that has not been mirrored yet
    still reads as its HGNC id rather than as nothing.
    """
    if not genes:
        return "none yet"
    return ", ".join(
        f"{html.escape(symbols.get(gene, gene))} ({html.escape(gene)})" for gene in genes
    )


def _render(
    corpus: Corpus,
    symbols: Mapping[str, str],
    validity: Mapping[str, GeneValidity],
    published: Collection[str],
) -> str:
    assertion_count = len(corpus.assertions)
    # The very set `gene_facts` keys `genes/index.json` on, not a recomputation
    # of it. Deliberately not `len(validity)`, which counts every gene either
    # mirror curates whether or not a panel calls it definitive — that count
    # belongs to "Where this data comes from", not to what this build publishes.
    published_gene_count = len(published)
    mirrored_gene_count = len(validity)
    # A third population, and the narrowest: the genes this atlas has curated
    # itself. Derived once and handed to `_asserted_genes`, so the count and the
    # sentence naming those genes cannot drift apart.
    curated_genes = sorted({assertion.gene for assertion in corpus.assertions})
    asserted = _asserted_genes(curated_genes, symbols)

    body = f"""  <h1>CHD Atlas</h1>
  <p class="tagline">Curated evidence linking genes, variants and proteins to congenital
    heart disease.</p>

  <h2>What this is</h2>
  <p>A curator reads the primary literature and mirrors upstream reference sources —
    ClinGen, GenCC, HPO — to build a static, queryable record of the evidence behind a
    gene's association with congenital heart disease. There is no server: every page
    below is built ahead of time, and every payload behind one is a JSON file a client
    fetches and filters itself.</p>

  <h2>What's published</h2>
  <dl>
    <dt>Curated gene-disease assertions</dt>
    <dd>{assertion_count}</dd>
    <dt>Genes published</dt>
    <dd>{published_gene_count}</dd>
    <dt>Genes the atlas has curated</dt>
    <dd>{len(curated_genes)}</dd>
    <dt>Functional evidence records</dt>
    <dd>{len(corpus.functional)}</dd>
    <dt>Publications cited</dt>
    <dd>{len(corpus.publications)}</dd>
    <dt>Phenotype terms (HPO)</dt>
    <dd>{len(corpus.phenotypes)}</dd>
    <dt>Featured manuscripts</dt>
    <dd>{len(corpus.featured)}</dd>
    <dt>Omics datasets</dt>
    <dd>{len(corpus.datasets)}</dd>
  </dl>
  <p>{_plural(assertion_count, "curated gene-disease assertion")} so far: {asserted}.</p>

  <h2>Where this data comes from</h2>
  <p>Gene-disease validity classifications are mirrored from <strong>ClinGen</strong> and
    <strong>GenCC</strong>; phenotype terms and their synonyms are mirrored from the
    <strong>Human Phenotype Ontology</strong>. The atlas authors no validity classification
    of its own — every classification on this site is an upstream panel's or submitter's,
    republished with its provenance intact. <a href="sources.json">sources.json</a> carries
    the licence and attribution terms each of these is mirrored under.</p>
  <dl>
    <dt>{_MIRRORED_ROW_LABEL}</dt>
    <dd>{mirrored_gene_count}</dd>
  </dl>

  <h2>Browse the data</h2>
  <ul class="links">
    <li><a href="genes/index.html"><code>genes/index.html</code></a> — every gene the atlas
      publishes, with filters.</li>
    <li><a href="genes/index.json"><code>genes/index.json</code></a> — every gene the atlas
      publishes, browse-row summaries and mirrored validity included.</li>
    <li><a href="manifest.json"><code>manifest.json</code></a> — what this build produced,
      and a checksum for every file in it.</li>
    <li><a href="sources.json"><code>sources.json</code></a> — what the atlas mirrors, and
      on whose terms.</li>
    <li><a href="{_REPOSITORY_URL}">{html.escape(_REPOSITORY_URL)}</a> — source, curation
      workflow and the full published API documentation.</li>
  </ul>

  <footer>
    <p>Read <a href="{_REPOSITORY_URL}/blob/main/docs/data-api.md">the published data API
      documentation</a> before writing against these files — in particular the note on
      contested genes, which must never be displayed as settled.</p>
  </footer>"""

    return document(title="CHD Atlas", root="", body=body)


def build_landing(
    corpus: Corpus,
    symbols: Mapping[str, str],
    validity: Mapping[str, GeneValidity],
    published: Collection[str],
    emitter: Emitter,
) -> None:
    """Emit `index.html`.

    `symbols` is the same `{gene: labels.symbol for gene, labels in genes.items()}`
    projection `build_genes` takes, not a second read of the mirror. `validity`
    is `build.validity.gene_validity()`'s return, the same mapping `build_genes`
    reduces through `derive.gene_facts` — read directly here rather than through
    that reduction, because this page wants only `len()`, not the per-gene
    fields `gene_facts` computes for a browse row.

    `published` is the same object `build_genes` keys `genes/index.json` on, for
    the reason that file exists to be counted rather than re-derived: the "Genes
    published" figure on the front page is a claim about that payload, and the
    only way it cannot contradict it is to count the set that produced it. Taken
    as a `Collection` and counted with `len`, which is exact because
    `build_genes` runs first in `build_site` and `Emitter` refuses to write one
    bundle path twice — a `published` carrying a duplicate fails the build there
    rather than inflating a figure here.
    """
    emitter.write_text(LANDING, _render(corpus, symbols, validity, published))
