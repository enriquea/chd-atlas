# src/chd_atlas/build/landing.py
"""`index.html`: the site's front door, and the one page a person opens directly.

Everything else under `build/` is JSON for a program to fetch. Before this
module existed the live site measured `200` on `manifest.json` and
`genes/index.json` and `404` on `/` — the data API worked and nothing told a
visitor what they were looking at, so this is not a cosmetic gap: a reader who
follows the repository link lands nowhere.

Self-contained by construction, like every other artifact this build writes:
no `<link>`, `<script src>` or `@import`, only an inline `<style>` block, so
the page renders with no network request beyond the one that fetched it.

Deterministic for the reason every other builder here is. No build timestamp;
every number on the page is read from what `build_site` already assembled for
`build_genes`, never hardcoded. That guarantees each figure matches what this
module was *handed* — it does not, on its own, guarantee that figure matches
what another builder chose to *publish* from the same input. `len(validity)` is
every gene either mirror curates, full stop; `genes/index.json` carries one row
per member of `published`, D21's narrower population (a ClinGen expert panel
calls the gene definitive for an in-scope disease). The two counts genuinely
diverge — 154 mirrored genes against 23 published ones in the committed corpus
— so this page takes its published-gene figure from the *same* `published` set
`build_genes` was handed, not from `len(validity)` and not from a second
derivation of its own: one object, so the front page and the browse payload
cannot state different numbers. The mirrored-validity count is kept, but
labelled and placed so it cannot be read as coverage this site browses.
`tests/unit/test_build_landing.py` recomputes the same numbers from the same
fixture and compares them against the rendered text rather than against a
literal, for exactly that reason — a hardcoded expectation in the test would
only prove the two hardcodings agree with each other — and checks the
published-gene figure against a real build's `genes/index.json` besides.

The count of curated assertions, and the sentence naming the genes that carry
them, stay keyed on `corpus.assertions`. They are a different claim from
"published": 22 of the 23 genes this site publishes carry no assertion at all,
and a page that conflated the two would read as though every published gene had
been curated here.

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
from collections.abc import Collection, Mapping

from chd_atlas.build.emit import Emitter
from chd_atlas.build.paths import LANDING
from chd_atlas.build.validity import GeneValidity
from chd_atlas.corpus import Corpus

# Read-only, and not a claim this build can verify — unlike every count on the
# page, which is derived from the corpus this module was handed. Kept as a
# literal because there is no data source for it to drift out of sync with; a
# wrong repository link is a code review problem, not a silent-evidence-loss
# one.
_REPOSITORY_URL = "https://github.com/enriquea/chd-atlas"


def _plural(count: int, noun: str) -> str:
    """The only pluralisation this page needs: `1 assertion` vs `2 assertions`."""
    return f"{count} {noun}" if count == 1 else f"{count} {noun}s"


def _asserted_genes(corpus: Corpus, symbols: Mapping[str, str]) -> str:
    """Every gene carrying a curated assertion, escaped and comma-joined.

    Sorted by HGNC id, the same order `genes/index.json` publishes in, so this
    sentence and that file agree on which gene is named first without either
    reading the other. `symbols.get(gene, gene)` is the same registry fallback
    `build_genes` and `search.py` already apply to a gene absent from
    `mirrors/genes.tsv` — a symbol that has not been mirrored yet still reads
    as its HGNC id rather than as nothing.
    """
    genes = sorted({assertion.gene for assertion in corpus.assertions})
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
    asserted = _asserted_genes(corpus, symbols)

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>CHD Atlas</title>
<style>
  :root {{
    color-scheme: light dark;
    --fg: #1a1a1a;
    --bg: #ffffff;
    --muted: #555555;
    --border: #d0d0d0;
    --notice-bg: #fff4e5;
    --notice-border: #b45309;
    --link: #0b5fa5;
  }}
  @media (prefers-color-scheme: dark) {{
    :root {{
      --fg: #e8e8e8;
      --bg: #14161a;
      --muted: #a8a8a8;
      --border: #3a3d42;
      --notice-bg: #3a2a12;
      --notice-border: #d99a3d;
      --link: #6cb6f5;
    }}
  }}
  * {{ box-sizing: border-box; }}
  body {{
    margin: 0;
    padding: 2rem 1.25rem 4rem;
    background: var(--bg);
    color: var(--fg);
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial,
      sans-serif;
    line-height: 1.5;
  }}
  main {{ max-width: 42rem; margin: 0 auto; }}
  h1 {{ margin-bottom: 0.25rem; }}
  .tagline {{ color: var(--muted); margin-top: 0; }}
  h2 {{ margin-top: 2.25rem; border-bottom: 1px solid var(--border); padding-bottom: 0.25rem; }}
  .notice {{
    border: 1px solid var(--notice-border);
    background: var(--notice-bg);
    border-radius: 6px;
    padding: 1rem 1.25rem;
    margin: 1.5rem 0;
  }}
  .notice h2 {{ margin: 0 0 0.5rem; border: none; padding: 0; font-size: 1.05rem; }}
  dl {{ display: grid; grid-template-columns: auto 1fr; gap: 0.35rem 1rem; margin: 0; }}
  dt {{ color: var(--muted); }}
  dd {{ margin: 0; }}
  a {{ color: var(--link); }}
  code {{ background: rgba(127, 127, 127, 0.15); padding: 0.1rem 0.3rem; border-radius: 3px; }}
  ul.links {{ list-style: none; padding: 0; margin: 0; }}
  ul.links li {{ margin: 0.35rem 0; }}
  footer {{ margin-top: 3rem; color: var(--muted); font-size: 0.9rem; }}
</style>
</head>
<body>
<main>
  <h1>CHD Atlas</h1>
  <p class="tagline">Curated evidence linking genes, variants and proteins to congenital
    heart disease.</p>

  <section class="notice">
    <h2>Development status &amp; research use</h2>
    <p>This atlas is under active development. It is
      <strong>not a clinical decision-support tool</strong> and must not be used to make
      or guide a diagnostic, treatment or any other clinical decision.</p>
  </section>

  <h2>What this is</h2>
  <p>A curator reads the primary literature and mirrors upstream reference sources —
    ClinGen, GenCC, HPO — to build a static, queryable record of the evidence behind a
    gene's association with congenital heart disease. There is no server: every page
    below is a JSON file a client fetches and filters itself.</p>

  <h2>What's published</h2>
  <dl>
    <dt>Curated gene-disease assertions</dt>
    <dd>{assertion_count}</dd>
    <dt>Genes published</dt>
    <dd>{published_gene_count}</dd>
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
    <dt>Genes with mirrored validity in CHD scope (browsable once ClinGen grades it definitive)</dt>
    <dd>{mirrored_gene_count}</dd>
  </dl>

  <h2>Browse the data</h2>
  <ul class="links">
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
  </footer>
</main>
</body>
</html>
"""


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
