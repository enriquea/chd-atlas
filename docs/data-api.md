# Published data API

The atlas is served as static JSON. There is no server, no query language and no
pagination: a consumer fetches whole files and filters them client-side.

Build it with `uv run chd-atlas build --root . --out dist`. The command refuses
to write anything if `chd-atlas validate` would report an error, so a published
site is always one that passed its own gate. Every example below is copied from
a real build of the committed corpus rather than written by hand.

Two builds of one commit are byte-identical. Nothing in the output carries a
timestamp — a consumer wanting a publication date should read the commit named
in `manifest.json`.

## Reading this API

**Never construct a path.** Every payload that refers to another file carries
the path to fetch. `genes/index.json` gives each gene its `bundle`, and the
search index gives each record its `path`. The slug rule that turns `HGNC:11604`
into `HGNC_11604` is an implementation detail and is not part of this contract.

**`.json.gz` is not transparently decompressed.** GitHub Pages serves a `.json.gz`
as `application/gzip` and sets no `Content-Encoding`, so the browser hands you
raw compressed bytes. Decompress them yourself:

```js
const response = await fetch("search/index.json.gz");
const stream = response.body.pipeThrough(new DecompressionStream("gzip"));
const { records } = await new Response(stream).json();
```

Plain `.json` files are gzipped in transit by Pages and need no such handling.

**A contested gene must never be displayed as settled.** See the note at the end
of this document; it is the one consumer obligation this API imposes.

---

## `manifest.json`

What the build produced, and a checksum for every file in it.

```json
{
  "counts": {
    "assertions": 1, "datasets": 0, "featured": 1,
    "functional": 0, "phenotypes": 2, "publications": 1
  },
  "files": {
    "genes/index.json": "sha256:<64 hex>",
    "publications.json": "sha256:<64 hex>"
  },
  "schema_version": "1.1",
  "source_commit": "<40-hex commit sha, or null outside a git checkout>"
}
```

The structure above is copied from a real build; the checksum and commit values
are shown as placeholders on purpose. Both are derived from content that changes
with every commit, so pinning real ones here would guarantee this document is
wrong by the next one.

- `files` maps a relative URL to the sha256 of **the bytes actually served** at
  it. For a `.json.gz` that is the digest of the compressed file, which is what
  a consumer can verify against what it downloaded.
- `manifest.json` is absent from its own `files`. It cannot be present: the
  value would have to be computed from bytes containing it.
- `source_commit` is `null` when the build was made outside a git checkout — an
  unpacked tarball still produces a complete site, just one that cannot state
  its provenance.
- `counts` counts curated records, which is not the same as files. A gene with
  no assertion contributes to no count and gets no bundle.
- `schema_version` is `major.minor`. **Minor** rises when a field is added and
  nothing is removed or repurposed, so a consumer written against an earlier
  version keeps working; **major** rises when a field changes shape or leaves.
  `1.1` added `genes` to omics shard rows and `conflicting_lesion_groups` to
  gene index rows.

## `genes/index.json`

The browse payload. Downloaded by every visitor before they pick a gene, so it
carries what ranks or filters a row and the path to fetch the rest — and none of
the evidence itself.

```json
{
  "genes": [
    {
      "assertion_count": 1,
      "bundle": "genes/HGNC_11604.json",
      "confidence_by_lesion_group": { "septal": "definitive" },
      "evidence_counts": { "genetic_case": 1 },
      "functional_count": 0,
      "gene": "HGNC:11604",
      "has_conflicting_evidence": false,
      "headline_confidence": "definitive",
      "lesion_groups": ["septal"],
      "symbol": "TBX5",
      "variant_count": 0
    }
  ]
}
```

- The array is ordered by HGNC id, not by symbol. JSON arrays keep their order,
  so this is part of the contract.
- `bundle` is the path to fetch for the detail page. Do not build it yourself.
- `symbol` falls back to the HGNC id for a gene not yet in `mirrors/genes.tsv`,
  so it is always a non-empty string you can render and search.
- `confidence_by_lesion_group` may disagree with `headline_confidence` by
  design: a gene can be definitive for septal disease and refuted for
  conotruncal.
- The three counts describe what the bundle contains, so a browse row never
  promises more than the page delivers.

## `genes/<slug>.json`

One gene's whole detail page, in one fetch.

```json
{
  "gene": "HGNC:11604",
  "symbol": "TBX5",
  "headline_confidence": "definitive",
  "has_conflicting_evidence": false,
  "lesion_groups": ["septal"],
  "publications": ["PMID:8988165"],
  "assertions": [ { "id": "CHDA:AST:0000001", "classification": "definitive", "evidence": [ … ] } ],
  "functional": [],
  "variants": [],
  "omics": {}
}
```

- `assertions` carry their full `evidence` array, including each item's
  `locator`, `strength` and `summary` — the record a curator is judged on.
  Assertion fields: `classification`, `curated_on`, `curator`, `evidence`,
  `extracardiac_features`, `gene`, `id`, `inheritance`, `last_reviewed`,
  `lesion_groups`, `mechanism`, `notes`, `phenotypes`, `source_tier`,
  `syndromic`.
- `functional` holds **every** functional record about the gene, not only those
  an assertion cites.
- `omics` and `variants` are always present and may be empty. Read them without
  guarding for a missing key.
- **`variants` are embedded; omics rows are linked.** That asymmetry is a
  curation policy, not a property of the data — this atlas curates variants by
  hand, so the count per gene is bounded by effort. Omics tables are not, so a
  bundle carries per-modality summaries with `shards` to fetch. See the omics
  section below for the shape of those summaries, and for a limitation on
  selecting a gene's rows out of a `phospho` shard.
- `publications` lists the PMIDs the gene's assertion evidence cites, in lexical
  order. It does not include PMIDs cited only by its functional records.
- `assertions` and `functional` are ordered by id.

## `omics/<modality>/<accession>.json`

A gene bundle's `omics` maps a modality — `expression`, `profiles`, `proteomics`
or `phospho` — to a summary of that gene's rows:

```json
{
  "omics": {
    "expression": {
      "count": 412,
      "shards": ["omics/expression/GSE1000.json"],
      "top": [ { "dataset": "GSE1000", "gene": "HGNC:11604", "log2fc": 2.1, "fdr": 0.001 } ]
    }
  }
}
```

- `count` is every row about the gene, across every shard listed.
- `shards` are the files holding them. Each shard is
  `{"table": "<modality>", "rows": [ … ]}` — the mirror rows, each with one field
  added by the build (see below).
- **`top` is capped at 25 rows.** It is a preview, ranked by significance, not a
  page of results. `count` is frequently larger, and the cap is not carried in
  the payload, so do not infer completeness from `len(top)`.

**To get the rows a bundle counted, filter the shard on `genes`.**

Every shard row carries `genes`, a list of the HGNC ids that row is evidence
about:

```json
{ "dataset": "PXD012345", "protein": "Q99593", "position": 100, "genes": ["HGNC:11604"] }
```

```js
const shard = await (await fetch(summary.shards[0])).json();
const mine = shard.rows.filter(row => row.genes.includes("HGNC:11604"));
// mine.length === summary.count, when the gene has one shard
```

That equality is the point of the field, and it holds by construction: `count` is
derived from the same attribution the rows publish, computed once. It is a
**list** because one protein accession can belong to several genes — a histone
cluster, for instance — and a single-valued field would silently drop all but one.

`genes` is present on every modality, including `expression` and `profiles` whose
rows already carry their own `gene` column, so a consumer filters one way
everywhere. Where the two exist side by side they agree; `genes` is the one
`count` is built from.

`mirrors/ptm_sites.tsv` is a validated mirror table that this site does not
publish at all. It is reference data about modification sites rather than
evidence about a gene, and no bundle links it.

## `variants/index.json` and `variants/<chrom>.json.gz`

```json
{ "shards": ["variants/1.json.gz", "variants/X.json.gz"] }
```

- The index exists so a consumer can enumerate chromosomes without probing for
  404s, and is emitted even when empty — as it is in the committed corpus today.
- `shards` is in karyotype order (1…22, X, Y, MT), not lexical, so it can drive
  a chromosome picker directly.
- Each shard is `{"chrom": "12", "rows": [ … ]}` and holds only rows on the
  chromosome it is named for.
- Gzipped: see the decompression note above.

## `publications.json`

```json
{
  "publications": [
    {
      "id": "PMID:8988165",
      "title": "Mutations in human TBX5 [corrected] cause limb and cardiac malformation in Holt-Oram syndrome.",
      "journal": "Nature genetics",
      "year": 1997,
      "authors": ["Basson CT", "Bachinsky DR", "…"],
      "study_type": "family_linkage",
      "doi": "10.1038/ng0197-30",
      "pmcid": null,
      "own_lab": false,
      "cohort_size": null,
      "ancestry": []
    }
  ]
}
```

Ordered by PMID lexically, not numerically — so `PMID:10` precedes `PMID:9`. A
PMID is issued at indexing time, so numeric order ranks by nothing a reader
asked for; pages should rank by year or by the curated featured list.

## `featured.json`

The landing page's manuscripts, in curator-chosen `order`.

```json
{
  "featured": [
    {
      "order": 1,
      "topic": "…",
      "blurb": "The founding demonstration that TBX5 haploinsufficiency causes …",
      "publication": { "id": "PMID:8988165", "title": "…", "journal": "Nature genetics", "…": "…" }
    }
  ]
}
```

`publication` is the **resolved object**, not a PMID string, so the landing page
renders without a second fetch. It is never a bare string.

## `phenotypes.json`

```json
{
  "phenotypes": [
    { "id": "HP:0001629", "label": "Ventricular septal defect",
      "lesion_group": "septal", "synonyms": ["VSD"] }
  ]
}
```

`lesion_group` is the facet the gene index's `lesion_groups` and
`confidence_by_lesion_group` key on, which is what lets a phenotype filter drive
a gene list.

## `datasets.json`

```json
{ "datasets": [] }
```

One record per omics dataset: accession, archive, technology, tissue, stage,
organism, sample count, licence and its contrasts. This is what an omics row's
`dataset` column resolves against, the way `publications.json` resolves a PMID.
Empty in the committed corpus today.

## `search/index.json.gz`

A flat array of records over genes, publications and phenotypes. Deliberately
not an inverted index: at this corpus size a client filters the whole array in a
fraction of a frame.

```json
{
  "records": [
    { "kind": "gene", "id": "HGNC:11604", "label": "TBX5",
      "path": "genes/HGNC_11604.json",
      "terms": ["TBX5", "HGNC:11604", "T-box transcription factor 5"] },
    { "kind": "phenotype", "id": "HP:0001631", "label": "Atrial septal defect",
      "path": "phenotypes.json",
      "terms": ["Atrial septal defect", "HP:0001631", "ASD"] }
  ]
}
```

- `terms` is the haystack: the strings a visitor might type, deduplicated,
  including each record's own identifier. **Matching is the client's job** —
  this file ships no scoring, no stemming and no ranking.
- `label` is what a result row displays; `id` identifies the thing.
- `path` is the payload that answers the query.
- Variants and datasets are not indexed. The variant space grows without bound
  and would dominate the size of the file every visitor downloads; a dataset has
  no title or description of its own, so its accession is the only string naming
  it, and that already resolves through `datasets.json`.
- Genes come from the assertion set, so a gene with omics or variant evidence
  but no curated assertion is not searchable — the atlas browses curated claims.

---

## Contested genes: the one consumer obligation

`headline_confidence` is the strongest classification asserted for a gene on a
single linear scale, where `definitive` outranks `refuted`. ClinGen treats
disputed and refuted as a **separate axis** rather than weaker rungs of the same
ladder, so a gene carrying both a definitive and a refuted assertion resolves to
`definitive` and the refutation is invisible in that field alone.

`has_conflicting_evidence` is the other half of that pair. It appears in both
the browse row and the bundle, and is always written alongside
`headline_confidence`.

**A consumer must pair the two and present a contested gene distinctly** — a
badge, a different colour, an explicit note. Rendering `headline_confidence`
alone would tell a reader that a gene the field disputes is settled science,
which is the one failure this resource exists to prevent.

`confidence_by_lesion_group` is the finer-grained view, and it collapses the same
way: each group's value is `strongest()` over that group's classifications, so a
group that is both `definitive` and `refuted` resolves to `definitive` exactly as
the headline does.

**`conflicting_lesion_groups` is its flag**, and the same obligation applies one
level down. It lists the groups that are themselves contested:

```json
{
  "has_conflicting_evidence": true,
  "confidence_by_lesion_group": { "septal": "definitive", "conotruncal": "moderate" },
  "conflicting_lesion_groups": ["septal"]
}
```

Read together, those say: the gene is contested, and specifically about septal
disease — the conotruncal association is not in dispute. `has_conflicting_evidence`
alone cannot make that distinction, since it is computed over the gene's whole
classification set.

The list is always present and may be empty. Every group it names is a key of
`confidence_by_lesion_group`, so the two join directly.
