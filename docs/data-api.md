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
of this document; it is the one consumer obligation this API imposes, and
`atlas_curation` extends it to the genes this atlas has not yet curated — 22 of
the 23 published today.

---

## `manifest.json`

What the build produced, and a checksum for every file in it.

```json
{
  "counts": {
    "assertions": 1, "datasets": 0, "featured": 1,
    "functional": 0, "phenotypes": 3, "publications": 1
  },
  "files": {
    "genes/index.json": "sha256:<64 hex>",
    "publications.json": "sha256:<64 hex>"
  },
  "schema_version": "2.2",
  "source_commit": "<40-hex commit sha, or null outside a git checkout>",
  "status": "in-development"
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
- `counts` counts curated records, which is not the same as files. A gene the
  atlas has not curated contributes to no count and still gets a bundle — 22 of
  the 23 genes published today are in exactly that position, published on an
  expert panel's classification rather than on curation done here. Read
  `atlas_curation` on the browse row to tell the two apart.
- `schema_version` is `major.minor`. **Minor** rises when a field is added and
  nothing is removed or repurposed, so a consumer written against an earlier
  version keeps working; **major** rises when a field changes shape or leaves.
  `1.1` added `genes` to omics shard rows and `conflicting_lesion_groups` to
  gene index rows. `2.0` removed `classification` and `source_tier` from the
  curated assertion — the atlas no longer authors a gene-disease validity call
  of its own — and added the gene bundle's `validity` object in their place.
  The removal is what makes it major: a 1.x reader looking for a
  classification on the assertion now finds none. `2.1` added `status`,
  purely additively. `2.2` added `atlas_curation` to every gene index row and
  every gene bundle, also additively — in the release that widened the
  published gene set from the genes this atlas has curated to the genes a
  ClinGen expert panel calls definitive. More rows of an unchanged shape is not
  a schema change; the new field is, and it is what tells the two kinds of row
  apart.
- `status` is the atlas's own readiness, so a program can read it without
  scraping `index.html`'s prose. Today it is always `"in-development"` — one
  curated gene-disease assertion alongside mirrored ClinGen/GenCC validity for
  many more genes than that. This is a research resource, not a clinical
  decision-support tool; see `index.html` at the site root for the statement
  in full, and the note on [contested genes](#contested-genes-the-one-consumer-obligation)
  below for what the mirrored validity fields do and do not assert.

## The site root: `index.html`

The page a person opens directly rather than fetches as JSON, and the entry
point to the other 24. It states what the atlas is, the same development-status
and research-use statement `status` above is the machine-readable half of, and
the real counts behind it — curated assertions, genes with mirrored
ClinGen/GenCC validity, and the rest of `counts` — read from the same build
that produced this document's other examples rather than written by hand. It
links to `genes/index.html`, `genes/index.json`, `manifest.json`,
`sources.json` and the repository. Self-contained: no external request, no
build timestamp, byte-identical between two builds of one commit like
everything else here.

**Every page is checksummed in `manifest.json` exactly like a payload.** The
build behind this document publishes 25 of them — this one, `genes/index.html`
and one per published gene — and each has an entry in `files` giving the sha256
of the bytes served at it. A page is published output, so it is verifiable
output.

## `genes/index.json`

The browse payload. Downloaded by every visitor before they pick a gene, so it
carries what ranks or filters a row and the path to fetch the rest — and none of
the evidence itself.

**Which genes are listed.** One row per gene a ClinGen expert panel classifies
`Definitive` for a disease in this atlas's CHD scope (`curation/chd_scope.yaml`)
— 23 genes today. That is a mirrored decision, not a curated one: 22 of the 23
carry no assertion authored here, and `atlas_curation` is the field that says
which. GenCC is not a route in: it aggregates submissions rather than
adjudicating them, and the five in-scope genes it alone calls definitive
include one whose submissions run from `Definitive` to `No Known Disease
Relationship`.

**Read `atlas_curation` before presenting a row as curated content.** A
`headline_confidence` of `definitive` is an upstream expert panel's call and
never this atlas's, and on 22 of the 23 rows there is no assessment by this
atlas behind it at all.

```json
{
  "genes": [
    {
      "assertion_count": 1,
      "atlas_curation": "curated",
      "bundle": "genes/HGNC_11604.json",
      "confidence_by_lesion_group": { "septal": "definitive" },
      "conflicting_lesion_groups": [],
      "evidence_counts": { "genetic_case": 1 },
      "functional_count": 0,
      "gene": "HGNC:11604",
      "has_conflicting_evidence": false,
      "has_source_discordance": false,
      "headline_confidence": "definitive",
      "lesion_groups": ["septal"],
      "symbol": "TBX5",
      "validity_state": "expert_curated",
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
- `headline_confidence` and the rest of this row come from the mirrored
  ClinGen/GenCC validity records for the gene, never from a curated assertion —
  this atlas mirrors gene-disease validity, it does not author it. The example
  above is TBX5's real row: an in-scope ClinGen Congenital Heart Disease Gene
  Curation Expert Panel record makes it `"expert_curated"` with a `"definitive"`
  headline. `headline_confidence` is `null` for a gene no authority has
  assessed. It is never `"no_known_association"` for that case: that
  classification is itself an assessed verdict ("a panel looked and found
  nothing"), and asserting it for a gene nobody has assessed would state a
  conclusion no authority reached.
- `atlas_curation` is `"curated"` when the atlas holds at least one curated
  assertion for the gene and `"not_yet_curated"` when it holds none — 22 of the
  23 rows published today. It is about *this* resource's own work and says
  nothing about the gene's validity: an uncurated row still carries a
  `"definitive"` `headline_confidence`, because that is an expert panel's call
  and is exactly why the gene is listed. Read this rather than testing
  `assertion_count == 0`, so a browse filter does not have to reimplement the
  rule. It appears on the browse row and in the bundle, written in one place so
  the two cannot disagree.
- `validity_state` says how well curated the gene is: `"expert_curated"` (an
  in-scope ClinGen record exists), `"submitter_curated"` (only GenCC has
  assessed it) or `"uncurated"` (neither mirror has). `headline_confidence` is
  `null` **iff** the gene is `"uncurated"` *or* every in-scope record maps to
  no rung on this atlas's scale — the two are not the same condition, and
  `validity_state` is how a consumer tells them apart. Six genes in the
  committed mirrors (HGNC:24595, HGNC:4317, HGNC:6188, HGNC:7881, HGNC:9380,
  HGNC:9381) publish `null` while `"submitter_curated"`: each carries only an
  Orphanet `Supportive` submission, which `vocab.GENCC_CLASSIFICATIONS` maps
  to `None` because the submitter asserted an association without grading its
  evidence, not because nobody looked. Do not infer "no authority has assessed
  this gene" from `headline_confidence: null` alone — check `validity_state`
  for that.
- `has_source_discordance` is `true` when one mirrored source contests the gene
  while the *other* supports it. It is narrower than `has_conflicting_evidence`,
  which also fires when a single source is internally split across diseases or
  panels — see [Contested genes](#contested-genes-the-one-consumer-obligation).
- `confidence_by_lesion_group` applies the gene's mirrored `headline_confidence`
  to every lesion group its curated assertions declare — ClinGen and GenCC
  classify a gene against a disease, never against a specific lesion, so there
  is no finer-grained signal to divide the groups with. It is empty exactly
  when `headline_confidence` is `null`.
- `conflicting_lesion_groups` names every group in `confidence_by_lesion_group`
  when the gene is contested, and none when it is not — the same "no per-group
  signal" reasoning applies, so this can never name a proper subset of the
  gene's declared groups. It is the per-group counterpart of
  `has_conflicting_evidence` — see [Contested genes](#contested-genes-the-one-consumer-obligation),
  which is the one obligation this API places on a consumer.
- **The two lesion-group fields appear here and nowhere else.** The gene bundle
  carries `headline_confidence`, `validity_state`, `has_conflicting_evidence` and
  `has_source_discordance`, but neither of them, so a detail page that needs
  group-level confidence must carry it over from the browse row it was opened
  from. It cannot be derived from the bundle: that would mean reimplementing
  the classification ranking and the contested test, neither of which is
  published.
- The three counts describe what the bundle contains, so a browse row never
  promises more than the page delivers.

## `genes/index.html`

The browse page: the same 23 rows `genes/index.json` publishes, rendered as a
table a person can read and filter. Each row carries the HGNC id — linked to
that gene's page — the symbol, `headline_confidence`, `validity_state`,
`atlas_curation` and the gene's lesion groups. Above the table sit a text box
matching id or symbol and four menus (lesion group, confidence, validity state,
atlas curation), whose options are the values actually present in the build
rather than every value the vocabulary allows.

**Every row is rendered by the build, and the inline script only hides rows.**
There is no empty `<tbody>` filled in by a fetch, so `curl`, a crawler and a
reader with JavaScript disabled all get the complete table of 23 genes,
unfiltered — which is what the page shows before anyone touches a control in
any case. Nothing is loaded from anywhere: the stylesheet and the script are
inline, and the page makes no external request.

**`genes/index.json` is the machine-readable contract; this page is not.** The
column set, the markup and the `data-` attributes the filter reads may change
without a `schema_version` bump. Write against the payload and leave the page
to people.

## `genes/<slug>.json`

One gene's whole detail page, in one fetch.

```json
{
  "gene": "HGNC:11604",
  "symbol": "TBX5",
  "headline_confidence": "definitive",
  "validity_state": "expert_curated",
  "atlas_curation": "curated",
  "has_conflicting_evidence": false,
  "has_source_discordance": false,
  "lesion_groups": ["septal"],
  "validity": { "state": "expert_curated", "has_source_discordance": false, "records": [ … ] },
  "publications": ["PMID:8988165"],
  "assertions": [ { "id": "CHDA:AST:0000001", "lesion_groups": ["septal"], "evidence": [ … ] } ],
  "functional": [],
  "variants": [],
  "omics": {}
}
```

- `assertions` carry their full `evidence` array, including each item's
  `locator`, `strength` and `summary` — the record a curator is judged on.
  Assertion fields: `curated_on`, `curator`, `evidence`, `extracardiac_features`,
  `gene`, `id`, `inheritance`, `last_reviewed`, `lesion_groups`, `mechanism`,
  `notes`, `phenotypes`, `syndromic`. **No `classification` or `source_tier`**:
  the atlas mirrors gene-disease validity rather than curating its own, so an
  assertion says only what the curator is the authority for -- which lesions a
  gene is claimed for, and on what evidence. `headline_confidence` above is
  where the mirrored classification lives.
- `functional` holds **every** functional record about the gene, not only those
  an assertion cites.
- `atlas_curation` reads the same here as on the browse row. On the 22 genes
  published today without curation here it is `"not_yet_curated"`, and
  `assertions`, `publications` and `functional` are then empty arrays: the page
  is the panel's classification plus whatever this atlas has recorded, which
  may be nothing at all. An empty `assertions` array is a curation gap, never a
  fetch that failed.
- `omics` and `variants` are always present and may be empty. Read them without
  guarding for a missing key.
- **`variants` are embedded; omics rows are linked.** That asymmetry is a
  curation policy, not a property of the data — this atlas curates variants by
  hand, so the count per gene is bounded by effort. Omics tables are not, so a
  bundle carries per-modality summaries with `shards` to fetch. The omics section
  below gives their shape and how to select a gene's rows out of a fetched
  shard.
- `publications` lists the PMIDs the gene's assertion evidence cites, in lexical
  order. It does not include PMIDs cited only by its functional records.
- `assertions` and `functional` are ordered by id.

### The bundle's `validity` object: mirrored, attributed, never authored here

**This atlas publishes no gene-disease validity classification of its own.**
`headline_confidence`, `validity_state`, `has_conflicting_evidence` and
`has_source_discordance` at the top of the bundle come entirely from ClinGen
and GenCC, mirrored and attributed rather than asserted by a curator — see
[Contested genes](#contested-genes-the-one-consumer-obligation) for how to
read them safely. The bundle's `validity` object is where the mirrored records
behind those fields live:

- `state` repeats `validity_state`. `has_source_discordance` repeats the
  top-level field of the same name. Both are published twice on purpose: the
  flat fields are what the browse row (`genes/index.json`) and the bundle
  publish identically, and this object is the self-contained provenance
  record for a consumer that wants only "who curated this gene, and what did
  each of them say" without cross-referencing the fields beside it.
- `records` is one entry per in-scope mirrored classification, ClinGen's and
  GenCC's alike, in no particular order a consumer should rely on beyond what
  is published. Every record carries the same key set regardless of source:

  ```json
  {
    "source": "clingen",
    "classification_term": "Definitive",
    "classification": "definitive",
    "disease": "MONDO:0007732",
    "disease_label": "Holt-Oram syndrome",
    "moi": "AD",
    "sop": "SOP11",
    "classification_date": "2025-03-25T16:00:00.000Z",
    "gcep": "Syndromic Disorders Gene Curation Expert Panel",
    "report_url": "https://search.clinicalgenome.org/kb/gene-validity/CGGV:assertion_24e6c85a-33cf-4248-be1f-6431c7c6b1e5-2025-03-25T160000.000Z",
    "submitter": null
  }
  ```

  - `classification_term` is the authority's own word, verbatim — "Definitive",
    "Disputed Evidence", "Supportive". `classification` maps that term onto
    this atlas's `Classification` scale (the same values `headline_confidence`
    publishes), or `null` where the term is not a rung on it at all — GenCC's
    `Supportive`, a submitter asserting an association without grading its
    strength. A consumer that wants to render exactly what the authority said
    reads `classification_term`; one that wants to filter or rank reads
    `classification` and must handle `null`.
  - `sop`, `classification_date` and `gcep` are populated only on a ClinGen
    record; `submitter` only on a GenCC record. The field a record's own
    source does not carry is published as `null` rather than omitted, so
    every object in `records` has the same shape and a consumer never has to
    check `source` before it can look a key up.
  - `sop` is published because ClinGen's committed mirror spans SOP4 through
    SOP12 with no crosswalk between framework versions published anywhere. A
    classification attributed to ClinGen without its SOP version is an
    unqualified claim — the same rule under a different framework applies
    equally to a rung on the classification ladder itself.
  - `moi` and `disease`/`disease_label` are the authority's own mode of
    inheritance and disease term for that specific record — ClinGen and GenCC
    do not always agree on either, which is part of why two records for one
    gene can differ.
- `state` is one of three values, the same ones `validity_state` publishes:
  `"expert_curated"` (at least one in-scope ClinGen record exists),
  `"submitter_curated"` (only GenCC has assessed the gene) or `"uncurated"`
  (neither mirror has). An uncurated gene publishes `records: []` and
  `has_source_discordance: false` — the explicit empty shape, not an absent
  `validity` key, so a consumer can tell "no authority has assessed this
  gene" from "the field is missing" without guessing.

`sources.json` (below) carries the licence terms this atlas mirrors ClinGen
and GenCC under, the same way it does for HPO.

## `genes/<slug>.html`

One page per published gene, rendering that gene's bundle for a reader. A
summary column carries the headline classification as a chip — with a
`conflicting evidence` and a `sources disagree` chip beside it when those flags
are set — then `validity_state`, `atlas_curation`, the lesion groups, the
assertion, functional-record and publication counts, and a link to the gene's
own JSON. Beside it, the mirrored validity table: one row per record, giving
the source, the panel or submitter, the disease, the mode of inheritance, the
authority's own `classification_term` **verbatim** rather than the rung this
atlas maps it onto, the SOP, the date, and a link to the upstream report where
one is published.

Then either the curated evidence — each assertion, its evidence items with
their class, strength and summary, and the publications that evidence cites —
or, for the 22 genes published today with no curation here, a paragraph saying
exactly that:

> The atlas has **not yet curated** a lesion assertion for this gene. The
> classification above is an expert panel's, mirrored with its provenance
> intact; no classification on this page is the atlas's own assessment.

That paragraph is there instead of the section simply being absent. A missing
evidence section is indistinguishable from "the atlas looked and found
nothing", and a reader deciding what a gene means clinically must not have to
infer which.

Read both halves of that sentence precisely, because each is narrower than it
first appears. **A lesion assertion** is what is absent, not evidence in
general: `atlas_curation` is derived from curated `LesionAssertion` records
alone, so a gene can carry functional-evidence records the atlas curated and
still report `not_yet_curated`. Such a page adds a second paragraph naming
those records as the atlas's own work, rather than leaving the first to deny
them. And **no classification** is the atlas's own — not "nothing on this
page", which would have denied that same curated work one column away from the
rail counting it.

A gene page carries no script at all — every value on it was rendered at build
time — and, like the browse page, makes no external request. The bundle remains
the machine-readable contract: the page's markup is not one.

**No payload carries the path to a gene's page.** A browse row's `bundle` is
the JSON; the route to the page is the link on `genes/index.html`, which is why
the slug rule stays an implementation detail (see
[Reading this API](#reading-this-api)).

## `omics/<modality>/<accession>.json`

A gene bundle's `omics` maps a modality — `expression`, `profiles`, `proteomics`
or `phospho` — to a summary of that gene's rows:

```json
{
  "omics": {
    "expression": {
      "count": 412,
      "shards": ["omics/expression/GSE1000.json"],
      "top": [
        { "dataset": "GSE1000", "gene": "HGNC:11604", "log2fc": 2.1, "fdr": 0.001,
          "genes": ["HGNC:11604"] }
      ]
    }
  }
}
```

- `count` is every row about the gene, across every shard listed.
- `top` holds the same row objects the shard does, `genes` included.
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
      "blurb": "One of two back-to-back 1997 reports identifying TBX5 mutations …",
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
a gene list. It is `null` for a term that is not itself a cardiac lesion — an
extracardiac feature (e.g. a limb malformation cited by a syndromic assertion)
registered here only so its label is checked against the pinned HPO release,
the same guarantee every cardiac term gets.

## `datasets.json`

```json
{ "datasets": [] }
```

One record per omics dataset: accession, archive, technology, tissue, stage,
organism, sample count, licence and its contrasts. This is what an omics row's
`dataset` column resolves against, the way `publications.json` resolves a PMID.
Empty in the committed corpus today.

## `sources.json`

What the atlas mirrors, and on whose terms.

```json
{
  "sources": [
    {
      "id": "hpo",
      "name": "Human Phenotype Ontology",
      "version": "hp/releases/2026-06-23",
      "retrieved_on": "2026-07-31",
      "url": "https://hpo.jax.org/",
      "licence": "https://hpo.jax.org/app/license",
      "redistribution": "permitted_with_attribution",
      "ontology_prefix": "HP",
      "ontology_file": "ontologies/hp-2026-06-23.obo"
    }
  ]
}
```

**Read this before redistributing anything from this site.** Phenotype labels
and synonyms in `phenotypes.json` and in the search index are transcribed from
the pinned HPO release, whose terms are `permitted_with_attribution` — so a
consumer republishing them carries the same obligation, and this file is where
the attribution to satisfy it comes from.

The repository's `LICENSE` (Apache-2.0) covers the **code**. It does not govern
mirrored third-party content, whose terms are the ones recorded here.

`version` is the upstream release identifier and `retrieved_on` the date it was
taken, so a claim can be traced to the exact release it rests on.

The registry also carries `clingen` and `gencc` — the two sources every
mirrored record in a gene bundle's `validity.records` is attributed to.
Both are recorded `"redistribution": "permitted"` (CC0-1.0), unlike HPO's
`"permitted_with_attribution"`: neither authority's licence requires an
attribution notice to redistribute their content, though GenCC's terms
*request* one. This file is where that distinction is recorded, rather than
assumed from the presence of a `validity` object.

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

`headline_confidence` is the strongest classification the mirrored ClinGen and
GenCC records assert for a gene, on a single linear scale where `definitive`
outranks `refuted`. ClinGen treats disputed and refuted as a **separate axis**
rather than weaker rungs of the same ladder, so a gene whose mirrored records
carry both a definitive and a refuted classification resolves to `definitive`
and the refutation is invisible in that field alone.

**`headline_confidence` is `null` for a gene no authority has assessed, and
must never be rendered as `"no_known_association"`.** The two are not
interchangeable: `no_known_association` is itself an assessed verdict — a
panel looked and found nothing. `null` is not that specific claim, but it is
not always "no panel looked" either: a GenCC submitter can assert an
association under a term this atlas maps to no rung at all (`Supportive`,
which declines to grade the evidence) and still publish `headline_confidence:
null` — six genes in the committed mirrors do exactly this, each via a single
Orphanet `Supportive` submission (see `validity_state` above). `null` means
either nobody has assessed the gene, or nobody who did assessed it on a scale
this atlas can rank; check `validity_state` to tell the two apart. Coercing
either case's `null` to `"no_known_association"` would fabricate a conclusion
nobody reached. This atlas publishes no gene-disease validity classification
of its own; every value `headline_confidence` can take, including the absence
of one, is mirrored and attributed from ClinGen or GenCC — see the bundle's
`validity` object, documented under [`genes/<slug>.json`](#genesslugjson), for
the records behind it and `sources.json` for the licence terms those two
mirrors are republished under.

`has_conflicting_evidence` is the other half of that pair. It appears in both
the browse row and the bundle, and is always written alongside
`headline_confidence`. `has_source_discordance` is a narrower relative: it is
`true` only when the contesting and the supporting classification come from
*different* mirrored sources. A single source split against itself would set
`has_conflicting_evidence` without setting this one — that is why the two
fields are not redundant — but that split does not currently occur among the
154 genes these mirrors curate within CHD scope.

**No row in `genes/index.json` sets either flag today, and that is not a reason
to skip implementing them.** Measured against a real build (2026-08-04): of the
154 mirrored genes, exactly one — LEFTY2, HGNC:3122, where ClinGen's own
`Disputed` call sits alongside GenCC's supportive one — sets
`has_conflicting_evidence`, and it is the same gene that sets
`has_source_discordance`. It is not among the 23 published, because a gene is
published on a ClinGen `Definitive` call and LEFTY2 does not have one. So a
consumer testing their rendering against the live data will find nothing
contested to look at. The first gene ClinGen both grades definitive for one
in-scope disease and disputes for another will appear as an ordinary
`"definitive"` row to any client that did not implement this, which is the
failure this whole section exists to prevent. The divergence
`has_source_discordance` catches is likewise real but unrealised: ninety genes
in the full, pre-scope ClinGen mirror carry both a supportive and a contesting
call, but none of the diseases those ninety concern is in CHD scope
(`build/validity.py`'s `_has_source_discordance` docstring has the
measurement).

**A consumer must pair `headline_confidence` with `has_conflicting_evidence`
and present a contested gene distinctly** — a badge, a different colour, an
explicit note. Rendering `headline_confidence` alone would tell a reader that a
gene the field disputes is settled science, which is the one failure this
resource exists to prevent.

`confidence_by_lesion_group` is **not** a finer-grained view of the same
question. ClinGen and GenCC classify a gene against a disease, never against a
specific lesion, so the mirrored records carry no per-group information at
all — every lesion group a curated assertion names for the gene publishes the
*identical* `strongest()` of the gene's mirrored records. It differs from
`headline_confidence` only in shape, as a map over the gene's declared groups
for a consumer already filtering by lesion, never in value.

**`conflicting_lesion_groups` is `has_conflicting_evidence`'s per-group
counterpart**, and for the same reason it cannot single out which of a
contested gene's groups is the disputed one: it lists *every* group in
`confidence_by_lesion_group` when the gene is contested, and none when it is
not.

```json
{
  "has_conflicting_evidence": true,
  "confidence_by_lesion_group": { "conotruncal": "definitive", "septal": "definitive" },
  "conflicting_lesion_groups": ["conotruncal", "septal"]
}
```

Read together, those say: the gene is contested, and every lesion group it is
curated for inherits that contest equally. There is no mirrored signal that
could clear one group while leaving another disputed.

**Where these two fields live.** `confidence_by_lesion_group` and
`conflicting_lesion_groups` appear in `genes/index.json` and **nowhere else**.
The gene bundle carries `headline_confidence`, `validity_state`,
`has_conflicting_evidence` and `has_source_discordance`, but neither of these
two, so a detail page needing group-level confidence must carry it over from
the browse row it was opened from. It cannot be recovered from the bundle:
that would mean reimplementing the classification ranking and the contested
test against the mirrored records, and neither rule is published.

**The two states the pair can express:**

| `has_conflicting_evidence` | `conflicting_lesion_groups` | what it means |
| --- | --- | --- |
| `false` | `[]` | nothing about this gene is disputed |
| `true` | every group in `confidence_by_lesion_group` | the gene is disputed, and the dispute applies equally to every lesion group it is curated for |

A gene disputed about only *some* of its lesion groups is not a state this API
can express: the mirrors classify by disease, not by lesion, so there is
nothing in the source data to divide the groups on. `conflicting_lesion_groups`
is therefore always either every group the gene declares or none of them —
never a proper subset.

The list is always present and may be empty. Every group it names is a key of
`confidence_by_lesion_group`, so the two join directly.

### The same failure in a second shape: a gene the atlas has not curated

`headline_confidence` is published for every gene, curated here or not, so a
consumer that renders it alone presents 22 of the 23 published genes as this
atlas's assessment of them. It is not: it is a ClinGen expert panel's, and this
atlas has recorded nothing about those genes beyond republishing it with its
provenance.

**A consumer must not present a `not_yet_curated` gene as one the atlas has
assessed.** `atlas_curation` appears on the browse row and in the bundle,
written in one place so the two cannot disagree, and it is what tells the two
kinds of row apart — read it rather than testing `assertion_count == 0`, so a
filter does not have to reimplement the rule. An empty `assertions` array is a
curation gap, never a verdict and never a fetch that failed, the same way
`headline_confidence: null` is not `no_known_association`. This atlas's own
gene pages state the gap in prose for exactly that reason.
