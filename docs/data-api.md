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
    "functional": 0, "phenotypes": 3, "publications": 1
  },
  "files": {
    "genes/index.json": "sha256:<64 hex>",
    "publications.json": "sha256:<64 hex>"
  },
  "schema_version": "2.1",
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
- `counts` counts curated records, which is not the same as files. A gene with
  no assertion contributes to no count and gets no bundle.
- `schema_version` is `major.minor`. **Minor** rises when a field is added and
  nothing is removed or repurposed, so a consumer written against an earlier
  version keeps working; **major** rises when a field changes shape or leaves.
  `1.1` added `genes` to omics shard rows and `conflicting_lesion_groups` to
  gene index rows. `2.0` removed `classification` and `source_tier` from the
  curated assertion — the atlas no longer authors a gene-disease validity call
  of its own — and added the gene bundle's `validity` object in their place.
  The removal is what makes it major: a 1.x reader looking for a
  classification on the assertion now finds none. `2.1` added `status`,
  purely additively.
- `status` is the atlas's own readiness, so a program can read it without
  scraping `index.html`'s prose. Today it is always `"in-development"` — one
  curated gene-disease assertion alongside mirrored ClinGen/GenCC validity for
  many more genes than that. This is a research resource, not a clinical
  decision-support tool; see `index.html` at the site root for the statement
  in full, and the note on [contested genes](#contested-genes-the-one-consumer-obligation)
  below for what the mirrored validity fields do and do not assert.

## The site root: `index.html`

The one page a person opens directly rather than fetches as JSON. It states
what the atlas is, the same development-status and research-use statement
`status` above is the machine-readable half of, and the real counts behind it
— curated assertions, genes with mirrored ClinGen/GenCC validity, and the
rest of `counts` — read from the same build that produced this document's
other examples rather than written by hand. It links to `manifest.json`,
`genes/index.json`, `sources.json` and the repository. Self-contained: no
external request, no build timestamp, byte-identical between two builds of one
commit like everything else here.

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

## `genes/<slug>.json`

One gene's whole detail page, in one fetch.

```json
{
  "gene": "HGNC:11604",
  "symbol": "TBX5",
  "headline_confidence": "definitive",
  "validity_state": "expert_curated",
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
154 genes these mirrors curate within CHD scope. Measured against a real
build: exactly one gene sets `has_conflicting_evidence` — LEFTY2 (HGNC:3122),
where ClinGen's own `Disputed` call sits alongside GenCC's supportive one —
and it is the same gene that sets `has_source_discordance`. The divergence
this field exists to catch is real but currently unrealised in published
data: ninety genes in the full, pre-scope ClinGen mirror carry both a
supportive and a contesting call, but none of the diseases those ninety
concern is in CHD scope (`build/validity.py`'s `_has_source_discordance`
docstring has the measurement).

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
