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
`atlas_curation` extends it to the genes this atlas has not yet curated — 91 of
the 92 published today.

---

## `manifest.json`

What the build produced, and a checksum for every file in it.

```json
{
  "counts": {
    "assertions": 1, "burden_rows": 915, "cohort_families": 3,
    "cohorts": 13, "datasets": 1, "featured": 1, "functional": 0,
    "genes": 92, "phenotypes": 3, "profile_datasets": 1, "profile_genes": 92,
    "publications": 5
  },
  "files": {
    "genes/index.json": "sha256:<64 hex>",
    "publications.json": "sha256:<64 hex>"
  },
  "schema_version": "2.12",
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
- `counts` is a census of **what this build published**, which is not the same as
  a count of files and — since `2.7` — no longer only a count of curated records.
  Six of its keys count a curated collection: `assertions`, `datasets` (omics),
  `featured`, `functional`, `phenotypes`, `publications`. A gene the atlas has
  not curated contributes to none of them and still gets a bundle — 91 of the 92
  genes published today are in exactly that position, published on an external
  authority's classification rather than on curation done here. Read
  `atlas_curation` on the browse row to tell the two apart.

  Three keys count the build instead, and were added because the object read as
  a census of the atlas while describing only its curation: through `2.6` it
  published `assertions: 1, datasets: 0, functional: 0` for a site serving 23
  genes and 290 burden statistics — the corpus as it stood that day — and named
  neither.

  - `genes` is the published population — the number of gene bundles, and the
    length of `genes/index.json`'s `genes` array.
  - `burden_rows` is the number of burden statistics reaching a bundle, which is
    exactly the sum of every browse row's `burden_row_count`. It counts what a
    consumer can fetch, **not** what `mirrors/burden.tsv` holds: that file is
    deliberately wider than the publication gate, and today 59 of its 150 genes
    publish no page.
  - `cohort_families` is the number of independent cohort families, the length
    of any bundle's `independent_datasets.families` array — the same for every
    gene by construction. It is deliberately not named `independent_datasets`,
    because `datasets` above already means omics datasets and two adjacent keys
    reading `datasets: 0` and `independent_datasets: 3` invite a misreading.

  `index.html` states the same three figures in prose. They are derived once in
  the build and handed to both, so the page and this object cannot disagree.
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
  apart. `2.3` added `burden` to every gene bundle and `burden_row_count` to
  every gene index row, again additively; both keys are always present, so a
  2.2 reader keeps working and a 2.3 reader never guards for a missing one.
  `2.4` added `pvalue_adjusted` and `pvalue_adjustment` to every burden object,
  and `2.5` added `count_unit`. Both additive, both always present. Both carry a
  *display* obligation heavier than the parsing one, described with the
  [`burden` array](#the-bundles-burden-array-per-study-never-pooled) below: a raw
  and a corrected p can point opposite ways, and two rows counting different
  units are not comparable at all. `2.6` added `independent_datasets`, also
  additive and always present, and its display obligation is the heaviest of the
  three: it is a **count of datasets, not a validity call**, and rendering it as
  a verdict beside a mirrored ClinGen `definitive` tells a clinician something
  the data do not say. See
  [`independent_datasets`](#independent_datasets-a-count-of-datasets-never-a-verdict).
  `2.7` added `genes`, `burden_rows` and `cohort_families` to `counts`, described
  above. Additive: every `2.6` key is present and unchanged, so the only thing a
  `2.6` reader misses is the half of the census that was never there.
  `2.8` added [`admitted_by` and `asserted_by`](#admitted_by-and-asserted_by-why-this-gene-is-here)
  to every gene bundle, and **widened the population from 23 genes to 92**. Both
  keys are additive and always present, so a 2.7 parser keeps working — but the
  population change carries the heaviest display obligation on this list, and it
  is not additive. Until 2.8 every published gene was `definitive`, so a
  consumer could render `headline_confidence` verbatim and be right every time.
  That is now 23 `definitive`, 1 `strong`, 9 `moderate`, 43 `limited` and 16
  `null`. A consumer that renders the chip without the grade's meaning turns
  `limited` — a panel saying the case is **not yet made** — into a weak yes, and
  one that treats presence in `genes/index.json` as "a panel called this a CHD
  gene" is now wrong for 16 genes. This is a population change inside an
  unchanged shape, which is why it is MINOR by the rule above; the obligation it
  creates is real regardless of the version letter.
  `2.10` added `has_no_association_report` and `no_association_reported_by` to
  every gene bundle and browse row, described above. Additive, and a third axis
  rather than a widening of `has_conflicting_evidence`.
  `2.9` added [`cohorts.json`](#cohortsjson) and a `cohorts` count. Additive: no
  existing key changes meaning and no payload loses one. It closes a gap open
  since `2.3` — every burden row has named its sample collections by bare id
  since `burden` was published, and nothing resolved those ids to anything, so a
  consumer had the numbers and none of the caveats that qualify them.
  `2.11` added `profile_datasets` and `profile_genes` to `counts`, restricted
  to `published` for the same reason `genes` and `burden_rows` already are,
  plus a "Developmental expression" pair of cards on `index.html`. Additive:
  both keys are always present, and both were `0` when the version was minted,
  because no `profiles` mirror had been committed yet. `E-MTAB-6814` has
  landed since, and the example above carries the figures a build produces
  today. See
  [`expression_profile`](#the-bundles-expression_profile-object-a-developmental-transcriptome-never-a-contrast)
  for the whole layer these two counts summarise.
  `2.12` added `order` to every entry of a dataset record's
  [`stages`](#datasetsjson) array, and **changed the order of the
  `expression_profile.datasets[].stages[]` array in every gene bundle from
  alphabetical to chronological**. The added key is additive and always
  present, so a 2.11 parser keeps working; the reordering is a correction to
  data that was already published, and it is the release. Through `2.11` a
  gene's developmental series came back in dictionary order — `4 week post
  conception` eighth, *after* `19 week post conception`, and `elderly` second
  of the eight post-natal stages. It was deterministic and reproducible and it
  was not a chronology. **If you plotted `stages` in array order you drew an
  axis in the wrong sequence, and you will now draw a different picture; that
  is the fix arriving, not a regression.** If you worked around it by sorting
  the stage tokens yourself, remove the workaround — sorting the tokens is
  exactly what the bug did. Like `2.2`'s and `2.8`'s population changes this is
  MINOR by the rule above, because no field is added, removed or reshaped and
  every entry still carries every key `2.11` published; the display obligation
  is real regardless of the letter.
- `status` is the atlas's own readiness, so a program can read it without
  scraping `index.html`'s prose. Today it is always `"in-development"` — one
  curated gene-disease assertion alongside mirrored ClinGen/GenCC validity for
  many more genes than that. This is a research resource, not a clinical
  decision-support tool; see `index.html` at the site root for the statement
  in full, and the note on [contested genes](#contested-genes-the-one-consumer-obligation)
  below for what the mirrored validity fields do and do not assert.

## The site root: `index.html`

The page a person opens directly rather than fetches as JSON, and the entry
point to the other 93. It states what the atlas is, the same development-status
and research-use statement `status` above is the machine-readable half of, and
the real counts behind it — every key of `counts`, plus three figures that are
not manifest keys: the genes carrying mirrored ClinGen/GenCC validity, the genes
carrying burden evidence, and the genes this atlas has itself curated. All are
read from the same build that produced this document's other examples rather
than written by hand. It links to
`genes/index.html`, `genes/index.json`, `manifest.json`, `sources.json` and the
repository. Self-contained: no external request, no build timestamp,
byte-identical between two builds of one commit like everything else here.

It also carries a key to the four evidence glyphs, built from the same constant
as the browse page's strip legend and the gene page's matrix legend. The three
cannot come to describe the same four states in different words.

**Every page is checksummed in `manifest.json` exactly like a payload.** The
build behind this document publishes 94 of them — this one, `genes/index.html`
and one per published gene — and each has an entry in `files` giving the sha256
of the bytes served at it. A page is published output, so it is verifiable
output.

## `genes/index.json`

The browse payload. Downloaded by every visitor before they pick a gene, so it
carries what ranks or filters a row and the path to fetch the rest — and none of
the evidence itself.

**Which genes are listed.** One row per gene an external authority already
treats as a congenital heart disease gene — 92 genes today, on either of two
warrants:

- a **ClinGen** record classifying it `Limited` or better for a disease that an
  external authority treats as congenital heart disease (76 genes), or
- **two or more Gene Curation Coalition submitters** independently asserting it
  at `Limited` or better, where no ClinGen panel disputes it (16 genes).

**No disease is in scope on this atlas's own judgement, and neither is any
gene.** `curation/chd_scope.yaml` records, for every term, which authority uses
it, and a validator (`SCP005`) checks that attribution against the mirrors
rather than taking it on trust. Every bundle carries `admitted_by`, naming the
one warrant that admitted that gene. That makes the whole chain a mirrored
decision rather than a curated one: 91 of the 92 genes carry no assertion
authored here, and `atlas_curation` is the field that says which.

**The two warrants are not symmetric, deliberately.** ClinGen admits alone
because it is the only source here with chartered expert panels and published
SOPs; GenCC harmonises submissions rather than adjudicating between them, so it
needs two. ClinGen's own GenCC submissions are excluded from that count — it is
GenCC's largest in-scope submitter (111 rows over 109 genes, measured
2026-08-06), and counting them would let one body vote twice. A ClinGen
`Disputed` or `Refuted` record vetoes the GenCC route entirely: a source trusted
to admit a gene alone must be trusted to refuse one.

**Read `atlas_curation` before presenting a row as curated content.** Every
`headline_confidence` is an upstream expert panel's call and never this atlas's,
and on 91 of the 92 rows there is no assessment by this atlas behind it at all.

**This population widened on 2026-08-06, from 23 genes to 92.** A consumer
written against the earlier data will have seen only `"definitive"` in
`headline_confidence` and only `"expert_curated"` in `validity_state`. Both
assumptions are now false — see the next section.

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
      "has_no_association_report": false,
      "has_source_discordance": false,
      "headline_confidence": "definitive",
      "no_association_reported_by": [],
      "lesion_groups": ["septal"],
      "symbol": "TBX5",
      "validity_state": "expert_curated",
      "variant_count": 0,
      "burden_row_count": 14,
      "independent_datasets": {
        "tested": 2,
        "enriched": 2,
        "corrected": 1,
        "families": [
          { "studies": ["PMID:34324492"], "state": "not_tested" },
          { "studies": ["PMID:40127276"], "state": "corrected" },
          { "studies": ["PMID:42230622"], "state": "nominal" }
        ]
      }
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
  above is TBX5's real row: its one in-scope ClinGen record — the **Syndromic
  Disorders** Gene Curation Expert Panel, for **Holt-Oram syndrome** — makes it
  `"expert_curated"` with a `"definitive"` headline. It is not a Congenital
  Heart Disease GCEP record, which this paragraph claimed until 2026-08-06 and
  which the mirror has never held for TBX5; `admitted_by.panel` and
  `admitted_by.disease_label` in the bundle now name both, so a reader can check
  this sentence rather than trust it. That a gene's headline can come from a
  panel whose remit is a syndrome, for a disease whose label names no cardiac
  feature, is the ordinary case here and not an exception: 10 of the 23
  originally published genes were graded by some panel other than the CHD GCEP. `headline_confidence` is `null` for a gene no authority has
  assessed. It is never `"no_known_association"` for that case: that
  classification is itself an assessed verdict ("a panel looked and found
  nothing"), and asserting it for a gene nobody has assessed would state a
  conclusion no authority reached.
- `atlas_curation` is `"curated"` when the atlas holds at least one curated
  assertion for the gene and `"not_yet_curated"` when it holds none — 91 of the
  92 rows published today. It is about *this* resource's own work and says
  nothing about the gene's validity: an uncurated row still carries whatever
  `headline_confidence` an expert panel assigned, which is exactly why the gene
  is listed. Read this rather than testing
  `assertion_count == 0`, so a browse filter does not have to reimplement the
  rule. It appears on the browse row and in the bundle, written in one place so
  the two cannot disagree.
- `validity_state` says how well curated the gene is: `"expert_curated"` (an
  in-scope ClinGen record exists), `"submitter_curated"` (only GenCC has
  assessed it) or `"uncurated"` (neither mirror has). `headline_confidence` is
  `null` **iff `validity_state` is not `"expert_curated"`** — the headline is
  the grade of the ClinGen record that admitted the gene, so a gene no ClinGen
  panel has graded in scope has no headline, however many GenCC submitters
  graded it and however highly. Measured 2026-08-06: the 16 `null` rows are
  exactly the 16 `"submitter_curated"` rows, as sets and not merely as counts.
  This rule read "`null` iff the gene is `"uncurated"` *or* every in-scope
  record maps to no rung" until 2026-08-06, which is false for every one of
  those 16 — GDF1 carries a G2P `Definitive` and a Labcorp `Strong`, and
  headlines `null` because neither is ClinGen's. Six genes in the
  committed mirrors (HGNC:24595, HGNC:4317, HGNC:6188, HGNC:7881, HGNC:9380,
  HGNC:9381) would resolve to `null` while `"submitter_curated"`: each carries
  exactly one in-scope record, an Orphanet `Supportive` submission, which
  `vocab.GENCC_CLASSIFICATIONS` maps to `None` because the submitter asserted an
  association without grading its evidence, not because nobody looked. Do not
  infer "no authority has assessed this gene" from `headline_confidence: null`
  alone — check `validity_state` for that.
- **`headline_confidence: null` and `"submitter_curated"` are now in the
  payload, and a consumer written before 2026-08-06 will not have seen them.**
  Measured against a real build (2026-08-06): `headline_confidence` is
  `"definitive"` on 23 rows, `"limited"` on 43, `"moderate"` on 9, `"strong"` on
  1 and `null` on 16; `validity_state` is `"expert_curated"` on 76 and
  `"submitter_curated"` on 16. Every `null` headline is on a
  `"submitter_curated"` row — no panel graded the gene, and `admitted_by` names
  the submitters instead. Until 2026-08-05 all 23 published rows carried
  `"definitive"` and `"expert_curated"` and those were the only values either
  field took, so any rendering tested against that data is now incomplete.
  `"uncurated"` still does not occur, and the six `Supportive`-only genes above
  are still in the mirrors rather than the output.

- **`limited` is not a weak yes.** It is a ClinGen expert panel saying the
  gene-disease case is *not yet made* — few probands, or evidence that is
  suggestive but not compelling. 43 of the 92 rows carry it, so a consumer that
  renders the string with no gloss, or that colours it on a good-to-bad ramp
  without saying what the rungs mean, presents the plurality of this atlas as
  better supported than the grading authority does. The site's own pages carry
  ClinGen's definitions beside every chip for this reason.
- `has_source_discordance` is `true` when one mirrored source contests the gene
  while the *other* supports it. It is narrower than `has_conflicting_evidence`,
  which also fires when a single source is internally split across diseases or
  panels — see [Contested genes](#contested-genes-the-one-consumer-obligation).
- `has_no_association_report` is **a third axis, not a third rung**, added in
  schema `2.10`. It is `true` when some authority reported *no known
  association* in scope while another asserted one, and
  `no_association_reported_by` names who — sorted, and empty exactly when the
  flag is `false`, so the two can never say different things.

  It is deliberately **not** folded into `has_conflicting_evidence`, which means
  exactly `disputed`/`refuted`. "A panel looked and found no reported evidence"
  is not "a panel disputes this", and merging them would give one laboratory's
  null result the weight of a chartered panel's refutation. ClinGen treats the
  assertion the same way — a distinct verdict rather than a rung of the
  definitive-to-limited ladder.

  Before `2.10` that disagreement reached no published byte at all. **GDF1 is
  the live case and currently the only one:** G2P grades it `Definitive`,
  Labcorp `Strong`, and Illumina reports `No Known Disease Relationship`. It
  published `has_conflicting_evidence: false` with nothing beside it, so a
  consumer was told the evidence did not conflict while two bodies disagreed
  about whether an association exists. A gene whose *only* in-scope record is
  `no_known_association` does **not** set this flag — it is not in disagreement
  with anything, and there are 9 such genes in the mirrors, none published.
- `confidence_by_lesion_group` applies the gene's mirrored `headline_confidence`
  to every lesion group its curated assertions declare — ClinGen and GenCC
  classify a gene against a disease, never against a specific lesion, so there
  is no finer-grained signal to divide the groups with. Its keys are exactly the
  gene's `lesion_groups`, so it is empty whenever that array is empty — which is
  every gene this atlas has not curated yet — and empty as well when
  `headline_confidence` is `null`.
- **An empty `confidence_by_lesion_group` does not mean "not assessed."** It
  means this atlas has recorded no lesion group for the gene, which is a fact
  about this resource's curation queue and says nothing about the gene's
  validity. Measured against a real build (2026-08-06), 75 of the 92 published
  rows carry a non-null `headline_confidence` together with
  `"confidence_by_lesion_group": {}` — an expert panel graded each of them, and
  no curator here has written an assertion for them yet. TBX5 is the only row
  with a non-empty map today (`{"septal": "definitive"}`). A consumer that
  renders an empty map as "not assessed" therefore understates 75 expert-panel
  calls as unassessed while still counting them among the genes it presents. Read `validity_state`
  for whether any authority has assessed the gene, and `atlas_curation` for
  whether this atlas has curated it; this field answers neither question.
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

The browse page: the same 92 rows `genes/index.json` publishes, rendered as a
table a person can read and filter. Each row carries the HGNC id — linked to
that gene's page — the symbol, `headline_confidence`, `validity_state`,
`atlas_curation` and the gene's lesion groups. Above the table sit a text box
matching id or symbol and five menus (lesion group, confidence, validity state,
atlas curation, burden evidence), whose options are the values actually present
in the build rather than every value the vocabulary allows. The confidence menu
now includes genes with no headline at all: 16 of the 92 are admitted on
submitter agreement and carry `null`, and the browse page renders those with a
`not classified` chip rather than an empty cell.

**Every row is rendered by the build, and the inline script only hides rows.**
There is no empty `<tbody>` filled in by a fetch, so `curl`, a crawler and a
reader with JavaScript disabled all get the complete table of 92 genes,
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
  "has_no_association_report": false,
  "no_association_reported_by": [],
  "lesion_groups": ["septal"],
  "validity": { "state": "expert_curated", "has_source_discordance": false, "records": [ … ] },
  "publications": ["PMID:8988165"],
  "assertions": [ { "id": "CHDA:AST:0000001", "lesion_groups": ["septal"], "evidence": [ … ] } ],
  "functional": [],
  "variants": [],
  "omics": {},
  "burden": [ { "study": "PMID:42230622", "cohort_stratum": "all", … } ]
  "independent_datasets": {
    "tested": 2, "enriched": 2, "corrected": 1,
    "families": [
      { "studies": ["PMID:34324492"], "state": "not_tested" },
      { "studies": ["PMID:40127276"], "state": "corrected" },
      { "studies": ["PMID:42230622"], "state": "nominal" }
    ]
  }
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
- `atlas_curation` reads the same here as on the browse row. On the 91 genes
  published today without curation here it is `"not_yet_curated"`, and
  `assertions`, `publications` and `functional` are then empty arrays: the page
  is the panel's classification plus whatever this atlas has recorded, which
  may be nothing at all. An empty `assertions` array is a curation gap, never a
  fetch that failed.
- `omics`, `variants` and `burden` are always present and may be empty. Read
  them without guarding for a missing key.
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
  - **A GenCC record has no date and no submission identifier, and this is a
    deliberate omission rather than a missing join.** `classification_date` is
    ClinGen's date, so it is `null` on every GenCC record — but GenCC does
    publish a date of its own, and an identifier for the submission itself, and
    this atlas mirrors both and publishes neither.
    `mirrors/gencc_submissions.tsv` carries `submitted_on` and `sgc_id`; no
    published byte carries either. Measured against a real build (2026-08-04):
    of the 344 validity records across the 92 gene bundles, 266 are GenCC, and
    all 266 have a non-null `submitted_on` and a non-null `sgc_id` in the
    mirror (re-measured 2026-08-06). So a consumer rendering a date column has
    nothing to show for 266 of 344 records, and cannot key a GenCC record by
    anything more stable than the `(source, disease, moi, submitter)` tuple that
    distinguishes it.
  - **`report_url` does not always lead back to the submission.** Of those 266
    GenCC records, 109 publish `report_url: null` and 48 more point at
    `https://panelapp-aus.org` — a homepage, not a submission — so 157 of the
    266 offer a reader no published route back to the record behind the
    classification (measured 2026-08-06). Attribute these to their `submitter`
    and present them as classifications this atlas mirrors; do not render them
    as though every one resolves to a citable report.
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
  (neither mirror has). An uncurated gene carries `records: []` and
  `has_source_discordance: false` — the explicit empty shape, not an absent
  `validity` key, so a consumer can tell "no authority has assessed this
  gene" from "the field is missing" without guessing. **Two of the three are
  reachable in the payload since 2026-08-06**: measured on a real build, 76
  bundles say `"expert_curated"` and 16 say `"submitter_curated"`. Before the
  gate widened only the first occurred, because a gene was published on an
  in-scope ClinGen `Definitive` call and that made it `"expert_curated"` by
  construction. `"uncurated"` remains unreachable — a published gene has at
  least one in-scope record by definition — so no bundle carries an empty
  `records` array, and a consumer that treats `records: []` as a failed fetch
  would still break on the first gene that had one.

`sources.json` (below) carries the licence terms this atlas mirrors ClinGen
and GenCC under, the same way it does for HPO.

### `admitted_by` and `asserted_by`: why this gene is here

Added in schema 2.8. The gate is a claim about external authorities, and these
two keys are what let a consumer **check** it rather than trust it.

`admitted_by` is the single warrant that cleared the gate. Exactly one warrant
admits a gene, and ClinGen is checked first, so a gene with both is published on
ClinGen's word:

```json
{ "authority": "clingen", "classification": "definitive",
  "disease": "MONDO:0007732", "disease_label": "Holt-Oram syndrome",
  "panel": "Syndromic Disorders Gene Curation Expert Panel", "submitters": [] }
```

```json
{ "authority": "gencc_agreement", "classification": null, "disease": null,
  "disease_label": null, "panel": null,
  "submitters": ["Ambry Genetics", "Labcorp Genetics (formerly Invitae)"] }
```

`authority` is `"clingen"` on 76 of the 92 genes and `"gencc_agreement"` on 16.
**Every key is present on every gene** — `submitters` is `[]` rather than absent
where ClinGen admitted the gene, and `classification`/`disease`/`disease_label`/
`panel` are `null` rather than absent where GenCC did — because an object whose
shape varies is a trap for a consumer reading a field off one gene and expecting
it on the next.

`admitted_by.classification` is the same value as `headline_confidence`, by
construction rather than by coincidence: both are read from the one function
that decides which record admitted the gene. A `"gencc_agreement"` gene has
`null` for both, however highly its submitters graded it.

`asserted_by` is every distinct institution **asserting** the gene in scope,
deduped by institution, with a count:

```json
{ "count": 5, "institutions": ["Ambry Genetics", "G2P",
  "Labcorp Genetics (formerly Invitae)", "Orphanet", "PanelApp Australia"] }
```

Two things it deliberately is not. It is **not a count of records** — ClinGen
submits to GenCC under its own name, so counting `gcep` and `submitter` values
naively overstates by exactly one per gene (134 against 111 over the 23 genes
published before the widening). And a **dissent is not an assertion**: a record
of `no_known_association`, `disputed` or `refuted` is excluded, because counting
the institutions that disagree as institutions that agree is the opposite of
what the field's name says. GDF1 is the case — it published `"count": 6`
including Illumina's `No Known Disease Relationship` until 2026-08-06, and GDF1
is the gene that motivates requiring two submitters in the first place. GenCC's
`Supportive` **is** counted: it asserts an association without grading its
evidence, which is an assertion, merely ungraded.

**`count` is not a score, and must never be rendered as one.** Eight authorities
asserting a gene is not evidence it is eight times better supported than a gene
with one; it frequently means it sits on more commercial test panels. The atlas
publishes no validity call of its own, and a rank derived from this count would
be exactly that.

Neither key appears in `genes/index.json` — the browse rows carry
`validity_state` instead, which separates the same two populations
(`"expert_curated"` 76, `"submitter_curated"` 16). Fetch the bundle for the
warrant itself.

### The bundle's `burden` array: per study, never pooled

Published rare-variant burden statistics for the gene, one object per
(study, cohort stratum, variant class, consequence class, frequency threshold).
**915 rows reach the API**, across 91 of the 92 published genes, from
3 studies (PMID:34324492, PMID:40127276, PMID:42230622), contributing
61, 150, 704 rows respectively. `mirrors/burden.tsv` holds 1,475 rows for 150 genes; the other 560,
covering 59 genes, are for genes the site does not publish and reach no
bundle and no page. They are held
so that widening the publication gate later needs no re-mirroring — which is
exactly what happened on 2026-08-06, when the gate widened from 23 genes to 92
and 625 already-mirrored rows became fetchable with no change to the mirror.
The same reason `mirrors/genes.tsv` registers 154 genes and 92 publish.

**One published gene carries no burden rows at all** and its `burden` array is
`[]`. That is not a fetch failure: the publication gate is an authority's
classification, not the presence of burden evidence, and a gene no burden study
covered is exactly what a widened gate admits.
`burden_row_count` on the browse row is this array's length.

Every count in this section is asserted against a real build by
`tests/test_docs_match_the_build.py`, so it fails rather than rots when a study
lands. It rotted once — this paragraph described two studies and 200 rows for
one commit after the third study shipped.

```json
{
  "study": "PMID:42230622",
  "cohort_stratum": "syndromic",
  "lesion_group": null,
  "variant_class": "snv_indel",
  "consequence_class": "lof",
  "origin": "any",
  "maf_max": 0.001,
  "count_unit": "individuals",
  "n_case_carriers": 5,
  "n_cases": 1471,
  "comparator": "control_cohort",
  "n_control_carriers": 0,
  "n_controls": 45082,
  "expected_count": null,
  "effect": null,
  "effect_measure": "odds_ratio",
  "effect_bound": "unbounded_above",
  "ci_low": 28.1,
  "ci_high": null,
  "pvalue": 3.13e-08,
  "pvalue_test": "fisher_exact",
  "pvalue_adjusted": null,
  "pvalue_adjustment": null,
  "case_cohorts": ["cnchd", "ddd", "nottingham"],
  "control_cohorts": ["ukbb"],
  "method_note": null,
  "source": "audain2026_sd3"
}
```

Every key is present on every object, `null` where the row's comparator does
not populate it, so a consumer never has to guard for a missing one.

**`count_unit` qualifies all four count keys, and they are not comparable
without it.** The four studies this schema was designed against do not count the
same thing:

| `count_unit` | numerator | denominator |
| --- | --- | --- |
| `individuals` | people carrying at least one qualifying variant | people sequenced |
| `alleles` | qualifying alleles observed | alleles called — roughly twice the people, and varying per gene with coverage |
| `de_novo_mutations` | de novo mutations observed | **trios**, not alleles and not people |

A consumer dividing `n_case_carriers` by `n_cases` across studies without
reading this key is doing arithmetic on three different quantities. Someone
carrying two qualifying variants counts once under `individuals` and twice under
`alleles`.

**`comparator` is the field the whole array turns on.** The published burden
literature answers one question — is this gene hit more often than expected? —
and differs only in what "expected" meant:

| `comparator` | expectation from | populated | `effect_measure` |
| --- | --- | --- | --- |
| `control_cohort` | `n_control_carriers` / `n_controls` | both control fields | `odds_ratio`, `rate_ratio` |
| `mutation_model` | `expected_count` | `expected_count` | `enrichment_ratio` |
| `none` | nothing (a case series) | neither | none, and no `pvalue` |

**Both `control_cohort` (840 rows) and `mutation_model` (75 rows) appear today;
`none` does not.** PMID:40127276 contributes the `mutation_model` rows — de novo
mutations in 3,887 trios against a mutability-based expectation — and they carry
null control fields, a non-null `expected_count`, and an `enrichment_ratio`
rather than an odds ratio. A consumer that reads only `n_control_carriers` will
see `null` on 75 of the 915 rows, several of which are the strongest results the
atlas publishes. The build refuses to publish a row whose statistic contradicts
its comparator.

**`consequence_class` is not a partition: `damaging` is the union of `lof` and
`missense_damaging`.** PMID:40127276 reports its primary analysis over
loss-of-function and damaging missense together, and publishes that composite
alongside its two components — 50 of the 915 rows. Its `n_case_carriers` is
exactly the sum of the two component rows for the same gene, stratum and
comparator, verified on every one. **Summing `n_case_carriers` across
`consequence_class` therefore double-counts those variants.** Group by
`consequence_class` and pick, or exclude `damaging`; do not aggregate over it.
The composite is published rather than dropped because its p-value is not a
function of its components' and is the statistic that study defines its results
by.

Five obligations, each of which is a wrong claim if you get it wrong:

1. **Never render `effect` without `effect_measure` beside it.** One column
   holds odds ratios and de novo enrichments alike, and an odds ratio of 3.1 and
   an enrichment of 3.1 are different claims. A cell reading `3.1` under a header
   reading "effect" equates them.
2. **`effect: null` with `effect_bound: "unbounded_above"` is the strongest
   result in the data, not a missing one.** Fisher's exact test returns an
   infinite odds ratio where no control carries; `Infinity` is accepted by
   `JSON.parse` nowhere, so the number cannot be published and `ci_low` carries
   the finding instead — "at least 28.1". 23 published rows today. Rendering it as
   blank or as "not tested" discards the clearest signals in the study.
3. **An absent (stratum, consequence) cell is not a null result — but why it is
   absent is the study's rule, not the atlas's.** Read it per study rather than
   as one law. For `PMID:42230622`, zero of its 1,192 mirror rows have no case
   carrier *and* no control carrier, so a 2×2 of all zeros supports no test and
   that study emitted none. The other two do not follow the same rule:
   `PMID:34324492` tests one consequence class by construction (CNV deletions),
   and `PMID:40127276` **observed 14,364 synonymous variants and still published
   no synonymous row**, because its gene-level table reports only damaging
   classes. So an absent cell means "this study did not report one here", and
   only for `PMID:42230622` does it additionally mean "no carrier on either
   side". One row in the mirror now does carry zeros on both sides.
4. **Do not filter out `consequence_class: "synonymous"`.** It is a study's own
   negative control — synonymous variants should show no enrichment — and where
   one is significant, that gene's comparison is poorly calibrated, which a
   reader can only see if you show it. **It is not available everywhere:** 261
   of the 915 published rows are synonymous and all of them come from
   `PMID:42230622`. The other two studies publish no synonymous row at all, so
   their results have no negative control on this page and must be read without
   one. (This obligation once called synonymous "the most numerous class, 435
   rows against 345 loss-of-function", measured against a one-study corpus; then
   "92 loss-of-function and 69 synonymous", measured against the 23-gene
   population. Both went stale the same way. Today loss-of-function is 306 and
   synonymous 261, beside 298 damaging missense and 50 damaging composite —
   re-derive these rather than citing them.)
5. **Read `pvalue_adjusted` where it is present.** A raw and a corrected
   p-value can point opposite ways: CHD7 in `PMID:34324492` is `0.0068` raw and
   `0.991` after the study's own family-wise permutation correction, so a
   consumer rendering `pvalue` alone shows as significant a result the study
   reported as null. `pvalue_adjustment` names the correction, because a
   family-wise permutation correction and a Bonferroni factor are different
   claims. Both keys are `null` where the study published none — and there the
   study's own comparison count, `tests_reported` in `publications.json`, is
   what a reader has instead. **The atlas computes no correction itself.**
6. **Do not pool across studies.** The CHD literature reuses cohorts, so a
   combined p-value counts the same people twice. `case_cohorts` and
   `control_cohorts` name the collections each row drew on precisely so overlap
   is visible.

   **These ids resolve against [`cohorts.json`](#cohortsjson), and you should
   resolve them.** `case_cohorts` and `control_cohorts` publish bare strings like
   `"taa_cases"`; that file maps each to a name, a URL where one exists, and the
   caveats that qualify every number drawn from it. `taa_cases` is 777 thoracic
   aortic aneurysm probands who **do not have congenital heart disease** — a row
   citing it is not a plain CHD case set, and only the description says so.

   Until schema `2.9` these ids resolved to nothing at all: the descriptions
   lived in the repository and reached a reader only through the gene *page*,
   under "About these cohorts". A consumer reading the JSON got the numbers
   without any of the qualifications, which is exactly the failure the cohort
   columns exist to prevent.

`n_cases` and `n_controls` are the row's own denominators, and they are what the
statistic beside them was computed from. They may differ from the figures a
paper's abstract reports — for `PMID:42230622` they are 3,876 cases (1,471
syndromic + 2,405 non-syndromic) against 45,082 controls, while the abstract
gives 4,747 and 52,881. The paper does not reconcile the two. Use the row's own
denominators: they are the ones its statistic was computed from.

### `independent_datasets`: a count of datasets, never a verdict

On every gene bundle and every `genes/index.json` row. It answers one question —
**how many independent datasets have looked at this gene, and what did they
find** — and it answers it by counting, never by combining. No pooled statistic
is computed here; see the previous section for why.

```json
"independent_datasets": {
  "tested": 2,
  "enriched": 2,
  "corrected": 1,
  "families": [
    { "studies": ["PMID:34324492"], "state": "not_tested" },
    { "studies": ["PMID:40127276"], "state": "corrected" },
    { "studies": ["PMID:42230622"], "state": "nominal" }
  ]
}
```

**A family is not a study.** It is a set of studies that share at least one
sample collection, walked transitively — so two papers drawing on the same
cohort describe the same people and count **once**. Today the three curated
studies draw on disjoint collections, so there are three singleton families; the
grouping exists so that the day a fourth study reuses DDD or PCGC, nobody is
told a reused cohort is independent evidence.

| `state` | meaning |
| --- | --- |
| `corrected` | enriched, and survives that study's own published correction |
| `nominal` | enriched, but nominal only — or the study published no correction at all |
| `no_enrichment` | that dataset tested the gene and detected nothing |
| `not_tested` | that dataset did not test this gene |

A family counts as enriched only if a non-synonymous row has `pvalue < 0.05`
**and points upward**. Direction is read from what the study published — the
`unbounded_above` flag, then `effect` against 1, then the published rates where
a study reports no effect at all. A significantly *depleted* row is not
agreement, and a synonymous row is never support: it is the negative control.

**Three obligations, and the first one matters most.**

1. **This is not a validity call, and it must never be rendered as one.**
   `headline_confidence` beside it is a mirrored ClinGen classification. A
   consumer that renders "0 of 2" as a verdict next to a green `definitive` chip
   tells a clinician the data contradict the classification. They do not:
   **KDM6A causes Kabuki syndrome and shows nothing in either dataset that
   tested it**, because burden tests at these cohort sizes routinely detect
   nothing for genes with overwhelming family and functional evidence. If you
   display this, display alongside it that no enrichment here is not evidence
   against a gene.
2. **Read `tested` as the denominator, never `len(families)`.** A gene absent
   from a study's panel was not examined and found wanting; it was not examined.
3. **`not_tested` must render distinguishably from `no_enrichment`.** Collapsing
   them is what turns "nobody looked" into "somebody looked and found nothing".

`families` is always present, ordered deterministically, and **always carries one
entry per cohort family in the whole corpus** — three today. A gene no study
reported gets three `not_tested` entries, not an empty array; `[]` occurs only if
the corpus has no burden data at all. That is deliberate: every gene shows the
same slots, so a dataset that did not test this gene is visible as an absence
rather than as a shorter list.

**The browse page heads this column "burden across studies", not
`independent_datasets`.** The key names what is counted for a program; the
header names what it is for a reader. Both avoid a verdict word on purpose —
"replicated in" was the shorter candidate and was rejected, because for a gene
showing 0 of 2 it reads as "not replicated", which is a claim the data do not
make.

### The bundle's `expression_profile` object: a developmental transcriptome, never a contrast

On every gene bundle. Always present, and empty (`{"datasets": []}`) for a
gene no profile dataset's mirror covers — the same rule `omics` and `burden`
above already keep, so "no profile dataset says anything about this gene"
cannot be confused with "the build dropped it".

```json
{ "expression_profile": { "datasets": [] } }
```

**That is still the real shape for a gene no profile dataset's rows mention.**
`mirrors/profiles/E-MTAB-6814.tsv` is committed, and `manifest.json`'s
`counts.profile_genes`/`counts.profile_datasets` say how far that reaches:
92 genes carry a developmental expression profile, across 1 dataset. Most
published genes carry the populated shape today, not the empty one above.
The nested shape below is illustrative — restructured for
readability (a real gene's `stages` runs to a dozen or more entries, one per
token the dataset declares) rather than copied verbatim from one bundle — but
every value in its `phase` block is real: it is `assign_phase`'s actual,
reproducible answer for a 7-elapsed-week stage against the boundaries
`curation/cardiac_phases.yaml` curates today. The `stage` token, the `unit`
and `n_genes` are real too: stage tokens are the dataset's own long-form
strings, this dataset reports `tpm`, and its grids were built from a 58,735-gene
transcriptome. The `tau`, the abundances and the percentiles are the invented
part.

```json
{
  "expression_profile": {
    "datasets": [
      {
        "dataset": "E-MTAB-6814",
        "quantile_shard": "omics/profile_quantiles/E-MTAB-6814.json",
        "stages": [
          {
            "stage": "7 week post conception",
            "phase": { "outcome": "matched",
                       "phase_ids": ["ventricular_septum_morphogenesis",
                                     "heart_valve_morphogenesis"],
                       "reason": null },
            "specificity": { "tau": 0.62, "scale": "log2(x+1)", "method": "…",
                              "tissues": ["heart", "kidney", "liver"],
                              "n_tissues": 3, "highest_in": "heart",
                              "medians": { "heart": 42.0, "kidney": 8.0, "liver": 6.0 } },
            "specificity_unavailable_reason": null,
            "tissues": [
              { "tissue": "heart", "median_abundance": 42.0, "unit": "tpm",
                "n_samples": 3,
                "placement": { "q25_percentile": 44, "median_percentile": 50,
                                "q75_percentile": 56, "median_abundance": 42.0,
                                "unit": "tpm", "n_samples": 3, "n_genes": 58735,
                                "method": "lowest percentile of a tied breakpoint (bisect_left)" },
                "not_placed_reason": null }
            ]
          }
        ]
      }
    ]
  }
}
```

**Two mirror tables feed this, and each reaches a consumer a different way.**
`mirrors/profiles/<accession>.tsv` (`dataset`, `gene`, `tissue`, `stage`,
`median_abundance`, `unit`, `q25`, `q75`, `n_samples`) is one of the four
tables [`omics/<modality>/<accession>.json`](#omicsmodalityaccessionjson)
above shards like any other — `profiles` is one of that section's four
modalities — so its raw rows reach `omics/profiles/<accession>.json`
regardless of gene publication, and `median_abundance`/`unit`/`n_samples`
below are read straight from it. `mirrors/profile_quantiles/<accession>.tsv`
never goes through that mechanism at all: it has no gene column, so
`build_omics` skips it outright, and it is published only as
[`omics/profile_quantiles/<accession>.json`](#omicsprofile_quantilesaccessionjson)
(below), one full shard per dataset.

- `datasets` is one entry per profile dataset whose mirror mentions this
  gene. `quantile_shard` names that accession's
  `omics/profile_quantiles/<accession>.json` payload, or `null` if the build
  wrote no shard for it.
- `stages` is one entry per stage token the mirror's rows use for this
  (gene, dataset) pair — `stage` is `null` for a measurement with no
  developmental stage recorded at all.

  **The array is in chronological order, and was not before schema `2.12`.**
  It is sorted on the dataset record's own curated
  [`order`](#datasetsjson) — not on `wpc`, which is `null` for every
  post-natal stage, and not on the stage token. A token the dataset does not
  declare sorts after every one it does; the `null` stage sorts last. The
  entries themselves do not carry `order`; the array's sequence is what
  publishes it, so read the series as it arrives.

  Through schema `2.11` this array came back **alphabetically**: `4 week post
  conception` published eighth, after `19 week post conception`, and `elderly`
  second of the eight post-natal stages. That was deterministic and
  reproducible between builds — and it was not a chronology. If you plotted
  these entries in array order against an earlier release you drew a
  developmental trajectory with its axis shuffled, and the same code now draws
  a different, correct picture. If you sorted the stage tokens yourself to
  work around it, stop: sorting the tokens is what the defect did.
- `phase` places the stage in the curated cardiac morphogenetic window
  (`curation/cardiac_phases.yaml`):

  | `phase.outcome` | meaning |
  | --- | --- |
  | `matched` | the stage's own developmental age (`wpc`) falls inside one or more curated cardiac phases; `phase_ids` names all of them |
  | `outside_window` | a real age exists, but no curated phase covers it |
  | `post_natal` | the stage has no developmental age at all |
  | `undeclared` | the dataset's own record does not declare this stage token |
  | `null` | the row itself carries no stage token; `reason` says so |

  **`phase_ids` may name more than one phase, and often does.** Human cardiac
  morphogenesis runs several processes concurrently — at 6 elapsed weeks post
  conception, atrial septation, ventricular septation and outflow tract
  septation are all underway at once — so a stage legitimately matches every
  phase whose window contains its `wpc`, never just the nearest or the first
  declared. `curation/cardiac_phases.yaml` curates six phases transcribed from
  Buijtendijk et al. 2020 (PMID:32048790), each carrying the Gene Ontology
  term and the Carnegie stage(s)/HsapDv id(s) its boundaries were read from —
  not published in this bundle field, only in the curation file itself, so a
  consumer auditing a boundary starts there. `phase_ids` is empty exactly
  when `outcome` is not `matched`.

  **One of the six, `heart_looping`, never appears in `phase_ids` at all.**
  The source states when it starts but never states when it ends, and this
  atlas will not guess a cutoff nothing supports: `end_basis` on that phase
  is `not_stated` rather than `stated`, and a phase in that state is excluded
  from every stage's match, at every `wpc`, not merely "too far" past its own
  start. It is still curated — with its real, sourced start — so a consumer
  reading `curation/cardiac_phases.yaml` directly sees it; a consumer reading
  only bundle JSON never does, because this atlas would rather publish
  nothing for that phase than assert it is still running at a stage the
  source gives no basis for.
- `specificity` is Yanai's τ (tau) over every organ a dataset sampled at one
  stage, or `null` with `specificity_unavailable_reason` naming why:

  | `specificity_unavailable_reason` | meaning |
  | --- | --- |
  | `dataset_not_registered` | the mirror names a dataset accession with no curated record at all |
  | `detection_floor_undeclared` | the dataset record exists but declares no detection floor |
  | `one_organ_sampled` | fewer than two organs were sampled at this stage — τ's denominator is `n − 1` |
  | `peak_below_detection_floor` | every sampled organ's median is below the dataset's own floor |
  | `undefined` | a residual this atlas cannot itself further diagnose |

  `tissue`-level entries carry a **related but not identical** vocabulary
  through `not_placed_reason` — the two share `dataset_not_registered` and
  `detection_floor_undeclared`, but a single tissue can never be short of
  organs or undefined the way a whole stage's τ can, and gains two reasons of
  its own instead: `no_quantile_grid` (no complete 101-point breakpoint grid
  exists for this cell) and `below_detection_floor` (a grid and a floor both
  exist; the median itself falls below it). Four possible values, never the
  stage-level five.
- **τ is published bare — no adjective, no band.** A value of 0.71 is not
  glossed "intermediate" or "specific": choosing a threshold would be a
  classification this atlas authors about how tissue-specific a gene is,
  which its governing rule forbids for every derived figure in this layer, the
  same way D12 forbids the atlas authoring a gene-disease validity call. If
  you display τ, display the number.
- **τ is blind to *where* a gene peaks — read `highest_in`, never the number
  alone.** τ measures concentration, not location. Measured: a gene at heart
  20, liver 200 and five other organs at 5 (rpkm; seven organs sampled) scores
  τ = 0.963 computed on the raw values and 0.623 on the `log2(x+1)` scale this
  atlas actually publishes — and in both, the peak organ is **liver**, not
  heart. A consumer who reads a moderately high τ as "heart-preferential" for
  a gene concentrated in liver states the opposite of the truth. `highest_in`
  is `null` when two or more organs tie exactly at the peak, where naming one
  of them would be arbitrary.
- `scale` (`"log2(x+1)"`) and `method` travel beside every τ so a consumer can
  re-derive it without reading this atlas's source rather than trusting a
  paraphrase of it. `method` reads exactly: tau (Yanai et al. 2005): mean over
  organs of (1 - x_i/x_max), x = log2(median+1); a negative median is clamped
  to 0 before the transform, and every organ's raw median is used even below
  the dataset's detection floor. `medians` is τ's own input, published exactly
  as measured — including any organ a gene bundle's own `omics` preview may
  have already dropped from `top` — so τ's inputs stay reachable from the one
  payload that publishes τ (design decision D39(b)).
- `placement` is one gene's percentile band in one (dataset, tissue, stage)
  cell: `q25_percentile`/`median_percentile`/`q75_percentile`, the outer two
  `null` below `n_samples = 3`, where a quartile of two points is not a
  quartile. `method`
  (`"lowest percentile of a tied breakpoint (bisect_left)"`) names the
  tie-break — a value tied with a run of identical breakpoints reports the
  *lowest* percentile in that run, so an unexpressed gene never reads above
  the median it shares with every other unexpressed one. `n_genes` is the
  size of the *source* transcriptome that cell's grid was built from — tens
  of thousands of genes, never the 92 this atlas publishes — because "top 4%"
  names no population without its denominator.

**Percentiles are not comparable across organs, and — more importantly — not
comparable across developmental stages.** Each organ transcribes a different
fraction of the whole gene universe, so "top 10% in liver" and "top 10% in
heart" describe two different reference distributions, not one ranking read
two ways. And a reference distribution's own shape changes as the heart
matures, so a gene with flat, unchanging absolute abundance across
development can still show a *moving* percentile from one stage to the next,
purely because everything else in the transcriptome is moving around it. The
cross-stage caveat matters more here, because this layer's headline
question — is this gene's cardiac expression distinctive, and does that
change over development — is asked *within* one organ, across its own
stages: reading a percentile trend without this caveat can report
development where none occurred, or miss it where it did.

**Why the quantile grid is published at all, and exactly how far that goes.**
A percentile's input is the whole source transcriptome — tens of thousands of
genes — and this atlas publishes 92 of them; design decision D32 forbids
re-hosting the matrix a percentile was read against. Without a published grid,
`median_percentile` would be a number nobody outside this atlas could check.
[`omics/profile_quantiles/<accession>.json`](#omicsprofile_quantilesaccessionjson)
(below) publishes exactly the 101 breakpoints instead — which is a
**one-level guarantee, not an unbounded one**: a gene's percentile is
re-derivable from its own `median_abundance` and the published breakpoints;
the breakpoints themselves are not re-derivable, because they were read off
the whole matrix D32 keeps unpublished. State it as the trade it is —
auditable one level down, not provably correct all the way back to the source
data.

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
or, for the 91 genes published today with no curation here, a paragraph saying
exactly that:

> The atlas has **not yet curated** a lesion assertion for this gene. The
> classification above is an expert panel's, mirrored with its provenance
> intact; no classification on this page is the atlas's own assessment.

That wording is for a gene a ClinGen panel graded — 76 of the 92. The 16
admitted on submitter agreement have no panel classification to describe, and
say so instead:

> The atlas has **not yet curated** a lesion assertion for this gene, and **no
> ClinGen expert panel has graded it** for a disease an external authority
> treats as congenital heart disease. It is published because the Gene Curation
> Coalition submitters named above independently assert it; their
> classifications are mirrored with their provenance intact, and no
> classification on this page is the atlas's own assessment.

`validity_state` selects between them, not `headline_confidence is null`: the
two agree on every gene published today, and two figures that are equal are one
figure to every test.

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
- **`top` is capped at 25 rows.** For `expression`, `proteomics` and `phospho`
  it is a preview, ranked by significance (ascending FDR, ties broken by the
  table's own sort order), not a page of results. `count` is frequently
  larger, and the cap is not carried in the payload, so do not infer
  completeness from `len(top)`.
- **`profiles` ranks `top` differently, because that table has no
  significance column at all** — it is an abundance table, not a contrast.
  The slice is *stratified*: the cardiac series leads, ranked by the derived
  `placement.median_percentile`
  (see [`expression_profile`](#the-bundles-expression_profile-object-a-developmental-transcriptome-never-a-contrast)
  above), but slots are **reserved for the non-cardiac tissues**, up to half of `top`.

  With `top` capped at 25 that reservation is at most 12, so every other
  tissue present is represented only while there are 12 or fewer of them —
  measured, and true of this dataset, which has six. Past that the cardiac
  series can consume its budget and the remaining tissues share what is left,
  so do not read the reservation as a guarantee that every tissue appears.

  **The reservation exists because of τ.** A gene's `tau` figure (also
  documented there) is computed over every organ a dataset sampled at one
  developmental stage, and design decision D39(b) requires τ's own inputs to
  stay reachable from the same payload that publishes τ. Without the
  reservation, a bulk developmental atlas sampling many stages per organ
  fills every one of the 25 rows with the cardiac series alone — measured:
  ranking `profiles` the way the other three tables are ranked put 0 of 14
  heart rows in a 25-row slice at 14 stages per organ (every row ties at "no
  FDR", so the tie-break falls back to an alphabetical tissue order that
  heart loses), and ranking cardiac-first with no reservation at all
  over-corrected, leaving 0 of 6 comparison organs in a 25-row slice at 23
  stages per organ. **A consumer who does not know this will read the
  reserved rows as a ranking bug**: they are exactly the organs a gene's own
  `tau` was computed over, shown so that figure is auditable from the same
  payload that publishes it.

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

`mirrors/genes.tsv` is the mirror table this site publishes least of, and the
one most likely to be assumed present. It holds 154 rows — exactly the genes
ClinGen or GenCC curates within the scope `curation/chd_scope.yaml` declares, of
which 92 clear the publication gate — across nine columns, and **five of those
nine reach no published byte at all**: `ensembl_gene`, `ncbi_gene`, `locus`,
`uniprot` and `mane_select`, each populated on all 154 rows (measured
2026-08-04). Only `hgnc_id`, `symbol`, `name` and `aliases` are published, the
last two as `terms` in `search/index.json.gz`.

So **this API carries no cross-reference to any other identifier space.** A
consumer that needs an Ensembl gene id, a MANE Select transcript, a UniProt
accession or a cytogenetic band for a published gene must resolve `hgnc_id`
against HGNC itself; there is no key to look for and no bundle field that will
appear later without a `schema_version` bump. The atlas identifies a gene by its
HGNC id and leaves the mapping to the authority that maintains it.

## `omics/profile_quantiles/<accession>.json`

The percentile grid a `placement`
(see [`expression_profile`](#the-bundles-expression_profile-object-a-developmental-transcriptome-never-a-contrast)
above) was read against — the other half of the trade design decision D39(b)
makes. This is **not** one of the four modalities in the section above: it has
no gene column at all — there is no gene to attribute one breakpoint row to —
so `build_omics` skips it outright (that function shards only the tables its
own `_GENE_COLUMN` map names). A separate builder
(`profiles.build_profile_quantiles`) emits it instead, to this separate,
sibling path under `omics/`.

```json
{
  "table": "profile_quantiles",
  "rows": [
    { "dataset": "E-MTAB-6814", "tissue": "heart",
      "stage": "7 week post conception", "percentile": 50,
      "value": 0.0, "unit": "tpm", "n_genes": 58735 }
  ]
}
```

Exactly the `mirrors/profile_quantiles/<accession>.tsv` rows for that dataset,
republished verbatim — `dataset`, `tissue`, `stage`, `percentile`, `value`,
`unit`, `n_genes` — sorted by `(dataset, tissue, stage, percentile)`. A
complete grid is 101 rows for one (dataset, tissue, stage) cell, `percentile`
running 0 through 100; `n_genes` is the size of the source transcriptome that
grid was built from, the same figure a `placement` in that cell publishes
under its own `n_genes`.

A gene bundle reaches this file through
`expression_profile.datasets[].quantile_shard`, never by constructing the path
itself — the "never construct a path" rule this document opens with. Without
this file, `median_percentile` would be an unauditable number: the whole point
of publishing it is so a consumer can look up the median a gene reported and
confirm which percentile it lands in, the same lookup this atlas's own build
performed once.

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
      "ancestry": [],
      "tests_reported": null
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
{
  "datasets": [
    {
      "id": "E-MTAB-6814", "archive": "arrayexpress", "design": "profile",
      "technology": "bulk_rnaseq", "organism": "NCBITaxon:9606",
      "tissue": "forebrain, heart, hindbrain, kidney, liver, ovary, testis",
      "developmental_stage": "4 weeks post conception to elderly adult (13 prenatal, 8 post-natal stages)",
      "n_samples": 287, "licence": "EMBL-EBI Terms of Use",
      "publication": "PMID:31243369", "contrasts": [],
      "cardiac_tissues": ["heart"],
      "detection_floor": 1.0, "floor_source": "…", "quantile_estimator": "linear",
      "stages": [
        { "token": "4 week post conception", "wpc": 4.0, "order": 1 },
        { "token": "19 week post conception", "wpc": 19.0, "order": 13 },
        { "token": "neonate", "wpc": null, "order": 14 },
        { "token": "elderly", "wpc": null, "order": 21 }
      ]
    }
  ]
}
```

`stages` and `floor_source` are abridged above — the real record declares all
21 stage tokens and spells the floor's provenance out in full. Every other
value is copied from a build of the committed corpus.

One record per omics dataset, serialised generically from the curated model
(`build/literature.py`'s `_dump`, the same function `publications.json`,
`featured.json` and `phenotypes.json` use) — every field it declares reaches
this file, with nothing filtered out: `id`, `archive`, `technology`, `tissue`,
`developmental_stage`, `organism`, `n_samples`, `licence`, `contrasts`,
`design`, `cardiac_tissues`, `stages`, `detection_floor`, `floor_source`,
`quantile_estimator` and `publication`. This is what an omics row's `dataset`
column resolves against, the way `publications.json` resolves a PMID.

**The last six (`design` onward) are what distinguishes a bulk developmental
expression dataset from every other kind mirrored here.** `design` is either
`"contrast"` (a differential-expression comparison — `expression`,
`proteomics`, `phospho`) or `"profile"` (a raw abundance series across
development — `profiles`, see
[`expression_profile`](#the-bundles-expression_profile-object-a-developmental-transcriptome-never-a-contrast)
above). `cardiac_tissues` and `stages` are `[]` and
`detection_floor`/`floor_source`/`quantile_estimator` are `null` on a
`"contrast"` record — never omitted — so every record has one shape whichever
design it declares. `cardiac_tissues` names which of a profile dataset's own
tissue tokens this atlas reads as "the heart" for that dataset.

`stages` is one entry per developmental stage the dataset declares. Each
carries the stage's own `token` — the literal string the mirror's rows use,
which is what an `expression_profile` stage joins against — and two positions
on the developmental axis:

- `wpc` is the stage's age in weeks post conception, and is `null` for every
  post-natal stage — 8 of E-MTAB-6814's 21. That null is not a gap in the
  curation; a stage after birth has no such age.
- `order` (added in schema `2.12`) is the **curated chronological position**:
  1-based, ascending, and unique within a dataset. It covers the post-natal
  stages, which is the whole reason it exists — `wpc` is `null` for all of
  them, so those tokens carry no other chronology, and the only remaining way
  to sequence them is the stage *string*, which puts `elderly` before
  `infant`. It is what orders the developmental series published on every gene
  bundle: see
  [`expression_profile`](#the-bundles-expression_profile-object-a-developmental-transcriptome-never-a-contrast)
  above.

**Sort on `order` rather than trusting this array's own sequence.** The
records in `datasets.json` are serialised straight from the curated file, so
`stages` here arrives in the order a curator declared it. That matches `order`
today and nothing enforces it to — the atlas checks that `order` is unique
(`PRF011`) and does not contradict a stage's own `wpc` (`PRF012`), not that
the block was typed in sequence. The gene bundle's stage array *is* sorted by
`order`; this one is a curated record printed as written.

Which is not academic: **this array's own sequence also moved in `2.12`.** The
same commit that numbered the stages re-sequenced the eight post-natal ones in
the curated file, so they now arrive `neonate` first rather than `adolescent`
first. Every token and every `wpc` is unchanged — only the order is different,
and only for the post-natal block.

The committed corpus holds one dataset today, the `"profile"`-design
`E-MTAB-6814` shown above. No `"contrast"`-design dataset has been curated
yet, so every `contrasts` array a consumer meets here is empty.

## `cohorts.json`

The sample collections every burden row names, and the caveats that qualify
them. Added in schema `2.9`.

```json
{
  "cohorts": [
    {
      "id": "taa_cases",
      "name": "Sporadic thoracic aortic aneurysm case series (PMID:34324492)",
      "description": "777 individuals with sporadic thoracic aortic aneurysm, included in the CNV case set of PMID:34324492 alongside the CHD cases. TAA IS NOT CONGENITAL HEART DISEASE and is explicitly out of this atlas's scope …",
      "url": null
    }
  ]
}
```

**This file is the reason `case_cohorts` and `control_cohorts` are worth
reading.** A burden row names its collections by bare id; without this table
those ids resolve to nothing, and a consumer computing anything from the row is
working with numbers stripped of the sentences that say what they count.

- `id` is exactly the string a burden row's `case_cohorts` / `control_cohorts`
  array carries. Every id in every published row appears here — the build
  publishes the whole curated registry rather than the subset the published rows
  happen to cite, and a validator (`BUR009`) refuses a row naming a collection
  this file does not carry, so resolution cannot fail.
- `name` is the collection's full name, always a non-empty string.
- `description` carries the caveats, and it is the field this file exists for.
  Examples from the current registry: `ukbb`'s controls are adults recruited at
  40–69 while the cases were largely enrolled in childhood; `ddd` ascertains on
  developmental disorder rather than heart disease, so its contribution is
  enriched for **syndromic** CHD; `gnomad_controls` is not screened for
  congenital heart disease. None of these is expressible as a column.
- `url` is the collection's public page, or **`null`** where it has none (4 of
  the 13 today). Always present, never omitted.

The array is sorted by `id`. There is no `role` field: a collection is cases in
one study and could be controls in another, so which it was is a property of the
burden row that cites it — the column it appears in — not of the collection.

`counts.cohorts` in `manifest.json` is the length of this array. It is **not**
`counts.cohort_families`, which counts how many independent collections the
burden evidence groups into (3 today) and is what
[`independent_datasets`](#independent_datasets-a-count-of-datasets-never-a-verdict)
is derived from.

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
- Genes come from the **published** set — the same 92 genes `genes/index.json`
  lists and `genes/<slug>.json` serves, each one admitted by an external
  authority under the rule above. Deliberately not the assertion set: keyed on
  that, the index would hold one record while the site published 92 bundles, so
  a visitor typing GATA4 would get nothing while `genes/HGNC_4173.json` was
  being served. Deliberately not the 154-gene mirror registry either, which
  would offer 62 results whose `path` names a bundle no builder wrote. A gene outside the
  published set is not searchable because it has no page to open.
- **A search record carries no curation or confidence signal.** `kind`, `id`,
  `label`, `path` and `terms` are the whole record — there is no
  `headline_confidence`, no `validity_state`, no `atlas_curation` and no
  `has_conflicting_evidence` in it. The obligations this document places on a
  consumer — read `atlas_curation` before presenting a row as curated content,
  and pair `headline_confidence` with `has_conflicting_evidence` — therefore
  **cannot be met from this file alone**. A typeahead built on it must join each
  hit back to `genes/index.json`, or fetch the bundle at `path`, before it
  labels a result as anything beyond a name that matched. This is not a
  hypothetical gap: of the 92 genes indexed today, 91 are `"not_yet_curated"`
  (measured 2026-08-06), and 16 carry no expert-panel grade at all, so a UI that
  presents a search hit as atlas-curated content would be wrong about nearly
  every one of them.

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
to skip implementing them.** Measured against a real build (2026-08-06): of the
154 mirrored genes, exactly one — LEFTY2, HGNC:3122, where ClinGen's own
`Disputed` call sits alongside GenCC's supportive one — sets
`has_conflicting_evidence`, and it is the same gene that sets
`has_source_discordance`. It is not among the 92 published, and **the reason it
is not changed on 2026-08-06.** Under the old gate it simply lacked a ClinGen
`Definitive` call. Under the widened one it has two GenCC submitters at
`Limited` (G2P and PanelApp Australia) and would qualify on submitter agreement
— so the gate now refuses it explicitly, because ClinGen's Congenital Heart
Disease GCEP disputes it and a source trusted to admit a gene alone must be
trusted to refuse one. Without that veto LEFTY2 would be published, and a
headline taken over every mirrored record would have rendered it `limited`,
burying the dispute in a payload no consumer had reason to inspect. So a
consumer testing their rendering against the live data will still find nothing
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
consumer that renders it alone presents 91 of the 92 published genes as this
atlas's assessment of them. It is not: it is a ClinGen expert panel's, and this
atlas has recorded nothing about those genes beyond republishing it with its
provenance. On 16 of the 92 it is `null` — no panel graded the gene, and
`admitted_by` names the GenCC submitters whose agreement admitted it instead.

**A consumer must not present a `not_yet_curated` gene as one the atlas has
assessed.** `atlas_curation` appears on the browse row and in the bundle,
written in one place so the two cannot disagree, and it is what tells the two
kinds of row apart — read it rather than testing `assertion_count == 0`, so a
filter does not have to reimplement the rule. An empty `assertions` array is a
curation gap, never a verdict and never a fetch that failed, the same way
`headline_confidence: null` is not `no_known_association`. This atlas's own
gene pages state the gap in prose for exactly that reason.
