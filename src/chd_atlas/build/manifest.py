# src/chd_atlas/build/manifest.py
"""The manifest: what was built, from which commit, and a checksum per file.

This is what a consumer verifies a download against and what a cache keys on,
so `files` maps each published URL to the sha256 of the exact bytes served at
it — the compressed bytes for a `.gz`, not the JSON inside. A checksum that
disagrees with the file is worse than no manifest at all: a consumer verifying
a download rejects a file that is correct.

Nothing here varies between two builds of the same commit. A build timestamp
was considered and rejected: it would make the output irreproducible, so no
test could assert byte-identical rebuilds and the checksums would stop being a
pure function of the inputs. A consumer wanting a publication date reads the
commit's.

The manifest cannot list itself — its own checksum would have to be computed
from bytes containing that checksum. The exclusion is not a filter but an
ordering: `write_json` encodes the payload before it records the new entry, so
`manifest.json` enters `checksums` a moment after the bytes naming every other
file are already fixed.
"""

from __future__ import annotations

import os
import subprocess
from collections.abc import Mapping
from pathlib import Path
from typing import Final

from chd_atlas.build.emit import Emitter
from chd_atlas.corpus import Corpus

# What a consumer can rely on across two builds, in the usual major.minor sense
# — spelled out here because the field was published for three releases without
# anyone saying what it meant.
#
# MINOR rises when a field is added and nothing is removed or repurposed: a
# consumer written against the previous version keeps working, and one that
# wants the new field can test for the version rather than for the key.
# MAJOR rises when a field changes shape or leaves, which breaks such a reader.
#
# 1.1 added `genes` to every omics shard row (issue #3) and
# `conflicting_lesion_groups` to every gene index row (issue #4). Both are
# additive, so 1.0 readers are unaffected.
#
# 2.0 removed `classification` and `source_tier` from the curated assertion and
# moved gene-disease validity to mirrored, attributed records (design decision
# D12: the atlas no longer authors a validity call of its own). A 1.x reader
# looking for a classification on the assertion finds none, which is a
# breaking change and so MAJOR rather than MINOR, even though the release's
# other change — the gene bundle's new `validity` object, carrying every
# mirrored record and its provenance — is purely additive on its own.
#
# 2.1 added `status`, so a programmatic consumer can read the atlas's own
# readiness without scraping `index.html`'s prose for it. Nothing existing
# changed shape or left, so this is MINOR despite following a MAJOR release —
# each release's letter is decided against what it changes, never against the
# rung before it.
#
# 2.2 added `atlas_curation` to every gene index row and every gene bundle
# (design decision D21/D22). The same release widened the published population
# from the genes carrying a curated assertion to the genes a ClinGen expert
# panel calls definitive — 1 row to 23 — which is not a schema change at all:
# more rows of an unchanged shape. The new field is what stops that widening
# being a silent one, since 22 of the 23 rows now describe a gene this atlas has
# not curated. Additive, so MINOR: a 2.1 reader keeps working.
#
# That widening is also where the rule above had to be extended, because the
# rule speaks only about fields and a release can change which *rows* appear
# without touching one. A population change inside an unchanged shape stays
# MINOR: a consumer's parsing is unaffected, since every row still carries every
# field 2.1 published and none of them has changed shape. A consumer's
# *display* is very much affected — 22 of the 23 rows describe a gene a 2.1
# reader had never seen and would render as this atlas's own curated content —
# so the warning belongs in `docs/data-api.md`, where a reader meets the file,
# rather than in a version letter that cannot carry it. Read the letter for
# "will my parser still work"; read the document for "does my page still say
# something true".
#
# 2.2 also added `genes/index.html` and `genes/<slug>.html`, taking the site
# from the single page 2.1 shipped — `index.html`, this build's first HTML
# artifact — to 25 (measured 2026-08-04 against a build of the committed
# corpus). Additive in a stronger sense than a field is: no consumer of the data
# API is affected by a page appearing beside a payload, and every page is
# checksummed into `files` like any other published byte.
#
# 2.3 adds `burden` to every gene bundle and `burden_row_count` to every
# `genes/index.json` row: the per-study rare-variant burden statistics from
# `mirrors/burden.tsv` -- 1,192 rows over 145 genes *at that version*, from one
# study. The mirror has grown twice since (1,475 rows, 150 genes, three studies);
# these figures describe what 2.3 shipped and are deliberately not updated,
# because a changelog entry that tracks the present tells a reader nothing about
# the version it names. Current counts live in `docs/data-api.md`, where
# `tests/test_docs_match_the_build.py` measures them against a real build.
#
# Additive, so MINOR — a 2.2 reader keeps working, and both keys are always
# present (an empty array and a zero for a gene no study reported), so a
# consumer never has to guard for a missing one.
#
# The display warning that accompanied 2.2's population change applies here in a
# sharper form, and `docs/data-api.md` carries it: `effect` must never be
# rendered without `effect_measure` beside it, because an odds ratio of 3.1 and
# a de novo enrichment of 3.1 are different claims sharing one column; and
# `effect` is `null` on 34 rows whose odds ratio the study published as
# infinite, where `effect_bound` is `unbounded_above` and `ci_low` carries the
# whole finding. A consumer treating that null as "not tested" would drop the
# strongest results in the data.
#
# 2.4 adds `pvalue_adjusted` and `pvalue_adjustment` to every burden object,
# and `mirrors/burden.tsv` gains its second study (PMID:34324492, CNV deletion
# burden by permutation). Additive, so MINOR: both keys are always present,
# `null` where the study published no correction.
#
# The display warning matters more than the parsing one, again. A raw p and a
# corrected p can point opposite ways -- CHD7 in PMID:34324492 is 0.0068 raw and
# 0.991 family-wise corrected -- so a consumer rendering `pvalue` alone will
# show a result as significant that the study reported as null. Where
# `pvalue_adjusted` is present it is the number to read.
# 2.5 adds `count_unit` to every burden object. Additive, so MINOR: the key is
# always present, and every row published before it carried `individuals`.
#
# It is the third display warning in a row, and the most load-bearing of them,
# because unlike the two above it changes what the *existing* keys mean.
# `n_case_carriers`, `n_cases`, `n_control_carriers` and `n_controls` are not
# comparable across rows with different `count_unit`: `individuals` counts
# people, `alleles` counts alleles -- so a person with two qualifying variants
# counts twice and the denominator is roughly twice the sample size -- and
# `de_novo_mutations` counts mutations against a denominator of *trios*. A
# consumer that divides `n_case_carriers` by `n_cases` across studies without
# reading this key is comparing three different quantities. That is not
# hypothetical: PMID:42230622 and PMID:40127276 both publish case-control rows
# for the same genes, in different units.
# 2.6 adds `independent_datasets` to every gene bundle and every
# `genes/index.json` row. Additive, so MINOR: the key is always present, and
# `families` always carries one entry per cohort family in the corpus -- a gene
# no study reported gets `not_tested` entries, never a shorter list.
#
# **It is a count of datasets, not a validity call, and the display obligation
# is that it must never be rendered as one.** `headline_confidence` is a
# mirrored ClinGen classification; this object sits beside it and says how many
# independent cohort families showed enrichment. A consumer that renders
# "0 of 2" as a verdict beside a green `definitive` chip will tell a clinician
# the data contradict the classification. They do not -- KDM6A causes Kabuki
# syndrome and shows nothing in either dataset that tested it, because burden
# tests at these cohort sizes routinely miss genes with overwhelming family
# evidence.
#
# Two consequences for any consumer: read `tested` as the denominator, never
# `len(families)`; and render `not_tested` distinguishably from
# `no_enrichment`, because collapsing them is what turns "nobody looked" into
# "somebody looked and found nothing".
#
# 2.7 adds `genes`, `burden_rows` and `cohort_families` to `counts`. Additive,
# so MINOR: every 2.6 key is present and unchanged.
#
# It is a correctness fix, not a convenience. Every `counts` key through 2.6
# counted a `corpus` collection, so the object described *the curation* while
# reading as a census of the build -- measured 2026-08-05 it published
# `assertions: 1, datasets: 0, functional: 0` for a site publishing 23 genes and
# 290 burden statistics, and named neither. A consumer summarising the atlas from
# this object reported an empty one. `index.html` had the same omission and is
# fixed in the same release, because `docs/data-api.md` describes the page as
# stating "every key of `counts`" -- they are one census in two formats and must
# not be able to disagree.
#
# So `counts` is now a census of what the build *published*, not of what the
# curation holds. The distinction matters for reading a diff against 2.6:
# `genes` is D21's published population (23), not `len(corpus.assertions)`'s
# genes (1) and not the 154 genes either validity mirror covers.
#
# `burden_rows` counts rows reaching a published bundle -- 290 -- not the 1,475
# in `mirrors/burden.tsv`. The mirror is deliberately wider than the publication
# gate, so 127 of its 150 genes publish no page; a count of the mirror here would
# advertise evidence no consumer of this API can fetch. It is exactly the sum of
# every `genes/index.json` row's `burden_row_count`, and a consumer may check it
# that way.
#
# `cohort_families`, and deliberately **not** `independent_datasets` even though
# that is the key every gene bundle carries. `counts` already has `datasets`,
# which means *omics* datasets and is 0; two adjacent keys reading `datasets: 0`
# and `independent_datasets: 3` is a misreading waiting to happen. This is the
# length of any bundle's `independent_datasets.families` array, which is the same
# for every gene by construction.
SCHEMA_VERSION: Final = "2.7"

# What `status` publishes today. A literal rather than something derived from
# the corpus, unlike every field in `counts`: there is no measurement of "is
# this atlas ready for clinical use" to compute, only a decision to record, and
# recording it as a constant is what stops a future release forgetting to flip
# it — the two-word status and the paragraph on `index.html` explaining what it
# concretely means are edited in the same commit, or neither is trustworthy.
STATUS: Final = "in-development"


def source_commit(root: Path) -> str | None:
    """The commit the build was made from, or None when there is not one.

    None rather than a raise: building from an unpacked tarball is legitimate
    and should produce a complete site, just one that cannot state its
    provenance. The same answer covers a `root` that does not exist and a
    checkout with no commits yet — `git -C` failing to change directory records
    exactly as much provenance as a missing `.git` does, so neither earns a
    branch of its own here.

    The `returncode` check is load-bearing rather than defensive, which is not
    obvious: in a repository with no commits `git rev-parse HEAD` exits 128 but
    prints the literal string `HEAD` on *stdout* (measured, git 2.50). Without
    the check this would publish `"source_commit": "HEAD"` — a provenance claim
    resolving to nothing, which is worse than the null a consumer can test for.

    `check=False` is not a second decision: `check=True` raises
    `CalledProcessError`, a `SubprocessError`, which the handler below already
    catches, so the two spellings behave identically. `timeout` is a real one —
    without it a git that blocks on a dead network mount or a credential prompt
    hangs the build with nothing to read.

    The toplevel comparison is what makes "not a checkout" mean it. `git -C X`
    searches *upward* from X, so a `--root` that is merely unpacked inside some
    unrelated repository resolves that repository's HEAD, and this published a
    real, verifiable sha naming a commit containing none of the built data —
    strictly worse than the `"HEAD"` the returncode check exists to prevent,
    because a false provenance claim that resolves is one nobody thinks to
    check. Two environment variables do the same thing more quietly: `GIT_DIR`
    and `GIT_WORK_TREE` override `-C` outright, so they are cleared rather than
    trusted. Reproduced before fixing, and pinned by
    `test_a_root_inside_an_unrelated_checkout_has_no_provenance`.
    """
    # `-C` still points git at `root`; `--show-toplevel` then reports which
    # repository it actually found, which is the answer that has to be checked.
    environment = {key: value for key, value in os.environ.items() if not key.startswith("GIT_")}
    try:
        located = subprocess.run(
            ["git", "-C", str(root), "rev-parse", "--show-toplevel", "HEAD"],
            capture_output=True,
            text=True,
            check=False,
            timeout=10,
            env=environment,
        )
    # OSError is git absent from PATH, which a build must survive for the same
    # reason a tarball build must.
    except (OSError, subprocess.SubprocessError):
        return None
    if located.returncode != 0:
        return None

    # Line-based, not whitespace-based. `--show-toplevel` prints a path, and a
    # path may contain spaces — `~/My Projects/atlas` is ordinary on macOS.
    # Splitting on whitespace made that three tokens instead of two and returned
    # `None`, so a perfectly valid checkout published no provenance and did it
    # by failing into the same answer that means "not a checkout at all".
    # Reported by review on #7 and reproduced before fixing.
    #
    # The commit is the last line and never contains whitespace; everything
    # before it is the path, rejoined, so even a newline in a directory name
    # compares correctly rather than silently answering `None`.
    lines = located.stdout.splitlines()
    if len(lines) < 2:
        return None
    toplevel, commit = "\n".join(lines[:-1]), lines[-1]
    try:
        if Path(toplevel).resolve() != root.resolve():
            return None
    except OSError:
        return None
    # `or None`, because "" is falsy but still reads as a recorded commit to a
    # consumer testing the key against null.
    return commit.strip() or None


def write_manifest(
    corpus: Corpus,
    emitter: Emitter,
    commit: str | None,
    counts: Mapping[str, int],
) -> None:
    """Write `manifest.json` and seal the emitter. The build ends here.

    `emitter.checksums` is complete only once every other builder has run, so
    this must be the last thing a build emits. Sealing is what makes that a rule
    rather than an intention — see `Emitter.seal`.

    `counts` carries the figures this module cannot derive from a `Corpus`: the
    published gene population, the burden rows reaching a bundle, and the number
    of cohort families. They are passed in rather than recomputed for the reason
    `runner.py` threads `published` to three builders — the front page, the
    browse payload and this manifest must state one census, and the only way
    they cannot disagree is to count the same objects. Merged over the corpus
    counts rather than nested under a key of their own, because a consumer
    reading `counts` wants the census, not a lesson about which half of it came
    from where.
    """
    corpus_counts = {
        "assertions": len(corpus.assertions),
        "datasets": len(corpus.datasets),
        "featured": len(corpus.featured),
        "functional": len(corpus.functional),
        "phenotypes": len(corpus.phenotypes),
        "publications": len(corpus.publications),
    }
    # A build count silently shadowing a corpus count would republish one figure
    # under another's name, and `sort_keys` would leave the output looking
    # perfectly ordinary. `datasets` is the live hazard: it means *omics*
    # datasets here and a caller reaching for the cohort-family count is one
    # plausible keystroke from overwriting it. `raise`, never `assert` -- `-O`
    # strips `assert`, and this is a guard on a mistake no test would otherwise
    # see.
    collisions = sorted(set(counts) & set(corpus_counts))
    if collisions:
        raise ValueError(
            f"build counts {collisions} would shadow the corpus counts of the same name; "
            f"rename them rather than overwriting a figure derived from the corpus"
        )

    emitter.write_json(
        "manifest.json",
        {
            "schema_version": SCHEMA_VERSION,
            "source_commit": commit,
            "status": STATUS,
            "counts": {**corpus_counts, **counts},
            # A snapshot, not an ordering. Sorting here was specified and then
            # dropped: `encode_json` passes `sort_keys=True`, which orders keys
            # at every level, so an unsorted mapping encodes to the same bytes.
            # Measured rather than reasoned — restoring the sort changes no
            # published byte and no test can see it — because a sort that buys
            # nothing has already been documented twice in this build as the
            # thing that keeps it reproducible.
            #
            # The copy is equally undetectable and is kept as hygiene, which is
            # the honest description: the payload is encoded before `write_json`
            # records anything, so aliasing `checksums` here would not put
            # `manifest.json` into its own `files` today. It would if that
            # ordering ever changed, and a mapping handed out is one a caller
            # can mutate.
            "files": dict(emitter.checksums),
        },
    )
    emitter.seal()
