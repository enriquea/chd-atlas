# src/chd_atlas/validate/scope.py
"""Checks the curated CHD scope list against the mirrored validity data.

`curation/chd_scope.yaml` is the atlas's own editorial claim about what counts
as congenital heart disease (see `models/scope.py` for why neither mechanical
rule -- expert panel, MONDO `is_a` closure -- works). MONDO is deliberately
*not* pinned as an ontology the way HPO is: the full release would dwarf the
11 MB HPO pin for a resource this module only ever needs to spot-check a
handful of terms in. So a scope term cannot be resolved the way a phenotype
term is in `validate/ontology.py`; instead it is checked against the mirror
-- `mirrors/clingen_gene_validity.tsv` and `mirrors/gencc_submissions.tsv`,
which already carry a `disease`/`disease_label` pair for every term either
authority curates against.

That substitution has a known, accepted blind spot: it does not catch a term
that ClinGen or GenCC still cites but that MONDO has since obsoleted upstream.
Checking against a resolver would catch that; checking against the mirror
cannot, because the mirror is exactly what an obsoleted-but-still-cited term
would still agree with. This is accepted rather than fixed, because the atlas
publishes each scope term with the label its authority gave it -- not with
MONDO's current name for it -- so the claim stays attributable to that
authority either way, obsolete upstream or not.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping
from typing import Final

from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.models.scope import ScopeEntry

# Deliberately wide. Measured 2026-08-03 against the committed ClinGen mirror:
# of 13 flagship CHD genes' Definitive/Strong disease labels, 12 of 13 --
# Holt-Oram syndrome (TBX5), CHARGE syndrome (CHD7), Kabuki syndrome 1 and 2,
# NOTCH1-related AOS spectrum disorder, cutis laxa AD 1 (ELN),
# blepharocheilodontic syndrome 2 (CTNND1), lymphatic malformation 1 (FLT4),
# Noonan syndrome, Mowat-Wilson syndrome, Cornelia de Lange syndrome -- fail
# every keyword here. This pattern alone cannot find them; see
# `scope_candidates`'s docstring for what covers that gap and what still does
# not.
_CARDIAC_KEYWORDS: Final = re.compile(
    r"heart|cardi|septal|aort|fallot|ventricul|atrial|truncus|valve|"
    r"hypoplastic left|heterotax|situs|coarctation|ductus|pulmonary atresia|"
    r"conotruncal|vascular ring",
    re.IGNORECASE,
)


def validate_scope_terms(
    scope: Iterable[ScopeEntry],
    mirrored: Mapping[str, str] | None,
    location: str,
) -> list[ValidationIssue]:
    """Check every scope term resolves in the mirror, under the label given.

    `mirrored` maps a MONDO term to the label its authority (ClinGen or GenCC)
    published for it; the caller builds this once from both mirror tables
    rather than this function reading them, so it stays a pure check over
    already-loaded data, matching `validate_source_references` and
    `validate_terms` elsewhere in this package.

    `mirrored is None` means neither mirror table could be read at all, not
    that scope has no terms to check against. Checking every entry against an
    empty mapping would report SCP001 once per scope term -- 68 failures
    naming the symptom, none naming the cause -- exactly the cascade
    REF000/SRC000/ONT000 exist to prevent elsewhere in this package. So a
    missing mirror is reported once, as SCP000, and nothing else runs.

    SCP001 (term absent from the mirror) alone cannot catch a copy-pasted
    wrong MONDO id that happens to exist in the mirror under an unrelated
    label -- it would resolve, and the wrong disease would be silently
    admitted. SCP002 exists for that: the label is the only field a human
    reviewing the scope file actually reads, so it is what gets checked
    against the mirror's own label for the same id.
    """
    if mirrored is None:
        return [
            ValidationIssue(
                "SCP000",
                Severity.WARNING,
                location,
                "skipped scope checks: no validity mirror could be read",
            )
        ]

    issues: list[ValidationIssue] = []
    for entry in scope:
        authority_label = mirrored.get(entry.id)
        if authority_label is None:
            issues.append(
                ValidationIssue(
                    "SCP001",
                    Severity.ERROR,
                    location,
                    f"scope term {entry.id} is not in the mirrored ClinGen/GenCC validity data",
                )
            )
            continue
        if entry.label != authority_label:
            issues.append(
                ValidationIssue(
                    "SCP002",
                    Severity.ERROR,
                    location,
                    f"scope term {entry.id} is labelled '{entry.label}' here, "
                    f"but the mirror calls it '{authority_label}'",
                )
            )
    return issues


def _describe(term: str, mirrored: Mapping[str, str]) -> str:
    label = mirrored.get(term)
    return f"{term} ({label})" if label else term


def scope_candidates(
    scope_terms: set[str],
    mirrored: Mapping[str, str],
    diseases_by_gene: Mapping[str, set[str]],
    in_scope_genes: set[str],
    location: str,
) -> list[ValidationIssue]:
    """Report disease terms that look like CHD but are not on the scope list.

    Reported only -- never auto-admitted. `curation/chd_scope.yaml` is a
    curator's editorial claim (see `models/scope.py`); this function exists so
    scope drift becomes a reviewed queue instead of silent staleness, not so
    it can decide scope on anyone's behalf.

    Two independent nets, because one alone is measurably blind to the class
    of entity this check exists to catch. A label-keyword net cannot surface
    a syndromic CHD entity, because a syndrome's own name usually carries no
    cardiac word at all -- see `_CARDIAC_KEYWORDS` above for the 2026-08-03
    measurement that motivated adding the second net.

    Net 1 ("label keyword") matches every `mirrored` label not already in
    `scope_terms` against `_CARDIAC_KEYWORDS`. It errs wide on purpose: a
    false positive costs a curator one rejection, a false negative costs a
    gene that is never reviewed at all. Measured 2026-08-03 by unioning
    `mirrors/clingen_gene_validity.tsv` and `mirrors/gencc_submissions.tsv`
    into one `disease -> disease_label` mapping (8,154 distinct terms) and
    applying this pattern against everything not among the 68 admitted scope
    terms: 268 candidates. The exact count depends on how the caller
    deduplicates the two mirrors, which is this function's caller's job, not
    this function's -- but it is hundreds either way, which is why the return
    is a digest and not one issue per hit.

    Net 2 ("gene bridge") looks instead at genes already admitted: for every
    gene in `in_scope_genes`, every *other* disease term it is curated under
    (via `diseases_by_gene`) that is not already in `scope_terms`. This is the
    drift case that actually occurs in practice -- an already-admitted gene
    gains a new curation. Measured example: GATA6 is Definitive for
    GATA6-related CHD (in scope) and Limited for dilated cardiomyopathy (not);
    the latter surfaces here and would be correctly rejected on review.

    Neither net discovers a *new* syndromic entity whose label carries no
    cardiac term and whose gene is not already in scope -- that is this
    function's blind spot, and it is a real one, not a hedge. Holt-Oram
    syndrome is the worked example: before TBX5 was ever admitted, "Holt-Oram
    syndrome" would have passed net 1 silently (the label matches no keyword
    in `_CARDIAC_KEYWORDS`) and net 2 silently (TBX5 was not yet an in-scope
    gene for anything to bridge from). Nothing in this module closes that gap;
    it is why the scope list itself, not this check, is what a curator must
    still read the literature to grow.

    Returns at most two issues, one per net, each carrying its own count and
    its own candidates in sorted order -- one issue per candidate would bury
    every real error behind the other ~90 net-1 hits measured on the
    committed mirror, and the report would stop being a digest.

    WARNING, not ERROR, and deliberately so: `ValidationReport.ok` ignores
    warnings, so an unreviewed new entity never blocks a deploy on its own.
    This is a *gap* warning in the sense `REF013` establishes -- the check
    ran and found real, reportable data that a human still has to act on --
    not a *skip* warning like SCP000, where a check could not run at all.
    `validate/runner.py` documents why the two must not be conflated: a skip
    warning is never allowed to arrive alone, but a gap warning is expected to
    arrive alone by design, because the corpus it describes is still sound and
    publishable exactly as curated.
    """
    net1 = sorted(
        term
        for term, label in mirrored.items()
        if term not in scope_terms and _CARDIAC_KEYWORDS.search(label)
    )
    net2 = sorted(
        {
            disease
            for gene in in_scope_genes
            for disease in diseases_by_gene.get(gene, set())
            if disease not in scope_terms
        }
    )

    issues: list[ValidationIssue] = []
    if net1:
        issues.append(
            ValidationIssue(
                "SCP003",
                Severity.WARNING,
                location,
                f"{len(net1)} disease term(s) not in scope have a cardiac-sounding "
                f"label (net 1, label keyword): "
                f"{', '.join(_describe(term, mirrored) for term in net1)}",
            )
        )
    if net2:
        issues.append(
            ValidationIssue(
                "SCP003",
                Severity.WARNING,
                location,
                f"{len(net2)} disease term(s) not in scope are curated for a gene "
                f"already in scope (net 2, gene bridge): "
                f"{', '.join(_describe(term, mirrored) for term in net2)}",
            )
        )
    return issues


def validate_curation_is_in_scope(
    curated_genes: Mapping[str, str],
    in_scope_genes: set[str],
) -> list[ValidationIssue]:
    """Refuse a curated gene that no scope term admits.

    `curated_genes` maps an HGNC id to the file that curates it -- an
    assertion's file, typically, built by the caller from the corpus rather
    than read here, matching this module's other two checks.

    This is the silent-evidence-loss guard for the scope backbone: without
    it, a curator can write a complete, referentially sound assertion about a
    gene `curation/chd_scope.yaml` never admitted, and it reaches no
    published page while `chd-atlas validate` reports 0 errors and every
    checksum verifies -- this project's characteristic failure, described in
    `CLAUDE.md` and reproduced in exactly this shape by the phospho-join and
    empty-omics-shard defects that motivated writing that description down.

    ERROR, not WARNING, unlike SCP003's unreviewed candidate. SCP003 warns
    about a term nobody has judged yet; this is curator judgement that
    already exists -- the assertion was written, reviewed, and committed --
    and would be silently discarded by the build rather than merely
    unreviewed. The message tells the curator the two ways to resolve it:
    admit the gene's disease term into scope, or remove the record, because
    only a human can decide which of those is correct.
    """
    return [
        ValidationIssue(
            "SCP004",
            Severity.ERROR,
            curated_genes[gene],
            f"gene {gene} is curated but no scope term admits it; admit its "
            f"disease term in curation/chd_scope.yaml, or remove the record",
        )
        for gene in sorted(curated_genes)
        if gene not in in_scope_genes
    ]
