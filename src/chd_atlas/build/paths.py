# src/chd_atlas/build/paths.py
"""Turning identifiers into path segments, in exactly one place.

The index and the bundle writer must agree on every filename or the site links to
files that were never written, so neither constructs a path itself.
"""

from __future__ import annotations

import re
from typing import Final

# Colon is illegal in a Windows filename and reserved in a URL authority; a slash
# would silently become a directory separator rather than fail, and
# `Zygosity.NOT_APPLICABLE` is the literal value "n/a". The dot is kept: it is
# legal in a filename, and `Emitter._write` is what rejects a path segment that
# is exactly "." or "..".
_UNSAFE: Final = re.compile(r"[^A-Za-z0-9.-]")


def slug(value: str) -> str:
    r"""A path-safe form of one identifier.

    Injective over the *ASCII* values of every grammar in `identifiers.py` but
    one, so two distinct identifiers cannot land on a single filename. Within
    ASCII the only characters this rewrites are the colon and the underscore,
    and no two values differ solely by swapping the two: the underscore appears
    only in `ContrastId`, which admits no capital, while every grammar bearing a
    colon starts with one. The hyphen and the dot are left alone, which is what
    keeps that list of two complete.

    The excepted grammar is `Doi` (`^10\.\d{4,9}/\S+$`), whose suffix admits
    every non-whitespace character, so "10.1000/ab/cd", "10.1000/ab:cd" and
    "10.1000/ab_cd" all reduce to one slug. Nothing builds a path from a DOI, and
    anything that starts to needs a rule that encodes rather than substitutes.

    "Within ASCII" is load-bearing rather than pedantic. `\d` is Unicode-aware in
    both Python's `re` and the rust engine pydantic validates with, so
    "HGNC:1160" suffixed with an Arabic-indic one and the same with an
    Arabic-indic two are two distinct values that `HgncId` accepts and that share
    the slug "HGNC_1160_". `\d` appears in every grammar but `ContrastId` —
    alone in spelling its alphabet out as `[a-z0-9]` — so this is not a corner of
    the schema: it reaches `HgncId`, the one grammar `gene_bundle_path` exists to
    serve. Two things keep it from quietly costing a gene page, neither of them
    this function: no curated file contains a non-ASCII digit today, and
    `Emitter._write` raises on the second write to a path rather than replacing
    the first, so a collision fails the build instead of publishing one page
    where two were expected. Note that the first of those is an observation
    about the data, not an invariant anything enforces.

    `test_build_paths.py` pins all three claims, the first of them against at
    least one ASCII value per grammar. That sample is maintained by hand and
    cannot notice a grammar it holds no value for, so a change to
    `identifiers.py` means rechecking the argument above and extending the
    sample, not trusting a green suite.
    """
    return _UNSAFE.sub("_", value)


def gene_bundle_path(gene: str) -> str:
    """Relative URL of one gene's bundle. Forward-slashed: consumers see a URL."""
    return f"genes/{slug(gene)}.json"
