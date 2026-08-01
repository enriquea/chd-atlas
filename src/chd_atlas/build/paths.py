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

    Injective over every grammar in `identifiers.py` but one, so two distinct
    identifiers cannot land on a single filename. Across those grammars the only
    characters this rewrites are the colon and the underscore, and no two of
    their values differ solely by swapping the two: the underscore appears only
    in `ContrastId`, which admits no capital, while every grammar bearing a colon
    starts with one.

    The exception is `Doi` (`^10\.\d{4,9}/\S+$`), whose suffix admits every
    non-whitespace character, so "10.1000/ab/cd", "10.1000/ab:cd" and
    "10.1000/ab_cd" all reduce to one slug. Nothing builds a path from a DOI, and
    anything that starts to needs a rule that encodes rather than substitutes.

    `test_build_paths.py` pins both halves against one value per grammar. That
    sample is maintained by hand and cannot notice a grammar being widened, so
    adding one to `identifiers.py` — or letting an existing one admit a capital
    or a new punctuation mark — means rechecking the argument above rather than
    trusting a green suite to have caught it.
    """
    return _UNSAFE.sub("_", value)


def gene_bundle_path(gene: str) -> str:
    """Relative URL of one gene's bundle. Forward-slashed: consumers see a URL."""
    return f"genes/{slug(gene)}.json"
