# src/chd_atlas/validate/ontology.py
"""Resolves ontology CURIEs against releases pinned in ``mirrors/sources.yaml``.

Pinning matters: an unpinned ontology upgrade could silently invalidate
curation, or silently start accepting terms that did not exist when a claim
was made.

NCBITaxon is deliberately not loaded as an ontology. The full release is
hundreds of megabytes to validate a handful of model organisms, so taxa
resolve against the ``MODEL_ORGANISMS`` allowlist in ``vocab.py`` instead.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path

import pronto

from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.vocab import MODEL_ORGANISMS

_TAXON_PREFIX = "NCBITaxon"


@dataclass
class OntologyRegistry:
    ontologies: dict[str, pronto.Ontology] = field(default_factory=dict)
    load_issues: list[ValidationIssue] = field(default_factory=list)

    @classmethod
    def from_files(cls, root: Path, pinned: dict[str, str]) -> OntologyRegistry:
        registry = cls()
        for prefix, relative in sorted(pinned.items()):
            path = root / relative
            if not path.is_file():
                registry.load_issues.append(
                    ValidationIssue(
                        "ONT004",
                        Severity.ERROR,
                        str(path),
                        f"pinned ontology for prefix '{prefix}' not found",
                    )
                )
                continue
            registry.ontologies[prefix] = pronto.Ontology(str(path))
        return registry

    def has_prefix(self, prefix: str) -> bool:
        return prefix == _TAXON_PREFIX or prefix in self.ontologies


def validate_terms(
    curies: Iterable[str], registry: OntologyRegistry, location: str
) -> list[ValidationIssue]:
    issues: list[ValidationIssue] = []
    for curie in curies:
        prefix, _, _ = curie.partition(":")

        if prefix == _TAXON_PREFIX:
            if curie not in MODEL_ORGANISMS:
                issues.append(
                    ValidationIssue(
                        "ONT001",
                        Severity.ERROR,
                        location,
                        f"{curie} is not an allowed model organism",
                    )
                )
            continue

        ontology = registry.ontologies.get(prefix)
        if ontology is None:
            issues.append(
                ValidationIssue(
                    "ONT003",
                    Severity.ERROR,
                    location,
                    f"no pinned ontology for prefix '{prefix}' (term {curie})",
                )
            )
            continue

        term = ontology.get_term(curie) if curie in ontology else None
        if term is None:
            issues.append(
                ValidationIssue(
                    "ONT001",
                    Severity.ERROR,
                    location,
                    f"{curie} not found in the pinned {prefix} release",
                )
            )
        elif term.obsolete:
            issues.append(
                ValidationIssue(
                    "ONT002",
                    Severity.ERROR,
                    location,
                    f"{curie} is obsolete in the pinned {prefix} release",
                )
            )
    return issues
