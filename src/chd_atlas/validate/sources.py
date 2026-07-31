"""The source registry: which external resource, which release, which licence.

Licence provenance is a correctness requirement rather than documentation,
because the atlas redistributes third-party data.
"""

from __future__ import annotations

from datetime import date
from pathlib import Path

from pydantic import BaseModel, ConfigDict, Field, ValidationError, model_validator
from ruamel.yaml import YAML
from ruamel.yaml.error import YAMLError

from chd_atlas.issues import Severity, ValidationIssue


class Source(BaseModel):
    model_config = ConfigDict(extra="forbid")

    id: str = Field(min_length=1)
    name: str = Field(min_length=1)
    version: str = Field(min_length=1)
    retrieved_on: date
    url: str = Field(min_length=1)
    licence: str = Field(min_length=1)
    redistribution: str = Field(min_length=1)
    ontology_prefix: str | None = None
    ontology_file: str | None = None

    @model_validator(mode="after")
    def ontology_prefix_is_pinned_to_a_file(self) -> Source:
        if self.ontology_prefix is not None and self.ontology_file is None:
            raise ValueError(
                f"source '{self.id}' declares ontology_prefix "
                f"'{self.ontology_prefix}' but no ontology_file"
            )
        return self


class SourceRegistry(BaseModel):
    model_config = ConfigDict(extra="forbid")

    sources: list[Source] = Field(min_length=1)

    @model_validator(mode="after")
    def ids_are_unique(self) -> SourceRegistry:
        seen: set[str] = set()
        for source in self.sources:
            if source.id in seen:
                raise ValueError(f"duplicate source '{source.id}'")
            seen.add(source.id)
        return self

    @property
    def by_id(self) -> dict[str, Source]:
        return {source.id: source for source in self.sources}

    def ontology_files(self) -> dict[str, str]:
        """Map ontology prefix (``HP``) to its pinned OBO path."""
        return {
            source.ontology_prefix: source.ontology_file
            for source in self.sources
            if source.ontology_prefix is not None and source.ontology_file is not None
        }


_EMPTY = SourceRegistry.model_construct(sources=[])


def load_sources(root: Path) -> tuple[SourceRegistry, list[ValidationIssue]]:
    path = root / "mirrors" / "sources.yaml"
    if not path.is_file():
        return _EMPTY, [
            ValidationIssue("SRC001", Severity.ERROR, str(path), "source registry not found")
        ]

    yaml = YAML(typ="safe")
    try:
        raw = yaml.load(path.read_text())
    except YAMLError as exc:
        return _EMPTY, [
            ValidationIssue("YAML001", Severity.ERROR, str(path), f"could not parse YAML: {exc}")
        ]

    try:
        return SourceRegistry.model_validate(raw), []
    except ValidationError as exc:
        issues = [
            ValidationIssue(
                "SRC002",
                Severity.ERROR,
                str(path),
                f"{'.'.join(str(part) for part in error['loc'])}: {error['msg']}",
            )
            for error in exc.errors()
        ]
        return _EMPTY, issues


def validate_source_references(
    used: dict[str, set[str]], registry: SourceRegistry
) -> list[ValidationIssue]:
    """Check every ``source`` value used by a mirror table has a registry entry.

    A record whose source is unregistered has no recorded licence, which is a
    redistribution risk rather than a cosmetic omission.
    """
    known = set(registry.by_id)
    issues: list[ValidationIssue] = []
    for location, sources in sorted(used.items()):
        for source in sorted(sources):
            if source not in known:
                issues.append(
                    ValidationIssue(
                        "SRC003",
                        Severity.ERROR,
                        location,
                        f"source '{source}' has no entry in mirrors/sources.yaml "
                        f"(licence unknown)",
                    )
                )
    return issues
