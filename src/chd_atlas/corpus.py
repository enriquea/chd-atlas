# src/chd_atlas/corpus.py
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from pydantic import BaseModel, ValidationError
from ruamel.yaml import YAML
from ruamel.yaml.error import YAMLError

from chd_atlas.issues import Severity, ValidationIssue
from chd_atlas.models.assertion import AssertionFile, GeneDiseaseAssertion
from chd_atlas.models.dataset import Dataset
from chd_atlas.models.functional import FunctionalEvidence, FunctionalFile
from chd_atlas.models.literature import (
    FeaturedFile,
    FeaturedManuscript,
    PhenotypeFile,
    PhenotypeTerm,
    Publication,
    PublicationFile,
)


@dataclass(frozen=True)
class Corpus:
    """Every interpretive record loaded from ``curation/``."""

    root: Path
    assertions: tuple[GeneDiseaseAssertion, ...] = ()
    functional: tuple[FunctionalEvidence, ...] = ()
    publications: tuple[Publication, ...] = ()
    featured: tuple[FeaturedManuscript, ...] = ()
    phenotypes: tuple[PhenotypeTerm, ...] = ()
    datasets: tuple[Dataset, ...] = ()


@dataclass
class _Accumulator:
    issues: list[ValidationIssue] = field(default_factory=list)

    def error(self, code: str, location: Path | str, message: str) -> None:
        self.issues.append(
            ValidationIssue(
                code=code,
                severity=Severity.ERROR,
                location=str(location),
                message=message,
            )
        )


def _read_yaml(path: Path, acc: _Accumulator) -> Any | None:
    yaml = YAML(typ="safe")
    try:
        return yaml.load(path.read_text())
    except YAMLError as exc:
        acc.error("YAML001", path, f"could not parse YAML: {exc}")
        return None


def _parse[ModelT: BaseModel](model: type[ModelT], path: Path, acc: _Accumulator) -> ModelT | None:
    raw = _read_yaml(path, acc)
    if raw is None:
        return None
    try:
        return model.model_validate(raw)
    except ValidationError as exc:
        for error in exc.errors():
            location = ".".join(str(part) for part in error["loc"])
            acc.error("SCHEMA001", path, f"{location}: {error['msg']}")
        return None


def load_curation(root: Path) -> tuple[Corpus, list[ValidationIssue]]:
    """Load every YAML record under ``root/curation``.

    Errors are accumulated rather than raised so that a single run reports
    every problem in the corpus, not just the first one encountered.
    """
    acc = _Accumulator()
    curation = root / "curation"
    if not curation.is_dir():
        acc.error("CORPUS001", curation, "curation directory not found")
        return Corpus(root=root), acc.issues

    assertions: list[GeneDiseaseAssertion] = []
    for path in sorted((curation / "assertions").glob("*.yaml")):
        parsed = _parse(AssertionFile, path, acc)
        if parsed is not None:
            assertions.extend(parsed.assertions)

    functional: list[FunctionalEvidence] = []
    for path in sorted((curation / "functional").glob("*.yaml")):
        parsed_functional = _parse(FunctionalFile, path, acc)
        if parsed_functional is not None:
            functional.extend(parsed_functional.functional_evidence)

    datasets: list[Dataset] = []
    for path in sorted((curation / "datasets").glob("*.yaml")):
        parsed_dataset = _parse(Dataset, path, acc)
        if parsed_dataset is not None:
            datasets.append(parsed_dataset)

    publications: tuple[Publication, ...] = ()
    publications_path = curation / "publications.yaml"
    if publications_path.is_file():
        parsed_publications = _parse(PublicationFile, publications_path, acc)
        if parsed_publications is not None:
            publications = tuple(parsed_publications.publications)

    featured: tuple[FeaturedManuscript, ...] = ()
    featured_path = curation / "featured.yaml"
    if featured_path.is_file():
        parsed_featured = _parse(FeaturedFile, featured_path, acc)
        if parsed_featured is not None:
            featured = tuple(parsed_featured.featured)

    phenotypes: tuple[PhenotypeTerm, ...] = ()
    phenotypes_path = curation / "phenotypes.yaml"
    if phenotypes_path.is_file():
        parsed_phenotypes = _parse(PhenotypeFile, phenotypes_path, acc)
        if parsed_phenotypes is not None:
            phenotypes = tuple(parsed_phenotypes.phenotypes)

    corpus = Corpus(
        root=root,
        assertions=tuple(assertions),
        functional=tuple(functional),
        publications=publications,
        featured=featured,
        phenotypes=phenotypes,
        datasets=tuple(datasets),
    )
    return corpus, acc.issues
