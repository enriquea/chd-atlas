"""Deterministic JSON Schema emission for the YAML record models."""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from typing import Final

from pydantic import BaseModel

from chd_atlas.models.assertion import AssertionFile
from chd_atlas.models.dataset import Dataset
from chd_atlas.models.functional import FunctionalFile
from chd_atlas.models.literature import FeaturedFile, PhenotypeFile, PublicationFile

EXPORTED_MODELS: Final[dict[str, type[BaseModel]]] = {
    "assertion_file": AssertionFile,
    "functional_file": FunctionalFile,
    "publication_file": PublicationFile,
    "featured_file": FeaturedFile,
    "phenotype_file": PhenotypeFile,
    "dataset": Dataset,
}


def _write_atomically(path: Path, payload: str) -> None:
    """Replace ``path`` with ``payload`` in a single step.

    Opening the destination for writing truncates it before the payload is
    written, so a write that failed part-way left one schema truncated while
    the rest of the directory kept its previous content. The drift test
    compares committed bytes against freshly generated ones, so it would report
    that corrupted file as merely stale and send a curator to re-run the export
    rather than telling them the file on disk is not a schema at all.

    Writing a sibling temporary file and renaming it into place leaves the
    target as either its old content or the complete new content.
    """
    # The temporary file is a sibling so that `os.replace` stays within one
    # filesystem, where it is atomic.
    descriptor, name = tempfile.mkstemp(dir=path.parent, prefix=f"{path.name}.", suffix=".tmp")
    temporary = Path(name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(payload)
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def export_schemas(target: Path) -> list[Path]:
    """Write one ``*.schema.json`` per model. Output is stable across runs."""
    target.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for name, model in sorted(EXPORTED_MODELS.items()):
        schema = model.model_json_schema()
        schema["title"] = model.__name__
        path = target / f"{name}.schema.json"
        _write_atomically(path, json.dumps(schema, indent=2, sort_keys=True) + "\n")
        written.append(path)
    return written
