from __future__ import annotations

from dataclasses import dataclass
from enum import StrEnum


class Severity(StrEnum):
    ERROR = "error"
    WARNING = "warning"


@dataclass(frozen=True, order=True)
class ValidationIssue:
    """A single problem found in the curation corpus.

    Ordering is by (code, severity, location, message) so that reports are
    byte-identical across runs regardless of validator execution order.

    `location` should begin with a repo-relative path whenever the issue is
    attributable to a file (e.g. "curation/assertions/TBX5.yaml" or
    "mirrors/variants/12.tsv:row 42"), so that consumers can group issues by
    file. When an issue is not attributable to a file, `location` is free
    text describing the entity instead (e.g. "assertion CHDA:AST:0000001").
    """

    code: str
    severity: Severity
    location: str
    message: str

    def format(self) -> str:
        return f"{self.severity.value.upper()} [{self.code}] {self.location}: {self.message}"

    def __str__(self) -> str:
        return self.format()
