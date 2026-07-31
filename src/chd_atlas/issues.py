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
    """

    code: str
    severity: Severity
    location: str
    message: str

    def format(self) -> str:
        return f"{self.severity.value.upper()} [{self.code}] {self.location}: {self.message}"
