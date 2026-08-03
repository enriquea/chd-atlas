"""Pins that the scope checks (`validate/scope.py`) are actually wired into
`validate_repository`, not merely importable and unit-tested in isolation.

This module exists separately from `test_runner.py` because the plan that
introduced it (`docs/superpowers/plans/2026-08-03-validity-backbone.md`, Task
10) named it explicitly; `test_runner.py` already covers every other branch of
`validate_repository` against synthetic `tmp_path` repositories.
"""

from __future__ import annotations

from pathlib import Path

from chd_atlas.validate.runner import validate_repository

REPO_ROOT = Path(__file__).parent.parent.parent


def test_scope_checks_run_against_the_real_repository() -> None:
    """The committed repository is in scope and reports no scope error.

    A smoke test in the sense CLAUDE.md permits: it is the only thing that
    exercises the wiring from `validate_repository` down through the mirror
    read into `validate_scope_terms`, `scope_candidates` and
    `validate_curation_is_in_scope`.

    Asserting only `report.error_count == 0` is not enough here -- that
    passes identically whether the scope checks ran or were never called.
    Measured directly: with the scope-check call site in `validate_repository`
    reverted entirely, `.venv/bin/python -m pytest` on a naive version of this
    test asserting only `report.error_count == 0` still passed (1 passed) --
    it cannot tell a live wire from a dead one. So this test also asserts that
    a scope check actually produced output.

    SCP003 (a mirrored cardiac-looking disease term that is not on the scope
    list) is emitted whenever such a term exists, and 268 (net 1, label
    keyword) plus 223 (net 2, gene bridge) were measured on 2026-08-03 against
    the committed 68-term `curation/chd_scope.yaml` and the committed
    ClinGen/GenCC mirrors, so it fires today. If the scope list ever grows to
    cover every candidate, replace the `"SCP003" in codes` assertion with one
    that the mirror-reading helper (`validate/runner.py::_mirrored_validity`)
    returned a non-empty result -- do not delete the guard, or this test
    degrades back into the worthless version described above.
    """
    report = validate_repository(REPO_ROOT)

    codes = {issue.code for issue in report.issues}
    assert "SCP003" in codes, "scope checks did not run"
    assert report.error_count == 0, report.render()
