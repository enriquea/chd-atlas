from __future__ import annotations

from pathlib import Path

import typer

from chd_atlas.build.runner import BuildRefused, build_site
from chd_atlas.schema_export import export_schemas
from chd_atlas.validate.runner import validate_repository

app = typer.Typer(help="CHD Atlas curation tooling.", no_args_is_help=True)
schemas_app = typer.Typer(help="JSON Schema utilities.", no_args_is_help=True)
app.add_typer(schemas_app, name="schemas")


@app.command()
def validate(
    root: Path = typer.Option(Path("."), help="Repository root to validate."),
) -> None:
    """Validate the curation corpus and mirror tables."""
    # Every validator reports its own target as missing, so a --root that does
    # not exist produces output byte-identical to a real but empty repository:
    # a confident list of content errors about a repository that was never
    # read. Exit 2 rather than 1 so CI can distinguish "you pointed me at the
    # wrong place" from "this repository has errors".
    if not root.is_dir():
        typer.echo(f"error: --root {root} is not a directory")
        raise typer.Exit(code=2)

    report = validate_repository(root)
    typer.echo(report.render())
    raise typer.Exit(code=0 if report.ok else 1)


@app.command()
def build(
    root: Path = typer.Option(Path("."), help="Repository root to build from."),
    out: Path = typer.Option(Path("dist"), help="Directory to write the site into."),
) -> None:
    """Build the published data API into a directory."""
    # Checked before the build, for the same reason `validate` checks it: a
    # --root that does not exist otherwise reaches the validator, which reports
    # every target as missing and refuses — so a typo produces a confident list
    # of content errors about a repository that was never read. Exit 2 rather
    # than 1 so CI can tell "you pointed me at the wrong place" from "this
    # repository has errors"; the first is a pipeline bug, the second a
    # curation one, and they are fixed by different people.
    if not root.is_dir():
        typer.echo(f"error: --root {root} is not a directory")
        raise typer.Exit(code=2)

    # `--out` may legitimately not exist, since the build creates it — but if
    # something is already there it has to be a directory. `is_symlink()` is
    # tested first because `exists()` follows the link, so a dangling symlink
    # answers False to both `exists()` and `is_dir()` and would slip past a
    # check written the obvious way.
    #
    # Checked here rather than left to fail during the build for the same
    # reason `--root` is: an argument that cannot work should not cost a second
    # and a half of validation before saying so.
    if (out.is_symlink() or out.exists()) and not out.is_dir():
        typer.echo(f"error: --out {out} exists and is not a directory")
        raise typer.Exit(code=2)

    try:
        written = build_site(root, out)
    except BuildRefused as exc:
        # The rendered report, not just the fact of failure: a curator who ran
        # `build` needs to know what to fix without being sent to a second
        # command to find out. `build_site` writes nothing when it refuses, so
        # there is no partial `dist/` here for a deploy step to pick up.
        typer.echo(str(exc))
        raise typer.Exit(code=1) from exc
    except OSError as exc:
        # A destination that cannot be written is a bad argument, not a defect,
        # and the same exit code `schemas export` uses for the same reason. The
        # guard above catches the shapes that are visible before the build; this
        # catches the rest — a read-only parent, a full disk, a path that
        # becomes unwritable while the build runs. Everything that reads the
        # repository already reports its own failures as validation issues and
        # was gated on above, so an OSError reaching here is a write.
        typer.echo(f"error: could not write the site to {out}: {exc}")
        raise typer.Exit(code=2) from exc

    typer.echo(f"wrote {len(written)} file(s) to {out}")


@schemas_app.command("export")
def schemas_export(
    target: Path = typer.Option(Path("schemas"), help="Directory to write schemas into."),
) -> None:
    """Regenerate committed JSON Schema files from the Pydantic models."""
    # A target that is an existing file, or a directory that cannot be written,
    # is a bad argument rather than a defect: report the path and exit 2 rather
    # than showing the user a traceback. FileExistsError and PermissionError
    # are both OSError.
    try:
        written = export_schemas(target)
    except OSError as exc:
        typer.echo(f"error: could not write schemas to {target}: {exc}")
        raise typer.Exit(code=2) from exc

    typer.echo(f"wrote {len(written)} schema(s) to {target}")
