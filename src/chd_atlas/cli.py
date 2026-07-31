from __future__ import annotations

from pathlib import Path

import typer

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
    report = validate_repository(root)
    typer.echo(report.render())
    raise typer.Exit(code=0 if report.ok else 1)


@schemas_app.command("export")
def schemas_export(
    target: Path = typer.Option(Path("schemas"), help="Directory to write schemas into."),
) -> None:
    """Regenerate committed JSON Schema files from the Pydantic models."""
    written = export_schemas(target)
    typer.echo(f"wrote {len(written)} schema(s) to {target}")
