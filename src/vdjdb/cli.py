"""The ``vdjdb`` command.

Subcommands land as their ROADMAP phase does; unimplemented ones exit 2 with the phase that
delivers them rather than a traceback, so ``vdjdb --help`` is an honest map of what works today.
"""
from __future__ import annotations

from pathlib import Path

import typer

from . import __version__
from .config import Paths

app = typer.Typer(
    name="vdjdb",
    help="Build, proofread and release the VDJdb database.",
    no_args_is_help=True,
    add_completion=False,
)

# subcommand -> the ROADMAP phase that implements it
_PENDING = {
    "build": "4 (feature/pipeline-core)",
    "motifs": "10-11 (feature/motifs-*)",
    "summary": "12 (feature/summary)",
    "release": "14 (feature/release-tooling)",
    "make": "6 (feature/new-format)",
    "convert": "7 (feature/airr)",
    "refs": "12 (feature/summary)",
    "changelog": "14 (feature/release-tooling)",
}


def _pending(name: str) -> None:
    phase = _PENDING[name]
    typer.secho(
        f"`vdjdb {name}` is not implemented yet — ROADMAP phase {phase}.",
        fg=typer.colors.YELLOW,
        err=True,
    )
    raise typer.Exit(2)


@app.command()
def version() -> None:
    """Print the package version and the resolved repo root."""
    typer.echo(f"vdjdb-db {__version__}")
    typer.echo(f"root {Paths.discover().root}")


@app.command()
def qc(
    paths: list[Path] = typer.Argument(None, help="Chunk files; default: every chunk in chunks/."),
    strict: bool = typer.Option(True, help="Exit non-zero on any QC error."),
    report: Path | None = typer.Option(None, help="Write a TSV error report here."),
) -> None:
    """Validate chunk files."""
    from .qc.runner import run_qc

    raise typer.Exit(run_qc(paths or None, strict=strict, report=report))


@app.command()
def schema(
    table: str = typer.Option("vdjdb", help="vdjdb, vdjdb-web, slim, full, "
                                            "cluster_members or motif_pwms."),
    format: str = typer.Option("meta", help="meta, header or json."),
) -> None:
    """Render a table's metadata, header or JSON schema from the field registry."""
    import json as _json

    from .schema import TABLES, fields, header, render_meta, render_slim_meta

    if table not in TABLES:
        typer.secho(f"unknown table {table!r}; known: {', '.join(sorted(TABLES))}",
                    fg=typer.colors.RED, err=True)
        raise typer.Exit(2)
    if format == "header":
        typer.echo(header(table))
    elif format == "meta":
        # slim's metadata is two columns; every other table's is the eight-column form.
        renderer = render_slim_meta if table == "slim" else render_meta
        typer.echo(renderer(table), nl=False)
    elif format == "json":
        typer.echo(_json.dumps([f.__dict__ if hasattr(f, "__dict__") else
                                {k: getattr(f, k) for k in f.__slots__}
                                for f in fields(table)], indent=2))
    else:
        typer.secho(f"unknown format {format!r}; known: meta, header, json",
                    fg=typer.colors.RED, err=True)
        raise typer.Exit(2)


@app.command()
def build(out: Path = typer.Option(Path("out"))) -> None:
    """Assemble the database into the three output formats."""
    _pending("build")


@app.command()
def motifs(out: Path = typer.Option(Path("out"))) -> None:
    """Infer TCRNET and TCREMP motifs."""
    _pending("motifs")


@app.command()
def summary(out: Path = typer.Option(Path("out"))) -> None:
    """Render the static and interactive dashboards."""
    _pending("summary")


@app.command()
def diff(
    reference: Path = typer.Argument(..., help="Reference release zip or directory."),
    candidate: Path = typer.Argument(..., help="Candidate build directory or zip."),
    rules: Path = typer.Option(Path("rules/expected_diffs.toml"),
                               help="Declared expected differences."),
    report: Path | None = typer.Option(None, help="Write the ledger here as Markdown."),
    only: str | None = typer.Option(None, help="Comma-separated members to compare; "
                                               "for deliberately partial candidates."),
) -> None:
    """Compare a candidate build against a released bundle and attribute every difference."""
    from .compare.diff import diff as run_diff
    from .compare.diff import render

    result = run_diff(reference, candidate, rules if rules.exists() else None,
                      only=[s.strip() for s in only.split(",")] if only else None)
    text = render(result)
    if report:
        report.parent.mkdir(parents=True, exist_ok=True)
        report.write_text(text)
    typer.echo(text, nl=False)
    raise typer.Exit(0 if result.ok else 1)


@app.command()
def release(tag: str = typer.Option(..., help="Release tag, e.g. v2026.09.0.")) -> None:
    """Assemble the release bundles, checksums and manifest."""
    _pending("release")


if __name__ == "__main__":
    app()
