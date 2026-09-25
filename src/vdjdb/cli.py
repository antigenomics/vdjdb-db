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
    "diff": "2 (feature/golden-harness)",
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
def diff(against: str = typer.Option(..., help="Release tag or path to a reference zip.")) -> None:
    """Produce the difference ledger against a released build."""
    _pending("diff")


@app.command()
def release(tag: str = typer.Option(..., help="Release tag, e.g. v2026.09.0.")) -> None:
    """Assemble the release bundles, checksums and manifest."""
    _pending("release")


if __name__ == "__main__":
    app()
