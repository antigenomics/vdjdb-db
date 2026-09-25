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
    "motifs": "10-11 (feature/motifs-*)",
    "summary": "12 (feature/summary)",
    "release": "14 (feature/release-tooling)",
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
    table: str = typer.Option("vdjdb", help="Any declared table: vdjdb, vdjdb-web, slim, full, "
                                            "records, chains, evidence, cluster_members, "
                                            "motif_pwms."),
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
def build(
    out: Path = typer.Option(Path("out"), help="Output directory."),
    chunks: Path | None = typer.Option(None, help="Chunk directory; default chunks/."),
    tables: bool = typer.Option(True, help="Write the definitive tables and the new format."),
    legacy: bool = typer.Option(True, help="Write the legacy projection."),
    airr: bool = typer.Option(True, help="Write the AIRR projection."),
    release: str = typer.Option("dev", help="Release tag recorded on new evidence rows."),
) -> None:
    """Assemble the database: the definitive tables, and every format projected from them."""
    from .assemble.master import build_master
    from .assemble.tables import build_tables
    from .emit.airr import from_tables as airr_frames
    from .emit.airr import write_all as write_airr
    from .emit.legacy import write_all as write_legacy
    from .emit.vdjdb3 import write_all as write_new
    from .io.chunks import chunk_files

    paths = chunk_files(chunks) if chunks else None
    built = build_tables(build_master(paths), release=release)
    out.mkdir(parents=True, exist_ok=True)

    if tables:
        for name, frame in built.items():
            typer.echo(f"{name:10} {frame.height:>8,} rows  {len(frame.columns):>3} cols")
        for name, path in write_new(built, out / "tables").items():
            typer.echo(f"{name:22} {path.stat().st_size:>12,} bytes")
    if legacy:
        for name, path in write_legacy(built, out / "legacy").items():
            typer.echo(f"{name:22} {path.stat().st_size:>12,} bytes")
    if airr:
        for name, path in write_airr(airr_frames(built), out / "airr").items():
            typer.echo(f"{name:26} {path.stat().st_size:>12,} bytes")


@app.command()
def make(
    what: str = typer.Argument(..., help="Which projection: legacy."),
    tables: Path = typer.Option(Path("out/tables"), help="A built new-format directory."),
    out: Path | None = typer.Option(None, help="Output directory; default <tables>/../<what>."),
) -> None:
    """Project a shipped format out of an already-built database.

    The point of the `make` split: the legacy export reads the **tables that shipped**, never
    `chunks/`, so it cannot drift into a parallel implementation of the build.
    """
    from .emit.vdjdb3 import read_tables

    if what != "legacy":
        typer.secho(f"unknown projection {what!r}; known: legacy", fg=typer.colors.RED, err=True)
        raise typer.Exit(2)
    from .emit.legacy import write_all as write_legacy

    dest = out or tables.parent / what
    for name, path in write_legacy(read_tables(tables), dest).items():
        typer.echo(f"{name:22} {path.stat().st_size:>12,} bytes")


@app.command()
def rules(
    chunks: Path | None = typer.Option(None, help="Chunk directory; default chunks/."),
    out: Path = typer.Option(Path("rules/expected_diffs.toml"), help="Ledger rule file."),
    report: Path | None = typer.Option(None, help="Also write the harmonisation report as TSV."),
) -> None:
    """Regenerate the ledger's declared renames from the nomenclature harmonisation.

    A nomenclature correction to an identity column removes a row and adds one, so it has no cell to
    attribute. The generated block tells the ledger to apply the same rewrite to the reference before
    keying; what a reviewer reads is this block's diff.
    """
    from .curate.nomenclature import (
        disambiguate_alleles,
        harmonise_segments,
        legacy_resolver,
        write_renames,
    )
    from .curate.patch import apply_antigen_patch
    from .io.chunks import chunk_files, read_chunks

    paths = chunk_files(chunks) if chunks else None
    harmonised, rep = harmonise_segments(apply_antigen_patch(read_chunks(paths)))
    _, alleles = disambiguate_alleles(harmonised)
    n = write_renames(rep, out, legacy_resolver(), alleles)
    typer.echo(f"{n} renames, {rep['rows'].sum():,} spelling + {alleles['rows'].sum():,} allele "
               f"records, written to {out}")
    if report:
        report.parent.mkdir(parents=True, exist_ok=True)
        rep.write_csv(report, separator="\t")
        alleles.write_csv(report.with_name("alleles.tsv"), separator="\t")
        typer.echo(f"report -> {report} and {report.with_name('alleles.tsv')}")


@app.command()
def convert(
    what: str = typer.Argument("airr", help="Target format: airr."),
    tables: Path | None = typer.Option(None, help="A built new-format directory."),
    legacy: Path | None = typer.Option(None, help="A legacy `vdjdb.txt` to convert instead."),
    out: Path = typer.Option(Path("out/airr"), help="Output directory."),
) -> None:
    """Convert the database to another standard, from the tables or from a legacy release.

    Both sources land on the same emitter: legacy `vdjdb.txt` already speaks VDJdb's column names,
    so there is no second implementation to drift. The legacy path exists for users holding an older
    release zip; it carries no `d_call` (the file has no D column) and, on the current corpus, 1,501
    fewer chains, which is exactly what the legacy build drops.
    """
    import polars as pl

    from .emit.airr import from_legacy, from_tables, write_all
    from .emit.vdjdb3 import read_tables

    if what != "airr":
        typer.secho(f"unknown target {what!r}; known: airr", fg=typer.colors.RED, err=True)
        raise typer.Exit(2)
    if (tables is None) == (legacy is None):
        typer.secho("give exactly one of --tables or --legacy", fg=typer.colors.RED, err=True)
        raise typer.Exit(2)

    frames = (from_tables(read_tables(tables)) if tables is not None
              else from_legacy(pl.read_csv(legacy, separator="\t", infer_schema=False,
                                           quote_char=None)))
    for name, path in write_all(frames, out).items():
        typer.echo(f"{name:26} {path.stat().st_size:>12,} bytes")


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
