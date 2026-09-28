"""The ``vdjdb`` command.

Every subcommand is implemented; ``vdjdb --help`` lists them all.
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
    engine: str = typer.Option("arda", help="CDR3 markup engine: arda or legacy."),
) -> None:
    """Assemble the database: the definitive tables, and every format projected from them."""
    from .assemble.master import build_master
    from .assemble.tables import build_tables
    from .emit.airr import from_tables as airr_frames
    from .emit.airr import write_all as write_airr
    from .emit.legacy import write_all as write_legacy
    from .emit.vdjdb3 import write_all as write_new
    from .io.chunks import chunk_files
    from .timing import report as timing_report
    from .timing import reset as timing_reset
    from .timing import stage
    from .timing import write as timing_write

    timing_reset()
    paths = chunk_files(chunks) if chunks else None
    with stage("assemble.master.build_master"):
        master = build_master(paths, engine=engine)
    built = build_tables(master, release=release)
    out.mkdir(parents=True, exist_ok=True)

    # Written every run, next to the other reports: a wall time nobody records is a wall time nobody
    # can regress against, which is how the 156 s in `add_junction_nt` went unmeasured until someone
    # profiled it by hand (ROADMAP_local section 49).
    timings = timing_write(out / "reports" / "build-timings.tsv", rows=built["records"].height)
    typer.echo(timing_report(timings))

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

    The legacy export reads the tables that shipped, never `chunks/`, so it cannot drift into a
    second implementation of the build.
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
    out: Path = typer.Option(Path("rules/expected_diffs.toml"), help="Expected-difference rules."),
    report: Path | None = typer.Option(None, help="Also write the harmonisation report as TSV."),
    reference: Path | None = typer.Option(None, help="Reference bundle the comparison runs against; "
                                                     "antigen renames are derived from it."),
) -> None:
    """Regenerate the declared renames from the nomenclature harmonisation.

    A nomenclature correction to an identity column removes a row and adds one, so it has no cell to
    attribute. The generated block tells `vdjdb diff` to apply the same rewrite to the reference
    before keying; a reviewer reads this block's diff.
    """
    import polars as pl

    from .curate.nomenclature import (
        disambiguate_alleles,
        harmonise_mhc,
        harmonise_references,
        harmonise_segments,
        legacy_resolver,
        write_renames,
    )
    from .curate.patch import apply_antigen_patch
    from .io.chunks import chunk_files, read_chunks

    paths = chunk_files(chunks) if chunks else None
    raw = read_chunks(paths)
    patched = apply_antigen_patch(raw)
    harmonised, rep = harmonise_segments(patched)
    allele_fixed, alleles = disambiguate_alleles(harmonised)
    mhc_fixed, mhc = harmonise_mhc(allele_fixed)
    _, refs = harmonise_references(mhc_fixed)
    if not refs.is_empty():
        mhc = pl.concat([mhc, refs.select(pl.lit("#347").alias("issue"),
                                         pl.lit("reference.id").alias("column"),
                                         "from", "to", "rows")], how="vertical")
    from .compare.diff import Bundle, _read_table
    from .curate.patch import render_patch_renames

    antigen_block = ""
    if reference is not None:
        ref_table = _read_table(Bundle(reference).read_bytes("vdjdb.txt"))
        antigen_block = render_patch_renames(ref_table, patched)
    n = write_renames(rep, out, legacy_resolver(), alleles, mhc, (antigen_block,))
    typer.echo(f"{n} renames, {rep['rows'].sum():,} spelling + {alleles['rows'].sum():,} allele "
               f"+ {mhc['rows'].sum():,} MHC records, written to {out}")
    if report:
        report.parent.mkdir(parents=True, exist_ok=True)
        rep.write_csv(report, separator="\t")
        alleles.write_csv(report.with_name("alleles.tsv"), separator="\t")
        mhc.write_csv(report.with_name("mhc.tsv"), separator="\t")
        typer.echo(f"reports -> {report.parent}/")


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
    release zip; it has no `d_call` (the file has no D column) and, on the current corpus, 1,501
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
def motifs(
    # `out/motifs`, not `out`: that is where `release/manifest.py`, `vdjdb motif-metrics` and
    # `tests/release/test_motif_reproduction.py` all read the motif files from, so the default used
    # to put them somewhere no consumer looked unless `--out` was passed.
    out: Path = typer.Option(Path("out/motifs"), help="Where to write the motif files."),
    reports: Path = typer.Option(Path("out/reports"),
                                 help="Where the two diagnostics go: the per-epitope breakdown and "
                                      "the stage timings. One reports directory for the whole "
                                      "build, per docs/outputs.md."),
    tables: Path = typer.Option(Path("out/tables"), help="The definitive tables to infer from."),
    p: float = typer.Option(None, help="Enrichment p threshold (default: per-chain tcrnet.TUNED)."),
    methods: str = typer.Option("tcrnet,tcremp", help="Which methods to run."),
) -> None:
    """Infer TCRNET and TCREMP motifs: cluster_members*.txt and motif_pwms*.txt."""
    from .motifs import run
    written = run(tables, out, reports, p=p,
                  methods=tuple(m.strip() for m in methods.split(",")))
    for name, rows in written.items():
        typer.echo(f"{name:24} {rows:>9,} rows")


@app.command(name="motif-metrics")
def motif_metrics(
    tables: Path = typer.Option(Path("out/tables"), help="The tables the cohort is read from."),
    motifs: Path = typer.Option(Path("out/motifs"), help="This build's motif files."),
    legacy: Path = typer.Option(..., help="The last legacy release: a zip or extracted directory."),
    latest: Path | None = typer.Option(
        None, help="The latest release. Defaults to --legacy, which is the truth until a release "
                   "ships after 2026-06-03."),
    out: Path = typer.Option(Path("out/reports"), help="Where the report is written."),
    baseline: Path = typer.Option(Path("rules/motif_metrics.tsv"),
                                  help="Committed baseline to gate against."),
    record: bool = typer.Option(
        False, help="Overwrite the baseline with this run instead of gating against it. For a "
                    "reviewed commit that explains the move -- never for a build."),
) -> None:
    """Score the motif clustering of three sources on this build's cohort, and gate both ways.

    The three are always the same: the last legacy release, the latest release, and this build. A
    fourth row is the partition that clusters nothing, because a bar it passes is not a bar.

    Exits 1 on an undeclared regression, naming the axis, both values and the delta.
    """
    import polars as pl

    from .compare.diff import Bundle
    from .validate import motif_bench as mb
    from .validate import motif_metrics as mmv

    chains = pl.read_parquet(tables / "chains.parquet")
    records = pl.read_parquet(tables / "records.parquet")

    def released(path: Path) -> pl.DataFrame:
        b = Bundle(path)
        dest = out / f"_{path.name}.cluster_members.txt"
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(b.read_bytes("cluster_members.txt"))
        return mb.read_members(dest)

    cur = {m: mb.read_members(motifs / f"cluster_members{suffix}.txt")
           for m, suffix in (("tcrnet", ""), ("tcremp", "_tcremp"))
           if (motifs / f"cluster_members{suffix}.txt").exists()}
    if not cur:
        typer.secho(f"no motif files in {motifs}; run `vdjdb motifs --out {motifs}`",
                    fg=typer.colors.RED, err=True)
        raise typer.Exit(2)

    measured = mmv.measure(chains, records,
                           mmv.sources(released(legacy), released(latest or legacy), cur))
    out.mkdir(parents=True, exist_ok=True)
    measured.write_csv(out / "motif-metrics.tsv", separator="\t")
    base = mmv.load_baseline(baseline)
    (out / "motif-metrics.md").write_text(mmv.report(measured, base) + "\n")
    typer.echo(f"{out / 'motif-metrics.tsv'}  {measured.height} rows")

    if record:
        baseline.parent.mkdir(parents=True, exist_ok=True)
        # Carry the declared reasons forward: re-recording must not silently undeclare a trade that
        # is still being made, which would turn the next run's pass into an unexplained pass.
        keep = base.select(*mmv.KEY, "reason") if "reason" in base.columns else None
        rec = (measured.join(keep, on=list(mmv.KEY), how="left") if keep is not None
               else measured.with_columns(pl.lit(None, dtype=pl.Utf8).alias("reason")))
        rec.select(*mmv.COLUMNS).write_csv(baseline, separator="\t")
        typer.secho(f"baseline rewritten: {baseline}. Say why in the commit message.",
                    fg=typer.colors.YELLOW)
        return

    bad = mmv.regressions_against_latest(measured, base)
    drift = mmv.regressions_against_baseline(measured, base)
    for title, frame in (("worse than the latest release", bad),
                         ("drifted from the committed baseline", drift)):
        if frame.height:
            typer.secho(f"\n{frame.height} gated axes {title}:", fg=typer.colors.RED, err=True)
            typer.echo(frame, err=True)
    if bad.height or drift.height:
        typer.secho("\nA drift that a corpus change explains is accepted by moving "
                    f"{baseline} in a commit that says why (--record).", err=True)
        raise typer.Exit(1)
    typer.secho("no gated axis regressed", fg=typer.colors.GREEN)


@app.command()
def summary(
    legacy: Path = typer.Option(Path("out/legacy"), help="Legacy projection this build produced."),
    reference: Path | None = typer.Option(None, help="A previous fragment, for the SSIM layer."),
    verbose: bool = typer.Option(False, help="Show knitr's chunk-by-chunk progress."),
    assets: Path | None = typer.Option(
        None, help="Write the figures here instead of inlining them. Needs a vdjdb-web change: "
                   "its Scala side finds images by matching data:image/png;base64."),
) -> None:
    """Render the release dashboard and verify the fragment `vdjdb-web` will serve."""
    from .summary import render as r

    r.render(legacy, quiet=not verbose)
    typer.echo(f"{r.extract(assets=assets):,} lines -> {r.FRAGMENT}")
    if code := r.check(reference=reference):
        raise typer.Exit(code)


@app.command()
def refs(
    tables: Path = typer.Option(Path("out/tables"), help="Definitive tables to read references from."),
    out: Path = typer.Option(None, help="Where to write the table; defaults to the committed path."),
) -> None:
    """Resolve the publication year of every reference and write `summary/reference_years.tsv`.

    Network-bound and not part of a build: the table is a committed, reviewed input that makes the
    dashboard render offline, refreshed by its own pull request (hard rule 9).

    Reads the built `records`, not `chunks/`, because the table has to key on the reference id the
    dashboard will see: #347 rewrites a DOI or a publisher URL into a PMID, so resolving the raw chunk
    values drops every reference the harmonisation converts. Measured: three of them, covering 669
    records.

    That leaves a staleness hazard, which the check below closes. Resolving against a build that
    predates a newly landed chunk writes the table without that chunk's reference, and the dashboard
    then fails in CI fifteen minutes later with "1 reference(s) have no year". `records` carries
    `chunk.file`, so the build states which chunks it saw and the comparison is exact rather than a
    guess from modification times, which a branch switch alone would invalidate.
    """
    import polars as pl

    from .io.chunks import chunk_files
    from .summary import references as refs_mod

    built = tables / "records.parquet"
    if not built.exists():
        typer.secho(f"no records table at {built}; run `vdjdb build` first",
                    fg=typer.colors.RED, err=True)
        raise typer.Exit(2)
    records = pl.read_parquet(built)
    seen = set(records["chunk.file"].unique().to_list())
    if unseen := sorted({p.name for p in chunk_files()} - seen):
        typer.secho(f"{built} was built without {len(unseen)} chunk(s) that exist now: "
                    f"{', '.join(unseen[:5])}. Their references would be resolved away silently. "
                    f"Re-run `vdjdb build` first.", fg=typer.colors.RED, err=True)
        raise typer.Exit(2)
    table = refs_mod.refresh(records, out or refs_mod.TABLE)
    for row in table.group_by("source").len().sort("source").iter_rows():
        typer.echo(f"{row[0]:14} {row[1]:>5,}")
    missing = refs_mod.unresolved(records, table)
    typer.echo(f"{'resolved':14} {table.height:>5,} of "
               f"{records['reference.id'].n_unique():,} distinct reference.id")
    if missing.height:
        typer.secho(f"{missing.height} unresolved, "
                    f"{int(missing['records'].sum()):,} records:", fg=typer.colors.YELLOW, err=True)
        for ref, n in missing.head(10).iter_rows():
            typer.echo(f"    {ref}  ({n:,} records)", err=True)


@app.command(name="release")
def release_cmd(
    tag: str = typer.Option(..., help="v<YYYY>.<MM>.<PATCH>, e.g. v2026.09.1."),
    build_dir: Path = typer.Option(Path("out"), "--build", help="What the build produced."),
    out: Path = typer.Option(Path("out/release"), help="Where the assets go."),
    previous: Path | None = typer.Option(None, help="Previous release zip, for the changelog."),
    previous_lifecycle: Path | None = typer.Option(None, help="Previous release's lifecycle TSV, so "
                                                             "retirements carry the release that "
                                                             "last held them."),
) -> None:
    """Assemble the release: three zips, `manifest.json`, `SHA256SUMS`, `latest-version.txt`.

    Does not publish. Steps 4-6 of the release - verify, publish and the commit-then-check of
    `latest-version.txt` - belong to `release.yml`, which has the repository write access and
    network access they need.
    """
    from .identity.checks import read_tables
    from .identity.lifecycle import advance, compare, present, read, write
    from .release import bundle as b
    from .release import changelog as cl

    version = b.version_of(tag)
    b.prepare_latest(tag)
    typer.echo(f"latest-version.txt line 1 -> {b.legacy_url(tag)}")
    b.stage(build_dir)
    # The lifecycle is written here and not by a build: a curation branch that adds a clonotype and
    # removes it again has retired nothing, so only a release moves these rows (`ROADMAP.md` 10.4).
    lifecycle_path = out / "identity-lifecycle.tsv"
    before, now = read(previous_lifecycle), present(read_tables(build_dir / "tables"))
    if not now.is_empty():
        out.mkdir(parents=True, exist_ok=True)
        write(advance(before, now, release=tag), lifecycle_path)
        typer.echo(compare(before, now, release=tag).as_markdown())
        typer.echo(f"identity-lifecycle.tsv -> {now.height:,} active ids")
    result = b.build(build_dir, out, tag, extras=(lifecycle_path,))
    for entry in result["bundles"]:
        typer.echo(f"{entry['role']:8} {entry['file']:34} {entry['bytes']:>13,} bytes  "
                   f"{len(entry['members']):>2} members")
    if previous:
        notes = cl.render(cl.diff(previous, build_dir / "legacy",
                                  years=Path("summary/reference_years.tsv")), tag=tag)
        (out / "RELEASE_NOTES.md").write_text(notes + "\n")
        typer.echo(f"release notes -> {out / 'RELEASE_NOTES.md'}")
    typer.echo(f"manifest -> {out / 'manifest.json'}   checksums -> {out / 'SHA256SUMS'}")
    typer.echo(f"version {version}")


@app.command()
def changelog(
    previous: Path = typer.Argument(..., help="Previous release zip or legacy directory."),
    current: Path = typer.Argument(Path("out/legacy"), help="This build's legacy projection."),
    tag: str = typer.Option("", help="Tag to title the notes with."),
) -> None:
    """The reference diff between two releases - which studies arrived, left, or were renamed."""
    from .release import changelog as cl

    d = cl.diff(previous, current, years=Path("summary/reference_years.tsv"))
    typer.echo(cl.render(d, tag=tag))


@app.command()
def diff(
    reference: Path = typer.Argument(..., help="Reference release zip or directory."),
    candidate: Path = typer.Argument(..., help="Candidate build directory or zip."),
    rules: Path = typer.Option(Path("rules/expected_diffs.toml"),
                               help="Declared expected differences."),
    report: Path | None = typer.Option(None, help="Write the comparison report here as Markdown."),
    json_out: Path | None = typer.Option(
        None, "--json", help="Also write the report as JSON, without the per-cell detail, so the "
                             "release tests can read it instead of running this again."),
    only: str | None = typer.Option(None, help="Comma-separated members to compare; "
                                               "for deliberately partial candidates."),
) -> None:
    """Compare a candidate build against a released bundle and attribute every difference."""
    from .compare.diff import diff as run_diff
    from .compare.diff import render, summary_json

    result = run_diff(reference, candidate, rules if rules.exists() else None,
                      only=[s.strip() for s in only.split(",")] if only else None)
    text = render(result)
    if report:
        report.parent.mkdir(parents=True, exist_ok=True)
        report.write_text(text)
    if json_out:
        json_out.parent.mkdir(parents=True, exist_ok=True)
        json_out.write_text(summary_json(result) + "\n")
    typer.echo(text, nl=False)
    raise typer.Exit(0 if result.ok else 1)


if __name__ == "__main__":
    app()


identity_app = typer.Typer(
    help="Stable identifiers: check them, resolve one, compare two releases, write the lifecycle.",
    no_args_is_help=True,
)
app.add_typer(identity_app, name="identity")


@identity_app.command("check")
def identity_check(
    tables: Path = typer.Option(Path("out/tables"), help="A built new-format directory."),
    previous: Path | None = typer.Option(None, help="The previous release's lifecycle TSV."),
    previous_tables: Path | None = typer.Option(None, help="The previous release's tables, for the "
                                                          "TCR_hash check."),
) -> None:
    """Assert the identity invariants on a built directory.

    Six of the seven in `ROADMAP.md` section 10.5. The seventh, that a permuted chunk order changes
    no id, is a property of two builds and is a unit test rather than a check on one directory.
    Checks needing history are skipped when `--previous` is absent, because a fork has none.
    """
    from .identity.checks import check, read_tables
    from .identity.lifecycle import read

    built = read_tables(tables)
    if not built:
        typer.secho(f"no tables under {tables}", fg=typer.colors.RED, err=True)
        raise typer.Exit(2)
    prev_chains = read_tables(previous_tables).get("chains") if previous_tables else None
    found = check(built, previous_lifecycle=read(previous), previous_chains=prev_chains)
    for name, frame in sorted(built.items()):
        typer.echo(f"{name:12} {frame.height:>10,} rows")
    if not found:
        typer.secho("every invariant holds", fg=typer.colors.GREEN)
        return
    for finding in found:
        typer.secho(str(finding), fg=typer.colors.RED, err=True)
    raise typer.Exit(1)


@identity_app.command("resolve")
def identity_resolve(
    identifier: str = typer.Argument(..., help="An id, e.g. CT3efa364e7f84aac9."),
    lifecycle: Path = typer.Option(..., help="A lifecycle TSV from a release."),
) -> None:
    """What happened to one id: its level, state, and the releases that carried it."""
    from .identity.levels import level_of
    from .identity.lifecycle import read, resolve

    level = level_of(identifier)
    if level is None:
        typer.secho(f"{identifier}: no VDJdb id has that prefix", fg=typer.colors.RED, err=True)
        raise typer.Exit(2)
    row = resolve(read(lifecycle), identifier)
    if row is None:
        typer.secho(f"{identifier}: a {level} id, never published", fg=typer.colors.YELLOW)
        raise typer.Exit(1)
    for key, value in row.items():
        typer.echo(f"{key:15} {value}")


@identity_app.command("diff")
def identity_diff(
    previous: Path = typer.Argument(..., help="The previous release's lifecycle TSV."),
    tables: Path = typer.Option(Path("out/tables"), help="A built new-format directory."),
    release: str = typer.Option("dev", help="Tag naming the build under comparison."),
) -> None:
    """Ids added, returned and retired since a release. Writes nothing."""
    from .identity.checks import read_tables
    from .identity.lifecycle import compare, present, read

    report = compare(read(previous), present(read_tables(tables)), release=release)
    typer.echo(report.as_markdown())


@identity_app.command("lifecycle")
def identity_lifecycle(
    release: str = typer.Option(..., help="The release tag these ids are first or last seen in."),
    tables: Path = typer.Option(Path("out/tables"), help="A built new-format directory."),
    previous: Path | None = typer.Option(None, help="The previous release's lifecycle TSV."),
    out: Path = typer.Option(Path("out/release/identity-lifecycle.tsv"), help="Where to write it."),
) -> None:
    """Write the lifecycle table for a release.

    Called by the release job and not by a build: a curation branch that adds a clonotype and removes
    it again has not retired anything, so only a release moves these rows.
    """
    from .identity.checks import read_tables
    from .identity.lifecycle import advance, compare, present, read, write

    before = read(previous)
    now = present(read_tables(tables))
    write(advance(before, now, release=release), out)
    typer.echo(compare(before, now, release=release).as_markdown())
    typer.echo(f"wrote {out} ({now.height:,} active ids)")


corpus_app = typer.Typer(
    help="The reference corpus: build it, refresh its PubMed input, search it, condition on it.",
    no_args_is_help=True,
)
app.add_typer(corpus_app, name="corpus")


@corpus_app.command("build")
def corpus_build(
    tables: Path = typer.Option(Path("out/tables"), help="A built new-format directory."),
    out: Path = typer.Option(Path("out/corpus"), help="Where the corpus goes."),
    k: int = typer.Option(3, help="k of the CDR3 and epitope k-mer families."),
    text: Path | None = typer.Option(None, help="corpus/text_terms.tsv; default the committed one."),
) -> None:
    """Build the corpus: documents, vocabulary, postings and the MHC dictionary.

    Recomputed from the tables every time; nothing is read back from a previous corpus (hard rule 9).
    The text family needs `corpus/text_terms.tsv`, which `vdjdb corpus refs` writes; without it the
    corpus builds from the receptor, antigen and MHC families and says so.
    """
    import polars as pl

    from .corpus import build as cb
    from .corpus import pubmed as cp

    need = {n: tables / f"{n}.parquet" for n in ("records", "chains", "restriction")}
    if missing := [str(p) for p in need.values() if not p.exists()]:
        typer.secho(f"missing: {', '.join(missing)}", fg=typer.colors.RED, err=True)
        raise typer.Exit(2)
    frames = {n: pl.read_parquet(p) for n, p in need.items()}
    terms = cp.load_terms(text)
    if terms.is_empty():
        typer.secho("no text_terms.tsv: building without the word family "
                    "(run `vdjdb corpus refs`)", fg=typer.colors.YELLOW)
    records = (pl.read_csv(cp.PUBMED_TABLE, separator="\t", infer_schema=False)
               if cp.PUBMED_TABLE.exists() else None)
    corpus = cb.build(frames["records"], frames["chains"], frames["restriction"],
                      text_terms=terms, pubmed_records=records, k=k)
    for name, path in cb.write(corpus, out).items():
        typer.echo(f"{name:10} {corpus[name].height:>9,} rows  {path}")
    families = (corpus["terms"].group_by("family").len().sort("len", descending=True))
    for row in families.iter_rows():
        typer.echo(f"  {row[0] or '(none)':16} {row[1]:>8,} terms")


@corpus_app.command("refs")
def corpus_refs(
    tables: Path = typer.Option(Path("out/tables"), help="A built new-format directory."),
) -> None:
    """Refresh the committed PubMed input: records and word counts, no running text.

    Hits the network, like `vdjdb refs`, and writes two reviewed inputs that a pull request carries.
    A build never runs this (hard rule 9): the tables it writes are inputs, not results.
    """
    import polars as pl

    from .corpus import pubmed as cp

    path = tables / "records.parquet"
    if not path.exists():
        typer.secho(f"no records table at {path}", fg=typer.colors.RED, err=True)
        raise typer.Exit(2)
    refs = sorted(set(pl.read_parquet(path, columns=["reference.id"])["reference.id"].to_list()))
    records, terms, missing = cp.build_tables(refs)
    for p in cp.write(records, terms):
        typer.echo(f"{p}  {p.stat().st_size:,} bytes")
    typer.echo(f"{records.height:,} PubMed records, {terms.height:,} term counts, "
               f"{terms['term'].n_unique():,} distinct words")
    if missing:
        typer.secho(f"{len(missing)} PMID(s) returned no record: {', '.join(missing[:10])}",
                    fg=typer.colors.YELLOW, err=True)


@corpus_app.command("query")
def corpus_query(
    cdr3: str = typer.Option("", help="Space-joined CDR3 motifs, as the refsearch client sends."),
    epitope: str = typer.Option("", help="Space-joined epitopes."),
    species: str = typer.Option("", help="Space-joined host species; default the client's three."),
    extra: str = typer.Option("", help="search_by_antigen and/or filter_stop_words."),
    corpus: Path = typer.Option(Path("out/corpus"), help="A built corpus."),
    limit: int = typer.Option(10, help="Rows to return; the refsearch client keeps 10."),
) -> None:
    """Rank references against a refsearch-shaped query."""
    from .corpus import build as cb
    from .corpus import query as cq

    built = cb.read(corpus)
    if not built:
        typer.secho(f"no corpus at {corpus}", fg=typer.colors.RED, err=True)
        raise typer.Exit(2)
    terms = cq.refsearch_query(cdr3, epitope, extra_parameters=extra, species_to_search=species)
    typer.echo(f"{len(terms)} query tokens")
    hits = cq.score(built, terms, limit=limit)
    if hits.is_empty():
        typer.secho("no document carries any of those tokens", fg=typer.colors.YELLOW)
        return
    for row in hits.iter_rows(named=True):
        typer.echo(f"{row['score']:.6f}  {row['matched']:>3} matched  {row['reference.id']}")


@corpus_app.command("lift")
def corpus_lift(
    term: str = typer.Argument(..., help="A token, e.g. k:CAS."),
    given: list[str] = typer.Option(None, "--given", help="Condition on this token; repeatable."),
    over: str = typer.Option("occurrences", help="documents or occurrences."),
    corpus: Path = typer.Option(Path("out/corpus"), help="A built corpus."),
) -> None:
    """How much more often a token appears among documents carrying the conditions.

    Every number it rests on is printed with it: a lift alone says nothing about how much evidence is
    under it.
    """
    from .corpus import build as cb
    from .corpus import query as cq

    built = cb.read(corpus)
    if not built:
        typer.secho(f"no corpus at {corpus}", fg=typer.colors.RED, err=True)
        raise typer.Exit(2)
    typer.echo(str(cq.lift(built, term, list(given or []), over=over)))
