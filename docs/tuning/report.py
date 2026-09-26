# 2026-09-26  Turn the sweeps into the committed tables the docs cite, and the data the figures
# read. Nothing is recomputed here; this only reshapes what `sweeps.py` measured (CLAUDE.md 0b).
#
#   uv run python docs/tuning/report.py
#
# Writes, all committed:
#   docs/tuning/scorecard.tsv          every cell, every instrument, one row per configuration
#   docs/tuning/scorecard_plot.tsv     the narrow frame tuning.gp reads
#   docs/tuning/per_epitope_ecdf.tsv   ECDFs behind per_epitope.gp
#   docs/tuning/per_epitope_summary.tsv
#   docs/tuning/pooled.gp              measured constants for the figure's rules
# and out/reports/tuning/tables.md, the markdown the docs paste from -- generated, so no number
# in docs/clustering.md is transcribed by hand.
from __future__ import annotations

from pathlib import Path

import polars as pl

SWEEPS = Path("out/reports/tuning")
DOCS = Path("docs/tuning")
PER_EPITOPE = Path("out/motifs_new/reports/motifs_per_epitope.tsv")
COLS = ["gene", "algo", "config", "retention", "lift", "lift_all", "f1", "obj_precision",
        "recall", "q", "h", "p", "purity", "precision", "epitopes", "perc_med", "ret_med",
        "cids", "raw_clusters", "admissible", "admissible_legacy_bar"]
SHOW = ["algo", "config", "lift", "f1", "q", "purity", "precision", "retention", "epitopes",
        "perc_med", "cids"]
REFS = ["legacy", "shipped-tcrnet", "shipped-tcremp", "trivial"]
FMT = {"retention": "{:.4f}", "lift": "{:.3f}", "f1": "{:.4f}", "q": "{:.4f}",
       "purity": "{:.4f}", "precision": "{:.4f}", "perc_med": "{:.3f}"}
BEST = {"lift": "max", "f1": "max", "q": "max", "purity": "max", "precision": "max",
        "epitopes": "max", "perc_med": "min"}


def config() -> pl.Expr:
    """One human-readable parameter string per row, whatever the sweep's own columns were."""
    return (pl.when(pl.col("algo").is_in(REFS)).then(pl.lit("--"))
            .when(pl.col("algo") == "dbscan").then("coef " + pl.col("coef").cast(pl.Utf8))
            .when(pl.col("algo").str.starts_with("hdbscan"))
            .then("mcs" + pl.col("min_cluster_size").cast(pl.Utf8)
                  + " ms" + pl.col("min_samples").cast(pl.Utf8))
            .otherwise("mcs" + pl.col("min_cluster_size").cast(pl.Utf8)
                       + " M" + pl.col("M").cast(pl.Utf8)))


def md(df: pl.DataFrame, fmt: dict | None = None, bold: dict | None = None) -> str:
    """A markdown table: one value per cell, best along the comparison axis bolded (CLAUDE.md 5)."""
    fmt = fmt or {}
    best = {c: (df[c].drop_nulls().max() if d == "max" else df[c].drop_nulls().min())
            for c, d in (bold or {}).items() if df[c].drop_nulls().len()}
    out = ["| " + " | ".join(df.columns) + " |",
           "|" + "|".join("---:" if df[c].dtype.is_numeric() else "---"
                          for c in df.columns) + "|"]
    for r in df.iter_rows(named=True):
        cells = []
        for c in df.columns:
            v = r[c]
            s = ("--" if v is None else fmt[c].format(v) if c in fmt
                 else f"{v:,}" if isinstance(v, int) else str(v))
            cells.append(f"**{s}**" if c in best and v == best[c] else s)
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


def scorecard() -> pl.DataFrame:
    parts = []
    for p in sorted(SWEEPS.glob("*.csv")):
        df = pl.read_csv(p, infer_schema_length=None)
        if "algo" not in df.columns:          # the HDBSCAN sweep keys on `method`, not `algo`
            df = df.with_columns(("hdbscan-" + pl.col("method")).alias("algo"))
        # polars evaluates every branch of when/then, so the parameter columns a given sweep
        # does not have must exist as nulls before config() is applied.
        for c in ("coef", "min_cluster_size", "min_samples", "M"):
            if c not in df.columns:
                df = df.with_columns(pl.lit(None, dtype=pl.Float64).alias(c))
        df = df.with_columns(config().alias("config"))
        for c in COLS:
            if c not in df.columns:
                df = df.with_columns(pl.lit(None).alias(c))
        parts.append(df.select(COLS))
        print(f"{p.name}: {df.height} rows")
    sc = (pl.concat(parts, how="vertical_relaxed")
          .with_columns(pl.col("admissible", "admissible_legacy_bar").fill_null(False))
          .unique(subset=["gene", "algo", "config"], keep="first", maintain_order=True)
          .sort("gene", "algo", "retention"))
    sc.write_csv(DOCS / "scorecard.tsv", separator="\t")
    (sc.select("gene", "algo", "config", "retention", "lift", "f1", "q", "purity", "precision",
               "epitopes", "perc_med", "cids", pl.col("admissible").cast(pl.Int8))
       .write_csv(DOCS / "scorecard_plot.tsv", separator="\t"))
    return sc


def tables(sc: pl.DataFrame) -> str:
    out = []
    for gene in ("TRA", "TRB"):
        g = sc.filter(pl.col("gene") == gene)
        ref = g.filter(pl.col("algo").is_in(REFS)).select(SHOW)
        fam = (g.filter(~pl.col("algo").is_in(REFS)).sort("f1", descending=True)
                .group_by("algo", maintain_order=True).first().select(SHOW))
        out += [f"\n### {gene} — the reference rows\n",
                md(ref, FMT, BEST),
                f"\n### {gene} — best cell per algorithm, ranked on F1\n",
                md(pl.concat([ref, fam]).sort("f1", descending=True), FMT, BEST)]

    rows = []
    for gene in ("TRA", "TRB"):
        g = sc.filter(pl.col("gene") == gene)
        tr = g.filter(pl.col("algo") == "trivial").row(0, named=True)
        lg = g.filter(pl.col("algo") == "legacy").row(0, named=True)
        for axis in ("q", "purity", "precision", "epitopes", "lift", "f1"):
            f = "{:.0f}" if axis == "epitopes" else "{:.4f}"
            rows.append({"gene": gene, "axis": axis,
                         "do-nothing partition": f.format(tr[axis]),
                         "legacy bar": f.format(lg[axis]),
                         "excludes it?": "yes" if tr[axis] < lg[axis] else "**no**"})
    out += ["\n### The instrument audit — which bars the do-nothing partition clears\n",
            md(pl.DataFrame(rows))]

    for gene in ("TRA", "TRB"):
        h = sc.filter((pl.col("gene") == gene) & pl.col("algo").str.starts_with("hdbscan"))
        tr = sc.filter((pl.col("gene") == gene) & (pl.col("algo") == "trivial")).row(0, named=True)
        out += [f"\n### {gene} HDBSCAN, full grid, by purity "
                f"(do-nothing purity = {tr['purity']:.4f})\n",
                md(h.sort("purity", descending=True).select([*SHOW, "admissible"]), FMT, BEST)]
    return "\n".join(out)


def per_epitope() -> None:
    if not PER_EPITOPE.exists():
        print(f"no {PER_EPITOPE}; run `uv run vdjdb motifs` first")
        return
    df = pl.read_csv(PER_EPITOPE, separator="\t").filter(pl.col("species") == "HomoSapiens")
    parts = []
    for stat in ("retention", "percolation"):
        for (gene, method), g in df.group_by(["gene", "method"], maintain_order=True):
            v = g.select(stat).drop_nulls().sort(stat)
            if v.height:
                parts.append(v.with_row_index("__i").select(
                    pl.lit(gene).alias("gene"), pl.lit(method).alias("method"),
                    pl.lit(stat).alias("stat"), pl.col(stat).alias("value"),
                    ((pl.col("__i") + 1) / v.height).round(5).alias("frac")))
    pl.concat(parts).write_csv(DOCS / "per_epitope_ecdf.tsv", separator="\t")

    summary = (df.group_by("gene", "method").agg(
        pl.len().alias("epitopes"),
        (pl.col("clustered").sum() / pl.col("clonotypes").sum()).round(4).alias("pooled_retention"),
        pl.col("retention").median().round(4).alias("median_epitope_retention"),
        pl.col("percolation").median().round(4).alias("median_percolation"),
        (pl.col("percolation") >= 0.9).sum().alias("epitopes_ge90pct_one_cluster"),
        (pl.col("clusters") > 0).sum().alias("covered")).sort("gene", "method"))
    summary.write_csv(DOCS / "per_epitope_summary.tsv", separator="\t")
    with (DOCS / "pooled.gp").open("w") as fh:
        fh.write("# generated by docs/tuning/report.py; do not edit\n")
        for r in summary.iter_rows(named=True):
            fh.write(f"pooled_{r['gene']}_{r['method']} = {r['pooled_retention']}\n")
            fh.write(f"median_{r['gene']}_{r['method']} = {r['median_epitope_retention']}\n")
    print(summary)


if __name__ == "__main__":
    DOCS.mkdir(parents=True, exist_ok=True)
    sc = scorecard()
    SWEEPS.mkdir(parents=True, exist_ok=True)
    (SWEEPS / "tables.md").write_text(tables(sc))
    per_epitope()
    print(f"\n{sc.height} rows, {sc['algo'].n_unique()} algorithms -> docs/tuning/scorecard.tsv")
    print(f"markdown -> {SWEEPS}/tables.md")
