# 2026-09-26  Shared scoring harness for the clustering bake-off (VDJdb motif stage).
#
# Every candidate algorithm is scored through ONE code path: produce a label vector in the
# tcremp contract (-1 noise, otherwise unique across epitopes), push it through clusters() so the
# legacy length-split and the min_cluster post-filter apply identically, then measure:
#
#   lift       independent-study lift on the NON-DISPLAY denominator (docs/denoising.md 6.1)
#   f1         F1 of `clustered` predicting `replicated` -- the ranking statistic. Lift falls to
#              1.0 mechanically as retention -> 1, so it cannot compare methods that sit at
#              different retentions; F1 has no such bias, and it is what section 11.1 ranks.
#   q,h,p      the homogeneity-parsimony trade-off (clustereval)
#   purity/precision/retention   metrics_lib, vendored
#   epitopes   epitopes with >=1 shipped cid -- the fourth admissibility axis
#   perc_med   median per-epitope largest-cluster share: the percolation diagnostic
#
# Nothing here reads TCRvdb.
from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import polars as pl

from vdjdb.assemble.evidence import support_counts
from vdjdb.motifs import tcremp
from vdjdb.motifs.tcremp import _objective
from vdjdb.validate import motif_bench as mb
from vdjdb.validate import qscore

T = Path("out/tables")
REF = Path("/tmp/vdjdbref/vdjdb-2026-06-03")
MIN_CLUSTER = 5          # the legacy post-filter, fixed so the algorithm is the only variable
SPECIES = "HomoSapiens"


def load() -> dict:
    chains = pl.read_parquet(T / "chains.parquet")
    records = pl.read_parquet(T / "records.parquet")
    legacy = mb.read_members(REF / "cluster_members.txt")
    support = support_counts(records, chains).select(
        "clonotype_id", "antigen.epitope", (pl.col("studies") > 1).alias("replicated"))
    disp = (chains.join(records.filter(pl.col("method.identification").str.to_lowercase()
                                       .str.contains(tcremp.DISPLAY)).select("record_id"),
                        on="record_id").select("clonotype_id").unique()
            .with_columns(pl.lit(True).alias("__d")))
    return {"chains": chains, "records": records, "legacy": legacy,
            "support": support, "disp": disp,
            "cohort": tcremp.cohort(chains, records)}


def context(db: dict, gene: str, *, n_components: int = tcremp.N_COMPONENTS) -> dict:
    """Everything a sweep needs for one chain: the embedding, the bench cohort, legacy's score."""
    grp = db["cohort"].filter((pl.col("species") == SPECIES) & (pl.col("gene") == gene))
    bench = mb.cohort(db["chains"], db["records"], gene=gene)
    leg = db["legacy"].filter((pl.col("gene") == gene) & (pl.col("species") == SPECIES))
    t = time.time()
    X = tcremp.embed_reduced(grp, SPECIES, gene, n_components=n_components)
    secs = time.time() - t
    lab_df = (grp.select("clonotype_id", "antigen.epitope", "junction_aa", "v_call", "j_call")
              .join(db["support"], on=["clonotype_id", "antigen.epitope"], how="left")
              .with_columns(pl.col("replicated").fill_null(False)))
    keep = lab_df.join(db["disp"], on="clonotype_id", how="left")["__d"].is_null().to_numpy()
    ctx = {"gene": gene, "grp": grp, "bench": bench, "X": X,
           "epi": grp["antigen.epitope"].to_numpy(), "lab_df": lab_df, "keep": keep,
           "is_rep": lab_df["replicated"].to_numpy(),
           "nn": tcremp.mean_nn_distance(X), "embed_secs": secs,
           "bench_epitopes": bench["antigen.epitope"].n_unique()}
    ctx["legacy"] = _measure(ctx, leg)
    return ctx


def _measure(ctx: dict, mm: pl.DataFrame) -> dict:
    """Every instrument, for any frame with ``(antigen.epitope, cdr3aa, v.segm, j.segm)``.

    Legacy's shipped file and a candidate's output both go through this function, so the two cannot
    be scored by different code. Scoring them separately made section 35's first TRA reading wrong
    by 1.6x.
    """
    s = mb.score(mb.assign(ctx["bench"], mm))
    q = qscore.score(ctx["bench"], mm)
    pe = mb.per_epitope(ctx["bench"], mm)
    covered = pe.filter(pl.col("clusters") > 0)
    member = (ctx["lab_df"].join(
        mm.select("antigen.epitope", pl.col("cdr3aa").alias("junction_aa"),
                  pl.col("v.segm").alias("v_call"), pl.col("j.segm").alias("j_call"),
                  pl.lit(True).alias("__c")).unique(),
        on=["antigen.epitope", "junction_aa", "v_call", "j_call"], how="left")
        ["__c"].fill_null(False).to_numpy())
    lab = np.where(member, 0, -1)
    o = _objective(lab[ctx["keep"]], ctx["is_rep"][ctx["keep"]])
    oa = _objective(lab, ctx["is_rep"])
    return {"lift": o["lift"], "lift_all": oa["lift"], "f1": o["f1"],
            "obj_precision": o["precision"], "recall": o["recall"],
            "q": q["q"], "h": q["h"], "p": q["p"], "cids": q["clusters"],
            "purity": s["purity"], "precision": s["precision"], "retention": s["retention"],
            "epitopes": covered.height,
            "perc_med": float(covered["percolation"].median() or 0.0),
            "ret_med": float(pe["retention"].median() or 0.0)}


def score_labels(ctx: dict, labels: np.ndarray, *, min_cluster: int = MIN_CLUSTER) -> dict | None:
    """One row of the scorecard for one label vector. ``None`` when nothing survives the filter."""
    m = tcremp.clusters(ctx["grp"], ctx["X"], labels=labels, min_cluster=min_cluster)
    if not m.height:
        return None
    mm = m.select("species", "gene", "antigen.epitope", pl.col("junction_aa").alias("cdr3aa"),
                  pl.col("v_call").alias("v.segm"), pl.col("j_call").alias("j.segm"), "cid")
    row = _measure(ctx, mm)
    row["raw_clusters"] = len(np.unique(labels[labels >= 0]))
    return row


def admissible(ctx: dict, row: dict, *, purity_floor: float | None = None) -> bool:
    """The four-axis rule of docs/denoising.md 7.1.

    ``purity_floor`` replaces the legacy-relative purity/precision bar with an absolute one when
    given: the 0.94 floor of docs/clustering.md section 8.1. Q and epitope coverage stay
    legacy-relative either way.
    """
    leg = ctx["legacy"]
    pbar = leg["purity"] if purity_floor is None else purity_floor
    xbar = leg["precision"] if purity_floor is None else purity_floor
    return bool(row["q"] >= leg["q"] and row["purity"] >= pbar
                and row["precision"] >= xbar and row["epitopes"] >= leg["epitopes"])


FMT = ("{tag:<26} lift {lift:6.3f} F1 {f1:.4f} | Q {q:.4f} (h {h:.3f} p {p:.3f}) | "
       "pur {purity:.4f} prec {precision:.4f} ret {retention:.4f} | "
       "ep {epitopes:3d} perc {perc_med:.3f} cids {cids:,}")


def line(tag: str, row: dict, flag: str = "") -> str:
    return FMT.format(tag=tag, **row) + (f" {flag}" if flag else "")
