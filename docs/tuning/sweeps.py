# 2026-09-26  Every candidate partition, scored through one harness.
#
#   uv run python docs/tuning/sweeps.py dbscan|hdbscan|lumbermark|hybrid|hybrid-len|all
#
# Writes one CSV per sweep into out/reports/tuning/; `report.py` merges them into the committed
# scorecard. Each sweep produces a label vector in the tcremp contract and hands it to
# clusterlab.score_labels, so no two methods can be scored by different code -- the mistake that
# made section 35's first TRA reading wrong by 1.6x.
#
# Nothing here reads TCRvdb. The objective is independent-study replication (ROADMAP section 11.1).
from __future__ import annotations

import itertools
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).parent))
import clusterlab as cl
from lumberlab import lumbermark_labels

warnings.filterwarnings("ignore", message="The number of clusters detected")

OUT = Path("out/reports/tuning")

#: The relaxed purity floor, per chain. ``None`` keeps the legacy-relative bar. It must sit ABOVE
#: the do-nothing partition's purity or it is not a bar at all -- measured at 0.9204 (TRA) and
#: 0.9343 (TRB), so 0.93 would admit a partition that clusters nothing. See docs/clustering.md 8.
BAR: dict[str, float | None] = {"TRA": None, "TRB": 0.94}

GRID = {
    "dbscan": {"coef": (0.4, 0.6, 0.9, 1.15, 1.3, 1.55, 1.8, 2.1, 2.5)},
    "hdbscan": {"mcs": (3, 5, 8, 10, 15, 20, 30, 50), "meth": ("eom", "leaf"),
                "ms": (2, 3, 5, 10)},
    "lumbermark": {"mcs": (3, 5, 10, 20, 50), "M": (1, 5, 10)},
    "hybrid": {"mcs": (3, 5, 10, 20, 50), "M": (1, 5, 10)},
    "hybrid-len": {"mcs": (3, 5, 8), "M": (1, 5)},
}


def _gate(ctx, scored, *, recruited: bool) -> np.ndarray:
    """The TCRNET noise model as a boolean mask over the cohort.

    ``recruited=False`` is the enrichment test alone; ``True`` adds each enriched clonotype's
    within-scope neighbours, which is the vertex set the shipped TCRNET actually partitions.
    """
    enr = (ctx["grp"].select("antigen.epitope", "junction_aa", "v_call", "j_call")
           .join(scored.filter((pl.col("gene") == ctx["gene"]) & pl.col("enriched"))
                 .select("antigen.epitope", "junction_aa", "v_call", "j_call")
                 .unique().with_columns(pl.lit(True).alias("__g")),
                 on=["antigen.epitope", "junction_aa", "v_call", "j_call"], how="left")
           ["__g"].fill_null(False).to_numpy())
    if not recruited:
        return enr
    from vdjdb.motifs.cluster import _recruited
    from vdjdb.motifs.tcrnet import TUNED

    scope = TUNED[ctx["gene"]]["scope"]
    seqs = ctx["grp"]["junction_aa"].to_list()
    out = np.zeros(len(seqs), dtype=bool)
    for ep in np.unique(ctx["epi"]):
        idx = np.flatnonzero(ctx["epi"] == ep)
        hits = [k for k, g in enumerate(enr[idx]) if g]
        if hits:
            out[idx[_recruited([seqs[i] for i in idx], hits, scope)]] = True
    return out


def _strat_labels(X, epitopes, lengths, mcs, m_smooth, gate):
    """Lumbermark inside each ``(epitope, CDR3 length)`` stratum -- the unit emit actually ships."""
    from lumbermark import Lumbermark

    out = np.full(len(X), -1, dtype=np.int64)
    offset = 0
    for ep, ln in sorted({(e, int(n)) for e, n in zip(epitopes, lengths, strict=True)}):
        idx = np.flatnonzero((epitopes == ep) & (lengths == ln) & gate)
        if len(idx) < max(2 * mcs, m_smooth + 2):
            continue
        lab = Lumbermark(n_clusters=len(idx) - 1, M=m_smooth, min_cluster_size=mcs,
                         min_cluster_factor=0.0).fit(
            np.ascontiguousarray(X[idx], dtype=np.float64)).labels_
        hit = lab >= 0
        out[idx[hit]] = lab[hit] + offset
        offset += int(lab.max()) + 1 if hit.any() else 0
    return out


def _enrichment(db):
    """One TCRNET pass for both chains. It builds a million-clonotype index per chain."""
    from vdjdb.motifs import tcrnet

    return tcrnet.enriched_clonotypes(
        db["chains"].join(db["records"].select("record_id", "species"), on="record_id")
        .filter(pl.col("species") == cl.SPECIES).drop("species"), db["records"])


def run(which: str) -> None:
    db = cl.load()
    scored = _enrichment(db) if which.startswith("hybrid") else None
    rows: list[dict] = []
    for gene in ("TRA", "TRB"):
        ctx = cl.context(db, gene)
        print(cl.line(f"[{gene}] legacy", ctx["legacy"]), flush=True)
        if which == "dbscan":                 # the reference rows ride along with the shipped path
            rows.append({"gene": gene, "algo": "legacy", "coef": 0.0, **ctx["legacy"],
                         "admissible": True})
            for name, path in (("tcrnet", "cluster_members.txt"),
                               ("tcremp", "cluster_members_tcremp.txt")):
                p = Path("out/motifs_new") / path
                if not p.exists():
                    continue
                mm = cl.mb.read_members(p).filter((pl.col("gene") == gene)
                                                  & (pl.col("species") == cl.SPECIES))
                r = cl._measure(ctx, mm)
                rows.append({"gene": gene, "algo": f"shipped-{name}", "coef": 0.0, **r,
                             "admissible": cl.admissible(ctx, r, purity_floor=BAR[gene])})
                print(cl.line(f"[{gene}] shipped {name}", r), flush=True)
            # One cluster per epitope, nothing excluded: the instrument's blind spot, measured
            # rather than argued (docs/clustering.md section 8).
            triv = np.unique(ctx["epi"], return_inverse=True)[1].astype(np.int64)
            r = cl.score_labels(ctx, triv, min_cluster=1)
            rows.append({"gene": gene, "algo": "trivial", "coef": 0.0, **r,
                         "admissible": cl.admissible(ctx, r, purity_floor=BAR[gene])})
            print(cl.line(f"[{gene}] trivial one-per-epitope", r), flush=True)

        gate = None
        if which.startswith("hybrid"):
            gate = _gate(ctx, scored, recruited=(which == "hybrid-recruited"))
            print(f"  TCRNET gate keeps {gate.sum():,} of {len(gate):,} ({gate.mean():.1%})",
                  flush=True)
        lengths = ctx["grp"]["junction_aa"].str.len_chars().to_numpy()
        g = GRID[which if which != "hybrid-recruited" else "hybrid"]

        cells = (list(itertools.product(g["coef"])) if which == "dbscan" else
                 list(itertools.product(g["mcs"], g["meth"], g["ms"])) if which == "hdbscan" else
                 list(itertools.product(g["mcs"], g["M"])))
        for cell in cells:
            t = time.time()
            if which == "dbscan":
                (coef,) = cell
                eps = cl.tcremp.chain_eps(ctx["X"], coef)
                labels = cl.tcremp.cluster_labels(ctx["X"], ctx["epi"], eps)
                tag, extra = f"dbscan coef{coef}", {"algo": "dbscan", "coef": coef}
            elif which == "hdbscan":
                mcs, meth, ms = cell
                if ms > mcs:
                    continue
                labels = cl.tcremp.hdbscan_labels(ctx["X"], ctx["epi"], min_cluster_size=mcs,
                                                  min_samples=ms, method=meth)
                tag = f"hdbscan-{meth} mcs{mcs} ms{ms}"
                extra = {"algo": f"hdbscan-{meth}", "min_cluster_size": mcs, "min_samples": ms,
                         "method": meth}
            else:
                mcs, m = cell
                labels = (_strat_labels(ctx["X"], ctx["epi"], lengths, mcs, m, gate)
                          if which == "hybrid-len" else
                          lumbermark_labels(ctx["X"], ctx["epi"], mcs, m, gate=gate))
                tag = f"{which} mcs{mcs} M{m}"
                extra = {"algo": which, "min_cluster_size": mcs, "M": m}
            row = cl.score_labels(ctx, labels)
            if row is None:
                print(f"[{gene}] {tag} no clusters", flush=True)
                continue
            adm = cl.admissible(ctx, row, purity_floor=BAR[gene])
            rows.append({"gene": gene, **extra, **row, "admissible": adm,
                         "admissible_legacy_bar": cl.admissible(ctx, row),
                         "secs": round(time.time() - t, 1)})
            print(cl.line(f"[{gene}] {tag}", row, "OK" if adm else ""), flush=True)
            OUT.mkdir(parents=True, exist_ok=True)
            pl.DataFrame(rows).write_csv(OUT / f"{which}.csv")
    print(f"wrote {OUT}/{which}.csv", flush=True)


if __name__ == "__main__":
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    for name in (["dbscan", "hdbscan", "lumbermark", "hybrid", "hybrid-recruited", "hybrid-len"]
                 if what == "all" else [what]):
        run(name)
