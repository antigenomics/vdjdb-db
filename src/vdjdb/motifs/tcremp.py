"""TCREMP: cluster CDR3s by where they land in a prototype-distance embedding.

The second motif method, run beside :mod:`vdjdb.motifs.tcrnet` rather than instead of it. TCRNET
asks whether a clonotype has more neighbours than a background repertoire explains, which only sees
sequences that are literally one substitution apart. TCREMP embeds each clonotype as its distance to
a fixed prototype panel and clusters in that space, so it groups receptors that are similar in a way
an edit distance does not express. Measured on the benchmark cohort, that is worth roughly 2.4x the
retention at comparable purity (ROADMAP section 8.2).

Four things here are measured decisions, not defaults.

**Chunking is mandatory, not an optimisation.** A single-shot ``embed`` -> ``StandardScaler`` ->
``PCA`` on human TRB peaks at **9.77 GB**, and 12.41 GB if the scaler is fitted on the full matrix;
that OOMs a 16 GB runner alongside the rest of the build. The scaler and the PCA are fitted on a
seeded :data:`FIT_SAMPLE` subsample and every chunk is reduced to 50 dimensions **as it is embedded**,
so the wide ``(n, 3 * n_prototypes)`` matrix never exists whole: peak is one chunk, ~480 MB at 2,000
prototypes, against 22 MB for the entire reduced result. Chunks run sequentially, one internally
threaded ``embed()`` each -- never a worker pool (CLAUDE.md hard rule 3, ROADMAP section 8.8).

**Kneedle is dead at production scale and is not in the operating path.** On pooled human TRB,
112,983 clonotypes, it returns knee 1 -- index 1 of 112,983, fraction 0.000 -- so the `floor_frac`
guard fires every time and the rule that actually runs is ``eps = coef * mean(1st-NN distance)``.
:func:`knee_eps_debug` keeps it as a cross-check; the method is not Kneedle-based and should not be
described as such (ROADMAP section 8.3).

**Pooled geometry, per-epitope scope.** ``eps`` is estimated once per chain from the pooled
k-distance curve and then DBSCAN runs per epitope. Re-estimating it on a per-epitope n of 30-300
is exactly the degenerate regime the knee fails in. ⚠ Purity is then 1.000 *by construction*,
because a cid cannot span epitopes -- that means cross-epitope contamination is impossible, not that
none exists, and the pooled clustering is the honest measurement (ROADMAP section 8.4).

**``coef`` is fitted on independent-study support, never on TCRvdb.** See :func:`fit_coef`.
"""
from __future__ import annotations

import numpy as np
import polars as pl

from ..config import SEED

#: Clonotypes embedded per call. One internally-threaded ``embed()`` each, run in sequence.
CHUNK = 20_000

#: Clonotypes the ``StandardScaler`` and ``PCA`` are fitted on, drawn with :data:`SEED`. Re-fitted
#: every build -- storing a fit would risk shipping one that no longer matches the data, and it
#: buys minutes (CLAUDE.md hard rule 9).
FIT_SAMPLE = 25_000

#: PCA dimensions. The benchmark's front-end, applied to every method so the representation is the
#: only variable.
N_COMPONENTS = 50

#: ``DBSCAN(min_samples=...)``. Two: the TCREMP paper's Table-1 value, and the smallest number that
#: can be a cluster at all.
MIN_SAMPLES = 2

#: Epitopes with fewer records than this are not clustered -- the benchmark cohort's floor, and
#: below it there is nothing for a density method to find.
MIN_RECORDS = 30

#: The species VDJdb infers motifs for, and their mirpy names. Human and mouse only.
SPECIES: dict[str, str] = {"HomoSapiens": "human", "MusMusculus": "mouse"}

#: ``method.identification`` substring marking a **display selection**. Those records are not
#: independent natural observations -- a library selected against one pMHC yields thousands of
#: receptors one substitution apart by construction -- and a display paper is a single
#: ``reference.id``, so every one of its clonotypes contributes zero independently-replicated pairs
#: while filling a quarter of the human TRB denominator. Measured: 29,688 of 192,753 records
#: (15.4 %), all SLLMWITQV / PMID:40498839. Held out of :func:`fit_coef` for that reason; whether
#: they should also be held out of the clustering is the author's call (ROADMAP section 30.6).
DISPLAY = "display"


def cohort(chains: pl.DataFrame, records: pl.DataFrame, *,
           min_records: int = MIN_RECORDS) -> pl.DataFrame:
    """The clustering cohort: one row per distinct clonotype-epitope pair, human and mouse.

    ``duplicate_count`` is how many records report that clonotype against that epitope. Epitopes
    below ``min_records`` **records** -- not clonotypes -- are dropped, which is the benchmark
    cohort's definition and what its per-chain counts are quoted over.
    """
    df = (chains
          .join(records.select("record_id", "species", "antigen.epitope"), on="record_id")
          .filter(pl.col("species").is_in(pl.Series(list(SPECIES)).implode())
                  & (pl.col("cdr3") != "") & pl.col("gene").is_in(pl.Series(["TRA", "TRB"]).implode())))
    big = (df.group_by(["species", "gene", "antigen.epitope"]).len()
             .filter(pl.col("len") >= min_records)
             .drop("len"))
    return (df.join(big, on=["species", "gene", "antigen.epitope"])
              .group_by(["species", "gene", "antigen.epitope",
                         pl.col("cdr3").alias("junction_aa"),
                         pl.col("v.segm").alias("v_call"), pl.col("j.segm").alias("j_call")],
                        maintain_order=True)
              .agg(pl.len().cast(pl.Int32).alias("duplicate_count"),
                   pl.col("clonotype_id").first())
              .sort("species", "gene", "antigen.epitope", "junction_aa", "v_call", "j_call"))


def embed_reduced(clonotypes: pl.DataFrame, species: str, locus: str, *,
                  chunk: int = CHUNK, fit_sample: int = FIT_SAMPLE,
                  n_components: int = N_COMPONENTS, seed: int = SEED) -> np.ndarray:
    """Embed and reduce in one pass. ``(n_clonotypes, n_components)`` float64.

    The scaler and the PCA are fitted on one seeded subsample, then every chunk is embedded and
    immediately projected, so the wide matrix is never held. Deterministic in ``seed`` and in the
    input order (CLAUDE.md hard rule 7).
    """
    from mir.embedding import TCREmp
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    cols = clonotypes.select("junction_aa", "v_call", "j_call")
    model = TCREmp.from_defaults(SPECIES[species], locus)

    n = cols.height
    rng = np.random.default_rng(seed)
    idx = np.sort(rng.choice(n, size=min(fit_sample, n), replace=False))
    fit = model.embed(cols[idx])
    scaler = StandardScaler().fit(fit)
    pca = PCA(n_components=min(n_components, fit.shape[1], len(idx) - 1),
              random_state=seed).fit(scaler.transform(fit))
    del fit

    out = np.empty((n, pca.n_components_), dtype=np.float64)
    for start in range(0, n, chunk):
        block = cols[start:start + chunk]
        out[start:start + block.height] = pca.transform(scaler.transform(model.embed(block)))
    return out


def mean_nn_distance(X: np.ndarray) -> float:
    """Mean distance to the 1st non-self neighbour, over all of ``X``.

    The reference implementation asks for 4 neighbours, sorts each column and reads column 1.
    Sorting a column does not change its mean, and columns 2 and 3 are never used, so 2 neighbours
    give the identical number for half the work.
    """
    from sklearn.neighbors import NearestNeighbors

    d, _ = NearestNeighbors(n_neighbors=2).fit(X).kneighbors(X)
    return float(d[:, 1].mean())


def chain_eps(X: np.ndarray, coef: float) -> float:
    """The chain-global DBSCAN radius: ``coef * mean(1st-NN distance)``, on the pooled geometry."""
    return round(mean_nn_distance(X) * coef, 3)


def knee_eps_debug(X: np.ndarray, coef: float, k: int = 4, floor_frac: float = 0.40) -> dict:
    """Kneedle, as a diagnostic only. Reports the knee and whether the floor guard fired.

    Kept so the claim in ROADMAP section 8.3 stays checkable rather than remembered. Nothing in the
    operating path calls this.
    """
    from kneed import KneeLocator
    from sklearn.neighbors import NearestNeighbors

    d, _ = NearestNeighbors(n_neighbors=k).fit(X).kneighbors(X)
    dist = np.sort(d, axis=0)[:, 1]
    n = len(dist)
    knee = KneeLocator(range(1, n + 1), dist, S=1.0, curve="concave",
                       interp_method="polynomial", polynomial_degree=10, online=True,
                       direction="increasing").knee
    fired = knee is None or int(knee) < floor_frac * n
    return {"knee": None if knee is None else int(knee), "n": n,
            "knee_frac": None if knee is None else int(knee) / n,
            "floor_guard_fired": bool(fired),
            "eps_knee": None if knee is None else round(float(dist[min(int(knee), n - 1)]) * coef, 3),
            "eps_mean": round(float(dist.mean()) * coef, 3)}


def cluster_labels(X: np.ndarray, epitopes: np.ndarray, eps: float, *,
                   min_samples: int = MIN_SAMPLES) -> np.ndarray:
    """DBSCAN per epitope at one shared ``eps``. ``-1`` is noise, as DBSCAN means it.

    Labels are made unique across epitopes by offsetting, so a label identifies a cluster and not
    just a cluster-within-an-epitope.
    """
    from sklearn.cluster import DBSCAN

    out = np.full(len(X), -1, dtype=np.int64)
    offset = 0
    for ep in np.unique(epitopes):               # np.unique sorts: the offset must not vary
        m = epitopes == ep
        lab = DBSCAN(eps=eps, min_samples=min_samples).fit_predict(X[m])
        hit = lab >= 0
        out[np.flatnonzero(m)[hit]] = lab[hit] + offset
        offset += int(lab.max()) + 1 if hit.any() else 0
    return out


# ---------------------------------------------------------------------------------------------
# Tuning
# ---------------------------------------------------------------------------------------------

#: ``coef`` values swept by :func:`fit_coef`. It has to reach well below the published 0.75: the
#: first sweep stopped at 0.5 and every chain picked the grid edge, which is not a fit.
COEF_GRID: tuple[float, ...] = (0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.8, 2.1)

#: Fitted per chain on the section 11.1 independent-study objective, human, phage-display epitopes
#: held out, 2026-09-25. Both chains take an **interior** optimum at 0.4:
#:
#: ===== ==== ====== ====== ==========
#: chain n    F1     lift   clustered
#: ===== ==== ====== ====== ==========
#: TRB   86k  0.2135 5.31x  4,659
#: TRA   58k  0.1723 4.24x  3,421
#: ===== ==== ====== ====== ==========
#:
#: ⚠ **The published 0.75 does not transfer** -- it was calibrated on standalone `tcremp` with
#: ~3,000 OLGA prototypes and Smith-Waterman, a different metric space (ROADMAP section 8.2). At
#: 0.75 human TRB scores F1 0.19 / lift 3.6x against 0.21 / 5.3x at 0.4.
#: ⚠ **Holding the display block out changes the answer**, it does not merely tidy it: with
#: PMID:40498839 in, human TRB fits coef **1.3** at lift **1.44x**; with it out, **0.4** at
#: **5.31x**. 29,698 display clonotypes with zero independent replication were paying for a wider
#: radius (ROADMAP section 30.6).
COEF: dict[str, float] = {"TRA": 0.4, "TRB": 0.4}


def replicated(records: pl.DataFrame, chains: pl.DataFrame) -> pl.DataFrame:
    """Clonotype-epitope pairs more than one publication reports. The tuning label.

    Reuses :func:`vdjdb.assemble.evidence.support_counts`, so the objective and the shipped
    ``independent_study`` evidence rows are the same computation and cannot drift (section 11.1).
    """
    from ..assemble.evidence import support_counts

    return (support_counts(records, chains)
            .select("clonotype_id", "antigen.epitope",
                    (pl.col("studies") > 1).alias("replicated")))


def _objective(labels: np.ndarray, is_replicated: np.ndarray) -> dict:
    """Score "clustered" as a prediction of "independently replicated". F1, with its parts.

    A clustering that finds real convergent selection should preferentially recover exactly the
    clonotypes a second laboratory saw. Clustering everything wins recall and loses precision;
    clustering nothing scores zero. F1 peaks where the clustering is being selective about the
    right thing -- which is the whole of the objective, and the only number the sweep ranks on.
    """
    clustered = labels >= 0
    tp = int((clustered & is_replicated).sum())
    fp = int((clustered & ~is_replicated).sum())
    fn = int((~clustered & is_replicated).sum())
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    base = float(is_replicated.mean())
    return {"f1": f1, "precision": precision, "recall": recall,
            # How much more often a clustered clonotype is independently replicated than one drawn
            # at random. F1 ranks the grid; lift says whether there is any signal to rank -- an F1
            # of 0.065 means nothing on its own when the base rate is 2.4 %.
            "lift": precision / base if base else 0.0, "base_rate": base,
            "clustered": int(clustered.sum()), "replicated": int(is_replicated.sum()),
            "n_clusters": len(np.unique(labels[clustered])) if clustered.any() else 0}


def fit_coef(cohort_chain: pl.DataFrame, X: np.ndarray, is_replicated: np.ndarray, *,
             grid: tuple[float, ...] = COEF_GRID,
             min_samples: int = MIN_SAMPLES) -> pl.DataFrame:
    """Sweep ``coef`` for one chain and score each value on :func:`_objective`.

    ⚠ **Never fitted against TCRvdb.** That set is held out, touched once, and only in aggregate
    (ROADMAP section 11.2). Tuning on it would destroy the only independent read this project has.

    The mean 1st-NN distance is computed once -- it does not depend on ``coef`` -- so the sweep
    costs one DBSCAN pass per grid point, not one nearest-neighbour search per point.
    """
    epitopes = cohort_chain["antigen.epitope"].to_numpy()
    base = mean_nn_distance(X)
    rows = []
    for coef in grid:
        labels = cluster_labels(X, epitopes, round(base * coef, 3), min_samples=min_samples)
        rows.append({"coef": coef, "eps": round(base * coef, 3),
                     **_objective(labels, is_replicated)})
    return pl.DataFrame(rows).sort("f1", descending=True)


# ---------------------------------------------------------------------------------------------
# Clustering
# ---------------------------------------------------------------------------------------------

def clusters(cohort_chain: pl.DataFrame, X: np.ndarray, eps: float, *,
             min_samples: int = MIN_SAMPLES, min_cluster: int = 5) -> pl.DataFrame:
    """Cluster one chain and label it in the shape :mod:`vdjdb.motifs.emit` expects.

    ``cid`` carries a ``L<len>`` suffix -- **one legacy cid per (cluster, CDR3 length)**. A DBSCAN
    cluster in embedding space may span lengths, where a PWM may not, and ``vdjdb-web`` splits every
    cid by ``len`` before building a cluster anyway; without the suffix two display clusters would
    share one ``clusterId`` and the reader's ``strict = true`` path would break. Measured on the
    REDCEA production files, this is the common case and not an edge case: **809 of 847 TRA and
    1,002 of 1,082 TRB cids span more than one length** (ROADMAP section 8.6).

    ``x``/``y`` are the first two principal components -- the embedding's own layout, so no separate
    force-directed pass is needed and the picture means something.
    """
    from .cluster import _INITIAL, _repr_allele

    epitopes = cohort_chain["antigen.epitope"].to_numpy()
    labels = cluster_labels(X, epitopes, eps, min_samples=min_samples)
    g = (cohort_chain
         .with_columns(pl.Series("__label", labels),
                       pl.Series("x", X[:, 0]), pl.Series("y", X[:, 1]),
                       pl.col("junction_aa").str.len_chars().alias("__len"))
         .filter(pl.col("__label") >= 0))
    if not g.height:
        return g.drop("__label", "__len")

    # The stratum, not the DBSCAN cluster, is what gets a cid -- so `csz` is unambiguous and
    # `freq` in the PWM is a distribution.
    sizes = g.group_by(["species", "gene", "antigen.epitope", "__label", "__len"]).agg(
        pl.len().alias("csz"), pl.col("junction_aa").min().alias("__first"))
    out = []
    for (species, gene, epitope), grp in sizes.group_by(
            ["species", "gene", "antigen.epitope"], maintain_order=True):
        # Size descending then the smallest member: a number that follows the content, so a cluster
        # whose membership is unchanged keeps its id without anything being stored (hard rule 9).
        order = (grp.filter(pl.col("csz") >= min_cluster)
                    .sort(["csz", "__first"], descending=[True, False])
                    .with_row_index("__n", offset=1))
        if order.height:
            prefix = f"{_INITIAL[species]}.{gene[-1]}.{epitope}"
            out.append(order.with_columns(
                (pl.lit(prefix) + "." + pl.col("__n").cast(pl.Utf8)
                 + "L" + pl.col("__len").cast(pl.Utf8)).alias("cid")).drop("__n", "__first"))
    if not out:
        return g.head(0).drop("__label", "__len")

    g = g.join(pl.concat(out, how="vertical"),
               on=["species", "gene", "antigen.epitope", "__label", "__len"], how="inner")
    g = g.join(g.group_by("cid").agg(
        pl.col("v_call").map_batches(_repr_allele, returns_scalar=True).alias("v.segm.repr"),
        pl.col("j_call").map_batches(_repr_allele, returns_scalar=True).alias("j.segm.repr"),
    ), on="cid")
    return (g.drop("__label", "__len")
             .sort("species", "gene", "antigen.epitope", "cid", "junction_aa", "v_call", "j_call"))


def display_epitopes(records: pl.DataFrame) -> pl.Series:
    """Epitopes whose records come from a display selection. See :data:`DISPLAY`."""
    return (records.filter(pl.col("method.identification").str.contains(f"(?i){DISPLAY}"))
                   ["antigen.epitope"].unique().sort())
