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
              .agg(pl.len().cast(pl.Int32).alias("duplicate_count"))
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
