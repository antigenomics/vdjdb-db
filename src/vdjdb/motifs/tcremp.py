"""TCREMP: cluster CDR3s by where they land in a prototype-distance embedding.

The second motif method, run beside :mod:`vdjdb.motifs.tcrnet` rather than instead of it. TCRNET
asks whether a clonotype has more neighbours than a background repertoire explains, which only sees
sequences that are literally one substitution apart. TCREMP embeds each clonotype as its distance to
a fixed prototype panel and clusters in that space, so it groups receptors that are similar in a way
an edit distance does not express. Measured on the benchmark cohort, that gives roughly 2.4x the
retention at comparable purity (ROADMAP section 8.2).

Four choices here are measured, not defaults.

Chunking is required, not an optimisation. A single-shot ``embed`` -> ``StandardScaler`` ->
``PCA`` on human TRB peaks at 9.77 GB, and 12.41 GB if the scaler is fitted on the full matrix;
that OOMs a 16 GB runner alongside the rest of the build. The scaler and the PCA are fitted on an
identity-priority :data:`FIT_SAMPLE` sample of unique receptor inputs; every chunk is reduced
to 50 dimensions as it is embedded,
so the wide ``(n, 3 * n_prototypes)`` matrix is never materialised in full. At 2,000
prototypes a 20,000-row raw float32 block is 480 MB; its float64 reduction matrix is
960 MB. Scaling and PCA use these owned buffers in place, released before the next
block. Chunks run sequentially, one
internally threaded ``embed()`` each -- never a worker pool (CLAUDE.md hard rule 3, ROADMAP
section 8.8).

Kneedle is degenerate at production scale and is not in the operating path. On pooled human TRB,
112,983 clonotypes, it returns knee 1 -- index 1 of 112,983, fraction 0.000 -- so the `floor_frac`
guard fires every time and the rule that runs is ``eps = coef * mean(1st-NN distance)``.
:func:`knee_eps_debug` keeps it as a cross-check; the method is not Kneedle-based and should not be
described as such (ROADMAP section 8.3).

Pooled geometry, per-epitope scope. ``eps`` is estimated once per chain from the pooled
k-distance curve and then DBSCAN runs per epitope. Re-estimating it on a per-epitope n of 30-300
is the degenerate regime the knee fails in. ⚠ Purity is then 1.000 by construction, because a cid
cannot span epitopes -- that means cross-epitope contamination is impossible, not that none exists,
and the pooled clustering is what is measured (ROADMAP section 8.4).

``coef`` is fitted on independent-study support, never on TCRvdb. See :func:`fit_coef`.
"""
from __future__ import annotations

import hashlib
import json
import sys
from contextlib import contextmanager
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import polars as pl

from ..config import SEED

if TYPE_CHECKING:
    from mir.embedding import TCREmp
    from sklearn.decomposition import PCA
    from sklearn.neighbors import NearestNeighbors
    from sklearn.preprocessing import StandardScaler

#: Clonotypes embedded per call. One internally-threaded ``embed()`` each, run in sequence.
CHUNK = 20_000

#: Unique receptors the ``StandardScaler`` and ``PCA`` are fitted on, ranked with :data:`SEED`. Re-fitted
#: every build: a stored fit could be one that no longer matches the data, and refitting costs
#: minutes (CLAUDE.md hard rule 9).
FIT_SAMPLE = 25_000

#: The embedding input identity, independent of peptide labels and observation counts.
INPUT_KEY = ("junction_aa", "v_call", "j_call")

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

#: ``method.identification`` substring marking a display selection. Those records are not
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
    below ``min_records`` records -- not clonotypes -- are dropped, which is the benchmark cohort's
    definition and what its per-chain counts are quoted over.
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


def receptor_inputs(clonotypes: pl.DataFrame) -> tuple[pl.DataFrame, np.ndarray]:
    """Sorted unique receptor inputs and their indices for restoring the labelled rows."""
    keys = clonotypes.select(INPUT_KEY).with_columns(pl.col("v_call", "j_call").fill_null(""))
    inputs = keys.unique().sort(INPUT_KEY)
    inverse = keys.join(inputs.with_row_index("__input"), on=INPUT_KEY, how="left",
                        validate="m:1", maintain_order="left")["__input"].to_numpy()
    return inputs, inverse


def fit_indices(inputs: pl.DataFrame, *, size: int = FIT_SAMPLE,
                seed: int = SEED) -> np.ndarray:
    """Sample unique receptors without replacement by stable seeded identity priority.

    BLAKE2b priorities avoid row-position dependence and Python's process-randomized hash.
    Canonical receptor order also fixes PCA's input order, including hash ties.
    """
    priorities = [hashlib.blake2b(
        json.dumps([seed, *row], separators=(",", ":"), ensure_ascii=True).encode("ascii"),
        digest_size=16).digest() for row in inputs.select(INPUT_KEY).iter_rows()]
    return (inputs.with_row_index("__index").with_columns(pl.Series("__priority", priorities))
            .sort("__priority", *INPUT_KEY).head(size).sort(INPUT_KEY)["__index"].to_numpy())


@contextmanager
def _single_threaded_blas():
    """Fix PCA's reduction order and restore the caller's numerical thread limits.

    Randomized PCA changes components with BLAS thread count. Accelerate
    needs its native per-thread API because threadpoolctl does not control it.
    """
    from threadpoolctl import threadpool_limits

    accelerate = None
    if sys.platform == "darwin":
        import ctypes

        accelerate = ctypes.CDLL("/System/Library/Frameworks/Accelerate.framework/Accelerate")
        if not all(hasattr(accelerate, name) for name in ("BLASGetThreading", "BLASSetThreading")):
            raise RuntimeError("Reproducible TCREmp PCA with Accelerate requires macOS 15 or newer")
        accelerate.BLASGetThreading.restype = ctypes.c_uint
        accelerate.BLASSetThreading.argtypes = [ctypes.c_uint]
        previous = accelerate.BLASGetThreading()
        if accelerate.BLASSetThreading(1) != 0:
            raise RuntimeError("Accelerate cannot provide single-threaded PCA")
    try:
        with threadpool_limits(limits=1, user_api="blas"):
            yield
    finally:
        if accelerate is not None:
            accelerate.BLASSetThreading(previous)


@dataclass(frozen=True)
class ReducedEmbedding:
    """In-memory reference prototypes, scaler and PCA; queries only transform.

    No fitted state is stored between builds. Missing V/J calls use the embedding
    producer's allele fallback. Optional species/gene columns must match this fit.
    """

    prototype: TCREmp
    scaler: StandardScaler
    pca: PCA
    species: str
    locus: str

    def transform(self, clonotypes: pl.DataFrame, *, chunk: int = CHUNK) -> np.ndarray:
        """Deduplicate, embed in bounded batches, and restore the original row order."""
        if chunk < 1:
            raise ValueError("chunk must be positive")
        if clonotypes["junction_aa"].is_null().any() or clonotypes["junction_aa"].eq("").any():
            raise ValueError("query junction_aa must not contain null or empty values")
        for column, expected in (("species", self.species), ("gene", self.locus)):
            if column in clonotypes.columns and not clonotypes[column].eq(expected).fill_null(False).all():
                raise ValueError(f"query {column} must match the frozen reference {expected}")
        cols, inverse = receptor_inputs(clonotypes)
        out = np.empty((cols.height, self.pca.n_components_), dtype=np.float64)
        for start in range(0, cols.height, chunk):
            block = cols[start:start + chunk]
            # OpenBLAS CBLAS SGEMM has a CPU-selected direct path outside
            # OPENBLAS_CORETYPE. Float64 keeps PCA on the pinned DGEMM path.
            embedded = self.prototype.embed(block).astype(np.float64)
            with _single_threaded_blas():
                out[start:start + block.height] = self.pca.transform(
                    self.scaler.transform(embedded))
            del embedded
        return out[inverse]


def fit_embedding(clonotypes: pl.DataFrame, species: str, locus: str, *,
                  chunk: int = CHUNK, fit_sample: int = FIT_SAMPLE,
                  n_components: int = N_COMPONENTS,
                  seed: int = SEED) -> tuple[np.ndarray, ReducedEmbedding]:
    """Fit once on reference receptors, returning their coordinates and reusable transform.

    Peptide duplication and input order cannot change the identity-priority reference
    sample. Queries never enter that sample, the scaler or the PCA fit.
    """
    from mir.embedding import TCREmp
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    cols, _ = receptor_inputs(clonotypes)
    model = TCREmp.from_defaults(SPECIES[species], locus)
    idx = fit_indices(cols, size=fit_sample, seed=seed)
    if len(idx) < 2:
        raise ValueError("TCREmp reduction requires at least two unique receptor inputs")
    fit = model.embed(cols[idx]).astype(np.float64)
    with _single_threaded_blas():
        scaler = StandardScaler(copy=False).fit(fit)
        pca = PCA(n_components=min(n_components, fit.shape[1], len(idx) - 1),
                  random_state=seed, copy=False).fit(scaler.transform(fit))
    del fit
    frozen = ReducedEmbedding(model, scaler, pca, species, locus)
    return frozen.transform(clonotypes, chunk=chunk), frozen


def embed_reduced(clonotypes: pl.DataFrame, species: str, locus: str, *,
                  chunk: int = CHUNK, fit_sample: int = FIT_SAMPLE,
                  n_components: int = N_COMPONENTS, seed: int = SEED) -> np.ndarray:
    """Embed reference rows with the unchanged scaler/PCA operation order.

    Use :func:`fit_embedding` when subsequent queries need the same frozen geometry.
    """
    return fit_embedding(clonotypes, species, locus, chunk=chunk, fit_sample=fit_sample,
                         n_components=n_components, seed=seed)[0]


def mean_nn_distance(X: np.ndarray, *, clonotypes: pl.DataFrame | None = None) -> float:
    """Mean distance to the 1st non-self neighbour, over unique receptor inputs when given.

    The reference implementation asks for 4 neighbours, sorts each column and reads column 1.
    Sorting a column does not change its mean, and columns 2 and 3 are never used, so 2 neighbours
    give the identical number for half the work.
    """
    from sklearn.neighbors import NearestNeighbors

    if clonotypes is not None:
        if len(X) != clonotypes.height:
            raise ValueError("coordinates and receptor rows must have equal lengths")
        _, inverse = receptor_inputs(clonotypes)
        X = X[np.unique(inverse, return_index=True)[1]]
    d, _ = NearestNeighbors(n_neighbors=2).fit(X).kneighbors(X)
    return float(d[:, 1].mean())


def chain_eps(X: np.ndarray, coef: float, *, clonotypes: pl.DataFrame | None = None) -> float:
    """The chain-global DBSCAN radius: ``coef * mean(1st-NN distance)``, on the pooled geometry."""
    return round(mean_nn_distance(X, clonotypes=clonotypes) * coef, 3)


def knee_eps_debug(X: np.ndarray, coef: float, k: int = 4, floor_frac: float = 0.40) -> dict:
    """Kneedle, as a diagnostic only. Reports the knee and whether the floor guard fired.

    Kept so the measurement in ROADMAP section 8.3 stays checkable. Nothing in the operating path
    calls this.
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


def knee(curve: np.ndarray, *, concave: bool = True, grid: int = 1000,
         min_strength: float = 0.05, floor_frac: float = 0.05,
         ceil_frac: float = 0.95) -> dict:
    r"""Kneedle on a sorted curve, with the degeneracies that killed the reference version guarded.

    The knee is the point of maximum deviation from the straight chord, in the unit square::

        x = linspace(0, 1, n)          y = (curve - min) / (max - min)
        knee = argmax(y - x)           strength = max(y - x)

    which is Kneedle's difference curve stated directly, with no polynomial fit. The reference
    implementation fits degree 10 and on pooled human TRB -- 112,983 clonotypes -- returns knee
    index 1, fraction 0.000. A degree-10 fit over 113k points oscillates, and the knee it reports is
    the first oscillation, not a feature of the data.

    Three protections, each against an observed failure:

    * **Resample onto ``grid`` points first.** The fit is then independent of ``n``, so the same
      curve shape gives the same knee at 300 points and at 300,000. This is the one that fixes the
      degree-10 oscillation.
    * **``strength`` is the guard, not a separate test.** ``max(y - x)`` is exactly 0 for a straight
      line and rises with how sharp the corner is, so a curve with no knee reports it in the same
      number that locates one. Below ``min_strength`` there is no knee to find, which is the answer
      for a near-linear k-distance curve; the reference version returned an index anyway.
    * **A knee pinned to either end is rejected** (``floor_frac`` / ``ceil_frac``). At the bottom it
      puts ``eps`` below the data and retention collapses to a few per cent; at the top it puts
      ``eps`` above every distance and the epitope percolates into one cluster.

    ``degenerate`` is True when any of the three fires, and ``reason`` names which. The caller
    should fall back rather than use a degenerate knee -- see :func:`chain_eps`.
    """
    n = len(curve)
    if n < 3:
        return {"index": None, "value": None, "strength": 0.0, "frac": None,
                "degenerate": True, "reason": "too-short"}
    c = np.sort(np.asarray(curve, dtype=float))
    lo, hi = float(c[0]), float(c[-1])
    if hi <= lo:
        return {"index": None, "value": None, "strength": 0.0, "frac": None,
                "degenerate": True, "reason": "flat"}

    g = min(grid, n)
    xs = np.linspace(0.0, 1.0, g)
    y = (np.interp(xs, np.linspace(0.0, 1.0, n), c) - lo) / (hi - lo)
    d = (y - xs) if concave else (xs - y)

    i = int(np.argmax(d))
    strength, frac = float(d[i]), i / (g - 1)
    reason = ("no-knee" if strength < min_strength
              else "at-floor" if frac < floor_frac
              else "at-ceiling" if frac > ceil_frac else "")
    idx = round(frac * (n - 1))
    return {"index": idx, "value": float(c[idx]), "strength": strength, "frac": frac,
            "degenerate": bool(reason), "reason": reason}


def _dbscan(X: np.ndarray, epitopes: np.ndarray, eps: float,
            min_samples: int) -> tuple[np.ndarray, np.ndarray]:
    from sklearn.cluster import DBSCAN

    out = np.full(len(X), -1, dtype=np.int64)
    core = np.zeros(len(X), dtype=bool)
    offset = 0
    for ep in np.unique(epitopes):
        rows = np.flatnonzero(epitopes == ep)
        fitted = DBSCAN(eps=eps, min_samples=min_samples).fit(X[rows])
        lab = fitted.labels_
        hit = lab >= 0
        out[rows[hit]] = lab[hit] + offset
        core[rows[fitted.core_sample_indices_]] = True
        offset += int(lab.max()) + 1 if hit.any() else 0
    return out, core


def cluster_labels(X: np.ndarray, epitopes: np.ndarray, eps: float, *,
                   min_samples: int = MIN_SAMPLES) -> np.ndarray:
    """DBSCAN per epitope at one shared ``eps``. ``-1`` is noise, as DBSCAN means it.

    Labels are made unique across epitopes by offsetting, so a label identifies a cluster and not
    just a cluster-within-an-epitope.
    """
    return _dbscan(X, epitopes, eps, min_samples)[0]


@dataclass(frozen=True)
class ClusterReference:
    """Frozen reference membership and per-peptide indexes of retained cluster cores.

    A retained geometric cluster has at least one native length stratum that passes
    ``clusters``' existing min_cluster filter. Its core points can support a query of
    another length without inventing a legacy PWM/display cid for that length.
    """

    members: pl.DataFrame
    labels: np.ndarray
    indexes: dict[str, tuple[NearestNeighbors, np.ndarray]]
    strata: pl.DataFrame
    cluster_names: pl.DataFrame
    eligible_epitopes: tuple[str, ...]
    dimensions: int
    scope: tuple[tuple[str, str], ...]

    def support(self, queries: pl.DataFrame, X: np.ndarray) -> pl.DataFrame:
        """Bounded nearest-core and second-core support in the frozen embedding.

        One native nearest-neighbour batch per peptide uses existing core indexes.
        For distance d and the fitted radius e, support is e/(e+d). The second
        neighbour distinguishes member density when an exact core self-hit has d=0.
        Eligible peptides without retained cores have zero support; others are missing.
        These indices are not calibrated binding probabilities.
        """
        assigned = self.assign(queries, X)  # Shared scope/shape/finite validation.
        X = np.asarray(X)
        near = np.full(queries.height, np.nan)
        density = near.copy()
        peptides = queries['antigen.epitope'].fill_null('').to_numpy()
        for ep in np.unique(peptides):
            rows = np.flatnonzero(peptides == ep)
            if ep not in self.eligible_epitopes:
                continue
            near[rows] = density[rows] = 0.
            if ep not in self.indexes:
                continue
            index, labels = self.indexes[ep]
            distances = index.kneighbors(X[rows], n_neighbors=min(2, len(labels)),
                                         return_distance=True)[0]
            near[rows] = index.radius / (index.radius + distances[:, 0])
            density[rows] = index.radius / (index.radius + distances[:, -1])
        return assigned.with_columns(pl.Series('core_support', near).fill_nan(None),
                                     pl.Series('core_density', density).fill_nan(None))

    def assign(self, queries: pl.DataFrame, X: np.ndarray) -> pl.DataFrame:
        """Return ordered queries with geometric support, ambiguity and legacy cid.

        Queries never fit or extend a cluster. Search only their peptide's reference
        cores; nearest core wins, with the native cluster's smallest cid breaking ties.
        ``matched_clusters`` counts distinct geometric clusters before that choice.
        ``eligible`` marks reference-cohort peptides, including scored nonmembers.
        """
        X = np.asarray(X)
        if X.ndim != 2 or X.shape != (queries.height, self.dimensions):
            raise ValueError("query coordinates must match the row count and reference dimensions")
        if not np.isfinite(X).all():
            raise ValueError("query coordinates must be finite")
        for column, expected in self.scope:
            if column in queries.columns and not queries[column].eq(expected).fill_null(False).all():
                raise ValueError(f"query {column} must match the frozen reference {expected}")
        key = ["antigen.epitope", *INPUT_KEY]
        keys = queries.select(key).with_columns(
            pl.col("v_call", "j_call", "antigen.epitope").fill_null(""))
        if keys["junction_aa"].is_null().any() or keys["junction_aa"].eq("").any():
            raise ValueError("query junction_aa must not contain null or empty values")
        inputs = keys.with_row_index("__row").unique(subset=key, keep="first").sort(key)
        inverse = keys.join(inputs.select(key).with_row_index("__query"), on=key,
                            how="left", validate="m:1", maintain_order="left")["__query"].to_numpy()
        coords = X[inputs["__row"].to_numpy()]
        parts = []
        peptides = inputs["antigen.epitope"].to_numpy()
        for ep in np.unique(peptides):
            if ep not in self.indexes:
                continue
            index, labels = self.indexes[ep]
            rows = np.flatnonzero(peptides == ep)
            edges = index.radius_neighbors_graph(coords[rows], mode="distance").tocoo()
            if edges.nnz:
                parts.append(pl.DataFrame({"__query": rows[edges.row],
                    "__label": labels[edges.col], "__distance": edges.data}))
        result = inputs.with_row_index("__query").with_columns(
            pl.col("antigen.epitope").is_in(self.eligible_epitopes).alias("eligible"),
            pl.col("junction_aa").str.len_chars().alias("__len"))
        if parts:
            matches = (pl.concat(parts).group_by("__query", "__label")
                       .agg(pl.col("__distance").min())
                       .join(self.cluster_names, on="__label", validate="m:1")
                       .sort("__query", "__distance", "__first_cid", "__label"))
            counts = matches.group_by("__query").len(name="matched_clusters")
            chosen = matches.unique(subset="__query", keep="first", maintain_order=True)
            result = (result.join(chosen.select("__query", "__label"), on="__query", how="left")
                      .join(counts, on="__query", how="left")
                      .join(self.strata, on=["antigen.epitope", "__label", "__len"],
                            how="left", validate="m:1"))
        else:
            result = result.with_columns(pl.lit(None, dtype=pl.Int64).alias("__label"),
                pl.lit(None, dtype=pl.UInt32).alias("matched_clusters"),
                pl.lit(None, dtype=pl.String).alias("cid"))
        result = result.sort("__query").with_columns(
            pl.col("__label").fill_null(-1).alias("geometric_label"),
            pl.col("matched_clusters").fill_null(0), pl.col("cid").fill_null(""))
        result = result.with_columns((pl.col("geometric_label") >= 0).alias("geometric_support"),
                                     (pl.col("matched_clusters") > 1).alias("ambiguous"))
        return queries.with_columns([pl.Series(c, result[c].to_numpy()[inverse]) for c in
            ["eligible", "geometric_label", "geometric_support", "matched_clusters", "ambiguous", "cid"]])


def fit_clusters(cohort_chain: pl.DataFrame, X: np.ndarray, eps: float, *,
                 min_samples: int = MIN_SAMPLES, min_cluster: int = 5) -> ClusterReference:
    """Fit reference DBSCAN once and retain its original members and query core indexes.

    Pass one eligible species/chain cohort, as returned by :func:`cohort`, and its
    frozen reference coordinates. Query peptide labels never affect this fit.
    """
    from sklearn.neighbors import NearestNeighbors

    X = np.asarray(X)
    if (X.ndim != 2 or len(X) != cohort_chain.height or not len(X)
            or X.shape[1] < 2 or not np.isfinite(X).all()):
        raise ValueError("reference coordinates must be finite, nonempty, match the rows "
                         "and have at least two dimensions")
    scope = tuple((c, cohort_chain[c].item(0)) for c in ("species", "gene")
                  if c in cohort_chain.columns and cohort_chain.height)
    if any(not cohort_chain[c].eq(value).fill_null(False).all() for c, value in scope):
        raise ValueError("fit_clusters requires one species and chain")
    epitopes = cohort_chain["antigen.epitope"].to_numpy()
    labels, core = _dbscan(X, epitopes, eps, min_samples)
    members = clusters(cohort_chain, X, labels=labels, min_cluster=min_cluster)
    key = ["antigen.epitope", *INPUT_KEY]
    member_ids = (members.select(*key, "cid") if "cid" in members.columns else
                  members.select(key).with_columns(pl.lit("", dtype=pl.String).alias("cid")))
    strata = (cohort_chain.select(key).with_columns(pl.Series("__label", labels),
                pl.col("junction_aa").str.len_chars().alias("__len"))
              .join(member_ids.unique(), on=key, how="inner", validate="m:1")
              .select("antigen.epitope", "__label", "__len", "cid").unique())
    names = strata.group_by("__label").agg(pl.col("cid").min().alias("__first_cid"))
    retained = names["__label"].to_numpy()
    indexes = {}
    for ep in np.unique(epitopes):
        pick = (epitopes == ep) & core & np.isin(labels, retained)
        if pick.any():
            points = np.ascontiguousarray(X[pick])
            points.setflags(write=False)
            core_labels = labels[pick].copy()
            core_labels.setflags(write=False)
            indexes[ep] = (NearestNeighbors(radius=eps, n_jobs=1).fit(points), core_labels)
    labels.setflags(write=False)
    return ClusterReference(members, labels, indexes, strata, names,
                            tuple(np.unique(epitopes)), X.shape[1], scope)


def hdbscan_labels(X: np.ndarray, epitopes: np.ndarray, min_cluster_size: int = 5, *,
                   min_samples: int = MIN_SAMPLES, selection_epsilon: float = 0.0,
                   method: str = "eom") -> np.ndarray:
    """HDBSCAN per epitope. Same contract as :func:`cluster_labels`: ``-1`` is noise, labels unique.

    ⚠ Measured, not enabled. The shipped path is :func:`cluster_labels`; this is here so the
    comparison is runnable. Judge it on ``Q`` (:mod:`vdjdb.validate.qscore`) and the section 11.1
    lift, never on lift alone -- density methods buy lift by shattering, and shattering is the
    failure mode ``Q`` exists to catch.

    DBSCAN commits to one radius for every epitope, and VDJdb's epitopes do not share a density: a
    display-selected library yields thousands of receptors one substitution apart by construction
    (:data:`DISPLAY`, 29,688 of 192,753 records) while a 30-record epitope is sparse. HDBSCAN
    condenses a hierarchy and selects by cluster stability instead, so there is no global radius to
    be wrong. It also folds ``min_samples`` and the ``min_cluster`` post-filter into the single
    ``min_cluster_size``.

    ``selection_epsilon`` is the anti-shatter knob -- subclusters closer than it are merged back --
    and ``method="leaf"`` is its opposite, taking every leaf of the condensed tree. Default ``"eom"``
    (excess of mass) is the parsimonious one.

    ``sklearn.cluster.HDBSCAN``, so no dependency beyond the one already used for PCA and DBSCAN.
    """
    from sklearn.cluster import HDBSCAN

    out = np.full(len(X), -1, dtype=np.int64)
    offset = 0
    for ep in np.unique(epitopes):               # np.unique sorts: the offset must not vary
        m = epitopes == ep
        sub = X[m]
        if len(sub) < min_cluster_size:
            continue
        lab = HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples,
                      cluster_selection_epsilon=selection_epsilon,
                      cluster_selection_method=method, copy=True).fit_predict(sub)
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

#: Fitted under the two-stage rule in ``docs/denoising.md`` section 7.1 -- admissible on ``Q``,
#: purity, precision and epitope coverage against the shipped annotation, then ranked on the
#: section 11.1 independent-study lift. Not the purity/retention criterion sections 31-34 used: that
#: pair trades against itself and had no interior optimum (section 36).
#:
#: The frontier is a single crossing, not a search. Lift falls monotonically as ``coef`` widens
#: while ``Q``, retention and epitope coverage all rise, so the answer is the smallest ``coef``
#: whose four admissibility quantities all clear the bar -- verified by probing the crossing at 0.05
#: resolution rather than trusting a grid point. ``min_cluster`` 5 dominates 3 on lift, purity and
#: precision at every ``coef`` measured, so it is not a trade either.
#:
#: Measured on the ``cluster_members_tcremp.txt`` the build writes, over human clonotype-epitope pairs
#: in epitopes with >= 30 records -- the same cohort TCRNET is scored on (ROADMAP_local section 37.1):
#:
#: ======  =========  ======  ======  ======  ======  ======  ======  =========  =========
#: chain   config     lift    legacy  Q       legacy  purity  legacy  epitopes   legacy
#: ======  =========  ======  ======  ======  ======  ======  ======  =========  =========
#: TRA     coef 1.8   1.855   1.703   0.1729  0.1691  0.9105  0.8658  104 / 118  103
#: TRB     coef 1.55  3.585   2.855   0.4947  0.4433  0.9829  0.9790  105 / 178  103
#: ======  =========  ======  ======  ======  ======  ======  ======  =========  =========
#:
#:
#: ⚠ **Measured on the corpus of 2026-09-26 and not restated since.** These numbers chose the
#: shipped parameter and the comparison they make is internally consistent, but the legacy
#: column moves whenever a chunk merges - the released file is fixed and the cohort it is
#: scored against is not. The current value of every axis, for the released files and for this
#: build, is measured on every build into ``out/reports/motif-metrics.tsv`` and gated against
#: ``rules/motif_metrics.tsv`` (``vdjdb motif-metrics``). Read those for today's figures.
#: Precision: TRA 0.8931 against legacy's 0.8567, TRB 0.9829 against 0.9756. Retention: TRA 0.2508
#: against 0.2105, TRB 0.3747 against 0.3218. Both chains improve on all five pooled axes and on
#: epitope coverage, with no cost to name here, unlike TCRNET's TRB cell.
#:
#: ⚠ TRB was ``coef`` 1.15 until the per-epitope breakdown was run. 1.15 maximises pooled lift
#: (4.490) and is admissible on ``Q``, purity and precision, but it covers 88 of 178 epitopes
#: against legacy's 103, and drops from 42 to 36 epitopes where the clustering beats local chance.
#: It clusters essentially the same total clonotypes (37,094 against 37,210) concentrated into
#: fifteen fewer epitopes and 746 clusters instead of 1,074. An epitope with no motif gets no
#: denoising, so that is a worse database bought with a better average. Coverage is an admissibility
#: axis for that reason (``docs/denoising.md`` section 7.1), and 1.55 is the smallest ``coef``
#: clearing all four. The cost is pooled lift, 4.490 -> 3.585, still +25.6 % on legacy.
#:
#: ⚠ Lift is on the non-display denominator (``docs/denoising.md`` section 6.1). The same TRB
#: clustering reads 3.585 there and a much lower figure on the full cohort; display-selected
#: records contribute zero independently-replicated pairs while filling 29,692 of 116,053 clonotype
#: slots. A lift figure needs its cohort named.
#:
#: ``n_components`` is 50 on both chains. Section 34 had TRB at 100; that came from the superseded
#: purity criterion, and 50 is the benchmark's own front-end applied to every method, so the
#: representation is the only variable.
#:
#: ⚠ The published 0.75 does not transfer in either direction -- different embedding (standalone
#: `tcremp`, ~3,000 OLGA prototypes, Smith-Waterman; ROADMAP section 8.2).
TUNED: dict[str, dict] = {
    "TRA": {"coef": 1.8, "min_cluster": 5, "n_components": 50},
    "TRB": {"coef": 1.55, "min_cluster": 5, "n_components": 50},
}

#: Back-compatible view of :data:`TUNED` for callers that only want the radius.
COEF: dict[str, float] = {g: c["coef"] for g, c in TUNED.items()}


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

    A clustering that finds convergent selection should preferentially recover exactly the
    clonotypes a second laboratory saw. Clustering everything wins recall and loses precision;
    clustering nothing scores zero. F1 peaks where the clustering is selective about the right
    thing, and F1 is the only number the sweep ranks on.
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

    ⚠ Never fitted against TCRvdb. That set is held out, touched once, and only in aggregate
    (ROADMAP section 11.2). Tuning on it would spend the project's one independent read.

    The mean 1st-NN distance is computed once -- it does not depend on ``coef`` -- so the sweep
    costs one DBSCAN pass per grid point, not one nearest-neighbour search per point.
    """
    epitopes = cohort_chain["antigen.epitope"].to_numpy()
    base = mean_nn_distance(X, clonotypes=cohort_chain)
    rows = []
    for coef in grid:
        labels = cluster_labels(X, epitopes, round(base * coef, 3), min_samples=min_samples)
        rows.append({"coef": coef, "eps": round(base * coef, 3),
                     **_objective(labels, is_replicated)})
    return pl.DataFrame(rows).sort("f1", descending=True)


# ---------------------------------------------------------------------------------------------
# Clustering
# ---------------------------------------------------------------------------------------------

def clusters(cohort_chain: pl.DataFrame, X: np.ndarray, eps: float | None = None, *,
             min_samples: int = MIN_SAMPLES, min_cluster: int = 5,
             labels: np.ndarray | None = None) -> pl.DataFrame:
    """Cluster one chain and label it in the shape :mod:`vdjdb.motifs.emit` expects.

    Takes either ``eps`` -- per-epitope DBSCAN at that radius, the shipped path -- or precomputed
    ``labels``. The second form is how an alternative algorithm is measured through this same cid
    machinery instead of growing a second copy of it: :func:`hdbscan_labels` and
    :func:`vdjdb.motifs.cluster._leiden` both produce ``labels`` in the one contract, ``-1`` for
    noise and otherwise unique across epitopes.

    ``cid`` takes an ``L<len>`` suffix: one legacy cid per (cluster, CDR3 length). A DBSCAN cluster
    in embedding space may span lengths, where a PWM may not, and ``vdjdb-web`` splits every cid by
    ``len`` before building a cluster anyway; without the suffix two display clusters would share one
    ``clusterId`` and the reader's ``strict = true`` path would break. Measured on the REDCEA
    production files, this is the common case: 809 of 847 TRA and 1,002 of 1,082 TRB cids span more
    than one length (ROADMAP section 8.6).

    ``x``/``y`` are the first two principal components -- the embedding's own layout, so no separate
    force-directed pass is needed and the coordinates are interpretable.
    """
    from .cluster import _INITIAL, _repr_allele

    if labels is None:
        if eps is None:
            raise ValueError("clusters() needs either eps or precomputed labels")
        labels = cluster_labels(X, cohort_chain["antigen.epitope"].to_numpy(), eps,
                               min_samples=min_samples)
    g = (cohort_chain
         .with_columns(pl.Series("__label", labels),
                       pl.Series("x", X[:, 0]), pl.Series("y", X[:, 1]),
                       pl.col("junction_aa").str.len_chars().alias("__len"))
         .filter(pl.col("__label") >= 0))
    if not g.height:
        return g.drop("__label", "__len")

    # The stratum, not the DBSCAN cluster, is what gets a cid -- so `csz` is unambiguous and
    # `freq` in the PWM is a distribution.
    sizes = g.group_by(["species", "gene", "antigen.epitope", "__label", "__len"],
                       maintain_order=True).agg(
        pl.len().alias("csz"), pl.col("junction_aa").min().alias("__first"))
    out = []
    for (species, gene, epitope), grp in sizes.group_by(
            ["species", "gene", "antigen.epitope"], maintain_order=True):
        # Size descending then the smallest member: a number that follows the content, so a cluster
        # whose membership is unchanged keeps its id without anything being stored (hard rule 9).
        #
        # `__len` and `__label` are the tiebreak, and they are what makes the key total: within one
        # epitope a stratum is identified by `(__label, __len)`, so no two rows here can tie on all
        # four. Without them the ordering was decided by whatever order the sort left equal keys in,
        # and that flipped between runs. Measured 2026-09-27 on the 2026-06-03 corpus: exactly one
        # tie per TRA chain (human `GILGFVFTL`, two strata at csz 47 both starting `CAAGGSQGNLIF`;
        # mouse `VEALYLVSG`, two at csz 5) and none on either TRB, which is why TCREMP's TRA files
        # differed between two runs of the same build and its TRB files did not.
        order = (grp.filter(pl.col("csz") >= min_cluster)
                    .sort(["csz", "__first", "__len", "__label"],
                          descending=[True, False, False, False])
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
