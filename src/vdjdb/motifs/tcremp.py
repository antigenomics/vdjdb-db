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


def knee(curve: np.ndarray, *, concave: bool = True, grid: int = 1000,
         min_strength: float = 0.05, floor_frac: float = 0.05,
         ceil_frac: float = 0.95) -> dict:
    r"""Kneedle on a sorted curve, with the degeneracies that killed the reference version guarded.

    The knee is the point of maximum deviation from the straight chord, in the unit square::

        x = linspace(0, 1, n)          y = (curve - min) / (max - min)
        knee = argmax(y - x)           strength = max(y - x)

    which is Kneedle's difference curve stated directly, with **no polynomial fit**. That matters:
    the reference implementation fits degree 10 and on pooled human TRB -- 112,983 clonotypes --
    returns knee index 1, fraction 0.000. A degree-10 fit over 113k points oscillates, and the knee
    it reports is the first oscillation, not a feature of the data.

    Three protections, each against an observed failure rather than an imagined one:

    * **Resample onto ``grid`` points first.** The fit is then independent of ``n``, so the same
      curve shape gives the same knee whether it carries 300 points or 300,000. This is the one that
      fixes the degree-10 oscillation.
    * **``strength`` is the guard, not a separate test.** ``max(y - x)`` is exactly 0 for a straight
      line and rises with how sharp the corner is, so a curve with no knee reports it in the same
      number that locates one. Below ``min_strength`` there is no knee to find -- which is the honest
      answer for a near-linear k-distance curve, where the reference version returned an index
      anyway.
    * **A knee pinned to either end is rejected** (``floor_frac`` / ``ceil_frac``). At the bottom it
      puts ``eps`` below the data and retention collapses to a few per cent; at the top it puts
      ``eps`` above every distance and the epitope percolates into one cluster.

    ``degenerate`` is True when any of the three fires, and ``reason`` names which. The caller is
    expected to fall back rather than to use a degenerate knee -- see :func:`chain_eps`.
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



def hdbscan_labels(X: np.ndarray, epitopes: np.ndarray, min_cluster_size: int = 5, *,
                   min_samples: int = MIN_SAMPLES, selection_epsilon: float = 0.0,
                   method: str = "eom") -> np.ndarray:
    """HDBSCAN per epitope. Same contract as :func:`cluster_labels`: ``-1`` is noise, labels unique.

    ⚠ **Measured, not enabled.** The shipped path is :func:`cluster_labels`; this is here so the
    comparison is runnable rather than argued. Judge it on ``Q``
    (:mod:`vdjdb.validate.qscore`) **and** the section 11.1 lift, never on lift alone -- density
    methods buy lift by shattering, and shattering is the failure mode ``Q`` exists to catch.

    The reason to have it at all is that DBSCAN commits to **one radius for every epitope**, and
    VDJdb's epitopes do not share a density: a display-selected library yields thousands of receptors
    one substitution apart by construction (:data:`DISPLAY`, 29,688 of 192,753 records) while a
    30-record epitope is sparse. HDBSCAN condenses a hierarchy and selects by cluster stability
    instead, so there is no global radius to be wrong. It also folds ``min_samples`` and the
    ``min_cluster`` post-filter into the single ``min_cluster_size``.

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

#: **Fitted under the two-stage rule in ``docs/denoising.md`` section 7.1** -- admissible on ``Q``,
#: purity and precision against the shipped annotation, then ranked on the section 11.1
#: independent-study lift. Not the purity/retention criterion sections 31-34 used: that pair trades
#: against itself and had no interior optimum (section 36).
#:
#: **The frontier is a single crossing, not a search.** Lift falls monotonically as ``coef`` widens
#: and ``Q`` rises, so the answer is the *smallest* ``coef`` whose ``Q`` clears the legacy bar --
#: verified by probing the crossing at 0.05 resolution rather than trusting a grid point.
#: ``min_cluster`` 5 dominates 3 on lift, purity **and** precision at every ``coef`` measured, so it
#: is not a trade either.
#:
#: Re-measured on the ``cluster_members_tcremp.txt`` the build actually writes, over human
#: clonotype-epitope pairs in epitopes with >= 30 records -- the same cohort TCRNET is scored on
#: (ROADMAP_local section 37.1):
#:
#: ======  =========  ======  ======  ======  ======  ======  ======  ======  ======
#: chain   config     lift    legacy  Q       legacy  purity  legacy  prec    legacy
#: ======  =========  ======  ======  ======  ======  ======  ======  ======  ======
#: TRA     coef 1.8   1.855   1.703   0.1728  0.1691  0.9104  0.8658  0.8931  0.8567
#: TRB     coef 1.15  4.490   2.855   0.4437  0.4433  0.9873  0.9790  0.9859  0.9756
#: ======  =========  ======  ======  ======  ======  ======  ======  ======  ======
#:
#: Retention: TRA 0.2508 against legacy's 0.2105, TRB 0.3232 against 0.3218.
#:
#: **Both chains improve on every one of the five axes** -- there is no cost to name here, unlike
#: TCRNET's TRB cell. TRB gains **+57 % lift** over the shipped annotation, TRA **+8.9 %**.
#:
#: ⚠ **TRB sits on the Q boundary**: 0.4437 against a bar of 0.4433, a margin of +0.0004. It clears
#: the stated criterion and is therefore the answer the rule gives, but the margin is rounding-level,
#: so re-check this cell whenever the corpus changes. ``coef`` 1.2 is the nearest cell with a real
#: margin (Q 0.4670) and costs 3.5 % of lift.
#:
#: ⚠ **Lift is on the non-display denominator** (``docs/denoising.md`` section 6.1). The same TRB
#: clustering reads **4.490** there and **1.272** on the full cohort; display-selected records
#: contribute zero independently-replicated pairs while filling 29,692 of 116,053 clonotype slots.
#: A lift figure without its cohort is not a number.
#:
#: ``n_components`` is 50 on both chains. Section 34 had TRB at 100; that came from the superseded
#: purity criterion, and 50 is the benchmark's own front-end applied to every method, so the
#: representation is the only variable.
#:
#: ⚠ The published **0.75 does not transfer** in either direction -- different embedding (standalone
#: `tcremp`, ~3,000 OLGA prototypes, Smith-Waterman; ROADMAP section 8.2).
TUNED: dict[str, dict] = {
    "TRA": {"coef": 1.8, "min_cluster": 5, "n_components": 50},
    "TRB": {"coef": 1.15, "min_cluster": 5, "n_components": 50},
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

def clusters(cohort_chain: pl.DataFrame, X: np.ndarray, eps: float | None = None, *,
             min_samples: int = MIN_SAMPLES, min_cluster: int = 5,
             labels: np.ndarray | None = None) -> pl.DataFrame:
    """Cluster one chain and label it in the shape :mod:`vdjdb.motifs.emit` expects.

    Takes **either** ``eps`` -- per-epitope DBSCAN at that radius, the shipped path -- **or**
    precomputed ``labels``. The second form is how an alternative algorithm is measured through this
    same cid machinery instead of growing a second copy of it: :func:`hdbscan_labels` and
    :func:`vdjdb.motifs.cluster._leiden` both produce ``labels`` in the one contract, ``-1`` for
    noise and otherwise unique across epitopes.

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
