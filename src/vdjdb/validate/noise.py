"""How much of a motif is chance, and how much of a lift is publicity.

Two quantities that `docs/denoising.md` argues from, computed rather than asserted.

**The chance-recruitment rate** says what fraction of a motif's members a *non-convergent* record
would supply for free. TCRNET already estimates the background neighbour probability per clonotype,
so this costs a binomial tail and nothing else -- and it is what turns "a wider ball recruits
bystanders" from a warning into a number.

**The publicity-controlled lift** says how much of the independent-study enrichment survives once
generation probability is held fixed. A public clonotype is more likely both to be seen by a second
laboratory *and* to have sequence neighbours, so the raw lift is confounded upward. The control is
the IMMREP25 audit's: permute the label within strata, so the null keeps the covariate and destroys
only the association.

⚠ **Nothing here reads TCRvdb.**
"""
from __future__ import annotations

import math

import numpy as np
import polars as pl

from ..config import SEED


def ball_volume(length: int, subs: int) -> int:
    r"""Number of peptides within ``subs`` substitutions of a length-``length`` CDR3.

    :math:`V_s(L) = \sum_{j \le s} \binom{L}{j} 19^j`. The point of having it is the **ratio**
    between two scopes: it is the factor by which the chance-neighbour probability inflates when the
    scope widens, and it is why scope is the most dangerous knob in the motif stage.
    """
    return sum(math.comb(length, j) * 19 ** j for j in range(subs + 1))


def chance_recruitment(n_sample: int, p_background: float, min_degree: int = 2) -> float:
    r"""Probability a non-convergent clonotype still clears the degree floor, by chance alone.

    A record that does not belong to the epitope has neighbours in the sample at the background
    rate, so its degree is :math:`D \sim \mathrm{Bin}(n-1, \hat p)` with
    :math:`\hat p = (n_{\mathrm{control}}+1)/(M+1)` -- TCRNET's own pseudocounted estimate, already
    computed per clonotype. It is recruited when :math:`D \ge d_{\min}`.

    Note the dependence on ``n_sample``: the same background rate recruits far more freely in a
    large epitope than a small one, so a single p-value threshold is not a single noise level.
    """
    from scipy.stats import binom

    return float(binom.sf(min_degree - 1, max(n_sample - 1, 0), p_background))


def _lift(clustered: np.ndarray, replicated: np.ndarray) -> float:
    base = replicated.mean()
    if not clustered.any() or base == 0:
        return 0.0
    return float(replicated[clustered].mean() / base)


def controlled_lift(df: pl.DataFrame, *, clustered: str = "clustered",
                    replicated: str = "replicated", stratum: str = "stratum",
                    n_perm: int = 1000, seed: int = SEED) -> dict:
    r"""Raw lift, the within-stratum permutation null, and the ratio of the two.

    ``stratum`` holds the covariate to hold fixed -- generation-probability decile, via
    :func:`pgen_stratum`. The null permutes ``replicated`` **within** each stratum, so a clonotype's
    publicity is preserved and only its association with the clustering is broken. ``ratio`` above 1
    is enrichment that publicity does not account for; a ratio near 1 means the raw lift was the
    covariate.

    **The permutation is drawn in closed form rather than performed.** Permuting within a stratum
    preserves that stratum's replicated count, so the lift's denominator -- the overall base rate --
    is identical in every draw, and its numerator is the clustered-and-replicated count, which is a
    sum of independent :math:`\mathrm{Hypergeometric}(n_g, k_g, c_g)` draws over strata
    (:math:`n_g` records, :math:`k_g` replicated, :math:`c_g` clustered). Sampling those directly is
    the same distribution exactly, not an approximation, and measured **300x faster** than shuffling
    the label vector -- 0.66 s to 2.2 ms at 2,000 draws over 20,000 rows in 10 strata.
    """
    rng = np.random.default_rng(seed)
    c = df[clustered].to_numpy().astype(bool)
    r = df[replicated].to_numpy().astype(bool)
    s = df[stratum].to_numpy()
    raw = _lift(c, r)

    n_clustered, base = int(c.sum()), float(r.mean())
    if not n_clustered or base == 0:
        return {"lift": raw, "null_mean": 0.0, "null_sd": 0.0, "ratio": 0.0, "p": 1.0,
                "n": int(df.height), "clustered": n_clustered, "replicated": int(r.sum())}

    v = np.unique(s)
    n_g = np.array([(s == x).sum() for x in v])
    k_g = np.array([(r & (s == x)).sum() for x in v])
    c_g = np.array([(c & (s == x)).sum() for x in v])
    # A stratum with nothing clustered contributes no draw; numpy requires nsample >= 1.
    m = c_g > 0
    tp = rng.hypergeometric(k_g[m], (n_g - k_g)[m], c_g[m], size=(n_perm, int(m.sum()))).sum(axis=1)
    null = (tp / n_clustered) / base

    mean = float(null.mean())
    return {"lift": raw, "null_mean": mean, "null_sd": float(null.std()),
            "ratio": raw / mean if mean else 0.0,
            "p": float((null >= raw).mean()), "n": int(df.height),
            "clustered": n_clustered, "replicated": int(r.sum())}


def pgen_stratum(df: pl.DataFrame, column: str = "cdr3nt.pgen", *, bins: int = 10) -> pl.Series:
    """Generation-probability decile, as the stratum to control publicity on.

    ``log10`` of the generation probability, cut into ``bins`` equal-count bins. Rows with no Pgen
    get their own stratum rather than being dropped -- missing a covariate is not a reason to lose
    a record, and a separate stratum absorbs whatever they have in common.
    """
    lg = (pl.col(column).cast(pl.Float64).log10()
          .replace([float("inf"), float("-inf")], None))
    q = df.select(lg.qcut(bins, labels=[str(i) for i in range(bins)], allow_duplicates=True)
                  .cast(pl.Utf8).fill_null("missing").alias("stratum"))
    return q["stratum"]
