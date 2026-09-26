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


def expected_false_members(n_sample: int, p_background: float, true_fraction: float, *,
                           min_degree: int = 2) -> float:
    r"""Expected count of noise clonotypes recruited into motifs:
    :math:`(1-\phi)\,n\,\alpha`, with :math:`\alpha` from :func:`chance_recruitment`.

    ``true_fraction`` is :math:`\phi`, the share of the epitope's records that really do bind it.
    It is not known -- that is the whole problem -- so this is used as a *sensitivity curve* over
    plausible :math:`\phi`, never as a point estimate.
    """
    return (1.0 - true_fraction) * n_sample * chance_recruitment(n_sample, p_background,
                                                                 min_degree)


def _lift(clustered: np.ndarray, replicated: np.ndarray) -> float:
    base = replicated.mean()
    if not clustered.any() or base == 0:
        return 0.0
    return float(replicated[clustered].mean() / base)


def controlled_lift(df: pl.DataFrame, *, clustered: str = "clustered",
                    replicated: str = "replicated", stratum: str = "stratum",
                    n_perm: int = 1000, seed: int = SEED) -> dict:
    """Raw lift, the within-stratum permutation null, and the ratio of the two.

    ``stratum`` holds the covariate to hold fixed -- generation-probability decile, via
    :func:`pgen_stratum`. The null permutes ``replicated`` **within** each stratum, so a clonotype's
    publicity is preserved and only its association with the clustering is broken. ``ratio`` above 1
    is enrichment that publicity does not account for; a ratio near 1 means the raw lift was the
    covariate.
    """
    rng = np.random.default_rng(seed)
    c = df[clustered].to_numpy().astype(bool)
    r = df[replicated].to_numpy().astype(bool)
    s = df[stratum].to_numpy()
    raw = _lift(c, r)

    idx = [np.flatnonzero(s == v) for v in np.unique(s)]
    null = np.empty(n_perm)
    for t in range(n_perm):
        perm = r.copy()
        for g in idx:
            perm[g] = rng.permutation(r[g])
        null[t] = _lift(c, perm)
    m = float(null.mean())
    return {"lift": raw, "null_mean": m, "null_sd": float(null.std()),
            "ratio": raw / m if m else 0.0,
            "p": float((null >= raw).mean()), "n": int(df.height),
            "clustered": int(c.sum()), "replicated": int(r.sum())}


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
