"""Per-cluster position weight matrices, and the background they are read against.

Stage three. A cluster is a set of CDR3s; its logo is one column per position, and each column is
read *against* what the background repertoire puts there, so that a residue the germline supplies
anyway does not look informative.

**The shipped files delete letters, and this is why they do not.** The legacy pipeline left-joined
each cluster column against a frozen background PWM and then applied ``filter(total.bg > 0)``. Any
``(pos, aa)`` the background had never seen at that ``(v, j, len)`` was dropped -- which is exactly
the rarest and most informative residues. Measured on the 2026-06-03 release: **1.00 % of letter
mass deleted across 253 of 589 clusters**, per-position ``freq`` then summing to less than 1 (minimum
observed 0.111), and ``I`` computed on a truncated distribution and therefore **inflated**. A further
**31 of the cids in `cluster_members.txt` lost every row**, so the web draws them in the tree with no
logo. ``need.impute`` was computed *after* the filter, so it was ``FALSE`` on all 13,456 rows and the
imputation machinery was dead code. ROADMAP section 8.5.

:func:`cluster_pwms` replaces the filter with a three-level cascade -- ``(v, j, len)`` -> ``(len)``
-> uniform -- taking the finest level that has any observation at that position, with a Laplace
pseudocount so no background frequency is ever zero. Nothing is dropped, ``freq`` sums to 1 by
construction, and ``need.impute`` records which level was actually used.

⚠ **Assert the sign, do not chase the shipped numbers.** At the 253 affected clusters the new
``sum(I)`` must come out *lower* than the shipped value, because the shipped one was computed over a
distribution missing its tail. A rewrite that reproduces the shipped numbers has reproduced the bug.
:func:`information_delta` is that check.
"""
from __future__ import annotations

import math

import numpy as np
import polars as pl

#: The productive 20, in the order a logo's columns are indexed.
ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
_CODE = np.full(256, -1, dtype=np.int8)
for _i, _a in enumerate(ALPHABET):
    _CODE[ord(_a)] = _i

#: ``log(20)`` -- information is reported on a 0..1 scale where 1 is a fully determined column.
_LOG_N = math.log(len(ALPHABET))

#: Laplace pseudocount added to every background cell before it becomes a frequency. One
#: observation per residue: enough that an unseen residue is rare rather than impossible, small
#: enough that it does not move a stratum with thousands of observations.
PSEUDOCOUNT = 1.0

#: Which level of the cascade supplied a column's background, recorded per row so the provenance is
#: in the file rather than in a comment. Replaces the legacy ``need.impute``, which was always FALSE.
LEVELS = ("vj_len", "len", "uniform")


def counts(seqs: list[str], length: int) -> np.ndarray:
    """``(length, 20)`` position x residue counts over equal-length ``seqs``.

    One ``frombuffer`` and one ``bincount`` per column -- no Python loop over sequences. Residues
    outside :data:`ALPHABET` are dropped by the ``-1`` code, which cannot happen on ``chunks/``
    (the unconventional-AA records are quarantined out of the build) but would silently corrupt a
    column if it did.
    """
    if not seqs:
        return np.zeros((length, len(ALPHABET)), dtype=np.int64)
    m = _CODE[np.frombuffer("".join(seqs).encode("ascii"), dtype=np.uint8)].reshape(-1, length)
    return np.stack([np.bincount(c[c >= 0], minlength=len(ALPHABET)) for c in m.T])


def background_pwms(background: pl.DataFrame, strata: pl.DataFrame) -> dict:
    """Background counts for every stratum a cluster will ask for, at both levels of the cascade.

    ``background`` is a control repertoire with ``cdr3aa``, ``v`` and ``j``; ``strata`` names the
    ``(v.segm.repr, j.segm.repr, len)`` triples in use. Returns
    ``{("vj_len", v, j, len): array, ("len", len): array}``.

    The background's calls are **gene-level** (``TRBV24-1``) while a cluster's representative is an
    allele (``TRBV24-1*01``), so the join is on the gene. Matching the strings as given would miss
    every triple and silently collapse the cascade to its ``len`` level.

    Only the strata actually in use are materialised -- a full ``(v, j, len, pos, aa)`` explode over
    a 15M-row control is ~225M rows for the handful of triples the clusters need (CLAUDE.md
    section 10, filter before materializing).
    """
    bg = background.select(
        pl.col("cdr3aa"),
        pl.col("v").str.split("*").list.first().alias("v.gene"),
        pl.col("j").str.split("*").list.first().alias("j.gene"),
        pl.col("cdr3aa").str.len_chars().alias("len"),
    )
    want = strata.select(
        pl.col("v.segm.repr").str.split("*").list.first().alias("v.gene"),
        pl.col("j.segm.repr").str.split("*").list.first().alias("j.gene"),
        pl.col("len"),
    ).unique()

    out: dict = {}
    for (length,), grp in bg.filter(
            pl.col("len").is_in(want["len"].unique().implode())).group_by(["len"], maintain_order=True):
        out[("len", length)] = counts(grp["cdr3aa"].to_list(), length)
        for (v, j), sub in grp.join(want, on=["v.gene", "j.gene", "len"]).group_by(
                ["v.gene", "j.gene"], maintain_order=True):
            out[("vj_len", v, j, length)] = counts(sub["cdr3aa"].to_list(), length)
    return out


def _background_column(bg: dict, v: str, j: str, length: int,
                       pos: int) -> tuple[np.ndarray, np.ndarray, str]:
    """The finest background column at ``pos``, the coarse one beside it, and the level used.

    The coarse column is the ``len`` level, which the legacy file carries as ``count.bg.i`` /
    ``total.bg.i`` -- there to be imputed *from*, though the legacy never did because it computed
    ``need.impute`` after the filter that would have set it.
    """
    coarse = bg.get(("len", length))
    if coarse is None or not coarse[pos].sum():
        coarse = np.ones((length, len(ALPHABET)))
    vj = bg.get(("vj_len", v.split("*")[0], j.split("*")[0], length))
    if vj is not None and vj[pos].sum():
        return vj[pos].astype(float), coarse[pos].astype(float), "vj_len"
    if bg.get(("len", length)) is not None and bg[("len", length)][pos].sum():
        return coarse[pos].astype(float), coarse[pos].astype(float), "len"
    return np.ones(len(ALPHABET)), coarse[pos].astype(float), "uniform"


def _information(freq: np.ndarray) -> float:
    """``1 + sum(p log p) / log 20`` -- 0 for a uniform column, 1 for a fully determined one."""
    nz = freq[freq > 0]
    return 1.0 + float((nz * np.log(nz)).sum()) / _LOG_N


def cluster_pwms(members: pl.DataFrame, background: pl.DataFrame) -> pl.DataFrame:
    """One row per ``(cid, len, pos, aa)`` with counts, frequencies and information.

    ``members`` is the output of :func:`vdjdb.motifs.cluster.clusters`. A cluster spanning several
    CDR3 lengths gets one **stratum** per length: the PWM is a fixed-width object and a column is
    only meaningful among sequences of the same length (ROADMAP section 8.6).

    ``I`` is the column's own information; ``I.norm`` is it net of the background column's, so a
    position the repertoire already fixes scores near zero. ``height.*`` are the logo letter heights,
    ``freq`` times the respective information.
    """
    m = members.with_columns(pl.col("junction_aa").str.len_chars().alias("len"))
    strata = m.select("v.segm.repr", "j.segm.repr", "len").unique()
    bg = background_pwms(background, strata)

    rows: list[dict] = []
    for keys, grp in m.group_by(["species", "gene", "antigen.epitope", "cid", "len"],
                                maintain_order=True):
        species, gene, epitope, cid, length = keys
        seqs = sorted(grp["junction_aa"].to_list())
        obs = counts(seqs, length)
        csz, v, j = grp.height, grp["v.segm.repr"][0], grp["j.segm.repr"][0]

        for pos in range(length):
            col = obs[pos].astype(float)
            freq = col / col.sum()
            bg_col, coarse_col, level = _background_column(bg, v, j, length, pos)
            total_bg = float(bg_col.sum())
            freq_bg = (bg_col + PSEUDOCOUNT) / (total_bg + PSEUDOCOUNT * len(ALPHABET))
            info = _information(freq)
            info_norm = info - _information(freq_bg)

            for k, aa in enumerate(ALPHABET):
                if not col[k]:
                    continue        # a residue the cluster never shows has no letter in the logo
                rows.append({
                    "species": species, "antigen.epitope": epitope, "gene": gene,
                    "aa": aa, "pos": pos, "len": length,
                    "v.segm.repr": v, "j.segm.repr": j, "cid": cid, "csz": csz,
                    "count": int(col[k]), "count.bg": int(bg_col[k]), "total.bg": int(total_bg),
                    "count.bg.i": int(coarse_col[k]), "total.bg.i": int(coarse_col.sum()),
                    "level.bg": level, "freq": float(freq[k]), "freq.bg": float(freq_bg[k]),
                    "I": info, "I.norm": info_norm,
                    "height.I": float(freq[k]) * info, "height.I.norm": float(freq[k]) * info_norm,
                })

    if not rows:
        return pl.DataFrame(schema={"cid": pl.Utf8, "pos": pl.Int64, "aa": pl.Utf8})
    return pl.DataFrame(rows).sort("species", "gene", "antigen.epitope", "cid", "len", "pos", "aa")


def information_delta(pwms: pl.DataFrame, shipped: pl.DataFrame) -> pl.DataFrame:
    """Per-cluster ``sum(I)``, ours against the shipped file's -- the sign assertion of section 8.5.

    At a cluster the legacy filter truncated, ours must be **lower**: the shipped ``I`` was computed
    over a distribution missing its rarest residues, and dropping mass from a distribution can only
    make it look more determined than it is. A cluster where ours is higher is a bug in this module,
    not an improvement.
    """
    def per_cid(df: pl.DataFrame, name: str) -> pl.DataFrame:
        return (df.group_by(["cid", "len", "pos"]).agg(pl.first("I"))
                  .group_by("cid").agg(pl.col("I").sum().alias(name),
                                       pl.col("I").len().alias(f"{name}.positions")))

    return (per_cid(pwms, "I.ours")
            .join(per_cid(shipped.with_columns(pl.col("I").cast(pl.Float64),
                                               pl.col("pos").cast(pl.Int64),
                                               pl.col("len").cast(pl.Int64)), "I.shipped"),
                  on="cid", how="inner")
            .with_columns((pl.col("I.ours") - pl.col("I.shipped")).alias("delta"))
            .sort("delta", descending=True))
