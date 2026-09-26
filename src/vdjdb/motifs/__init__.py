"""Motif inference: which CDR3s recur against an epitope more than the repertoire explains.

Two methods, run side by side and shipped as two pairs of files.

**TCRNET** asks whether a clonotype has more one-substitution neighbours than a matched background
repertoire accounts for, then clusters what survives on the Hamming-1 graph:
:mod:`~vdjdb.motifs.tcrnet` scores, :mod:`~vdjdb.motifs.cluster` clusters.

**TCREMP** embeds each clonotype as its distance to a fixed prototype panel and clusters by density
in that space, so it groups receptors that are similar in a way an edit distance cannot express:
:mod:`~vdjdb.motifs.tcremp`.

:mod:`~vdjdb.motifs.pwm` turns any cluster into a logo and is shared. :mod:`~vdjdb.motifs.emit`
projects either result into the positionally-parsed legacy files.

Nothing here is stored between builds. The background is fetched (an input) and the control index is
content-addressed by ``seqtree``; enrichment, embeddings, scalers, PCAs, clusters and PWMs are all
recomputed every run, and cluster numbers follow cluster content rather than a stored assignment, so
they are stable across releases without anything being remembered (CLAUDE.md hard rule 9).

**TCRNET is tuned** under the two-stage rule in ``docs/denoising.md`` section 7.1 -- admissible on
Q, purity and precision against the shipped annotation, then ranked on independent-study lift. Its
``TUNED`` carries the scorecard and the one cost it pays. TCREMP is being fitted the same way;
until it is, its defaults are correct and reproducible but not optimised. ROADMAP section 30 lists
what remains open, with the measurement that would settle each.

⚠ **A lift figure is meaningless without its denominator.** Display-selected records contribute no
independently-replicated pairs while filling a quarter of the human TRB cohort, so lift is scored on
the non-display subset and the same clustering reads 1.098 or 5.259 depending only on which cohort
is used (``docs/denoising.md`` section 6.1).
"""
from __future__ import annotations

from pathlib import Path

import polars as pl

from . import cluster, emit, pwm, tcremp, tcrnet

__all__ = ["cluster", "emit", "pwm", "run", "tcremp", "tcrnet"]


def _pwms(members: pl.DataFrame) -> pl.DataFrame:
    """PWMs for every ``(species, gene)`` present, against that chain's background repertoire."""
    parts = []
    for (species, gene), grp in members.group_by(["species", "gene"], maintain_order=True):
        background = tcrnet.background_frame(species, gene)
        if background is not None:
            parts.append(pwm.cluster_pwms(grp, background))
    return pl.concat(parts, how="vertical") if parts else pwm.cluster_pwms(members.head(0),
                                                                          pl.DataFrame())


def run_tcrnet(chains: pl.DataFrame, records: pl.DataFrame, out: Path, *,
               p: float | None = None, min_sample: int = tcrnet.MIN_SAMPLE,
               min_cluster: int | None = None) -> dict[str, int]:
    """TCRNET: enrichment, the two-stage neighbourhood graph, connected components.

    ``p``, the scope, the Leiden ``resolution`` and ``min_cluster`` all default to
    :data:`vdjdb.motifs.tcrnet.TUNED`. It is keyed per chain, but both chains currently hold the
    same configuration -- the per-chain split section 34 introduced did not survive being scored on
    the corrected objective.
    """
    scored = tcrnet.enriched_clonotypes(chains, records, p=p, min_sample=min_sample)
    members = cluster.clusters(scored, min_cluster=min_cluster)
    return emit.write(emit.cluster_members(members, chains, records),
                      emit.motif_pwms(_pwms(members), records), out)


def run_tcremp(chains: pl.DataFrame, records: pl.DataFrame, out: Path, *,
               coef: dict[str, float] | None = None,
               min_cluster: int | None = None) -> dict[str, int]:
    """TCREMP: chunked embedding, a chain-global radius, per-epitope DBSCAN.

    One embedding pass per ``(species, gene)``; the radius is estimated on that chain's pooled
    geometry and then DBSCAN runs per epitope, so it is never re-estimated on an n of 30-300
    (ROADMAP section 8.4).
    """
    cohort = tcremp.cohort(chains, records)
    parts = []
    for (species, gene), grp in cohort.group_by(["species", "gene"], maintain_order=True):
        tuned = tcremp.TUNED.get(gene, {})
        X = tcremp.embed_reduced(grp, species, gene,
                                 n_components=tuned.get("n_components", tcremp.N_COMPONENTS))
        eps = tcremp.chain_eps(X, (coef or tcremp.COEF)[gene])
        parts.append(tcremp.clusters(grp, X, eps, min_cluster=(
            min_cluster if min_cluster is not None else tuned.get("min_cluster", 5))))
    members = pl.concat([p for p in parts if p.height], how="vertical")
    return emit.write(emit.cluster_members(members, chains, records),
                      emit.motif_pwms(_pwms(members), records), out, suffix="_tcremp")


def run(tables: Path, out: Path, *, p: float | None = None,
        min_sample: int = tcrnet.MIN_SAMPLE,
        min_cluster: int | None = None,
        methods: tuple[str, ...] = ("tcrnet", "tcremp")) -> dict[str, int]:
    """Both methods from the definitive tables. Returns ``{filename: rows}``."""
    chains = pl.read_parquet(tables / "chains.parquet")
    records = pl.read_parquet(tables / "records.parquet")

    written: dict[str, int] = {}
    if "tcrnet" in methods:
        written |= run_tcrnet(chains, records, out, p=p, min_sample=min_sample,
                              min_cluster=min_cluster)
    if "tcremp" in methods:
        written |= run_tcremp(chains, records, out, min_cluster=min_cluster)
    return written
