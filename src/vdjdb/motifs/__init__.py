"""Motif inference: which CDR3s recur against an epitope more than the repertoire explains.

Three stages, one per module, each a pure function of the one before:

* :mod:`~vdjdb.motifs.tcrnet` -- per-epitope neighbourhood enrichment against a matched background.
* :mod:`~vdjdb.motifs.cluster` -- the Hamming-1 graph over what survives, its components, a layout.
* :mod:`~vdjdb.motifs.pwm` -- a position weight matrix per cluster stratum, read against the
  background repertoire.

:mod:`~vdjdb.motifs.emit` projects the result into the two positionally-parsed legacy files.
:func:`run` is the whole pipeline; ``vdjdb motifs`` is a thin wrapper over it.

Nothing here is stored between builds. The background is fetched (an input) and the control index is
content-addressed by ``seqtree``; everything computed -- enrichment, clusters, PWMs -- is recomputed
every run, and cluster numbers follow cluster content rather than a stored assignment, so they are
stable across releases without anything being remembered (CLAUDE.md hard rule 9).
"""
from __future__ import annotations

from pathlib import Path

import polars as pl

from . import cluster, emit, pwm, tcrnet

__all__ = ["cluster", "emit", "pwm", "run", "tcrnet"]


def run(tables: Path, out: Path, *, q: float = tcrnet.Q_THRESHOLD,
        min_sample: int = tcrnet.MIN_SAMPLE,
        min_cluster: int = cluster.MIN_CLUSTER) -> dict[str, int]:
    """Infer TCRNET motifs from the definitive tables and write both legacy files under ``out``.

    Returns ``{filename: rows}``. The background is loaded once per ``(species, gene)`` present --
    four indices and four frames at most -- and reused across that chain's epitopes.
    """
    chains = pl.read_parquet(tables / "chains.parquet")
    records = pl.read_parquet(tables / "records.parquet")

    enriched = tcrnet.enriched_clonotypes(chains, records, q=q, min_sample=min_sample)
    members = cluster.clusters(enriched, min_cluster=min_cluster)

    pwms = []
    for (species, gene), grp in members.group_by(["species", "gene"], maintain_order=True):
        background = tcrnet.background_frame(species, gene)
        if background is None:
            continue
        pwms.append(pwm.cluster_pwms(grp, background))
    all_pwms = pl.concat(pwms, how="vertical") if pwms else pwm.cluster_pwms(members.head(0),
                                                                            pl.DataFrame())

    return emit.write(emit.cluster_members(members, chains, records),
                      emit.motif_pwms(all_pwms, records), out)
