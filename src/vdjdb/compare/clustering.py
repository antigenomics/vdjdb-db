"""Do the shipped motif tables still put the same records together?

``vdjdb diff`` keys ``cluster_members.txt`` on ``cid``, and a cid is
``<species>.<chain>.<epitope>.<n>`` where ``n`` is a position in a sorted list - so renumbering one
cluster relabels every cluster after it and the row comparison reports the whole file as replaced.
That is why the CI comparison names its five members with ``--only`` and never the motif files: the
instrument had no answer for them, and the only thing measured about the clustering was
``vdjdb motif-metrics``, which scores purity, retention and F1 on a cohort. Those are statistics
about a clustering, not a comparison of the tables that ship.

The question a cid cannot answer is whether the **partition** survived. What a consumer of
``cluster_members.txt`` depends on is a record's cluster-mates: `vdjdb-web` lists a record's
neighbours, so a record keeping its neighbours under a new cid is no change to anyone, and a record
losing them is a change to everyone.

**The headline is per record, not per pair, and that is a measurement rather than a preference.**
The released TRB clustering puts **19,971 of its 36,906 clonotypes in one cluster**,
``H.B.SLLMWITQV.1``, and that single cluster holds **199,410,435 of the 210,634,576 co-clustered
pairs in the file, 94.7 %**. Any pair-weighted agreement is therefore a measurement of whether that
one blob was reproduced: the do-nothing partition, which clusters every clonotype of an epitope
together and is the bar a method has to beat before its numbers mean anything, scores 0.9991
pair-preservation and 0.9939 adjusted Rand against the release on TRB. So the pair figures are
recorded and the gate is :data:`neighbours_preserved`, which averages over clonotypes and gives the
19,971-member cluster the weight of its 19,971 records rather than of its 199 million pairs.

Everything here is contingency-table arithmetic in one polars pass - no pair is ever materialised,
which matters at 199 million of them.
"""
from __future__ import annotations

import polars as pl

from ..validate.motif_bench import KEY, members_map


def _pairs(df: pl.DataFrame, by: str | list[str]) -> int:
    """``sum C(size, 2)`` over the groups of ``by`` - the pairs those groups place together."""
    if not df.height:
        return 0
    n = df.group_by(by).len()
    return int(n.select((pl.col("len") * (pl.col("len") - 1) // 2).sum()).item() or 0)


def agreement(reference: pl.DataFrame, candidate: pl.DataFrame) -> dict[str, float]:
    """How much of the reference clustering the candidate still reports.

    Both sides go through :func:`vdjdb.validate.motif_bench.members_map`, so a clonotype in two
    clusters is resolved the same way the metrics resolve it and the two instruments cannot disagree
    about which cluster a clonotype is in.
    """
    ref, cand = members_map(reference), members_map(candidate)
    joined = ref.join(cand.rename({"cluster": "cluster.candidate"}), on=list(KEY), how="left")

    # Per reference-clustered clonotype: its reference cluster size, and how many of that cluster's
    # members the candidate keeps with it. A clonotype the candidate does not place at all keeps
    # none. Singletons are excluded - they have no cluster-mate to lose, and counting them as a
    # perfect score would let a candidate that clusters nothing read as agreement.
    sized = joined.with_columns(
        pl.len().over("cluster").alias("ref.size"),
        pl.when(pl.col("cluster.candidate").is_null()).then(pl.lit(0, dtype=pl.UInt32))
          .otherwise(pl.len().over("cluster", "cluster.candidate")).alias("kept"))
    with_mates = sized.filter(pl.col("ref.size") > 1)
    neighbours = with_mates.select(
        ((pl.when(pl.col("kept") > 0).then(pl.col("kept") - 1).otherwise(0))
         / (pl.col("ref.size") - 1)).mean()).item() if with_mates.height else 0.0

    shared = joined.filter(pl.col("cluster.candidate").is_not_null())
    ref_all = _pairs(ref, "cluster")
    ref_shared = _pairs(shared, "cluster")
    cand_shared = _pairs(shared, "cluster.candidate")
    together = _pairs(shared, ["cluster", "cluster.candidate"])

    n = shared.height
    total = n * (n - 1) // 2
    # Expected agreement under the permutation null of the adjusted Rand index. Not a control
    # (`CLAUDE.md` 0c) - it is the definition of the index, not a hypothesis test.
    expected = (ref_shared * cand_shared / total) if total else 0.0
    ceiling = (ref_shared + cand_shared) / 2
    return {
        "neighbours_preserved": float(neighbours or 0.0),
        "clonotypes.reference": float(ref.height),
        "clonotypes.candidate": float(cand.height),
        "clonotypes.shared": float(n),
        "pairs.reference": float(ref_all),
        "pairs.together": float(together),
        "pairs_preserved": (together / ref_all) if ref_all else 0.0,
        "ari": ((together - expected) / (ceiling - expected)) if ceiling != expected else 1.0,
    }
