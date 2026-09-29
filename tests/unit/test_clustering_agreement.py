"""Does the shipped clustering still give a record its cluster-mates?

Nothing compared the motif tables before this: `vdjdb diff` keys `cluster_members.txt` on `cid`, a
cid carries a position in a sorted list, so the CI comparison names its five legacy members with
`--only` and the motif files were judged by purity and retention alone - statistics about a
clustering rather than a comparison of the tables that ship.
"""
from __future__ import annotations

import polars as pl
import pytest

from vdjdb.compare.clustering import agreement

SCHEMA = {"cid": pl.Utf8, "species": pl.Utf8, "gene": pl.Utf8, "cdr3aa": pl.Utf8,
          "v.segm": pl.Utf8, "j.segm": pl.Utf8}


def members(*pairs: tuple[str, str]) -> pl.DataFrame:
    """``(cid, cdr3aa)`` rows in one species/chain, shaped like `cluster_members.txt`."""
    return pl.DataFrame([(c, "HomoSapiens", "TRB", s, "TRBV1*01", "TRBJ1-1*01")
                         for c, s in pairs], orient="row", schema=SCHEMA)


def test_a_relabelled_clustering_is_the_same_clustering() -> None:
    """The case `vdjdb diff` cannot see: every cid renamed, every cluster-mate kept."""
    a = members(("H.B.E.1", "CA"), ("H.B.E.1", "CB"), ("H.B.E.2", "CC"), ("H.B.E.2", "CD"))
    b = members(("H.B.E.7", "CA"), ("H.B.E.7", "CB"), ("H.B.E.9", "CC"), ("H.B.E.9", "CD"))
    got = agreement(a, b)
    assert got["neighbours_preserved"] == 1.0
    assert got["pairs_preserved"] == 1.0
    assert got["ari"] == 1.0


def test_a_split_cluster_loses_exactly_the_mates_it_lost() -> None:
    a = members(("c", "CA"), ("c", "CB"), ("c", "CC"))
    b = members(("x", "CA"), ("x", "CB"), ("y", "CC"))
    got = agreement(a, b)
    # CA and CB keep one of two mates, CC keeps none of two: (0.5 + 0.5 + 0) / 3.
    assert got["neighbours_preserved"] == pytest.approx(1 / 3)
    assert got["pairs_preserved"] == pytest.approx(1 / 3)   # 1 of 3 pairs


def test_a_clonotype_the_candidate_drops_keeps_no_mates() -> None:
    """Dropping a record is a loss, and the pair-level figures alone would call it agreement."""
    a = members(("c", "CA"), ("c", "CB"), ("c", "CC"))
    b = members(("x", "CA"), ("x", "CB"))
    got = agreement(a, b)
    assert got["clonotypes.shared"] == 2.0
    assert got["neighbours_preserved"] == pytest.approx(1 / 3)
    # On the two clonotypes it did place, the candidate agrees perfectly - which is why the adjusted
    # Rand index over the shared set is not the gate.
    assert got["ari"] == 1.0


def test_clustering_nothing_is_not_agreement() -> None:
    a = members(("c", "CA"), ("c", "CB"), ("c", "CC"))
    empty = members()
    got = agreement(a, empty)
    assert got["neighbours_preserved"] == 0.0
    assert got["pairs_preserved"] == 0.0


def test_reference_singletons_do_not_score() -> None:
    """A cluster of one has no mate to keep, so counting it would reward clustering nothing."""
    a = members(("c", "CA"), ("c", "CB"), ("solo", "CZ"))
    b = members(("x", "CA"), ("x", "CB"))
    assert agreement(a, b)["neighbours_preserved"] == 1.0


def test_a_merged_cluster_keeps_every_mate() -> None:
    """Merging is not a loss of cluster-mates. It is a loss of purity, which `purity` measures."""
    a = members(("c1", "CA"), ("c1", "CB"), ("c2", "CC"), ("c2", "CD"))
    b = members(("one", "CA"), ("one", "CB"), ("one", "CC"), ("one", "CD"))
    got = agreement(a, b)
    assert got["neighbours_preserved"] == 1.0
    assert got["pairs_preserved"] == 1.0
    assert got["ari"] < 1.0        # the partitions are not the same partition


def test_one_giant_cluster_dominates_the_pair_figures_and_not_the_gate() -> None:
    """Why the gate is per record. The released TRB clustering has 19,971 of 36,906 clonotypes in
    one cluster, 94.7 % of the file's pairs, so a pair-weighted score measures that blob alone."""
    big = [("blob", f"C{i}") for i in range(200)]
    small = [("s1", "CX"), ("s1", "CY"), ("s1", "CZ")]
    a = members(*big, *small)
    # The blob is reproduced; the three-member cluster is shattered.
    b = members(*big, ("a", "CX"), ("b", "CY"), ("c", "CZ"))
    got = agreement(a, b)
    assert got["pairs_preserved"] > 0.999          # 19,900 of 19,903 pairs
    assert got["neighbours_preserved"] == pytest.approx(200 / 203)
