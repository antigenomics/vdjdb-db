"""Motif inference: the statistic, the graph, and the PWM cascade that stopped deleting letters."""
from __future__ import annotations

import math

import numpy as np
import polars as pl
import pytest

from vdjdb.motifs import cluster as C
from vdjdb.motifs import emit as E
from vdjdb.motifs import pwm as P
from vdjdb.motifs import tcrnet as T


def test_the_legacy_statistic_never_returns_zero_where_the_shipped_one_does():
    """The whole reason ``tcrnet``'s own p-value is unusable (ROADMAP section 8.1).

    ``poisson.sf(k - 1, 0.0)`` is exactly 0.0 for every k >= 1, so a clonotype with neighbours and
    no background neighbour is handed certainty. The pseudocount is what removes that.
    """
    from scipy.stats import poisson
    assert poisson.sf(4, 0.0) == 0.0            # the defect, asserted so it cannot silently change

    df = pl.DataFrame({"d": [5], "nc": [0]})
    p = df.select(T.legacy_pvalue(pl.col("d"), pl.col("nc"), n_sample=1000, m_control=1_000_000))
    assert 0.0 < p.item() < 1.0


def test_the_statistic_is_the_groovy_formula():
    """Bit-faithful to ``DegreeStatisticsAnnotator.computePValue``, computed here by hand."""
    from scipy.stats import binom
    d, nc, n, m = 7, 3, 500, 1_000_000
    p = (nc + 1) / (m + 1)
    want = binom.sf(d - 1, n, p) / (1 - (1 - p) ** n)
    got = pl.DataFrame({"d": [d], "nc": [nc]}).select(
        T.legacy_pvalue(pl.col("d"), pl.col("nc"), n, m)).item()
    assert got == pytest.approx(want, rel=1e-12)


def test_more_background_neighbours_is_never_more_significant():
    """Monotone in the background count, which a pseudocount must not break."""
    df = pl.DataFrame({"d": [6, 6, 6], "nc": [0, 5, 50]})
    p = df.select(T.legacy_pvalue(pl.col("d"), pl.col("nc"), 800, 1_000_000)).to_series()
    assert p[0] < p[1] < p[2]


def test_bh_is_monotone_and_bounded():
    p = pl.DataFrame({"p": [0.001, 0.01, 0.2, 0.5, 0.9]})
    q = p.select(T._bh(pl.col("p"))).to_series().to_list()
    assert q == sorted(q)
    assert max(q) <= 1.0
    assert q[0] == pytest.approx(0.005)          # 0.001 * 5 / 1


def test_edges_are_the_hamming_one_ball_without_self_or_mirror():
    seqs = ["CASSF", "CASSY", "CAWWW", "CASSW"]
    edges = C._edges(seqs, "1,0,0,1")
    assert all(i < j for i, j in edges)
    assert set(edges) == {(0, 1), (0, 3), (1, 3)}    # CAWWW is two substitutions from every other


def test_components_label_the_connected_pieces():
    assert len(set(C._components(5, [(0, 1), (1, 2), (3, 4)]))) == 2
    assert len(set(C._components(3, []))) == 3


def test_a_neighbour_of_an_enriched_clonotype_joins_the_graph():
    """The legacy Rmd's two-stage construction, and phase 10's whole coverage gap.

    `CWWWWW` is not enriched but sits one substitution from `CWWWWA`, which is -- so it belongs to
    the motif. `CYYYYY` is neither and stays out.
    """
    seqs = ["CWWWWA", "CWWWWC", "CWWWWD", "CWWWWE", "CWWWWF", "CWWWWW", "CYYYYY"]
    enriched = [0, 1, 2, 3, 4]                       # not 5, not 6
    keep = C._recruited(seqs, enriched, "1,0,0,1")
    assert keep == [0, 1, 2, 3, 4, 5]
    assert keep == sorted(keep)                      # order must not follow hit order


def test_only_the_enriched_seed_the_graph():
    """A clonotype two substitutions from every enriched one is not recruited."""
    seqs = ["CWWWWA", "CWWWWC", "CYYYYY"]
    assert C._recruited(seqs, [0, 1], "1,0,0,1") == [0, 1]


def test_cluster_numbering_follows_content_not_vertex_order():
    """A counter would renumber on every reshuffle and break bookmarked motif URLs."""
    base = pl.DataFrame({
        "species": ["HomoSapiens"] * 12, "gene": ["TRB"] * 12, "antigen.epitope": ["EEE"] * 12,
        # six mutually-adjacent, then five, so the numbering has something to order
        "junction_aa": ["CASSA", "CASSC", "CASSD", "CASSE", "CASSF", "CASSG",
                        "CYYYA", "CYYYC", "CYYYD", "CYYYE", "CYYYF", "CWWWW"],
        "v_call": ["TRBV1*01"] * 12, "j_call": ["TRBJ1*01"] * 12,
        "enriched": [True] * 12,
    })
    a = C.clusters(base)
    b = C.clusters(base.sample(fraction=1.0, shuffle=True, seed=7))
    key = ["junction_aa", "cid"]
    assert a.sort("junction_aa").select(key).equals(b.sort("junction_aa").select(key))
    assert a.filter(pl.col("junction_aa") == "CWWWW").height == 0     # below MIN_CLUSTER


#: `H.B.ALSKGVHFV.1` position 8 of the 2026-06-03 release: counts N1 G1 D3 S4 over csz 9, with the
#: `freq.bg` the file ships beside them. Both information columns are pinned to it.
_SHIPPED_POSITION = (("N", 1, 0.0194423126119212), ("G", 1, 0.252408970751258),
                     ("D", 3, 0.0237912509593246), ("S", 4, 0.143898695318496))


def _shipped_column() -> tuple[np.ndarray, np.ndarray]:
    freq, bg = np.zeros(20), np.zeros(20)
    for aa, c, q in _SHIPPED_POSITION:
        freq[P.ALPHABET.index(aa)] = c
        bg[P.ALPHABET.index(aa)] = q
    return freq / freq.sum(), bg


def test_information_reproduces_the_shipped_value():
    freq, _ = _shipped_column()
    assert P._information(freq) == pytest.approx(0.594459870571867, rel=1e-12)


def test_normalised_information_is_the_halved_cross_entropy_the_file_ships():
    """Not ``I`` minus the background's own information -- that is the natural guess and is wrong.

    `vdjdb-motifs/scripts/compute_motif_pwms.py`: ``-sum(p log q) / log 20 / 2``.
    """
    freq, bg = _shipped_column()
    assert P._information_norm(freq, bg) == pytest.approx(0.450398191673246, rel=1e-12)
    assert P._information_norm(freq, bg) != pytest.approx(P._information(freq) - P._information(bg))


def test_information_is_zero_for_uniform_and_one_for_determined():
    assert P._information(np.full(20, 1 / 20)) == pytest.approx(0.0, abs=1e-12)
    one = np.zeros(20)
    one[3] = 1.0
    assert P._information(one) == pytest.approx(1.0)


def test_counts_are_position_by_residue():
    c = P.counts(["CASSF", "CASSY", "CAWWY"], 5)
    assert c.shape == (5, 20)
    assert c[0][P.ALPHABET.index("C")] == 3
    assert c[4][P.ALPHABET.index("Y")] == 2


def _tiny_members() -> pl.DataFrame:
    return pl.DataFrame({
        "species": ["HomoSapiens"] * 5, "gene": ["TRB"] * 5, "antigen.epitope": ["EEE"] * 5,
        "junction_aa": ["CASSF", "CASSY", "CASSW", "CASDF", "CASEF"],
        "v_call": ["TRBV1*01"] * 5, "j_call": ["TRBJ1*01"] * 5,
        "cid": ["H.B.EEE.1"] * 5, "csz": [5] * 5,
        "v.segm.repr": ["TRBV1*01"] * 5, "j.segm.repr": ["TRBJ1*01"] * 5,
    })


def test_every_position_keeps_all_its_mass():
    """The bug the cascade exists to remove: shipped ``freq`` sums to as little as 0.023."""
    bg = pl.DataFrame({"cdr3aa": ["CASSF", "CASSY", "CAAAA"], "v": ["TRBV1"] * 3,
                       "j": ["TRBJ1"] * 3})
    got = P.cluster_pwms(_tiny_members(), bg)
    mass = got.group_by(["cid", "len", "pos"]).agg(pl.col("freq").sum())
    assert mass["freq"].to_numpy() == pytest.approx(1.0)


def test_a_residue_the_background_never_saw_still_gets_a_letter():
    """``W`` at position 4 is absent from the background; the legacy filter deleted it."""
    bg = pl.DataFrame({"cdr3aa": ["CASSF", "CASSY"], "v": ["TRBV1"] * 2, "j": ["TRBJ1"] * 2})
    got = P.cluster_pwms(_tiny_members(), bg)
    w = got.filter((pl.col("pos") == 4) & (pl.col("aa") == "W"))
    assert w.height == 1
    assert w["count.bg"].item() == 0
    assert w["freq.bg"].item() > 0.0          # the pseudocount: rare, never impossible


def test_the_cascade_falls_back_and_says_so():
    """No ``(v, j, len)`` stratum at all -> the ``len`` level, recorded in ``level.bg``."""
    bg = pl.DataFrame({"cdr3aa": ["CQQQF", "CQQQY"], "v": ["TRBV9"] * 2, "j": ["TRBJ9"] * 2})
    got = P.cluster_pwms(_tiny_members(), bg)
    assert set(got["level.bg"]) == {"len"}


def test_the_legacy_column_orders_are_the_contract():
    """``Motifs.scala`` types both files positionally with no header check (hard rule 1)."""
    assert len(E.MEMBER_COLUMNS) == 19
    assert len(E.PWM_COLUMNS) == 27
    assert E.MEMBER_COLUMNS[:3] == ("species", "antigen.epitope", "antigen.gene")
    assert E.PWM_COLUMNS[3:6] == ("aa", "pos", "len")


def test_information_is_lower_once_a_column_keeps_its_tail():
    """ROADMAP section 8.5's sign, as an inequality on one column.

    Dropping a residue from a distribution can only make the rest look more determined, so the
    legacy's filtered ``I`` is an over-estimate and ours must sit below it.
    """
    full = np.array([4, 3, 2, 1] + [0] * 16, dtype=float)
    truncated = np.array([4, 3, 2, 0] + [0] * 16, dtype=float)
    assert P._information(full / full.sum()) < P._information(truncated / truncated.sum())
    assert math.isfinite(P._information(full / full.sum()))
