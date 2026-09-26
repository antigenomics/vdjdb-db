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


def test_leiden_at_resolution_zero_is_the_connected_components():
    """The two partitions are one knob, not two code paths -- so the sweep is continuous."""
    edges = [(0, 1), (1, 2), (0, 2), (2, 3), (3, 4), (4, 5), (3, 5), (7, 8)]
    def parts(labels):
        out = {}
        for i, lab in enumerate(labels):
            out.setdefault(lab, set()).add(i)
        return {frozenset(v) for v in out.values()}
    assert parts(C._leiden(9, edges, 0.0)) == parts(C._components(9, edges))


def test_leiden_splits_a_percolated_component_and_stays_deterministic():
    """Two triangles joined by one edge: a component, but two motifs. And Leiden is randomised,
    so the same graph must return the same partition every call (rule 7)."""
    edges = [(0, 1), (1, 2), (0, 2), (2, 3), (3, 4), (4, 5), (3, 5)]
    assert len(set(C._components(6, edges))) == 1
    labels = C._leiden(6, edges, 0.3)
    assert len(set(labels)) == 2
    assert labels[:3] == [labels[0]] * 3 and labels[3:] == [labels[3]] * 3
    assert C._leiden(6, edges, 0.3) == labels


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


def test_q_is_zero_at_both_degenerate_ends():
    """The property that makes Q one number instead of purity and retention: shattering and
    percolating are both scored 0, so neither end can be gamed (ROADMAP section 35)."""
    import clustereval
    import numpy as np

    y = np.array(["A"] * 10 + ["B"] * 10)
    def q(k):
        h = float(clustereval.homogeneity_score(y, k))
        p = float(clustereval.parsimony_score(y, k))
        return 0.0 if h + p == 0 else 2 * h * p / (h + p)

    assert q(np.array([0] * 10 + [1] * 10)) == pytest.approx(1.0)   # clusters == epitopes
    assert q(np.arange(20)) == pytest.approx(0.0)                    # every clonotype alone
    assert q(np.zeros(20, dtype=int)) == pytest.approx(0.0)          # one cluster for everything


def test_unclustered_clonotypes_are_singletons_not_dropped():
    """How Q charges for coverage: leaving a clonotype out costs parsimony, so retention does not
    need to be a separate axis."""
    from vdjdb.validate import qscore

    cohort = pl.DataFrame({"antigen.epitope": ["E1", "E1"], "species": ["HomoSapiens"] * 2,
                           "gene": ["TRB"] * 2, "cdr3aa": ["CASSA", "CASSB"],
                           "v.segm": ["V1"] * 2, "j.segm": ["J1"] * 2})
    members = cohort.head(1).with_columns(pl.lit("H.B.E1.1").alias("cid"))
    f = qscore.frame(cohort, members)
    assert f.height == 2
    assert f["cluster"].n_unique() == 2
    assert f["cluster"].str.starts_with("__singleton.").sum() == 1


def test_ball_volume_is_the_substitution_neighbourhood():
    """V_s(L) = sum_j C(L,j) 19^j. The j=0 term is the sequence itself."""
    from vdjdb.validate.noise import ball_volume

    assert ball_volume(14, 0) == 1
    assert ball_volume(14, 1) == 1 + 14 * 19
    assert ball_volume(14, 2) == 1 + 14 * 19 + math.comb(14, 2) * 361
    # The combinatorial bound the measured background rate must never be replaced by: 124x, where
    # the measured inflation of p_hat from scope 1 to scope 2 is 12x (docs/denoising.md section 5.1).
    assert round(ball_volume(14, 2) / ball_volume(14, 1)) == 124


def test_chance_recruitment_rises_with_the_ball_and_the_sample():
    """alpha is monotone in both arguments -- which is why one p-value is not one noise level."""
    from vdjdb.validate.noise import chance_recruitment

    assert chance_recruitment(1000, 1e-6) < chance_recruitment(1000, 1e-5)
    assert chance_recruitment(100, 1e-5) < chance_recruitment(10_000, 1e-5)
    assert chance_recruitment(1, 1e-5) == pytest.approx(0.0)


def test_publicity_control_finds_no_lift_when_there_is_none():
    """The null must sit at 1.0 when clustering and replication are unrelated, and the ratio must
    collapse to ~1 when the association is entirely explained by the stratum."""
    import numpy as np

    from vdjdb.validate.noise import controlled_lift

    rng = np.random.default_rng(0)
    n = 2000
    # replication independent of clustering: raw lift ~1, ratio ~1
    df = pl.DataFrame({"clustered": rng.random(n) < 0.3,
                       "replicated": rng.random(n) < 0.1,
                       "stratum": ["a"] * n})
    r = controlled_lift(df, n_perm=100)
    assert r["lift"] == pytest.approx(1.0, abs=0.25)
    assert r["ratio"] == pytest.approx(1.0, abs=0.25)

    # association driven entirely by the stratum: raw lift is high, the within-stratum null matches
    # it, so the ratio returns to ~1 and the enrichment is correctly attributed to the covariate.
    hot = np.arange(n) < n // 4
    df = pl.DataFrame({"clustered": hot, "replicated": hot,
                       "stratum": np.where(hot, "hi", "lo")})
    r = controlled_lift(df, n_perm=100)
    assert r["lift"] > 3.0
    assert r["ratio"] == pytest.approx(1.0, abs=0.05)


def test_the_guarded_knee_is_invariant_to_how_many_points_the_curve_carries():
    """The defect the guard exists for: a degree-10 fit over 112,983 points returns knee index 1.

    Resampling onto a fixed grid makes the knee a property of the curve's *shape*, so the same
    elbow is found at n = 300 and at n = 300,000 rather than dissolving into fit oscillation.
    """
    from vdjdb.motifs.tcremp import knee

    fracs = []
    for n in (300, 30_000, 300_000):
        x = np.linspace(0, 1, n)
        y = np.where(x < 0.8, x * 0.2, 0.16 + (x - 0.8) * 4.2)   # convex, elbow at 0.80
        k = knee(y, concave=False)
        assert not k["degenerate"]
        fracs.append(k["frac"])
    assert all(abs(f - 0.80) < 0.01 for f in fracs)
    assert max(fracs) - min(fracs) < 0.005      # n changes the answer by less than half a percent


def test_the_guarded_knee_declines_a_curve_that_has_no_knee():
    """A straight line has strength exactly 0 and must be reported as degenerate, not given an index.

    This is the case the reference implementation answered anyway, which is how ``eps`` ended up
    below the data and retention collapsed to a few per cent.
    """
    from vdjdb.motifs.tcremp import knee

    k = knee(np.linspace(1.0, 5.0, 50_000))
    assert k["strength"] == pytest.approx(0.0, abs=1e-9)
    assert k["degenerate"] and k["reason"] == "no-knee"

    assert knee(np.full(100, 3.0))["reason"] == "flat"           # no range at all
    assert knee(np.array([1.0, 2.0]))["reason"] == "too-short"

    # A knee pinned to the bottom of the curve is rejected rather than used.
    spike = np.concatenate([[0.0], np.linspace(5.0, 5.001, 9999)])
    assert knee(spike)["reason"] == "at-floor"


def test_the_guarded_knee_does_not_depend_on_the_units_of_the_curve():
    """Normalising to the unit square is what makes ``strength`` a comparable number across chains."""
    from vdjdb.motifs.tcremp import knee

    x = np.linspace(0, 1, 5000)
    y = np.where(x < 0.6, x * 0.1, 0.06 + (x - 0.6) * 3.0)
    a, b = knee(y, concave=False), knee(y * 1000.0 + 7.0, concave=False)
    assert a["frac"] == pytest.approx(b["frac"])
    assert a["strength"] == pytest.approx(b["strength"])


def test_hdbscan_labels_keeps_the_cluster_labels_contract():
    """Same contract as ``cluster_labels``: -1 is noise, and a label never spans two epitopes."""
    from vdjdb.motifs.tcremp import hdbscan_labels

    rng = np.random.default_rng(0)
    # two well-separated blobs per epitope, two epitopes
    blobs = np.vstack([rng.normal(c, 0.05, (30, 2)) for c in ([0, 0], [5, 5], [0, 0], [5, 5])])
    epitopes = np.array(["A"] * 60 + ["B"] * 60)

    labels = hdbscan_labels(blobs, epitopes, min_cluster_size=5)
    assert labels.min() >= -1
    clustered = labels >= 0
    assert clustered.sum() > 0
    for lab in np.unique(labels[clustered]):
        assert len(set(epitopes[labels == lab])) == 1        # never spans an epitope

    # An epitope with fewer records than min_cluster_size is left entirely unclustered.
    tiny = hdbscan_labels(blobs[:3], np.array(["A"] * 3), min_cluster_size=5)
    assert (tiny == -1).all()


def _members(rows):
    return pl.DataFrame(
        {"species": ["HomoSapiens"] * len(rows), "gene": ["TRB"] * len(rows),
         "cdr3aa": [r[0] for r in rows], "v.segm": ["TRBV1"] * len(rows),
         "j.segm": ["TRBJ1"] * len(rows), "cid": [r[1] for r in rows]})


def test_per_epitope_counts_clonotypes_not_records():
    """A clonotype reported by many papers must not weigh many times in its epitope's retention."""
    from vdjdb.validate.motif_bench import per_epitope

    cohort = pl.DataFrame({           # CASSA reported 3x, CASSB once
        "species": ["HomoSapiens"] * 4, "gene": ["TRB"] * 4,
        "antigen.epitope": ["E1"] * 4,
        "cdr3aa": ["CASSA", "CASSA", "CASSA", "CASSB"],
        "v.segm": ["TRBV1"] * 4, "j.segm": ["TRBJ1"] * 4})
    r = per_epitope(cohort, _members([("CASSA", "c1")]))
    assert r.height == 1
    assert r["clonotypes"][0] == 2 and r["clustered"][0] == 1
    assert r["retention"][0] == pytest.approx(0.5)


def test_per_epitope_percolation_is_one_when_an_epitope_collapses():
    """Percolation is the failure a lift figure cannot see: one cluster holding everything."""
    from vdjdb.validate.motif_bench import per_epitope

    seqs = ["CASS" + a for a in "ABCDEF"]
    cohort = pl.DataFrame({
        "species": ["HomoSapiens"] * 6, "gene": ["TRB"] * 6,
        "antigen.epitope": ["E1"] * 3 + ["E2"] * 3, "cdr3aa": seqs,
        "v.segm": ["TRBV1"] * 6, "j.segm": ["TRBJ1"] * 6})
    # E1 collapses into one cluster; E2 is split into three singletons
    m = _members([(s, "one") for s in seqs[:3]] + [(s, f"s{i}") for i, s in enumerate(seqs[3:])])
    r = per_epitope(cohort, m).sort("antigen.epitope")
    e1, e2 = r.row(0, named=True), r.row(1, named=True)
    assert e1["clusters"] == 1 and e1["percolation"] == pytest.approx(1.0)
    assert e2["clusters"] == 3 and e2["percolation"] == pytest.approx(1 / 3, abs=1e-4)
    assert e2["singleton_clusters"] == 3 and e1["singleton_clusters"] == 0


def test_per_epitope_lift_is_null_where_there_is_no_base_rate():
    """An epitope with no replicated clonotype has no lift -- not a lift of zero.

    Averaging a fabricated zero across such epitopes is how a per-epitope mean lift ends up lower
    than the pooled one for no reason.
    """
    from vdjdb.validate.motif_bench import per_epitope

    seqs = ["CASS" + a for a in "ABCD"]
    cohort = pl.DataFrame({
        "species": ["HomoSapiens"] * 4, "gene": ["TRB"] * 4,
        "antigen.epitope": ["E1", "E1", "E2", "E2"], "cdr3aa": seqs,
        "v.segm": ["TRBV1"] * 4, "j.segm": ["TRBJ1"] * 4})
    rep = cohort.with_columns(pl.Series("replicated", [True, False, False, False]))
    r = per_epitope(cohort, _members([(s, "c1") for s in seqs]), replicated=rep).sort("antigen.epitope")
    e1, e2 = r.row(0, named=True), r.row(1, named=True)
    assert e1["replicated"] == 1 and e1["lift"] == pytest.approx(1.0)   # clustered everything
    assert e2["replicated"] == 0 and e2["lift"] is None                 # no base rate
