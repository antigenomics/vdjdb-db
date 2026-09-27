"""The vendored clustering metrics, and the wrapper the acceptance rule reads them through.

`validate/metrics_lib.py` is copied verbatim from the benchmark this build has to stay comparable
with, so it must not be edited. It had 0 % coverage across 97 statements while producing the
purity, retention and F1 that `docs/clustering.md` ranks 252 configurations on. These pin its
behaviour on hand-built assignments where the answer can be worked out by hand, so an upstream
change or an accidental edit shows up as a failing test rather than as a moved scorecard.

`motif_bench.score` is the wrapper; it is checked against the same frames.
"""
from __future__ import annotations

import pytest

# `metrics_lib` is vendored verbatim and imports sklearn at module level, so the skip has to be
# here rather than in it. pandas arrives with the same extra.
pytest.importorskip("sklearn", reason="needs the `motifs` extra")

import pandas as pd

from vdjdb.validate.metrics_lib import (
    binominal_test,
    count_clstr_purity,
    metrics_from_assignments,
    per_epitope_metrics,
)


def assign(pairs):
    """`[(epitope, cluster), ...]` as the frame the metrics take."""
    return pd.DataFrame(pairs, columns=["antigen.epitope", "cluster"])


#: Two clusters, each entirely one epitope, nothing unclustered.
PERFECT = assign([("AAA", 0), ("AAA", 0), ("AAA", 0), ("BBB", 1), ("BBB", 1), ("BBB", 1)])

#: The same records with every cluster half one epitope and half the other.
MIXED = assign([("AAA", 0), ("AAA", 0), ("BBB", 0), ("BBB", 1), ("BBB", 1), ("AAA", 1)])


def test_a_perfect_clustering_scores_one_on_every_axis():
    m, _, _ = metrics_from_assignments(PERFECT)
    assert m["purity"] == 1.0
    assert m["retention"] == 1.0
    assert m["precision"] == 1.0
    assert m["recall"] == 1.0
    assert m["f1"] == 1.0
    assert m["n_records"] == 6
    assert m["n_epitopes"] == 2
    assert m["n_clusters_kept"] == 2


def test_purity_is_matched_over_total_across_kept_clusters():
    """Each cluster is 2 of 3 its majority epitope, so purity is 4/6."""
    m, _, _ = metrics_from_assignments(MIXED)
    assert m["purity"] == pytest.approx(4 / 6, abs=5e-5)


def test_noise_is_excluded_from_clusters_and_lowers_retention():
    """Cluster -1 is the noise label: never a cluster, and its records count against retention."""
    df = assign([("AAA", 0), ("AAA", 0), ("AAA", -1), ("AAA", -1)])
    m, _, _ = metrics_from_assignments(df)
    assert m["retention"] == 0.5
    assert m["n_clusters_kept"] == 1


def test_a_singleton_is_not_a_cluster():
    """`total_cluster > 1` is the rule, so a cluster of one never counts."""
    df = assign([("AAA", 0), ("AAA", 0), ("BBB", 7)])
    _, binom, _ = metrics_from_assignments(df)
    kept = binom[binom["is_cluster"] == 1]["cluster"].tolist()
    assert kept == [0]


def test_nothing_clustered_gives_zeros_rather_than_raising():
    """The degenerate case an admissibility sweep hits constantly."""
    df = assign([("AAA", -1), ("BBB", -1)])
    m, _, _ = metrics_from_assignments(df)
    assert m["purity"] == 0.0 and m["retention"] == 0.0 and m["f1"] == 0.0
    assert m["n_records"] == 2


def test_count_clstr_purity_returns_none_when_no_cluster_survives():
    """A None, not a zero, which is why the wrapper above special-cases the degenerate branch."""
    binom = binominal_test(assign([("AAA", -1), ("BBB", -1)]), "cluster", "antigen.epitope",
                           compute_pvalue=False)
    assert count_clstr_purity(binom) is None


def test_the_enriched_flag_needs_the_fraction_threshold():
    """`is_cluster` says a cluster exists; `enriched_clstr` says it is mostly one epitope."""
    binom = binominal_test(MIXED, "cluster", "antigen.epitope", compute_pvalue=False)
    assert set(binom["is_cluster"]) == {1}
    assert set(binom["enriched_clstr"]) == {0}, "2 of 3 is below the 0.7 default"


def test_the_threshold_is_a_parameter():
    binom = binominal_test(MIXED, "cluster", "antigen.epitope", threshold=0.6,
                           compute_pvalue=False)
    assert set(binom["enriched_clstr"]) == {1}


def test_p_values_are_computed_only_when_asked():
    off = binominal_test(PERFECT, "cluster", "antigen.epitope", compute_pvalue=False)
    on = binominal_test(PERFECT, "cluster", "antigen.epitope", compute_pvalue=True)
    assert set(off["p_value"]) == {0.0}
    assert all(0.0 < p <= 1.0 for p in on["p_value"])


def test_per_epitope_counts_unclustered_records_as_false_negatives():
    """Recall per epitope must see records the clustering dropped, or it reads far too high."""
    df = assign([("AAA", 0), ("AAA", 0), ("AAA", -1), ("AAA", -1)])
    _, _, data = metrics_from_assignments(df)
    out = per_epitope_metrics(data).set_index("antigen.epitope")
    assert out.loc["AAA", "precision"] == 1.0
    assert out.loc["AAA", "recall"] == 0.5
    assert out.loc["AAA", "retention"] == 0.5
    assert out.loc["AAA", "n_records"] == 4


def test_per_epitope_returns_one_row_per_epitope_in_sorted_order():
    _, _, data = metrics_from_assignments(MIXED)
    out = per_epitope_metrics(data)
    assert out["antigen.epitope"].tolist() == ["AAA", "BBB"]
