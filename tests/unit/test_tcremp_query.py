"""Queries use a reference fit and its retained DBSCAN cores without changing either."""
from __future__ import annotations

import numpy as np
import polars as pl
import pytest

from vdjdb.motifs import tcremp as TE


def rows(sequences, epitopes=None):
    n = len(sequences)
    return pl.DataFrame({"species": ["HomoSapiens"] * n, "gene": ["TRB"] * n,
                         "antigen.epitope": epitopes or ["PEP"] * n,
                         "junction_aa": sequences, "v_call": ["TRBV5-1*01"] * n,
                         "j_call": ["TRBJ1-1*01"] * n,
                         "clonotype_id": [str(i) for i in range(n)], "duplicate_count": [1] * n})


def test_frozen_embedding_preserves_original_fit_and_batches_queries(monkeypatch):
    from mir.embedding import TCREmp
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    calls = []
    instances = []

    class Prototype:
        def embed(self, frame):
            calls.append(frame.select(TE.INPUT_KEY))
            a = np.array([sum(map(ord, s)) for s in frame["junction_aa"]], dtype=np.float32)
            return np.column_stack((a, a * a, np.sin(a), np.cos(a))).astype(np.float32)

    def construct(*args):
        model = Prototype()
        instances.append(model)
        return model

    monkeypatch.setattr(TCREmp, "from_defaults", construct)
    reference = rows(["CASAF", "CASCF", "CASDF", "CASEF", "CASFF", "CASGF"])
    # The existing operation order, independently of the new wrapper.
    unique, inverse = TE.receptor_inputs(reference)
    sample = TE.fit_indices(unique, size=4)
    fit = Prototype().embed(unique[sample])
    scaler = StandardScaler().fit(fit)
    pca = PCA(n_components=2, random_state=TE.SEED).fit(scaler.transform(fit))
    want = np.vstack([pca.transform(scaler.transform(Prototype().embed(unique[i:i + 2])))
                      for i in range(0, unique.height, 2)])[inverse]
    calls.clear()
    X, frozen = TE.fit_embedding(reference, "HomoSapiens", "TRB", fit_sample=4, n_components=2, chunk=2)
    assert np.array_equal(X, want)
    assert len(instances) == 1 and frozen.prototype is instances[0]
    legacy = TE.embed_reduced(reference, "HomoSapiens", "TRB", fit_sample=4, n_components=2, chunk=2)
    assert np.array_equal(legacy, X)
    before = (frozen.scaler.mean_.copy(), frozen.pca.components_.copy())
    monkeypatch.setattr(frozen.scaler, "fit", lambda *a, **k: pytest.fail("query refitted scaler"))
    monkeypatch.setattr(frozen.pca, "fit", lambda *a, **k: pytest.fail("query refitted PCA"))
    query = pl.concat([reference.head(3), reference.head(1)]).reverse()
    calls.clear()
    count = len(instances)
    projected = frozen.transform(query, chunk=2)
    assert len(instances) == count and sum(c.height for c in calls) == 3
    np.testing.assert_allclose(projected, np.vstack([X[:3], X[:1]])[::-1], rtol=1e-6, atol=1e-6)
    assert np.array_equal(frozen.transform(query.reverse(), chunk=2)[::-1], projected)
    missing = reference.head(1).with_columns(pl.lit(None, dtype=pl.String).alias("v_call"))
    empty = missing.with_columns(pl.lit("").alias("v_call"))
    assert np.array_equal(frozen.transform(missing), frozen.transform(empty))
    assert frozen.transform(reference.head(0)).shape == (0, 2)
    assert np.array_equal(before[0], frozen.scaler.mean_)
    assert np.array_equal(before[1], frozen.pca.components_)
    for junction in (None, ""):
        bad = reference.head(1).with_columns(pl.lit(junction, dtype=pl.String).alias("junction_aa"))
        with monkeypatch.context() as guard:
            guard.setattr(TE, "receptor_inputs",
                          lambda *a: pytest.fail("invalid junction reached input join"))
            with pytest.raises(ValueError, match="junction_aa"):
                frozen.transform(bad)


def test_core_queries_preserve_reference_and_cannot_expand_through_a_border(monkeypatch):
    from sklearn.cluster import DBSCAN

    reference = rows(["CASAF", "CASCF", "CASDF", "CASEF", "CASFF"])
    X = np.column_stack(([0., .1, .2, .3, 1.1], np.zeros(5)))
    labels = TE.cluster_labels(X, reference["antigen.epitope"].to_numpy(), 1., min_samples=4)
    want = TE.clusters(reference, X, labels=labels, min_cluster=2)
    fitted = TE.fit_clusters(reference, X, 1., min_samples=4, min_cluster=2)
    assert fitted.members.equals(want) and np.array_equal(fitted.labels, labels)
    dropped = TE.fit_clusters(reference, X, 1., min_samples=4, min_cluster=6)
    assert dropped.members.height == 0
    assert not dropped.assign(reference, X)["geometric_support"].any()
    assert dropped.assign(reference, X)["eligible"].all()
    snapshot = fitted.members.clone()
    monkeypatch.setattr(DBSCAN, "fit", lambda *a, **k: pytest.fail("query reclustered reference"))
    query = rows(["CASGF", "CASSGGF", "CASHF", "CASGF"], ["PEP", "PEP", "OTHER", "PEP"])
    Q = np.array([[2., 0.], [.15, 0.], [.15, 0.], [2., 0.]])
    got = fitted.assign(query, Q)
    assert got["geometric_support"].to_list() == [False, True, False, False]
    assert got["eligible"].to_list() == [True, True, False, True]
    assert got["cid"].to_list() == ["", "", "", ""]  # unseen length cannot invent a legacy cid
    assert got["matched_clusters"].to_list() == [0, 1, 0, 0]
    assert not got["ambiguous"].any()
    restored = fitted.assign(reference, X)
    assert restored["cid"].to_list() == want["cid"].to_list()
    assert fitted.members.equals(snapshot)
    assert fitted.assign(query.reverse(), Q[::-1]).reverse().equals(got)
    assert fitted.assign(query.head(0), Q[:0]).height == 0
    for junction in (None, ""):
        bad = query.head(1).with_columns(pl.lit(junction, dtype=pl.String).alias("junction_aa"))
        with pytest.raises(ValueError, match="junction_aa"):
            fitted.assign(bad, Q[:1])


def test_ambiguous_query_keeps_two_frozen_clusters_and_searches_only_its_peptide(monkeypatch):
    from sklearn.neighbors import NearestNeighbors

    reference = rows(["CASAF", "CASCF", "CASDF", "CASEF", "CASFF", "CASGF", "CASHF", "CASIF",
                      "CASKF", "CASLF", "CASMF", "CASNF"], ["PEP"] * 8 + ["OTHER"] * 4)
    X = np.column_stack(([0., .1, .2, .3, 2., 2.1, 2.2, 2.3, 0., .1, .2, .3], np.zeros(12)))
    frozen = TE.fit_clusters(reference, X, 1., min_samples=4, min_cluster=2)
    before = frozen.labels.copy()
    sizes = []
    original = NearestNeighbors.radius_neighbors_graph

    def searched(self, query, *args, **kwargs):
        sizes.append(self.n_samples_fit_)
        return original(self, query, *args, **kwargs)

    monkeypatch.setattr(NearestNeighbors, "radius_neighbors_graph", searched)
    got = frozen.assign(rows(["CASPF", "CASPF"]), np.array([[1.15, 0.], [1.15, 0.]]))
    assert got["matched_clusters"].to_list() == [2, 2]
    assert got["ambiguous"].to_list() == [True, True]
    assert got["geometric_support"].all() and got["cid"].n_unique() == 1
    assert np.array_equal(before, frozen.labels)
    assert sizes == [8]  # deduplicated query; no cores belonging to OTHER
