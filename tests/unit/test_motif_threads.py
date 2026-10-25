"""Neighbour queries own one OpenMP allocation and restore the caller's pools."""
import numpy as np
import pytest
from threadpoolctl import threadpool_info, threadpool_limits

from vdjdb.motifs import tcremp


@pytest.mark.parametrize('operation', ['radius', 'dbscan'])
def test_neighbour_queries_own_allocation_and_restore_pools(monkeypatch, operation):
    import sklearn.cluster
    import sklearn.neighbors
    import vdjtools.cores

    monkeypatch.setattr(vdjtools.cores, 'kernel_threads', lambda threads=0: 2)

    class Neighbours:
        def __init__(self, **kwargs):
            pass

        def fit(self, X):
            pools = threadpool_info()
            assert all(p['num_threads'] == 2 for p in pools if p['user_api'] == 'openmp')
            assert all(p['num_threads'] == 1 for p in pools if p['user_api'] == 'blas')
            self.labels_ = np.zeros(len(X), dtype=int)
            self.core_sample_indices_ = np.arange(len(X))
            return self

        def kneighbors(self, X):
            return np.ones((len(X), 2)), np.zeros((len(X), 2), dtype=int)

    monkeypatch.setattr(sklearn.neighbors, 'NearestNeighbors', Neighbours)
    monkeypatch.setattr(sklearn.cluster, 'DBSCAN', Neighbours)
    with threadpool_limits(limits=1):
        X = np.zeros((2, 3))
        if operation == 'radius':
            assert tcremp.mean_nn_distance(X) == 1
        else:
            assert tcremp.cluster_labels(X, np.array(['p', 'p']), 1).tolist() == [0, 0]
        assert all(p['num_threads'] == 1 for p in threadpool_info())
