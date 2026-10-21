"""Peptide duplication must not reweight the shared TCREmp geometry."""
from __future__ import annotations

import os
import subprocess
import sys

import numpy as np
import polars as pl
import pytest

from vdjdb.motifs import tcremp as TE


def test_fit_sample_is_independent_of_process_hash_seed():
    code = """
import polars as pl
from vdjdb.motifs.tcremp import fit_indices
inputs = pl.DataFrame({'junction_aa': [f'CASS{i}F' for i in range(60)],
                       'v_call': ['TRBV5-1*01'] * 60, 'j_call': ['TRBJ1-1*01'] * 60})
print(inputs[fit_indices(inputs, size=12)].write_csv())
"""
    outputs = [subprocess.check_output([sys.executable, "-c", code],
               env={**os.environ, "PYTHONHASHSEED": seed}, timeout=15)
               for seed in ("1", "987")]
    assert outputs[0] == outputs[1]


def test_fit_sample_follows_receptor_identity_not_row_position():
    inputs = pl.DataFrame({
        "junction_aa": [f"CASS{i}F" for i in range(60)],
        "v_call": ["TRBV5-1*01"] * 60, "j_call": ["TRBJ1-1*01"] * 60,
    })
    chosen = inputs[TE.fit_indices(inputs, size=12)]
    shuffled = inputs.reverse()
    assert chosen.equals(shuffled[TE.fit_indices(shuffled, size=12)])
    # Removing receptors outside the sample must not replace the selected receptors.
    reduced = inputs.filter(pl.col("junction_aa").is_in(
        chosen["junction_aa"].to_list() + inputs["junction_aa"].head(6).to_list()))
    assert chosen.equals(reduced[TE.fit_indices(reduced, size=12)])
    assert chosen.height == chosen.unique().height == 12
    assert not chosen.equals(inputs[TE.fit_indices(inputs, size=12, seed=TE.SEED + 1)])


def test_radius_counts_unique_receptors_but_keeps_distinct_equal_vectors():
    inputs = pl.DataFrame({
        "junction_aa": ["CASAF", "CASAF", "CASAF", "CASCF", "CASDF"],
        "v_call": ["TRBV5-1*01"] * 5, "j_call": ["TRBJ1-1*01"] * 5,
    })
    X = np.array([[0.0], [0.0], [0.0], [0.0], [3.0]])
    # Three receptors at 0, 0, 3 have nearest-neighbour distances 0, 0, 3.
    assert TE.mean_nn_distance(X, clonotypes=inputs) == pytest.approx(1.0)
    assert TE.chain_eps(X, 1.5, clonotypes=inputs) == 1.5


def test_embedding_is_unchanged_by_peptide_duplication_and_input_order():
    inputs = pl.DataFrame({
        "junction_aa": ["CASSAGQETQYF", "CASSCGQETQYF", "CASSDGQETQYF",
                        "CASSEGQETQYF", "CASSFGQETQYF", "CASSGGQETQYF",
                        "CASSHGQETQYF", "CASSIGQETQYF", "CASSKGQETQYF",
                        "CASSLGQETQYF", "CASSMGQETQYF", "CASSNGQETQYF"],
        "v_call": ["TRBV5-1*01"] * 12, "j_call": ["TRBJ1-1*01"] * 12,
    })
    X = TE.embed_reduced(inputs, "HomoSapiens", "TRB", fit_sample=8, n_components=2, chunk=3)
    duplicated = pl.concat([inputs, inputs.head(8)]).reverse()
    Y = TE.embed_reduced(duplicated, "HomoSapiens", "TRB", fit_sample=8,
                         n_components=2, chunk=3)
    assert np.array_equal(Y, np.concatenate([X, X[:8]])[::-1])
    assert TE.chain_eps(X, 1.55, clonotypes=inputs) == TE.chain_eps(
        Y, 1.55, clonotypes=duplicated)


def test_pca_geometry_is_independent_of_numerical_thread_allocation():
    """Randomized PCA must not change coordinates with the caller's BLAS budget."""
    code = """
import hashlib
import numpy as np
import polars as pl
from mir.embedding import TCREmp
from vdjdb.motifs import tcremp as TE

matrix = np.random.default_rng(TE.SEED).normal(size=(1024, 6000)).astype(np.float32)
class Prototype:
    def embed(self, frame):
        return matrix[[int(s[1:-1]) for s in frame['junction_aa']]]
TCREmp.from_defaults = lambda *args: Prototype()
inputs = pl.DataFrame({'junction_aa': [f'C{i}F' for i in range(len(matrix))],
                       'v_call': ['TRBV5-1*01'] * len(matrix),
                       'j_call': ['TRBJ1-1*01'] * len(matrix)})
X, frozen = TE.fit_embedding(inputs, 'HomoSapiens', 'TRB', fit_sample=1024, chunk=512)
assert frozen.pca._fit_svd_solver == 'randomized'
assert frozen.pca.components_.dtype == np.float64
print(hashlib.sha256(X.tobytes()).hexdigest())
"""
    outputs = [subprocess.check_output([sys.executable, '-c', code], env={
        **os.environ, 'OPENBLAS_NUM_THREADS': str(threads), 'OMP_NUM_THREADS': str(threads),
        'MKL_NUM_THREADS': str(threads), 'VECLIB_MAXIMUM_THREADS': str(threads),
        'POLARS_MAX_THREADS': str(threads)}, timeout=30) for threads in (1, 4)]
    assert outputs[0] == outputs[1]


def test_deterministic_pca_restores_initialized_blas_limits():
    from threadpoolctl import threadpool_info, threadpool_limits

    with threadpool_limits(limits=2, user_api='blas'):
        before = threadpool_info()
        if sys.platform == 'darwin':
            import ctypes

            accelerate = ctypes.CDLL('/System/Library/Frameworks/Accelerate.framework/Accelerate')
            previous = accelerate.BLASGetThreading()
        with TE._single_threaded_blas():
            assert all(p['num_threads'] == 1 for p in threadpool_info() if p['user_api'] == 'blas')
            if sys.platform == 'darwin':
                assert accelerate.BLASGetThreading() == 1
        assert threadpool_info() == before
        if sys.platform == 'darwin':
            assert accelerate.BLASGetThreading() == previous



def test_release_and_reproducibility_use_the_same_blas_environment():
    from pathlib import Path

    import yaml

    root = Path(__file__).resolve().parents[2] / '.github' / 'workflows'
    release = yaml.safe_load((root / 'build.yml').read_text())['jobs']['build']['env']
    probe = yaml.safe_load((root / 'tcremp-reproducibility.yml').read_text())['jobs']['probe']['env']
    for name in ('OPENBLAS_CORETYPE', 'OPENBLAS_L2_SIZE'):
        assert release[name] == probe[name]
