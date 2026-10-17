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
