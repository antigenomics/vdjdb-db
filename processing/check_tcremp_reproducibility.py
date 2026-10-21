"""Check TCREmp stages across fresh processes and Linux CPU allocations.

Uses generated productive junctions, no assay data or stored fitted state. Each
invocation re-embeds and refits through the production entry point.
"""
from __future__ import annotations

import argparse
import faulthandler
import hashlib
import json
import os
import platform
from pathlib import Path


def main() -> None:
    faulthandler.enable()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--threads', type=int, required=True, choices=(1, 2, 4))
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    # Vary the caller's allocation before importing either BLAS provider.
    for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
                 'VECLIB_MAXIMUM_THREADS', 'POLARS_MAX_THREADS'):
        os.environ[name] = str(args.threads)

    import numpy as np
    import polars as pl
    from mir.embedding import TCREmp
    from threadpoolctl import threadpool_info
    from vdjdb.motifs import tcremp as TE

    def digest(array):
        return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()

    # Match the production 6,000-feature randomized-PCA path with a bounded fit.
    rng = np.random.default_rng(TE.SEED)
    aa = np.array(list('ACDEFGHIKLMNPQRSTVWY'))
    junctions = ['C' + ''.join(row) + 'F' for row in rng.choice(aa, (1024, 12))]
    inputs = pl.DataFrame({'junction_aa': junctions,
                           'v_call': ['TRBV5-1*01'] * len(junctions),
                           'j_call': ['TRBJ1-1*01'] * len(junctions)})
    raw = []
    original = TCREmp.embed

    def recording_embed(self, frame):
        matrix = original(self, frame)
        raw.append(digest(matrix))
        return matrix

    TCREmp.embed = recording_embed
    X, frozen = TE.fit_embedding(inputs, 'HomoSapiens', 'TRB', fit_sample=1024, chunk=512)
    assert frozen.pca.components_.dtype == np.float64
    eps = TE.chain_eps(X, TE.COEF['TRB'], clonotypes=inputs)
    labels = TE.cluster_labels(X, np.array(['PEPTIDE'] * len(X)), eps)
    pools = threadpool_info()
    expected_core = os.environ.get('OPENBLAS_CORETYPE')
    if platform.system() == 'Linux' and expected_core:
        blas = [pool for pool in pools if pool['internal_api'] == 'openblas']
        assert blas and all(pool['architecture'].lower() == expected_core.lower() for pool in blas), pools
    cpu = platform.processor()
    if platform.system() == 'Linux':
        cpu = next((line.split(':', 1)[1].strip()
                    for line in Path('/proc/cpuinfo').read_text().splitlines()
                    if line.startswith('model name')), cpu)
    result = {
        'threads': args.threads, 'machine': platform.machine(), 'cpu': cpu,
        'platform': platform.platform(), 'numpy': np.__version__, 'pools': pools,
        'solver': frozen.pca._fit_svd_solver, 'shape': list(X.shape),
        'prototype_features': int(frozen.pca.n_features_in_),
        'pca_dtype': str(frozen.pca.components_.dtype),
        'openblas_l2_size_kib': os.environ.get('OPENBLAS_L2_SIZE', ''),
        'clustered': int((labels >= 0).sum()),
        'clusters': int(len(np.unique(labels[labels >= 0]))),
        'stages': {'raw': raw, 'scale_mean': digest(frozen.scaler.mean_),
                   'scale_variance': digest(frozen.scaler.var_),
                   'pca_mean': digest(frozen.pca.mean_),
                   'pca_components': digest(frozen.pca.components_),
                   'coordinates': digest(X), 'labels': digest(labels), 'eps': eps},
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
