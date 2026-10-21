"""Probe TCREmp numerical stages across fresh processes and Linux CPU allocations.

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
    parser.add_argument('--primitive-only', action='store_true')
    args = parser.parse_args()
    # Set before importing either BLAS provider; the probe varies the caller's allocation.
    for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
                 'VECLIB_MAXIMUM_THREADS', 'POLARS_MAX_THREADS'):
        os.environ[name] = str(args.threads)

    import numpy as np
    import polars as pl
    from mir.embedding import TCREmp
    from threadpoolctl import threadpool_info
    from scipy import linalg
    from sklearn.preprocessing import StandardScaler
    from vdjdb.motifs import tcremp as TE

    def digest(array):
        return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()

    if args.primitive_only:
        print(json.dumps({'pools': threadpool_info(), 'cpu':
              Path('/proc/cpuinfo').read_text().split('\n\n')[0]
              if platform.system() == 'Linux' else platform.processor()}), flush=True)
        A = np.random.RandomState(TE.SEED).normal(size=(1024, 6000)).T
        Q = np.random.RandomState(TE.SEED).normal(size=(1024, 60))
        with TE._single_threaded_blas():
            for repeat in range(2):
                print(f'GEMM {repeat}: {A.shape} @ {Q.shape}', flush=True)
                print(digest(A @ Q), flush=True)
        return

    # Larger than sklearn's randomized-PCA switch, as in the production 25k fit.
    rng = np.random.default_rng(TE.SEED)
    aa = np.array(list('ACDEFGHIKLMNPQRSTVWY'))
    junctions = ['C' + ''.join(row) + 'F' for row in rng.choice(aa, (1024, 12))]
    inputs = pl.DataFrame({'junction_aa': junctions,
                           'v_call': ['TRBV5-1*01'] * len(junctions),
                           'j_call': ['TRBJ1-1*01'] * len(junctions)})
    raw = []
    numerical = []
    active_pools = []
    original = TCREmp.embed

    # Record boundaries inside randomized PCA without replacing its calculations.
    def record_function(owner, name):
        function = getattr(owner, name)

        def wrapped(*positional, **keywords):
            if not active_pools:
                active_pools.extend(threadpool_info())
            array = positional[-1] if name == 'transform' else positional[0]
            before = digest(array)
            print(f"TCREmp probe {name}: input {array.shape} {array.dtype}", flush=True)
            result = function(*positional, **keywords)
            arrays = result if isinstance(result, tuple) else (result,)
            numerical.append({'operation': name, 'input': before,
                              'outputs': [digest(a) for a in arrays]})
            return result

        setattr(owner, name, wrapped)

    for name in ('lu', 'qr', 'svd'):
        record_function(linalg, name)
    record_function(StandardScaler, 'transform')

    def recording_embed(self, frame):
        matrix = original(self, frame)
        raw.append(digest(matrix))
        return matrix

    TCREmp.embed = recording_embed
    print(json.dumps({'before_fit_pools': threadpool_info(),
                      'cpu': Path('/proc/cpuinfo').read_text().split('\n\n')[0]
                      if platform.system() == 'Linux' else platform.processor()}), flush=True)
    X, frozen = TE.fit_embedding(inputs, 'HomoSapiens', 'TRB', fit_sample=1024, chunk=512)
    eps = TE.chain_eps(X, TE.COEF['TRB'], clonotypes=inputs)
    labels = TE.cluster_labels(X, np.array(['PEPTIDE'] * len(X)), eps)
    pools = threadpool_info()
    if platform.system() == 'Linux' and os.environ.get('OPENBLAS_CORETYPE') == 'Haswell':
        blas = [pool for pool in pools if pool['internal_api'] == 'openblas']
        assert blas and all(pool['architecture'] == 'Haswell' for pool in blas), pools
    result = {
        'threads': args.threads, 'machine': platform.machine(),
        'processor': platform.processor(), 'platform': platform.platform(),
        'libc': platform.libc_ver(),
        'numpy_cpu_features': np._core._multiarray_umath.__cpu_features__,
        'numpy': np.__version__, 'pools': pools, 'active_pools': active_pools,
        'solver': frozen.pca._fit_svd_solver, 'shape': list(X.shape),
        'prototype_features': int(frozen.pca.n_features_in_),
        'pca_dtype': str(frozen.pca.components_.dtype),
        'clustered': int((labels >= 0).sum()),
        'clusters': int(len(np.unique(labels[labels >= 0]))),
        'stages': {'raw': raw, 'numerical': numerical,
                   'random_projection': digest(np.random.RandomState(TE.SEED).normal(
                       size=(1024, TE.N_COMPONENTS + 10)).astype(frozen.pca.components_.dtype)), 'scale_mean': digest(frozen.scaler.mean_),
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
