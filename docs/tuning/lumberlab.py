# 2026-09-26  Per-epitope Lumbermark in the tcremp label contract. Shared by the sweeps.
import numpy as np


def lumbermark_labels(X, epitopes, min_cluster_size, m_smooth, gate=None):
    """Per-epitope Lumbermark in the tcremp label contract: -1 noise, unique across epitopes."""
    from lumbermark import Lumbermark

    out = np.full(len(X), -1, dtype=np.int64)
    offset = 0
    for ep in np.unique(epitopes):                 # np.unique sorts: the offset must not vary
        m = epitopes == ep if gate is None else (epitopes == ep) & gate
        idx = np.flatnonzero(m)
        # Too few to split into two floors, or fewer points than the core-distance neighbourhood
        # needs: no claim is made for that epitope rather than one cluster covering all of it.
        if len(idx) < max(2 * min_cluster_size, m_smooth + 2):
            continue
        sub = np.ascontiguousarray(X[idx], dtype=np.float64)
        lab = Lumbermark(n_clusters=len(idx) - 1, M=m_smooth, min_cluster_size=min_cluster_size,
                         min_cluster_factor=0.0).fit(sub).labels_
        hit = lab >= 0
        out[idx[hit]] = lab[hit] + offset
        offset += int(lab.max()) + 1 if hit.any() else 0
    return out
