"""Optional continuous motif indices; hard producer assignments remain authoritative."""
from __future__ import annotations

import numpy as np


def separated(binary: np.ndarray, support: np.ndarray) -> np.ndarray:
    """Map support to [0,1/3] for noise and [2/3,1] for members; preserve missingness.

    This index preserves hard ordering and is not a calibrated binding probability.
    """
    b, u = np.asarray(binary, dtype=float), np.asarray(support, dtype=float)
    if b.shape != u.shape:
        raise ValueError('Membership and support must have identical shapes')
    if np.any(~np.isnan(b) & ~np.isin(b, [0, 1])):
        raise ValueError('Membership must be binary or missing')
    if np.any(~np.isnan(u) & (~np.isfinite(u) | (u < 0) | (u > 1))):
        raise ValueError('Support must be in [0,1] or missing')
    return (2 * b + u) / 3
