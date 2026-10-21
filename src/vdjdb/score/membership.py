"""Optional continuous motif indices; hard producer assignments remain authoritative."""
from __future__ import annotations

import numpy as np


def separated(binary: np.ndarray, support: np.ndarray, *, width: float = 1/3) -> np.ndarray:
    """Map noise to [0,width] and members to [1-width,1]; preserve missingness.

    This index preserves hard ordering and is not a calibrated binding probability.
    A width in [0,0.5) controls departure from hard evidence in soft inference.
    Zero recovers hard evidence where both inputs are present; the default is 1/3.
    """
    if not np.isfinite(width) or not 0 <= width < .5:
        raise ValueError('width must be finite and in [0,0.5)')
    b, u = np.asarray(binary, dtype=float), np.asarray(support, dtype=float)
    if b.shape != u.shape:
        raise ValueError('Membership and support must have identical shapes')
    if np.any(~np.isnan(b) & ~np.isin(b, [0, 1])):
        raise ValueError('Membership must be binary or missing')
    if np.any(~np.isnan(u) & (~np.isfinite(u) | (u < 0) | (u > 1))):
        raise ValueError('Support must be in [0,1] or missing')
    if width == 1/3:
        return (2 * b + u) / 3
    return (1-width)*b + width*u
