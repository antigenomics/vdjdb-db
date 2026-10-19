"""Metrics and training-only fitting for standalone motif research assessments.

No assay is read here. Callers own grouped splits and must pass training rows only
when fitting paired coefficients or thresholds. Partial AP retains tied score groups.
"""
from __future__ import annotations

import numpy as np


def _inputs(y, score):
    y, score = np.asarray(y), np.asarray(score, dtype=float)
    if y.ndim != 1 or y.shape != score.shape or not len(y):
        raise ValueError('Expected nonempty matching label and score vectors')
    if not np.isin(y, [0, 1]).all() or not np.isfinite(score).all():
        raise ValueError('Labels must be binary and scores finite')
    return y.astype(int), score


def ranking_metrics(y, score) -> dict:
    """AP, ROC-AUC, PR area and normalized partial AP at recall 0.1, 0.2 and 0.5.

    partial_ap_rR = integral_0^R precision(r) dr / R, using AP's constant
    precision over each threshold group's recall increment. Ties enter together.
    PR trapezoids are reported separately. ROC partial AUC uses sklearn's
    standardized max_fpr convention and is named by FPR, never by recall.
    """
    from sklearn.metrics import auc, average_precision_score, precision_recall_curve, roc_auc_score

    y, score = _inputs(y, score)
    positives = int(y.sum())
    precision, recall, _ = precision_recall_curve(y, score)
    out = {'n': len(y), 'positive': positives, 'negative': len(y) - positives,
           'ap': float(average_precision_score(y, score)) if positives else None,
           'pr_auc_trapezoid': float(auc(recall, precision)) if positives else None,
           'roc_auc': float(roc_auc_score(y, score)) if 0 < positives < len(y) else None,
           'roc_auc_fpr0.1': float(roc_auc_score(y, score, max_fpr=.1))
           if 0 < positives < len(y) else None}
    r, p = recall[::-1], precision[::-1]
    for limit in [.1, .2, .5]:
        width = np.maximum(0, np.minimum(r[1:], limit) - np.minimum(r[:-1], limit))
        out[f'partial_ap_r{limit}'] = float(width @ p[1:] / limit) if positives else None
        # Precision at the first attainable nonempty recall >= limit, keeping ties.
        ok = np.flatnonzero(r >= limit)
        out[f'precision_r{limit}'] = float(p[ok[0]]) if positives and len(ok) else None
    return out


def decision_metrics(y, score, threshold=.5) -> dict:
    """Binary decisions at a scalar or per-observation training-selected threshold."""
    y, score = _inputs(y, score)
    threshold = np.asarray(threshold, dtype=float)
    if np.isnan(threshold).any() or threshold.ndim > 1 or (threshold.ndim and threshold.shape != y.shape):
        raise ValueError('Threshold must be a scalar or matching vector without NaN')
    call = score >= threshold
    tp, fp = int(np.sum(call & (y == 1))), int(np.sum(call & (y == 0)))
    fn, tn = int(np.sum(~call & (y == 1))), int(np.sum(~call & (y == 0)))
    return {'tp': tp, 'fp': fp, 'fn': fn, 'tn': tn,
            'precision': tp/(tp+fp) if tp+fp else 0., 'recall': tp/(tp+fn) if tp+fn else 0.,
            'f1': 2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else 0.,
            'accuracy': (tp+tn)/len(y),
            'balanced_accuracy': .5*(tp/(tp+fn)+tn/(tn+fp)) if tp+fn and tn+fp else None}


def select_threshold(y, score) -> float:
    """Maximize training F1 over attainable decisions; ties nearest .5, then larger.

    Include .5 explicitly to retain the original decision when it is optimal.
    The threshold above the maximum score permits an all-negative call.
    """
    y, score = _inputs(y, score)
    order = np.argsort(-score, kind='stable')
    s, labels = score[order], y[order]
    end = np.r_[np.flatnonzero(s[:-1] != s[1:]), len(s)-1]
    tp = np.cumsum(labels)[end]
    f1 = 2*tp / np.maximum(y.sum() + end + 1, 1)
    thresholds = np.r_[s[end], np.nextafter(s[0], np.inf), .5]
    values = np.r_[f1, 0., decision_metrics(y, score, .5)['f1']]
    best = np.flatnonzero(values == values.max())
    chosen = np.lexsort((-thresholds[best], np.abs(thresholds[best]-.5)))[0]
    return float(thresholds[best[chosen]])


def fit_pair(x, y, *, interaction: bool = False) -> dict:
    """Fit a monotone logistic pair formula with fixed L2 penalty 1e-3.

    sigmoid(intercept + wa*a + wb*b [+ wab*a*b]); nonnegative weights.
    """
    from scipy.optimize import minimize
    from scipy.special import expit

    x, y = np.asarray(x, dtype=float), np.asarray(y)
    if (x.shape != (len(y), 2) or not len(y) or not np.isfinite(x).all()
            or not np.isin(y, [0, 1]).all() or np.any((x < 0) | (x > 1))):
        raise ValueError('Paired training needs finite [0,1] chain scores and binary labels')
    z = np.column_stack([np.ones(len(x)), x, x[:, 0]*x[:, 1]]) if interaction else np.c_[np.ones(len(x)), x]
    def objective(w):
        eta = z @ w
        penalty = .001 * np.sum(w[1:]**2) / 2
        loss = np.mean(np.logaddexp(0, eta)-y*eta) + penalty
        grad = z.T @ (expit(eta)-y) / len(y)
        grad[1:] += .001*w[1:]
        return loss, grad
    result = minimize(objective, np.zeros(z.shape[1]), jac=True, method='L-BFGS-B',
                      bounds=[(None, None)]+[(0, None)]*(z.shape[1]-1), options={'ftol': 1e-12})
    if not result.success:
        raise RuntimeError(f'Paired fit failed: {result.message}')
    return {'intercept': float(result.x[0]), 'weights': result.x[1:].tolist(),
            'interaction': interaction, 'l2': .001}


def predict_pair(model: dict, x) -> np.ndarray:
    """Evaluate the saved two-chain formula without refitting."""
    from scipy.special import expit

    x = np.asarray(x, dtype=float)
    w = np.asarray(model['weights'], dtype=float)
    if (x.ndim != 2 or x.shape[1] != 2 or not np.isfinite(x).all()
            or np.any((x < 0) | (x > 1)) or not np.isfinite(w).all() or np.any(w < 0)
            or w.shape != (3 if model['interaction'] else 2,)
            or not np.isfinite(model['intercept'])):
        raise ValueError('Invalid paired model or scores')
    z = np.c_[x, x[:, 0]*x[:, 1]] if model['interaction'] else x
    return expit(model['intercept'] + z @ w)


def fit_thresholds(y, score, peptides) -> dict:
    """Fit a global F1 cutoff and two-class peptide cutoffs on training data only.

    Peptides with a single observed class use the global cutoff, as do unseen
    peptides. No minimum sample count or class label is inferred for missing data.
    """
    y, score = _inputs(y, score)
    peptides = np.asarray(peptides, dtype=str)
    if peptides.shape != y.shape:
        raise ValueError('Peptides must match the training rows')
    cuts = {ep: select_threshold(y[peptides == ep], score[peptides == ep])
            for ep in np.unique(peptides) if ep and np.unique(y[peptides == ep]).size == 2}
    return {'global_threshold': select_threshold(y, score), 'epitope_thresholds': cuts,
            'fallback': 'global threshold for unseen or single-class training peptides'}


def apply_thresholds(policy: dict, peptides) -> np.ndarray:
    """Resolve peptide-specific cutoffs with the policy's training-derived fallback."""
    peptides = np.asarray(peptides, dtype=str)
    if peptides.ndim != 1:
        raise ValueError('Expected a peptide vector')
    cutoff = float(policy['global_threshold'])
    specific = policy['epitope_thresholds']
    values = np.array([specific.get(ep, cutoff) for ep in peptides], dtype=float)
    if not np.isfinite(cutoff) or not np.isfinite(values).all():
        raise ValueError('Threshold policy must contain finite cutoffs')
    return values
