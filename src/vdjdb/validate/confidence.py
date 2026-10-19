"""Matched binary/discrete/continuous score comparisons with explicit error denominators."""
from __future__ import annotations

import polars as pl


def metrics(predictions: pl.DataFrame) -> pl.DataFrame:
    """Summarize saved held-out forecasts; never fit, choose thresholds or drop unknown calls.

    Input columns: case, epitope, predictor, fold, y, score and predicted. Each
    predictor must contain the identical labelled cases/folds within each epitope.
    ``predicted`` is the original binary decision, not a newly optimized cutoff.
    """
    from sklearn.metrics import average_precision_score, balanced_accuracy_score, matthews_corrcoef

    key = ['epitope', 'case', 'predictor']
    if predictions.select(key).n_unique() != predictions.height:
        raise ValueError('Duplicate case/predictor')
    if predictions.select(pl.any_horizontal(
            pl.col('y', 'score', 'predicted', 'fold').is_null()).any()).item():
        raise ValueError('Metrics require complete labels, scores, folds and decisions')
    if predictions.filter(~pl.col('y').is_in([0, 1]) | ~pl.col('predicted').is_in([0, 1]) |
                          ~pl.col('score').is_finite()).height:
        raise ValueError('Invalid label, decision or score')
    rows = []
    for (epitope,), group in sorted(predictions.group_by('epitope'), key=lambda x: x[0]):
        reference = None
        for (name,), part in sorted(group.group_by('predictor'), key=lambda x: x[0]):
            part = part.sort('case')
            matched = part.select('case', 'fold', 'y')
            if reference is not None and not reference.equals(matched):
                raise ValueError('Comparisons require identical cases, folds and outcomes')
            reference = matched
            y, called = part['y'], part['predicted']
            if y.n_unique() != 2:
                raise ValueError('Both outcome classes are required')
            tp = int(((y == 1) & (called == 1)).sum())
            fp = int(((y == 0) & (called == 1)).sum())
            fn = int(((y == 1) & (called == 0)).sum())
            tn = int(((y == 0) & (called == 0)).sum())
            rows.append({'epitope': epitope, 'predictor': name, 'n': len(y),
                         'tp': tp, 'fp': fp, 'fn': fn, 'tn': tn,
                         'ap': average_precision_score(y, part['score']),
                         'accuracy': (tp + tn) / len(y),
                         'balanced_accuracy': balanced_accuracy_score(y, called),
                         'mcc': matthews_corrcoef(y, called),
                         'fpr': fp / (fp + tn),
                         'false_discovery_fraction': fp / (tp + fp) if tp + fp else None})
    return pl.DataFrame(rows).sort('epitope', 'predictor')


def check_nonregression(table: pl.DataFrame, *, joint: str = 'v3') -> None:
    """Require joint AP/accuracy/balanced accuracy at least as high as every supplied baseline."""
    for (_,), group in table.group_by('epitope'):
        selected = group.filter(pl.col('predictor') == joint)
        if selected.height != 1 or group.height < 2:
            raise ValueError('Each epitope requires one joint score and at least one baseline')
        for name in ['ap', 'accuracy', 'balanced_accuracy']:
            if group.filter(pl.col(name) > selected[name][0]).height:
                raise ValueError(f'Joint score regresses {name} for {group["epitope"][0]}')
