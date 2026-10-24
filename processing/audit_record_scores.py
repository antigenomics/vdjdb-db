"""Audit aligned score forecasts with paired per-epitope and macro differences.

Input TSV columns: endpoint, epitope, case, label, predictor, score.
Every predictor must cover identical cases/labels. Keep private case tables and
outputs outside Git. No fitting occurs here.
"""
import argparse
from pathlib import Path

import numpy as np
import polars as pl

from vdjdb.config import SEED


def weighted_metrics(y, score, weights):
    order = np.argsort(-score, kind='stable')
    ends = np.r_[np.flatnonzero(np.diff(score[order]) != 0), len(y)-1]
    tp = np.cumsum(weights[:, order] * y[order], axis=1)[:, ends]
    fp = np.cumsum(weights[:, order] * (1-y[order]), axis=1)[:, ends]
    recall, fpr = tp / tp[:, -1:], fp / fp[:, -1:]
    precision = np.divide(tp, tp+fp, out=np.ones_like(tp), where=(tp+fp) > 0)
    dr = np.diff(recall, prepend=0, axis=1)
    return np.column_stack([
        np.sum(dr * precision, axis=1),
        np.sum(dr * (precision + np.column_stack([np.ones(len(weights)), precision[:, :-1]])) / 2, axis=1),
        np.sum(np.diff(fpr, prepend=0, axis=1) *
               (recall + np.column_stack([np.zeros(len(weights)), recall[:, :-1]])) / 2, axis=1)])


def audit(cases, candidate, comparators, draws=10000):
    if draws < 2:
        raise ValueError('At least two bootstrap draws required')
    fields = ['endpoint', 'epitope', 'case', 'predictor']
    if cases.select(pl.struct(fields).n_unique()).item() != cases.height:
        raise ValueError('Duplicate forecast cases')
    if cases.select(pl.any_horizontal(pl.all().is_null()).any()).item():
        raise ValueError('Missing forecast fields')
    if cases.filter(~pl.col('label').is_in([0, 1]) | ~pl.col('score').is_finite()).height:
        raise ValueError('Binary labels and finite scores required; exclude unknown outcomes explicitly')
    names = [candidate, *comparators]
    if set(cases['predictor'].unique()) != set(names) or len(names) != len(set(names)):
        raise ValueError('Provide exactly the candidate and distinct comparators')
    points, differences = [], []

    def compare(endpoint, epitope, distributions, values):
        for comparator in comparators:
            intervals = np.quantile(distributions[candidate]-distributions[comparator], [.025, .975], axis=0)
            for j, metric in enumerate(['ap', 'pr_auc', 'roc_auc']):
                differences.append(dict(endpoint=endpoint, epitope=epitope, candidate=candidate,
                    comparator=comparator, metric=metric,
                    difference=float(values[candidate][j]-values[comparator][j]),
                    ci95_low=float(intervals[0, j]), ci95_high=float(intervals[1, j]),
                    regression=bool(intervals[1, j] < 0), draws=draws))

    for endpoint_index, endpoint in enumerate(sorted(cases['endpoint'].unique())):
        part = cases.filter(pl.col('endpoint') == endpoint)
        macro, point_macro = {n: [] for n in names}, {n: [] for n in names}
        for epitope_index, epitope in enumerate(sorted(part['epitope'].unique())):
            group = part.filter(pl.col('epitope') == epitope)
            reference = group.filter(pl.col('predictor') == comparators[0]).sort('case')
            y = reference['label'].to_numpy()
            if set(y) != {0, 1}:
                raise ValueError('Every evaluated epitope requires both classes')
            rng = np.random.default_rng(SEED + 100 * endpoint_index + epitope_index)
            weights = np.zeros((draws, len(y)))
            for cls in (0, 1):
                mask = y == cls
                weights[:, mask] = rng.multinomial(mask.sum(), np.full(mask.sum(), 1/mask.sum()), size=draws)
            distributions, values = {}, {}
            for name in names:
                rows = group.filter(pl.col('predictor') == name).sort('case')
                if not rows.select('case', 'label').equals(reference.select('case', 'label')):
                    raise ValueError('Forecasts must have identical cases and labels')
                score = rows['score'].to_numpy()
                values[name] = weighted_metrics(y, score, np.ones((1, len(y))))[0]
                distributions[name] = weighted_metrics(y, score, weights)
                macro[name].append(distributions[name])
                point_macro[name].append(values[name])
                for j, metric in enumerate(['ap', 'pr_auc', 'roc_auc']):
                    points.append(dict(endpoint=endpoint, epitope=epitope, model=name, metric=metric,
                                       value=float(values[name][j]), n=len(y), positive=int(y.sum())))
            compare(endpoint, epitope, distributions, values)
        compare(endpoint, 'macro', {n: np.mean(v, axis=0) for n, v in macro.items()},
                {n: np.mean(v, axis=0) for n, v in point_macro.items()})
    return pl.DataFrame(points), pl.DataFrame(differences)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases', type=Path, required=True)
    parser.add_argument('--candidate', required=True)
    parser.add_argument('--comparators', nargs='+', required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        parser.error('Use a new output directory')
    points, differences = audit(pl.read_csv(args.cases, separator='\t',
                                          schema_overrides={'case': pl.String}), args.candidate, args.comparators)
    args.out.mkdir(parents=True)
    points.write_csv(args.out/'metric_points.tsv', separator='\t')
    differences.write_csv(args.out/'paired_differences.tsv', separator='\t')
    if differences['regression'].any():
        raise SystemExit('Regression gate failed; inspect paired_differences.tsv')
    print(f'{differences.height} paired regression gates passed')


if __name__ == '__main__':
    main()
