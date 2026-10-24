"""Exact bounded inference for frozen discrete record-score networks."""
import numpy as np
import polars as pl

TARGET = "validated"


def predict_network(model, data, *, epitopes=None, max_cells=4_000_000, bounds=False):
    """Exact batched posteriors, summing joint tables over missing predictors.

    The small evidence vocabulary permits bounded joint enumeration. max_cells
    rejects larger networks before allocation; use variable elimination for such
    models rather than silently allocating an unbounded Cartesian product.
    """
    if not isinstance(data, pl.DataFrame):
        data = pl.from_pandas(data)
    if model.get('format') != 'vdjdb-learned-bn-1':
        raise ValueError('Unsupported model format')
    if 'epitope_models' in model:
        context = np.asarray(epitopes)
        if context.shape != (len(data),):
            raise ValueError('Aligned epitope context required for epitope-specific models')
        base = {k: v for k, v in model.items() if k != 'epitope_models'}
        result = np.asarray(predict_network(base, data, max_cells=max_cells, bounds=bounds))
        for epitope, local in sorted(model['epitope_models'].items()):
            mask = context == epitope
            if mask.any():
                result[mask] = predict_network(local, data.filter(mask),
                                               max_cells=max_cells, bounds=bounds)
        return result.tolist()
    states = model['states']
    names = list(states)
    if set(data.columns) - set(names):
        raise ValueError('Unknown prediction fields')
    shape = [2, *(len(states[n]) for n in names)]
    if np.prod(shape, dtype=np.int64) > max_cells:
        raise ValueError('Exact joint exceeds max_cells; use variable elimination')
    joint = np.broadcast_to(np.asarray(model['prior']).reshape([2] + [1] * len(names)), shape).copy()
    axes = {n: i + 1 for i, n in enumerate(names)}
    for factor in model['factors']:
        family = [TARGET, *factor['parents'], factor['name']]
        family_axes = [0, *(axes[n] for n in family[1:])]
        table = np.asarray(factor['values']).reshape([shape[i] for i in family_axes])
        order = np.argsort(family_axes)
        table = table.transpose(order)
        broadcast = [1] * len(shape)
        for axis in family_axes:
            broadcast[axis] = shape[axis]
        joint *= table.reshape(broadcast)
    if not np.isfinite(joint).all() or not np.isclose(joint.sum(), 1.):
        raise ValueError('Model is not a normalized joint distribution')
    encoded = np.full((len(data), len(names)), -1, dtype=np.int64)
    for i, name in enumerate(names):
        if name not in data:
            continue
        column = data[name]
        known = (column.is_not_null() & ~column.is_in(['', '?'])).to_numpy()
        mapped = column.replace_strict({v: j for j, v in enumerate(states[name])},
                                       default=-1, return_dtype=pl.Int64).to_numpy()
        if (mapped[known] < 0).any():
            raise ValueError(f'Unknown prediction state in {name}')
        encoded[known, i] = mapped[known].astype(int)
    result = np.empty((len(data), 3)) if bounds else np.empty(len(data))
    masks, assignments = np.unique(encoded >= 0, axis=0, return_inverse=True)
    for group, mask in enumerate(masks):
        rows = np.flatnonzero(assignments == group)
        marginal = joint.sum(axis=tuple(np.flatnonzero(~mask) + 1))
        indices = tuple(encoded[rows, i] for i in np.flatnonzero(mask))
        positive = marginal[(1, *indices)]
        total = marginal[(0, *indices)] + positive
        if bounds:
            # Bounds range over completions of absent predictors; they are not CIs.
            posterior = joint[1] / joint.sum(axis=0)
            lower, upper = posterior, posterior
            for axis in reversed(np.flatnonzero(~mask)):
                legal = [i for i, state in enumerate(states[names[axis]]) if state != '?']
                lower = np.take(lower, legal, axis=axis).min(axis=axis)
                upper = np.take(upper, legal, axis=axis).max(axis=axis)
            nominal = positive / total
            result[rows] = np.column_stack([nominal, np.minimum(nominal, lower[indices]),
                                           np.maximum(nominal, upper[indices])])
        else:
            result[rows] = positive / total
    return result.tolist()
