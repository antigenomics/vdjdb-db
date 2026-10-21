"""Small discrete Bayesian evidence classifier with exact missing-value inference.

Install the ``confidence`` extra for fitting. Missing fields are marginalized at prediction;
no evidence returns the fitted network prior. Fitted artifacts contain no receptor identities.
"""
from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from itertools import product

from .evidence import COUNT_FEATURES, FEATURES, categorical

STATES = {name: (['?', '0', '1', '2', '3+'] if name in COUNT_FEATURES else
                 ['?', 'low', 'mid', 'high'] if name == 'iptm' else
                 ['?', '0', '1', '2', '3'] if name == 'legacy_score' else ['?', '0', '1'])
          for name in FEATURES}

# Bound simultaneous nine-feature assignment dictionaries in mixtures and envelopes.
_ASSIGNMENT_BATCH = 4096


def fit(rows: Sequence[Mapping[str, str]], labels: Sequence[int], *,
        equivalent_sample_size: float | Mapping[str, float] = 10., structure: str = 'tan') -> dict:
    """Learn a discrete Bayesian graph and Dirichlet-smoothed conditional tables.

    ``structure='tan'`` learns one feature-dependence tree; ``'naive'`` omits those edges.
    ``'tan_leaf'`` learns the other features first, attaching ipTM as a leaf by conditional
    mutual information. ``'tan_independent_iptm'`` gives ipTM only the validation parent.
    Both keep the nonstructural graph independent of ipTM grading.
    A smoothing map specifies one prior size per node without changing the graph.
    Training uses explicit unreported states, avoiding complete-case exclusion. These states
    describe reporting, not nonbinding. At inference, unknown features are unobserved rather
    than evidence for a class. The resulting probabilities target the training assay endpoint.
    """
    import pandas as pd
    from pgmpy.estimators import BayesianEstimator, TreeSearch
    from pgmpy.models import DiscreteBayesianNetwork
    from sklearn.metrics import mutual_info_score

    if structure not in ('tan', 'tan_leaf', 'tan_independent_iptm', 'naive'):
        raise ValueError('Unsupported Bayesian structure')
    if len(rows) != len(labels) or set(labels) != {0, 1}:
        raise ValueError('Training requires aligned observations and both classes')
    if isinstance(equivalent_sample_size, Mapping):
        if set(equivalent_sample_size) != {*FEATURES, 'validated'}:
            raise ValueError('Smoothing map requires every feature and validated')
        equivalent_sample_size = dict(equivalent_sample_size)
        sizes = equivalent_sample_size.values()
    else:
        sizes = [equivalent_sample_size]
    if any(not math.isfinite(size) or size <= 0 for size in sizes):
        raise ValueError('Positive finite equivalent sample size required')
    if any(set(row) - set(FEATURES) for row in rows):
        raise ValueError('Unknown training predictor fields')
    data = pd.DataFrame([{n: r.get(n, '?') for n in FEATURES} for r in rows])
    for name in FEATURES:
        if not set(data[name]) <= set(STATES[name]):
            raise ValueError(f'Invalid training state for {name}')
    data['validated'] = [str(int(y)) for y in labels]
    data = data.astype('category')
    # Explicit lexical tie break makes constant/tied features deterministic.
    tree_features = [n for n in FEATURES
                     if structure not in ('tan_leaf', 'tan_independent_iptm') or n != 'iptm']
    root = min(tree_features, key=lambda n: (-mutual_info_score(data[n], data['validated']), n))
    if structure in ('tan', 'tan_leaf', 'tan_independent_iptm'):
        graph = TreeSearch(data if structure == 'tan' else data.drop(columns=['iptm']),
                           root_node=root, n_jobs=1).estimate(
            estimator_type='tan', class_node='validated', show_progress=False)
        edges = sorted(graph.edges())
        if structure == 'tan_leaf':
            def conditional_information(name):
                return sum((data['validated'] == cls).mean() * mutual_info_score(
                    data.loc[data['validated'] == cls, name],
                    data.loc[data['validated'] == cls, 'iptm']) for cls in ('0', '1'))

            parent = min(tree_features, key=lambda n: (-conditional_information(n), n))
            edges.extend([(parent, 'iptm'), ('validated', 'iptm')])
        elif structure == 'tan_independent_iptm':
            edges.append(('validated', 'iptm'))
    else:
        edges = [('validated', n) for n in FEATURES]
    model = DiscreteBayesianNetwork(edges)
    states = {**STATES, 'validated': ['0', '1']}
    model.fit(data, estimator=BayesianEstimator, prior_type='BDeu',
              equivalent_sample_size=equivalent_sample_size, state_names=states, n_jobs=1)
    if not model.check_model():
        raise ValueError('Invalid fitted Bayesian network')
    cpds = []
    for node in sorted(model.nodes()):
        cpd = model.get_cpds(node)
        cpds.append({'variable': node, 'variables': list(cpd.variables),
                     'cardinality': cpd.cardinality.tolist(), 'values': cpd.get_values().tolist(),
                     'states': cpd.state_names})
    return {'format': 'vdjdb-evidence-tan-1', 'features': FEATURES, 'edges': sorted(map(list, model.edges())),
            'cpds': cpds, 'observed_states': {n: data[n].value_counts().to_dict() for n in FEATURES},
            'structure': structure, 'root': root, 'n_training': len(labels),
            'equivalent_sample_size': equivalent_sample_size,
            'missing_fit': 'unreported category', 'missing_predict': 'exact marginalization'}


def predict(artifact: dict, rows: Sequence[Mapping[str, str | float | None]], *,
            soft_membership: bool = False, _state: dict | None = None) -> list[float]:
    """Score categorical evidence or uncertain binary motif membership.

    A numeric motif value q mixes hard nonmember/member posteriors with weights
    (1-q, q). Multiple uncertain memberships use an independent mixing distribution.
    These are caller-supplied membership probabilities, not binding probabilities.
    """
    if _state is None:
        _state = {}
    if any(isinstance(value, (int, float)) for row in rows for name, value in row.items()
           if name in FEATURES and name.endswith(('tcrnet', 'tcremp'))):
        patterns = [tuple(sorted(row.items())) for row in rows]
        unique = dict(zip(patterns, rows, strict=True))
        expanded, groups, values = [], [], {}

        def flush():
            probabilities = iter(predict(artifact, expanded, _state=_state))
            for pattern, weights in groups:
                values[pattern] = math.fsum(weight * next(probabilities) for weight in weights)
            expanded.clear()
            groups.clear()

        for pattern, row in unique.items():
            mixtures = [(dict(row), 1.)]
            for name, value in row.items():
                if name not in FEATURES or not name.endswith(('tcrnet', 'tcremp')):
                    continue
                if not isinstance(value, (int, float)):
                    continue
                if not math.isfinite(value) or not 0 <= value <= 1:
                    raise ValueError(f'Invalid membership probability for {name}')
                if value not in (0, 1) and not soft_membership:
                    raise ValueError('Fractional motif evidence requires soft_membership=True')
                mixtures = [({**r, name: str(state)}, weight * probability)
                            for r, weight in mixtures
                            for state, probability in ((0, 1-value), (1, value)) if probability > 0]
            if len(expanded) + len(mixtures) > _ASSIGNMENT_BATCH:
                flush()
            expanded.extend(r for r, _ in mixtures)
            groups.append((pattern, [weight for _, weight in mixtures]))
        if expanded:
            flush()
        return [values[pattern] for pattern in patterns]

    from pgmpy.factors.discrete import TabularCPD
    from pgmpy.inference import VariableElimination
    from pgmpy.models import DiscreteBayesianNetwork

    if artifact.get('format') != 'vdjdb-evidence-tan-1' or artifact.get('features') != FEATURES:
        raise ValueError('Unsupported classifier artifact')
    if 'infer' not in _state:
        model = DiscreteBayesianNetwork(artifact['edges'])
        model.add_nodes_from(cpd['variable'] for cpd in artifact['cpds'])
        for cpd in artifact['cpds']:
            variables, card = cpd['variables'], cpd['cardinality']
            model.add_cpds(TabularCPD(cpd['variable'], card[0], cpd['values'],
                                     evidence=variables[1:] or None,
                                     evidence_card=card[1:] or None, state_names=cpd['states']))
        model.check_model()
        _state['infer'] = VariableElimination(model)
    patterns = []
    for row in rows:
        unknown = set(row) - set(FEATURES)
        if unknown:
            raise ValueError(f'Unknown predictor fields: {sorted(unknown)}')
        observed = []
        for name in FEATURES:
            state = row.get(name)
            if state in (None, '', '?'):
                continue
            if state not in STATES[name]:
                raise ValueError(f'Invalid state {state!r} for {name}')
            support = artifact.get('observed_states', {}).get(name)
            if support is not None and support.get(state, 0) == 0:
                continue
            observed.append((name, state))
        patterns.append(tuple(observed))
    infer = _state['infer']
    # Reuse categorical results across bounded batches within this inference call.
    values = _state.setdefault('values', {})
    for pattern in sorted(set(patterns) - values.keys()):
        result = infer.query(['validated'], evidence=dict(pattern), show_progress=False)
        values[pattern] = float(result.get_value(validated='1'))
    return [values[p] for p in patterns]


def predict_marginals(artifact: dict, rows: Sequence[Mapping[str, str | float | None]], *,
                      soft_membership: bool = False
                      ) -> list[dict[str, float]]:
    """Full-interaction validation from paired, TRA-only and TRB-only evidence.

    Single-chain views omit paired legacy score, donor total and structure confidence.
    They condition on that chain's evidence, not on a single-chain binding outcome.
    """
    views = [list(rows), [{k: v for k, v in row.items() if k.startswith('a_')} for row in rows],
             [{k: v for k, v in row.items() if k.startswith('b_')} for row in rows]]
    values = predict(artifact, [row for view in views for row in view], soft_membership=soft_membership)
    return [{name: values[j * len(rows) + i] for j, name in enumerate(('paired', 'TRA', 'TRB'))}
            for i in range(len(rows))]


def _numeric_rows(artifact: dict, frame, *, soft_membership: bool = False) -> list[dict]:
    """Encode numeric evidence consistently for forecasts and missing-evidence ranges."""
    motifs = [name for name in FEATURES if name.endswith(('tcrnet', 'tcremp'))
              and name in frame.columns]
    names = [name for name in FEATURES if name in frame.columns and name not in motifs]
    rows = categorical(frame.select(names), iptm_bounds=tuple(
        artifact.get('iptm_bounds', (.5, .8)))).to_dicts() if names else [{}] * frame.height
    if motifs:
        import polars as pl

        if not soft_membership:
            categorical(frame.select(motifs))  # Reject fractional membership by default.
        membership = frame.select(pl.col(motifs).cast(pl.Float64).fill_nan(None)).to_dicts()
        rows = [{**row, **soft} for row, soft in zip(rows, membership, strict=True)]
    return rows


def predict_evidence(artifact: dict, frame, *, soft_membership: bool = False) -> list[dict[str, float]]:
    """Score numeric evidence using saved grades and one network for all three views.

    Missing columns are unobserved; extra identity columns are not predictors.
    """
    return predict_marginals(artifact, _numeric_rows(artifact, frame, soft_membership=soft_membership),
                             soft_membership=soft_membership)


def _paired_ranges(artifact: dict, rows: list[dict]) -> list[tuple[float, float]]:
    """Sensitivity to missing chain/ipTM evidence, not statistical confidence limits.

    Complete missing chain features and ipTM using states observed in training. Metadata
    remains marginalized. Include the nominal forecast, since explicit unreported states
    in the fitted network can put it outside the observed-state completion envelope.
    """
    patterns = [tuple(row.get(name) for name in FEATURES) for row in rows]
    unique = dict(zip(patterns, rows, strict=True))
    expanded, groups, bounds = [], [], {}
    inference_state = {}

    def flush():
        values = predict(artifact, expanded, soft_membership=True, _state=inference_state)
        offset = 0
        for pattern, size in groups:
            probabilities = values[offset:offset + size]
            bounds[pattern] = (min(probabilities), max(probabilities))
            offset += size
        expanded.clear()
        groups.clear()

    for pattern, row in unique.items():
        missing = {}
        for name in FEATURES:
            if not (name.startswith(('a_', 'b_')) or name == 'iptm'):
                continue
            support = artifact.get('observed_states', {}).get(name)
            value = row.get(name)
            if isinstance(value, (int, float)) and 0 < value < 1:
                continue
            state = str(int(value)) if isinstance(value, (int, float)) else value
            if state not in (None, '', '?') and (support is None or support.get(state, 0) > 0):
                continue
            choices = [state for state in STATES[name] if state != '?' and
                       (support is None or support.get(state, 0) > 0)]
            if choices:
                missing[name] = choices
        completions = [dict(row)]
        if missing:
            completions.extend({**row, **dict(zip(missing, values, strict=True))}
                               for values in product(*missing.values()))
        if len(expanded) + len(completions) > _ASSIGNMENT_BATCH:
            flush()
        groups.append((pattern, len(completions)))
        expanded.extend(completions)
    if expanded:
        flush()
    return [bounds[pattern] for pattern in patterns]


def score_records(artifact: dict, frame, *, soft_membership: bool = False
                  ) -> list[dict[str, float | bool | None]]:
    """TRA/TRB/paired forecasts and a range for missing chain or structural evidence.

    Record score retains the paired forecast even without direct pair evidence.
    Without direct pair evidence, include the chain-only forecasts in the range too.
    A call is positive only when the range is entirely above the model's decision
    threshold (default 0.5), negative when below, and otherwise unknown (None).
    The lower endpoint includes equality. This range is not a credible interval.
    Pair support is descriptive and never forces a score to zero.
    """
    import polars as pl

    threshold = artifact.get('decision_threshold', .5)
    if not isinstance(threshold, (int, float)) or not math.isfinite(threshold) or not 0 <= threshold <= 1:
        raise ValueError('Decision threshold must be finite and in [0,1]')
    rows = _numeric_rows(artifact, frame, soft_membership=soft_membership)
    probabilities = predict_marginals(artifact, rows, soft_membership=soft_membership)
    ranges = _paired_ranges(artifact, rows)
    def value(name):
        return pl.col(name).cast(pl.Float64).fill_nan(None) if name in frame.columns else pl.lit(None)

    chains = [pl.any_horizontal((value(name) > 0 if name.endswith('studies') else
                               value(name + '_hard' if name + '_hard' in frame.columns else name) == 1)
                               for name in FEATURES if name.startswith(prefix))
              for prefix in ('a_', 'b_')]
    direct = (value('legacy_score').is_not_null() | (value('donors') > 0) |
              value('iptm').is_not_null())
    pair = (pl.all_horizontal(pl.col('a_junction', 'b_junction').fill_null('') != '')
            if {'a_junction', 'b_junction'} <= set(frame.columns) else pl.lit(True))
    direct = direct & pair
    support = frame.with_columns((pair & (direct | (chains[0] & chains[1]))).fill_null(False)
                                 .alias('_pair_supported'), direct.fill_null(False).alias('_direct_pair'))
    result = []
    for p, supported, paired, (lower, upper) in zip(
            probabilities, support['_pair_supported'], support['_direct_pair'], ranges, strict=True):
        if not paired:
            lower = min(lower, p['TRA'], p['TRB'])
            upper = max(upper, p['TRA'], p['TRB'])
        result.append({**p, 'pair_supported': bool(supported), 'direct_pair_evidence': bool(paired),
                       'record_score': p['paired'], 'paired_lower': lower, 'paired_upper': upper,
                       'positive': True if lower >= threshold else False if upper < threshold else None})
    return result
