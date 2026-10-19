"""Small discrete Bayesian evidence classifier with exact missing-value inference.

Install the ``confidence`` extra for fitting. Missing fields are marginalized at prediction;
no evidence returns the fitted network prior. Fitted artifacts contain no receptor identities.
"""
from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from .evidence import COUNT_FEATURES, FEATURES, categorical

STATES = {name: (['?', '0', '1', '2', '3+'] if name in COUNT_FEATURES else
                 ['?', 'low', 'mid', 'high'] if name == 'iptm' else
                 ['?', '0', '1', '2', '3'] if name == 'legacy_score' else ['?', '0', '1'])
          for name in FEATURES}


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


def predict(artifact: dict, rows: Sequence[Mapping[str, str | float | None]]) -> list[float]:
    """Score categorical evidence or uncertain binary motif membership.

    A numeric motif value q mixes hard nonmember/member posteriors with weights
    (1-q, q). Multiple uncertain memberships use an independent mixing distribution.
    These are caller-supplied membership probabilities, not binding probabilities.
    """
    if any(isinstance(value, (int, float)) for row in rows for name, value in row.items()
           if name in FEATURES and name.endswith(('tcrnet', 'tcremp'))):
        expanded, weights = [], []
        for row in rows:
            mixtures = [(dict(row), 1.)]
            for name, value in row.items():
                if name not in FEATURES or not name.endswith(('tcrnet', 'tcremp')):
                    continue
                if not isinstance(value, (int, float)):
                    continue
                if not math.isfinite(value) or not 0 <= value <= 1:
                    raise ValueError(f'Invalid membership probability for {name}')
                mixtures = [({**r, name: str(state)}, weight * probability)
                            for r, weight in mixtures
                            for state, probability in ((0, 1-value), (1, value)) if probability > 0]
            expanded.extend(r for r, _ in mixtures)
            weights.append([weight for _, weight in mixtures])
        probabilities = iter(predict(artifact, expanded))
        return [math.fsum(weight * next(probabilities) for weight in group) for group in weights]

    from pgmpy.factors.discrete import TabularCPD
    from pgmpy.inference import VariableElimination
    from pgmpy.models import DiscreteBayesianNetwork

    if artifact.get('format') != 'vdjdb-evidence-tan-1' or artifact.get('features') != FEATURES:
        raise ValueError('Unsupported classifier artifact')
    model = DiscreteBayesianNetwork(artifact['edges'])
    model.add_nodes_from(cpd['variable'] for cpd in artifact['cpds'])
    for cpd in artifact['cpds']:
        variables, card = cpd['variables'], cpd['cardinality']
        model.add_cpds(TabularCPD(cpd['variable'], card[0], cpd['values'],
                                 evidence=variables[1:] or None,
                                 evidence_card=card[1:] or None, state_names=cpd['states']))
    model.check_model()
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
    infer = VariableElimination(model)
    # Within-call deduplication is not a persisted feature or inference cache.
    values = {}
    for pattern in sorted(set(patterns)):
        result = infer.query(['validated'], evidence=dict(pattern), show_progress=False)
        values[pattern] = float(result.get_value(validated='1'))
    return [values[p] for p in patterns]


def predict_marginals(artifact: dict, rows: Sequence[Mapping[str, str | float | None]]
                      ) -> list[dict[str, float]]:
    """Full-interaction validation from paired, TRA-only and TRB-only evidence.

    Single-chain views omit paired legacy score, donor total and structure confidence.
    They condition on that chain's evidence, not on a single-chain binding outcome.
    """
    views = [list(rows), [{k: v for k, v in row.items() if k.startswith('a_')} for row in rows],
             [{k: v for k, v in row.items() if k.startswith('b_')} for row in rows]]
    values = predict(artifact, [row for view in views for row in view])
    return [{name: values[j * len(rows) + i] for j, name in enumerate(('paired', 'TRA', 'TRB'))}
            for i in range(len(rows))]


def predict_evidence(artifact: dict, frame) -> list[dict[str, float]]:
    """Score numeric evidence using the saved grades and one network for all three views.

    Missing feature columns are unobserved. Extra columns, such as receptor identifiers,
    are preserved by the caller and are not predictors. No epitope name enters inference.
    """
    motifs = [name for name in FEATURES if name.endswith(('tcrnet', 'tcremp'))
              and name in frame.columns]
    names = [name for name in FEATURES if name in frame.columns and name not in motifs]
    rows = categorical(frame.select(names), iptm_bounds=tuple(
        artifact.get('iptm_bounds', (.5, .8)))).to_dicts() if names else [{}] * frame.height
    if motifs:
        import polars as pl

        membership = frame.select(pl.col(motifs).cast(pl.Float64).fill_nan(None)).to_dicts()
        rows = [{**row, **soft} for row, soft in zip(rows, membership, strict=True)]
    return predict_marginals(artifact, rows)


def score_records(artifact: dict, frame) -> list[dict[str, float | bool]]:
    """Record acceptance score, retaining the unchanged assay and chain probabilities.

    Positive support must concern both chains, or the exact pair (legacy annotation,
    donor support or a matching predicted structure). Unsupported pairs score zero;
    this is an evidence requirement, not a claim of zero binding probability.
    Use one global threshold, 0.5. Missing evidence cannot establish pair support.
    """
    import polars as pl

    probabilities = predict_evidence(artifact, frame)
    def value(name):
        return pl.col(name).cast(pl.Float64).fill_nan(None) if name in frame.columns else pl.lit(None)

    chains = [pl.any_horizontal((value(name) > 0 if name.endswith('studies') else value(name) == 1)
                               for name in FEATURES if name.startswith(prefix))
              for prefix in ('a_', 'b_')]
    direct = (value('legacy_score').is_not_null() | (value('donors') > 0) |
              value('iptm').is_not_null())
    supported = frame.with_columns((direct | (chains[0] & chains[1])).fill_null(False)
                                   .alias('_pair_supported'))['_pair_supported']
    return [{**p, 'pair_supported': bool(support),
             'record_score': p['paired'] if support else 0.,
             'positive': bool(support and p['paired'] >= .5)}
            for p, support in zip(probabilities, supported, strict=True)]
