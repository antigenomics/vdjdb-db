"""Small independent probability checks for fresh network fitting."""
import json

import pytest

pytest.importorskip('pgmpy')
import numpy as np
import pandas as pd

from vdjdb.score.learning import fit_network, learn_structure, predict_network


def test_supervised_counts_and_missing_marginal():
    data = pd.DataFrame({'x': ['0', '0', '1', '1'], 'z': ['0', '1', '0', '1']})
    model = fit_network(data, [0, 0, 1, 1], [('validated', 'x'), ('validated', 'z')],
                        {'x': ['0', '1'], 'z': ['0', '1']}, prior_mass=4)
    actual = predict_network(json.loads(json.dumps(model)), data)
    assert actual == pytest.approx([.25, .25, .75, .75])
    missing = pd.DataFrame({'x': ['?', '1'], 'z': ['?', '?']})
    assert predict_network(model, missing) == pytest.approx([.5, .75])
    assert model['converged']
    assert model['n_labeled'] == 4


def test_unlabeled_em_and_objective():
    data = pd.DataFrame({'x': ['0'] * 4 + ['1'] * 16})
    labels = [0, 0, 0, 1] + [1, 1, 1, 0] + [-1] * 12
    model = fit_network(data, labels, [('validated', 'x')], {'x': ['0', '1']},
                        prior_mass=1, max_iter=1000)
    assert model['n_labeled'] == 8 and model['n_unlabeled'] == 12
    assert model['converged']
    assert np.min(np.diff(model['objective_trace'])) >= -1e-8
    assert np.isfinite(predict_network(model, data)).all()
    assert model == fit_network(data, labels, [('validated', 'x')], {'x': ['0', '1']},
                                prior_mass=1, max_iter=1000)
    with pytest.raises(ValueError, match='labels'):
        fit_network(data, [2] * 20, [('validated', 'x')], {'x': ['0', '1']})


@pytest.mark.parametrize('method', ['hc-bic', 'hc-aic', 'aracne'])
def test_graph_learning_and_cycle_rejection(method):
    data = pd.DataFrame({'x': ['0', '0', '1', '1'] * 10, 'z': ['0', '1', '0', '1'] * 10})
    labels = [0, 0, 1, 1] * 10
    edges = learn_structure(data, labels, method=method)
    assert ('validated', 'x') in edges and ('validated', 'z') in edges
    assert all(v != 'validated' for u, v in edges)
    assert edges == learn_structure(data, labels, method=method)
    with pytest.raises(ValueError, match='acyclic'):
        fit_network(data, labels, [*edges, ('x', 'z'), ('z', 'x')],
                    {'x': ['0', '1'], 'z': ['0', '1']})


def test_dependent_missing_posterior_matches_enumeration():
    from itertools import product
    data = pd.DataFrame(list(product(['0', '1'], repeat=3)) * 3, columns=['a', 'b', 'c'])
    labels = [int(a == b or c == '1') for a, b, c in data.itertuples(index=False, name=None)]
    model = fit_network(data, labels, [('validated', n) for n in data] + [('c', 'a'), ('a', 'b')],
                        {n: ['0', '1'] for n in data})
    query = pd.DataFrame({'a': ['?', '0'], 'b': ['1', '?'], 'c': ['0', '?']})
    expected = []
    for row in query.to_dict('records'):
        sums = [0., 0.]
        for cls in (0, 1):
            for values in product(range(2), repeat=3):
                assignment = dict(zip(['a', 'b', 'c'], values, strict=True))
                if any(value != '?' and assignment[name] != int(value) for name, value in row.items()):
                    continue
                probability = model['prior'][cls]
                for factor in model['factors']:
                    index = 0
                    for parent in factor['parents']:
                        index = index * 2 + assignment[parent]
                    probability *= factor['values'][cls][index][assignment[factor['name']]]
                sums[cls] += probability
        expected.append(sums[1] / sum(sums))
    assert predict_network(model, query) == pytest.approx(expected)


def test_explicit_epitope_models_require_context_and_preserve_other_rows():
    data = pd.DataFrame({'x': ['0', '0', '1', '1']})
    edges, states = [('validated', 'x')], {'x': ['0', '1']}
    global_model = fit_network(data, [0, 0, 1, 1], edges, states)
    local_model = fit_network(data, [1, 1, 0, 0], edges, states)
    model = {**global_model, 'epitope_models': {'e': local_model}}
    with pytest.raises(ValueError, match='epitope'):
        predict_network(model, data)
    original = predict_network(global_model, data)
    local = predict_network(local_model, data)
    actual = predict_network(model, data, epitopes=['e', 'other', 'e', 'other'])
    assert actual == pytest.approx([local[0], original[1], local[2], original[3]])


def test_hill_climb_ties_are_stable_across_hash_seeds():
    import os
    import subprocess
    import sys
    code = """
import json
import pandas as pd
from vdjdb.score.learning import learn_structure
x = ['0', '0', '1', '1'] * 20
data = pd.DataFrame({name: x for name in ['a', 'b', 'c']})
print(json.dumps(learn_structure(data, [0, 1] * 40)))
"""
    outputs = [subprocess.check_output(
        [sys.executable, '-c', code],
        env={**os.environ, 'PYTHONHASHSEED': str(seed)}, text=True)
        for seed in (1, 2, 3, 4)]
    assert len(set(outputs)) == 1
