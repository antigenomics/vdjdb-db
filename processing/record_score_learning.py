"""Offline discrete graph learning and unknown-outcome EM for record scores."""
from itertools import combinations

import networkx as nx
import numpy as np
from scipy.special import logsumexp

TARGET = "validated"


def _labels(labels, n):
    values = np.asarray(labels)
    if values.shape != (n,) or not np.isin(values, [-1, 0, 1]).all():
        raise ValueError('Aligned labels must be -1 (unknown), 0 or 1')
    values = values.astype(np.int8)
    if set(values[values >= 0]) != {0, 1}:
        raise ValueError('Observed labels require both classes')
    return values


def learn_structure(data, labels, *, method='hc-bic', max_parents=3):
    """Class-root graphs: HC BIC/AIC, or conditional-MI ARACNE pruning.

    Every feature has the outcome as a parent. ARACNE removes the weakest edge
    in each triangle (lexical tie break), then orients surviving edges by
    descending feature/outcome mutual information. This deterministic orientation
    is a classifier restriction, not causal discovery. HC permits up to three
    total parents by default, including the outcome.
    """
    from sklearn.metrics import mutual_info_score

    y = _labels(labels, len(data))
    names = sorted(data.columns)
    if TARGET in names or not names:
        raise ValueError('Predictors must be nonempty and exclude validated')
    observed = data.loc[y >= 0, names].astype(str).copy()
    y = y[y >= 0]
    required = [(TARGET, n) for n in names]
    if method in ('hc-bic', 'hc-aic'):
        from pgmpy.estimators import ExpertKnowledge, HillClimbSearch

        class OrderedHillClimbSearch(HillClimbSearch):
            def _legal_operations(self, *args, **kwargs):
                return iter(sorted(super()._legal_operations(*args, **kwargs),
                                   key=lambda item: item[0]))

        observed[TARGET] = y.astype(str)
        graph = OrderedHillClimbSearch(observed, use_cache=False).estimate(
            scoring_method='bic-d' if method == 'hc-bic' else 'aic-d',
            expert_knowledge=ExpertKnowledge(required_edges=required,
                                             forbidden_edges=[(n, TARGET) for n in names]),
            max_indegree=max_parents, max_iter=10000, show_progress=False)
        return sorted(graph.edges())
    if method != 'aracne':
        raise ValueError('Unknown graph-learning method')
    information = {(a, b): sum(np.mean(y == cls) * mutual_info_score(
        observed.loc[y == cls, a], observed.loc[y == cls, b]) for cls in (0, 1))
        for a, b in combinations(names, 2)}
    removed = set()
    for triple in combinations(names, 3):
        edges = list(combinations(triple, 2))
        removed.add(min(edges, key=lambda edge: (information[edge], edge)))
    order = sorted(names, key=lambda name: (-mutual_info_score(observed[name], y), name))
    rank = {name: i for i, name in enumerate(order)}
    edges = required.copy()
    parents = dict.fromkeys(names, 1)
    for edge in sorted(information, key=lambda e: (-information[e], e)):
        if edge in removed or information[edge] <= 0:
            continue
        a, b = sorted(edge, key=rank.get)
        if parents[b] < max_parents:
            edges.append((a, b))
            parents[b] += 1
    return sorted(edges)


def _encode(data, states):
    columns = []
    for name, values in states.items():
        mapped = data[name].map({v: i for i, v in enumerate(values)})
        if mapped.isna().any():
            raise ValueError(f'Unknown fitting state in {name}')
        columns.append(mapped.to_numpy(dtype=np.int64))
    return np.column_stack(columns)


def fit_network(data, labels, edges, states, *, prior_mass=10., max_iter=500, tolerance=1e-8):
    """Fit fixed graph by MAP EM, clamping observed outcomes on every iteration.

    Uniform Dirichlet pseudocounts sum to prior_mass for each node. Convergence
    uses the relative change of log likelihood plus the log prior; raw observed
    likelihood and its AIC/BIC are reported separately. A capped fit is explicitly
    marked unconverged. NumPy bincount reduces weighted sufficient statistics.
    """
    y = _labels(labels, len(data))
    states = {n: list(states[n]) for n in sorted(states)}
    if set(data.columns) != set(states) or TARGET in states or not states:
        raise ValueError('Predictor columns and state declarations must match')
    if any(not v or len(v) != len(set(v)) for v in states.values()):
        raise ValueError('State declarations must be nonempty and unique')
    if not np.isfinite(prior_mass) or prior_mass <= 0 or max_iter < 1 or tolerance <= 0:
        raise ValueError('Positive fitting limits and prior mass required')
    graph = nx.DiGraph()
    graph.add_nodes_from([TARGET, *states])
    graph.add_edges_from(edges)
    if set(graph.nodes) != {TARGET, *states}:
        raise ValueError('Graph contains undeclared nodes')
    if not nx.is_directed_acyclic_graph(graph):
        raise ValueError('Graph must be acyclic')
    if graph.in_degree(TARGET) or any(not graph.has_edge(TARGET, n) for n in states):
        raise ValueError('Class must be a root parent of every predictor')
    x = _encode(data, states)
    n_labeled, n_unlabeled = int((y >= 0).sum()), int((y < 0).sum())
    patterns, multiplicity = np.unique(np.column_stack([y, x]), axis=0, return_counts=True)
    y, x = patterns[:, 0], patterns[:, 1:]
    names = list(states)
    positions = {n: i for i, n in enumerate(names)}
    cards = {n: len(v) for n, v in states.items()}
    families = []
    for name in names:
        parents = sorted(p for p in graph.predecessors(name) if p != TARGET)
        q = int(np.prod([cards[p] for p in parents], dtype=np.int64))
        codes = np.zeros(len(x), dtype=np.int64)
        for parent in parents:
            codes = codes * cards[parent] + x[:, positions[parent]]
        index = codes * cards[name] + x[:, positions[name]]
        families.append((name, parents, q, index))
    # Supervised initialization anchors outcome identity; no random pseudo-labels.
    weights = np.column_stack([y == 0, y == 1]).astype(float) * multiplicity[:, None]

    def maximize(w):
        prior = w.sum(axis=0) + prior_mass / 2
        prior /= prior.sum()
        tables = []
        for name, _parents, q, index in families:
            shape = (2, q, cards[name])
            counts = np.array([np.bincount(index, weights=w[:, cls], minlength=q * cards[name])
                               for cls in (0, 1)]).reshape(shape)
            counts += prior_mass / counts.size
            tables.append(counts / counts.sum(axis=2, keepdims=True))
        return prior, tables

    def expectation(prior, tables):
        joint = np.broadcast_to(np.log(prior), (len(x), 2)).copy()
        penalty = float(prior_mass / 2 * np.log(prior).sum())
        for (_, _, _, index), table in zip(families, tables, strict=True):
            joint += np.log(table.reshape(2, -1)[:, index].T)
            penalty += prior_mass / table.size * np.log(table).sum()
        normalizer = logsumexp(joint, axis=1)
        likelihood = float(normalizer[y < 0] @ multiplicity[y < 0] +
                           joint[y >= 0, y[y >= 0]] @ multiplicity[y >= 0])
        w = np.exp(joint - normalizer[:, None])
        w[y >= 0] = np.column_stack([y[y >= 0] == 0, y[y >= 0] == 1])
        return w * multiplicity[:, None], likelihood, float(likelihood + penalty)

    prior, tables = maximize(weights)
    trace = []
    converged = False
    for iteration in range(max_iter + 1):
        weights, likelihood, objective = expectation(prior, tables)
        trace.append(objective)
        if len(trace) > 1 and abs(trace[-1] - trace[-2]) <= tolerance * (1 + abs(trace[-2])):
            converged = True
            break
        if iteration == max_iter:
            break
        prior, tables = maximize(weights)
    parameters = 1 + sum(2 * q * (cards[name] - 1) for name, _, q, _ in families)
    return {'format': 'vdjdb-learned-bn-1', 'states': states, 'edges': sorted(map(list, graph.edges())),
            'prior': prior.tolist(), 'factors': [
                {'name': name, 'parents': parents, 'values': table.tolist()}
                for (name, parents, _, _), table in zip(families, tables, strict=True)],
            'n_labeled': n_labeled, 'n_unlabeled': n_unlabeled, 'n_patterns': len(x),
            'prior_mass': prior_mass, 'iterations': iteration, 'converged': converged,
            'objective_trace': trace, 'log_likelihood': likelihood, 'parameters': parameters,
            'aic': 2 * parameters - 2 * likelihood,
            'bic': float(parameters * np.log(len(data)) - 2 * likelihood),
            'missing_fit': 'explicit reporting state', 'missing_predict': 'exact marginalization'}
