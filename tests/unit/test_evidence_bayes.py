"""Optional pgmpy contract, run in the pinned confidence environment."""
import json

import pytest

pytest.importorskip('pgmpy')
from vdjdb.score.bayes import fit, predict, predict_marginals


def test_exact_missing_inference_roundtrip_and_deterministic_tree():
    rows = [{'a_studies': '2' if i % 2 else '0', 'b_tcrnet': '1' if i % 3 else '0'}
            for i in range(24)]
    model = fit(rows, [i % 2 for i in range(24)])
    restored = json.loads(json.dumps(model))
    p = predict(restored, [{}, {'a_studies': None}, {'a_studies': '2'}, {'a_studies': '0'}])
    assert p[0] == pytest.approx(.5)
    assert p[1] == p[0] and p[2] > p[0] > p[3]
    assert model == fit(rows, [i % 2 for i in range(24)])
    views = predict_marginals(model, [rows[1], {}])
    assert views[0]['TRA'] == pytest.approx(predict(model, [{'a_studies': rows[1]['a_studies']}])[0])
    assert views[0]['TRB'] == pytest.approx(predict(model, [{'b_tcrnet': rows[1]['b_tcrnet']}])[0])
    assert views[1] == {'paired': p[0], 'TRA': p[0], 'TRB': p[0]}
    with pytest.raises(ValueError):
        predict(model, [{'padj': '0'}])
    with pytest.raises(ValueError):
        predict(model, [{'iptm': 'ultra'}])


def test_unobserved_legal_evidence_states_do_not_create_support():
    rows = [{'a_studies': '2' if i % 2 else '0'} for i in range(24)]
    model = fit(rows, [int(i % 4 == 0) for i in range(24)])
    assert predict(model, [{'iptm': 'high', 'donors': '3+'}])[0] == pytest.approx(predict(model, [{}])[0])


def test_naive_network_keeps_exact_missing_inference_with_fewer_edges():
    rows = [{'a_studies': str(i % 2), 'b_tcrnet': str(int(i % 3 == 0))}
            for i in range(36)]
    model = fit(rows, [i % 2 for i in range(36)], structure='naive')
    assert len(model['edges']) == 9
    assert all(parent == 'validated' for parent, _ in model['edges'])
    assert predict(model, [{}, {'a_studies': None}]) == pytest.approx([.5, .5])
    assert predict(model, [{'a_studies': '1'}])[0] > .5
    with pytest.raises(ValueError):
        fit(rows, [i % 2 for i in range(36)], structure='unknown')


def test_unpredicted_structure_is_marginalized_in_observed_chain_context():
    rows = [{'a_studies': '1', 'iptm': 'high' if i % 2 else 'low'} for i in range(24)]
    for structure in ['tan', 'naive']:
        model = fit(rows, [i % 2 for i in range(24)], structure=structure)
        scores = predict(model, [
            {'a_studies': '1'}, {'a_studies': '1', 'iptm': None},
            {'a_studies': '1', 'iptm': '?'}, {'a_studies': '1', 'iptm': 'high'},
        ])
        assert scores[:3] == pytest.approx([scores[0]] * 3)
        assert scores[3] > scores[0]
