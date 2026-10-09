"""Execute the workflow selector with controlled Git inputs and external commands."""
import subprocess
from pathlib import Path

import pytest
import yaml

WORKFLOW = Path(__file__).parents[2] / '.github/workflows/assessment.yml'


def script():
    workflow = yaml.load(WORKFLOW.read_text(), Loader=yaml.BaseLoader)
    step = next(s for s in workflow['jobs']['assess']['steps'] if 'env' in s)
    return step['run'].split("<<'PY'\n", 1)[1].rsplit('\nPY', 1)[0]


@pytest.mark.parametrize('selected', ['', 'chunks/missing.tsv', '../outside.tsv', '$(bad)'])
def test_invalid_selection_stops_before_any_assessment(tmp_path, monkeypatch, selected):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv('SELECTED_CHUNKS', selected)
    monkeypatch.setattr(subprocess, 'check_output', lambda *a, **kw: 'chunks/valid.tsv\0')
    commands = []
    monkeypatch.setattr(subprocess, 'run', lambda *a, **kw: commands.append(a))
    with pytest.raises(SystemExit, match='Select existing tracked chunks/'):
        exec(script(), {})
    assert not commands


@pytest.mark.parametrize('predict', ['true', 'false'])
def test_selection_is_explicit_deduplicated_and_prediction_is_opt_in(tmp_path, monkeypatch, predict):
    monkeypatch.chdir(tmp_path)
    Path('chunks').mkdir()
    Path('chunks/valid.tsv').touch()
    monkeypatch.setenv('SELECTED_CHUNKS', 'chunks/valid.tsv chunks/valid.tsv')
    monkeypatch.setenv('PREDICT', predict)
    monkeypatch.setenv('GITHUB_STEP_SUMMARY', str(tmp_path / 'summary.md'))
    monkeypatch.setattr(subprocess, 'check_output', lambda *a, **kw: 'chunks/valid.tsv\0')
    commands = []

    def run(args, **kwargs):
        assert kwargs == {'check': True}
        commands.append(args)
        if 'submission' in args:
            out = Path(args[args.index('--out') + 1])
            out.parent.mkdir(parents=True)
            out.write_text('selected submission')

    monkeypatch.setattr(subprocess, 'run', run)
    exec(script(), {})
    assert ('epitope-reference' in [c[3] for c in commands]) == (predict == 'true')
    assessment = next(c for c in commands if c[3] == 'assess-epitopes')
    assert assessment.count('chunks/valid.tsv') == 1
    assert ('--pmhc-reference' in assessment) == (predict == 'true')
    assert not {'build', 'motifs', 'summary', 'release'} & {c[3] for c in commands}
    assert 'selected submission' in (tmp_path / 'summary.md').read_text()
