"""A missing QC report must not be presented as a clean submission."""
from pathlib import Path

import pytest
import yaml

WORKFLOW = Path(__file__).parents[2] / '.github/workflows/chunk-check.yml'


@pytest.mark.parametrize('report,selected', [(None, False), (None, True),
                                             ('', True), ('a.tsv\t2\tduplicate\t\n', True)])
def test_pr_comment_distinguishes_missing_clean_and_failed_qc(
    tmp_path, monkeypatch, capsys, report, selected,
):
    workflow = yaml.load(WORKFLOW.read_text(), Loader=yaml.BaseLoader)
    step = next(s for s in workflow['jobs']['check']['steps']
                if s.get('name') == 'Summarise for the pull request')
    script = step['run'].split('\n', 1)[1].split('\nPY', 1)[0]
    monkeypatch.chdir(tmp_path)
    if selected:
        Path('changed.txt').write_text('')
    if report is not None:
        Path('qc-report.tsv').write_text('file\trow\tcode\tdetail\n' + report)
    exec(script, {})
    comment = capsys.readouterr().out
    if report is None:
        assert 'QC did not complete' in comment
        assert 'No QC findings' not in comment and 'was validated' not in comment
    elif report:
        assert '1 QC finding' in comment and '`duplicate`' in comment
    else:
        assert 'No QC findings' in comment
