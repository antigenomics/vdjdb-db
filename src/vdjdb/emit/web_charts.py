"""Static legacy web scatter charts from shipped motif coordinates; no re-embedding."""
import re
from pathlib import Path

from ..io.text import read_csv, resolve_path


def write(motifs: Path, out: Path) -> int:
    import plotly.graph_objects as go
    from plotly.offline import get_plotlyjs

    out.mkdir(parents=True, exist_ok=True)
    (out / 'plotly.min.js').write_text(get_plotlyjs())
    # This directory is the generated chart artifact, never a source/input directory.
    for previous in out.rglob('*.html'):
        previous.unlink()
    count = 0
    names = set()
    for method, source in [('tcrnet', 'cluster_members.txt'),
                           ('tcremp', 'cluster_members_tcremp.txt')]:
        if method == 'tcremp' and not resolve_path(motifs / source).exists():
            continue
        data = read_csv(motifs / source, separator='\t', quote_char=None,
            columns=['species', 'gene', 'antigen.epitope', 'cid', 'cdr3aa', 'x', 'y'])
        folder = out / method
        folder.mkdir(exist_ok=True)
        for (species, gene, epitope), group in sorted(data.group_by('species', 'gene', 'antigen.epitope')):
            parts = (species, gene, epitope) if method == 'tcrnet' else (species, epitope, gene)
            name = '_'.join(re.sub('[^a-z0-9]', '', value.lower()) for value in parts) + '.html'
            target = folder / name
            if target in names:
                raise ValueError('motif chart URL collision')
            names.add(target)
            figure = go.Figure()
            for (cid,), cluster in sorted(group.group_by('cid')):
                figure.add_trace(go.Scattergl(x=cluster['x'].to_list(), y=cluster['y'].to_list(),
                    text=cluster['cdr3aa'].to_list(), name=cid, mode='markers',
                    marker={'size': 5}, hovertemplate='%{text}<extra>%{fullData.name}</extra>'))
            figure.update_layout(title=epitope, template='plotly_white', autosize=True,
                                 margin={'l': 30, 'r': 10, 't': 10, 'b': 30})
            figure.write_html(target, include_plotlyjs='../plotly.min.js', div_id='motif-plot')
            count += 1
    return count
