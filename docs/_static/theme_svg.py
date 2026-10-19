"""Apply light/dark display colours to generated documentation SVGs; no data changes.

Run after Graphviz or gnuplot. The embedding page supplies the colour scheme.
"""
import re
import sys
from pathlib import Path

STYLE = """<style id="figure-theme">
:root { fill: var(--ink); --ink: #253444; --line: #526375; --panel: #f4f6f8;
        --grid: #d0d7de; --blue: #0072b2; }
@media (prefers-color-scheme: dark) {
  :root { --ink: #e5e9ef; --line: #a8b5c5; --panel: #252d38; --grid: #465362; --blue: #56b4e9; }
}
</style>"""
COLOURS = {
    "black": "var(--ink)", "#555555": "var(--line)",
    "#666666": "var(--line)", "#465362": "var(--line)",
    "#bac2ca": "var(--grid)", "#f4f6f8": "var(--panel)",
    "rgb(221, 221, 221)": "var(--grid)",
    "#0072b2": "var(--blue)", "rgb(  0, 114, 178)": "var(--blue)",
}


def themed(svg: str) -> str:
    """Keep geometry and numbers intact while adapting neutral colours.

    >>> sample = '<svg><text fill="black">0.11111111</text></svg>'
    >>> '0.11111111' in themed(sample) and 'fill="var(--ink)"' in themed(sample)
    True
    >>> themed(themed(sample)) == themed(sample)
    True
    """
    svg = re.sub(r'<style id="figure-theme">.*?</style>', '', svg, flags=re.S)
    svg = re.sub(r'(<svg\b[^>]*>)', lambda m: m[0] + STYLE, svg, count=1)
    return re.sub(r"(fill|stroke|color)=([\"'])([^\"']+)\2",
                  lambda m: f'{m[1]}={m[2]}{COLOURS.get(m[3].lower(), m[3])}{m[2]}', svg)


if __name__ == "__main__":
    for name in sys.argv[1:]:
        path = Path(name)
        path.write_text(themed(path.read_text()))
