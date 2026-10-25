"""Review artwork from existing TCREmp memberships and inferred junction boundaries.
Run `python images/render_splash.py` with rsvg-convert installed.
Reads the committed display selection; no scientific analysis is run.
"""

import csv
import html
import math
import subprocess
from pathlib import Path

OUT = Path(__file__).resolve().parent
with (OUT / "vdjdb-splash-members.tsv").open() as source:
    rows = list(csv.DictReader(source, delimiter="\t"))
for row in rows:
    for key in ("v_end_nt", "d_start_nt", "d_end_nt", "j_start_nt"):
        row[key] = int(row[key])
    assert (
        0
        <= row["v_end_nt"]
        <= row["d_start_nt"]
        < row["d_end_nt"]
        <= row["j_start_nt"]
        <= 3 * len(row["cdr3"])
    )
peptides = list(dict.fromkeys(row["antigen.epitope"] for row in rows))
assert len(rows) == 105 and len(peptides) == 5
assert len({row["cid"] for row in rows}) == 15

P = {
    "light": dict(
        bg="#f8fafc",
        ink="#182b40",
        muted="#667789",
        line="#c4cfda",
        panel="#edf2f6",
        V="#0072b2",
        D="#aa6200",
        J="#9960a8",
        N="#6c7c8b",
        branch="#147d70",
    ),
    "dark": dict(
        bg="#101c2b",
        ink="#edf4fb",
        muted="#a4b6c8",
        line="#3b5268",
        panel="#1a2e43",
        V="#63bcec",
        D="#f1bd64",
        J="#d2a3e5",
        N="#a4b6c8",
        branch="#58c9b3",
    ),
}
C = (1000, 940)


def xy(r, deg):
    t = math.radians(deg)
    return C[0] + r * math.cos(t), C[1] + r * math.sin(t)


def label(x, y, s, size, color, anchor="middle", weight=400):
    return (
        f'<text x="{x:.2f}" y="{y:.2f}" text-anchor="{anchor}" fill="{color}"'
        f' font-family="Helvetica,Arial,sans-serif" font-size="{size}"'
        f' font-weight="{weight}">{html.escape(s)}</text>'
    )


def edge(r1, t1, r2, t2, color, width=1, opacity=1):
    x1, y1 = xy(r1, t1)
    x2, y2 = xy(r2, t2)
    a, b = xy((r1 + r2) / 2, t1)
    c, d = xy((r1 + r2) / 2, t2)
    return (
        f'<path d="M{x1},{y1} C{a},{b} {c},{d} {x2},{y2}" fill="none" stroke="{color}"'
        f' stroke-width="{width}" opacity="{opacity}"/>'
    )


def symbol(x, y, kind, color, background):
    if kind == "epitope":
        return f'<path d="M{x},{y - 8} l8,8 -8,8 -8,-8 Z" fill="{color}"/>'
    return (
        f'<circle cx="{x}" cy="{y}" r="8" fill="{background}"'
        f' stroke="{color}" stroke-width="2.5"/>'
    )


for theme, p in P.items():
    assert p["branch"] not in {p[k] for k in ("V", "D", "J", "N")}
    defs = []
    parts = []
    epitope_labels = []
    parts.append(f'<rect width="2000" height="1900" fill="{p["bg"]}"/>')
    parts.append(label(1000, 73, "VDJdb · a map of receptor recognition", 42, p["ink"], weight=600))
    parts.append(label(1000, 116, "Epitope → paratope motif → CDR3 sequence", 24, p["muted"]))
    # Decorative rings also expose the hierarchy without implying measured distances.
    for r in [258, 430, 598]:
        parts.append(
            f'<circle cx="1000" cy="940" r="{r}" fill="none" stroke="{p["line"]}"'
            ' stroke-width="1" opacity=".25"/>'
        )
    for ei, ep in enumerate(peptides):
        center = -90 + ei * 72
        color = p["branch"]
        parts.append(edge(91, center, 232, center, color, 3))
        x, y = xy(265, center)
        epitope_labels.append(
            f'<rect x="{x - 104}" y="{y - 23}" width="208" height="50" rx="25" fill="{p["bg"]}"'
            f' stroke="{color}" stroke-width="1.5"/>'
        )
        epitope_labels.append(symbol(x - 84, y + 1, "epitope", color, p["bg"]))
        epitope_labels.append(label(x + 12, y + 7, ep, 23, p["ink"], weight=600))
        group = [r for r in rows if r["antigen.epitope"] == ep]
        cids = list(dict.fromkeys(r["cid"] for r in group))
        for mi, cid in enumerate(cids):
            ma = center + (mi - 1) * 21
            parts.append(edge(292, center, 430, ma, color, 2, 0.8))
            mx, my = xy(430, ma)
            parts.append(symbol(mx, my, "paratope", color, p["bg"]))
            lx, ly = xy(401, ma)
            parts.append(label(lx, ly + 5, cid.rsplit(".", 1)[-1], 15, color))
            sample = [r for r in group if r["cid"] == cid]
            for si, row in enumerate(sample):
                angle = ma + (si - 3) * 2.6
                parts.append(edge(438, ma, 598, angle, color, 1, 0.55))
                sx, sy = xy(598, angle)
                parts.append(f'<circle cx="{sx}" cy="{sy}" r="2.7" fill="{color}"/>')
                # Each residue is split into thirds at mixed codons using stored nt boundaries.
                seq = row["cdr3"]
                left = math.cos(math.radians(angle)) < 0
                rotate = angle + 180 if left else angle
                x, y = xy(612, angle)
                advance = 9.6
                start = -len(seq) * advance if left else 0
                letters = []
                for j, ch in enumerate(seq):
                    kinds = []
                    for nt in range(j * 3, j * 3 + 3):
                        kinds.append(
                            "V"
                            if nt < row["v_end_nt"]
                            else "J"
                            if nt >= row["j_start_nt"]
                            else "D"
                            if row["d_start_nt"] <= nt < row["d_end_nt"]
                            else "N"
                        )
                    if len(set(kinds)) == 1:
                        fill = p[kinds[0]]
                    else:
                        ident = f"codon{ei}_{mi}_{si}_{j}"
                        stops = []
                        for k, kind in enumerate(kinds):
                            stops.extend(
                                [
                                    f'<stop offset="{k / 3:.6f}" stop-color="{p[kind]}"/>',
                                    f'<stop offset="{(k + 1) / 3:.6f}" stop-color="{p[kind]}"/>',
                                ]
                            )
                        defs.append(f'<linearGradient id="{ident}">' + "".join(stops) + "</linearGradient>")
                        fill = f"url(#{ident})"
                    letters.append(
                        f'<text x="{start + j * advance}" y="5"'
                        f' font-family="Menlo,DejaVu Sans Mono,monospace" font-size="16"'
                        f' fill="{fill}">{ch}</text>'
                    )
                parts.append(
                    f'<g transform="translate({x},{y}) rotate({rotate})">' + "".join(letters) + "</g>"
                )
    parts.extend(epitope_labels)
    parts.append(
        f'<circle cx="1000" cy="940" r="91" fill="{p["panel"]}" stroke="{p["V"]}" stroke-width="2.5"/>'
    )
    parts.append(label(1000, 940, "VDJdb", 36, p["ink"], weight=600))
    parts.append(label(1000, 969, "human TRB", 15, p["muted"]))
    for x, kind, title in [(810, "epitope", "Epitope"), (1045, "paratope", "Paratope motif (cid)")]:
        parts.append(symbol(x, 1737, kind, p["branch"], p["bg"]))
        parts.append(label(x + 18, 1744, title, 20, p["ink"], "start"))
    for x, k, t in [
        (610, "V", "V segment"),
        (845, "D", "D segment"),
        (1080, "J", "J segment"),
        (1315, "N", "Junction gap"),
    ]:
        parts.append(f'<circle cx="{x}" cy="1782" r="6" fill="{p[k]}"/>')
        parts.append(label(x + 16, 1789, t, 20, p["ink"], "start"))
    parts.append(
        label(1000, 1831, "105 displayed sequences · 15 TCREmp paratope motifs · 5 epitopes", 20, p["muted"])
    )
    parts.append(
        label(
            1000,
            1865,
            "Anchor-inclusive CDR3s; inferred V/D/J boundaries. Mixed codons use split colours. "
            "Branches show membership, not distance.",
            17,
            p["muted"],
        )
    )
    svg = (
        '<svg xmlns="http://www.w3.org/2000/svg" width="2000" height="1900" viewBox="0 0 2000 1900"'
        ' role="img"><title>VDJdb paratope motif atlas</title><desc>Five epitope groups,'
        " three existing TCREmp motifs per epitope, and seven sequences per motif, "
        "coloured using stored nucleotide boundaries. Teal diamonds mark epitopes; "
        "teal rings mark paratope motifs.</desc><defs>"
        + "".join(defs)
        + "</defs>"
        + "".join(parts)
        + "</svg>"
    )
    path = OUT / f"vdjdb-splash-{theme}.svg"
    path.write_text(svg)
    if theme == "light":
        subprocess.run(["rsvg-convert", "-o", str(OUT / "vdjdb-splash.png"), str(path)], check=True)
assert len(rows) == 105
print(OUT)
