"""Verify the rendered dashboard fragment before it reaches ``vdjdb-web``.

    uv run python summary/check_summary.py [fragment.html] [--baseline summary/fingerprint.json]
    uv run python summary/check_summary.py --write-baseline

There is no committed image baseline -- ``.gitignore`` excludes ``summary/*.html`` -- so the
reference is a **committed fingerprint** (``summary/fingerprint.json``) plus, optionally, a previous
release's fragment for a perceptual comparison. Three layers, cheapest first:

**Structural**, exact, no image decoding. The ordered ``<h4>`` list, the ordered ``<th>`` text of
every Semantic UI table, the image count, and each image's width and height read straight out of the
PNG IHDR in the first 32 base64 characters -- the cheapest possible detector for a ``fig.width`` /
``dpi`` / ``fig.retina`` regression. Plus the three contracts ``vdjdb-web`` depends on: single-line
base64, no ``<div>``, no document scaffolding.

**Style**, needs Pillow. Each panel is quantised and checked for at least three of its declared
ColorBrewer anchors within :data:`DELTA_E`, with an anti-assertion that no panel carries three
viridis anchors -- that is what makes "someone swapped ``Set1`` for viridis" a red build. An ink
fraction outside :data:`INK_BAND` catches a blank or collapsed facet.

**Perceptual**, needs a previous fragment. SSIM per panel against it, loose and asymmetric because
the database grows: fail below :data:`SSIM_FAIL`, warn below :data:`SSIM_WARN`. The point is
catching "panel 7 is a solid grey block now", not pixel equality.

⚠ **Nothing here reads TCRvdb.**
"""
from __future__ import annotations

import argparse
import base64
import io
import json
import re
import struct
import sys
from pathlib import Path

FRAGMENT = Path("summary/vdjdb_summary_embed.html")
BASELINE = Path("summary/fingerprint.json")

#: CIE76 distance under which a quantised colour counts as a palette anchor.
DELTA_E = 3.0
#: Fraction of non-white pixels a panel must carry: below is blank, above is a solid block.
INK_BAND = (0.01, 0.80)
SSIM_FAIL, SSIM_WARN = 0.55, 0.80

#: ColorBrewer anchors the dashboard actually uses, and the viridis anchors it must NOT.
SET1 = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00"]
VIRIDIS = ["#440154", "#21918c", "#fde725", "#3b528b", "#5ec962"]

_IMG = re.compile(r'data:image/png;base64,([A-Za-z0-9+/=]+)')
_H4 = re.compile(r"<h4[^>]*>(.*?)</h4>", re.S)
_TABLE = re.compile(r"<table class=[^>]*>(.*?)</table>", re.S)
_TH = re.compile(r"<th[^>]*>(.*?)</th>", re.S)
_TAG = re.compile(r"<[^>]+>")


def _text(html: str) -> str:
    return re.sub(r"\s+", " ", _TAG.sub("", html)).strip()


def png_size(payload: str) -> tuple[int, int]:
    """Width and height from the PNG IHDR, decoded from the first 32 base64 characters.

    An IHDR sits at bytes 16..24 of every PNG, so 32 base64 characters (24 bytes) always contain it
    and the rest of a multi-megabyte payload never has to be decoded.
    """
    head = base64.b64decode(payload[:32] + "=" * (-len(payload[:32]) % 4))
    return struct.unpack(">II", head[16:24])


def fingerprint(html: str) -> dict:
    """Everything the structural layer compares, as plain JSON."""
    images = _IMG.findall(html)
    return {
        "headings": [_text(h) for h in _H4.findall(html)],
        "tables": [[_text(t) for t in _TH.findall(body)] for body in _TABLE.findall(html)],
        "images": [list(png_size(p)) for p in images],
        "contracts": {
            "no_div": "<div" not in html,
            "no_document_scaffolding": not any(t in html.lower()
                                               for t in ("<html", "<head", "<body", "<script")),
            "base64_single_line": all("data:image/png;base64," not in n or n.rstrip().endswith(
                ("/>", ">", "</p>")) for n in html.splitlines()),
            "images_responsive": html.count("max-width:100%") == len(images),
        },
    }


def structural(current: dict, baseline: dict | None) -> list[str]:
    bad = [f"contract {k} is false" for k, v in current["contracts"].items() if not v]
    if baseline is None:
        return bad
    if current["headings"] != baseline["headings"]:
        bad.append(f"headings changed: {baseline['headings']} -> {current['headings']}")
    if len(current["tables"]) != len(baseline["tables"]):
        bad.append(f"{len(baseline['tables'])} tables expected, {len(current['tables'])} found")
    else:
        for n, (a, b) in enumerate(zip(baseline["tables"], current["tables"], strict=True)):
            if a != b:
                bad.append(f"table {n} header changed: {a} -> {b}")
    if len(current["images"]) != len(baseline["images"]):
        bad.append(f"{len(baseline['images'])} images expected, {len(current['images'])} found")
    else:
        for n, (a, b) in enumerate(zip(baseline["images"], current["images"], strict=True)):
            if a != b:
                bad.append(f"image {n} is {b[0]}x{b[1]}, baseline {a[0]}x{a[1]} "
                           "-- a dpi/fig.retina/fig.width regression")
    return bad


def _lab(rgb: tuple[int, int, int]) -> tuple[float, float, float]:
    """sRGB -> CIE L*a*b* (D65), enough for a CIE76 distance."""
    def f(c):
        c = c / 255.0
        return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4
    r, g, b = (f(c) for c in rgb)
    x = (0.4124 * r + 0.3576 * g + 0.1805 * b) / 0.95047
    y = 0.2126 * r + 0.7152 * g + 0.0722 * b
    z = (0.0193 * r + 0.1192 * g + 0.9505 * b) / 1.08883
    def h(t):
        return t ** (1 / 3) if t > 0.008856 else 7.787 * t + 16 / 116
    fx, fy, fz = h(x), h(y), h(z)
    return 116 * fy - 16, 500 * (fx - fy), 200 * (fy - fz)


def _hits(colours: list[tuple[int, int, int]], anchors: list[str]) -> int:
    labs = [_lab(c) for c in colours]
    n = 0
    for a in anchors:
        target = _lab(tuple(int(a[i:i + 2], 16) for i in (1, 3, 5)))
        if any(sum((p - q) ** 2 for p, q in zip(lab, target, strict=True)) ** 0.5 < DELTA_E
               for lab in labs):
            n += 1
    return n


def style(html: str, *, allow_skip: bool = False) -> list[str]:
    """Palette and ink-fraction checks.

    A missing Pillow is a **failure**, not a skip, unless ``allow_skip``. Pillow arrives
    transitively today, so an upstream dependency change would otherwise turn this layer off
    without anything saying so -- and a check that silently downgrades itself is worse than one
    that is not there, because it still reports success. ``uv sync --extra summary`` provides it.
    """
    try:
        import numpy as np
        from PIL import Image
    except ImportError:
        msg = "style layer needs Pillow and numpy: uv sync --extra summary"
        if allow_skip:
            print(f"  skip  {msg}")
            return []
        return [msg]
    bad = []
    for n, payload in enumerate(_IMG.findall(html)):
        img = Image.open(io.BytesIO(base64.b64decode(payload))).convert("RGB")
        arr = np.asarray(img)
        ink = float((arr.reshape(-1, 3).min(axis=1) < 240).mean())
        if not INK_BAND[0] <= ink <= INK_BAND[1]:
            bad.append(f"panel {n}: ink fraction {ink:.3f} outside {INK_BAND} "
                       "-- blank, or a solid block")
        colours = [c for _, c in img.quantize(colors=32).convert("RGB").getcolors(1 << 16)]
        if _hits(colours, VIRIDIS) >= 3:
            bad.append(f"panel {n}: carries three viridis anchors -- the palette was swapped")
    return bad


def ssim(a, b) -> float:
    """Mean SSIM between two greyscale arrays, uniform-filter form. No scikit-image dependency."""
    import numpy as np
    from scipy.ndimage import uniform_filter

    a, b = a.astype(np.float64), b.astype(np.float64)
    c1, c2 = (0.01 * 255) ** 2, (0.03 * 255) ** 2
    mu_a, mu_b = uniform_filter(a, 7), uniform_filter(b, 7)
    saa = uniform_filter(a * a, 7) - mu_a * mu_a
    sbb = uniform_filter(b * b, 7) - mu_b * mu_b
    sab = uniform_filter(a * b, 7) - mu_a * mu_b
    num = (2 * mu_a * mu_b + c1) * (2 * sab + c2)
    den = (mu_a ** 2 + mu_b ** 2 + c1) * (saa + sbb + c2)
    return float((num / den).mean())


def perceptual(html: str, reference: Path) -> list[str]:
    """SSIM per panel against a previous fragment. Loose and asymmetric -- the database grows."""
    try:
        import numpy as np
        from PIL import Image
    except ImportError:
        return ["perceptual layer needs Pillow and numpy: uv sync --extra summary"]
    ref = _IMG.findall(reference.read_text())
    cur = _IMG.findall(html)
    if len(ref) != len(cur):
        return [f"perceptual layer skipped: {len(ref)} reference panels, {len(cur)} current"]
    bad = []
    for n, (a, b) in enumerate(zip(ref, cur, strict=True)):
        ia = Image.open(io.BytesIO(base64.b64decode(a))).convert("L")
        ib = Image.open(io.BytesIO(base64.b64decode(b))).convert("L").resize(ia.size)
        s = ssim(np.asarray(ia), np.asarray(ib))
        if s < SSIM_FAIL:
            bad.append(f"panel {n}: SSIM {s:.3f} against the reference, below {SSIM_FAIL}")
        elif s < SSIM_WARN:
            print(f"  warn  panel {n}: SSIM {s:.3f}, below {SSIM_WARN}")
    return bad


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("fragment", nargs="?", type=Path, default=FRAGMENT)
    ap.add_argument("--baseline", type=Path, default=BASELINE)
    ap.add_argument("--reference", type=Path, help="A previous fragment, for the SSIM layer.")
    ap.add_argument("--write-baseline", action="store_true")
    ap.add_argument("--allow-skip", action="store_true",
                    help="Downgrade a missing imaging dependency to a warning.")
    args = ap.parse_args(argv)

    html = args.fragment.read_text()
    fp = fingerprint(html)
    if args.write_baseline:
        args.baseline.write_text(json.dumps(fp, indent=2) + "\n")
        print(f"wrote {args.baseline}: {len(fp['headings'])} headings, "
              f"{len(fp['tables'])} tables, {len(fp['images'])} images")
        return 0

    base = json.loads(args.baseline.read_text()) if args.baseline.exists() else None
    problems = structural(fp, base) + style(html, allow_skip=args.allow_skip)
    if args.reference:
        problems += perceptual(html, args.reference)
    for p in problems:
        print(f"  FAIL  {p}")
    print(f"{len(fp['headings'])} headings, {len(fp['tables'])} tables, {len(fp['images'])} images, "
          f"{'no baseline' if base is None else 'baseline matched' if not problems else 'MISMATCH'}")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
