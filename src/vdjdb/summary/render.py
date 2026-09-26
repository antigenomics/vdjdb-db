"""Render the release dashboard and verify the fragment `vdjdb-web` will serve.

Three steps, each able to fail the build on its own:

1. ``rmarkdown::render`` on ``summary/vdjdb_summary.Rmd``, pointed at the legacy projection of this
   build, with ``clean = FALSE`` so knitr's intermediate survives. The render makes **no network
   call** -- publication years come from the committed table :mod:`vdjdb.summary.references`
   writes -- so it is reproducible and runs offline.
2. one more `pandoc` pass over that intermediate with ``summary/embed.html`` and
   ``summary/embed.lua``, which **emit the publishable fragment directly**.
3. ``summary/check_summary.py`` asserts the structure, the palette and the three contracts
   ``vdjdb-web`` depends on, against the committed fingerprint.

**Step 2 replaced a post-processing script and that is the point.** ``MakeEmbedableHtml.py``
line-scanned pandoc's finished output for ``<div``, ``<pre class="r">`` and ``<table>`` -- three
guesses about markup that each silently blank ``/overview`` when they stop matching. Stating the
same transforms on the document tree costs one extra pandoc invocation (the R never runs twice) and
removes the guesses. Verified against the string-surgery output: identical line count, identical
image digests, and the only textual difference is the position of one attribute inside 8 ``<img>``
tags.

``-auto_identifiers`` is deliberate: without it pandoc puts an ``id`` on every ``<h4>``, which the
rmarkdown path did not because it hung ids on section ``<div>``s instead. Adding anchors to the
fragment is a reasonable thing to want and a separate decision from this refactor.

The paper figures (``summary/vdjdb_paper_figures.Rmd``) are not part of a release and are not
rendered here; they need ``maps`` and ``scatterpie``, which the release path deliberately does not.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

SUMMARY = Path("summary")
RMD = SUMMARY / "vdjdb_summary.Rmd"
RENDERED = SUMMARY / "vdjdb_summary.html"
INTERMEDIATE = SUMMARY / "vdjdb_summary.knit.md"
FIGURES = SUMMARY / "vdjdb_summary_files" / "figure-html"
FRAGMENT = SUMMARY / "vdjdb_summary_embed.html"
TEMPLATE = "embed.html"
FILTER = "embed.lua"

#: The reader pandoc must use on knitr's intermediate. The first three extensions are rmarkdown's
#: own; `-auto_identifiers` keeps the fragment's markup as it has always been (see the docstring).
READER = "markdown+autolink_bare_uris+tex_math_single_backslash-auto_identifiers"


def render(legacy: Path, *, quiet: bool = True) -> Path:
    """Run ``rmarkdown::render``, keeping knitr's intermediate for :func:`extract`."""
    if shutil.which("Rscript") is None:
        raise RuntimeError("Rscript is not on PATH; the dashboard needs R and rmarkdown.")
    # `legacy` is resolved here rather than inside the Rmd because knitr sets the working
    # directory to the document's own, so a relative path would mean something different there.
    expr = (f'rmarkdown::render("{RMD}", quiet={"TRUE" if quiet else "FALSE"}, clean = FALSE, '
            f'params = list(legacy = "{legacy.resolve()}"))')
    subprocess.run(["Rscript", "-e", expr], check=True)
    return RENDERED


def extract(intermediate: Path = INTERMEDIATE, fragment: Path = FRAGMENT, *,
            assets: Path | None = None, asset_prefix: str = "figures/") -> int:
    """One pandoc pass over knitr's intermediate that emits the fragment. Returns its line count.

    Run from ``summary/`` because the intermediate references its figures relatively and
    ``--embed-resources`` resolves them from the working directory.

    ``assets`` switches the figures from inlined base64 to files written there, referenced as
    ``asset_prefix + name``. Measured on the current dashboard: the fragment goes from **5.14 MB to
    94.8 KB, 55x smaller**, and the 3.78 MB of PNGs become separately cacheable rather than
    re-sent on every page load.

    **It is not the default, and cannot be until vdjdb-web changes.** The Scala side matches
    ``data:image/png;base64`` to find the images; pointed at a fragment with external ``src``
    attributes it would render eight broken images. The capability ships ready for that change.
    """
    if shutil.which("pandoc") is None:
        raise RuntimeError("pandoc is not on PATH; the fragment is produced by it.")
    if not intermediate.exists():
        raise FileNotFoundError(
            f"{intermediate} is missing -- render() must run with `clean = FALSE`, or knitr "
            "deletes the intermediate this pass reads.")
    # Built in one piece rather than inserted into positionally: an index into an argv is a
    # statement about a list that has already changed once.
    cmd = ["pandoc", intermediate.name,
           "--from", READER, "--to", "html4", "--standalone",
           "--syntax-highlighting", "none",
           "--template", TEMPLATE, "--lua-filter", FILTER,
           *(["--embed-resources"] if assets is None
             else ["-M", f"asset-prefix={asset_prefix}"]),
           # Absolute: pandoc runs in `summary/` so the figures resolve relatively, but the
           # fragment may be written anywhere -- a bare name would land beside the intermediate.
           "-o", str(fragment.resolve())]
    subprocess.run(cmd, cwd=intermediate.parent, check=True)
    if assets is not None:
        assets.mkdir(parents=True, exist_ok=True)
        for png in sorted(FIGURES.glob("*.png")):
            (assets / png.name).write_bytes(png.read_bytes())
    return len(fragment.read_text().splitlines())


def check(fragment: Path = FRAGMENT, reference: Path | None = None) -> int:
    cmd = [sys.executable, str(SUMMARY / "check_summary.py"), str(fragment)]
    if reference:
        cmd += ["--reference", str(reference)]
    return subprocess.run(cmd, check=False).returncode
