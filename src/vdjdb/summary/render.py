"""Render the release dashboard and verify the fragment `vdjdb-web` will serve.

Three steps, each of which can fail the build on its own:

1. ``rmarkdown::render`` on ``summary/vdjdb_summary.Rmd``, pointed at the legacy projection of this
   build. The render makes **no network call** -- publication years come from the committed table
   :mod:`vdjdb.summary.references` writes -- so it is reproducible and runs offline.
2. ``summary/MakeEmbedableHtml.py`` extracts the publishable fragment.
3. ``summary/check_summary.py`` asserts the structure, the palette and the three contracts
   ``vdjdb-web`` depends on, against the committed fingerprint.

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
FRAGMENT = SUMMARY / "vdjdb_summary_embed.html"


def render(legacy: Path, *, quiet: bool = True) -> Path:
    """Run ``rmarkdown::render`` against the legacy tables of this build."""
    if shutil.which("Rscript") is None:
        raise RuntimeError("Rscript is not on PATH; the dashboard needs R and rmarkdown.")
    # `legacy` is resolved here rather than inside the Rmd because knitr sets the working
    # directory to the document's own, so a relative path would mean something different there.
    expr = (f'rmarkdown::render("{RMD}", quiet={"TRUE" if quiet else "FALSE"}, '
            f'params = list(legacy = "{legacy.resolve()}"))')
    subprocess.run(["Rscript", "-e", expr], check=True)
    return RENDERED


def extract(rendered: Path = RENDERED, fragment: Path = FRAGMENT) -> int:
    proc = subprocess.run([sys.executable, str(SUMMARY / "MakeEmbedableHtml.py"),
                           str(rendered), str(fragment)],
                          capture_output=True, text=True, check=True)
    return int(proc.stdout.split()[0])


def check(fragment: Path = FRAGMENT, reference: Path | None = None) -> int:
    cmd = [sys.executable, str(SUMMARY / "check_summary.py"), str(fragment)]
    if reference:
        cmd += ["--reference", str(reference)]
    return subprocess.run(cmd, check=False).returncode
