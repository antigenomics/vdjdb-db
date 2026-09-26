"""Every file the build reads must actually be in the repository.

A source file that `.gitignore` swallows is invisible locally -- it is right there on disk, the code
finds it, the tests pass -- and missing everywhere else. `summary/embed.html` was ignored by
`summary/*.html`, a rule that exists for the *rendered* dashboard, and CI failed four steps deep
with "Could not find data file". `git add -A` says nothing when it skips an ignored file.

The repo has this trap in at least two places: `summary/*.html` and `summary/*.txt`, which is why
committed tables there are `.tsv`. So the check is on the *inputs*, not on the rule.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import pytest


def _tracked() -> set[str]:
    out = subprocess.run(["git", "ls-files"], capture_output=True, text=True, check=True)
    return set(out.stdout.split())


#: Files the build or the docs read at run time and cannot work without. Paths, not globs: a glob
#: would match whatever happens to be on disk, which is the failure mode being guarded against.
REQUIRED_INPUTS = [
    "summary/embed.tpl",
    "summary/embed.lua",
    "summary/install-packages.R",
    "summary/check_summary.py",
    "summary/panels.py",
    "summary/annotations.tsv",
    "summary/reference_years.tsv",
    "summary/fingerprint.json",
    "summary/vdjdb_summary.Rmd",
    "summary/preview/index.html",
    "docs/conf.py",
    "docs/_ext/vdjdb_schema.py",
    "docs/_static/dashboard.html",
    "docs/tuning/scorecard.tsv",
    "attic/BuildDatabase.groovy",
]


@pytest.mark.parametrize("path", REQUIRED_INPUTS)
def test_required_input_is_tracked(path):
    assert Path(path).exists(), f"{path} is missing from the working tree"
    assert path in _tracked(), (
        f"{path} exists on disk but is NOT tracked -- almost certainly swallowed by a .gitignore "
        "rule written for a generated file of the same shape. Rename it out of the pattern rather "
        "than negating the rule; that is what `summary/*.tsv` does.")


def test_the_render_modules_paths_are_all_tracked():
    """Whatever `render.py` points pandoc at, specifically -- it is the one that broke."""
    from vdjdb.summary import render

    tracked = _tracked()
    for p in (render.TEMPLATE, render.FILTER, render.RMD):
        rel = p.relative_to(Path.cwd()) if p.is_absolute() else p
        assert str(rel) in tracked, f"{rel} is not tracked"
