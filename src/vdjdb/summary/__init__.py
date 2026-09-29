"""The release dashboard: its offline inputs, and the driver that renders it.

The dashboard is an R/rmarkdown document (``summary/vdjdb_summary.Rmd``) whose published fragment
``vdjdb-web`` injects into ``/overview``. Everything in this package exists to make that render
offline and deterministic: at the time of writing it made a live NCBI call at render time and
fell back to a fourteen-entry table of non-PubMed references last updated in 2021, so a reference
added since then dropped out of the by-year plots with no error.

:mod:`~vdjdb.summary.references` replaces that call with a committed, reviewed table. It is refreshed
by its own pull request and never written by a build -- the explicit carve-out in hard rule 9 for
a derived table that exists to make the build offline.
"""
from __future__ import annotations

from . import interactive, references, render

__all__ = ["interactive", "references", "render"]
