"""The legacy CDR3 fixer, vendored **unchanged**, as a bridge.

``py_src/Cdr3Fixer.py`` and its two helpers are copied here verbatim so the assembly can move off
``py_src/`` without changing a single output value. They are not ported and not improved: ROADMAP
phase 5 replaces the whole thing with ``arda.cdr3fix``, and porting 400 lines of k-mer scanning to
polars only to delete it next phase would be waste.

The one thing that *is* different is how they are called -- once per distinct
``(species, cdr3, v, j)`` rather than once per row (CLAUDE.md rule 4). That is a pure speed-up: the
functions are deterministic in their arguments, so deduplicating and joining back cannot change a
value.

**Delete this package in phase 5**, together with ``res/segments.txt`` and
``res/segments.aaparts.txt``.
"""
from __future__ import annotations

import sys
from pathlib import Path

# The vendored modules import each other by bare module name, as they did in py_src/.
sys.path.insert(0, str(Path(__file__).parent))

from Cdr3Fixer import Cdr3Fixer  # noqa: E402

__all__ = ["Cdr3Fixer"]
