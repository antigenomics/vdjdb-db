"""The legacy V/J guesser, isolated so phase 8 can delete exactly this file.

``arda.cdr3fix`` repairs a CDR3 against a named segment; it does not propose one. VDJdb chunks may
leave ``v.alpha`` or ``j.beta`` blank -- the README allows it explicitly -- and the legacy build
then guessed from the CDR3 by longest k-mer hit against the germline parts.

Keeping that guess while swapping the repair makes phase 5 a single, attributable deviation.
Issue #462 replaces the scan with OLGA Pgen scoring of the candidate set.
"""
from __future__ import annotations

import polars as pl

from .config import Paths


def guess_segments(keys: pl.DataFrame, gene: str | None = None) -> pl.DataFrame:
    """Add ``__gv`` / ``__gj``: the given segments, with blanks filled by the k-mer guess.

    ``v`` and ``j`` are left untouched -- they are the join key back to the table.
    """
    from .annotate._legacy_fixer import Cdr3Fixer

    res = Paths.discover().res
    fx = Cdr3Fixer(str(res / "segments.txt"), str(res / "segments.aaparts.txt"))

    vs, js = [], []
    for species, cdr3, v, j in keys.select("species", "cdr3", "v", "j").iter_rows():
        # The legacy guesser is keyed on the chain word, not the locus; when the caller does not
        # know it, read it off the segment that *is* present.
        g = gene or ("alpha" if (v or j).startswith("TRA") else "beta")
        vs.append(v or fx.guess_id(cdr3, species, g, True) or "")
        js.append(j or fx.guess_id(cdr3, species, g, False) or "")
    return keys.with_columns(pl.Series("__gv", vs, dtype=pl.Utf8),
                             pl.Series("__gj", js, dtype=pl.Utf8))
