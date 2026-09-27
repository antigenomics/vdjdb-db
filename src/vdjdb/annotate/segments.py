"""A V or J call for the chains that have none (#462).

The README lets a chunk leave ``v.alpha`` or ``j.beta`` blank, and 3,712 of 286,047 chains do
(1.3 %): 711 with no V, 3,521 with no J. The legacy build guessed from the CDR3 by longest k-mer hit
against the germline parts. Measured against the curated calls on 1,200 human chains per locus, with
the true call hidden:

======== ====================== ======================
locus    k-mer scan             Pgen (``infer_nt``)
======== ====================== ======================
TRB V    0 of 1,189 (0.0 %)      283 (23.8 %)
TRA V    1 of 1,111 (0.1 %)      557 (50.1 %)
TRB J    1,133 (95.9 %)          1,153 (97.5 %)
TRA J    882 (79.3 %)            1,065 (95.8 %)
======== ====================== ======================

The V column reflects a bug rather than a comparison: ``Cdr3Fixer.guess_id`` puts ``return ""``
inside the five-prime loop, so it tries one prefix length and gives up -- 3 non-empty guesses in
4,000 sequences. The J branch has the same statement correctly in a ``for...else``. So VDJdb's V
guesser has never worked, and the 711 chains with no V fail the legacy build's "a CDR3 needs a V and
a J" test, so their records are dropped from `vdjdb.txt`.

What ships here is a new-format column, not a change to that behaviour: ``v.inferred`` and
``j.inferred`` are filled only where the curated call is missing, and the curated columns are
never touched. Whether the legacy build should start keeping those records is a decision about
shipped data, and it belongs with #327 in phase 9.

Read the accuracy above before using ``v.inferred``: a V call recovered from the junction alone is
right about a quarter of the time for TRB, because TRBV contributes only a few residues to it. The J
side is reliable.
"""
from __future__ import annotations

import polars as pl

from .junction import MODELS

#: The columns this stage adds to ``chains``.
INFERRED_COLUMNS: tuple[str, ...] = ("v.inferred", "j.inferred")

_KEY = ("species", "gene", "cdr3", "v.segm", "j.segm")


def add_inferred_segments(chains: pl.DataFrame, records: pl.DataFrame) -> pl.DataFrame:
    """Add :data:`INFERRED_COLUMNS`, empty wherever the curator already named the segment."""
    from vdjtools.model import infer_nt, load_bundled

    from .junction import _resolver

    keyed = chains.join(records.select("record_id", "species"), on="record_id", how="left")
    target = keyed.filter(
        (pl.col("cdr3") != "")
        & ((pl.col("v.segm") == "") | (pl.col("j.segm") == ""))
        & pl.col("species").is_in(list(MODELS))
    )
    blank = (pl.lit("").alias("v.inferred"), pl.lit("").alias("j.inferred"))
    if target.is_empty():
        return keyed.with_columns(*blank).drop("species").sort("record_id", "gene")

    keys = target.select(_KEY).unique().sort(_KEY)
    groups, vs, js = [], [], []
    # Grouped so each model loads once; sorted so the order never depends on group iteration.
    for (sp, gene), group in sorted(keys.group_by("species", "gene", maintain_order=True),
                                    key=lambda kv: kv[0]):
        source, organism = MODELS[sp]
        model = load_bundled(gene, source, organism=organism)
        vmap = _resolver(model, "genes_v", "v_allele")
        jmap = _resolver(model, "genes_j", "j_allele")
        groups.append(group)
        for _, _, cdr3, v, j in group.iter_rows():
            # The missing side goes in as None, so the model chooses it rather than echoing
            # ours back.
            s = infer_nt(model, cdr3, v=vmap.get(v) if v else None, j=jmap.get(j) if j else None)
            vs.append("" if v or s is None else (s.v_call or ""))
            js.append("" if j or s is None else (s.j_call or ""))

    ordered = pl.concat(groups, how="vertical")
    lookup = ordered.with_columns(pl.Series("v.inferred", vs, dtype=pl.Utf8),
                                  pl.Series("j.inferred", js, dtype=pl.Utf8))
    return (keyed.join(lookup, on=_KEY, how="left")
            .with_columns(pl.col("v.inferred", "j.inferred").fill_null(""))
            .drop("species")
            .sort("record_id", "gene"))
