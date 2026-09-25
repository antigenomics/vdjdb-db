"""How much to believe the D call, from ``arda.dpost``.

The D segment is short, heavily trimmed, and there are only two of them in TRB, so assigning one from
a junction is often close to a coin flip. A point estimate that does not say so is misleading, which
is what this module fixes: :mod:`vdjdb.annotate.junction` says *which* D and *where*, and this says
*how sure*.

Measured on 2,000 distinct human TRB keys: the posterior for the winning gene has median **0.791**
and falls below 0.6 on **21.8 %**; the entropy over the posterior has median **0.740** and exceeds
0.9 -- essentially undecidable between TRBD1 and TRBD2 -- on **28.9 %**. Independently, ``infer_nt``
and ``arda.dpost`` name the same D gene on only **78.5 %**, and the inferred call matches the curated
``d.segm`` at gene level on 76.5 %. All three numbers say the same thing, so ``d.posterior`` is not a
decoration: a consumer that filters on it is doing the only correct thing with a D call.

The posterior reported is for the gene ``d.inferred`` names, **not** for arda's own winner. That is
deliberate: the two disagree on a fifth of chains, and the number beside a call must be the
probability of *that* call. When they disagree the posterior is low, which is exactly the signal.
"""
from __future__ import annotations

import polars as pl

#: VDJdb species -> arda's name. Only loci with a D and a shipped model appear; TRA never does.
SPECIES: dict[str, str] = {"HomoSapiens": "human", "MusMusculus": "mouse"}

#: The columns this stage adds to ``chains``.
D_COLUMNS: tuple[str, ...] = ("d.posterior", "d.entropy")

_KEY = ("species", "cdr3", "v.segm", "j.segm", "d.inferred")


def add_d_posterior(chains: pl.DataFrame, records: pl.DataFrame) -> pl.DataFrame:
    """Add :data:`D_COLUMNS` to ``chains``, null wherever no D was inferred."""
    from arda.dpost import posterior_d

    from .cdr3fix import ensure_reference

    keyed = chains.join(records.select("record_id", "species"), on="record_id", how="left")
    target = keyed.filter((pl.col("d.inferred") != "")
                          & pl.col("species").is_in(list(SPECIES)))
    blank = (pl.lit(None, pl.Float64).alias("d.posterior"),
             pl.lit(None, pl.Float64).alias("d.entropy"))
    if target.is_empty():
        return keyed.with_columns(*blank).drop("species").sort("record_id", "gene")

    # Without this, arda mistakes this repository for its own checkout and every call returns None
    # -- silently, which is how it went unnoticed for a whole phase. See annotate/cdr3fix.py.
    ensure_reference()

    keys = target.select(_KEY).unique().sort(_KEY)
    post, ent = [], []
    for sp, cdr3, v, j, d in keys.iter_rows():
        p = posterior_d(cdr3, v, j, SPECIES[sp])
        # by_gene is keyed on the gene; `d.inferred` is an allele. A gene arda did not consider has
        # probability 0 under its model, which is a real answer, not a missing one.
        post.append(None if p is None else p.by_gene.get(d.split("*")[0], 0.0))
        ent.append(None if p is None else p.entropy)

    lookup = keys.with_columns(pl.Series("d.posterior", post, dtype=pl.Float64),
                               pl.Series("d.entropy", ent, dtype=pl.Float64))
    return (keyed.join(lookup, on=_KEY, how="left")
            .drop("species")
            .sort("record_id", "gene"))
