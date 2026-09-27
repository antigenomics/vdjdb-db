"""Most-plausible junction nucleotide sequences (#461).

VDJdb records amino-acid junctions; a great deal of downstream work -- generation probability,
recombination markup, full-contig synthesis -- needs nucleotides. ``vdjtools.model.infer_nt`` asks
the recombination model for the single most likely nucleotide junction behind an amino-acid one, so
what this produces is inferred, not observed, and the table says so: ``cdr3nt.pgen`` is its
generation probability and ``cdr3nt.margin`` how far it beat the runner-up.

Measured: on 600 distinct human TRB keys the OLGA and arda models agree on 7.2 % of the nucleotide
sequences they both return (293 both-resolved). The models disagree about which of many synonymous
nucleotide histories is most likely, not about the protein. Treat ``cdr3nt`` as a plausible
representative, never as evidence.

Nothing is cached (CLAUDE.md rule 9). The call is deterministic in ``(species, gene, cdr3, v, j)``,
so it runs once per distinct key within the build and joins back -- that is rule 4's deduplication,
not a stored result -- and all 230 chunks cost minutes.
"""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor

import polars as pl

#: VDJdb species -> (``vdjtools`` model source, organism). OLGA is human-only and is the community
#: reference for human; arda is the only source with mouse. Species absent here get no ``cdr3nt``:
#: there is no model, and a guess from the wrong organism's marginals would be worse than an empty
#: cell. The corpus also has 1,402 MacacaMulatta chains, which is why this is a lookup and not an
#: assertion.
MODELS: dict[str, tuple[str, str]] = {
    "HomoSapiens": ("olga", "human"),
    "MusMusculus": ("arda", "mouse"),
}

#: The columns this stage adds to ``chains``. The D geometry comes from the same scenario that
#: produced ``cdr3nt``, which is why it belongs in this module: ``d.start`` and ``d.end`` index that
#: nucleotide sequence, so taking them from a second model would ship coordinates that do not point
#: at the sequence beside them. ``arda.dpost`` supplies how much to believe the call
#: (:mod:`vdjdb.annotate.dgene`), not where it sits.
NT_COLUMNS: tuple[str, ...] = ("cdr3nt", "cdr3nt.pgen", "cdr3nt.margin",
                               "d.inferred", "d.start", "d.end")

#: Contiguous slices of the sorted key set, one per worker. Never a pool of per-record tasks: the
#: native call releases the GIL only partly (measured 2.11x on 4 threads), and dispatch on 114k
#: one-row tasks would cost more than it saves. Reassembled in slice order, so the worker count
#: cannot change the answer (CLAUDE.md rule 7).
SLICES = 4

#: ``infer_nt`` takes what it calls ``cdr3_aa``, but it means the junction -- Cys104..Phe/Trp118
#: inclusive, which is what VDJdb's ``cdr3`` column holds. Passing a true CDR3 would infer a
#: junction two codons short, with no error. See :mod:`vdjdb.convert.coords`.
_JUNCTION_IS_WHAT_IT_WANTS = True


def _resolver(model: object, table: str, column: str) -> dict[str, str | None]:
    """A call -> the allele to pass the model, or ``None`` to marginalise over that segment.

    The model is keyed by allele and raises on anything else, so every call is resolved before the
    run rather than by catching an exception per sequence. Three outcomes, none of which invents a
    call:

    * a known allele -- used as given;
    * a gene name whose model has exactly one allele -- that allele, because there is no choice to
      make. VDJdb has 2,371 V and 1,784 J calls with no allele at all (#389);
    * anything else -- absent from the mapping, so ``.get()`` yields ``None`` and the model
      marginalises over that segment. (``.get`` is the accessor; a missing key is the third case,
      not an error.) That is 3,500 of
      284,764 chains (1.2 %), and the reasons read as a catalogue for phase 9: ``TRBJ1.2`` and
      ``TRBJ 2-7`` (a dot and a space where a dash belongs), ``TRAJ01-1*01`` (zero-padded),
      ``TRAV21-DV12`` against the model's ``TRAV21/DV12*01``, alleles no model lists
      (``TRBV19*02``, 107 chains) and pseudogenes (``TRBV21-1``, 303).
    """
    genes = model.genomic[table]                                   # type: ignore[attr-defined]
    alleles = {a: a for a in genes[column]}
    sole = dict(genes.join(genes.group_by("gene").len().filter(pl.col("len") == 1),
                           on="gene", how="inner").select("gene", column).iter_rows())
    return {"": None, **{g: a for g, a in sole.items() if g not in alleles}, **alleles}


def _infer_slice(rows: list[tuple[str | None, str | None, str]], model: object) -> list[tuple]:
    from vdjtools.model import infer_nt

    out = []
    for v, j, cdr3 in rows:
        s = infer_nt(model, cdr3, v=v, j=j)
        out.append(("", None, None, "", None, None) if s is None else
                   (s.cdr3_nt, s.pgen,
                    s.pgen / s.runner_up_pgen if s.runner_up_pgen else float("inf"),
                    s.d_call or "", s.d_start, s.d_end))
    return out


def infer(keys: pl.DataFrame, species: str, gene: str, *, workers: int = SLICES) -> pl.DataFrame:
    """``(cdr3, v.segm, j.segm)`` -> the same frame plus :data:`NT_COLUMNS`.

    ``keys`` must already be distinct and sorted; an unsupported species returns empty columns
    rather than raising, because a corpus is allowed to contain species no model covers.
    """
    from vdjtools.model import load_bundled

    blank = keys.with_columns(pl.lit("").alias("cdr3nt"),
                              pl.lit(None, pl.Float64).alias("cdr3nt.pgen"),
                              pl.lit(None, pl.Float64).alias("cdr3nt.margin"),
                              pl.lit("").alias("d.inferred"),
                              pl.lit(None, pl.Int64).alias("d.start"),
                              pl.lit(None, pl.Int64).alias("d.end"))
    if species not in MODELS or keys.is_empty():
        return blank
    source, organism = MODELS[species]
    model = load_bundled(gene, source, organism=organism)
    vmap = _resolver(model, "genes_v", "v_allele")
    jmap = _resolver(model, "genes_j", "j_allele")

    rows = [(vmap.get(v), jmap.get(j), cdr3)
            for cdr3, v, j in keys.select("cdr3", "v.segm", "j.segm").iter_rows()]
    n = max(1, min(workers, len(rows)))
    bounds = [(i * len(rows) // n, (i + 1) * len(rows) // n) for i in range(n)]
    with ThreadPoolExecutor(n) as pool:
        parts = list(pool.map(lambda b: _infer_slice(rows[b[0]:b[1]], model), bounds))
    flat = [r for part in parts for r in part]          # slice order, never completion order

    return keys.with_columns(
        pl.Series("cdr3nt", [r[0] for r in flat], dtype=pl.Utf8),
        pl.Series("cdr3nt.pgen", [r[1] for r in flat], dtype=pl.Float64),
        pl.Series("cdr3nt.margin", [r[2] for r in flat], dtype=pl.Float64),
        # 0-based half-open, in the coordinate space of `cdr3nt` above -- vdjtools' Scenario space
        # (CLAUDE.md). TRA has no D, so these stay empty there by construction.
        pl.Series("d.inferred", [r[3] for r in flat], dtype=pl.Utf8),
        pl.Series("d.start", [r[4] for r in flat], dtype=pl.Int64),
        pl.Series("d.end", [r[5] for r in flat], dtype=pl.Int64),
    )


def add_junction_nt(chains: pl.DataFrame, records: pl.DataFrame,
                    *, workers: int = SLICES) -> pl.DataFrame:
    """Add :data:`NT_COLUMNS` to ``chains``, one model load per (species, locus)."""
    species = records.select("record_id", "species")
    keyed = chains.join(species, on="record_id", how="left")
    resolvable = ((pl.col("cdr3") != "") & (pl.col("v.segm") != "") & (pl.col("j.segm") != ""))

    parts = []
    for (sp, gene), group in sorted(
            keyed.filter(resolvable).group_by("species", "gene", maintain_order=True),
            key=lambda kv: kv[0]):          # sorted: the join order must not vary by group order
        keys = group.select("cdr3", "v.segm", "j.segm").unique().sort("cdr3", "v.segm", "j.segm")
        parts.append(infer(keys, sp, gene, workers=workers)
                     .with_columns(pl.lit(sp).alias("species"), pl.lit(gene).alias("gene")))

    lookup = (pl.concat(parts, how="vertical") if parts else
              keyed.head(0).select("cdr3", "v.segm", "j.segm", "species", "gene",
                                   pl.lit("").alias("cdr3nt"),
                                   pl.lit(None, pl.Float64).alias("cdr3nt.pgen"),
                                   pl.lit(None, pl.Float64).alias("cdr3nt.margin"),
                                   pl.lit("").alias("d.inferred"),
                                   pl.lit(None, pl.Int64).alias("d.start"),
                                   pl.lit(None, pl.Int64).alias("d.end")))
    return (keyed.join(lookup, on=["species", "gene", "cdr3", "v.segm", "j.segm"], how="left")
            .with_columns(pl.col("cdr3nt", "d.inferred").fill_null(""))   # rule 6
            .drop("species")
            .sort("record_id", "gene"))
