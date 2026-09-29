"""Most-plausible junction nucleotide sequences (#461).

VDJdb records amino-acid junctions; a great deal of downstream work -- generation probability,
recombination markup, full-contig synthesis -- needs nucleotides. ``vdjtools.model.infer_nt_batch``
asks the recombination model for the single most likely nucleotide junction behind each amino-acid
one, so what this produces is inferred, not observed, and the table says so: ``cdr3nt.pgen`` is its
generation probability and ``cdr3nt.margin`` how far it beat the runner-up.

Measured: on 600 distinct human TRB keys the OLGA and arda models agree on 7.2 % of the nucleotide
sequences they both return (293 both-resolved). The models disagree about which of many synonymous
nucleotide histories is most likely, not about the protein. Treat ``cdr3nt`` as a plausible
representative, never as evidence.

**One batched call per (species, locus), and nothing around it.** ``infer_nt_batch`` releases the GIL
and partitions the batch across its own kernel threads, so a pool of our own would oversubscribe the
machine and read as "batching did not help" -- its docstring says so, and hard rule 3 says it here.
This stage ran as four ``vdjdb infer-nt`` processes over contiguous slices until `vdjtools` 4.5
published the batch entry point (`antigenomics/vdjtools#181`, opened from this build's profile);
measured on 3,000 distinct human TRB keys, **1.115 ms/key serial against 0.106 ms batched, 10.5x, and
all 3,000 nucleotide sequences identical**.

Nothing is cached (CLAUDE.md rule 9). The call is deterministic in ``(species, gene, cdr3, v, j)``,
so it runs once per distinct key within the build and joins back -- that is rule 4's deduplication,
not a stored result.
"""
from __future__ import annotations

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
                               "d.inferred", "d.start", "d.end",
                               "v.end.inferred", "j.start.inferred")

#: ``infer_nt_batch`` takes what it calls ``cdr3_aas``, but it means junctions -- Cys104..Phe/Trp118
#: inclusive, which is what VDJdb's ``cdr3`` column holds. Passing true CDR3s would infer junctions
#: two codons short, with no error. See :mod:`vdjdb.convert.coords`.
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


#: What an unmapped V/J boundary reads as. Already the value ``v.end``/``j.start`` carry where the
#: markup engine declines, so the fallback columns use it too rather than inventing a second marker.
UNMAPPED = -1


def blank_columns(keys: pl.DataFrame) -> pl.DataFrame:
    """``keys`` plus :data:`NT_COLUMNS`, all empty. The answer for a species with no model."""
    return keys.with_columns(pl.lit("").alias("cdr3nt"),
                             pl.lit(None, pl.Float64).alias("cdr3nt.pgen"),
                             pl.lit(None, pl.Float64).alias("cdr3nt.margin"),
                             pl.lit("").alias("d.inferred"),
                             pl.lit(None, pl.Int64).alias("d.start"),
                             pl.lit(None, pl.Int64).alias("d.end"),
                             pl.lit(UNMAPPED, pl.Int64).alias("v.end.inferred"),
                             pl.lit(UNMAPPED, pl.Int64).alias("j.start.inferred"))


def infer(keys: pl.DataFrame, species: str, gene: str) -> pl.DataFrame:
    """``(cdr3, v.segm, j.segm)`` -> the same frame plus :data:`NT_COLUMNS`.

    ``keys`` must already be distinct and sorted; an unsupported species returns empty columns rather
    than raising, because a corpus is allowed to contain species no model covers.

    ``infer_nt_batch`` returns one row per input row in input order, with nulls where the model
    cannot explain a junction, so the result is joined positionally and the worker count cannot
    change the answer -- there is no worker count (CLAUDE.md rule 7).
    """
    from vdjtools.model import germline_boundary, infer_nt_batch, load_bundled

    from ..convert.coords import nt_to_aa_boundary_expr

    if species not in MODELS or keys.is_empty():
        return blank_columns(keys)
    source, organism = MODELS[species]
    model = load_bundled(gene, source, organism=organism)
    vmap = _resolver(model, "genes_v", "v_allele")
    jmap = _resolver(model, "genes_j", "j_allele")
    cdr3, v, j = (keys["cdr3"].to_list(),
                  [vmap.get(x) for x in keys["v.segm"]],
                  [jmap.get(x) for x in keys["j.segm"]])
    got = infer_nt_batch(model, cdr3, v=v, j=j)
    # A second batched call, for the boundary only. `Scenario.v_end` is a property of the argmax
    # recombination history, and maximising P(sequence) explains N-region nucleotides as templated
    # whenever it can, so it credits the germline too far (`antigenomics/vdjtools#182`, opened from
    # this column's measurement). `germline_boundary` asks the germline alignment instead. Measured
    # against `isalgo/airr_control`'s `human.trb.ntvj` on 8,132 VDJdb human TRB junctions with
    # unambiguous nucleotide truth: `v.end` exact 7,554 of 8,132 (92.89 %) against the history's
    # 7,244 (89.08 %), `j.start` 7,966 (97.96 %) against 7,483 (92.02 %). It costs 2.67 us/row
    # against `infer_nt_batch`'s 82.81, so 3.2 % for four points of accuracy on one coordinate and
    # six on the other.
    #
    # It needs a germline to align against, so where our call was `None` and the DP marginalised, the
    # allele the DP settled on goes in - the same allele `v.inferred` records. Passing `None` through
    # instead would decline on 2,221 chains the scenario answered, and a boundary conditional on a
    # named allele is what both columns then mean.
    bounds = germline_boundary(
        model, cdr3,
        v=[a or b for a, b in zip(v, got["v_call"].to_list(), strict=True)],
        j=[a or b for a, b in zip(j, got["j_call"].to_list(), strict=True)])
    return keys.with_columns(
        got["cdr3_nt"].fill_null("").alias("cdr3nt"),                       # rule 6
        got["pgen"].alias("cdr3nt.pgen"),
        # How far the winner beat the runner-up. No runner-up means nothing else was in contention,
        # which is an unbounded margin rather than a missing one.
        pl.when(got["runner_up_pgen"].is_null() | (got["runner_up_pgen"] == 0))
          .then(pl.when(got["pgen"].is_null()).then(None).otherwise(float("inf")))
          .otherwise(got["pgen"] / got["runner_up_pgen"]).alias("cdr3nt.margin"),
        # 0-based half-open, in the coordinate space of `cdr3nt` above -- vdjtools' Scenario space
        # (CLAUDE.md). TRA has no D, so these stay empty there by construction.
        got["d_call"].fill_null("").alias("d.inferred"),
        got["d_start"].cast(pl.Int64).alias("d.start"),
        got["d_end"].cast(pl.Int64).alias("d.end"),
        # The germline alignment's V/J boundary, converted out of vdjtools' nucleotide space into
        # the `v.end`/`j.start` residue space by the fitted ceiling (`convert.coords`). Null where
        # there is no germline to align against - no call, or a call this model does not carry -
        # which `nt_to_aa_boundary_expr` reads as UNMAPPED. `add_junction_nt` masks it again against
        # the markup engine's answer, so it only ever fills where that declined.
        nt_to_aa_boundary_expr(bounds["v_end"], unmapped=UNMAPPED).alias("v.end.inferred"),
        nt_to_aa_boundary_expr(bounds["j_start"], unmapped=UNMAPPED).alias("j.start.inferred"))


def add_junction_nt(chains: pl.DataFrame, records: pl.DataFrame) -> pl.DataFrame:
    """Add :data:`NT_COLUMNS` to ``chains``, one model load and one batched call per (species, locus)."""
    species = records.select("record_id", "species")
    keyed = chains.join(species, on="record_id", how="left")
    resolvable = ((pl.col("cdr3") != "") & (pl.col("v.segm") != "") & (pl.col("j.segm") != ""))

    parts = []
    for (sp, gene), group in sorted(
            keyed.filter(resolvable).group_by("species", "gene", maintain_order=True),
            key=lambda kv: kv[0]):          # sorted: the join order must not vary by group order
        keys = group.select("cdr3", "v.segm", "j.segm").unique().sort("cdr3", "v.segm", "j.segm")
        parts.append(infer(keys, sp, gene)
                     .with_columns(pl.lit(sp).alias("species"), pl.lit(gene).alias("gene")))

    lookup = (pl.concat(parts, how="vertical") if parts else
              keyed.head(0).select("cdr3", "v.segm", "j.segm", "species", "gene",
                                   pl.lit("").alias("cdr3nt"),
                                   pl.lit(None, pl.Float64).alias("cdr3nt.pgen"),
                                   pl.lit(None, pl.Float64).alias("cdr3nt.margin"),
                                   pl.lit("").alias("d.inferred"),
                                   pl.lit(None, pl.Int64).alias("d.start"),
                                   pl.lit(None, pl.Int64).alias("d.end"),
                                   pl.lit(UNMAPPED, pl.Int64).alias("v.end.inferred"),
                                   pl.lit(UNMAPPED, pl.Int64).alias("j.start.inferred")))
    return (keyed.join(lookup, on=["species", "gene", "cdr3", "v.segm", "j.segm"], how="left")
            .with_columns(pl.col("cdr3nt", "d.inferred").fill_null(""))   # rule 6
            # A fallback, never an override (#631). The model's boundary survives only where the
            # markup engine declined; where the alignment answered, the column reads UNMAPPED and a
            # consumer coalescing the two cannot overwrite an alignment answer even by accident.
            # That is stronger than filling the shipped column, and it is what keeps the legacy
            # export - which reads `v.end` and `j.start` - byte-identical.
            .with_columns(*[
                pl.when(pl.col(shipped) == UNMAPPED)
                  .then(pl.col(fallback).fill_null(UNMAPPED))
                  .otherwise(pl.lit(UNMAPPED, pl.Int64)).alias(fallback)
                for shipped, fallback in (("v.end", "v.end.inferred"),
                                          ("j.start", "j.start.inferred"))])
            .drop("species")
            .sort("record_id", "gene"))
