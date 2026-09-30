"""Junction nucleotides, D geometry and how sure the D call is - one library call (#461).

VDJdb records amino-acid junctions; a great deal of downstream work - generation probability,
recombination markup, full-contig synthesis - needs nucleotides, and the D segment needs naming and
placing before anything can draw a rearrangement. ``vdjtools.model.annotate_junctions`` answers all
of it in **one batched call**, so what this produces is inferred, not observed, and the table says
so: ``cdr3nt.pgen`` is the junction's generation probability, ``cdr3nt.margin`` how far it beat the
runner-up, and ``d.posterior`` the probability of the gene ``d.inferred`` names.

Treat ``cdr3nt`` as a plausible representative, never as evidence. Two bundled human models, one
bootstrapped from OLGA and one an EM fit over arda's IMGT allele namespace, agree on 7.2 % of the
nucleotide sequences they both return (600 distinct human TRB keys, 293 both-resolved): they
disagree about which of many synonymous nucleotide histories is most likely, not about the protein.

**Naming the D and placing it are separate questions, and one estimator answers each.** The
recombination model names the gene from its own scenario weights; the aligner places it greedily and
ungated. Measured on 4,000 real human TRB rearrangements whose D and coordinates come from the
nucleotide sequence:

===================================================  =====================  ================
model names, aligner places                          74.35 % gene correct   99.80 % placed
gated alignment, model posterior where it declines    71.40 %                55.75 %
the retired length-and-prior posterior                69.67 %                -
E-value-gated alignment chooses and places            47.93 %                55.75 %
===================================================  =====================  ================

That retired posterior was ``arda.dpost``, which this module used to read through a per-key Python
loop for 15.96 s of a 38.73 s build - 35.8 %, the largest stage. It is **gone from arda since
2.33.0 and does not reappear in vdjtools**: the number now comes from the model's own scenario
weights normalised over D genes, which is both better and free, because the pipeline already makes
that call. `annotate/dgene.py` was deleted rather than re-pointed.

⚠ The one thing that got worse is per-row positional precision: ``d.start`` is exact on 64.71 % of
correctly-called rows against 66.67 % under the E-value gate. It is exact on **1,922 rows rather
than 1,278**, because it answers 3,992 rather than 2,230. For a view drawing V/N/D/N/J that is the
trade to take.

⚠ **Every accuracy percentage above is human TRB.** The truth set (``isalgo/airr_control``) carries
TRA and TRB and no immunoglobulin, so none of them may be quoted for IGH.

**The model set is chosen per row, not configured.** ``model_source="auto"`` runs OLGA's bundled fit
first and arda's on whatever it left unexplained, and its last rung is a germline scaffold rather
than a fitted model - which is what makes the non-human records answer at all. A fitted model exists
for human and mouse and for nothing else, so before this every rhesus record came back empty: 1,457
keys answered zero times, invisible inside a single corpus-wide total. **Break any coverage check
down by species** (:mod:`vdjdb.annotate` has no configuration knob to get this wrong, but a test
can).

**One batched call, and nothing around it.** ``annotate_junctions`` takes ``species`` per row, loads
each model once, releases the GIL and partitions across its own kernel threads; every stage inside
it is already batched. A pool of our own would re-import the libraries and re-load the models per
worker, and read as "batching did not help" - hard rule 3.

Nothing is cached (CLAUDE.md rule 9). The call is deterministic in ``(species, cdr3, v, j)``, so it
runs once per distinct key within the build and joins back - that is rule 4's deduplication, not a
stored result. Measured on this repository's own chunk files: 190,902 distinct keys in **23.3 s, one
process, 122 us per key**, against the 15.96 s + 13.16 s + 7.71 s of the three stages it replaces.
"""
from __future__ import annotations

import polars as pl

#: VDJdb species names are passed through unchanged to :func:`infer`: ``annotate_junctions`` resolves
#: them itself, and its last rung is a germline scaffold rather than a fitted model, so a species with
#: no published model still answers. Nothing in this stage decides which model to use.
_SPECIES_IS_PASSED_THROUGH = True

#: VDJdb species -> (``vdjtools`` model source, organism), for the one caller that still needs a
#: **model handle** rather than an annotation: :mod:`vdjdb.annotate.contig` stitches a full variable
#: domain and has to hand ``stitch_contig`` a loaded model. Both sources are ``vdjtools``' own fits
#: and differ in how they were fit - ``olga`` is the OLGA bootstrap, human-only and the community
#: reference for human, and ``arda`` is an EM fit over real non-functional reads in arda's IMGT
#: allele namespace, the only bundled set covering mouse. There is no "arda model": arda is the
#: aligner and the germline namespace. A species absent here gets no stitched contig, because a
#: guess from the wrong organism's marginals would be worse than an empty cell.
MODELS: dict[str, tuple[str, str]] = {
    "HomoSapiens": ("olga", "human"),
    "MusMusculus": ("arda", "mouse"),
}


def _resolver(model: object, table: str, column: str) -> dict[str, str | None]:
    """A call -> the allele to pass the model, or ``None`` to marginalise over that segment.

    The model is keyed by allele and raises on anything else, so every call is resolved before the
    run rather than by catching an exception per sequence. Three outcomes, none of which invents a
    call: a known allele is used as given; a gene name whose model has exactly one allele resolves
    to that allele, because there is no choice to make; anything else is absent from the mapping, so
    ``.get()`` yields ``None`` and the model marginalises over that segment.
    """
    genes = model.genomic[table]                                   # type: ignore[attr-defined]
    alleles = {a: a for a in genes[column]}
    sole = dict(genes.join(genes.group_by("gene").len().filter(pl.col("len") == 1),
                           on="gene", how="inner").select("gene", column).iter_rows())
    return {"": None, **{g: a for g, a in sole.items() if g not in alleles}, **alleles}

#: The columns this stage adds to ``chains``. ``d.posterior`` belongs here rather than in a module of
#: its own because it comes from the same scenario that produced ``d.inferred``: the number beside a
#: call has to be the probability of *that* call, and taking it from a second estimator is what the
#: retired `annotate/dgene.py` did (78.5 % agreement with the call it annotated).
NT_COLUMNS: tuple[str, ...] = ("cdr3nt", "cdr3nt.pgen", "cdr3nt.margin",
                               "d.inferred", "d.start", "d.end", "d.posterior",
                               "v.end.inferred", "j.start.inferred")

#: ``annotate_junctions`` takes what it calls ``junction_aas``, and means it - Cys104..Phe/Trp118
#: inclusive, which is what VDJdb's ``cdr3`` column holds. Passing true CDR3s would infer junctions
#: two codons short, with no error. See :mod:`vdjdb.convert.coords`.
_JUNCTION_IS_WHAT_IT_WANTS = True

#: What an unmapped V/J boundary reads as, the value ``v.end``/``j.start`` already carry where the
#: markup engine declines, so the fallback columns use it too rather than inventing a second marker.
UNMAPPED = -1

#: The key the annotation is deduplicated on and joined back by.
KEY: tuple[str, ...] = ("species", "cdr3", "v.segm", "j.segm")


def blank_columns(keys: pl.DataFrame) -> pl.DataFrame:
    """``keys`` plus :data:`NT_COLUMNS`, all empty. The answer for an empty corpus."""
    return keys.with_columns(pl.lit("").alias("cdr3nt"),
                             pl.lit(None, pl.Float64).alias("cdr3nt.pgen"),
                             pl.lit(None, pl.Float64).alias("cdr3nt.margin"),
                             pl.lit("").alias("d.inferred"),
                             pl.lit(None, pl.Int64).alias("d.start"),
                             pl.lit(None, pl.Int64).alias("d.end"),
                             pl.lit(None, pl.Float64).alias("d.posterior"),
                             pl.lit(UNMAPPED, pl.Int64).alias("v.end.inferred"),
                             pl.lit(UNMAPPED, pl.Int64).alias("j.start.inferred"))


def infer(keys: pl.DataFrame) -> pl.DataFrame:
    """``(species, cdr3, v.segm, j.segm)`` -> the same frame plus :data:`NT_COLUMNS`.

    ``keys`` must already be distinct and sorted. One row comes back per row in, in input order,
    with nulls where the model cannot explain a junction, so the result is joined **positionally**
    and the worker count cannot change the answer - there is no worker count (CLAUDE.md rule 7).
    """
    from vdjtools.model import annotate_junctions

    from .cdr3fix import ensure_reference

    if keys.is_empty():
        return blank_columns(keys)
    # Without this, arda mistakes this repository for its own checkout, every germline lookup
    # returns nothing and the whole corpus comes back `impossible` with no error. See
    # `annotate/cdr3fix.py::ensure_reference`, and note that it is silent: a misconfiguration
    # becomes a database of wrong annotations.
    ensure_reference()
    got = annotate_junctions(keys["cdr3"].to_list(),
                             keys["v.segm"].to_list(),
                             keys["j.segm"].to_list(),
                             species=keys["species"].to_list())
    # **The inferred nucleotides must encode the junction that ships, and on 2,214 keys they do
    # not.** `cdr3_nt` comes back translating to a *substituted* sequence - `CASSSRAGGEQYF` ships
    # and `TGT...` translates to `CASSTRAGGEQYF`, one residue different - because the library repairs
    # with arda's default `max_replace=1` where `annotate/cdr3fix.py` deliberately uses 0 (#327:
    # substituting destroys the evidence that fixes an allele call). `cdr3_repaired` does **not**
    # report it: it equals the input on all but 117 of these rows, so a proxy check misses them.
    # Filed as `antigenomics/vdjtools#187`.
    #
    # So the gate is the direct check, the same one `tests/release/test_tables_contract.py` makes:
    # everything read off the nucleotide sequence is dropped where the sequence does not translate
    # to its own junction. A coordinate into a sequence nobody can see is worse than an empty cell.
    from vdjtools.model import translate
    encodes = pl.Series("__ok", [nt is not None and nt != "" and translate(nt) == aa
                                 for nt, aa in zip(got["cdr3_nt"], keys["cdr3"], strict=True)],
                        dtype=pl.Boolean)

    def keep(col: pl.Expr | pl.Series) -> pl.Expr:
        """``col`` where the nucleotides encode their own junction, null where they do not."""
        return pl.when(encodes).then(col).otherwise(None)

    return keys.with_columns(
        keep(got["cdr3_nt"]).fill_null("").alias("cdr3nt"),                 # rule 6
        keep(got["pgen"]).alias("cdr3nt.pgen"),
        # How far the winner beat the runner-up. No runner-up means nothing else was in contention,
        # which is an unbounded margin rather than a missing one.
        keep(pl.when(got["runner_up_pgen"].is_null() | (got["runner_up_pgen"] == 0))
               .then(pl.when(got["pgen"].is_null()).then(None).otherwise(float("inf")))
               .otherwise(got["pgen"] / got["runner_up_pgen"])).alias("cdr3nt.margin"),
        keep(got["d_call"]).fill_null("").alias("d.inferred"),
        # **The published contract is 0-based half-open in `cdr3nt`'s space** (`docs/outputs.md`),
        # and the library reports 1-based closed in junction space, so the start loses one and the
        # end is already the half-open bound. Converting here rather than changing the contract keeps
        # every consumer of the shipped table working; `convert.coords` records the four spaces.
        keep(got["d_start_nt"].cast(pl.Int64) - 1).alias("d.start"),
        keep(got["d_end_nt"].cast(pl.Int64)).alias("d.end"),
        # The probability of the gene `d.inferred` names, from the same scenario weights that named
        # it. Null where the locus has no D, which is an answer rather than a gap.
        keep(got["d_posterior"]).alias("d.posterior"),
        # The nucleotide-derived V/J boundary, in residues. The library reads it off the inferred
        # nucleotide junction and recomputes the residue bound from it, so there is no `ceil(nt/3)`
        # conversion in this repository any more. `add_junction_nt` masks these against the markup
        # engine's own answer, so they only ever fill where that declined.
        *_fallback_boundaries(keys, got))


def _fallback_boundaries(keys: pl.DataFrame, got: pl.DataFrame) -> list[pl.Expr]:
    """``v.end.inferred`` / ``j.start.inferred``: an **ungated germline alignment**, per locus.

    This is deliberately a *second, different* computation, and that is the whole reason the columns
    exist. Where the markup declines a boundary it declines it as ``impossible`` - a junction whose
    anchors cannot be satisfied - and measured on this corpus it declines 3,076 V and 1,104 J keys
    that way, on **every one of which** the annotation's own ``v_end_nt`` / ``j_start_nt`` is also
    absent. ``germline_boundary`` asks the germline alignment instead and answers anyway, which is
    what filled 2,446 ``v.end.inferred`` and 488 ``j.start.inferred`` cells before this stage was
    rewritten; taking the boundary from the annotation alone lost all of them.

    It needs a germline to align against, so where our call is unresolvable the allele the model
    settled on goes in - the same allele ``v.inferred`` records. Passing ``None`` through instead
    would decline on chains the scenario answered.

    `Scenario.v_end` is **not** what is asked: maximising P(sequence) explains N-region nucleotides
    as templated whenever it can, so it credits the germline too far
    (`antigenomics/vdjtools#182`, opened from this column's measurement). Measured against
    `isalgo/airr_control`'s `human.trb.ntvj` on 8,132 VDJdb human TRB junctions with unambiguous
    nucleotide truth: `v.end` exact 7,554 of 8,132 (92.89 %) against the history's 7,244 (89.08 %),
    `j.start` 7,966 (97.96 %) against 7,483 (92.02 %). It costs 2.67 us/row.
    """
    from vdjtools.model import germline_boundary, load_bundled

    from ..convert.coords import nt_to_aa_boundary_expr

    v_out = [UNMAPPED] * keys.height
    j_out = [UNMAPPED] * keys.height
    frame = keys.with_row_index("__i").with_columns(
        got["v_call_nt"].alias("__vm"), got["j_call_nt"].alias("__jm"),
        # arda resolves the locus from the junction, which is the only thing here that knows it: a
        # record naming neither V nor J has no locus in its own columns.
        got["locus"].fill_null("").alias("locus"))
    # Per (species, locus), because `germline_boundary` takes one loaded model. A species with no
    # bundled fit gets no fallback, which is what it got before: the primary columns answer it
    # through the annotation's germline-scaffold rung, and a boundary from the wrong organism's
    # marginals would be worse than UNMAPPED.
    for (sp, locus), grp in sorted(frame.filter(pl.col("species").is_in(list(MODELS)))
                                        .group_by("species", "locus", maintain_order=True),
                                   key=lambda kv: kv[0]):   # sorted: the order must not vary
        if not locus:
            continue
        source, organism = MODELS[sp]
        try:
            model = load_bundled(locus, source, organism=organism)
        except (FileNotFoundError, ValueError):
            continue
        vmap = _resolver(model, "genes_v", "v_allele")
        jmap = _resolver(model, "genes_j", "j_allele")
        v = [vmap.get(a) or m for a, m in zip(grp["v.segm"], grp["__vm"], strict=True)]
        j = [jmap.get(a) or m for a, m in zip(grp["j.segm"], grp["__jm"], strict=True)]
        bounds = germline_boundary(model, grp["cdr3"].to_list(), v=v, j=j)
        conv = grp.select(
            nt_to_aa_boundary_expr(bounds["v_end"], unmapped=UNMAPPED).alias("v"),
            nt_to_aa_boundary_expr(bounds["j_start"], unmapped=UNMAPPED).alias("j"))
        for i, vv, jj in zip(grp["__i"], conv["v"], conv["j"], strict=True):
            v_out[i], j_out[i] = (vv if vv is not None else UNMAPPED), \
                                 (jj if jj is not None else UNMAPPED)
    return [pl.Series("v.end.inferred", v_out, dtype=pl.Int64),
            pl.Series("j.start.inferred", j_out, dtype=pl.Int64)]


def add_junction_nt(chains: pl.DataFrame, records: pl.DataFrame) -> pl.DataFrame:
    """Add :data:`NT_COLUMNS` to ``chains``. One batched call over the distinct key set."""
    species = records.select("record_id", "species")
    keyed = chains.join(species, on="record_id", how="left")
    # **A blank V or J is passed through, not filtered out.** `arda.cdr3fix` proposes a call for a
    # side the submission left blank (2.34) and the locus too (2.36), so those rows are exactly the
    # ones the inference is for. Requiring both calls cost `v.end.inferred` 2,418 of the 2,442 cells
    # it was measured on, because a chain with no V was never offered to the model at all.
    keys = keyed.filter(pl.col("cdr3") != "").select(KEY).unique().sort(KEY)
    lookup = infer(keys)
    return (keyed.join(lookup, on=list(KEY), how="left")
            .with_columns(pl.col("cdr3nt", "d.inferred").fill_null(""))   # rule 6
            # **A VJ locus has no D.** arda resolves the locus from the junction, so a junction
            # filed in the alpha column that reads as TRD comes back with a TRD D gene - 20 chains,
            # every one `CxxD...KLx`, TRDV2's own anchor. The D is a real reading of the sequence and
            # a wrong column to ship it in, so it is dropped and `d.segm.arda` is where the
            # disagreement can be read. `chains` promises only beta chains carry a D
            # (`tests/release/test_tables_contract.py`).
            .with_columns(*[pl.when(pl.col("gene") == "TRA").then(blank).otherwise(pl.col(col))
                              .alias(col)
                            for col, blank in (("d.inferred", pl.lit("")),
                                               ("d.start", pl.lit(None, pl.Int64)),
                                               ("d.end", pl.lit(None, pl.Int64)),
                                               ("d.posterior", pl.lit(None, pl.Float64)))])
            # A fallback, never an override (#631). The nucleotide boundary survives only where the
            # markup engine declined; where the markup answered, the column reads UNMAPPED and a
            # consumer coalescing the two cannot overwrite a markup answer even by accident. That is
            # stronger than filling the shipped column, and it is what keeps the legacy export -
            # which reads `v.end` and `j.start` - byte-identical.
            .with_columns(*[
                pl.when(pl.col(shipped) == UNMAPPED)
                  .then(pl.col(fallback).fill_null(UNMAPPED))
                  .otherwise(pl.lit(UNMAPPED, pl.Int64)).alias(fallback)
                for shipped, fallback in (("v.end", "v.end.inferred"),
                                          ("j.start", "j.start.inferred"))])
            .drop("species")
            .sort("record_id", "gene"))
