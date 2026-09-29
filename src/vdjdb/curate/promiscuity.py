"""Which alleles can present each epitope VDJdb records. Issue #372.

**This does not second-guess the allele a publication reported.** The restriction on a record is the
author's finding and the authority for it is the paper. What a presentation predictor adds is the
*other* alleles: given a donor's HLA type, which of VDJdb's epitopes could that donor present,
including pairings no publication has reported. That is a filter for a query, not a correction to a
record, and it is what #372 asks for - "see if our epitopes, linked to, say A*02, can be presented by
other alleles".

The one judgement the pipeline does make about `mhc.a` is whether the allele exists. A blank or
spurious allele is a transcription error rather than a finding, and that check belongs in `vdjdb qc`
against `proofreading/mhc_alleles.tsv.gz`, not here.

**The presentation axis, not the specificity axis.** :meth:`mhcmatch.Store.restriction` ranks how
strongly one allele prefers a peptide *relative to other alleles*, and its own docstring warns it can
band a canonical, widely-shared ligand as "weak" while that allele is the correct restriction -
measured on `FLRGRAYGL` it put `HLA-A*02:01` ahead of `HLA-B*08:01` on a neighbour vote of 0.73
against 0.03, contradicting the paper and 80 references. So the scorer here is
:func:`mhcmatch.predict.build_scorer` with ``background="proteome"``, the NetMHCpan ``%Rank_EL``
analogue, which asks "is this presented at all".

One scorer build serves every epitope and every allele: it is memoised on the store and depends only
on the panel, so the whole table is one pass of ``model.score`` per (epitope, allele) - never one
:func:`mhcmatch.predict.predict_windows` call per pair, which re-resolves the panel and re-derives the
cut every time.

The output is a **committed, reviewed input**, not a build artifact: `mhcmatch` fetches its reference
data from HuggingFace on first use, so scoring inside a build would put a network call in the critical
path and make the answer depend on the day. Same treatment as ``summary/reference_years.tsv`` -
refreshed by its own pull request through ``vdjdb promiscuity``, never written by ``vdjdb build``
(hard rule 9).
"""

from __future__ import annotations

from pathlib import Path
from typing import NamedTuple

import polars as pl

from ..config import Paths
from ..schema import RESTRICTION_COLUMNS

#: The model that produced a row, recorded in the table so `restriction` can carry it per row.
VERSION_COLUMN = "mhcmatch.version"

#: Written to `proofreading/epitope_promiscuity.tsv`, in this order.
COLUMNS = ("epitope", "mhc_a", "percent_rank", "p_present", "band", "recorded", VERSION_COLUMN)

#: Class I binding cores. An epitope outside this range is not a class I ligand, so the class I
#: scorer has nothing to say about it and it is skipped rather than scored badly.
LENGTHS = range(8, 12)


class Presented(NamedTuple):
    """One epitope scored against one allele of the panel."""

    epitope: str
    mhc_a: str
    percent_rank: float
    p_present: float
    band: str
    #: 1 when VDJdb records this epitope under this allele, 0 when only the predictor pairs them.
    #: A recorded pair is kept whatever it scores - the author reported it, and a reader filtering by
    #: donor type needs to see the record, not the prediction's opinion of it.
    recorded: int
    #: The `mhcmatch` that scored it. On the row rather than in a header because `restriction`
    #: carries it per row (ROADMAP §10.6) and a table refreshed in two passes could hold two.
    mhcmatch_version: str = ""


def presented(recorded: dict[str, set[str]], *, species: str = "human",
              tier: str = "shortlist") -> tuple[list[Presented], list[str]]:
    """Score every epitope in ``recorded`` against the whole class I panel.

    ``recorded`` maps an epitope to the alleles VDJdb records it under, in VDJdb's spelling. Returns
    ``(rows, unmatched)``: the rows a reader can filter by donor allele, and the VDJdb allele
    spellings the panel has no entry for. An unmatched allele is a gap in the panel, not a claim
    about the allele - whether an allele exists is IPD-IMGT/HLA's answer, not a predictor's.

    Rows are every allele that presents the epitope (``band`` above non-binder) plus every allele
    VDJdb records for it. Sorted by epitope then ascending %rank, ties by allele, so the file is
    byte-identical across runs (hard rule 7).

    Network-bound on first call: `mhcmatch` fetches its reference data from HuggingFace. That is why
    this is not in a build.
    """
    import mhcmatch
    from mhcmatch.predict import band_for, build_scorer

    version = getattr(mhcmatch, "__version__", "")

    store = mhcmatch.Store.from_pmhc(tier=tier, species=species)
    model, cal, _ = build_scorer(store, "mhc1", background="proteome")
    panel = sorted(store.panel_alleles("mhc1", "all"))

    # VDJdb's spelling and the panel's are not the same vocabulary; `panel_alleles` is the one place
    # that folds them together, so translate once rather than comparing strings per row.
    spelling: dict[str, str] = {}
    unmatched: list[str] = []
    for allele in sorted({a for alleles in recorded.values() for a in alleles}):
        hit = store.panel_alleles("mhc1", [allele])
        if hit:
            spelling[allele] = hit[0]
        else:
            unmatched.append(allele)

    rows: list[Presented] = []
    for epitope in sorted(recorded):
        if len(epitope) not in LENGTHS:
            continue
        mine = {spelling[a] for a in recorded[epitope] if a in spelling}
        scored: list[Presented] = []
        for allele in panel:
            s = model.score(epitope, allele)
            if s == float("-inf"):
                continue
            rank = cal.percent_rank(allele, s)
            if rank != rank:                      # nan: this allele has no background
                continue
            band = band_for(rank, "mhc1")
            if band == "non-binder" and allele not in mine:
                continue
            scored.append(Presented(epitope, allele, round(rank, 3),
                                    round(cal.p_present(allele, s), 4), band,
                                    1 if allele in mine else 0, version))
        rows.extend(sorted(scored, key=lambda p: (p.percent_rank, p.mhc_a)))
    return rows, unmatched


def _table(root: Path | None = None) -> pl.DataFrame:
    """The committed scores, or an empty frame with the right schema when the file is absent."""
    schema = {"epitope": pl.String, "mhc_a": pl.String, "percent_rank": pl.Float64,
              "p_present": pl.Float64, "band": pl.String, "recorded": pl.Int64}
    path = (root or Paths.discover().root) / "proofreading" / "epitope_promiscuity.tsv"
    if not path.exists():
        return pl.DataFrame(schema=schema)
    got = pl.read_csv(path, separator="\t", schema_overrides=schema)
    if VERSION_COLUMN not in got.columns:
        # Added after the table was first written. An older file carries no version rather than a
        # guessed one: a number whose model nobody can name is worse than a blank.
        got = got.with_columns(pl.lit("").alias(VERSION_COLUMN))
    return got


def annotate(restriction: pl.DataFrame, root: Path | None = None) -> pl.DataFrame:
    """Add :data:`vdjdb.schema.PROMISCUITY_COLUMNS` to ``restriction``. ROADMAP §10.6, phase 16.9.

    A join against a committed, reviewed input and nothing else. The model that produced
    ``proofreading/epitope_promiscuity.tsv`` is never run here: it fetches its reference data from
    HuggingFace, and a build that downloads a model is not offline or deterministic (hard rule 9).
    Re-running it is `vdjdb promiscuity`, which opens a pull request against that table.

    ``alleles.reported`` is the one column that is curation rather than prediction - how many
    distinct ``mhc.a`` values VDJdb itself records for the peptide - so it is counted from
    ``restriction`` and is present on class II rows, where the other four are blank. That is the
    column #372 asks about; the rest say how the prediction ranks what the curators wrote.

    **A curated allele the prediction outranks is not a defect.** Most epitopes here are promiscuous
    - the panel puts several alleles in the strong band - and the curated one is the allele the
    publication typed a donor for, which a proteome-background ranking has no access to. The columns
    are an annotation to sort and filter by, which is why §10.6 keeps them out of every key.
    """
    scores = _table(root)
    reported = (restriction.filter(pl.col("mhc.a") != "")
                           .group_by("antigen.epitope")
                           .agg(pl.col("mhc.a").n_unique().cast(pl.UInt32)
                                  .alias("alleles.reported")))
    if scores.is_empty():
        per_epitope = pl.DataFrame(schema={"antigen.epitope": pl.String, "mhc.a.top": pl.String,
                                           "promiscuity": pl.UInt32,
                                           "mhcmatch.version": pl.String})
        per_allele = pl.DataFrame(schema={"antigen.epitope": pl.String, "mhc.a": pl.String,
                                          "mhc.a.rank": pl.UInt32,
                                          "mhc.a.percentile": pl.Float64})
    else:
        ranked = scores.sort("percent_rank", "mhc_a").with_columns(
            pl.int_range(1, pl.len() + 1).over("epitope").cast(pl.UInt32).alias("rank"))
        per_epitope = (ranked.group_by("epitope")
                       .agg(pl.col("mhc_a").first().alias("mhc.a.top"),
                            (pl.col("band") == "strong").sum().cast(pl.UInt32).alias("promiscuity"),
                            pl.col(VERSION_COLUMN).first().alias(VERSION_COLUMN))
                       .rename({"epitope": "antigen.epitope"}))
        per_allele = (ranked.select(pl.col("epitope").alias("antigen.epitope"),
                                    pl.col("mhc_a").alias("mhc.a.trimmed"),
                                    pl.col("rank").alias("mhc.a.rank"),
                                    pl.col("percent_rank").alias("mhc.a.percentile")))
    # Both sides at two fields. `mhcmatch`'s panel is two-field and there is no deeper groove, so
    # `HLA-A*02:01:48` is scored as the `HLA-A*02:01` molecule it is a third-field allele of - the
    # best available answer, and the alternative is a blank on 74 pairs. Nothing else joins on it.
    return (restriction.with_columns(_trimmed().alias("mhc.a.trimmed"))
                       .join(reported, on="antigen.epitope", how="left")
                       .join(per_epitope, on="antigen.epitope", how="left")
                       .join(per_allele, on=["antigen.epitope", "mhc.a.trimmed"], how="left")
                       .with_columns(pl.col("alleles.reported").fill_null(0),
                                     pl.col("mhc.a.top").fill_null(""),
                                     pl.col(VERSION_COLUMN).fill_null(""))
                       .select(RESTRICTION_COLUMNS))


def _trimmed() -> pl.Expr:
    """``mhc.a`` at two fields, which is the resolution `mhcmatch`'s panel is named at."""
    from mhcmatch.pseudoseq import trim_allele

    return pl.col("mhc.a").map_elements(trim_allele, return_dtype=pl.String)
