"""TCRNET: which clonotypes sit in a denser neighbourhood than the background explains.

Stage one of the motif pipeline. Per ``(species, gene, epitope)`` it counts, for every unique CDR3,
how many neighbours it has **within the epitope-specific sample** and how many it has in a matched
background repertoire, and keeps the ones whose within-sample degree the background cannot account
for. Stage two (:mod:`vdjdb.motifs.cluster`) turns those into clusters; stage three
(:mod:`vdjdb.motifs.pwm`) turns each cluster into a logo.

Two things here are deliberate and neither is obvious from the call site.

**``vdjtools.overlap.tcrnet`` is used as a neighbour counter, not as a test.** Its ``_score``
computes ``E = (n_target / max(m_control, 1)) * n_control`` with no pseudocount, and
``scipy.stats.poisson.sf(k - 1, 0.0)`` is **exactly 0.0** for every ``k >= 1``. So every clonotype
with a within-sample neighbour and no background neighbour gets ``p_enrichment == 0.0``,
``q_value == 0.0``, and sorts to the top. Measured on human TRB, 111,407 unique CDR3s against the
bundled 250k control: ``n_control == 0`` for 81.0 % of queries and 44,010 rows (39.5 %) get exactly
zero. A bigger background does not fix it -- mean ``n_control`` scales linearly with ``M``, so ``E``
is ``M``-invariant in expectation and ``M`` only controls the zero-inflation. Reproduced here on
GILGFVFTL: 4,718 of 6,637. ROADMAP section 8.1; filed upstream against ``vdjtools``.

What is used instead is :func:`legacy_pvalue`, which is bit-faithful to the Groovy
``DegreeStatisticsAnnotator.computePValue`` and carries the pseudocount that removes the pathology.

**The background is always passed in.** Left to itself ``tcrnet()`` calls
``vdjmatch.evalue.background(locus, species)`` with no ``size`` and indexes the entire table.
:func:`control_for` loads a bounded, seeded, uniformly reservoir-sampled control instead, and
:func:`enrich` reads ``M`` back off the index it was given.
"""
from __future__ import annotations

import polars as pl

from ..config import SEED

#: Unique clonotypes **requested** for every background, matching the scale of the 1M
#: ``vdjdb-web-control`` subsamples the legacy pipeline used. Asking for one number rather than
#: "the whole table" is what stops ``M`` being whatever each source happens to hold -- human TRA
#: alone is 2,266,274. Where the source is smaller the whole of it is used, so the realised ``M``
#: is, measured 2026-09-25 after the productive filter: human TRA 1,000,000 · human TRB 1,000,000 ·
#: mouse TRB 694,241 · mouse TRA 272,827. :func:`enrich` reads it from the index rather than
#: assuming this constant, because ``M`` enters the statistic directly.
CONTROL_SIZE = 1_000_000

#: ``(species, gene)`` -> the ``seqtree.control`` name. Amino-acid builds, never ``ntvj``: the query
#: here is an epitope sample deduplicated to unique CDR3aa, and an nt-level control inserts the same
#: ``cdr3aa`` once per nt variant -- measured at 1.63x inflation for VDJdb-matching clonotypes
#: against 1.04x overall, concentrated on exactly the germline-proximal public sequences under test.
#: That would inflate the background degree and systematically *under*-call public motifs.
#: ROADMAP section 8.7. Species with no background get no motifs, as in the legacy pipeline.
CONTROLS: dict[tuple[str, str], str] = {
    ("HomoSapiens", "TRA"): "human_tra_aa",
    ("HomoSapiens", "TRB"): "human_trb_aa",
    ("MusMusculus", "TRA"): "mouse_tra_aa",
    ("MusMusculus", "TRB"): "mouse_trb_aa",
}

#: ``(species, gene)`` -> the ``isalgo/airr_control`` asset the **background PWM** is read from.
#: The same files the controls above are sampled from; :func:`background_frame` needs the ``v`` and
#: ``j`` calls, which a built ``seqtree.Index`` no longer carries.
ASSETS: dict[tuple[str, str], str] = {
    ("HomoSapiens", "TRA"): "human.tra.aa.vdjtools.tsv.gz",
    ("HomoSapiens", "TRB"): "human.trb.aa.vdjtools.tsv.gz",
    ("MusMusculus", "TRA"): "mouse.tra.aa.vdjtools.tsv.gz",
    ("MusMusculus", "TRB"): "mouse.trb.aa.vdjtools.tsv.gz",
}

#: vdjmatch edit-distance scope ``"subs,ins,dels,total"``. The legacy's is one substitution; see
#: :data:`TUNED`, where TRA takes two.
SCOPE = "1,0,0,1"

#: Threshold on the enrichment p-value. ⚠ **Applied to the raw p, not a BH-adjusted q**, because
#: that is what the legacy does: `compute_vdjdb_motifs.Rmd` writes `mutate(p.adj = p.value.g)`,
#: which is an identity -- the name says adjusted and the code adjusts nothing. Measured, the
#: difference is not small in the direction anyone expects: over 185,738 scored clonotypes BH
#: `q <= 0.05` calls **53,609** enriched where raw `p <= 0.05` calls **48,418**, because the
#: p-value distribution is bottom-heavy enough that the step-up threshold rises well above 0.05.
#: `q.legacy` is computed and carried regardless, so switching is a one-line change and a
#: measurement (ROADMAP section 30.3).
P_THRESHOLD = 0.05

#: Within-sample neighbours a clonotype needs before it can be called enriched, independent of its
#: p-value: `enriched = degree.s >= 2 & p.adj < 0.05`. A single neighbour is not a neighbourhood.
MIN_DEGREE = 2

#: An epitope sample smaller than this many **unique clonotypes** is not scored. The legacy's
#: `dt.epi.count %>% filter(total >= 10)`. The benchmark cohort's floor of 30 is a different
#: number for a different purpose and using it here cost 13 epitopes their motifs.
MIN_SAMPLE = 10

#: **Per-chain parameters, fitted rather than inherited.** Chosen by maximising retention subject to
#: purity **and** precision staying at or above the shipped TCRNET's, scored in
#: :mod:`vdjdb.validate.motif_bench` on one cohort -- the stated acceptance criterion applied
#: mechanically over `scope` x `p` x `min_degree` x `min_cluster` (ROADMAP section 34).
#:
#: ===== ================================ ========= ========= ========= =========
#: chain config                           retention legacy    purity    precision
#: ===== ================================ ========= ========= ========= =========
#: TRA   scope 2,0,0,2 · p .05 · mc 3     **0.4149**  0.2105    0.8905    0.8890
#: TRB   scope 1,0,0,1 · p .01 · mc 5     **0.3337**  0.3218    0.9790    0.9761
#: ===== ================================ ========= ========= ========= =========
#:
#: TRA gains **+97 % retention with purity and precision both up** -- a strict win on every axis.
#: TRB's bar is tight (the shipped file's purity is 0.9790) and the winner trades a little retention
#: for it: ⚠ the previous default `p = 0.05` reached 0.3382 but at purity **0.9781**, which *fails*
#: the criterion by 0.0009. Tightening `p` to 0.01 is what buys the bar back.
#:
#: ⚠ TRA's two-substitution neighbourhood is a **different definition of neighbour**, not a tuned
#: threshold -- it is here because it dominates on every measured axis, including removing 95 % of
#: the reproduction regression (842 lost clonotypes -> 44, ROADMAP section 33.2).
TUNED: dict[str, dict] = {
    "TRA": {"scope": "2,0,0,2", "p": 0.05, "min_degree": 2, "min_cluster": 3},
    "TRB": {"scope": "1,0,0,1", "p": 0.01, "min_degree": 2, "min_cluster": 5},
}


def control_for(species: str, gene: str, size: int = CONTROL_SIZE, seed: int = SEED):
    """The background index for one ``(species, gene)``, or ``None`` if there is no background.

    Delegates to ``seqtree.control.load_control``, which streams from ``isalgo/airr_control``,
    filters to the productive 20 and reservoir-samples **uniformly over unique clonotypes** rather
    than taking the abundance-sorted head. What it stores is the fetched table, which is an input,
    not a computed result (CLAUDE.md hard rule 9); the sample it draws is deterministic in the seed.
    """
    name = CONTROLS.get((species, gene))
    if name is None:
        return None
    from seqtree.control import load_control
    return load_control(name, size=size, seed=seed)


def legacy_pvalue(degree: pl.Expr, n_control: pl.Expr, n_sample: int, m_control: int) -> pl.Expr:
    """The legacy enrichment p-value, as an expression over a tcrnet result frame.

    Bit-faithful to Groovy ``DegreeStatisticsAnnotator.computePValue``::

        p = (n_control + 1) / (M + 1)
        P = binom.sf(degree - 1, N, p) / (1 - (1 - p) ** N)

    ``degree`` is the clonotype's within-sample neighbour count, ``n_control`` its background
    neighbour count, ``N`` (``n_sample``) the unique clonotypes in this sample and ``M``
    (``m_control``) the unique clonotypes in the background.

    The ``+1`` numerator is the whole point: it is what stops a clonotype with no background
    neighbour from being handed a probability of zero. The denominator conditions on the clonotype
    having at least one neighbour, which is the event that made it a candidate at all.
    """
    from scipy.stats import binom

    def _sf(d: pl.Series, nc: pl.Series) -> pl.Series:
        p = (nc.to_numpy() + 1.0) / (m_control + 1.0)
        tail = binom.sf(d.to_numpy() - 1, n_sample, p)
        return pl.Series(tail / (1.0 - (1.0 - p) ** n_sample))

    return pl.map_batches([degree, n_control], lambda s: _sf(s[0], s[1]), return_dtype=pl.Float64)


def _bh(p: pl.Expr) -> pl.Expr:
    """Benjamini-Hochberg step-up, as a monotone expression -- no ``statsmodels`` for six lines."""
    n = pl.len()
    rank = p.rank("ordinal")
    return (p * n / rank).reverse().cum_min().reverse().clip(upper_bound=1.0)


def enrich(sample: pl.DataFrame, control, *, scope: str = SCOPE) -> pl.DataFrame:
    """Score one epitope sample against ``control``. One row per unique clonotype.

    ``sample`` is the canonical clonotype frame (``junction_aa``, ``v_call``, ``j_call``,
    ``duplicate_count``). Adds ``p.legacy`` and ``q.legacy`` beside the columns ``tcrnet()``
    returns; ``p_enrichment`` and ``q_value`` are carried through **unused**, for the deviation
    report to quote.
    """
    from vdjtools.overlap import tcrnet

    # One call for the whole sample. Never once per clonotype, and never a pool around it: the
    # native search threads internally (CLAUDE.md hard rule 3).
    scored = tcrnet(sample, control=control, scope=scope)
    n, m = scored.height, len(control)
    return scored.with_columns(
        legacy_pvalue(pl.col("n_neighbors"), pl.col("n_control"), n, m).alias("p.legacy"),
    ).with_columns(
        _bh(pl.col("p.legacy")).alias("q.legacy"),
    )


def enriched_clonotypes(chains: pl.DataFrame, records: pl.DataFrame, *,
                        scope: str | None = None, p: float | None = None,
                        min_degree: int | None = None, min_sample: int = MIN_SAMPLE,
                        control_size: int = CONTROL_SIZE,
                        control_seed: int = SEED) -> pl.DataFrame:
    """Every scored clonotype, flagged ``enriched`` where the background cannot explain its degree.

    Groups by ``(species, gene, antigen.epitope)`` -- the scope the legacy pipeline used, which
    ``CalcDegreeStats.groovy`` reached with ``-g dummy`` rather than the V/VJ/VJL grouping the plan
    assumed, so the grouping-free ``tcrnet()`` is an exact match and not a deviation
    (ROADMAP section 8.2).

    One background index is loaded per ``(species, gene)`` and reused across that chain's epitopes;
    loading it per epitope would re-index a million clonotypes a few hundred times.
    """
    df = _samples(chains, records)
    out: list[pl.DataFrame] = []
    for (species, gene), chain in df.group_by(["species", "gene"], maintain_order=True):
        control = control_for(species, gene, control_size, control_seed)
        if control is None:
            continue
        # Per-chain scope and threshold: the two chains do not have the same optimum, and pinning
        # them to one number costs TRA 97 % of its retention (TUNED).
        tuned = TUNED.get(gene, {})
        chain_scope = tuned.get("scope", scope) if scope is None else scope
        chain_p = tuned.get("p", p) if p is None else p
        chain_degree = tuned.get("min_degree", min_degree) if min_degree is None else min_degree
        for (epitope,), grp in chain.group_by(["antigen.epitope"], maintain_order=True):
            sample = grp.select("junction_aa", "v_call", "j_call", "duplicate_count")
            if sample.height < min_sample:
                continue
            scored = enrich(sample, control, scope=chain_scope)
            # Every scored clonotype is returned, flagged -- not only the ones that pass. The graph
            # stage recruits a clonotype that is a neighbour of an enriched one even when it is not
            # itself enriched, which is what the legacy Rmd's two-stage `compute_edges` does.
            out.append(scored.with_columns(
                ((pl.col("p.legacy") <= chain_p) & (pl.col("n_neighbors") >= chain_degree))
                .alias("enriched"),
                pl.lit(chain_scope).alias("scope"),
                pl.lit(species).alias("species"), pl.lit(gene).alias("gene"),
                pl.lit(epitope).alias("antigen.epitope")))
    if not out:
        return pl.DataFrame(schema={"species": pl.Utf8, "gene": pl.Utf8,
                                    "antigen.epitope": pl.Utf8, "junction_aa": pl.Utf8,
                                    "enriched": pl.Boolean, "scope": pl.Utf8})
    # Sorted, because group_by order is the frame's and the caller keys on this (hard rule 7).
    return pl.concat(out, how="vertical").sort(
        "species", "gene", "antigen.epitope", "junction_aa", "v_call", "j_call")


def _samples(chains: pl.DataFrame, records: pl.DataFrame) -> pl.DataFrame:
    """The per-epitope clonotype samples, deduplicated on ``(cdr3, v, j)``.

    ``duplicate_count`` is the number of **records** reporting that clonotype against that epitope,
    which is what makes a public clonotype heavier than a singleton in the cluster size. VDJdb's
    ``cdr3`` is junction space (Cys104..Phe/Trp118 inclusive), which is exactly what vdjtools calls
    ``junction_aa`` -- no conversion, and none wanted (CLAUDE.md, the coordinate table).
    """
    return (
        chains
        .join(records.select("record_id", "species", "antigen.epitope"), on="record_id")
        # A clonotype with no V or J call is not a clonotype for this purpose: the legacy requires
        # `v.segm != "", j.segm != "", !is.na(j.segm)` before it builds the sample, and the cluster's
        # representative alleles are read off these columns.
        .filter((pl.col("cdr3") != "") & (pl.col("v.segm") != "") & (pl.col("j.segm") != ""))
        .group_by(["species", "gene", "antigen.epitope",
                   pl.col("cdr3").alias("junction_aa"),
                   pl.col("v.segm").alias("v_call"), pl.col("j.segm").alias("j_call")],
                  maintain_order=True)
        .agg(pl.len().cast(pl.Int32).alias("duplicate_count"))
        .sort("species", "gene", "antigen.epitope", "junction_aa", "v_call", "j_call")
    )


def background_frame(species: str, gene: str) -> pl.DataFrame | None:
    """The background repertoire as ``(cdr3aa, v, j)``, for the PWM prior. ``None`` if there is none.

    A **download**, not a cache: the asset is an input that arrives over the network and is
    content-addressed by ``huggingface_hub`` (CLAUDE.md hard rule 9). Only the three columns the
    PWM needs are read -- the tables carry fourteen and run to 424 MB compressed.

    Unlike the enrichment control this is **not** subsampled. The two are different statistics: the
    control fixes ``M`` in a tail probability, where a bounded, seeded draw is what makes the number
    comparable across chains; the PWM prior is a frequency estimate, where every row helps and none
    of them enters a p-value.
    """
    asset = ASSETS.get((species, gene))
    if asset is None:
        return None
    from huggingface_hub import hf_hub_download

    path = hf_hub_download("isalgo/airr_control", asset, repo_type="dataset")
    return pl.read_csv(path, separator="\t", columns=["cdr3aa", "v", "j"], infer_schema_length=0)
