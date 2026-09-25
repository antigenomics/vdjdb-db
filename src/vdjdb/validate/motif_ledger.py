"""What a motif rebuild preserved, lost and gained, clonotype by clonotype.

The aggregate metrics in :mod:`vdjdb.validate.motif_bench` say whether a clustering is *better*.
They do not say what happened to any particular motif, and "better on average" is not an answer to
"where did that cluster go". This is the difference ledger of section 6 applied to motifs.

**Lost and new are each two different things, and conflating them is the mistake this module
exists to avoid.** The reference release was built from an older corpus, so a clonotype can be
absent from our clustering because the clustering dropped it *or* because curation removed or
respelled it, and present in ours because the clustering found it *or* because the chunk is new.
:func:`ledger` separates those four cases by intersecting with each side's **corpus**, not only with
each side's clustered set:

===================  ==========================================================================
category             meaning
===================  ==========================================================================
``preserved``        clustered by both
``lost.unclustered`` in our corpus, clustered in the reference, **not** clustered by us
``lost.absent``      clustered in the reference, not in our corpus at all -- curation, not motifs
``new.recovered``    in the reference corpus and unclustered there, clustered by us
``new.data``         not in the reference corpus -- chunks added since
===================  ==========================================================================

``lost.unclustered`` is the only category that is a regression. ``lost.absent`` is a curation
change and belongs in the chunk ledger, not here.
"""
from __future__ import annotations

import polars as pl

#: The clonotype identity a motif file names. Matches
#: :data:`vdjdb.validate.motif_bench.KEY` so the two instruments agree about what a row is.
KEY: tuple[str, ...] = ("species", "gene", "cdr3aa", "v.segm", "j.segm")


def _clonotypes(members: pl.DataFrame) -> pl.DataFrame:
    return members.select(*KEY).unique(maintain_order=True)


def corpus(chains: pl.DataFrame, records: pl.DataFrame) -> pl.DataFrame:
    """Every clonotype the build could have clustered, as :data:`KEY`."""
    return (chains.join(records.select("record_id", "species"), on="record_id")
            .filter((pl.col("cdr3") != "") & (pl.col("v.segm") != "") & (pl.col("j.segm") != ""))
            .select("species", "gene", pl.col("cdr3").alias("cdr3aa"), "v.segm", "j.segm")
            .unique(maintain_order=True))


def ledger(reference: pl.DataFrame, candidate: pl.DataFrame,
           reference_corpus: pl.DataFrame, candidate_corpus: pl.DataFrame) -> pl.DataFrame:
    """One row per clonotype in either clustering, labelled with its category.

    ``reference`` and ``candidate`` are ``cluster_members``-shaped frames; the two corpora say what
    each side could have clustered. Returns the union with ``category`` and both cluster ids, so a
    caller can group by epitope, by size, or by anything else in the file.
    """
    ref, cand = _clonotypes(reference), _clonotypes(candidate)
    rc = reference_corpus.select(*KEY).unique().with_columns(pl.lit(True).alias("__in_ref_corpus"))
    cc = candidate_corpus.select(*KEY).unique().with_columns(pl.lit(True).alias("__in_cand_corpus"))

    both = (pl.concat([ref, cand], how="vertical").unique(maintain_order=True)
            .join(ref.with_columns(pl.lit(True).alias("__ref")), on=KEY, how="left")
            .join(cand.with_columns(pl.lit(True).alias("__cand")), on=KEY, how="left")
            .join(rc, on=KEY, how="left").join(cc, on=KEY, how="left")
            .with_columns(pl.col("__ref", "__cand", "__in_ref_corpus", "__in_cand_corpus")
                          .fill_null(False)))

    return (both.with_columns(
        pl.when(pl.col("__ref") & pl.col("__cand")).then(pl.lit("preserved"))
        .when(pl.col("__ref") & ~pl.col("__in_cand_corpus")).then(pl.lit("lost.absent"))
        .when(pl.col("__ref")).then(pl.lit("lost.unclustered"))
        .when(~pl.col("__in_ref_corpus")).then(pl.lit("new.data"))
        .otherwise(pl.lit("new.recovered")).alias("category"))
        .drop("__ref", "__cand", "__in_ref_corpus", "__in_cand_corpus")
        .sort(*KEY))


def summarise(led: pl.DataFrame) -> pl.DataFrame:
    """Counts per ``(species, gene, category)``, with the share of each side's clustered total."""
    return (led.group_by(["species", "gene", "category"]).len().rename({"len": "clonotypes"})
               .sort("species", "gene", "category"))


def cluster_agreement(reference: pl.DataFrame, candidate: pl.DataFrame) -> pl.DataFrame:
    """For clonotypes clustered by both, how well the two partitions agree.

    Adjusted Rand index and adjusted mutual information over the shared clonotypes, per
    ``(species, gene)``. A high ARI on a small overlap means the motifs that survived kept their
    shape; a low one means they were re-cut.
    """
    from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score

    r = reference.select(*KEY, pl.col("cid").alias("cid.ref")).unique(subset=KEY, keep="first")
    c = candidate.select(*KEY, pl.col("cid").alias("cid.cand")).unique(subset=KEY, keep="first")
    shared = r.join(c, on=KEY, how="inner")
    rows = []
    for (species, gene), grp in shared.group_by(["species", "gene"], maintain_order=True):
        if grp.height < 2:
            continue
        a, b = grp["cid.ref"].to_list(), grp["cid.cand"].to_list()
        rows.append({"species": species, "gene": gene, "shared_clonotypes": grp.height,
                     "ari": round(adjusted_rand_score(a, b), 4),
                     "ami": round(adjusted_mutual_info_score(a, b), 4),
                     "clusters_ref": grp["cid.ref"].n_unique(),
                     "clusters_cand": grp["cid.cand"].n_unique()})
    return pl.DataFrame(rows).sort("species", "gene")
