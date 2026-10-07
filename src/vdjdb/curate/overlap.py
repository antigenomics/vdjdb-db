"""Distinct-junction overlap screens for tracing experimental provenance."""
from __future__ import annotations

import polars as pl

STRATUM = ["species", "antigen.epitope", "mhc.a", "mhc.b", "mhc.class"]
GROUP = ["chunk.file", "reference.id"]
MIN_SHARED = 5
MIN_CONTAINMENT = 0.20


def overlaps(records: pl.DataFrame, *, by_sample: bool = False,
             submitted: bool = False) -> pl.DataFrame:
    """All group pairs, including zero matches, within each species/pMHC and chain mode.

    Sizes count distinct junctions, not rows or cells. Blank donor IDs are unknown;
    identifiers from different references are never treated as the same donor.
    """
    group = [*GROUP, *(["meta.subject.id"] if by_sample else [])]
    df = records.with_columns(*(pl.lit("").alias(c) for c in [*STRATUM, *group]
                               if c not in records.columns))
    results = []
    for mode, chains in [("alpha", ["alpha"]), ("beta", ["beta"]),
                         ("paired", ["alpha", "beta"]),
                         ("paired.pmhc", ["alpha", "beta"])]:
        junctions = [f"cdr3.{chain}" for chain in chains]
        stratum = ["species"] if mode == "paired.pmhc" else STRATUM
        keys = [*junctions, *STRATUM[1:]] if mode == "paired.pmhc" else junctions
        if submitted:
            df_mode = df.with_columns(*(
                pl.col(f"__cdr3old.{chain}").alias(f"cdr3.{chain}")
                for chain in chains if f"__cdr3old.{chain}" in df.columns))
        else:
            df_mode = df
        unique = (df_mode.filter(pl.all_horizontal(pl.col(c) != "" for c in junctions))
                  .select(*stratum, *group, *keys).unique())
        counts = unique.group_by(*stratum, *group).len(name="n")
        # Stable lexical group IDs avoid both self-pairs and reverse pairs.
        groups = unique.select(group).unique().sort(group).with_row_index("group")
        unique = unique.join(groups, on=group, validate="m:1")
        counts = counts.join(groups, on=group, validate="m:1")
        pairs = counts.join(counts, on=stratum, suffix=".other").filter(
            pl.col("group") < pl.col("group.other"))
        matches = (unique.select(*stratum, *keys, "group")
                   .join(unique.select(*stratum, *keys, "group"),
                         on=[*stratum, *keys], suffix=".other")
                   .filter(pl.col("group") < pl.col("group.other"))
                   .group_by(*stratum, "group", "group.other").len(name="shared"))
        pairs = (pairs.join(matches, on=[*stratum, "group", "group.other"], how="left")
                 .with_columns(pl.col("shared").fill_null(0), pl.lit(mode).alias("mode"))
                 .with_columns(
                     (pl.col("n").cast(pl.UInt64) * pl.col("n.other").cast(pl.UInt64))
                     .alias("size.product"),
                     (pl.col("shared") / pl.min_horizontal("n", "n.other"))
                     .alias("containment"),
                     (pl.col("reference.id") == pl.col("reference.id.other"))
                     .alias("same.reference"))
                 .with_columns(
                     ((pl.col("shared") + 1) / (pl.col("size.product") + 1))
                     .alias("ratio.plus1"),
                     ((pl.col("shared") >= MIN_SHARED)
                      & (pl.col("containment") >= MIN_CONTAINMENT))
                     .alias("review.provenance"))
                 .drop("group", "group.other"))
        if mode == "paired.pmhc":
            pairs = pairs.with_columns(*(pl.lit("*").alias(c) for c in STRATUM[1:]))
        results.append(pairs.select(*STRATUM, *group, "n", *[f"{c}.other" for c in group],
                                   "n.other", "shared", "mode", "size.product", "containment",
                                   "same.reference", "ratio.plus1", "review.provenance"))
    return pl.concat(results).sort(*STRATUM, "mode", *group, *[f"{c}.other" for c in group])


def report(pairs: pl.DataFrame, files: list[str]) -> str:
    """Compact candidate list; the complete TSV includes the zero-overlap pairs."""
    flagged = pairs.filter(pl.col("review.provenance") &
                           (pl.col("chunk.file").is_in(files)
                            | pl.col("chunk.file.other").is_in(files)))
    out = ["### Experimental provenance overlap", "",
           "Distinct junction sets within species and peptide-MHC (paired.pmhc compares distinct "
           "receptor-pMHC observations across peptides). High overlap means trace the "
           "source experiment; it does not establish independent validation or justify deletion.", "",
           "Review candidates have at least 5 shared junctions and 20% containment of the smaller set. "
           "The smoothed ratio is (shared + 1)/(n1*n2 + 1), a descriptive measure, not a p-value.", "",
           "| Reference 1 | Reference 2 | Peptide | Mode | n1 | n2 | Shared | Ratio +1 |",
           "|---|---|---|---|---:|---:|---:|---:|"]
    for row in flagged.iter_rows(named=True):
        out.append(f"| {row['reference.id']} | {row['reference.id.other']} | "
                   f"{row['antigen.epitope']} | {row['mode']} | {row['n']} | "
                   f"{row['n.other']} | {row['shared']} | {row['ratio.plus1']:.3g} |")
    if flagged.is_empty():
        out += ["", "No group exceeds the screening threshold. Small reused datasets can still escape it."]
    return "\n".join([*out, ""])
