"""The evidence table: one row per piece of support for a record.

Long, not wide. A record may have any number of pieces of evidence of any number of kinds -- a
motif cluster per method, a structure per PDB entry, a supporting publication per replication -- and
the wide form would be a table of mostly-empty boolean columns that grows a column per producer.
``vdjdb.parquet`` pivots it back to that shape at export (:mod:`vdjdb.emit.vdjdb3`), in the derived
direction rather than the authored one.

The first producer is independent replication, and it is the same computation as the motif tuning
objective (ROADMAP §11.1). "Two papers found this receptor against this epitope" is both the
strongest evidence the database holds about a record and the signal a clustering has to recover, so
there must not be two implementations of it that can disagree.

It is also the reason deduplication is within a chunk only. A chunk is one paper; collapsing two
chunks' matching rows into one record would delete this signal before it could be counted.

No held-out validation data is ever an evidence row. TCRvdb is proprietary, and it is read once at
the end to validate -- never to tune, never to annotate (`docs/outputs.md` §7).
"""
from __future__ import annotations

from pathlib import Path

import polars as pl

from ..schema import EVIDENCE_TABLE_COLUMNS

#: The support count is per receptor chain against one epitope: ``clonotype_id`` already covers
#: species, gene, CDR3, V and J.
SUPPORT_KEY: tuple[str, ...] = ("clonotype_id", "antigen.epitope")

INDEPENDENT_STUDY = "independent_study"


def support_counts(records: pl.DataFrame, chains: pl.DataFrame) -> pl.DataFrame:
    """Distinct reporting publications per clonotype-epitope pair, and which they are.

    One row per pair, with ``references`` the sorted distinct ``reference.id`` list and ``studies``
    its length. Phase 11 fits the clustering ``coef`` against ``studies > 1``; the evidence rows
    below are the same numbers reshaped, so the two cannot drift.

    Measured on the current corpus, human only: 4,129 of 187,238 clonotype-epitope pairs (2.21 %)
    have >= 2 distinct references when counted on `chunks/` as submitted, and 4,974 of 184,660
    when counted here, after CDR3 repair -- repair merges sequences, so pairs fall and replication
    rises. Both are the same definition at two stages of the pipeline; the tuning objective uses
    this one, because this is what ships.
    """
    return (
        chains.select("record_id", "clonotype_id")
        .join(records.select("record_id", "antigen.epitope", "reference.id"),
              on="record_id", how="left")
        .group_by(SUPPORT_KEY)
        .agg(pl.col("reference.id").unique().sort().alias("references"))
        .with_columns(pl.col("references").list.len().alias("studies"))
        .sort(SUPPORT_KEY)
    )


def independent_study(records: pl.DataFrame, chains: pl.DataFrame,
                      *, release: str = "dev") -> pl.DataFrame:
    """One row per other reporting publication for each replicated chain-epitope observation.

    ``evidence_value`` is one other reference, so a reader gets "who else saw this" rather than a
    list containing the paper they are already reading. ``evidence_score`` is the number of distinct
    references, the pair's total.
    """
    supported = support_counts(records, chains).filter(pl.col("studies") > 1)
    if supported.is_empty():
        return empty()
    return (
        chains.select("record_id", "gene", "clonotype_id")
        .join(records.select("record_id", "antigen.epitope", "reference.id"),
              on="record_id", how="left")
        .join(supported, on=SUPPORT_KEY, how="inner")
        .with_columns(pl.col("references").list.set_difference(pl.concat_list("reference.id"))
                      .list.sort().alias("__reference"))
        .explode("__reference", empty_as_null=True)
        .select(
            "record_id", "gene",
            # `(record_id, evidence_id)` is the key, so the id only has to separate a record's own
            # evidence. A readable natural value beats a hash nobody can trace back.
            pl.concat_str(pl.lit(INDEPENDENT_STUDY), pl.col("gene"), pl.col("__reference"), separator=":")
              .alias("evidence_id"),
            pl.lit(INDEPENDENT_STUDY).alias("evidence_type"),
            # Each source citation is an individual value, never a joined list.
            pl.col("__reference").alias("evidence_source"),
            pl.col("__reference").alias("evidence_value"),
            pl.col("studies").cast(pl.Float64).alias("evidence_score"),
            pl.lit(release).alias("first_seen_release"),
        )
        .sort("record_id", "evidence_id")
    )


def empty() -> pl.DataFrame:
    """An evidence table with no rows but the right schema, so a producer may legitimately find none."""
    return pl.DataFrame(schema={
        "record_id": pl.String, "gene": pl.String, "evidence_id": pl.String,
        "evidence_type": pl.String, "evidence_source": pl.String, "evidence_value": pl.String,
        "evidence_score": pl.Float64, "first_seen_release": pl.String,
    }).select(EVIDENCE_TABLE_COLUMNS)


def build_evidence(records: pl.DataFrame, chains: pl.DataFrame,
                   *, release: str = "dev", structure_metadata: pl.DataFrame | None = None,
                   structure_coordinates: Path | None = None,
                   structure_full_germline: bool = False) -> pl.DataFrame:
    """Every producer's rows, concatenated.

    Publication support and optional structure metadata are assembled here. The motif stage
    appends its freshly inferred rows to the same table.
    """
    if structure_coordinates is not None and structure_metadata is None:
        raise ValueError("structure coordinates require structure metadata")
    if structure_full_germline and structure_coordinates is None:
        raise ValueError("full structure germline requires structure coordinates")
    parts = [independent_study(records, chains, release=release)]
    if structure_metadata is not None:
        parts.append(structures(chains, structure_metadata, release=release, records=records,
                                coordinates=structure_coordinates, full_germline=structure_full_germline))
    return (pl.concat([p for p in parts if not p.is_empty()] or [empty()], how="vertical")
            .select(EVIDENCE_TABLE_COLUMNS)
            .sort("record_id", "evidence_id"))


STRUCTURE_METRICS = ("num_contacts", "ranking_confidence", "plddt", "ptm", "iptm",
                     "tcr_pmhc_iptm", "scanning_angle", "pitch_angle")


def structures(chains: pl.DataFrame, metadata: pl.DataFrame, *, release: str = "dev",
               records: pl.DataFrame | None = None, coordinates: Path | None = None,
               full_germline: bool = False) -> pl.DataFrame:
    """Link historical hashes and paired-junction/peptide model associations.

    Repeated metadata observations of one model must agree. Each predicted metric is a separate
    value, and native metrics are excluded. Absent metrics produce no row. Supplying
    ``records`` enables context-preserving paired-junction/peptide association. Coordinates
    additionally enable complete-domain proof using explicit alleles.
    """
    if full_germline and coordinates is None:
        raise ValueError("full structure germline requires structure coordinates")
    required = {"tcr_pmhc_hash", "is_native"}
    if not required.issubset(metadata.columns):
        raise ValueError(f"structure metadata requires {sorted(required)}")
    m = metadata.with_columns(pl.col("is_native").cast(pl.String).str.to_lowercase())
    if m.filter(~pl.col("is_native").is_in(["true", "false"]) | pl.col("is_native").is_null()).height:
        raise ValueError("is_native must be true or false")
    keys = chains.select(pl.col("TCR_hash").alias("tcr_pmhc_hash"), "record_id").unique()
    keys = keys.filter(pl.col("tcr_pmhc_hash").is_not_null() & (pl.col("tcr_pmhc_hash") != ""))
    if keys.group_by("record_id").len().filter(pl.col("len") > 1).height:
        raise ValueError("one record cannot identify different structural complexes")
    if records is not None:
        for column in ('species', 'mhc.class'):
            if column in metadata.columns:
                context = keys.join(metadata.select('tcr_pmhc_hash', column).unique(),
                    on='tcr_pmhc_hash', how='inner').join(records.select('record_id', column),
                    on='record_id', suffix='_current', validate='m:1')
                if context.filter(pl.col(column).is_null() |
                                  (pl.col(column) != pl.col(column + '_current'))).height:
                    raise ValueError(f'structure hash disagrees on {column}')
    predicted_keys = keys
    native_keys = None
    if records is not None:
        from .structures import associate_links
        predicted_extra, native_extra = associate_links(records, chains, m, keys)
        predicted_keys = pl.concat([keys, predicted_extra]).unique()
        native_keys = native_extra
    if coordinates is not None:
        if records is None:
            raise ValueError('coordinate recovery requires records')
        from .structures import recover_links
        keys = recover_links(records, chains, metadata, coordinates, keys, full_germline=full_germline)
        predicted_keys = pl.concat([keys, predicted_keys]).unique()
    parts = []
    native = m.filter(pl.col("is_native") == "true")
    if native.height:
        if "meta.structure.id" not in native.columns or native.filter(
            ~pl.col("meta.structure.id").fill_null("").str.contains(r"^[0-9][A-Za-z0-9]{3}$")
        ).height:
            raise ValueError("native metadata requires one PDB id per row")
        native_parts = native.select("tcr_pmhc_hash",
            pl.lit("structure_native").alias("evidence_type"),
            pl.lit("PDB").alias("evidence_source"),
            pl.col("meta.structure.id").str.to_lowercase().alias("evidence_value"),
            pl.lit(None, dtype=pl.Float64).alias("evidence_score"))
        native_linked = native_parts.join(keys, on="tcr_pmhc_hash", how="inner")
        if native_keys is not None:
            native_linked = pl.concat([native_linked, native_parts.join(native_keys.rename(
                {"pdb_id": "evidence_value"}), on=["tcr_pmhc_hash", "evidence_value"], how="inner")])
        parts.append(native_linked)
    metrics = [c for c in STRUCTURE_METRICS if c in m.columns]
    if metrics:
        predicted = (m.filter(pl.col("is_native") == "false")
            .select("tcr_pmhc_hash", *metrics)
            .unpivot(index="tcr_pmhc_hash", variable_name="evidence_value", value_name="__value")
            .with_columns(pl.col("__value").cast(pl.String).replace(["", "NA", "NaN", "nan"], None)
                          .cast(pl.Float64).alias("evidence_score"))
            .drop_nulls("evidence_score"))
        if predicted.filter(~pl.col("evidence_score").is_finite()).height:
            raise ValueError("structure metrics must be finite or absent")
        if predicted.group_by("tcr_pmhc_hash", "evidence_value").agg(
            pl.col("evidence_score").n_unique().alias("n")).filter(pl.col("n") > 1).height:
            raise ValueError("conflicting metrics for one structure hash; supply distinct model identities")
        parts.append(predicted.select("tcr_pmhc_hash",
            pl.lit("structure_model").alias("evidence_type"),
            pl.col("tcr_pmhc_hash").alias("evidence_source"), "evidence_value", "evidence_score")
            .join(predicted_keys, on="tcr_pmhc_hash", how="inner"))
    if not parts:
        return empty()
    return (pl.concat(parts).unique()
        .with_columns(pl.lit("").alias("gene"), pl.lit(release).alias("first_seen_release"),
            pl.when(pl.col("evidence_type") == "structure_native").then(
                pl.concat_str("evidence_type", "evidence_source", "evidence_value",
                              "tcr_pmhc_hash", separator=":"))
            .otherwise(pl.concat_str("evidence_type", "evidence_source", "evidence_value", separator=":"))
            .alias("evidence_id"))
        .select(EVIDENCE_TABLE_COLUMNS).unique().sort("record_id", "evidence_id"))


def motifs(records: pl.DataFrame, chains: pl.DataFrame, cohort: pl.DataFrame,
           members: pl.DataFrame, *, method: str, release: str = "dev") -> pl.DataFrame:
    """Producer-scored membership: positive cluster size, zero nonmember, absent unscored.

    TCRNET evidence follows junction identity without V/J restrictions; TCREmp uses its complete
    embedding input. Both are conditioned on species, chain and peptide, not an inferred MHC.
    """
    if method not in ("tcrnet", "tcremp"):
        raise ValueError(f"unknown motif method: {method}")
    key = ["species", "gene", "antigen.epitope", "junction_aa"]
    if method == "tcremp":
        key += ["v_call", "j_call"]
    eligible = cohort.select(key).unique()
    hits = members.select(*key, "cid", "csz").unique()
    if hits.group_by(key).len().filter(pl.col("len") > 1).height:
        raise ValueError("conflicting motif membership for one producer input")
    if hits.join(eligible, on=key, how="anti").height:
        raise ValueError("motif member absent from producer cohort")
    query = (chains.select("record_id", "gene", pl.col("cdr3").alias("junction_aa"),
                           pl.col("v.segm").alias("v_call"), pl.col("j.segm").alias("j_call"))
        .join(records.select("record_id", "species", "antigen.epitope"), on="record_id", validate="m:1"))
    return (query.join(eligible.join(hits, on=key, how="left", validate="1:1"),
                       on=key, how="inner", validate="m:1")
        .select("record_id", "gene",
            pl.concat_str(pl.lit(f"motif_{method}"), "gene", separator=":").alias("evidence_id"),
            pl.when(pl.col("csz").is_not_null()).then(pl.lit(f"motif_{method}"))
              .otherwise(pl.lit(f"motif_{method}_nonmember")).alias("evidence_type"),
            pl.lit(release).alias("evidence_source"),
            pl.col("cid").fill_null("").alias("evidence_value"),
            pl.col("csz").fill_null(0).cast(pl.Float64).alias("evidence_score"),
            pl.lit(release).alias("first_seen_release"))
        .sort("record_id", "evidence_id"))
