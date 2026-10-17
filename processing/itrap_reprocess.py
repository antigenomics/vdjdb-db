# 2026-10-16
# Reprocess deposited ITRAP count tables with the dCODE chain-graph caller.
# The negative-binomial model is not refit or retuned. These loaders only
# reshape deposited GEM or 10x matrices into the frames that caller expects.
# ITRAP1 raw.csv has no empty droplets and no negative-control dextramer.
# The viral matrix uses the deposited barcode-to-peptide key. Unmapped
# antibody features are not given a peptide name.

from __future__ import annotations

import csv
import gzip
import importlib.util
import math
from pathlib import Path

import numpy as np
import polars as pl

_SCRIPT = Path(__file__).resolve().parent / "dcode_chain_graph.py"
_spec = importlib.util.spec_from_file_location("dcode_chain_graph", _SCRIPT)
assert _spec is not None and _spec.loader is not None
dcode = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dcode)

AA = set("ACDEFGHIKLMNPQRSTVWY")
CONTIG_SCHEMA = {
    "barcode": pl.String, "is_cell": pl.String, "chain": pl.String,
    "v_gene": pl.String, "j_gene": pl.String, "productive": pl.String,
    "cdr3": pl.String, "umis": pl.Float64, "donor": pl.Int64,
}


def _parts(value) -> list[str]:
    if value is None:
        return []
    text = str(value).strip()
    if text.lower() in {"", "na", "nan", "none", "null"}:
        return []
    return text.split("|")


def _vj(field: str) -> tuple[str, str]:
    v_gene = ""
    j_gene = ""
    for token in field.split(";"):
        name = token.strip().split("*", 1)[0]
        if name.startswith("TRAV") or name.startswith("TRBV"):
            v_gene = name
        elif name.startswith("TRAJ") or name.startswith("TRBJ"):
            j_gene = name
    return v_gene, j_gene


def _donor_id(value) -> int | None:
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return None


def itrap1_tables(frame: pl.DataFrame) -> tuple[pl.DataFrame, pl.DataFrame, dict]:
    """GEM lists to a contig table and a barcode by peptide count matrix.

    A chain whose CDR3 list and UMI list differ in length is skipped. The
    returned donor column is the cell-hash sample id.
    """
    contig_rows: list[dict] = []
    count_rows: list[dict] = []
    peptides: set[str] = set()
    skipped_chains = 0
    for gem, sample, genes_a, cdr_a, umi_a, genes_b, cdr_b, umi_b, pep, umi_m in frame.select(
        "gem", "sample_id", "genes_lst_TRA", "cdr3_lst_TRA", "umi_count_lst_TRA",
        "genes_lst_TRB", "cdr3_lst_TRB", "umi_count_lst_TRB",
        "peptide_HLA_lst", "umi_count_lst_mhc",
    ).iter_rows():
        donor = _donor_id(sample)
        if donor is None:
            continue
        for locus, genes, cdrs, umis in (
            ("TRA", _parts(genes_a), _parts(cdr_a), _parts(umi_a)),
            ("TRB", _parts(genes_b), _parts(cdr_b), _parts(umi_b)),
        ):
            if len(cdrs) != len(umis) or (genes and len(genes) != len(cdrs)):
                skipped_chains += 1
                continue
            for index, cdr3 in enumerate(cdrs):
                sequence = cdr3.strip().upper()
                if len(sequence) < 8 or any(char not in AA for char in sequence):
                    continue
                v_gene, j_gene = _vj(genes[index] if index < len(genes) else "")
                try:
                    umi = float(umis[index])
                except ValueError:
                    continue
                contig_rows.append({
                    "barcode": gem,
                    "is_cell": "True",
                    "chain": locus,
                    "v_gene": v_gene,
                    "j_gene": j_gene,
                    "productive": "True",
                    "cdr3": sequence,
                    "umis": umi,
                    "donor": donor,
                })
        pep_parts = _parts(pep)
        umi_parts = _parts(umi_m)
        if len(pep_parts) != len(umi_parts):
            continue
        counts: dict[str, float] = {}
        for token, raw_umi in zip(pep_parts, umi_parts, strict=True):
            pieces = token.split()
            if len(pieces) != 2:
                continue
            peptide, allele = pieces[0].upper(), pieces[1].upper()
            if not dcode.ALLELE_RE.match(allele):
                continue
            if not peptide or any(char not in AA for char in peptide):
                continue
            column = f"{allele}_{peptide}"
            try:
                counts[column] = counts.get(column, 0.0) + float(raw_umi)
            except ValueError:
                continue
            peptides.add(column)
        counts["barcode"] = gem
        counts["donor"] = donor
        count_rows.append(counts)
    contigs = pl.DataFrame(contig_rows, schema=CONTIG_SCHEMA)
    peptide_cols = sorted(peptides)
    columns = ["barcode", "donor", *peptide_cols]
    filled = []
    for row in count_rows:
        full = {col: 0.0 for col in peptide_cols}
        full["barcode"] = row["barcode"]
        full["donor"] = row["donor"]
        for col in peptide_cols:
            if col in row:
                full[col] = row[col]
        filled.append(full)
    schema = {"barcode": pl.String, "donor": pl.Int64, **{col: pl.Float64 for col in peptide_cols}}
    features = pl.DataFrame(filled, schema=schema) if filled else pl.DataFrame(schema=schema)
    features = features.select(columns)
    stats = {
        "gems": frame.height,
        "contigs": contigs.height,
        "skipped_chain_lists": skipped_chains,
        "peptides": len(peptides),
        "donors": contigs["donor"].n_unique() if contigs.height else 0,
    }
    return contigs, features, stats


def run_frames(contigs: pl.DataFrame, features: pl.DataFrame) -> dict[str, pl.DataFrame]:
    """One donor at a time, through the committed chain graph and specificity model."""
    keys = ("chains", "edges", "clones", "calls", "checks")
    parts: dict[str, list[pl.DataFrame]] = {key: [] for key in keys}
    if contigs.is_empty():
        return {key: pl.DataFrame() for key in keys}
    for donor in contigs["donor"].unique().sort().to_list():
        donor_contigs = contigs.filter(pl.col("donor") == donor).drop("donor")
        barcodes = donor_contigs["barcode"].unique().to_list()
        donor_features = features.filter(
            (pl.col("donor") == donor) & pl.col("barcode").is_in(barcodes)
        ).drop("donor")
        chains = dcode.chain_table(donor_contigs, int(donor))
        edges = dcode.edge_table(chains)
        endogenous = dcode.mark_endogenous(chains)
        clones = dcode.clone_table(endogenous)
        dex_cols = dcode.dextramer_columns(donor_features.columns)
        nc_cols = dcode.negative_control_columns(donor_features.columns)
        clones = dcode.drop_doublets(clones, donor_features, dex_cols + nc_cols)
        calls, checks = dcode.assign_specificities(
            clones, donor_features, dex_cols, donor_features, nc_cols
        )
        parts["chains"].append(endogenous)
        parts["edges"].append(edges)
        parts["clones"].append(clones)
        parts["calls"].append(calls)
        parts["checks"].append(checks)
    out = {}
    for key, frames in parts.items():
        frames = [frame for frame in frames if frame.height]
        out[key] = pl.concat(frames, how="diagonal_relaxed") if frames else pl.DataFrame()
        if out[key].height:
            out[key] = out[key].sort(dcode.TABLE_SORT_COLUMNS[key])
    return out


def _gene_name(value) -> str:
    text = "" if value is None else str(value).strip()
    if text.lower() in {"", "na", "nan", "none"}:
        return ""
    return text.split("*", 1)[0]


def _matrix_antibody(
    matrix_path: Path,
    barcodes_path: Path,
    features_path: Path,
    peptide_of_feature: dict[str, str],
) -> pl.DataFrame:
    """10x matrix market file to a barcode by mapped-peptide count table."""
    with gzip.open(features_path, "rt") as handle:
        feature_names = [line.split("\t")[0] for line in handle if line.strip()]
    wanted: dict[int, str] = {}
    for index, name in enumerate(feature_names, start=1):
        key = name.replace("_", "-")
        column = peptide_of_feature.get(key)
        if column:
            wanted[index] = column
    with gzip.open(barcodes_path, "rt") as handle:
        barcodes = [line.strip() for line in handle if line.strip()]
    columns = sorted(set(wanted.values()))
    column_index = {name: i for i, name in enumerate(columns)}
    matrix = np.zeros((len(barcodes), len(columns)), dtype=np.float64)
    with gzip.open(matrix_path, "rt") as handle:
        expected_entries = None
        entries = 0
        for line in handle:
            line = line.strip()
            if not line or line.startswith("%"):
                continue
            fields = line.split()
            if len(fields) != 3:
                raise ValueError("MatrixMarket rows must have three fields")
            if expected_entries is None:
                n_features, n_barcodes, expected_entries = map(int, fields)
                if (n_features, n_barcodes) != (len(feature_names), len(barcodes)) or expected_entries < 0:
                    raise ValueError("MatrixMarket dimensions do not match the feature and barcode files")
                continue
            feature_i, barcode_i = map(int, fields[:2])
            value = float(fields[2])
            if not (1 <= feature_i <= len(feature_names) and 1 <= barcode_i <= len(barcodes)):
                raise ValueError("MatrixMarket coordinate is outside the declared dimensions")
            if not math.isfinite(value) or value < 0.0 or not value.is_integer():
                raise ValueError("MatrixMarket UMI counts must be finite nonnegative integers")
            entries += 1
            column = wanted.get(feature_i)
            if column is None:
                continue
            matrix[barcode_i - 1, column_index[column]] += value
        if expected_entries is None or entries != expected_entries:
            raise ValueError("MatrixMarket entry count does not match its declaration")
    frame = pl.DataFrame({"barcode": barcodes, **{name: matrix[:, i] for i, name in enumerate(columns)}})
    return frame


def load_peptide_key(path: Path) -> dict[str, str]:
    frame = pl.read_csv(path, separator="\t", infer_schema_length=0)
    out = {}
    for barcode, sequence, hla in frame.select("Barcode", "Sequence", "HLA").iter_rows():
        allele = str(hla).upper()
        peptide = str(sequence).upper()
        if dcode.ALLELE_RE.match(allele) and peptide and all(char in AA for char in peptide):
            out[str(barcode)] = f"{allele}_{peptide}"
    return out


def viral_tables(
    alpha_path: Path,
    beta_path: Path,
    meta_path: Path,
    matrix_path: Path,
    barcodes_path: Path,
    features_path: Path,
    peptide_key_path: Path,
) -> tuple[pl.DataFrame, pl.DataFrame, dict]:
    """Deposited viral contigs and antibody counts, restricted to hashed donors."""
    peptide_of_feature = load_peptide_key(peptide_key_path)
    features = _matrix_antibody(matrix_path, barcodes_path, features_path, peptide_of_feature)
    meta = pl.read_csv(meta_path, separator="\t", infer_schema_length=0)
    name_col = meta.columns[0]
    donor_of: dict[str, str] = {}
    for cell, donor in meta.select(name_col, "BC").iter_rows():
        if donor and str(donor).lower() not in {"na", "nan", "none"}:
            donor_of[str(cell)] = str(donor)
    labels = {label: index + 1 for index, label in enumerate(sorted(set(donor_of.values())))}
    contig_rows = []
    seen_cells = set()
    dropped_unhashed = 0
    skipped_ragged = 0
    for path, locus in ((alpha_path, "TRA"), (beta_path, "TRB")):
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            width = len(reader.fieldnames or [])
            for raw in reader:
                if raw is None or len(raw) != width or None in raw:
                    skipped_ragged += 1
                    continue
                cell = raw["cell_id"]
                productive = raw["productive"]
                v_call = raw["v_call"]
                j_call = raw["j_call"]
                junction = raw["junction_aa"]
                umi = raw["umi_count"]
                if str(productive).upper() not in {"T", "TRUE"}:
                    continue
                sequence = "" if junction is None else str(junction).strip().upper()
                if len(sequence) < 8 or any(char not in AA for char in sequence):
                    continue
                bare = str(cell).split("-", 1)[0]
                key = f"1OS_{bare}"
                label = donor_of.get(key)
                if label is None:
                    dropped_unhashed += 1
                    continue
                seen_cells.add(str(cell))
                try:
                    umi_value = float(umi)
                except (TypeError, ValueError):
                    continue
                contig_rows.append({
                    "barcode": str(cell),
                    "is_cell": "True",
                    "chain": locus,
                    "v_gene": _gene_name(v_call),
                    "j_gene": _gene_name(j_call),
                    "productive": "True",
                    "cdr3": sequence,
                    "umis": umi_value,
                    "donor": labels[label],
                })
    contigs = pl.DataFrame(contig_rows, schema=CONTIG_SCHEMA)
    keep = set(contigs["barcode"].unique().to_list()) if contigs.height else set()
    features = features.filter(pl.col("barcode").is_in(list(keep))).with_columns(pl.lit(0).alias("donor"))
    # Donor on a feature row follows its contig donor. A barcode is one donor.
    donor_frame = contigs.select("barcode", "donor").unique(subset=["barcode"])
    features = features.drop("donor").join(donor_frame, on="barcode", how="inner")
    stats = {
        "mapped_peptides": len(set(peptide_of_feature.values())),
        "matrix_barcodes": features.height,
        "contigs": contigs.height,
        "donors": contigs["donor"].n_unique() if contigs.height else 0,
        "dropped_unhashed_contigs": dropped_unhashed,
        "skipped_ragged_rows": skipped_ragged,
        "cells": len(seen_cells),
    }
    return contigs, features, stats


def write_result(result: dict[str, pl.DataFrame], out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    for key, frame in result.items():
        frame.write_csv(out / f"{key}.tsv", separator="\t")
