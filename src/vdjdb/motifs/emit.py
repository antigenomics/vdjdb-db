"""Project the motif tables into the two files ``vdjdb-web`` parses.

⚠ **Both files are parsed positionally.** ``app/backend/server/motifs/Motifs.scala`` hands Tablesaw a
fixed ``Array[ColumnType]`` with no header check -- 27 entries for ``motif_pwms.txt``, 19 for
``cluster_members.txt``. An inserted, removed or reordered column silently mistypes or shifts the
whole table rather than failing. :data:`PWM_COLUMNS` and :data:`MEMBER_COLUMNS` are that contract,
and the order in them is the only order either file may ever be written in (CLAUDE.md hard rule 1).

Everything this module adds is a join and a rename. The motif tables carry the cluster; the record
columns (``antigen.gene``, the MHC) and the per-chain geometry (``v.end``, ``j.start``) come back
from the definitive tables, exactly as the legacy exporter is a projection of the new build rather
than a second pipeline (ROADMAP section 6).
"""
from __future__ import annotations

from pathlib import Path

import polars as pl

#: ``cluster_members.txt``, 19 columns, in the order ``Motifs.scala`` types them.
MEMBER_COLUMNS: tuple[str, ...] = (
    "species", "antigen.epitope", "antigen.gene", "antigen.species", "mhc.a", "mhc.b", "mhc.class",
    "gene", "cdr3aa", "x", "y", "cid", "csz", "v.segm", "j.segm", "v.end", "j.start",
    "v.segm.repr", "j.segm.repr",
)

#: ``motif_pwms.txt``, 27 columns, likewise.
PWM_COLUMNS: tuple[str, ...] = (
    "species", "antigen.epitope", "gene", "aa", "pos", "len", "v.segm.repr", "j.segm.repr",
    "cid", "csz", "count", "count.bg", "total.bg", "count.bg.i", "total.bg.i", "need.impute",
    "freq", "freq.bg", "I", "I.norm", "height.I", "height.I.norm",
    "antigen.gene", "antigen.species", "mhc.a", "mhc.b", "mhc.class",
)

#: The record columns both files carry, resolved once per ``(species, epitope)``.
_ANNOTATION = ("antigen.gene", "antigen.species", "mhc.a", "mhc.b", "mhc.class")


def _annotation(records: pl.DataFrame) -> pl.DataFrame:
    """One annotation row per ``(species, antigen.epitope)``: the modal value of each column.

    An epitope can appear under more than one MHC across publications -- HLA promiscuity is real and
    is issue #372 -- but the legacy files carry one value per cluster, so the modal one is what a
    positional reader can be given. Ties break lexicographically, never on frame order (hard rule 7).
    """
    out = records.select("species", "antigen.epitope", *_ANNOTATION)
    for col in _ANNOTATION:
        modal = (out.group_by(["species", "antigen.epitope", col]).len()
                    .sort(["len", col], descending=[True, False])
                    .group_by(["species", "antigen.epitope"], maintain_order=True).first()
                    .select("species", "antigen.epitope", col))
        out = out.drop(col).unique(subset=["species", "antigen.epitope"], maintain_order=True) \
                 .join(modal, on=["species", "antigen.epitope"], how="left")
    return out


def cluster_members(members: pl.DataFrame, chains: pl.DataFrame,
                    records: pl.DataFrame) -> pl.DataFrame:
    """``cluster_members.txt`` as a frame, 19 columns in :data:`MEMBER_COLUMNS` order.

    ``v.end`` and ``j.start`` are 0-based amino-acid indices in junction space -- the same space
    ``cdr3aa`` is in, because VDJdb's ``cdr3`` includes both anchors. They come from the chain that
    reported the clonotype; where several did and they disagree, the first by sorted record order.
    """
    geometry = (
        chains.join(records.select("record_id", "species", "antigen.epitope"), on="record_id")
        .select("species", "antigen.epitope", "gene",
                pl.col("cdr3").alias("cdr3aa"), pl.col("v.segm"), pl.col("j.segm"),
                pl.col("v.end"), pl.col("j.start"), pl.col("record_id"))
        .sort("record_id")
        .unique(subset=["species", "antigen.epitope", "gene", "cdr3aa", "v.segm", "j.segm"],
                keep="first", maintain_order=True)
        .drop("record_id")
    )
    return (
        members
        .select("species", "antigen.epitope", "gene",
                pl.col("junction_aa").alias("cdr3aa"),
                pl.col("v_call").alias("v.segm"), pl.col("j_call").alias("j.segm"),
                "x", "y", "cid", "csz", "v.segm.repr", "j.segm.repr")
        .join(geometry, on=["species", "antigen.epitope", "gene", "cdr3aa", "v.segm", "j.segm"],
              how="left")
        .join(_annotation(records), on=["species", "antigen.epitope"], how="left")
        .with_columns(pl.col("v.end", "j.start").fill_null(-1))
        .select(*MEMBER_COLUMNS)
        .sort("species", "gene", "antigen.epitope", "cid", "cdr3aa")
    )


def motif_pwms(pwms: pl.DataFrame, records: pl.DataFrame) -> pl.DataFrame:
    """``motif_pwms.txt`` as a frame, 27 columns in :data:`PWM_COLUMNS` order.

    ``need.impute`` becomes a real flag: ``TRUE`` wherever the ``(v, j, len)`` stratum had nothing
    at that position and the column was read off the coarser background instead. In the shipped file
    it is ``FALSE`` on all 13,456 rows, because it was computed after the filter that would have set
    it (ROADMAP section 8.5).
    """
    return (
        pwms
        .with_columns(
            pl.when(pl.col("level.bg") == "vj_len").then(pl.lit("FALSE")).otherwise(pl.lit("TRUE"))
            .alias("need.impute"))
        .join(_annotation(records), on=["species", "antigen.epitope"], how="left")
        .select(*PWM_COLUMNS)
        .sort("species", "gene", "antigen.epitope", "cid", "len", "pos", "aa")
    )


def write(members: pl.DataFrame, pwms: pl.DataFrame, out: Path, *,
          suffix: str = "") -> dict[str, int]:
    """Write both files under ``out``. Returns ``{filename: rows}``.

    ``suffix`` names the method's variant: TCREMP ships alongside TCRNET as
    ``cluster_members_tcremp.txt`` / ``motif_pwms_tcremp.txt``, which ``vdjdb-web`` reads from the
    same directory and parses with the same fixed column types.

    Tab-separated, never quoted, ``\\n`` terminated: neither file contains a ``"`` and the Scala
    reader would take one literally.
    """
    out.mkdir(parents=True, exist_ok=True)
    written = {}
    for name, df, cols in ((f"cluster_members{suffix}.txt", members, MEMBER_COLUMNS),
                           (f"motif_pwms{suffix}.txt", pwms, PWM_COLUMNS)):
        assert tuple(df.columns) == cols, f"{name}: column order is a contract, got {df.columns}"
        df.write_csv(out / name, separator="\t", quote_style="never", line_terminator="\n")
        written[name] = df.height
    return written
