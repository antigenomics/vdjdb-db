"""The reference diff between two releases: the database changelog (#432).

Not a code changelog. A reader of a VDJdb release wants to know which studies arrived, how many
records each brought, and what the totals now are.

The publication years are already a committed input: :mod:`vdjdb.summary.references` resolved them
once, so dating the new studies costs a join rather than 600 network calls.

⚠ Nothing here reads TCRvdb.
"""
from __future__ import annotations

import zipfile
from pathlib import Path

import polars as pl

MEMBER = "vdjdb.slim.txt"


def _references(source: Path) -> pl.DataFrame:
    """``(reference.id, rows, content)`` from a release zip or a built legacy directory.

    ``vdjdb.slim.txt`` holds a comma-joined ``reference.id`` per row -- one slim row can come from
    several studies -- so the column is split before counting. Counting the joined strings instead
    would invent a distinct "reference" for every combination that occurs. ``rows`` therefore sums
    to more than the table's row count: it is a reference-row pair count, not a record count.

    ``content`` is a digest of the reference's own ``(cdr3, epitope)`` multiset. It distinguishes a
    re-identification from a study appearing and another disappearing: issue #347 replaced DOI and
    preprint URLs with PubMed ids, and without this the release notes would report three studies
    withdrawn and three added when nothing changed but a name.
    """
    if source.is_dir():
        text = (source / MEMBER).read_text()
    else:
        with zipfile.ZipFile(source) as zf:
            name = next(n for n in zf.namelist() if Path(n).name == MEMBER)
            text = zf.read(name).decode()
    df = pl.read_csv(text.encode(), separator="\t", infer_schema_length=0, quote_char=None)
    return (df.select("cdr3", "antigen.epitope",
                      pl.col("reference.id").str.split(",").alias("ref"))
              # Pinned rather than left to the default, which Polars 2.0 flips. It makes no
              # difference to this data -- `str.split` on an empty cell yields `[""]`, never an
              # empty list, so the `!= ""` filter below catches it either way -- and pinning stops
              # an upgrade from changing a release note with no error.
              .explode("ref", empty_as_null=False)
              .with_columns(pl.col("ref").str.strip_chars().alias("reference.id"))
              .filter(pl.col("reference.id") != "")
              .with_columns((pl.col("cdr3") + "|" + pl.col("antigen.epitope")).alias("__k"))
              .group_by("reference.id").agg(
                  pl.len().alias("rows"),
                  pl.col("__k").sort().str.join(";").hash().alias("content"))
              .sort("rows", descending=True))


def diff(previous: Path, current: Path, *, years: Path | None = None) -> dict:
    """What changed between two releases, as data. ``years`` is the committed reference table."""
    before, after = _references(previous), _references(current)
    added = after.join(before.select("reference.id"), on="reference.id", how="anti")
    removed = before.join(after.select("reference.id"), on="reference.id", how="anti")

    # A reference that vanished and one that appeared with the SAME rows is one study renamed, not
    # two events. Matched on content, never on row count alone -- two unrelated studies can
    # contribute the same number of rows.
    renamed = (removed.join(added, on="content", how="inner", suffix="_new")
               .select(pl.col("reference.id").alias("from"),
                       pl.col("reference.id_new").alias("to"),
                       pl.col("rows")))
    added = added.join(renamed.select(pl.col("to").alias("reference.id")),
                       on="reference.id", how="anti")
    removed = removed.join(renamed.select(pl.col("from").alias("reference.id")),
                           on="reference.id", how="anti")
    if years is not None and years.exists():
        y = pl.read_csv(years, separator="\t")
        added = added.join(y.select("reference.id", "year"), on="reference.id", how="left")
    return {
        "references_before": before.height,
        "references_after": after.height,
        "rows_before": int(before["rows"].sum()),
        "rows_after": int(after["rows"].sum()),
        "added": added.sort("rows", descending=True),
        "removed": removed.sort("rows", descending=True),
        "renamed": renamed.sort("rows", descending=True),
    }


def render(d: dict, *, tag: str = "") -> str:
    """The release-notes body, ready to paste without editing."""
    head = f"## VDJdb {tag}".rstrip() if tag else "## This release"
    dr = d["rows_after"] - d["rows_before"]
    ds = d["references_after"] - d["references_before"]
    out = [head, "",
           f"**{d['references_after']:,} studies**, contributing {d['rows_after']:,} "
           f"study-row pairs to `vdjdb.slim.txt` ({ds:+,} studies, {dr:+,} pairs).", ""]
    if d["added"].height:
        out += [f"### {d['added'].height} new studies", "",
                "| Study | Rows | Year |", "|---|---:|---:|"]
        has_year = "year" in d["added"].columns
        for row in d["added"].iter_rows(named=True):
            year = row.get("year") if has_year else None
            out.append(f"| {row['reference.id']} | {row['rows']:,} | {year or '--'} |")
        out.append("")
    if d["removed"].height:
        # A study leaving a release is a curation decision, not a routine event -- it means a
        # retraction, a re-attribution or a withdrawn chunk. Never summarise it as a count alone.
        out += [f"### {d['removed'].height} studies no longer present", "",
                "| Study | Rows it carried |", "|---|---:|"]
        for row in d["removed"].iter_rows(named=True):
            out.append(f"| {row['reference.id']} | {row['rows']:,} |")
        out.append("")
    if d["renamed"].height:
        out += [f"### {d['renamed'].height} studies re-identified", "",
                "Same rows, new identifier -- a nomenclature fix, not a data change.", "",
                "| Was | Now | Rows |", "|---|---|---:|"]
        for row in d["renamed"].iter_rows(named=True):
            out.append(f"| {row['from']} | {row['to']} | {row['rows']:,} |")
        out.append("")
    return "\n".join(out)
