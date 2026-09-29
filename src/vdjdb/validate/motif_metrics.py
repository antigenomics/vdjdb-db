"""Every motif metric, for three fixed sources, recorded and gated on every build.

The motif stage is the one place where the release cannot be reproduced byte for byte, so its
acceptance criterion is a set of measurements rather than a digest (``docs/clustering.md``). Those
measurements were being computed inside test assertions and thrown away, which has a specific cost:
when the corpus changes the numbers move, the assertions still pass because they carry slack, and
**nobody learns the new value**. That happened on 2026-09-28 -- three chunk merges moved the legacy
comparator's ``Q`` by 0.0042 and its retention by 0.0059, well inside the slack, and the drift had
to be reconstructed by hand from a log afterwards.

**Three sources are always measured and always compared** (:data:`SOURCES`):

``legacy``
    the last legacy release's ``cluster_members.txt`` -- 2026-06-03, the file every number in
    ``docs/clustering.md`` is stated against.
``latest``
    the latest release's. Today it is the same file, because no release has shipped since; the two
    rows are kept separate anyway, so the day they diverge the report shows it instead of silently
    changing what "the comparator" means.
``current``
    what this build produced, one row per method (``tcrnet``, ``tcremp``).

plus ``trivial``, the do-nothing partition, because it clears three of the four admissibility axes
on TRA and a bar it passes is not a bar (``docs/clustering.md`` 8).

All four are scored **on this build's cohort**, through one code path, so no two rows can be
measured by different code. Two gates follow, and they answer different questions:

* :func:`regressions_against_latest` -- did *this build* get worse than what is published? This is
  the build gate, and a failure is a defect in the code that just changed.
* :func:`regressions_against_baseline` -- did any source drift from its last recorded value? A
  fixed released file scored on a changed cohort moves, so this catches a *corpus* change: rewriting
  a ``v.segm`` cell or adding records gives the released clustering a key it has never seen. That is
  explainable and sometimes expected, which is why :data:`BASELINE` is editable -- but moving it is
  a reviewed diff whose commit message says why.

:data:`BASELINE` is a committed, reviewed input and is never written by a build (hard rule 9).
"""
from __future__ import annotations

from pathlib import Path

import polars as pl

from ..compare import clustering
from . import motif_bench as mb
from . import qscore

#: The committed baseline: one row per source per gated axis, so it is greppable and diffable.
BASELINE = Path("rules/motif_metrics.tsv")

#: The three fixed comparison points, plus the partition that clusters nothing. ``current`` expands
#: to one row per method the build produced.
SOURCES = ("legacy", "latest", "current-tcrnet", "current-tcremp", "trivial")

#: ``axis -> (direction, tolerance)``. ``up`` means a fall is the regression, ``down`` means a rise
#: is, ``record`` means the axis is reported and never gated.
#:
#: ``f1_assignment`` is named in full on purpose. It is ``metrics_lib``'s F1 of the cluster
#: assignment -- TRA 0.86 for the released clustering -- and it is **not** the objective F1 of
#: ``docs/denoising.md`` 6.1, which is ``clustered`` predicting independent replication and reads
#: 0.1086 for the same file. Two different numbers called ``f1`` in one repository is the confusion
#: this report exists to remove.
#:
#: The [0, 1] axes take 0.005, the slack ``tests/release/test_motif_reproduction.py`` already used:
#: the build is deterministic, so this is room for the corpus to grow between a recorded value and a
#: run, not for numerical noise.
#:
#: ``epitopes`` takes 0. Epitope coverage is the fourth admissibility axis (``docs/denoising.md``
#: 7.1) and it counts epitopes that receive any denoising at all, so losing one is a loss to the
#: database rather than rounding.
#:
#: ``percolation_excess`` is this source's percolation median **minus** ``latest``'s, both computed on
#: the epitopes the two of them share. It exists because the own-set median is not comparable across
#: sources that cover different epitopes, and it is a difference rather than a level because the
#: comparable pair depends on which source is being compared -- one ``latest`` row cannot serve all of
#: them, and pairing a source's restricted median against ``latest``'s unrestricted one is the mistake
#: that hid this on the first attempt.
#:
#: Measured on TRA: our TCREMP reads 0.4762 against the release's 0.4583 on their own sets, a gap of
#: 0.0179, but 0.4387 against 0.3765 on the 98 epitopes both cover -- an excess of 0.0622, three and a
#: half times larger. The seven epitopes only TCREMP covers are small and necessarily concentrated
#: (percolation 0.5 to 1.0) and pull its own-set median up; the release's five unique epitopes do the
#: same to its own. So the own-set figure understates the difference. ``latest``'s own excess is 0.0 by
#: construction.
#:
#: ``clusters`` and ``clonotypes`` are recorded and not gated: more clusters is better for coverage
#: and worse for parsimony, so neither direction is a regression on its own, and ``Q`` already
#: charges for the trade.
#:
#: ``partition_*`` is the one family measured **against the shipped file** rather than against the
#: cohort: :func:`vdjdb.compare.clustering.agreement` asks whether this source still gives a record
#: the cluster-mates ``latest`` gave it. Purity and retention are statistics about a clustering;
#: these say whether the table a consumer downloads still groups the records the last one grouped,
#: which is what `vdjdb-web` shows when it lists a record's neighbours.
#:
#: They exist because nothing compared those tables. ``vdjdb diff`` keys ``cluster_members.txt`` on
#: ``cid``, and a cid carries a position in a sorted list, so one renumbered cluster reads as the
#: whole file replaced - which is why the CI comparison names its five members with ``--only`` and
#: skips the motif files entirely.
#:
#: Only ``partition_neighbours_preserved`` is gated, and the reason is measured: 19,971 of the
#: released TRB clustering's 36,906 clonotypes sit in one cluster, which holds 94.7 % of the file's
#: co-clustered pairs, so ``partition_pairs_preserved`` and ``partition_ari`` are measurements of
#: whether that one blob was reproduced - the do-nothing partition scores 0.9991 and 0.9939 on them.
#: They are recorded because they are the exact quantity and because the blob is worth watching.
#:
#: Direction ``baseline``: gated against the committed baseline only, never against ``latest``. These
#: are agreement *with* ``latest``, so ``latest``'s own value is 1.0 by construction and comparing a
#: candidate to it would fail every candidate that is not byte-identical.
AXES: dict[str, tuple[str, float]] = {
    "purity": ("up", 0.005),
    "precision": ("up", 0.005),
    "recall": ("up", 0.005),
    "f1_assignment": ("up", 0.005),
    "retention": ("up", 0.005),
    "q": ("up", 0.005),
    "h": ("up", 0.005),
    "p": ("up", 0.005),
    "epitopes": ("up", 0.0),
    "percolation_median": ("down", 0.005),
    "percolation_excess": ("down", 0.005),
    "clusters": ("record", 0.0),
    "clonotypes": ("record", 0.0),
    "partition_neighbours_preserved": ("baseline", 0.005),
    "partition_pairs_preserved": ("record", 0.0),
    "partition_ari": ("record", 0.0),
    "partition_clonotypes_reference": ("record", 0.0),
    "partition_clonotypes_shared": ("record", 0.0),
}

KEY = ("species", "gene", "source", "axis")


def sources(legacy: pl.DataFrame, latest: pl.DataFrame,
            current: dict[str, pl.DataFrame]) -> dict[str, pl.DataFrame]:
    """The source name -> members frame mapping :func:`measure` scores, in :data:`SOURCES` order.

    ``trivial`` is not here: it is built inside :func:`measure` from each chain's own cohort, because
    ``motif_bench.cohort`` takes one gene at a time and a single partition passed in from outside
    would score one chain against an empty frame -- measured, it read as TRA purity 0.000.
    """
    out = {"legacy": legacy, "latest": latest}
    for method, members in current.items():
        out[f"current-{method}"] = members
    return out


def measure(chains: pl.DataFrame, records: pl.DataFrame,
            members_by_source: dict[str, pl.DataFrame], *, species: str = "HomoSapiens",
            genes: tuple[str, ...] = ("TRA", "TRB")) -> pl.DataFrame:
    """The metric table: one row per ``(species, gene, source, axis)``, one value per cell.

    Long rather than wide, because the gate is per axis and a long frame joins to the baseline on
    :data:`KEY` with no column-name arithmetic.
    """
    rows: list[dict] = []
    for gene in genes:
        cohort = mb.cohort(chains, records, species=species, gene=gene)
        per_gene = {**members_by_source, "trivial": mb.trivial_members(cohort)}
        ref = per_gene.get("latest")
        ref_chain = (ref.filter((pl.col("gene") == gene) & (pl.col("species") == species))
                     if ref is not None else None)
        reference_epitopes: list[str] = []
        reference_pe: pl.DataFrame | None = None
        if ref is not None:
            rpe = mb.per_epitope(cohort, ref.filter((pl.col("gene") == gene)
                                                    & (pl.col("species") == species)))
            reference_pe = rpe.filter(pl.col("clusters") > 0)
            reference_epitopes = sorted(reference_pe["antigen.epitope"])
        for source, members in per_gene.items():
            mm = members.filter((pl.col("gene") == gene) & (pl.col("species") == species))
            s = mb.score(mb.assign(cohort, mm))
            q = qscore.score(cohort, mm)
            pe = mb.per_epitope(cohort, mm)
            covered = pe.filter(pl.col("clusters") > 0)
            # Both sides on the SAME epitopes: the ones this source and `latest` both cover.
            shared = sorted(set(covered["antigen.epitope"]) & set(reference_epitopes))
            mine = covered.filter(pl.col("antigen.epitope").is_in(shared))["percolation"].median()
            theirs = (reference_pe.filter(pl.col("antigen.epitope").is_in(shared))["percolation"]
                      .median()) if reference_pe is not None else None
            vals = {
                "purity": s["purity"], "precision": s["precision"], "recall": s["recall"],
                "f1_assignment": s["f1"], "retention": s["retention"],
                "q": q["q"], "h": q["h"], "p": q["p"],
                "epitopes": float(covered.height),
                "percolation_median": float(covered["percolation"].median() or 0.0),
                "percolation_excess": (0.0 if mine is None or theirs is None
                                       else float(mine) - float(theirs)),
                "clusters": float(q["clusters"]), "clonotypes": float(q["clustered"]),
            }
            # Against the shipped table, not the cohort. `ref_chain` is `latest` restricted to this
            # chain, so a source is compared with the file it would replace.
            agr = (clustering.agreement(ref_chain, mm) if ref_chain is not None
                   else dict.fromkeys(("neighbours_preserved", "pairs_preserved", "ari",
                                       "clonotypes.reference", "clonotypes.shared"), 0.0))
            vals |= {"partition_neighbours_preserved": agr["neighbours_preserved"],
                     "partition_pairs_preserved": agr["pairs_preserved"],
                     "partition_ari": agr["ari"],
                     "partition_clonotypes_reference": agr["clonotypes.reference"],
                     "partition_clonotypes_shared": agr["clonotypes.shared"]}
            for axis, value in vals.items():
                rows.append({"species": species, "gene": gene, "source": source, "axis": axis,
                             "direction": AXES[axis][0], "tolerance": AXES[axis][1],
                             "value": round(float(value), 6)})
    return pl.DataFrame(rows).sort("species", "gene", "source", "axis")


#: The baseline's columns. ``reason`` is what makes a regression against ``latest`` acceptable: a
#: non-empty reason declares the trade, and the declared ``value`` caps it, so the axis fails again
#: the moment it gets worse than what was declared. Same shape as ``rules/expected_diffs.toml``,
#: where a rule carries both a cause and a measured count.
COLUMNS = {"species": pl.Utf8, "gene": pl.Utf8, "source": pl.Utf8, "axis": pl.Utf8,
           "direction": pl.Utf8, "tolerance": pl.Float64, "value": pl.Float64, "reason": pl.Utf8}


def load_baseline(path: Path = BASELINE) -> pl.DataFrame:
    """The committed baseline, or an empty frame shaped like one when it does not exist yet."""
    return pl.read_csv(path, separator="\t", schema=COLUMNS) if path.exists() \
        else pl.DataFrame(schema=COLUMNS)


def _worse(now: pl.Expr, was: pl.Expr, direction: pl.Expr, tol: pl.Expr) -> pl.Expr:
    """``baseline`` behaves like ``up`` here and is excluded from the ``latest`` comparison."""
    return (pl.when(now.is_null()).then(pl.lit(True))
            .when(direction.is_in(["up", "baseline"])).then(now < was - tol)
            .when(direction == "down").then(now > was + tol)
            .otherwise(pl.lit(False)))


def regressions_against_baseline(measured: pl.DataFrame,
                                 baseline: pl.DataFrame) -> pl.DataFrame:
    """Gated axes that moved outside tolerance since the baseline was recorded, any source.

    A source in the baseline but missing from ``measured`` counts as a regression: a clustering that
    stopped being produced scores nothing rather than scoring badly, and a left join would hide it.
    """
    j = baseline.join(measured.select(*KEY, pl.col("value").alias("now")),
                      on=list(KEY), how="left")
    return (j.filter(_worse(pl.col("now"), pl.col("value"), pl.col("direction"),
                            pl.col("tolerance")))
            .with_columns((pl.col("now") - pl.col("value")).alias("delta"))
            .select(*KEY, "direction", "tolerance", pl.col("value").alias("baseline"), "now",
                    "delta")
            .sort("gene", "source", "axis"))


def regressions_against_latest(measured: pl.DataFrame,
                               baseline: pl.DataFrame | None = None) -> pl.DataFrame:
    """Gated axes where a ``current-*`` source is worse than the ``latest`` release, undeclared.

    Both sides come from the same run on the same cohort, so there is no corpus-drift term here: a
    failure is either a defect in the code that just changed, or a trade somebody has to declare.

    A trade is declared by a non-empty ``reason`` on the baseline row for that exact axis, and the
    declared ``value`` caps it -- our TRA TCREMP percolates more than the release and buys retention,
    purity and two epitopes for it, which is a trade; the same axis drifting further is not, and
    fails. Without ``baseline`` nothing is declared and every regression is reported.
    """
    latest = (measured.filter(pl.col("source") == "latest")
              .select("species", "gene", "axis", pl.col("value").alias("latest")))
    # `baseline`-direction axes are agreement *with* `latest`, so its own value is 1.0 by
    # construction and this comparison would fail every candidate that is not byte-identical to it.
    cur = measured.filter(pl.col("source").str.starts_with("current-")
                          & (pl.col("direction") != "baseline"))
    j = cur.join(latest, on=["species", "gene", "axis"], how="inner")
    bad = (j.filter(_worse(pl.col("value"), pl.col("latest"), pl.col("direction"),
                           pl.col("tolerance")))
           .with_columns((pl.col("value") - pl.col("latest")).alias("delta")))
    if baseline is not None and baseline.height:
        dec = (baseline.filter(pl.col("reason").is_not_null() & (pl.col("reason") != ""))
               .select(*KEY, pl.col("value").alias("declared"), "reason"))
        bad = (bad.join(dec, on=list(KEY), how="left")
               # Declared and no worse than the declaration: accepted. Declared and worse: fails,
               # and `declared` in the output says what it was measured at when it was accepted.
               .filter(pl.col("declared").is_null()
                       | _worse(pl.col("value"), pl.col("declared"), pl.col("direction"),
                                pl.col("tolerance"))))
    return (bad.select(*KEY, "direction", "tolerance", "latest", pl.col("value").alias("now"),
                       "delta")
            .sort("gene", "source", "axis"))


def report(measured: pl.DataFrame, baseline: pl.DataFrame) -> str:
    """The table as markdown: every axis, every source, against both the baseline and ``latest``.

    Written into the CI step summary so a reviewer reads the numbers without downloading an
    artifact. That is the point of a report rather than a pass/fail line: the value that did *not*
    trip the gate is the one nobody would otherwise ever see.
    """
    latest = (measured.filter(pl.col("source") == "latest")
              .select("species", "gene", "axis", pl.col("value").alias("vs_latest")))
    j = (measured
         .join(baseline.select(*KEY, pl.col("value").alias("was")), on=list(KEY), how="left")
         .join(latest, on=["species", "gene", "axis"], how="left")
         .with_columns((pl.col("value") - pl.col("was")).alias("d_base"),
                       (pl.col("value") - pl.col("vs_latest")).alias("d_latest")))
    out = ["| gene | source | axis | value | baseline | vs baseline | vs latest | gate |",
           "|---|---|---|--:|--:|--:|--:|---|"]
    for r in j.sort("gene", "source", "axis").iter_rows(named=True):
        was = "-" if r["was"] is None else f"{r['was']:.4f}"
        db = "-" if r["d_base"] is None else f"{r['d_base']:+.4f}"
        dl = "-" if r["d_latest"] is None or r["source"] == "latest" else f"{r['d_latest']:+.4f}"
        gate = "recorded" if r["direction"] == "record" \
            else f"{r['direction']}, tol {r['tolerance']:g}"
        out.append(f"| {r['gene']} | {r['source']} | `{r['axis']}` | {r['value']:.4f} | {was} | "
                   f"{db} | {dl} | {gate} |")
    return "\n".join(out)
