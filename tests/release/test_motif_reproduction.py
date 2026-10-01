"""Is the rewritten motif stage equivalent to or better than the one that shipped?

Marked ``release``: needs ``VDJDB_REFERENCE_ZIP``, built tables (``VDJDB_TABLES``, default
``out/tables``) and built motif files (``VDJDB_MOTIFS``, default ``out/motifs``).

The motif stage is the one place where byte reproduction is neither possible nor wanted: TCRNET is
re-implemented on `vdjtools` with a BH `q_value` the legacy never had, and TCREMP did not exist in
the released pipeline at all. So the acceptance criterion is a comparison rather than a digest, and
this module is that criterion:

* the released clustering, scored through our harness, reproduces the bar `docs/clustering.md`
  publishes. That is what makes every other number here comparable to the documented ones;
* our clustering is not worse than the released one on any scored axis;
* TCRNET, the method that was ported rather than replaced, still clusters the clonotypes the
  released file clustered;
* the structure `vdjdb-web` assumes -- one epitope, one V, one J, one CDR3 length per cluster --
  holds in both files.

The scorecards that chose the configuration are `docs/tuning/scorecard.tsv`, pinned to the prose by
`tests/unit/test_tuning_docs.py`. Nothing here re-runs a sweep.
"""
from __future__ import annotations

import os
from pathlib import Path

import polars as pl
import pytest

from vdjdb.compare.diff import Bundle
from vdjdb.emit.vdjdb3 import read_table
from vdjdb.validate import motif_bench as mb
from vdjdb.validate import motif_metrics as mmv

pytestmark = pytest.mark.release

GENES = ("TRA", "TRB")
SPECIES = "HomoSapiens"

#: The released ``cluster_members.txt`` of 2026-06-03. Fixed: it is a published file.
RELEASE_MEMBER_ROWS = 55_636
RELEASE_MEMBER_CIDS = 1_928
RELEASE_PWM_CIDS = 1_791

#: The bar `docs/clustering.md` sections 4.1 and 4.2 publish for the released clustering. Asserted to
#: two decimals only: the exact current values live in `rules/motif_metrics.tsv`, which the build
#: gates against per axis, and duplicating them here is how the pair went stale -- these read 0.8658
#: and 0.8567 on TRA until the 2026-09-28 chunk merges moved them to 0.8641 and 0.8548, and the
#: `abs=0.01` tolerance meant nothing noticed.
PUBLISHED_BAR = {"TRA": {"purity": 0.8658, "precision": 0.8567},
                 "TRB": {"purity": 0.9790, "precision": 0.9756}}

#: Fraction of the clonotypes the released TCRNET clustered that ours clusters too. Measured
#: 2026-09-27: TRA 12,177 of 13,327 (0.914), TRB 35,867 of 36,416 (0.985). The floor is what a
#: reimplementation has to clear to be called the same method; the TRB figure is the interesting one,
#: because TRB is where the released clustering is strongest.
OVERLAP_FLOOR = {"TRA": 0.85, "TRB": 0.95}

#: Which axes must not fall, in which direction, and by how much is now one table:
#: `vdjdb.validate.motif_metrics.AXES`. It includes `retention` on purpose -- a method can buy purity
#: by clustering almost nothing, and that is the trade the released TCRNET made (0.2105 on TRA).

#: Slack on the published-bar comparison only; the per-axis gate carries its own tolerance.
SLACK = 0.005


@pytest.fixture(scope="module")
def reference() -> Bundle:
    env = os.environ.get("VDJDB_REFERENCE_ZIP")
    if not env or not Path(env).exists():
        pytest.skip("set VDJDB_REFERENCE_ZIP to a release zip or extracted directory")
    return Bundle(Path(env))


@pytest.fixture(scope="module")
def tables() -> tuple[pl.DataFrame, pl.DataFrame]:
    d = Path(os.environ.get("VDJDB_TABLES", "out/tables"))
    if not (d / "chains.parquet").exists():
        pytest.skip(f"no built tables at {d}; run `vdjdb build --out out/`")
    return read_table(d, "chains"), read_table(d, "records")


@pytest.fixture(scope="module")
def motifs() -> Path:
    """The motif files, refused if they predate the tables they would be scored against.

    Both directories default independently, so a fresh `vdjdb build --out out/` beside a motif
    directory from an earlier tuning scores one build's clustering on another build's cohort. That
    happened on 2026-09-27: motif files from 2026-09-25, 71,138 member rows against the current
    53,505, read as a single TCREMP TRB metric falling below its bar. Nothing was wrong with either
    artifact, and nothing in the assertions could say so.

    Modification time rather than a recorded provenance, because there is none to read: these are
    build outputs in an output directory, never a checkout, and in CI both come from one run. It
    catches the stale pair and says which; it cannot prove a matched pair came from one build.
    """
    d = Path(os.environ.get("VDJDB_MOTIFS", "out/motifs"))
    members = d / "cluster_members.txt"
    if not members.exists():
        pytest.skip(f"no built motif files at {d}; run `vdjdb motifs --out {d}`")
    chains = Path(os.environ.get("VDJDB_TABLES", "out/tables")) / "chains.parquet"
    if chains.exists() and members.stat().st_mtime < chains.stat().st_mtime:
        pytest.fail(f"{members} is older than {chains}: the clustering and the cohort are from two "
                    f"builds. Re-run `vdjdb motifs --tables {chains.parent} --out {d}`.")
    return d


@pytest.fixture(scope="module")
def released(reference: Bundle, tmp_path_factory) -> pl.DataFrame:
    p = tmp_path_factory.mktemp("released") / "cluster_members.txt"
    p.write_bytes(reference.read_bytes("cluster_members.txt"))
    return mb.read_members(p)


@pytest.fixture(scope="module")
def ours(motifs: Path) -> dict[str, pl.DataFrame]:
    return {"tcrnet": mb.read_members(motifs / "cluster_members.txt"),
            "tcremp": mb.read_members(motifs / "cluster_members_tcremp.txt")}


@pytest.fixture(scope="module")
def cohorts(tables) -> dict[str, pl.DataFrame]:
    chains, records = tables
    return {g: mb.cohort(chains, records, species=SPECIES, gene=g) for g in GENES}


def _chain(members: pl.DataFrame, gene: str) -> pl.DataFrame:
    return members.filter((pl.col("gene") == gene) & (pl.col("species") == SPECIES))


def _clustered(cohort: pl.DataFrame, members: pl.DataFrame) -> pl.DataFrame:
    """The distinct clonotypes of the cohort this clustering put in a cluster."""
    return (mb.assign(cohort, members).filter(pl.col("cluster").is_not_null())
            .select(*mb.KEY).unique())


# --------------------------------------------------------------------------------------------
# The released files, as published
# --------------------------------------------------------------------------------------------

def test_the_released_clustering_is_the_one_we_measured_against(released: pl.DataFrame) -> None:
    assert released.height == RELEASE_MEMBER_ROWS
    assert released["cid"].n_unique() == RELEASE_MEMBER_CIDS


def test_the_release_ships_clusters_with_no_motif(reference: Bundle,
                                                  released: pl.DataFrame) -> None:
    """A join of the two released files on ``cid`` leaves 137 clusters without a PWM.

    The legacy PWM builder dropped a whole cluster when a residue at some position was absent from
    the background, rather than giving it a pseudocount, so 1,879 member rows (1,584 TRA, 295 TRB)
    reach `vdjdb-web` with no motif to draw. Ours imputes and records it in ``need.impute``
    (`test_motifs.py::test_a_residue_the_background_never_saw_still_gets_a_letter`), which is why
    our two files agree on their cluster set and the released pair does not.

    Asserted so the release's own shape stays on the record: this is a property of the file we
    compare against, not a target to reproduce.
    """
    pwms = pl.read_csv(reference.read_bytes("motif_pwms.txt"), separator="\t",
                       infer_schema_length=0, quote_char=None)
    members, motifs = set(released["cid"].unique()), set(pwms["cid"].unique())
    assert len(motifs) == RELEASE_PWM_CIDS
    assert not motifs - members, "a motif for a cluster with no members"
    assert len(members - motifs) == RELEASE_MEMBER_CIDS - RELEASE_PWM_CIDS


@pytest.mark.parametrize("gene", GENES)
def test_the_harness_reproduces_the_published_legacy_bar(gene: str, released: pl.DataFrame,
                                                         cohorts) -> None:
    """Every comparison below is against these numbers, so they are checked first."""
    got = mb.score(mb.assign(cohorts[gene], _chain(released, gene)))
    for axis, want in PUBLISHED_BAR[gene].items():
        assert got[axis] == pytest.approx(want, abs=0.01), (
            f"{gene} {axis}: {got[axis]:.4f}, docs/clustering.md publishes {want}")


@pytest.fixture(scope="module")
def measured(tables, ours, released) -> pl.DataFrame:
    """The metric table: read from the build's own report when there is a fresh one, else computed.

    `vdjdb motif-metrics` already writes this on every build and already runs both gates, so
    recomputing it here cost 24 s of a 186 s suite to assert what the CLI had just exited 1 on
    (ROADMAP_local section 57.1). One fixture, shared by both tests, reusing the artifact.

    **Reused only when it is newer than both inputs.** A report from an earlier build scored on this
    build's cohort is the exact mistake the `motifs` fixture above documents happening on 2026-09-27,
    and it reads as a metric change rather than as a stale file.
    """
    chains, records = tables
    report = Path(os.environ.get("VDJDB_MOTIF_METRICS", "out/reports/motif-metrics.tsv"))
    newer_than = [Path(os.environ.get("VDJDB_TABLES", "out/tables")) / "chains.parquet",
                  Path(os.environ.get("VDJDB_MOTIFS", "out/motifs")) / "cluster_members.txt"]
    if report.exists() and all(
            p.exists() and report.stat().st_mtime >= p.stat().st_mtime for p in newer_than):
        return pl.read_csv(report, separator="\t")
    return mmv.measure(chains, records, mmv.sources(released, released, ours))


def test_no_source_drifted_from_the_committed_baseline(measured) -> None:
    """The recorded numbers, per axis, for all four sources -- including the released files' own.

    A fixed released file scored on a changed cohort moves, so this is the gate that catches a
    *corpus* change rather than a code one. It is the check that was missing on 2026-09-28.
    """
    drift = mmv.regressions_against_baseline(measured, mmv.load_baseline())
    assert drift.height == 0, (
        f"{drift.height} axes drifted; explain each and re-record with "
        f"`vdjdb motif-metrics --record`:\n{drift}")


# --------------------------------------------------------------------------------------------
# Ours against the release
# --------------------------------------------------------------------------------------------

def test_our_clustering_is_not_worse_than_the_latest_release_on_any_gated_axis(
        measured) -> None:
    """Both chains, both methods, every gated axis, through the same code the build's gate uses.

    An axis that is worse on purpose is declared with its reason in `rules/motif_metrics.tsv` and
    capped at the declared value, so the trade stays visible and cannot quietly get worse. Percolation
    is not among the gated axes: it depends on how many prominent motifs the epitope has.
    """
    bad = mmv.regressions_against_latest(measured, mmv.load_baseline())
    assert bad.height == 0, f"{bad.height} undeclared regressions against the release:\n{bad}"


@pytest.mark.parametrize("gene", GENES)
def test_tcrnet_still_clusters_what_the_released_tcrnet_clustered(gene: str, ours, released,
                                                                  cohorts) -> None:
    """TCRNET was ported, not replaced, so coverage is the property that says it is the same method.

    The clonotypes it no longer clusters are the price of the BH `q_value` the legacy never applied;
    the ones it clusters and the release did not are the enrichment the legacy p threshold missed.
    Both are reported, and only the shortfall is gated.
    """
    cohort = cohorts[gene]
    legacy = _clustered(cohort, _chain(released, gene))
    mine = _clustered(cohort, _chain(ours["tcrnet"], gene))
    both = legacy.join(mine, on=mb.KEY, how="inner").height
    assert legacy.height, "the released clustering covers none of this cohort"
    overlap = both / legacy.height
    assert overlap >= OVERLAP_FLOOR[gene], (
        f"{gene}: {both} of {legacy.height} released clonotypes still clustered "
        f"({overlap:.3f}), floor {OVERLAP_FLOOR[gene]}")


# --------------------------------------------------------------------------------------------
# The structure vdjdb-web assumes
# --------------------------------------------------------------------------------------------

@pytest.mark.parametrize("source", ["released", "tcrnet", "tcremp"])
def test_a_cluster_is_pinned_to_one_epitope_one_v_one_j_and_one_length(source, ours,
                                                                       released) -> None:
    """The legacy ``motif_pwms.txt`` schema carries a single ``len`` and V/J per cluster, so a
    cluster spanning two of any of them cannot be projected into it without being split
    (ROADMAP section 2, the mixed-V problem). Measured on all three files: none does.

    Purity also depends on it. `motif_bench.KEY` deliberately excludes the epitope, so a cid that
    spanned two epitopes would be scored as an impure cluster rather than as a schema violation.
    """
    members = released if source == "released" else ours[source]
    spread = members.group_by("cid").agg(
        pl.col("antigen.epitope").n_unique().alias("epitopes"),
        pl.col("v.segm.repr").n_unique().alias("v"),
        pl.col("j.segm.repr").n_unique().alias("j"),
        pl.col("gene").n_unique().alias("genes"),
        pl.col("species").n_unique().alias("species"),
        pl.col("cdr3aa").str.len_chars().n_unique().alias("lengths"))
    for axis in ("epitopes", "v", "j", "genes", "species", "lengths"):
        bad = spread.filter(pl.col(axis) > 1)
        assert bad.height == 0, f"{source}: {bad.height} cids span several {axis}: " \
                                f"{bad['cid'].head(3).to_list()}"


@pytest.mark.parametrize("source", ["released", "tcrnet", "tcremp"])
def test_every_cid_names_its_species_chain_and_epitope(source, ours, released) -> None:
    """`<species initial>.<chain initial>.<epitope>.<n>`, which is what makes a cid readable in
    `vdjdb-web` and unique across chains without a composite key."""
    members = released if source == "released" else ours[source]
    parts = members.select(pl.col("cid").str.split(".")).to_series().to_list()
    assert all(len(p) == 4 for p in parts[:1000])
    tagged = members.with_columns(pl.col("cid").str.split(".").alias("p"))
    assert tagged.filter(pl.col("p").list.get(1) != pl.col("gene").str.slice(-1)).height == 0
    assert tagged.filter(pl.col("p").list.get(2) != pl.col("antigen.epitope")).height == 0
