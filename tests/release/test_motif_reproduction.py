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
from vdjdb.validate import motif_bench as mb

pytestmark = pytest.mark.release

GENES = ("TRA", "TRB")
SPECIES = "HomoSapiens"

#: The released ``cluster_members.txt`` of 2026-06-03. Fixed: it is a published file.
RELEASE_MEMBER_ROWS = 55_636
RELEASE_MEMBER_CIDS = 1_928
RELEASE_PWM_CIDS = 1_791

#: The bar `docs/clustering.md` sections 4.1 and 4.2 publish for the released clustering, measured
#: through this harness on 2026-09-27: TRA purity 0.8658 / precision 0.8567, TRB 0.9790 / 0.9756.
#: Asserted to two decimals because the cohort is the current corpus and grows with every chunk,
#: while a point of purity is a real change in the instrument.
LEGACY_BAR = {"TRA": {"purity": 0.8658, "precision": 0.8567},
              "TRB": {"purity": 0.9790, "precision": 0.9756}}

#: Fraction of the clonotypes the released TCRNET clustered that ours clusters too. Measured
#: 2026-09-27: TRA 12,177 of 13,327 (0.914), TRB 35,867 of 36,416 (0.985). The floor is what a
#: reimplementation has to clear to be called the same method; the TRB figure is the interesting one,
#: because TRB is where the released clustering is strongest.
OVERLAP_FLOOR = {"TRA": 0.85, "TRB": 0.95}

#: Axes where our clustering must not fall below the released one. `retention` is included on
#: purpose: a method can buy purity by clustering almost nothing, and that trade is what the
#: released TCRNET made (0.2105 on TRA).
NOT_WORSE = ("purity", "precision", "recall", "f1", "retention")

#: Slack on those comparisons. The build is deterministic, so this is room for the corpus to grow
#: between the measurement and the run, not for numerical noise.
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
    return pl.read_parquet(d / "chains.parquet"), pl.read_parquet(d / "records.parquet")


@pytest.fixture(scope="module")
def motifs() -> Path:
    d = Path(os.environ.get("VDJDB_MOTIFS", "out/motifs"))
    if not (d / "cluster_members.txt").exists():
        pytest.skip(f"no built motif files at {d}; run `vdjdb motifs --out {d}`")
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
    for axis, want in LEGACY_BAR[gene].items():
        assert got[axis] == pytest.approx(want, abs=0.01), (
            f"{gene} {axis}: {got[axis]:.4f}, docs/clustering.md publishes {want}")


# --------------------------------------------------------------------------------------------
# Ours against the release
# --------------------------------------------------------------------------------------------

@pytest.mark.parametrize("method", ["tcrnet", "tcremp"])
@pytest.mark.parametrize("gene", GENES)
def test_our_clustering_is_not_worse_on_any_scored_axis(method: str, gene: str,
                                                        ours, released, cohorts) -> None:
    legacy = mb.score(mb.assign(cohorts[gene], _chain(released, gene)))
    got = mb.score(mb.assign(cohorts[gene], _chain(ours[method], gene)))
    worse = {a: (round(got[a], 4), round(legacy[a], 4))
             for a in NOT_WORSE if got[a] < legacy[a] - SLACK}
    assert not worse, f"{method} {gene} below the released clustering on {worse} (got, released)"


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
