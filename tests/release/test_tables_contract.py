"""The definitive tables' contracts, on a real build.

Marked ``release``: needs a built directory (``VDJDB_TABLES``, default ``out/tables``), because the
properties worth asserting here -- a key is unique, a hash does not collide, a foreign key resolves
-- are invisible on the handful of synthetic rows the unit tests use.
"""
from __future__ import annotations

import os
from pathlib import Path

import polars as pl
import pytest

from vdjdb.assemble.tables import CLONOTYPE_KEY
from vdjdb.emit.vdjdb3 import read_table
from vdjdb.schema import CHAIN_COLUMNS, EVIDENCE_TABLE_COLUMNS, RECORD_COLUMNS

pytestmark = pytest.mark.release

#: One chunk row is one record after within-chunk deduplication, so this is the row count the build
#: reads. 192,753 until `PMID_18025130` landed its 40 (#161), then 192,793 until the #646 anchor
#: substitutions.
#:
#: **A repair can lower it, and lowering it is the repair working.** A base-call error at a conserved
#: anchor makes one clone look like two: correcting it makes 33 rows duplicate a sibling inside their
#: own chunk, and `CHUNK_DEDUP_KEY` carries the CDR3, so deduplication collapses them. 192,793 ->
#: 192,763. The retired identifiers keep their amendment trail in `registry/records.tsv` and the
#: surviving record carries the same observation.
#:
#: 192,763 -> 192,641 for a second repair of the same shape, and a larger one. The B16 chunk (#397)
#: recorded a clonotype's count by repeating its row, and a spreadsheet autofill had incremented the
#: gene `Eef2` down the column into `Eef2`..`Eef188` - one label per row, all of them in
#: `CHUNK_DEDUP_KEY` - so 187 rows that are 65 clonotypes never deduplicated. Undoing the autofill
#: removes **122 records that were never real**. The count is now `count/sample total` in
#: `method.frequency`, which is the convention 3,548 values in the corpus already use (#696).
#:
#: 192,641 -> 192,623 for #390, and this one is not a repair of a chunk: `CHUNK_DEDUP_KEY` contains
#: `reference.id`, so a group of it spanning two chunk files is one publication reporting one clone
#: twice, and two rows of one paper are not two independent reports. 19 such groups exist over 38
#: rows; **18 merge** - the paper's own chunk is the base and the other fills its blanks - and 1 is
#: left alone, because `PDB_Database.tsv` and `PMID_34433824.tsv` give one clone `structural` and
#: `tetramer-sort`, which is a solved complex and the sort that found it. The chunk was a proxy for
#: the publication; where the two come apart, the publication is what deduplication is about.
#:
#: 192,623 -> 192,609 for #625: `PMID_39286976.tsv` reported seven paired TCRs as 18 single chains
#: per epitope, and pairing them from the paper's supplementary figure 10 takes 46 rows to 32.
EXPECTED_RECORDS = 192_609


@pytest.fixture(scope="module")
def tables() -> dict[str, pl.DataFrame]:
    d = Path(os.environ.get("VDJDB_TABLES", "out/tables"))
    if not (d / "records.parquet").exists():
        pytest.skip(f"no built tables at {d}; run `vdjdb build --out out/`")
    return {n: read_table(d, n) for n in ("records", "chains", "evidence")}


def test_column_orders_are_the_declared_ones(tables):
    assert tuple(tables["records"].columns) == RECORD_COLUMNS
    assert tuple(tables["chains"].columns) == CHAIN_COLUMNS
    assert tuple(tables["evidence"].columns) == EVIDENCE_TABLE_COLUMNS


def test_record_id_is_a_key_and_there_is_one_per_curated_line(tables):
    records = tables["records"]
    assert records.height == EXPECTED_RECORDS
    assert records["record_id"].n_unique() == records.height
    # chunk.file + chunk.row is the other name for the same line
    assert records.select("chunk.file", "chunk.row").n_unique() == records.height


def test_chains_are_keyed_on_record_and_gene_and_every_one_has_a_record(tables):
    chains, records = tables["chains"], tables["records"]
    assert chains.select("record_id", "gene").n_unique() == chains.height
    assert set(chains["gene"].unique()) == {"TRA", "TRB"}
    assert chains.join(records.select("record_id"), on="record_id", how="anti").is_empty()


def test_clonotype_id_neither_collides_nor_splits(tables):
    """A collision would silently merge two receptors; a split would lose the replication signal."""
    keyed = tables["chains"].join(tables["records"].select("record_id", "species"),
                                 on="record_id", how="left")
    assert keyed.select(pl.struct(CLONOTYPE_KEY)).n_unique() == keyed["clonotype_id"].n_unique()


def test_evidence_is_keyed_and_resolves_to_a_chain(tables):
    ev, chains = tables["evidence"], tables["chains"]
    assert ev.select("record_id", "evidence_id").n_unique() == ev.height
    assert ev.join(chains.select("record_id", "gene"), on=["record_id", "gene"],
                   how="anti").is_empty()
    assert (ev["evidence_score"] >= 2).all(), "independent support means at least two studies"


def test_the_d_geometry_indexes_the_nucleotide_sequence_beside_it(tables):
    """`d.start`/`d.end` come from the same scenario as `cdr3nt`, which is the reason they are
    trustworthy at all: a second model's coordinates would point at a different sequence."""
    d = tables["chains"].filter(pl.col("d.inferred") != "")
    assert d.filter(pl.col("d.end") > pl.col("cdr3nt").str.len_chars()).is_empty()
    assert d.filter(pl.col("d.start") >= pl.col("d.end")).is_empty()
    assert d.filter(pl.col("cdr3nt") == "").is_empty()


def test_only_beta_chains_have_a_d(tables):
    assert tables["chains"].filter((pl.col("gene") == "TRA")
                                   & (pl.col("d.inferred") != "")).is_empty()


def test_the_d_posterior_is_a_probability_and_is_often_low(tables):
    """Not a decoration: on this corpus 34.6 % of beta chains fall below 0.6, which is what a short,
    heavily trimmed segment with two candidates actually looks like."""
    d = tables["chains"].filter(pl.col("d.inferred") != "").drop_nulls("d.posterior")
    assert d.filter((pl.col("d.posterior") < 0) | (pl.col("d.posterior") > 1)).is_empty()
    assert (d["d.posterior"] < 0.6).sum() > 0


def test_an_inferred_segment_never_sits_beside_a_curated_one(tables):
    """#462 fills a gap; it does not second-guess a curator.

    Keyed on the **submitted** call, which is the one a curator wrote. The shipped `j.segm` is not the
    right test any more: the J proposal reaches it, so `j.inferred` and `j.segm` are the same value on
    3,272 chains by design (#658). `v.segm` stays blank on all 745 of its own, so there the two tests
    coincide - which is exactly the asymmetry `vdjdb.annotate.cdr3fix.markup` measured and states.
    """
    chains = tables["chains"]
    for side in ("v", "j"):
        beside = chains.filter((pl.col(f"{side}.segm.submitted") != "")
                               & (pl.col(f"{side}.inferred") != ""))
        shown = beside.select("record_id", "gene", f"{side}.segm.submitted",
                              f"{side}.inferred").head(5)
        assert beside.is_empty(), (
            f"{beside.height} chains carry a proposed {side} beside the one the publication "
            f"reported:\n{shown}")
    assert chains.filter((pl.col("v.segm") == "") & (pl.col("v.inferred") != "")).height > 0
    # The V proposal reaches no shipped column, so `v.segm` is blank wherever `v.inferred` is filled.
    assert chains.filter((pl.col("v.segm") != "") & (pl.col("v.inferred") != "")).is_empty()


def test_no_string_column_is_ever_null(tables):
    """ROADMAP.md rule 6: empty string is the only missing marker.

    Scoped to string columns, which is where the rule bites -- the pandas `None`/`NaN`/`""` three-way
    ambiguity that shipped more than one bug. A *numeric* column has no empty string to use, and NaN
    would be worse than null because it compares unequal to itself and propagates silently. So
    `cdr3nt.pgen` and `cdr3nt.margin` are null for the 24,950 chains with no inferred nucleotide
    junction, and that is the honest representation.
    """
    for name, frame in tables.items():
        nulls = {c: n for c, n in frame.null_count().row(0, named=True).items()
                 if n and frame.schema[c] == pl.String}
        assert not nulls, f"{name}: {nulls}"


def test_a_numeric_column_is_null_only_where_the_quantity_does_not_exist(tables):
    chains = tables["chains"]
    absent = chains.filter(pl.col("cdr3nt") == "")
    assert absent["cdr3nt.pgen"].null_count() == absent.height
    assert chains.filter(pl.col("cdr3nt") != "")["cdr3nt.pgen"].null_count() == 0


# ---------------------------------------------------------------------------------------------
# The inferred junction nucleotides encode the junction they were inferred from
# ---------------------------------------------------------------------------------------------

def test_every_inferred_cdr3nt_back_translates_to_its_own_junction(tables):
    """#461's acceptance criterion, on the corpus rather than on four fixture rows.

    It **is** guaranteed by construction, and this is what holds it to that. `infer_nt_batch`
    enumerates `(V, delV) x (J, delJ) x (D, delD, position)` and picks the best codon assignment
    *within* each scenario, so a scenario that cannot spell the given residues has probability zero
    and is never a candidate. Anything the model cannot encode comes back null rather than wrong:
    probed on human TRB, a stop codon, an `X`, a `Z`, a one- or two-residue junction, an empty string
    and a true CDR3 with its anchors stripped are all declined.

    The one input that would read as a mismatch is a **lower-case** junction, where the nucleotides
    are right and the comparison is case-sensitive. `vdjdb qc` rejects a residue outside the 20
    upper-case letters and zero chains in the built corpus carry one, so that is closed upstream
    rather than tolerated here.

    `vdjtools._core.translate_junctions` is the entry point for exactly this column - it is
    `to_unified_cdr3aa(translate(nt))`, the treatment a *junction* gets, not a generic translate. That
    is why it is the right instrument here and not merely the fast one: an out-of-frame junction is
    translated inward from both ends with the untranslatable middle collapsed to `_`, so it reads as a
    loud mismatch below, where a plain translate would silently drop a trailing partial codon.

    It also does the whole column in one threaded native call: measured on 263,437 sequences, **0.023 s
    against 0.284 s** for `vdjtools.model.translate` in a Python loop, 12.3x, and identical on every
    row. Rule 4 and the reach order both - one batched call into existing C++, never a per-row call and
    never a codon table written here.

    Measured 2026-09-29: 263,437 of 285,989 chains carry an inferred `cdr3nt`, **0 mismatches and 0
    whose length is not exactly three times the junction's**. The unit test covers four rows, which
    cannot see a residue class the corpus has and a fixture does not.
    """
    from vdjtools._core import translate_junctions

    got = tables["chains"].filter(pl.col("cdr3nt") != "").select("cdr3", "cdr3nt")
    assert got.height > 250_000, \
        f"only {got.height:,} chains carry an inferred cdr3nt; the stage produced almost nothing"

    ragged = got.filter(pl.col("cdr3nt").str.len_chars() != 3 * pl.col("cdr3").str.len_chars())
    assert ragged.height == 0, \
        f"{ragged.height} inferred sequences are not three nucleotides per residue:\n{ragged.head(5)}"

    wrong = (got.with_columns(pl.Series("back", translate_junctions(got["cdr3nt"].to_list())))
             .filter(pl.col("back") != pl.col("cdr3")))
    assert wrong.height == 0, (
        f"{wrong.height} inferred sequences do not encode their own junction:\n"
        f"{wrong.head(5)}")


# ---------------------------------------------------------------------------------------------
# The model V/J boundary is a fallback, never an override (#631)
# ---------------------------------------------------------------------------------------------

#: Chains the fallback fills, measured 2026-09-29 after the arda 2.31 / vdjtools 4.7 bump. The markup
#: engine declines `v.end` on 4,163 and `j.start` on 1,164; a germline alignment answers on these.
#: #631 stated 4,060 and 347, and the run before the bump was 3,484 and 346.
#:
#: **Both moves are gains, which is the case the slack below exists for.** A cell the alignment answers
#: is strictly better than one a fallback proposes, so the fallback shrinking is the outcome to want.
#:
#: `v.end.inferred` 3,484 -> 2,442. Decomposed on the same pair of builds, one row per
#: `(record_id, gene)`: 1,127 cells left because arda 2.31 now places the boundary itself on alleles
#: whose IMGT germline record is truncated (`antigenomics/arda#135`), so the fallback correctly stands
#: down; **0** left because the new source declined; 85 arrived. Of the 2,357 cells both sources
#: filled, they agree on 2,116.
#:
#: `j.start.inferred` 346 -> 488, all 142 of them arrivals, none lost.
#:
#: The source also changed, which is why the second move is up: the boundary now comes from
#: `vdjtools.model.germline_boundary` rather than from the argmax recombination history
#: (`antigenomics/vdjtools#182`). Against nucleotide truth the history is 89.08 % exact on `v.end` and
#: the germline alignment 92.89 %; 92.02 % against 97.96 % on `j.start`.
FALLBACK_FILLED = {"v.end.inferred": 2_442, "j.start.inferred": 488}

#: How far below the recorded count a run may sit before it reads as a regression rather than a
#: corpus change. Loose on purpose: a curator naming a V that arda can then align is a *gain*, and it
#: takes a row out of this population.
FALLBACK_SLACK = 200


@pytest.mark.parametrize(("shipped", "fallback"),
                         [("v.end", "v.end.inferred"), ("j.start", "j.start.inferred")])
def test_the_model_boundary_never_lands_where_the_alignment_answered(tables, shipped, fallback):
    """The safety property, on the corpus. `v.end` and `j.start` are what `vdjdb-web` reads out of
    the `cdr3fix` JSON and what the legacy tables carry, so a model boundary reaching one of those
    cells is a data change belonging to a curation decision, not to a build improvement.

    Separate columns rather than a filled one is what makes this hold by construction: the fallback
    reads -1 wherever the alignment answered, so a consumer coalescing the two cannot overwrite an
    alignment answer even by accident.
    """
    chains = tables["chains"]
    over = chains.filter((pl.col(shipped) >= 0) & (pl.col(fallback) >= 0))
    assert over.height == 0, (
        f"{over.height} chains carry a model {fallback} where the alignment already answered "
        f"{shipped}:\n{over.select('record_id', 'gene', shipped, fallback).head(5)}")


@pytest.mark.parametrize("fallback", sorted(FALLBACK_FILLED))
def test_the_fallback_still_fills_the_cells_it_was_measured_on(tables, fallback):
    """A fallback that quietly stopped filling anything is indistinguishable from one that works."""
    got = tables["chains"].filter(pl.col(fallback) >= 0).height
    want = FALLBACK_FILLED[fallback]
    assert got >= want - FALLBACK_SLACK, (
        f"{fallback} fills {got:,} chains against {want:,} recorded; the fallback has stopped "
        f"reaching cells it used to")


def test_a_filled_boundary_is_inside_the_junction_it_describes(tables):
    """A residue index outside its own junction is a coordinate-space error, which is the failure
    mode a conversion between four spaces has (`CLAUDE.md`). -1 is the declared missing marker.
    """
    chains = tables["chains"].with_columns(pl.col("cdr3").str.len_chars().alias("n"))
    for col in FALLBACK_FILLED:
        bad = chains.filter((pl.col(col) >= 0)
                            & ((pl.col(col) > pl.col("n")) | (pl.col("n") == 0)))
        assert bad.height == 0, \
            f"{bad.height} chains put {col} outside their junction:\n{bad.head(5)}"
        assert chains[col].null_count() == 0, f"{col} must use -1, never null (rule 6)"
