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
from vdjdb.schema import CHAIN_COLUMNS, EVIDENCE_TABLE_COLUMNS, RECORD_COLUMNS

pytestmark = pytest.mark.release

#: One chunk row is one record, so this is the row count the build reads. 192,753 until
#: `PMID_18025130` landed its 40 (#161).
EXPECTED_RECORDS = 192_793


@pytest.fixture(scope="module")
def tables() -> dict[str, pl.DataFrame]:
    d = Path(os.environ.get("VDJDB_TABLES", "out/tables"))
    if not (d / "records.parquet").exists():
        pytest.skip(f"no built tables at {d}; run `vdjdb build --out out/`")
    return {n: pl.read_parquet(d / f"{n}.parquet") for n in ("records", "chains", "evidence")}


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
    """#462 fills a gap; it does not second-guess a curator. 686 of the 711 chains with no V get a
    call, which is material for the phase 9 decision, not a change to what ships today."""
    chains = tables["chains"]
    assert chains.filter((pl.col("v.segm") != "") & (pl.col("v.inferred") != "")).is_empty()
    assert chains.filter((pl.col("j.segm") != "") & (pl.col("j.inferred") != "")).is_empty()
    assert chains.filter((pl.col("v.segm") == "") & (pl.col("v.inferred") != "")).height > 0


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
