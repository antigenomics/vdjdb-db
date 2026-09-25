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

#: One chunk row is one record, so this is the row count the build reads.
EXPECTED_RECORDS = 192_753


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


def test_empty_string_is_the_only_missing_marker(tables):
    """CLAUDE.md rule 6. A null in a shipped table is the pandas three-way ambiguity coming back."""
    for name, frame in tables.items():
        nulls = {c: n for c, n in frame.null_count().row(0, named=True).items() if n}
        assert not nulls, f"{name}: {nulls}"
