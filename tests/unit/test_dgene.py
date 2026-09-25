"""D-segment confidence: the posterior belongs to the call we ship, not to arda's own winner."""
from __future__ import annotations

import polars as pl
import pytest

from vdjdb.annotate import dgene

CHAINS = pl.DataFrame({
    "record_id": ["r1", "r2", "r3"],
    "gene": ["TRB", "TRA", "TRB"],
    "cdr3": ["CASSIRSSYEQYF", "CAVSDLEPNSSASKIIF", "CASSPGQGAYEQYF"],
    "v.segm": ["TRBV10-3*01", "TRAV12-2*01", "TRBV5-1*01"],
    "j.segm": ["TRBJ2-7*01", "TRAJ3*01", "TRBJ2-7*01"],
    "d.inferred": ["TRBD2*01", "", "TRBD1*01"],
})
RECORDS = pl.DataFrame({"record_id": ["r1", "r2", "r3"],
                        "species": ["HomoSapiens"] * 3})


@pytest.fixture(scope="module")
def scored() -> pl.DataFrame:
    return dgene.add_d_posterior(CHAINS, RECORDS)


def test_a_chain_with_no_inferred_d_gets_no_posterior(scored):
    """TRA has no D at all; a null here is the quantity not existing, not a failed lookup."""
    tra = scored.filter(pl.col("gene") == "TRA")
    assert tra["d.posterior"].null_count() == tra.height
    assert tra["d.entropy"].null_count() == tra.height


def test_a_beta_chain_gets_a_posterior_and_an_entropy(scored):
    trb = scored.filter(pl.col("gene") == "TRB")
    assert trb["d.posterior"].null_count() == 0
    assert all(0.0 <= p <= 1.0 for p in trb["d.posterior"])
    assert all(e >= 0.0 for e in trb["d.entropy"])


def test_the_posterior_is_for_the_gene_we_shipped(scored):
    """arda names its own winner; we report the probability of *our* call, so the two disagreeing
    shows up as a low number rather than as a confident wrong one."""
    from arda.dpost import posterior_d

    from vdjdb.annotate.cdr3fix import ensure_reference
    ensure_reference()

    row = scored.filter(pl.col("record_id") == "r1").row(0, named=True)
    p = posterior_d("CASSIRSSYEQYF", "TRBV10-3*01", "TRBJ2-7*01", "human")
    assert row["d.posterior"] == pytest.approx(p.by_gene["TRBD2"])
    assert row["d.entropy"] == pytest.approx(p.entropy)


def test_no_row_is_dropped_and_the_order_is_the_key(scored):
    assert scored.height == CHAINS.height
    assert scored["record_id"].to_list() == sorted(CHAINS["record_id"])
    assert "species" not in scored.columns


def test_an_unsupported_species_gets_nulls_without_calling_arda():
    macaque = pl.DataFrame({"record_id": ["r1"], "species": ["MacacaMulatta"]})
    got = dgene.add_d_posterior(CHAINS.head(1), macaque)
    assert got["d.posterior"].null_count() == 1
