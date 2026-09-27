"""The three stages of a segment call stay distinguishable, and the engine's call is the shipped one.

`chains` reports `v.segm.submitted` (what the publication said), `v.segm` (what ships) and
`v.segm.arda` (what the markup engine aligned against). Collapsing any two of them loses the only
thing that tells a curation decision from a markup one.
"""
from __future__ import annotations

import polars as pl
import pytest

from vdjdb.schema import CHAIN_COLUMNS

pytest.importorskip("arda", reason="the markup engine")


@pytest.fixture(scope="module")
def marked() -> pl.DataFrame:
    from vdjdb.annotate.cdr3fix import markup

    keys = pl.DataFrame({
        "species": ["HomoSapiens"] * 3 + ["Bos taurus"],
        # a bare gene the engine resolves to an allele; a specific gene; a blank J; a species
        # with no reference at all
        "cdr3": ["CASSIRSSYEQYF", "CASSIGPALNTEAFF", "CASRGQGFSYEQYF", "CASSQQF"],
        "v": ["TRBV10-3", "TRBV27*01", "TRBV6-5", "TRBV1"],
        "j": ["TRBJ2-7", "TRBJ1-1", "", "TRBJ1-1"],
    })
    return markup(keys, "beta")


def test_the_engines_own_call_is_reported_beside_the_shipped_one(marked):
    row = marked.filter(pl.col("cdr3") == "CASSIRSSYEQYF").row(0, named=True)
    # The submitted call names no allele; the engine names the one it aligned against, and that is
    # what ships -- every release to date does the same, on 99.2 % of `v.segm` cells.
    assert row["v"] == "TRBV10-3"
    assert row["__varda"] == "TRBV10-3*01"
    assert row["__v"] == "TRBV10-3*01"


def test_a_multi_call_uses_the_comma_the_legacy_format_expects(marked):
    """arda separates an ambiguous call with `;`. Every VDJdb column uses `,` (README STYLE)."""
    from vdjdb.annotate.cdr3fix import markup

    keys = pl.DataFrame({"species": ["HomoSapiens"], "cdr3": ["CASSLGGQGAFYNEQFF"],
                         "v": ["TRBV12-3,TRBV12-4"], "j": ["TRBJ2-1"]})
    out = markup(keys, "beta")
    for col in ("__v", "__varda", "__j", "__jarda"):
        assert ";" not in out[col][0], f"{col} leaked arda's `;` separator: {out[col][0]!r}"


def test_a_species_with_no_reference_still_produces_every_column(marked):
    """An unmappable species is marked unmapped, not dropped, and its schema matches the rest."""
    row = marked.filter(pl.col("species") == "Bos taurus").row(0, named=True)
    assert row["__varda"] == "" and row["__jarda"] == ""
    assert row["__vend"] == -1 and row["__jstart"] == -1
    # the call the curator gave survives: there is no evidence against it, only no reference for it
    assert row["__v"] == "TRBV1"


def test_the_chain_table_declares_all_three_stages():
    for col in ("v.segm.submitted", "j.segm.submitted", "d.segm.submitted",
                "v.segm.arda", "j.segm.arda", "v.segm", "j.segm", "d.segm"):
        assert col in CHAIN_COLUMNS, f"{col} is not a declared chain column"
