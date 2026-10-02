"""The reference-year resolver: the parsers that need no network, and the committed table itself.

The network resolvers (PubMed, RCSB, GitHub) are exercised by ``vdjdb refs``, not by the suite --
a test that fails when NCBI is slow is a test that gets disabled. What is asserted here is the part
that has actually gone wrong: identifier parsing, and whether the committed table still covers the
database it is supposed to describe.
"""
from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from vdjdb.emit.vdjdb3 import read_table
from vdjdb.summary import references as R


def test_offline_resolvers_cover_arxiv_biorxiv_and_the_two_literals():
    df = R.resolve([
        "https://arxiv.org/abs/2503.00648",
        "https://doi.org/10.1101/2020.05.04.20085779",
        "http://mediatum.ub.tum.de/doc/1136748",
        "",           # the 854-record curation gap: dropped, never counted as a lookup failure
        "   ",
    ])
    assert dict(zip(df["reference.id"], df["year"], strict=True)) == {
        "https://arxiv.org/abs/2503.00648": 2025,
        "https://doi.org/10.1101/2020.05.04.20085779": 2020,
        "http://mediatum.ub.tum.de/doc/1136748": 2014,
    }


def test_pmid_pattern_tolerates_the_stray_space_but_keeps_the_literal_id():
    # 22 records carry `PMID: 34433824`; the same paper is also present without the space. Both
    # must resolve, and both must key on what the database literally says.
    assert R._PMID.match("PMID: 34433824")[1] == "34433824"
    assert R._PMID.match("PMID:34433824")[1] == "34433824"
    assert R._PMID.match("PMID:not-a-number") is None


def test_patent_publication_years_resolve_offline():
    expected = {"US20220324939A1": 2022, "US20230060095A1": 2023,
                "WO2017048593A1": 2017, "WO2024163935A2": 2024}
    prefix = "https://patents.google.com/patent/"
    result = R.resolve([prefix + name for name in expected] + [prefix + "US12345678B2"])
    assert dict(zip(result["reference.id"], result["year"], strict=True)) == {
        prefix + name: year for name, year in expected.items()}
    assert result["source"].unique().to_list() == ["patent-publication-id"]


def test_pdb_and_issue_patterns():
    assert R._PDB.match("https://www.rcsb.org/structure/9WBD")[1] == "9WBD"
    assert R._PDB.match("https://www.rcsb.org/structure/9WBD/")[1] == "9WBD"
    assert R._ISSUE.match("https://github.com/antigenomics/vdjdb-db/issues/397").groups() == (
        "antigenomics", "vdjdb-db", "397")


def test_unresolved_reports_records_not_references():
    records = pl.DataFrame({"reference.id": ["PMID:1", "PMID:1", "PMID:2", ""]})
    table = pl.DataFrame({"reference.id": ["PMID:1"], "year": [2020], "source": ["pubmed"]})
    out = R.unresolved(records, table)
    assert out.to_dicts() == [{"reference.id": "PMID:2", "records": 1}]


@pytest.mark.skipif(not Path("out/tables/records.parquet").exists(),
                    reason="needs a build; run `uv run vdjdb build --out out/`")
def test_committed_table_still_covers_every_reference_in_the_database():
    """The table is only useful if it has not gone stale -- which is exactly how it failed before."""
    records = read_table(Path("out/tables"), "records")
    missing = R.unresolved(records, R.load())
    assert missing.is_empty(), (
        f"{missing.height} references have no year, covering "
        f"{int(missing['records'].sum()):,} records. Run `uv run vdjdb refs`.\n{missing}")
