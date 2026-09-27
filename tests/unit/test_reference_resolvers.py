"""The publication-year resolvers, with the network stubbed.

`summary/references.py` was 57 % covered: the offline resolvers had tests, the three that call an
API did not, so nothing checked that a PubMed, RCSB or GitHub payload is read correctly. The
dashboard stops rather than plotting an incomplete year axis, so a parse that quietly returns
nothing turns into a failed render with no explanation.

Every response here is a recorded shape, not a live call. Nothing in this file reaches the network.
"""
from __future__ import annotations

import json

import pytest

from vdjdb.summary import references as r


@pytest.fixture
def no_network(monkeypatch):
    """Fail loudly if anything tries to actually fetch."""
    def boom(url, **kw):
        raise AssertionError(f"unexpected network call: {url}")
    monkeypatch.setattr(r, "_get", boom)


def stub_get(monkeypatch, payload, record=None):
    def fake(url, **kw):
        if record is not None:
            record.append(url)
        return json.dumps(payload).encode()
    monkeypatch.setattr(r, "_get", fake)


def test_pubmed_reads_the_leading_year_of_every_pubdate_form(monkeypatch):
    stub_get(monkeypatch, {"result": {
        "uids": ["1", "2", "3", "4"],
        "1": {"pubdate": "2019 Dec 12"},
        "2": {"pubdate": "2020"},
        "3": {"pubdate": "2021 Jan-Feb"},
        "4": {"pubdate": ""},              # unresolvable, and must not become a guess
    }})
    assert r._pubmed(["1", "2", "3", "4"]) == {"1": 2019, "2": 2020, "3": 2021}


def test_pubmed_batches_rather_than_calling_per_id(monkeypatch):
    """One request per PMID_CHUNK ids. A call per id is the mistake hard rule 3 forbids."""
    calls = []
    stub_get(monkeypatch, {"result": {"uids": []}}, record=calls)
    r._pubmed([str(i) for i in range(r.PMID_CHUNK * 2 + 1)])
    assert len(calls) == 3, f"{len(calls)} requests for {r.PMID_CHUNK * 2 + 1} ids"


def test_pubmed_skips_the_uids_key_and_non_dict_entries(monkeypatch):
    stub_get(monkeypatch, {"result": {"uids": ["1"], "1": "not a record"}})
    assert r._pubmed(["1"]) == {}


def test_rcsb_reads_the_initial_release_date(monkeypatch):
    stub_get(monkeypatch, {"data": {"entries": [
        {"rcsb_id": "1ao7", "rcsb_accession_info": {"initial_release_date": "1997-06-16T00:00:00Z"}},
        {"rcsb_id": "5D2N", "rcsb_accession_info": {"initial_release_date": "2016-01-20T00:00:00Z"}},
    ]}})
    assert r._rcsb(["1AO7", "5D2N"]) == {"1AO7": 1997, "5D2N": 2016}


def test_rcsb_makes_one_call_for_every_id(monkeypatch):
    calls = []
    stub_get(monkeypatch, {"data": {"entries": []}}, record=calls)
    r._rcsb([f"{i}XYZ" for i in range(50)])
    assert len(calls) == 1


def test_rcsb_tolerates_a_missing_accession_block(monkeypatch):
    stub_get(monkeypatch, {"data": {"entries": [{"rcsb_id": "1AO7", "rcsb_accession_info": None}]}})
    assert r._rcsb(["1AO7"]) == {}


def test_rcsb_with_no_ids_does_not_call_at_all(no_network):
    assert r._rcsb([]) == {}


def test_issues_with_no_refs_does_not_shell_out(no_network, monkeypatch):
    def boom(*a, **kw):
        raise AssertionError("subprocess should not run")
    monkeypatch.setattr(r.subprocess, "run", boom)
    assert r._issues([]) == {}


def test_issues_builds_one_query_with_an_alias_per_issue(monkeypatch, no_network):
    seen = {}

    class Result:
        returncode = 0
        stdout = json.dumps({"data": {
            "i0": {"issue": {"createdAt": "2021-03-04T10:00:00Z"}},
            "i1": {"issue": {"createdAt": "2023-11-02T10:00:00Z"}},
        }})

    def fake_run(argv, **kw):
        seen["argv"] = argv
        return Result()

    monkeypatch.setattr(r.subprocess, "run", fake_run)
    got = r._issues([("antigenomics", "vdjdb-db", "12"), ("antigenomics", "vdjdb-db", "34")])
    assert got == {"https://github.com/antigenomics/vdjdb-db/issues/12": 2021,
                   "https://github.com/antigenomics/vdjdb-db/issues/34": 2023}
    query = seen["argv"][-1]
    assert query.count("repository(") == 2, "one query, one alias per issue"
    assert seen["argv"][:2] == ["gh", "api"]


def test_issues_returns_nothing_when_gh_fails(monkeypatch, no_network):
    class Failed:
        returncode = 1
        stdout = ""
    monkeypatch.setattr(r.subprocess, "run", lambda *a, **kw: Failed())
    assert r._issues([("o", "r", "1")]) == {}


def test_resolve_keeps_both_spellings_of_one_pmid(monkeypatch):
    """`PMID:34433824` and `PMID: 34433824` both occur, on 22 records.

    A dict keyed on the id would keep one and drop the other from the year plots, which is the
    failure this module exists to end.
    """
    stub_get(monkeypatch, {"result": {"uids": ["34433824"], "34433824": {"pubdate": "2021 Aug"}}})
    monkeypatch.setattr(r, "_issues", lambda refs: {})
    out = r.resolve(["PMID:34433824", "PMID: 34433824"])
    assert set(out["reference.id"]) == {"PMID:34433824", "PMID: 34433824"}
    assert set(out["year"]) == {2021}
    assert set(out["source"]) == {"pubmed"}


def test_resolve_drops_blank_reference_ids(monkeypatch, no_network):
    """854 records carry no reference at all. That is a curation gap, not a lookup failure."""
    monkeypatch.setattr(r, "_pubmed", lambda ids: {})
    monkeypatch.setattr(r, "_issues", lambda refs: {})
    out = r.resolve(["", "   ", None])
    assert out.height == 0
    assert out.columns == ["reference.id", "year", "source"]


def test_resolve_is_sorted_and_unique(monkeypatch, no_network):
    monkeypatch.setattr(r, "_pubmed", lambda ids: {})
    monkeypatch.setattr(r, "_issues", lambda refs: {})
    out = r.resolve(["https://arxiv.org/abs/2104.05918",
                     "https://doi.org/10.1101/2022.03.04.483006",
                     "https://arxiv.org/abs/2104.05918"])
    assert out["reference.id"].to_list() == sorted(out["reference.id"].to_list())
    assert out["reference.id"].n_unique() == out.height


def test_refresh_writes_the_table(tmp_path, monkeypatch, no_network):
    import polars as pl
    monkeypatch.setattr(r, "_pubmed", lambda ids: {})
    monkeypatch.setattr(r, "_issues", lambda refs: {})
    records = pl.DataFrame({"reference.id": ["https://arxiv.org/abs/2104.05918", ""]})
    out = tmp_path / "nested" / "reference_years.tsv"
    got = r.refresh(records, out)
    assert out.exists()
    assert out.read_text().startswith("reference.id\tyear\tsource")
    assert got["year"].to_list() == [2021]


def test_load_says_what_to_run_when_the_table_is_missing(tmp_path):
    """A silent fallback to a hardcoded table is how the years came to be four years stale."""
    with pytest.raises(FileNotFoundError, match="vdjdb refs"):
        r.load(tmp_path / "absent.tsv")
