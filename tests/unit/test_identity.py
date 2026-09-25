"""Record identity: allocation, amendment tracing, retirement, persistence."""
from __future__ import annotations

import polars as pl
import pytest

from vdjdb.identity import IdentityRegistry, canonical_content_hash, natural_key, reconcile
from vdjdb.identity.ids import NATURAL_KEY, RecordState, format_id, parse_id
from vdjdb.schema import ALL_COLUMNS

BASE = {
    "species": "HomoSapiens", "cdr3.alpha": "CAVSDLEPNSSASKIIF", "v.alpha": "TRAV12-2*01",
    "j.alpha": "TRAJ3*01", "cdr3.beta": "CASSIRSSYEQYF", "v.beta": "TRBV10-3*01",
    "d.beta": "", "j.beta": "TRBJ2-7*01", "mhc.a": "HLA-A*02:01", "mhc.b": "B2M",
    "mhc.class": "MHCI", "antigen.epitope": "GILGFVFTL", "antigen.gene": "M",
    "antigen.species": "InfluenzaA", "reference.id": "PMID:28629751",
}


def frame(*rows: dict) -> pl.DataFrame:
    full = []
    for i, r in enumerate(rows):
        d = {c: "" for c in ALL_COLUMNS}
        d.update(BASE)
        d.update(r)
        d["chunk.file"] = r.get("chunk.file", "PMID_28629751.txt")
        d["chunk.row"] = r.get("chunk.row", i)
        full.append(d)
    return pl.DataFrame(full)


def test_ids_are_allocated_monotonically_and_formatted():
    df = frame({}, {"cdr3.beta": "CASSLLLGGF"}, {"cdr3.beta": "CASSQQQGGF"})
    out, reg, rep = reconcile(df, IdentityRegistry(), release="v1")
    ids = out["record_id"].to_list()
    assert ids == [format_id(1), format_id(2), format_id(3)]
    assert [parse_id(i) for i in ids] == [1, 2, 3]
    assert len(rep.added) == 3 and rep.unchanged == 0
    assert len(reg) == 3


def test_rebuilding_unchanged_data_is_stable():
    df = frame({}, {"cdr3.beta": "CASSLLLGGF"})
    out1, reg, _ = reconcile(df, IdentityRegistry(), release="v1")
    out2, reg, rep = reconcile(df, reg, release="v2")
    assert out1["record_id"].to_list() == out2["record_id"].to_list()
    assert rep.unchanged == 2 and not rep.added and not rep.retired


def test_annotation_change_keeps_the_id_but_records_a_new_content_hash():
    """Re-annotating how a record was assayed does not make it a different record."""
    out1, reg, _ = reconcile(frame({}), IdentityRegistry(), release="v1")
    before = reg.to_frame().row(0, named=True)["content_hash"]

    out2, reg, rep = reconcile(frame({"method.identification": "tetramer-sort"}), reg, release="v2")
    after = reg.to_frame().row(0, named=True)["content_hash"]

    assert out1["record_id"][0] == out2["record_id"][0]
    assert before != after
    assert rep.annotated == 1 and not rep.added and not rep.retired


def test_typo_fix_is_traced_as_an_amendment_not_a_delete_plus_insert():
    """The case the whole design exists for."""
    out1, reg, _ = reconcile(frame({}), IdentityRegistry(), release="v1")
    original = out1["record_id"][0]

    fixed = frame({"cdr3.beta": "CASSIRSSYEQYFF"})  # one trailing F added
    out2, reg, rep = reconcile(fixed, reg, release="v2")

    assert out2["record_id"][0] == original, "a typo fix must not mint a new id"
    assert not rep.added and not rep.retired
    assert len(rep.amended) == 1
    rid, fld, old, new = rep.amended[0]
    assert rid == original
    assert fld == "cdr3.beta"
    assert old == "CASSIRSSYEQYF" and new == "CASSIRSSYEQYFF"
    assert reg.to_frame().row(0, named=True)["amendment_count"] == "1"


def test_two_field_change_is_a_new_record_not_an_amendment():
    """Amendment is deliberately conservative: one field, or it is a different record."""
    out1, reg, _ = reconcile(frame({}), IdentityRegistry(), release="v1")
    out2, reg, rep = reconcile(
        frame({"cdr3.beta": "CASSIRSSYEQYFF", "v.beta": "TRBV19*01"}), reg, release="v2")
    assert out2["record_id"][0] != out1["record_id"][0]
    assert len(rep.added) == 1 and len(rep.retired) == 1 and not rep.amended


def test_ambiguous_amendment_is_refused():
    """Two equally good candidates must not be guessed between -- a wrong link beats no link."""
    df = frame({"cdr3.beta": "CASSAAAAAF"}, {"cdr3.beta": "CASSBBBBBF"})
    _, reg, _ = reconcile(df, IdentityRegistry(), release="v1")
    # CASSCCCCCF differs from BOTH registry entries in exactly one field, so neither is a safe
    # match. The record must get a fresh id and both old ones must retire.
    _, reg, rep = reconcile(frame({"cdr3.beta": "CASSCCCCCF"}), reg, release="v2")
    assert not rep.amended
    assert len(rep.added) == 1 and len(rep.retired) == 2


def test_unambiguous_amendment_is_taken_even_with_other_records_present():
    """The sibling of the test above: one candidate, in a chunk that holds others."""
    df = frame({"cdr3.beta": "CASSAAAAAF"}, {"cdr3.beta": "CASSBBBBBF", "v.beta": "TRBV28*01"})
    _, reg, _ = reconcile(df, IdentityRegistry(), release="v1")
    # Only the first entry is one field away; the second also differs in v.beta.
    _out, reg, rep = reconcile(
        frame({"cdr3.beta": "CASSAAAAAFF"}, {"cdr3.beta": "CASSBBBBBF", "v.beta": "TRBV28*01"}),
        reg, release="v2")
    assert len(rep.amended) == 1 and rep.amended[0][1] == "cdr3.beta"
    assert not rep.added and not rep.retired


def test_removed_record_is_retired_not_deleted():
    df = frame({}, {"cdr3.beta": "CASSLLLGGF"})
    out1, reg, _ = reconcile(df, IdentityRegistry(), release="v1")
    kept = out1["record_id"][0]

    out2, reg, rep = reconcile(frame({}), reg, release="v2")
    assert out2["record_id"].to_list() == [kept]
    assert len(rep.retired) == 1
    states = dict(zip(reg.to_frame()["record_id"], reg.to_frame()["state"], strict=False))
    assert states[rep.retired[0]] == RecordState.RETIRED
    assert len(reg) == 2, "a retired record stays in the registry"


def test_retired_ids_are_never_reused():
    out1, reg, _ = reconcile(frame({}), IdentityRegistry(), release="v1")
    _, reg, _ = reconcile(frame({"cdr3.beta": "CASSZZZZZF", "v.beta": "TRBV28*01"}), reg, release="v2")
    _, reg, rep3 = reconcile(frame({"cdr3.beta": "CASSYYYYYF", "v.beta": "TRBV29-1*01"}), reg, release="v3")
    all_ids = reg.to_frame()["record_id"].to_list()
    assert len(all_ids) == len(set(all_ids)) == 3
    assert out1["record_id"][0] not in rep3.added


def test_registry_round_trips_through_tsv(tmp_path):
    _, reg, _ = reconcile(frame({}, {"cdr3.beta": "CASSLLLGGF"}), IdentityRegistry(), release="v1")
    p = tmp_path / "records.registry.tsv"
    reg.save(p)
    reloaded = IdentityRegistry.load(p)
    assert reloaded.to_frame().equals(reg.to_frame())

    # and it keeps allocating from where it left off
    _, reloaded, rep = reconcile(
        frame({}, {"cdr3.beta": "CASSLLLGGF"}, {"cdr3.beta": "CASSNEWNEWF"}), reloaded, release="v2")
    assert rep.added == [format_id(3)]


def test_registry_tsv_is_sorted_so_diffs_are_reviewable(tmp_path):
    rows = [{"cdr3.beta": f"CASS{i:04d}F"} for i in range(12)]
    _, reg, _ = reconcile(frame(*rows), IdentityRegistry(), release="v1")
    p = tmp_path / "r.tsv"
    reg.save(p)
    ids = pl.read_csv(p, separator="\t", infer_schema=False)["record_id"].to_list()
    assert ids == sorted(ids)


def test_natural_key_ignores_assay_annotation_but_content_hash_does_not():
    """How a record was assayed is annotation, not identity. Which donor it came from is identity."""
    a = dict.fromkeys(ALL_COLUMNS, "")
    a.update(BASE)
    annotated = dict(a, **{
        "method.identification": "tetramer-sort",
        "method.sequencing": "sanger",
        "meta.donor.MHC": "A02,B07",
        "meta.structure.id": "1AO7",
    })
    assert natural_key(a) == natural_key(annotated), "re-annotation must not change identity"
    assert canonical_content_hash(a) != canonical_content_hash(annotated), "but must be detected"

    other_donor = dict(a, **{"meta.subject.id": "donor2"})
    assert natural_key(a) != natural_key(other_donor), "a different donor is a different record"


@pytest.mark.parametrize("fld", NATURAL_KEY)
def test_every_natural_key_field_changes_the_key(fld):
    a = dict.fromkeys(ALL_COLUMNS, "")
    a.update(BASE)
    b = dict(a)
    b[fld] = (b.get(fld) or "") + "X"
    assert natural_key(a) != natural_key(b), f"{fld} is in NATURAL_KEY but does not affect it"


def test_natural_key_is_not_confused_by_a_tab_in_a_field():
    a = dict.fromkeys(ALL_COLUMNS, "")
    a.update(BASE, **{"antigen.gene": "A\tB"})
    b = dict.fromkeys(ALL_COLUMNS, "")
    b.update(BASE, **{"antigen.gene": "A", "antigen.species": "B"})
    assert canonical_content_hash(a) != canonical_content_hash(b)


def test_natural_key_is_the_chunk_dedup_key():
    """Identity and deduplication must agree on what "the same record" means.

    They drifted once: a narrower key that stopped at reference.id collided on 20,769 of the
    192,753 real records, because one paper reporting the same TCR against the same epitope in
    several donors is several records.
    """
    from vdjdb.schema import CHUNK_DEDUP_KEY

    assert NATURAL_KEY == CHUNK_DEDUP_KEY


def test_records_differing_only_in_donor_are_distinct_records():
    df = frame({"meta.subject.id": "donor1"}, {"meta.subject.id": "donor2"})
    out, reg, rep = reconcile(df, IdentityRegistry(), release="v1")
    assert out["record_id"].n_unique() == 2
    assert len(reg) == 2 and len(rep.added) == 2


def test_the_same_record_in_two_chunks_gets_one_id():
    """19 records in the real corpus are submitted in two chunks. One record, one id.

    Allocating a second id for the duplicate would also make it unstable: it would be re-added on
    every build, because the registry can only hold one entry per natural key.
    """
    df = frame({"chunk.file": "PMID_A.txt"}, {"chunk.file": "PMID_B.txt"})
    out, reg, rep = reconcile(df, IdentityRegistry(), release="v1")
    assert out["record_id"].n_unique() == 1
    assert len(reg) == 1 and len(rep.added) == 1
    assert rep.duplicated == [(out["record_id"][0], "PMID_B.txt")]


def test_duplicates_stay_stable_across_rebuilds():
    df = frame({"chunk.file": "PMID_A.txt"}, {"chunk.file": "PMID_B.txt"})
    out1, reg, _ = reconcile(df, IdentityRegistry(), release="v1")
    out2, reg, rep = reconcile(df, reg, release="v2")
    assert out1["record_id"].to_list() == out2["record_id"].to_list()
    assert not rep.added and not rep.retired
