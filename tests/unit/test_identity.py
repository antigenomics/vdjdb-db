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


def test_a_reference_id_correction_amends_rather_than_retiring_the_record():
    """#685. The pass bucketed candidates on `(chunk_file, reference.id)`, building the bucket from
    the previous build's reference and reading it with the new one, so a change to `reference.id`
    itself found an empty bucket and was retired instead of amended -- silently. Closing the space
    in `PMID: 34433824` cost 22 published ids that way."""
    out1, reg, _ = reconcile(frame({"reference.id": "PMID: 34433824"}), IdentityRegistry(),
                             release="v1")
    original = out1["record_id"][0]
    out2, reg, rep = reconcile(frame({"reference.id": "PMID:34433824"}), reg, release="v2")

    assert out2["record_id"][0] == original, "a reference respelling must not mint a new id"
    assert not rep.added and not rep.retired
    assert len(rep.amended) == 1
    _rid, field, old, new = rep.amended[0]
    assert (field, old, new) == ("reference.id", "PMID: 34433824", "PMID:34433824")


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


def test_two_records_one_field_apart_amend_by_row_when_no_line_moved():
    """The bijection that makes `chunk.row` exact rather than a guess.

    Both rows are one field from both registry entries, so the field's value cannot say which is
    which - and `test_ambiguous_amendment_is_refused` is right that a row number alone cannot either,
    because deleting a line shifts every number after it. What settles it is that the two unmatched
    rows and the two candidates occupy the *same* pair of row numbers: nothing moved, so the pairing
    is forced.

    The case: `menon_etal_2024.txt` rows 26 and 27 carry `TRBV5-3;TRBV5-5;TRBV5-8` and
    `TRBV5-3;TRBV5-8`, and normalising `;` to `,` moved both keys by that one field. Without this,
    `VDJDB0000187889` and `...890` retire and two fresh ids are minted - two published identifiers
    lost to a separator.
    """
    before = frame({"v.beta": "TRBV5-3;TRBV5-5;TRBV5-8", "chunk.row": 26},
                   {"v.beta": "TRBV5-3;TRBV5-8", "chunk.row": 27})
    out, reg, _ = reconcile(before, IdentityRegistry(), release="v1")
    was = dict(zip(before["chunk.row"].to_list(), out["record_id"].to_list(), strict=True))

    after = frame({"v.beta": "TRBV5-3,TRBV5-5,TRBV5-8", "chunk.row": 26},
                  {"v.beta": "TRBV5-3,TRBV5-8", "chunk.row": 27})
    out, reg, rep = reconcile(after, reg, release="v2")
    assert not rep.added and not rep.retired, rep
    assert len(rep.amended) == 2 and {a[1] for a in rep.amended} == {"v.beta"}
    assert dict(zip(after["chunk.row"].to_list(), out["record_id"].to_list(), strict=True)) == was


def test_a_bulk_repair_amends_every_row_it_touches():
    """The case the first version of the row tie-break got wrong, and the reason it cannot be narrower.

    Repairing the 4,324 junctions of vdjdb-db#646 moves that many keys at once, so one chunk's bucket
    holds hundreds of unmatched rows while each row has only its own two or three candidates. A test
    that asked for a bijection between the bucket's rows and *one row's* candidates never held, and
    every one of them fell through to a fresh id: 2,066 published identifiers retired and re-minted.

    What makes the row number an identity is that the file still has a row at every number the
    candidates sit on, which an edit in place always does.
    """
    rows = [{"cdr3.beta": f"CASS{a}{b}F", "chunk.row": i}
            for i, (a, b) in enumerate([("A", "A"), ("A", "B"), ("B", "A"), ("B", "B")])]
    before = frame(*rows)
    out, reg, _ = reconcile(before, IdentityRegistry(), release="v1")
    was = dict(zip(before["chunk.row"].to_list(), out["record_id"].to_list(), strict=True))

    # Every row gains the trailing anchor its germline encodes - one field, every row, in place.
    after = frame(*[{**r, "cdr3.beta": r["cdr3.beta"] + "F"} for r in rows])
    out, _, rep = reconcile(after, reg, release="v2")
    assert not rep.added and not rep.retired, rep
    assert len(rep.amended) == 4 and {a[1] for a in rep.amended} == {"cdr3.beta"}
    assert dict(zip(after["chunk.row"].to_list(), out["record_id"].to_list(), strict=True)) == was


def test_the_row_tiebreak_is_refused_when_a_line_moved():
    """The other half: three rows become two, so the row numbers no longer pair up and the ambiguity
    stands. A shifted line must not be read as an amendment of whatever now sits on its number.
    """
    before = frame({"v.beta": "TRBV5-3;TRBV5-5", "chunk.row": 0},
                   {"v.beta": "TRBV5-3;TRBV5-8", "chunk.row": 1},
                   {"v.beta": "TRBV5-3;TRBV5-9", "chunk.row": 2})
    _, reg, _ = reconcile(before, IdentityRegistry(), release="v1")
    after = frame({"v.beta": "TRBV5-3,TRBV5-5", "chunk.row": 0},
                  {"v.beta": "TRBV5-3,TRBV5-8", "chunk.row": 1})
    _, _, rep = reconcile(after, reg, release="v2")
    assert len(rep.amended) < 2, rep.amended


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

    assert list(NATURAL_KEY) == [*CHUNK_DEDUP_KEY, "chunk.file"], (
        "identity is the dedup key plus the chunk: a chunk is one paper, so two "
        "chunks are two independent reports")


def test_records_differing_only_in_donor_are_distinct_records():
    df = frame({"meta.subject.id": "donor1"}, {"meta.subject.id": "donor2"})
    out, reg, rep = reconcile(df, IdentityRegistry(), release="v1")
    assert out["record_id"].n_unique() == 2
    assert len(reg) == 2 and len(rep.added) == 2


def test_matching_rows_in_two_chunks_are_two_records():
    """A chunk is one paper, so two chunks are two independent reports -- never one record.

    19 pairs in the real corpus match field for field across chunks. Collapsing them would delete
    exactly the independent-replication signal phase 11 fits motif clustering against.
    """
    df = frame({"chunk.file": "PMID_A.txt"}, {"chunk.file": "PMID_B.txt"})
    out, reg, rep = reconcile(df, IdentityRegistry(), release="v1")
    assert out["record_id"].n_unique() == 2
    assert len(reg) == 2 and len(rep.added) == 2


def test_independent_reports_stay_stable_across_rebuilds():
    df = frame({"chunk.file": "PMID_A.txt"}, {"chunk.file": "PMID_B.txt"})
    out1, reg, _ = reconcile(df, IdentityRegistry(), release="v1")
    out2, reg, rep = reconcile(df, reg, release="v2")
    assert out1["record_id"].to_list() == out2["record_id"].to_list()
    assert not rep.added and not rep.retired


# ---------------------------------------------------------------------------------------------
# A build with no registry has to say so
# ---------------------------------------------------------------------------------------------

def test_a_missing_registry_warns_that_the_build_is_not_id_stable(tmp_path) -> None:
    """`ROADMAP.md` 10.4 claims the build reports this. Nothing did, and it is the state today.

    Silent is the failure mode that matters: ids allocated from 1 look exactly like ids reconciled
    against a registry, and the difference only appears when someone compares two releases.
    """
    import warnings

    import polars as pl

    from vdjdb.assemble.master import add_record_ids
    from vdjdb.schema import ALL_COLUMNS

    df = pl.DataFrame({**{c: [""] for c in ALL_COLUMNS}, "chunk.file": ["a.txt"], "chunk.row": [1]})
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        add_record_ids(df, registry=tmp_path / "absent.tsv")
    messages = [str(w.message) for w in caught if issubclass(w.category, UserWarning)]
    assert any("not id-stable" in m for m in messages), messages


def test_an_existing_registry_does_not_warn(tmp_path) -> None:
    import warnings

    import polars as pl

    from vdjdb.assemble.master import add_record_ids
    from vdjdb.identity.ids import REGISTRY_COLUMNS
    from vdjdb.schema import ALL_COLUMNS

    registry = tmp_path / "records.tsv"
    registry.write_text("\t".join(REGISTRY_COLUMNS) + "\n")
    df = pl.DataFrame({**{c: [""] for c in ALL_COLUMNS}, "chunk.file": ["a.txt"], "chunk.row": [1]})
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        add_record_ids(df, registry=registry)
    assert not [w for w in caught if "id-stable" in str(w.message)]


# ---------------------------------------------------------------------------------------------
# A registry a build wrote has to be a registry a build can read
# ---------------------------------------------------------------------------------------------

def _frame(**over):
    import polars as pl

    from vdjdb.schema import ALL_COLUMNS

    base = {c: [""] for c in ALL_COLUMNS}
    base.update({k: [v] for k, v in over.items()})
    return pl.DataFrame({**base, "chunk.file": ["a.txt"], "chunk.row": [1]})


def test_a_registry_written_by_a_build_is_matched_by_the_next_build(tmp_path) -> None:
    """The round trip. Writing it from anywhere else in the pipeline silently breaks every key.

    `NATURAL_KEY` carries `cdr3.alpha` and `cdr3.beta`; `add_record_ids` runs before `fix_cdr3` so a
    repair does not mint a new id. A registry seeded from `build_master`'s *return* value is therefore
    keyed on repaired sequences, nothing matches pass 1, and the amendment pass compares every row
    against a bucket of leftovers - measured at minutes against three seconds on the real corpus.
    """
    from vdjdb.assemble.master import add_record_ids

    path = tmp_path / "records.tsv"
    df = _frame(**{"cdr3.beta": "CASSLAPGATNEKLFF", "antigen.epitope": "GILGFVFTL"})
    first = add_record_ids(df, registry=path, write=path)
    assert path.exists()
    again = add_record_ids(df, registry=path)
    assert again["record_id"].to_list() == first["record_id"].to_list()


def test_a_registry_keyed_on_repaired_sequences_records_an_amendment_nobody_made(tmp_path) -> None:
    """Why the write has to happen inside the build, stated as the observable consequence.

    The amendment pass does recover the id - a one-field change in the same chunk and reference is an
    amendment by design - so a registry seeded from the post-repair frame does not corrupt ids. What it
    does is route every record through that pass, which costs minutes instead of seconds on the real
    corpus and records an amendment for a record no curator touched.
    """
    from vdjdb.assemble.master import add_record_ids
    from vdjdb.identity.ids import IdentityRegistry

    path = tmp_path / "records.tsv"
    raw = _frame(**{"cdr3.beta": "ASSLAPGATNEKLF"})           # as submitted, anchors trimmed
    repaired = _frame(**{"cdr3.beta": "CASSLAPGATNEKLFF"})    # what fix_cdr3 would produce
    add_record_ids(repaired, registry=path, write=path)       # seeded from the wrong frame
    kept = add_record_ids(raw, registry=path, write=path)     # what the next build actually passes
    entry = next(iter(IdentityRegistry.load(path)._by_key.values()))
    assert kept["record_id"][0] == entry.record_id            # the id survives, via the amendment pass
    assert entry.amendment_count == 1                         # and an amendment is recorded regardless


def test_a_correctly_seeded_registry_records_no_amendment(tmp_path) -> None:
    """The control for the test above: same two calls, same frame both times."""
    from vdjdb.assemble.master import add_record_ids
    from vdjdb.identity.ids import IdentityRegistry

    path = tmp_path / "records.tsv"
    raw = _frame(**{"cdr3.beta": "ASSLAPGATNEKLF"})
    add_record_ids(raw, registry=path, write=path)
    add_record_ids(raw, registry=path, write=path)
    assert next(iter(IdentityRegistry.load(path)._by_key.values())).amendment_count == 0


def _entry(rid: str, key: tuple[str, ...], *, row: int, chunk: str = "c.txt",
           state: str = "active", replaced_by: str = ""):
    from vdjdb.identity.ids import _Entry, _hash, _pack_note
    return _Entry(record_id=rid, state=state, natural_key_hash=_hash(list(key)),
                  content_hash=_hash(list(key)), chunk_file=chunk, chunk_row=row,
                  first_seen_release="dev", first_seen_commit="", last_seen_release="dev",
                  last_modified_release="dev", last_modified_commit="", amendment_count=0,
                  amended_from_key_hash="", replaced_by=replaced_by, note=_pack_note(key))


def _key(**over) -> tuple[str, ...]:
    from vdjdb.identity.ids import NATURAL_KEY
    base = dict.fromkeys(NATURAL_KEY, "")
    base["cdr3.alpha"], base["cdr3.beta"] = "CAVX", "CASSY"
    base["chunk.file"] = "c.txt"
    base.update(over)
    return tuple(base[c] for c in NATURAL_KEY)


def test_a_retirement_that_moved_two_key_fields_points_at_the_id_that_took_over():
    """#693, and `ROADMAP.md` §10.4: an id that vanishes is the failure nobody can diagnose.

    Two key fields moving is a new record by the rule the amendment pass is built on, so the
    retirement is right. What was missing is the pointer from the old id to the new one.
    """
    from vdjdb.identity.ids import IdentityRegistry, RecordState, reconcile
    from vdjdb.schema import ALL_COLUMNS

    before = _key(**{"antigen.gene": "Plod1", "meta.subject.cohort": "C1"})
    registry = IdentityRegistry([_entry("VDJDB0000000001", before, row=1)])

    row = dict.fromkeys(ALL_COLUMNS, "")
    row |= {"cdr3.alpha": "CAVX", "cdr3.beta": "CASSY", "antigen.gene": "Plod2",
            "meta.subject.cohort": "", "chunk.file": "c.txt", "chunk.row": 1}
    _, registry, report = reconcile(pl.DataFrame([row]), registry, release="dev")

    assert report.amended == [], "two fields moved, so this is not an amendment"
    assert len(report.added) == 1 and len(report.retired) == 1
    old = next(e for e in registry._by_key.values() if e.record_id == "VDJDB0000000001")
    assert old.state == RecordState.RETIRED
    assert old.replaced_by == report.added[0]


def test_a_successor_is_refused_when_the_receptor_is_not_the_same_one():
    """A new record at the line an old one left is a coincidence until the CDR3s agree."""
    from vdjdb.identity.ids import IdentityRegistry, reconcile
    from vdjdb.schema import ALL_COLUMNS

    before = _key(**{"antigen.gene": "Plod1", "meta.subject.cohort": "C1"})
    registry = IdentityRegistry([_entry("VDJDB0000000001", before, row=1)])

    row = dict.fromkeys(ALL_COLUMNS, "")
    row |= {"cdr3.alpha": "CAVZZZ", "cdr3.beta": "CASSZZZ", "antigen.gene": "Plod2",
            "meta.subject.cohort": "", "chunk.file": "c.txt", "chunk.row": 1}
    _, registry, _report = reconcile(pl.DataFrame([row]), registry, release="dev")

    old = next(e for e in registry._by_key.values() if e.record_id == "VDJDB0000000001")
    assert old.replaced_by == "", "a different TCR at that line says nothing about the old id"


def test_two_retirements_from_one_line_link_neither():
    """One in, one out. With two on a side the pairing has no evidence, so nothing is written."""
    from collections import Counter

    from vdjdb.identity.ids import _successor

    a = _entry("VDJDB0000000001", _key(**{"antigen.gene": "Plod1"}), row=1)
    allocated = {("c.txt", 1): [("VDJDB0000000009", _key(**{"antigen.gene": "Plod2"}))]}
    assert _successor(a, allocated, Counter({("c.txt", 1): 1})) == "VDJDB0000000009"
    assert _successor(a, allocated, Counter({("c.txt", 1): 2})) == ""
    two = {("c.txt", 1): [("VDJDB0000000009", _key()), ("VDJDB0000000010", _key())]}
    assert _successor(a, two, Counter({("c.txt", 1): 1})) == ""

