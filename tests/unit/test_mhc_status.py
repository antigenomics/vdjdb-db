"""Every MHC call resolves against one of two authorities, or the build stops.

`assert_mhc_resolves` is the gate the corpus already satisfies: measured over 192,793 records, no
blank MHC cell and no unresolved name in either column. These tests fix that as a contract rather
than a property of today's data.
"""
from __future__ import annotations

import polars as pl
import pytest

from vdjdb.assemble.epitopes import (
    _hla_prefixes,
    _nonhuman_names,
    assert_mhc_resolves,
    build_restriction,
    mhc_status,
)
from vdjdb.config import Paths

ROOT = Paths.discover().root


def _status(*values: str) -> list[str]:
    df = pl.DataFrame({"mhc.a": list(values)})
    return df.with_columns(mhc_status("mhc.a"))["mhc.a.status"].to_list()


def test_an_expression_suffix_sits_on_the_terminal_field() -> None:
    """Both spellings of a null allele resolve.

    Splitting `HLA-A*24:09N` on `:` yields `HLA-A*24` and `HLA-A*24:09N`, never `HLA-A*24:09`, so a
    prefix set built that way reports the two-field name of a null allele as naming nothing. 2,442
    prefixes are reachable only through this path.
    """
    assert "HLA-A*24:09" in _hla_prefixes(ROOT)
    assert "HLA-A*24:09N" in _hla_prefixes(ROOT)
    assert _status("HLA-A*24:09", "HLA-A*24:09N") == ["known", "known"]


def test_the_three_field_depths_all_resolve() -> None:
    """A VDJdb call is one, two or three fields; the authority stores four."""
    assert _status("HLA-A*02", "HLA-A*02:01", "HLA-B*35:42:01") == ["known"] * 3


def test_a_name_outside_the_hla_database_needs_a_declaration() -> None:
    assert _status("H2-Kb", "B2M", "Mamu-A*01") == ["declared"] * 3
    # Same shape, not declared: a murine haplotype nobody has recorded.
    assert _status("H2-Zz9") == ["unknown"]


def test_a_well_formed_allele_that_does_not_exist_is_unknown() -> None:
    """The gate's whole purpose: spelling that parses is not spelling that names something."""
    assert _status("HLA-A*99:99", "HLA-A*08:01", "HLA-B*12") == ["unknown"] * 3


def test_the_declared_vocabulary_covers_every_non_hla_name_in_the_corpus() -> None:
    """A row that names nothing in the corpus is as much a defect as a corpus name with no row."""
    from vdjdb.curate.nomenclature import harmonise_mhc
    from vdjdb.io.chunks import read_chunks

    harmonised, _ = harmonise_mhc(read_chunks())
    used = {v for c in ("mhc.a", "mhc.b") for v in harmonised[c].unique()
            if v and not v.startswith("HLA-")}
    declared = set(_nonhuman_names(ROOT))
    assert used - declared == set(), "corpus names with no declaration"
    assert declared - used == set(), "declarations naming nothing in the corpus"


def test_the_gate_names_the_value_the_column_and_the_chunk() -> None:
    df = pl.DataFrame({"mhc.a": ["HLA-A*02:01", "HLA-A*99:99", ""],
                       "mhc.b": ["B2M"] * 3,
                       "chunk.file": ["PMID_1.txt", "PMID_2.txt", "PMID_3.txt"]})
    with pytest.raises(ValueError) as e:
        assert_mhc_resolves(df)
    message = str(e.value)
    assert "2 MHC value(s)" in message
    for expected in ("HLA-A*99:99", "PMID_2.txt", "<blank>", "PMID_3.txt", "mhc.a",
                     "patches/mhc.dict", "proofreading/mhc_nonhuman.tsv"):
        assert expected in message
    assert "HLA-A*02:01" not in message and "PMID_1.txt" not in message


def test_build_restriction_refuses_to_assemble_an_unresolved_call() -> None:
    """The gate runs where the status is computed, so no path reaches a release around it."""
    df = pl.DataFrame({"antigen.epitope": ["GILGFVFTL"], "antigen.species": ["InfluenzaA"],
                       "mhc.a": ["HLA-A*08:01"], "mhc.b": ["B2M"], "mhc.class": ["MHCI"],
                       "reference.id": ["PMID:1"], "chunk.file": ["PMID_1.txt"]})
    with pytest.raises(ValueError, match="HLA-A\\*08:01"):
        build_restriction(df)


# -- IPD's Confirmed / Unconfirmed, as a fourth status (#634) ------------------------------------

def test_a_call_with_no_confirmed_allele_under_it_reads_unconfirmed() -> None:
    """The second of two questions. #624's gate asks whether IPD carries the name at any field depth
    and all four of these pass it; this asks whether anyone other than the submitter believes the
    allele exists, which for a two-field call in a specificity database is the more useful one.
    """
    assert _status("HLA-A*02:266", "HLA-B*57:06", "HLA-DQA1*01:11") == ["unconfirmed"] * 3


def test_a_call_with_a_confirmed_allele_under_it_stays_known() -> None:
    """Prefix-wise and *any*: `HLA-A*02` names thousands of alleles and the question is whether any
    is confirmed, not whether all are.
    """
    assert _status("HLA-A*02", "HLA-A*02:01", "HLA-DRB1*15:01") == ["known"] * 3


def test_unconfirmed_is_a_refinement_of_known_and_never_fatal() -> None:
    """4 calls over 5 records. A gate here would reject almost nothing at a cost to those records,
    and an unconfirmed allele is a real name (#634).
    """
    from vdjdb.assemble.epitopes import assert_mhc_resolves

    records = pl.DataFrame({"mhc.a": ["HLA-A*02:266"], "mhc.b": ["B2M"],
                            "chunk.file": ["PMID_1.txt"]})
    assert_mhc_resolves(records)        # must not raise


def test_a_class_two_allele_written_at_four_fields_still_reaches_its_groove():
    """ROADMAP phase 9e, and the workaround for `antigenomics/mhcmatch#3`.

    `mhcmatch`'s `resolve_allele` trims to two fields and `class2_key` does not, so the same molecule
    resolves or does not depending on how deep the submitter's spelling went. 29 of the corpus's 354
    class II pairs are only reachable after the trim. Delete this with the workaround when the
    upstream fix ships.
    """
    from vdjdb.curate.presentation import resolve

    shallow = resolve("HLA-DRA*01:01", "HLA-DRB1*11:01", "MHCII")
    deep = resolve("HLA-DRA*01:02:03", "HLA-DRB1*11:01:02", "MHCII")
    assert shallow == ("DRB1_1101", "exact")
    assert deep == shallow, "a deeper spelling of one molecule is the same groove"


def test_a_murine_class_two_molecule_has_no_hla_pseudosequence_and_says_so():
    """Not a defect in the corpus. `mhcmatch`'s class II pseudosequences are HLA.

    Reported rather than passed over, because a reader needs to know which pairs a presentation
    model cannot score at all - that is a coverage statement, not a finding against the record.
    """
    from vdjdb.curate.presentation import resolve

    assert resolve("H2-IAb", "H2-IAb", "MHCII") == ("", "none")


def test_one_molecule_filed_under_two_classes_names_which_side_is_the_outlier():
    """The check that needs no authority at all, and the one that found `H2-IAb` on 77 records."""
    import polars as pl

    from vdjdb.curate.presentation import report

    restriction = pl.DataFrame({
        "antigen.epitope": ["QVYSLIRPNENPAH", "AAAAAAAAAAAAAAA", "SIINFEKL"],
        "antigen.species": ["", "", ""],
        "mhc.a": ["H2-IAb", "H2-IAb", "H2-Kb"],
        "mhc.b": ["H2-IAb", "H2-IAb", "B2M"],
        "mhc.class": ["MHCI", "MHCII", "MHCI"],
        "records": [77, 537, 10], "references": [1, 20, 1],
    })
    flagged = report(restriction)
    outlier = flagged.filter(pl.col("mhc.class") == "MHCI").filter(pl.col("mhc.a") == "H2-IAb")
    assert outlier.height == 1
    assert "molecule is also MHCII on 1 of 2" in outlier["finding"].item()
    assert "epitope is 14 residues" in outlier["finding"].item()
    assert "H2-Kb" not in flagged["mhc.a"].to_list(), "a consistent class I pair is not a finding"


@pytest.mark.parametrize("a,b,cls", [
    ("HLA-DRA*01:01", "HLA-DRB1*01:01", "MHCI"),
    ("H2-IAb", "H2-IAb", "MHCI"),
    ("HLA-A*02:01", "B2M", "MHCII"),
])
def test_mhc_class_matches_both_harmonised_chains(a, b, cls):
    from vdjdb.assemble.epitopes import assert_mhc_class

    records = pl.DataFrame({"mhc.a": [a], "mhc.b": [b], "mhc.class": [cls],
                            "chunk.file": ["PMID_1.tsv"]})
    with pytest.raises(ValueError, match="MHC class disagrees"):
        assert_mhc_class(records)
    corrected = "MHCII" if cls == "MHCI" else "MHCI"
    assert_mhc_class(records.with_columns(pl.lit(corrected).alias("mhc.class")))
