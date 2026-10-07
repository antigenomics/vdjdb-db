"""The report a chunk pull request gets, which is about the corpus and not about the file.

Every assertion here is that a number is relative: a score that only exists because another chunk
reports the same clonotype, a value that is "new" only because no other chunk carries it. A report
computed from the file alone could not make any of these statements, which is the reason this reads
assembled records.
"""
from __future__ import annotations

import polars as pl
import pytest

from vdjdb.curate.submission import report

COLUMNS = ("chunk.file", "vdjdb.score", "cdr3.alpha", "cdr3.beta", "antigen.epitope",
           "antigen.species", "antigen.gene", "mhc.a", "mhc.b", "species", "reference.id")


def rows(*specs: tuple) -> pl.DataFrame:
    return pl.DataFrame(list(specs), schema=list(COLUMNS), orient="row")


NEW = ("new.txt", 0, "CAAA", "CBBB", "SLLMWITQC", "HomoSapiens", "CTAG1B",
       "HLA-A*02:01", "", "HomoSapiens", "PMID:1")
OLD = ("old.txt", 1, "CCCC", "CDDD", "GILGFVFTL", "InfluenzaA", "M",
       "HLA-A*02:01", "", "HomoSapiens", "PMID:2")


def test_a_value_no_other_chunk_carries_is_reported_as_new() -> None:
    text = report(["new.txt"], rows(NEW, OLD))
    assert "`SLLMWITQC`" in text and "`CTAG1B`" in text
    assert "`GILGFVFTL`" not in text          # the other chunk's epitope is not this chunk's finding
    assert "`HLA-A*02:01`" not in text        # shared, so not new


def test_a_value_the_corpus_already_has_is_not_new() -> None:
    shared = (*NEW[:4], "GILGFVFTL", "InfluenzaA", "M", *NEW[7:])
    text = report(["new.txt"], rows(shared, OLD))
    assert "| `antigen.epitope` | 0 |" in text


def test_the_score_histogram_is_per_chunk() -> None:
    text = report(["new.txt", "old.txt"], rows(NEW, OLD))
    assert "| `new.txt` | 1 | 1 | 0 | 0 | 0 |" in text
    assert "| `old.txt` | 1 | 0 | 1 | 0 | 0 |" in text


def test_a_recurrence_requires_source_evidence_for_independence() -> None:
    """Publication labels alone do not establish independent experiments."""
    echo = ("new.txt", 0, OLD[2], OLD[3], OLD[4], OLD[5], OLD[6], OLD[7], "", "HomoSapiens", "PMID:3")
    text = report(["new.txt"], rows(echo, OLD))
    assert "**1 record(s) repeat a clonotype/pMHC" in text
    assert "Check the source experiment" in text
    assert "independent replication, not duplication" not in text


def test_no_record_from_the_changed_file_says_so_rather_than_printing_zeroes() -> None:
    """A renamed or excluded file, which would otherwise render an empty table reading as success."""
    text = report(["absent.txt"], rows(OLD))
    assert "No assembled records come from the changed files" in text


def test_a_path_is_matched_on_its_basename() -> None:
    """The workflow passes `chunks/PMID_1.txt`; `chunk.file` holds the bare name."""
    assert "| `new.txt` |" in report(["chunks/new.txt"], rows(NEW, OLD))


@pytest.mark.parametrize("column", ["antigen.gene", "reference.id"])
def test_a_column_missing_from_the_frame_is_skipped_not_fatal(column: str) -> None:
    """`build_master` is the only caller today, but a report must not crash on a narrower frame."""
    text = report(["new.txt"], rows(NEW, OLD).drop(column))
    assert f"`{column}`" not in text and "Values new to VDJdb" in text


# ---------------------------------------------------------------------------------------------
# Look-alike values: an alert, never a gate
# ---------------------------------------------------------------------------------------------

def test_a_value_differing_only_in_case_is_flagged_against_the_one_it_should_have_joined() -> None:
    near = ("new.txt", 0, "CAAA", "CBBB", "SLLMWITQC", "InfluenzaA", "m", *NEW[7:])
    text = report(["new.txt"], rows(near, OLD))
    assert "`antigen.gene`: **`m`** against `M`, already in VDJdb" in text
    assert "blocks nothing" in text.lower()


def test_a_separator_difference_is_the_same_finding() -> None:
    old = (*OLD[:6], "MAGE-A3", *OLD[7:])
    near = ("new.txt", 0, "CAAA", "CBBB", "SLLMWITQC", "HomoSapiens", "MAGEA3", *NEW[7:])
    assert "**`MAGEA3`** against `MAGE-A3`" in report(["new.txt"], rows(near, old))


def test_plus_and_minus_are_never_folded_together() -> None:
    """Measured false positives: `CD8+CD95-` and `CD8+CD95+` are two populations, not one typo.

    They are not in `NOVELTY_COLUMNS`, so the guard is on the folding function itself - the next column
    added to that tuple inherits it.
    """
    from vdjdb.curate.submission import _fold, collisions

    assert _fold("CD8+CD27+CD45RA+CD95-") != _fold("CD8+CD27+CD45RA+CD95+")
    assert collisions(["CD8+CD95-", "CD8+CD95+"]) == {}


def test_a_genuinely_new_value_is_not_reported_as_a_look_alike() -> None:
    text = report(["new.txt"], rows(NEW, OLD))
    assert "Values new to VDJdb" in text and "probably a spelling" not in text.lower()


def test_the_report_never_signals_failure() -> None:
    """It returns text. There is no exit code, no raise, and no caller that could read one."""
    near = ("new.txt", 0, "CAAA", "CBBB", "SLLMWITQC", "InfluenzaA", "m", *NEW[7:])
    assert isinstance(report(["new.txt"], rows(near, OLD)), str)


# ---------------------------------------------------------------------------------------------
# The corpus-wide scan
# ---------------------------------------------------------------------------------------------

def test_the_scan_reports_one_row_per_spelling_with_its_record_count() -> None:
    from vdjdb.curate.submission import lookalikes

    frame = rows(NEW, (*OLD[:6], "m", *OLD[7:]), (*OLD[:6], "M", *OLD[7:]))
    got = lookalikes(frame).filter(pl.col("column") == "antigen.gene")
    assert sorted(got["spelling"].to_list()) == ["M", "m"]
    assert got["records"].sum() == 2


def test_the_governing_species_for_a_gene_is_the_antigens_not_the_donors() -> None:
    """Reading `species` here inverts the answer, which is what it did before this was fixed.

    `Nef` and `NEF` sit on one donor `species` and on two `antigen.species` - HIV-1 and InfluenzaA -
    and only the second says the two spellings might be two different genes.
    """
    from vdjdb.curate.submission import lookalikes

    frame = rows(("a.txt", 0, "CA", "CB", "EP1", "HIV-1", "Nef", *NEW[7:]),
                 ("b.txt", 0, "CC", "CD", "EP2", "InfluenzaA", "NEF", *NEW[7:]))
    assert frame["species"].n_unique() == 1                      # one donor species
    assert not lookalikes(frame)["same.species"].any()           # two antigen species


def test_one_antigen_species_makes_it_a_within_species_finding() -> None:
    from vdjdb.curate.submission import lookalikes

    frame = rows(("a.txt", 0, "CA", "CB", "EP1", "CMV", "IE1", *NEW[7:]),
                 ("b.txt", 0, "CC", "CD", "EP2", "CMV", "IE-1", *NEW[7:]))
    assert lookalikes(frame)["same.species"].all()


def test_a_corpus_with_no_look_alikes_gives_an_empty_frame_with_the_schema() -> None:
    """`build` writes this file every run, so an empty one still has to be readable."""
    from vdjdb.curate.submission import lookalikes

    got = lookalikes(rows(NEW, OLD))
    assert got.is_empty()
    assert got.columns == ["column", "folded", "spelling", "records", "spellings", "same.species"]


def test_a_peptide_under_two_species_is_reported_with_its_source_count() -> None:
    """#633. `epitopes` is keyed on `(epitope, species)`, so this is two rows there by design."""
    from vdjdb.curate.submission import epitope_sources

    frame = rows(("a.txt", 0, "CA", "CB", "VEALYLVCG", "HomoSapiens", "INS", *NEW[7:]),
                 ("b.txt", 0, "CC", "CD", "VEALYLVCG", "MusMusculus", "Ins2", *NEW[7:]))
    got = epitope_sources(frame)
    assert got.height == 2
    assert got["sources"].to_list() == [2, 2]
    assert got["genes"].to_list() == [1, 1]
    assert got["conserved"].all(), "VEALYLVCG is in human INS and mouse Ins2"


def test_two_gene_labels_for_one_peptide_and_species_are_reported_too() -> None:
    """`build_epitopes` keeps the modal label and said a function listed the rest. It did not exist,
    so the discarded labels were reported nowhere -- including 187 on one peptide."""
    from vdjdb.curate.submission import epitope_sources

    frame = rows(("a.txt", 0, "CA", "CB", "ESDPIVAQY", "HomoSapiens", "TTN", *NEW[7:]),
                 ("b.txt", 0, "CC", "CD", "ESDPIVAQY", "HomoSapiens", "TTN", *NEW[7:]),
                 ("c.txt", 0, "CE", "CF", "ESDPIVAQY", "HomoSapiens", "TITIN", *NEW[7:]))
    got = epitope_sources(frame)
    assert got.height == 1                       # one (epitope, species)
    assert got["sources"][0] == 1 and got["genes"][0] == 2
    assert got["antigen.gene"][0] == "TTN", "the modal label, matching what the catalogue keeps"
    assert not got["conserved"][0]


def test_a_single_sourced_peptide_is_not_reported() -> None:
    from vdjdb.curate.submission import epitope_sources

    got = epitope_sources(rows(NEW, OLD))
    assert got.is_empty()
    assert got.columns == ["antigen.epitope", "antigen.species", "antigen.gene",
                           "sources", "genes", "records", "conserved"]
