"""The antigen nomenclature patch.

`curate/patch.py` was 45 % covered. Three of its behaviours were each chosen to reproduce a
measured legacy result, and none of them had a test: `NA` is a gene symbol rather than a missing
value, duplicate epitopes resolve to the last entry, and an empty patch cell falls back to the
curated value instead of blanking it.
"""
from __future__ import annotations

import polars as pl
import pytest

from vdjdb.curate.patch import (
    CONFLICTING_EPITOPES,
    apply_antigen_patch,
    conflicts,
    load_antigen_patch,
    render_patch_renames,
)

HEADER = "antigen.epitope\tantigen.gene\tantigen.species\n"


def write_patch(tmp_path, body):
    p = tmp_path / "antigen.dict"
    p.write_text(HEADER + body)
    return p


def records(rows):
    return pl.DataFrame(rows, schema={"antigen.epitope": pl.Utf8, "antigen.gene": pl.Utf8,
                                      "antigen.species": pl.Utf8}, orient="row")


def test_na_is_a_gene_symbol_not_a_missing_value(tmp_path):
    """Influenza neuraminidase. Read as missing, it blanked the gene on 14 rows."""
    p = write_patch(tmp_path, "GILGFVFTL\tNA\tInfluenzaA\n")
    patch = load_antigen_patch(p)
    assert patch["__gene"].to_list() == ["NA"]
    out = apply_antigen_patch(records([("GILGFVFTL", "M", "Influenza")]), patch)
    assert out["antigen.gene"].to_list() == ["NA"]


def test_a_duplicate_epitope_resolves_to_the_last_entry(tmp_path):
    """pandas' `.to_dict()` kept the last; taking the first moved 1,158 rows."""
    p = write_patch(tmp_path, "AAA\tfirst\tSpeciesA\nAAA\tlast\tSpeciesB\n")
    patch = load_antigen_patch(p)
    assert patch.height == 1
    assert patch["__gene"].to_list() == ["last"]
    assert patch["__species"].to_list() == ["SpeciesB"]


def test_an_empty_patch_cell_keeps_the_curated_value(tmp_path):
    """A patch entry asserts something about the epitope, not about every column of the record."""
    p = write_patch(tmp_path, "AAA\t\tSpeciesB\n")
    out = apply_antigen_patch(records([("AAA", "curated-gene", "curated-species")]),
                              load_antigen_patch(p))
    assert out["antigen.gene"].to_list() == ["curated-gene"]
    assert out["antigen.species"].to_list() == ["SpeciesB"]


def test_an_epitope_the_patch_does_not_cover_is_untouched(tmp_path):
    p = write_patch(tmp_path, "AAA\tgene\tspecies\n")
    out = apply_antigen_patch(records([("BBB", "own-gene", "own-species")]),
                              load_antigen_patch(p))
    assert out["antigen.gene"].to_list() == ["own-gene"]
    assert out["antigen.species"].to_list() == ["own-species"]


def test_apply_leaves_no_working_columns_and_keeps_the_row_count(tmp_path):
    """A left join that fanned out would silently multiply records."""
    p = write_patch(tmp_path, "AAA\tgene\tspecies\nAAA\tgene2\tspecies2\n")
    df = records([("AAA", "x", "y"), ("AAA", "x", "y"), ("BBB", "x", "y")])
    out = apply_antigen_patch(df, load_antigen_patch(p))
    assert out.height == df.height
    assert [c for c in out.columns if c.startswith("__")] == []


def test_conflicts_finds_an_epitope_answered_two_ways(tmp_path):
    p = write_patch(tmp_path, "AAA\tGag\tHIV-1\nAAA\tPol\tHIV-1\nBBB\tM\tInfluenzaA\n")
    got = conflicts(p)
    assert got["antigen.epitope"].to_list() == ["AAA"]
    assert got["answers"][0].to_list() == ["HIV-1 Gag", "HIV-1 Pol"]


def test_conflicts_is_empty_when_every_epitope_has_one_answer(tmp_path):
    p = write_patch(tmp_path, "AAA\tGag\tHIV-1\nBBB\tM\tInfluenzaA\n")
    assert conflicts(p).height == 0


def test_the_committed_patch_conflicts_are_the_declared_set():
    """The named set may shrink as a curator settles one; it must not grow unnoticed."""
    got = conflicts()
    assert set(got["antigen.epitope"]) <= set(CONFLICTING_EPITOPES), (
        "a new epitope is answered two ways; settle it or add it to CONFLICTING_EPITOPES")
    for epitope, answers in got.iter_rows():
        assert tuple(answers) == CONFLICTING_EPITOPES[epitope]


def test_renames_are_scoped_to_slim_and_keyed_on_equality():
    """A 9-mer is a substring of a 10-mer, so a contains predicate rewrites the wrong rows."""
    reference = records([("SPRWYFYYL", "S", "SARS-CoV-2")])
    patched = records([("SPRWYFYYL", "Spike", "SARS-CoV-2")])
    out = render_patch_renames(reference, patched)
    assert 'when_equals = "SPRWYFYYL"' in out
    assert "when_contains" not in out
    assert 'files = "vdjdb.slim.txt"' in out
    assert 'from = "S"' in out and 'to = "Spike"' in out


def test_no_rename_is_emitted_when_nothing_changed():
    same = records([("AAA", "gene", "species")])
    assert render_patch_renames(same, same) == ""


def test_an_epitope_the_candidate_answers_two_ways_yields_no_rename():
    """The rename needs one target; several would make slim's grouping ambiguous."""
    reference = records([("AAA", "old", "species")])
    patched = records([("AAA", "one", "species"), ("AAA", "two", "species")])
    assert render_patch_renames(reference, patched) == ""


@pytest.mark.parametrize("epitope", sorted(CONFLICTING_EPITOPES))
def test_every_declared_conflict_is_still_in_the_patch_file(epitope):
    """A name left behind after the file changed is a comment that has stopped being true."""
    patch = pl.read_csv("patches/antigen_epitope_species_gene.dict", separator="\t",
                        quote_char=None, infer_schema_length=0)
    assert epitope in set(patch["antigen.epitope"])
