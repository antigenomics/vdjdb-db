"""IMGT segment nomenclature (#389): mechanical respellings only, never a guess."""
from __future__ import annotations

import polars as pl
import pytest

from vdjdb.curate import nomenclature as N


@pytest.mark.parametrize("call,expected", [
    # the /DV genes IMGT names once for both loci -- 1,377 chains write the short form
    ("TRAV14", "TRAV14/DV4"),
    ("TRAV29", "TRAV29/DV5"),
    ("TRAV29DV5", "TRAV29/DV5"),        # the slash dropped entirely
    ("TRBJ1.2", "TRBJ1-2"),             # a dot where a dash belongs
    ("TRBJ 2-7", "TRBJ2-7"),            # a space
    ("TRBD1-1*01", "TRBD1*01"),         # a D gene with a suffix IMGT does not use
    ("TCRBD2*02", "TRBD2*02"),          # the pre-IMGT TCR prefix
])
def test_mechanical_respellings(call, expected):
    assert N.normalise_call(call, "HomoSapiens") == expected


@pytest.mark.parametrize("call", [
    "TRBV10-3*01",      # already IMGT
    "TRBV8",            # several candidates: TRBV8-1, TRBV8-2 -- a curation question, not a fix
    "TRAJ16.5",         # no IMGT counterpart at all
    "",
])
def test_a_call_that_is_correct_or_undecidable_is_left_alone(call):
    assert N.normalise_call(call, "HomoSapiens") is None


@pytest.mark.parametrize("species,call,expected", [
    # IMGT names the gene `TRBV19`; nothing is named `TRBV19-1`, and the family has no second member
    ("HomoSapiens", "TRBV19-1", "TRBV19"),
    ("HomoSapiens", "TRBV19-1*01", "TRBV19*01"),
    ("HomoSapiens", "TRBV28-1", "TRBV28"),
    ("MacacaMulatta", "TRBV19-1", "TRBV19"),
    ("HomoSapiens", "TRAJ24-1", "TRAJ24"),
    ("HomoSapiens", "TRAV14-1", "TRAV14/DV4"),   # the suffix goes, then the /DV gene is named
])
def test_a_spurious_minus_one_is_dropped(species, call, expected):
    assert N.normalise_call(call, species) == expected


@pytest.mark.parametrize("species,call", [
    ("HomoSapiens", "TRBV6-1"),      # a real IMGT gene: TRBV6-2 exists, so -1 is not spurious
    ("HomoSapiens", "TRBV20-1"),     # a real IMGT gene, and the only member of its family
    ("HomoSapiens", "TRAJ37-2"),     # not a -1: dropping it would discard the curator's distinction
    ("HomoSapiens", "TRBV13-6"),     # likewise -- `TRBV13` exists, so a blind drop would "resolve"
    ("MusMusculus", "TRAV13-1"),     # real in mouse, where the family runs to TRAV13-4/DV7
    ("MusMusculus", "TRAV14D-1"),
])
def test_a_suffix_that_is_not_a_spurious_minus_one_survives(species, call):
    assert N.normalise_call(call, species) is None


def test_a_species_with_no_authority_is_never_touched():
    """An unchecked rewrite is worse than an odd spelling."""
    assert N.normalise_call("TRAV14", "GallusGallus") is None


def test_a_multi_call_is_normalised_member_by_member_and_sorted():
    """`TRBD2,TRBD1` and `TRBD1,TRBD2` are two spellings of one fact -- 89 chains' worth."""
    assert N.normalise_call("TRBD2,TRBD1", "HomoSapiens") == "TRBD1,TRBD2"
    assert N.normalise_call("TRBD2*01 or TCRBD2*02", "HomoSapiens") == "TRBD2*01,TRBD2*02"


def test_a_multi_call_with_one_unresolvable_member_is_refused_whole():
    """Half a correction is not a correction."""
    assert N.normalise_call("TRBD1,TRBV8", "HomoSapiens") is None


def test_the_frame_is_rewritten_and_the_report_accounts_for_every_change():
    df = pl.DataFrame({
        "species": ["HomoSapiens", "HomoSapiens", "MusMusculus"],
        "v.alpha": ["TRAV14", "TRAV14", "TRAV21-DV12"],
        "j.alpha": ["TRAJ3*01", "TRAJ3*01", "TRAJ3*01"],
        "v.beta": ["", "", ""], "d.beta": ["", "", ""], "j.beta": ["", "", ""],
    })
    out, report = N.harmonise_segments(df)
    assert out["v.alpha"].to_list() == ["TRAV14/DV4", "TRAV14/DV4", "TRAV21/DV12"]
    assert report["rows"].sum() == 3
    assert set(report["from"]) == {"TRAV14", "TRAV21-DV12"}
    assert report.filter(pl.col("from") == "TRAV14")["rows"][0] == 2


def test_nothing_to_do_yields_an_empty_report_not_an_error():
    df = pl.DataFrame({"species": ["HomoSapiens"], "v.alpha": ["TRAV12-2*01"],
                       "j.alpha": [""], "v.beta": [""], "d.beta": [""], "j.beta": [""]})
    out, report = N.harmonise_segments(df)
    assert report.is_empty() and out["v.alpha"][0] == "TRAV12-2*01"


# -- the generated rename block ----------------------------------------------------------------

def test_a_rename_is_declared_only_when_the_fixer_leaves_the_old_spelling_alone():
    """The injectivity condition, and the reason the comparison can trust the block.

    `get_closest_id` simplifies `TRAV6-7-DV9` to `TRAV6` and then tries `TRAV6-1*01`, `TRAV6-2*01`,
    ... taking the first hit -- so the reference ships `TRAV6-1*01` and is indistinguishable from
    records that genuinely are TRAV6-1. Declaring that rename rewrites both.
    """
    report = pl.DataFrame({"column": ["v.alpha", "v.alpha"], "species": ["MusMusculus"] * 2,
                           "from": ["TRAV6-7-DV9", "TRAV14"], "to": ["TRAV6-7/DV9", "TRAV14/DV4"],
                           "rows": pl.Series([15, 2], dtype=pl.UInt32)})
    resolve = lambda sp, call: {                     # noqa: E731
        "TRAV6-7-DV9": "TRAV6-1*01", "TRAV6-7/DV9": "TRAV6-7/DV9*01",
        "TRAV14": "TRAV14", "TRAV14/DV4": "TRAV14/DV4*01"}.get(call, call)
    block = N.render_renames(report, resolve)
    assert "TRAV6-7-DV9" not in block and "TRAV6-1" not in block
    assert 'from = "TRAV14"' in block and 'to = "TRAV14/DV4*01"' in block


def test_a_multi_call_is_never_declared_as_a_rename():
    """`fix_both` splits it and keeps the best member, so what ships is a selection, not a rename."""
    report = pl.DataFrame({"column": ["d.beta"], "species": ["HomoSapiens"],
                           "from": ["TRBD2,TRBD1"], "to": ["TRBD1,TRBD2"],
                           "rows": pl.Series([43], dtype=pl.UInt32)})
    assert "[[rename]]" not in N.render_renames(report, lambda sp, c: c)


def test_the_generated_block_is_replaced_in_place_not_appended(tmp_path):
    path = tmp_path / "rules.toml"
    path.write_text('[[rule]]\nid = "keep-me"\n')
    report = pl.DataFrame({"column": ["v.alpha"], "species": ["HomoSapiens"],
                           "from": ["TRAV14"], "to": ["TRAV14/DV4"],
                           "rows": pl.Series([2], dtype=pl.UInt32)})
    N.write_renames(report, path, lambda sp, c: c)
    N.write_renames(report, path, lambda sp, c: c)
    text = path.read_text()
    assert text.count("[[rename]]") == 1 and "keep-me" in text


# -- allele disambiguation from the CDR3 (#327) --------------------------------------------------

def _traj24(calls, cdr3s, species=None):
    n = len(calls)
    return pl.DataFrame({
        "species": species or ["HomoSapiens"] * n,
        "j.alpha": calls, "cdr3.alpha": cdr3s,
        "v.alpha": [""] * n, "v.beta": [""] * n, "d.beta": [""] * n, "j.beta": [""] * n,
    })


def test_the_cdr3_decides_the_allele_whatever_the_submitter_wrote():
    """73 of the 111 records explicitly called *01 carry the *02 signature; none carry *01's."""
    df = _traj24(["TRAJ24", "TRAJ24*01", "TRAJ24*02"], ["CAWGKLQF"] * 3)
    out, report = N.disambiguate_alleles(df)
    assert out["j.alpha"].to_list() == ["TRAJ24*02"] * 3
    assert report["rows"].sum() == 2          # the one already correct is not a change


def test_the_other_allele_signature_is_honoured_symmetrically():
    """`WGKFEF` appears zero times in the corpus today; the rule must still be the right one."""
    out, _ = N.disambiguate_alleles(_traj24(["TRAJ24*02"], ["CAWGKFEF"]))
    assert out["j.alpha"].to_list() == ["TRAJ24*01"]


def test_a_cdr3_with_no_signature_is_left_alone():
    """364 records have a CDR3 trimmed short of the anchor. No evidence, no correction."""
    out, report = N.disambiguate_alleles(_traj24(["TRAJ24", "TRAJ24*01"], ["CAVSDLE", "CAVSDLE"]))
    assert out["j.alpha"].to_list() == ["TRAJ24", "TRAJ24*01"]
    assert report.is_empty()


def test_a_cdr3_carrying_both_signatures_contradicts_itself_and_is_refused():
    out, _ = N.disambiguate_alleles(_traj24(["TRAJ24"], ["CAWGKLQFWGKFEF"]))
    assert out["j.alpha"].to_list() == ["TRAJ24"]


def test_another_species_is_not_touched_by_a_human_rule():
    out, _ = N.disambiguate_alleles(_traj24(["TRAJ24"], ["CAWGKLQF"], ["MusMusculus"]))
    assert out["j.alpha"].to_list() == ["TRAJ24"]


def test_a_different_gene_is_not_touched():
    out, _ = N.disambiguate_alleles(_traj24(["TRAJ42*01"], ["CAWGKLQF"]))
    assert out["j.alpha"].to_list() == ["TRAJ42*01"]


def test_the_allele_rename_carries_its_evidence_into_the_rules():
    """Not injective on value alone: the fixer resolves a bare `TRAJ24` to `*01`, so the reference
    ships the same value for the records the CDR3 corrects and the ones it does not."""
    report = pl.DataFrame({"issue": ["#327"] * 2, "column": ["j.alpha"] * 2,
                           "species": ["HomoSapiens"] * 2, "from": ["TRAJ24", "TRAJ24*01"],
                           "to": ["TRAJ24*02"] * 2, "signature": ["WGKLQF"] * 2,
                           "rows": pl.Series([974, 73], dtype=pl.UInt32)})
    block = N.render_allele_renames(report, lambda sp, c: "TRAJ24*01" if c == "TRAJ24" else c)
    assert block.count("[[rename]]") == 1, "both rows resolve to one reference value"
    assert 'from = "TRAJ24*01"' in block and 'to = "TRAJ24*02"' in block
    assert 'when_contains = "WGKLQF"' in block
    assert "records = 1047" in block


# -- MHC ----------------------------------------------------------------------------------------

def _mhc(a, b, species=None):
    n = len(a)
    return pl.DataFrame({"species": species or ["MusMusculus"] * n, "mhc.a": a, "mhc.b": b})


def test_murine_class_two_spellings_collapse_onto_one_molecule():
    """Three spellings of I-A(b) split one molecule's records three ways in the motif grouping.

    They collapse onto `H2-IAb`, the MGI-prefixed form the other fourteen murine strings already
    use, not onto the `I-Ab` that happened to carry the most records.
    """
    out, report = N.harmonise_mhc(_mhc(["H2-IAb", "I-Ab"], ["H2-IAb", "I-Ab"]))
    assert out["mhc.a"].to_list() == ["H2-IAb", "H2-IAb"]
    assert out["mhc.b"].to_list() == ["H2-IAb", "H2-IAb"]
    assert report.filter(pl.col("issue") == "mhc.dict")["rows"].sum() == 2


@pytest.mark.parametrize("old,new", [("H-2Aa", "H2-Aa"), ("H-2Eb1", "H2-Eb1"),
                                     ("H2-Ag7", "H2-IAg7"), ("H2-Ed", "H2-IEd")])
def test_the_other_murine_respellings(old, new):
    out, _ = N.harmonise_mhc(_mhc([""], [old]))
    assert out["mhc.b"].to_list() == [new]


def test_an_allele_absent_from_imgt_is_corrected_to_the_one_the_paper_reported():
    """#467: 0 rows for `A*24:01` in IPD-IMGT/HLA against 342 for `A*24:02`, and all 80 records come
    from one reference that reports testing in `*24:02`."""
    out, report = N.harmonise_mhc(_mhc(["HLA-A*24:01"], ["B2M"], ["HomoSapiens"]))
    assert out["mhc.a"].to_list() == ["HLA-A*24:02"]
    assert report.filter(pl.col("issue") == "mhc.dict")["rows"].sum() == 1


def test_the_class_two_chains_are_put_in_order():
    """`mhc.a` is the first chain. The gene symbol says which chain it is, so this needs no
    judgement -- 149 records had them the wrong way round."""
    out, report = N.harmonise_mhc(
        _mhc(["HLA-DRB1*01:01"], ["HLA-DRA*01:01"], ["HomoSapiens"]))
    assert out["mhc.a"].to_list() == ["HLA-DRA*01:01"]
    assert out["mhc.b"].to_list() == ["HLA-DRB1*01:01"]
    assert report.filter(pl.col("issue") == "mhc-chain-order")["rows"].sum() == 1


def test_a_pair_already_in_order_is_not_swapped_back():
    out, report = N.harmonise_mhc(_mhc(["HLA-DRA*01:01"], ["HLA-DRB1*01:01"], ["HomoSapiens"]))
    assert out["mhc.a"].to_list() == ["HLA-DRA*01:01"]
    assert report.filter(pl.col("issue") == "mhc-chain-order").is_empty()


def test_a_class_one_pair_is_left_alone():
    out, report = N.harmonise_mhc(_mhc(["HLA-A*02:01"], ["B2M"], ["HomoSapiens"]))
    assert out.row(0) == ("HomoSapiens", "HLA-A*02:01", "B2M")
    assert report.is_empty()


def test_the_chain_swap_is_not_declared_as_a_rename():
    """It rewrites two columns at once; the comparison takes it as a declared row delta."""
    report = pl.DataFrame({"issue": ["mhc-chain-order"], "column": ["mhc.a,mhc.b"],
                           "from": ["beta,alpha"], "to": ["alpha,beta"], "rows": [149]})
    assert N.render_mhc_renames(report) == ""


# -- references (#347) ---------------------------------------------------------------------------

def test_a_reference_with_a_pubmed_id_gets_it(tmp_path):
    root = tmp_path
    (root / "proofreading").mkdir()
    (root / "proofreading" / "reference_ids.tsv").write_text(
        "# a comment\nreference.id\tpmid\tnote\ndoi:10.1\tPMID:1\tx\n")
    df = pl.DataFrame({"reference.id": ["doi:10.1", "PMID:9", "https://www.10xgenomics.com/x"]})
    out, report = N.harmonise_references(df, root)
    assert out["reference.id"].to_list() == ["PMID:1", "PMID:9", "https://www.10xgenomics.com/x"]
    assert report["rows"].to_list() == [1]


def test_a_reference_with_no_pubmed_id_is_left_alone(tmp_path):
    """Most of #347's 30,977 non-PMID records cannot be resolved and should not be: a 10x
    application note, a PDB entry and a direct submission's own issue are not papers."""
    (tmp_path / "proofreading").mkdir()
    (tmp_path / "proofreading" / "reference_ids.tsv").write_text("reference.id\tpmid\tnote\n")
    df = pl.DataFrame({"reference.id": ["https://github.com/antigenomics/vdjdb-db/issues/193"]})
    out, report = N.harmonise_references(df, tmp_path)
    assert out["reference.id"][0].endswith("/193")
    assert report.is_empty()


def test_the_committed_table_resolves_what_it_claims():
    """The table is a committed input, so its content is part of the build's correctness."""
    from vdjdb.config import Paths

    table = pl.read_csv(Paths.discover().root / "proofreading" / "reference_ids.tsv",
                        separator="\t", infer_schema=False, comment_prefix="#")
    assert table.height >= 3
    assert all(p.startswith("PMID:") and p[5:].isdigit() for p in table["pmid"])
    assert table["reference.id"].n_unique() == table.height


def test_a_correction_scoped_to_one_reference_does_not_leak(tmp_path):
    """`HLA-A*24:09` is corrected only for PMID:39286976, whose title states A*24:02. Elsewhere it
    is left alone, because nothing has been checked about it there."""
    (tmp_path / "patches").mkdir()
    (tmp_path / "patches" / "mhc.dict").write_text(
        "mhc\treplacement\treference.id\tnote\n"
        "HLA-A*24:09\tHLA-A*24:02\tPMID:39286976\tthe paper says *24:02\n")
    df = pl.DataFrame({"species": ["HomoSapiens"] * 2, "mhc.a": ["HLA-A*24:09"] * 2,
                       "mhc.b": ["B2M"] * 2,
                       "reference.id": ["PMID:39286976", "PMID:99999999"]})
    out, report = N.harmonise_mhc(df, tmp_path)
    assert out["mhc.a"].to_list() == ["HLA-A*24:02", "HLA-A*24:09"]
    assert report["rows"].to_list() == [1]


def test_the_shipped_mhc_patch_is_well_formed():
    """A committed patch is part of the build's correctness. Four entries in the antigen patch had a
    space where the tab belongs and had therefore never matched a single record -- including
    `ALAGIGILTV` (MLANA), one of the most-studied human epitopes."""
    from vdjdb.config import Paths

    root = Paths.discover().root
    for name, n_cols in (("mhc.dict", 4), ("antigen_epitope_species_gene.dict", 3)):
        path = root / "patches" / name
        body = [ln for ln in path.read_text().splitlines()[1:]
                if ln and not ln.startswith("#")]
        bad = [ln for ln in body if ln.count("\t") != n_cols - 1]
        assert not bad, f"{name}: {bad[:3]}"
        keys = [ln.split("\t")[0] for ln in body]
        dupes = {k for k in keys if keys.count(k) > 1}
        if name == "antigen_epitope_species_gene.dict":
            from vdjdb.curate.patch import CONFLICTING_EPITOPES, conflicts
            # A duplicate that gives the *same* answer twice is harmless; one that gives two is
            # resolved by file order, silently. Those are declared, so the set cannot grow.
            assert set(conflicts()["antigen.epitope"]) == set(CONFLICTING_EPITOPES)
            assert set(CONFLICTING_EPITOPES) <= dupes


# --- #136 and #402: the separator, the roman numeral, and the one-digit allele ------------------

@pytest.mark.parametrize("call,want", [
    # `/` read as "or", in the three shapes a curator writes it
    ("TRBV12-3/TRBV12-4", "TRBV12-3,TRBV12-4"),
    ("TRBV6-2/TRBV6-3", "TRBV6-2,TRBV6-3"),
    ("TRBV19*01/02", "TRBV19*01,TRBV19*02"),
    ("TRBV20-1*01/04/05", "TRBV20-1*01,TRBV20-1*04,TRBV20-1*05"),
    ("TRBV12-3/4*01", "TRBV12-3*01,TRBV12-4*01"),
    # a one-digit allele: IMGT writes two
    ("TRBJ2-7*1", "TRBJ2-7*01"),
    ("TRAV20*2", "TRAV20*02"),
    ("TRBV7-9*3", "TRBV7-9*03"),
    # the Arden roman numeral, which the Arden table then maps
    ("TRBVIS1", "TRBV9"),
])
def test_a_separator_or_a_mistyped_digit_is_a_spelling_and_is_fixed(call, want):
    assert N.normalise_call(call, "HomoSapiens") == want


@pytest.mark.parametrize("call", [
    # `/` that is part of the gene name, not a separator: IMGT names this gene once, for two loci
    "TRAV14/DV4",
    "TRAV14/DV4*01",
    "TRAV23/DV6",
    # already IMGT
    "TRBV27",
    "TRBJ2-7*01",
])
def test_a_slash_inside_an_imgt_gene_name_is_left_alone(call):
    assert N.normalise_call(call, "HomoSapiens") is None


def test_an_ambiguous_slash_is_a_curation_decision_and_is_refused():
    """`TRBV11/2` could be `TRBV11-2` or two genes. Guessing which is not a spelling fix."""
    assert N.normalise_call("TRBV11/2", "HomoSapiens") is None


def test_the_mouse_traj47_rule_reads_the_anchor_residue_itself():
    """#647. *01 templates `HYANKMIC` and is ORF; *02 templates `DYANKMIF`. All 95 mouse chains
    read the Phe, so all 95 are *02 -- 81 written with no allele and 14 naming *01."""
    out, report = N.disambiguate_alleles(_traj24(["TRAJ47", "TRAJ47*01"], ["CAANKMIF"] * 2,
                                                 ["MusMusculus"] * 2))
    assert out["j.alpha"].to_list() == ["TRAJ47*02"] * 2
    assert report["rows"].sum() == 2


def test_human_traj47_is_out_of_reach_of_the_mouse_rule():
    """Both human alleles template `EYGNKLVF`, so the sequence cannot decide and must not try. The
    three human `TRAJ47-1` spellings are a name IMGT has not got (#389) and the species scope is
    what keeps them from being rewritten into a real gene here."""
    calls = ["TRAJ47", "TRAJ47*01", "TRAJ47-1*01"]
    out, report = N.disambiguate_alleles(_traj24(calls, ["CAANKMIF"] * 3))
    assert out["j.alpha"].to_list() == calls
    assert report.is_empty()
