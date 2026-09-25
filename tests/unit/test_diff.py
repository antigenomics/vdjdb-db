"""The difference ledger, on fixtures small enough to reason about.

The ledger is the acceptance gate for every later phase, so its own failure modes matter more than
most: a ledger that under-reports would wave a regression through, and one that over-reports would
be switched off.
"""
from __future__ import annotations

import zipfile
from pathlib import Path

import pytest

from vdjdb.compare.diff import Bundle, Rule, diff, load_rules, render

HEADER = ("complex.id\tgene\tcdr3\tv.segm\tj.segm\tspecies\tmhc.a\tmhc.b\tmhc.class\t"
          "antigen.epitope\tantigen.gene\tantigen.species\treference.id\tvdjdb.score\tTCR_hash\t"
          "method\tmeta\tcdr3fix\tweb.method\tweb.method.seq\tweb.cdr3fix.nc\tweb.cdr3fix.unmp")


def _row(cdr3: str, *, ref: str = "PMID:1", unmp: str = "no", score: str = "0") -> str:
    return "\t".join(["0", "TRB", cdr3, "TRBV1*01", "TRBJ1*01", "HomoSapiens",
                      "HLA-A*02:01", "B2M", "MHCI", "GILGFVFTL", "M", "InfluenzaA", ref,
                      score, "h", "{}", "{}", "{}", "sort", "sanger", "no", unmp])


def _bundle(tmp: Path, name: str, rows: list[str]) -> Path:
    d = tmp / name
    d.mkdir()
    (d / "vdjdb.txt").write_text(HEADER + "\n" + "\n".join(rows) + "\n")
    return d


@pytest.fixture
def base() -> list[str]:
    return [_row("CASSA"), _row("CASSB"), _row("CASSC")]


# --------------------------------------------------------------------------------------------
# The zero case
# --------------------------------------------------------------------------------------------

def test_identical_bundles_diff_to_zero(tmp_path: Path, base: list[str]) -> None:
    r = diff(_bundle(tmp_path, "a", base), _bundle(tmp_path, "b", base))
    assert r.ok
    assert not r.unattributed
    assert all(f.raw_equal and f.canonical_equal for f in r.files)


def test_row_order_alone_breaks_raw_but_not_canonical(tmp_path: Path, base: list[str]) -> None:
    """The reproduction contract: canonical equality gates, raw equality only informs.

    ``runBuidDatabase.py`` iterates ``os.listdir``, so chunk order is host-dependent and raw
    equality is unrecoverable for existing releases.
    """
    r = diff(_bundle(tmp_path, "a", base), _bundle(tmp_path, "b", list(reversed(base))))
    f = r.files[0]
    assert not f.raw_equal
    assert f.canonical_equal
    assert r.ok, "a pure reordering must not fail the gate"


# --------------------------------------------------------------------------------------------
# Detection
# --------------------------------------------------------------------------------------------

def test_a_changed_cell_is_reported_with_both_values(tmp_path: Path, base: list[str]) -> None:
    changed = [_row("CASSA", unmp="yes"), *base[1:]]
    r = diff(_bundle(tmp_path, "a", base), _bundle(tmp_path, "b", changed))
    assert not r.ok
    assert len(r.unattributed) == 1
    cell = r.unattributed[0]
    assert (cell.column, cell.old, cell.new) == ("web.cdr3fix.unmp", "no", "yes")


def test_two_changed_columns_in_one_row_are_two_cells(tmp_path: Path, base: list[str]) -> None:
    r = diff(_bundle(tmp_path, "a", base),
             _bundle(tmp_path, "b", [_row("CASSA", unmp="yes", score="2"), *base[1:]]))
    assert {c.column for c in r.unattributed} == {"web.cdr3fix.unmp", "vdjdb.score"}
    assert r.files[0].changed_rows == 1, "one row changed, whatever the cell count"


def test_added_and_removed_rows_are_counted_separately(tmp_path: Path, base: list[str]) -> None:
    r = diff(_bundle(tmp_path, "a", base), _bundle(tmp_path, "b", [*base[:2], _row("CASSZ")]))
    f = r.files[0]
    assert (f.only_in_reference, f.only_in_candidate) == (1, 1)
    assert f.changed_rows == 0, "different keys are not a change, they are a removal plus an add"


def test_rows_sharing_a_key_are_compared_as_a_multiset(tmp_path: Path) -> None:
    """One paper reporting the same TCR in several donors is several rows with one key."""
    ref = [_row("CASSA", score="0"), _row("CASSA", score="1")]
    cand = [_row("CASSA", score="1"), _row("CASSA", score="0")]   # same multiset, swapped
    r = diff(_bundle(tmp_path, "a", ref), _bundle(tmp_path, "b", cand))
    assert r.files[0].changed_rows == 0
    assert r.ok


def test_one_of_two_rows_sharing_a_key_changing_is_one_change(tmp_path: Path) -> None:
    ref = [_row("CASSA", score="0"), _row("CASSA", score="1")]
    cand = [_row("CASSA", score="0"), _row("CASSA", score="2")]
    r = diff(_bundle(tmp_path, "a", ref), _bundle(tmp_path, "b", cand))
    assert r.files[0].changed_rows == 1
    assert [(c.old, c.new) for c in r.unattributed] == [("1", "2")]


def test_a_missing_file_fails(tmp_path: Path, base: list[str]) -> None:
    a = _bundle(tmp_path, "a", base)
    b = tmp_path / "b"
    b.mkdir()
    r = diff(a, b)
    assert r.missing == ["vdjdb.txt"]
    assert not r.ok


def test_a_reordered_column_is_caught_before_any_row_comparison(tmp_path: Path) -> None:
    """``vdjdb-web`` parses positionally, so a reordering silently mistypes the whole table."""
    a = _bundle(tmp_path, "a", [_row("CASSA")])
    b = tmp_path / "b"
    b.mkdir()
    cols = HEADER.split("\t")
    swapped = [*cols[:13], cols[14], cols[13], *cols[15:]]
    (b / "vdjdb.txt").write_text("\t".join(swapped) + "\n" + _row("CASSA") + "\n")
    r = diff(a, b)
    assert "reordered" in r.files[0].note
    assert not r.files[0].canonical_equal
    assert not r.ok


# --------------------------------------------------------------------------------------------
# Rule attribution -- both halves of the gate
# --------------------------------------------------------------------------------------------

def _rules(tmp: Path, body: str) -> Path:
    p = tmp / "rules.toml"
    p.write_text(body)
    return p


def test_a_declared_rule_with_the_right_count_passes(tmp_path: Path, base: list[str]) -> None:
    rules = _rules(tmp_path, '[[rule]]\nid = "unmp"\nfile = "vdjdb.txt"\n'
                             'column = "web.cdr3fix.unmp"\nfrom = "no"\nto = "yes"\nrows = 1\n')
    r = diff(_bundle(tmp_path, "a", base),
             _bundle(tmp_path, "b", [_row("CASSA", unmp="yes"), *base[1:]]), rules)
    assert r.rule_counts == {"unmp": 1}
    assert not r.unattributed
    assert r.ok


def test_a_rule_firing_the_wrong_number_of_times_fails(tmp_path: Path, base: list[str]) -> None:
    """The half that turns a rule from a description into a measurement."""
    rules = _rules(tmp_path, '[[rule]]\nid = "unmp"\ncolumn = "web.cdr3fix.unmp"\nrows = 7998\n')
    r = diff(_bundle(tmp_path, "a", base),
             _bundle(tmp_path, "b", [_row("CASSA", unmp="yes"), *base[1:]]), rules)
    assert r.miscounted_rules == ["unmp"]
    assert not r.ok, "a declared column is not a licence to change any number of cells in it"


def test_a_rule_that_never_fires_fails_too(tmp_path: Path, base: list[str]) -> None:
    rules = _rules(tmp_path, '[[rule]]\nid = "unmp"\ncolumn = "web.cdr3fix.unmp"\nrows = 1\n')
    r = diff(_bundle(tmp_path, "a", base), _bundle(tmp_path, "b", base), rules)
    assert r.miscounted_rules == ["unmp"]
    assert not r.ok, "a stale rule means the fix it describes stopped happening"


def test_a_rule_without_a_count_attributes_but_does_not_gate(tmp_path: Path, base: list[str]) -> None:
    rules = _rules(tmp_path, '[[rule]]\nid = "wip"\ncolumn = "web.cdr3fix.unmp"\n')
    r = diff(_bundle(tmp_path, "a", base),
             _bundle(tmp_path, "b", [_row("CASSA", unmp="yes"), *base[1:]]), rules)
    assert r.rule_counts == {"wip": 1}
    assert r.ok


def test_a_rule_does_not_cover_a_different_column(tmp_path: Path, base: list[str]) -> None:
    rules = _rules(tmp_path, '[[rule]]\nid = "unmp"\ncolumn = "web.cdr3fix.unmp"\n')
    r = diff(_bundle(tmp_path, "a", base),
             _bundle(tmp_path, "b", [_row("CASSA", score="3"), *base[1:]]), rules)
    assert r.unattributed and r.unattributed[0].column == "vdjdb.score"
    assert not r.ok


def test_from_and_to_narrow_a_rule(tmp_path: Path, base: list[str]) -> None:
    rules = _rules(tmp_path, '[[rule]]\nid = "unmp"\ncolumn = "web.cdr3fix.unmp"\n'
                             'from = "no"\nto = "yes"\n')
    r = diff(_bundle(tmp_path, "a", [_row("CASSA", unmp="yes"), *base[1:]]),
             _bundle(tmp_path, "b", base), rules)          # the reverse direction
    assert r.unattributed, "yes -> no is not the declared no -> yes"


def test_rules_file_absent_is_not_an_error(tmp_path: Path) -> None:
    assert load_rules(tmp_path / "nope.toml") == []


def test_first_matching_rule_wins(tmp_path: Path, base: list[str]) -> None:
    rules = _rules(tmp_path, '[[rule]]\nid = "specific"\ncolumn = "web.cdr3fix.unmp"\n'
                             'from = "no"\n\n[[rule]]\nid = "broad"\nfile = "vdjdb.txt"\n')
    r = diff(_bundle(tmp_path, "a", base),
             _bundle(tmp_path, "b", [_row("CASSA", unmp="yes"), *base[1:]]), rules)
    assert r.rule_counts == {"specific": 1}


# --------------------------------------------------------------------------------------------
# Bundles and rendering
# --------------------------------------------------------------------------------------------

def test_a_zip_and_a_directory_compare_equal(tmp_path: Path, base: list[str]) -> None:
    """The release ships a zip with a ``vdjdb-<date>/`` prefix; a build produces a flat directory."""
    d = _bundle(tmp_path, "a", base)
    z = tmp_path / "ref.zip"
    with zipfile.ZipFile(z, "w") as zf:
        zf.writestr("vdjdb-2026-06-03/vdjdb.txt", (d / "vdjdb.txt").read_text())
    assert Bundle(z).names == Bundle(d).names == ["vdjdb.txt"]
    assert diff(z, d).ok


def test_render_states_the_verdict(tmp_path: Path, base: list[str]) -> None:
    text = render(diff(_bundle(tmp_path, "a", base), _bundle(tmp_path, "b", base)))
    assert "**Verdict: PASS**" in text
    text = render(diff(_bundle(tmp_path, "a2", base),
                       _bundle(tmp_path, "b2", [_row("CASSA", score="3"), *base[1:]])))
    assert "**Verdict: FAIL**" in text
    assert "vdjdb.score" in text


def test_rule_matches_on_every_declared_facet() -> None:
    from vdjdb.compare.diff import CellDiff
    c = CellDiff("vdjdb.txt", "web.cdr3fix.unmp", "no", "yes", "k")
    assert Rule("r").matches(c)
    assert Rule("r", file="vdjdb.txt", column="web.cdr3fix.unmp", from_="no", to="yes").matches(c)
    assert not Rule("r", file="other.txt").matches(c)
    assert not Rule("r", to="maybe").matches(c)
