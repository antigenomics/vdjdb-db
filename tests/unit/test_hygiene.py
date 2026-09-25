"""The reader strips surrounding whitespace, because it forks one value into two."""
from __future__ import annotations

import polars as pl

from vdjdb.io.chunks import read_chunk
from vdjdb.schema import ALL_COLUMNS

#: A real J-gene call in the corpus carried one of these. Named, because a literal is unreadable.
NBSP = chr(0xA0)


def test_surrounding_whitespace_is_stripped(tmp_path):
    """`tetramer-sort ` beside `tetramer-sort` is 103 records' worth of one value looking like two."""
    header = "\t".join(ALL_COLUMNS)
    row = "\t".join("tetramer-sort " if c == "method.identification" else
                    " HLA-DRB1*15 " if c == "meta.donor.MHC" else
                    NBSP + "TRAJ12*01" if c == "j.alpha" else "x" for c in ALL_COLUMNS)
    p = tmp_path / "PMID_1.txt"
    p.write_text(f"{header}\n{row}\n")
    df = read_chunk(p)
    assert df["method.identification"][0] == "tetramer-sort"
    assert df["meta.donor.MHC"][0] == "HLA-DRB1*15"
    assert df["j.alpha"][0] == "TRAJ12*01", "a non-breaking space counts as whitespace"


def test_an_all_whitespace_cell_becomes_the_missing_marker(tmp_path):
    header = "\t".join(ALL_COLUMNS)
    row = "\t".join("   " if c == "antigen.gene" else "x" for c in ALL_COLUMNS)
    p = tmp_path / "PMID_2.txt"
    p.write_text(f"{header}\n{row}\n")
    assert read_chunk(p)["antigen.gene"][0] == ""


def test_the_identical_chain_rule_is_advisory_and_finds_the_copy(tmp_path):
    """#561: the beta CDR3 copied into the alpha field, V and J left correct. 99 records."""
    from vdjdb.qc.rules import check
    from vdjdb.qc.runner import ADVISORY

    assert "alpha and beta cdr3 identical" in ADVISORY
    df = pl.DataFrame({c: ["" for _ in range(3)] for c in ALL_COLUMNS}).with_columns(
        pl.Series("cdr3.alpha", ["CASSIRSSYEQYF", "CAVSDLEPNSSASKIIF", ""]),
        pl.Series("cdr3.beta", ["CASSIRSSYEQYF", "CASSIRSSYEQYF", "CASSIRSSYEQYF"]),
        pl.Series("chunk.file", ["a", "a", "a"]),
        pl.Series("chunk.row", [1, 2, 3]),
    )
    found = check(df).filter(pl.col("rule") == "alpha and beta cdr3 identical")
    assert found.height == 1 and found["chunk.row"].to_list() == [1]
