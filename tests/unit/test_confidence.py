"""The VDJdb confidence score.

`score/confidence.py` was 39 % covered while deciding the number every consumer filters on. These
pin the two behaviours the released scores depend on - the `n/m` frequency parse where `n` is a
cell count, and the maximum over the sample signature - plus each branch of the two sub-scores.
"""
from __future__ import annotations

import polars as pl
import pytest

from vdjdb.score.confidence import (
    SCORE_SIGNATURE,
    add_score,
    cell_count,
    frequency,
    sequencing_score,
)

#: Every column `add_score` touches, with values that score 0 unless a test overrides them.
BLANK = {
    "method.frequency": "", "method.singlecell": "", "method.sequencing": "",
    "method.identification": "", "method.verification": "", "meta.structure.id": "",
    **{c: "" for c in SCORE_SIGNATURE},
}


def row(**over):
    return pl.DataFrame([{**BLANK, **over}])


def score(**over):
    return add_score(row(**over))["vdjdb.score"][0]


@pytest.mark.parametrize("cell, expected", [
    ("", 0.0),
    ("0.5", 0.5),
    ("50%", 0.5),
    ("2/47", 2 / 47),
    ("2//47", 2 / 47),        # occurs in the data
    ("nonsense", 0.0),
    ("1/0", 0.0),             # division by zero must not produce inf or raise
])
def test_frequency_parses_every_form_that_occurs(cell, expected):
    got = row(**{"method.frequency": cell}).select(frequency().alias("f"))["f"][0]
    assert got == pytest.approx(expected), cell


@pytest.mark.parametrize("cell, expected", [
    ("2/47", 2),
    ("2//47", 2),
    ("50%", 0),               # a percentage carries no cell count
    ("0.5", 0),
    ("", 0),
])
def test_cell_count_is_the_numerator_only(cell, expected):
    got = row(**{"method.frequency": cell}).select(cell_count().alias("n"))["n"][0]
    assert got == expected, cell


def test_a_structure_outranks_everything():
    assert score(**{"meta.structure.id": "1AO7"}) == 3


def test_single_cell_gives_full_sequencing_confidence():
    """Sequencing confidence 3, but the score is the minimum with specificity, which is 0 here."""
    assert score(**{"method.singlecell": "yes", "method.verification": "direct"}) == 3
    assert score(**{"method.singlecell": "yes"}) == 0


def test_singlecell_no_is_not_single_cell():
    assert score(**{"method.singlecell": "no", "method.verification": "direct"}) == 3


@pytest.mark.parametrize("verification, expected", [
    ("direct", 3),
    ("antigen-loaded targets", 2),
    ("cell sorting", 1),
    ("", 0),
])
def test_verification_sets_the_specificity_ceiling(verification, expected):
    """With `direct` verification the sequence is trusted too, so the minimum is the ceiling."""
    got = score(**{"method.verification": verification, "method.singlecell": "yes"})
    assert got == expected, verification


@pytest.mark.parametrize("identification, freq, expected", [
    ("cell culture", "0.6", 1),     # culture is judged hardest: >= 0.5
    ("cell culture", "0.4", 0),
    ("cell sorting", "0.06", 1),    # sort-based: >= 0.05
    ("cell sorting", "0.04", 0),
    ("antigen-loaded targets", "0.3", 1),   # stimulation: >= 0.25
    ("antigen-loaded targets", "0.2", 0),
])
def test_the_moderate_point_needs_the_assay_threshold(identification, freq, expected):
    got = score(**{"method.identification": identification, "method.frequency": freq,
                   "method.singlecell": "yes"})
    assert got == expected, (identification, freq)


def test_amplicon_depth_decides_sequencing_confidence():
    high = score(**{"method.sequencing": "amplicon-seq", "method.frequency": "0.02",
                    "method.identification": "cell culture", "method.verification": "direct"})
    low = row(**{"method.sequencing": "amplicon-seq", "method.frequency": "0.001",
                 "method.identification": "", "method.verification": ""})
    assert high == 3
    assert add_score(low)["vdjdb.score"][0] == 0


@pytest.mark.parametrize("cell, count, expected", [
    ("2/100", 2, 3),       # 0.02 of the reads, 2 reads
    ("1/100", 1, 1),       # frequency clears 0.01, one read does not
    ("2/1000", 2, 1),      # 2 reads, frequency 0.002 is under 0.01
    ("0.02", None, 3),     # a float carries no count: the frequency decides alone
    ("2%", None, 3),       # so does a percentage
    ("0.002", None, 1),    # and a float under the threshold still fails it
    ("0.02", 1, 1),        # a submitted count is used when there is one, whatever the float says
    ("0.02", 5, 3),
])
def test_amplicon_needs_two_reads_only_when_the_count_is_known(cell, count, expected):
    df = row(**{"method.sequencing": "amplicon-seq", "method.frequency": cell}).with_columns(
        pl.lit(count, dtype=pl.Int64).alias("method.frequency.count"))
    got = df.select(sequencing_score(frequency(), cell_count(), pl.col("method.frequency.count")))
    assert got.item() == expected, (cell, count)


def test_the_count_column_reaches_the_score():
    """The branch is reached through `add_score`, not only through the helper."""
    base = {"method.sequencing": "amplicon-seq", "method.identification": "cell culture",
            "method.verification": "", "method.frequency": "0.9", "method.singlecell": ""}
    one = row(**base).with_columns(pl.lit(1, dtype=pl.Int64).alias("method.frequency.count"))
    many = row(**base).with_columns(pl.lit(9, dtype=pl.Int64).alias("method.frequency.count"))
    none = row(**base).with_columns(pl.lit(None, dtype=pl.Int64).alias("method.frequency.count"))
    assert [add_score(f)["vdjdb.score"][0] for f in (one, many, none)] == [1, 1, 1]
    # sequencing 3 against specificity 1 gives 1 either way, so look at the term itself
    got = [f.with_columns(frequency().alias("__freq"), cell_count().alias("__count"),
                          pl.col("method.frequency.count").alias("__reads"))
            .select(sequencing_score(pl.col("__freq"), pl.col("__count"), pl.col("__reads"))).item()
           for f in (one, many, none)]
    assert got == [1, 3, 3]


def test_sanger_needs_two_cells_for_full_confidence():
    two = score(**{"method.sequencing": "sanger", "method.frequency": "2/47",
                   "method.verification": "antigen-loaded targets"})
    one = score(**{"method.sequencing": "sanger", "method.frequency": "1/47",
                   "method.verification": "antigen-loaded targets"})
    assert two == 2 and one == 2, "verification raises sequencing confidence to 3 either way"


def test_the_score_is_a_maximum_over_the_signature():
    """The same clonotype assayed twice takes the better reading, in both rows."""
    shared = {c: "x" for c in SCORE_SIGNATURE}
    df = pl.DataFrame([
        {**BLANK, **shared, "meta.structure.id": "1AO7"},
        {**BLANK, **shared},
    ])
    assert add_score(df)["vdjdb.score"].to_list() == [3, 3]


def test_a_different_signature_does_not_share_a_score():
    df = pl.DataFrame([
        {**BLANK, **{c: "x" for c in SCORE_SIGNATURE}, "meta.structure.id": "1AO7"},
        {**BLANK, **{c: "y" for c in SCORE_SIGNATURE}},
    ])
    assert add_score(df)["vdjdb.score"].to_list() == [3, 0]


def test_add_score_leaves_no_working_columns_behind():
    out = add_score(row())
    assert [c for c in out.columns if c.startswith("__")] == []
    assert "vdjdb.score" in out.columns
