"""``method.frequency`` as three independent columns (#696).

The column held at least three different things - a ratio, a bare float and a percentage - and the
one thing the confidence score wants, the read count, was reachable in only two thirds of the
non-blank cells and unreachable in all of them without parsing free text.
"""
from __future__ import annotations

import polars as pl
import pytest

from vdjdb.curate.frequency import discordant, split_frequency

COUNT, TOTAL = "method.frequency.count", "method.frequency.total"


def _frame(values: list[str], **extra: list[str]) -> pl.DataFrame:
    return pl.DataFrame({"method.frequency": values, **extra})


@pytest.mark.parametrize(("submitted", "count", "total"), [
    ("7/30", 7, 30),
    ("5719/33921", 5719, 33921),
    ("12 / 163", 12, 163),
    # 297 records over 2 chunks double the slash. A doubled separator carries no second meaning.
    ("1//13", 1, 13),
    ("0/62", 0, 62),
])
def test_an_unambiguous_ratio_gives_a_count_and_a_total(submitted, count, total):
    got, _ = split_frequency(_frame([submitted]))
    assert (got[COUNT][0], got[TOTAL][0]) == (count, total)


@pytest.mark.parametrize("submitted", [
    "", "0.052115583", "0.163%", "18.06%", "99%", "1e-04", "2e-05",
    # Not parsed on purpose: a partial read of a value nobody checked is worse than no read.
    "1/2/3", "12/163 of 200", "~7/30", "7/30 (23%)",
])
def test_anything_else_carries_no_count(submitted):
    got, _ = split_frequency(_frame([submitted]))
    assert got[COUNT][0] is None and got[TOTAL][0] is None


def test_the_submitted_string_is_never_rewritten():
    """The parse is derived; `method.frequency` is what the curator typed, including the typo."""
    values = ["1//13", "7/30", "0.163%"]
    got, _ = split_frequency(_frame(values))
    assert got["method.frequency"].to_list() == values


def test_a_submitted_count_wins_over_a_parsed_one():
    """The two are chunk columns, so a curator may write them - and then nothing is parsed.

    Always trust the submission first: `3/40` beside a submitted count of 5 keeps the 5 and reports
    the disagreement, rather than quietly replacing it with what the string implies.
    """
    got, report = split_frequency(_frame(["3/40"], **{COUNT: ["5"], TOTAL: ["40"]}))
    assert (got[COUNT][0], got[TOTAL][0]) == (5, 40)
    assert report.filter(pl.col("shape") == "submitted count/total")["rows"][0] == 1


def test_a_missing_half_of_a_submitted_pair_is_filled_from_the_string():
    got, _ = split_frequency(_frame(["3/40"], **{COUNT: ["5"], TOTAL: [None]}))
    assert (got[COUNT][0], got[TOTAL][0]) == (5, 40)


def test_the_report_names_the_shape_of_every_cell():
    _, report = split_frequency(_frame(["7/30", "0.163%", "0.05", "", "1e-04"]))
    assert dict(zip(report["shape"], report["rows"], strict=True)) == {
        "parsed count/total": 1, "percentage": 1, "float": 2, "blank": 1}


def test_a_frequency_its_own_count_contradicts_is_reported():
    """Live only because the count and total are chunk columns: a `method.frequency` cell is one
    shape or the other, so before #696 no record could report all three."""
    # 0.1005 is inside the 1 % tolerance and 0.50 is five times the ratio. A value *at* the
    # tolerance is deliberately not tested: which side of it float arithmetic lands on is not a
    # property worth pinning.
    frame = _frame(["0.50", "0.10", "0.1005", "10%", "0.163%"],
                   **{COUNT: ["1", "1", "1", "1", None], TOTAL: ["10", "10", "10", "10", None]})
    assert frame.select(discordant())["method.frequency"].to_list() == [
        True, False, False, False, False]


def test_a_zero_total_is_not_a_division():
    """`0/0` is meaningless and must not raise or report."""
    frame = _frame(["0.5"], **{COUNT: ["0"], TOTAL: ["0"]})
    assert frame.select(discordant())["method.frequency"].to_list() == [False]


def test_the_check_runs_on_an_all_string_frame():
    """`vdjdb qc` reads chunks with every cell a string, where a division would raise."""
    frame = pl.DataFrame({"method.frequency": ["0.50"], COUNT: ["1"], TOTAL: ["10"]},
                         schema={"method.frequency": pl.String, COUNT: pl.String, TOTAL: pl.String})
    assert frame.select(discordant())["method.frequency"].to_list() == [True]


def test_a_frame_without_the_column_is_returned_unchanged():
    frame = pl.DataFrame({"record_id": ["VDJDB1"]})
    got, report = split_frequency(frame)
    assert got.equals(frame) and report.is_empty()
