"""The promiscuity table: the contract of the committed file, not the predictor behind it.

`mhcmatch` fetches its reference data from HuggingFace, so nothing here scores anything -- a test that
fails when HuggingFace is slow is a test that gets disabled. What is asserted is the part that goes
wrong without anybody noticing: a reviewed input edited by hand, a writer that changes its column
order, or a sort that stops being deterministic and makes the file churn on every refresh.
"""
from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from vdjdb.curate import promiscuity as P

TABLE = Path(__file__).resolve().parents[2] / "proofreading" / "epitope_promiscuity.tsv"
BANDS = {"strong", "weak", "non-binder"}


@pytest.fixture(scope="module")
def table() -> pl.DataFrame:
    return pl.read_csv(TABLE, separator="\t", infer_schema=False)


def test_the_columns_are_the_ones_the_writer_declares(table: pl.DataFrame):
    assert tuple(table.columns) == P.COLUMNS


def test_every_row_is_a_class_i_ligand_length_and_a_known_band(table: pl.DataFrame):
    lengths = set(table["epitope"].str.len_chars().to_list())
    assert lengths <= set(P.LENGTHS), f"epitope lengths outside {P.LENGTHS}: {sorted(lengths)}"
    assert set(table["band"].to_list()) <= BANDS
    assert set(table["recorded"].to_list()) <= {"0", "1"}


def test_the_file_is_sorted_so_a_refresh_does_not_churn(table: pl.DataFrame):
    """Epitope, then ascending %rank, ties by allele. Hard rule 7: same inputs, same bytes."""
    typed = table.with_columns(pl.col("percent_rank").cast(pl.Float64))
    assert typed.equals(typed.sort("epitope", "percent_rank", "mhc_a"))


def test_a_non_binder_is_only_kept_when_vdjdb_records_the_pair(table: pl.DataFrame):
    """The filter the module documents: predictions are the presenting alleles, records are kept as
    they stand. A non-binder with ``recorded == 0`` would mean the cut stopped being applied."""
    assert table.filter((pl.col("band") == "non-binder") & (pl.col("recorded") == "0")).is_empty()


def test_the_table_answers_the_donor_filter_question(table: pl.DataFrame):
    """#372's use: given a donor allele, which epitopes could they present. The answer must be
    larger than the set VDJdb records under that allele, or the table adds nothing."""
    a2 = table.filter(pl.col("mhc_a") == "HLA-A*02:01")
    recorded = a2.filter(pl.col("recorded") == "1").height
    assert recorded > 0
    assert a2.height > recorded


def test_hpvtkyim_is_presented_by_b0801_which_is_what_resolved_its_record(table: pl.DataFrame):
    """The HCV NS3 8-mer `HPVTKYIM` is recorded in `chunks/` under `HLA-A*08:01`, an allele that does
    not exist at any resolution -- `HLA-A` has no `*08` group. Its only strong presenter in the
    107-allele panel is `HLA-B*08:01`, one locus letter away, which is the same defect class as #597
    (`HLA-A*80:01` for the `HLA-B*08:01` epitope `FLRGRAYGL`). Pinned here because it is the evidence
    the correction rests on.
    """
    rows = table.filter(pl.col("epitope") == "HPVTKYIM")
    strong = rows.filter(pl.col("band") == "strong")
    assert strong["mhc_a"].to_list() == ["HLA-B*08:01"]
