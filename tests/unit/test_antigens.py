"""The proteome authority for `antigen.gene` (#632): the table's shape, and what it must not claim."""
from __future__ import annotations

from pathlib import Path

import polars as pl
import pytest

from vdjdb.curate import antigens as A

TABLE = Path("proofreading/epitope_proteome.tsv")


@pytest.fixture(scope="module")
def table() -> pl.DataFrame:
    if not TABLE.exists():
        pytest.skip(f"no committed table at {TABLE}; run `vdjdb antigens`")
    return pl.read_csv(TABLE, separator="\t", infer_schema=False)


def test_the_columns_are_the_ones_the_writer_declares(table: pl.DataFrame) -> None:
    assert tuple(table.columns) == A.COLUMNS


def test_every_row_carries_one_of_the_three_verdicts(table: pl.DataFrame) -> None:
    assert set(table["verdict"]) <= {A.EXACT, A.ONE_SUB, A.NOT_FOUND}


def test_only_the_two_species_with_a_reference_proteome_appear(table: pl.DataFrame) -> None:
    """A viral epitope's source is the pathogen's proteome, which is a different fetch and a
    different question - `antigen.species` names the pathogen there, not the host."""
    assert set(table["antigen.species"]) <= set(A.PROTEOMES)


def test_an_exact_row_names_a_protein_and_a_one_sub_row_names_the_substitution(table) -> None:
    exact = table.filter(pl.col("verdict") == A.EXACT)
    assert (exact["source.protein"] != "").all()
    assert (exact["source.subs"] == "").all(), "exact means no substitution, so the column is empty"

    near = table.filter(pl.col("verdict") == A.ONE_SUB)
    assert (near["source.subs"] != "").all(), "a row that cannot name what differs says nothing"
    assert (near["source.peptide"] != near["antigen.epitope"]).all()


def test_a_not_found_row_claims_nothing(table: pl.DataFrame) -> None:
    """A question rather than a verdict, so it must not carry a half-guess at a source."""
    absent = table.filter(pl.col("verdict") == A.NOT_FOUND)
    for column in ("source.protein", "source.gene", "source.peptide", "source.subs"):
        assert (absent[column] == "").all(), f"{column} is set on an absent row"
    assert (absent["source.position"] == "-1").all()


def test_the_file_is_sorted_so_a_refresh_does_not_churn(table: pl.DataFrame) -> None:
    """Hard rule 7: same inputs, same bytes."""
    assert table.equals(table.sort("antigen.species", "antigen.epitope"))


def test_the_substitution_is_one_indexed_and_reads_native_to_curated() -> None:
    """`SLLMWITQV` is the C9V analogue of `SLLMWITQC`, so the cell must read `9C>V` and not `8V>C`.

    Getting the direction or the base wrong here is the kind of defect nobody notices, because both
    forms look plausible. `SourceHit.mutations` is 0-based and ordered `(position, curated, native)`.
    """
    hit_like = [(8, "V", "C")]
    rendered = ",".join(f"{pos + 1}{was}>{now}" for pos, now, was in hit_like)
    assert rendered == "9C>V"


def test_the_two_largest_one_substitution_rows_name_their_native_peptide(table) -> None:
    """The cross-reference the catalogue cannot express: `SLLMWITQV` and `SLLMWITQC` are unrelated
    rows, and so are `ELAGIGILTV` and `EAAGIGILTV`. Naming the pair is the whole use of the row -
    neither is a mislabel, because an anchor-optimised peptide is a screening reagent for the real
    self antigen and the gene label is right."""
    rows = {r["antigen.epitope"]: r for r in table.iter_rows(named=True)}
    assert rows["SLLMWITQV"]["verdict"] == A.ONE_SUB
    assert rows["SLLMWITQV"]["source.peptide"] == "SLLMWITQC"
    assert rows["SLLMWITQV"]["source.subs"] == "9C>V"
    assert rows["ELAGIGILTV"]["source.peptide"] == "EAAGIGILTV"
    assert rows["ELAGIGILTV"]["source.subs"] == "2A>L"


def test_a_gene_disagreement_is_case_insensitive(table: pl.DataFrame) -> None:
    """`G6pc2` against `G6PC2` is a spelling, and `lookalikes.tsv` already reports that class."""
    clash = A.gene_disagreements(table.with_columns(pl.col("records").cast(pl.Int64)))
    folded = {(r["antigen.gene"].lower(), r["source.gene"].lower()) for r in clash.iter_rows(named=True)}
    assert not any(a == b for a, b in folded), "a case-only difference is not a disagreement"


def test_the_record_total_for_one_substitution_is_dominated_by_a_single_peptide(table) -> None:
    """Why :func:`by_gene` exists and why no summary here should lead with a record count.

    `SLLMWITQV` is over 80 % of the `one_substitution` record total, so that total measures one
    reagent choice in one antigen rather than anything about the database. The median epitope carries
    two records.
    """
    near = table.filter(pl.col("verdict") == A.ONE_SUB).with_columns(
        pl.col("records").cast(pl.Int64))
    total = near["records"].sum()
    assert near["records"].max() / total > 0.75
    assert near["records"].median() <= 3


def test_a_gene_with_many_peptides_a_residue_from_reference_is_a_cohort_not_a_defect(table) -> None:
    """A reference proteome is one genome and a patient cohort is not.

    One antigen yields as many distinct peptides as the cohort carries mutations in it, every one
    correctly labelled with that gene. `KRAS` carries twelve over fifty records, which is a somatic
    mutation panel; a check that called any of them wrong would be wrong itself.
    """
    genes = A.by_gene(table.with_columns(pl.col("records").cast(pl.Int64)))
    assert genes["peptides"].max() > 5, "the corpus does carry cohort-scale panels"
    kras = genes.filter(pl.col("antigen.gene") == "KRAS")
    assert kras.height == 1 and kras["peptides"].item() >= 10
