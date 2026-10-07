"""Provenance screens count clonotypes and preserve zero-match comparisons."""
import polars as pl

from vdjdb.curate.overlap import overlaps


def example() -> pl.DataFrame:
    return pl.DataFrame([
        {"chunk.file": f"{ref}.tsv", "reference.id": ref, "species": "HomoSapiens",
         "antigen.epitope": "NLVPMVATV", "mhc.a": "HLA-A*02:01", "mhc.b": "B2M",
         "mhc.class": "MHCI", "cdr3.alpha": alpha, "cdr3.beta": beta,
         "meta.subject.id": donor}
        for ref, alpha, beta, donor in [
            ("1", "CAAF", "CBBF", "a"), ("1", "CAAF", "CBBF", "a"),
            ("1", "CAAF", "CBBF", "b"), ("2", "CAAF", "CDDF", "a"),
            ("3", "CEEF", "CFFF", "a")]])


def test_modes_denominators_and_zero_pairs() -> None:
    result = overlaps(example())
    assert result.height == 12
    ab = result.filter((pl.col("reference.id") == "1") & (pl.col("reference.id.other") == "2"))
    assert dict(ab.select("mode", "shared").iter_rows()) == {
        "alpha": 1, "beta": 0, "paired": 0, "paired.pmhc": 0}
    assert ab["n"].to_list() == [1, 1, 1, 1]
    assert dict(ab.select("mode", "ratio.plus1").iter_rows()) == {
        "alpha": 1., "beta": .5, "paired": .5, "paired.pmhc": .5}
    assert not result["review.provenance"].any()


def test_donors_are_scoped_to_reference_and_input_order_is_irrelevant() -> None:
    result = overlaps(example(), by_sample=True)
    within = result.filter(pl.col("same.reference"))
    assert within.height == 4
    assert within["shared"].to_list() == [1, 1, 1, 1]
    assert result.equals(overlaps(example().reverse(), by_sample=True))
    # Reused donor labels across publications do not prove a shared sample.
    assert not result.filter(pl.col("reference.id") != pl.col("reference.id.other"))["same.reference"].any()


def test_strata_and_missing_chains_do_not_match() -> None:
    source = example().with_columns(
        pl.when(pl.col("reference.id") == "2").then(pl.lit("HLA-A*02:02"))
        .otherwise(pl.col("mhc.a")).alias("mhc.a"),
        pl.when(pl.col("reference.id") == "3").then(pl.lit(""))
        .otherwise(pl.col("cdr3.alpha")).alias("cdr3.alpha"))
    result = overlaps(source)
    assert result.height == 2
    assert set(result["mode"]) == {"beta", "paired.pmhc"}
    assert result["shared"][0] == 0


def test_submitted_sequences_and_high_overlap_screen() -> None:
    rows = []
    for ref in ["1", "2"]:
        for i in range(5):
            row = example().row(0, named=True)
            row.update({"chunk.file": ref, "reference.id": ref,
                        "cdr3.alpha": f"C{i}F", "cdr3.beta": f"C{i}F",
                        "__cdr3old.alpha": f"{ref}-{i}", "__cdr3old.beta": f"{ref}-{i}"})
            rows.append(row)
    source = pl.DataFrame(rows)
    repaired = overlaps(source)
    assert repaired["review.provenance"].all()
    assert repaired["size.product"].to_list() == [25, 25, 25, 25]
    assert repaired["ratio.plus1"].to_list() == [6/26]*4
    assert overlaps(source, submitted=True)["shared"].sum() == 0


def test_reuse_spread_across_peptide_variants_is_detected() -> None:
    source = example().filter(pl.col("reference.id") == "1").head(1)
    records = pl.concat([source.with_columns(pl.lit(ref).alias("reference.id"),
                                           pl.lit(ref).alias("chunk.file"),
                                           pl.lit(str(i)).alias("antigen.epitope"))
                         for ref in ["1", "2"] for i in range(10)])
    result = overlaps(records)
    portfolio = result.filter(pl.col("mode") == "paired.pmhc")
    assert portfolio["shared"].to_list() == [10]
    assert portfolio["review.provenance"].all()
    assert not result.filter(pl.col("mode") == "paired")["review.provenance"].any()
