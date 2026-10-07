"""Author independence has strict boundaries and cannot use missing metadata."""
import polars as pl

from vdjdb.curate.authors import author_pairs, parse_authors


def fixture() -> pl.DataFrame:
    rows = []
    for ref, names in [("1", ["a", "b", "senior"]), ("2", ["a", "d", "other"]),
                       ("3", ["f", "g", "h", "new"]), ("4", ["z", "senior"])]:
        rows += [(ref, i, name, name, True, True) for i, name in enumerate(names)]
    return pl.DataFrame(rows, orient="row", schema={
        "reference.id": pl.String, "ordinal": pl.Int64, "name": pl.String,
        "author.key": pl.String, "complete": pl.Boolean, "individual": pl.Boolean})


def test_strict_one_third_senior_overlap_and_unknown() -> None:
    pairs = pl.DataFrame({"reference.id": ["1"]*5, "reference.id.other": ["2", "3", "4", "5", "1"]})
    out = author_pairs(pairs, fixture())
    assert out["author.independence"].to_list() == [
        "same-or-overlapping-laboratory", "independent", "same-or-overlapping-laboratory",
        "unknown", "same-or-overlapping-laboratory"]
    assert out["authors.overlap"][0] == 1/3
    assert out.equals(author_pairs(pairs, fixture().reverse()))


def test_parser_order_names_truncation_and_initial_variants() -> None:
    xml = b'''<PubmedArticleSet><PubmedArticle><MedlineCitation><PMID>42</PMID><Article>
    <AuthorList CompleteYN="N"><Author><LastName>Meyer</LastName><ForeName>Hannah V</ForeName></Author>
    <Author><CollectiveName>Group</CollectiveName></Author></AuthorList>
    </Article></MedlineCitation></PubmedArticle></PubmedArticleSet>'''
    out = parse_authors(xml)
    assert out["ordinal"].to_list() == [1, 2]
    assert out["author.key"].to_list() == ["meyerh", "group"]
    assert not out["complete"].any()
    pairs = pl.DataFrame({"reference.id": ["PMID:42"], "reference.id.other": ["PMID:42"]})
    assert author_pairs(pairs, out)["author.independence"][0] == "unknown"


def test_four_classes_and_reuse_warning_do_not_confuse_assay_with_replication() -> None:
    from vdjdb.curate.authors import categories
    from vdjdb.curate.overlap import overlaps
    from vdjdb.score.confidence import SCORE_SIGNATURE

    rows = []
    for serial, ref, clone, replica, verification in [
        (1, "1", "one", "", ""), (2, "1", "two", "", "direct"),
        (3, "1", "three", "a", ""), (4, "1", "three", "b", ""),
        (5, "1", "four", "", ""), (6, "3", "four", "", ""),
        (7, "5", "one", "", "")]:
        row = dict.fromkeys(SCORE_SIGNATURE, "")
        row.update({"record_id": str(serial), "chunk.file": ref, "reference.id": ref,
                    "cdr3.alpha": clone, "cdr3.beta": clone, "species": "HomoSapiens",
                    "antigen.epitope": "NLVPMVATV", "mhc.a": "HLA-A*02:01",
                    "meta.replica.id": replica, "method.verification": verification})
        rows.append(row)
    records = pl.DataFrame(rows)
    result = categories(records, fixture(), overlaps(records))
    assert result["validation.category"].to_list() == [
        "1_observed_only", "2_observed_and_validated", "3_repeated_within_study",
        "3_repeated_within_study", "4_independent_corroboration", "4_independent_corroboration",
        "1_observed_only"]
    assert result["cross.reference.alpha"][0]  # Missing authors are not independent evidence.
    suspect = overlaps(records).with_columns(
        ((pl.col("reference.id") == "1") & (pl.col("reference.id.other") == "3"))
        .alias("review.provenance"))
    held = categories(records, fixture(), suspect)
    assert not held["independent.support"].any()
    assert held.filter(pl.col("provenance.warning"))["record_id"].to_list() == ["5", "6"]
    assert held.filter(pl.col("record_id") == "2")["assay.validated"].item()


def test_reviewed_input_covers_corpus_pmids_and_preserves_order() -> None:
    from vdjdb.corpus.pubmed import PMID
    from vdjdb.curate.authors import AUTHOR_TABLE
    from vdjdb.curate.nomenclature import harmonise_references
    from vdjdb.io.chunks import read_chunks

    records, _ = harmonise_references(read_chunks())
    needed = {ref for ref in records['reference.id'].unique() if PMID.fullmatch(ref)}
    authors = pl.read_csv(AUTHOR_TABLE, separator="\t")
    assert needed <= set(authors['reference.id'])
    assert authors.select('reference.id', 'ordinal').n_unique() == authors.height
    assert authors['author.key'].ne('').all()
    order = authors.group_by('reference.id').agg(pl.col('ordinal').min().alias('first'),
                                               pl.col('ordinal').max().alias('last'),
                                               pl.len().alias('n'))
    assert order['first'].eq(1).all()
    assert order['last'].equals(order['n'].cast(pl.Int64))
