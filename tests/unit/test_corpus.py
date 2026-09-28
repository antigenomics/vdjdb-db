"""The reference corpus: the MHC dictionary, the token families, the weighting, and the two queries.

The weighting is checked against ``sklearn.feature_extraction.text.TfidfVectorizer`` rather than only
against itself, which is the reason :mod:`vdjdb.corpus.build` uses that library's conventions: an
implementation agreeing with a reference is evidence, and one agreeing with its own last run is not.
"""
from __future__ import annotations

import math

import polars as pl
import pytest

from vdjdb.corpus import build, mhc, pubmed, query, tokens

# --------------------------------------------------------------------------------------------
# Fixtures: two papers, three receptors, two epitopes, one class I and one class II
# --------------------------------------------------------------------------------------------

RECORDS = pl.DataFrame([
    {"record_id": "VDJDB0000000001", "reference.id": "PMID:1", "species": "HomoSapiens",
     "antigen.epitope": "GILGFVFTL", "antigen.gene": "M", "antigen.species": "InfluenzaA",
     "mhc.a": "HLA-A*02:01", "mhc.b": "B2M", "mhc.class": "MHCI"},
    {"record_id": "VDJDB0000000002", "reference.id": "PMID:1", "species": "HomoSapiens",
     "antigen.epitope": "GILGFVFTL", "antigen.gene": "M", "antigen.species": "InfluenzaA",
     "mhc.a": "HLA-A*02:01", "mhc.b": "B2M", "mhc.class": "MHCI"},
    {"record_id": "VDJDB0000000003", "reference.id": "PMID:2", "species": "MusMusculus",
     "antigen.epitope": "ASNENMETM", "antigen.gene": "NP", "antigen.species": "InfluenzaA",
     "mhc.a": "H2-Db", "mhc.b": "B2M", "mhc.class": "MHCI"},
])

CHAINS = pl.DataFrame([
    {"record_id": "VDJDB0000000001", "gene": "TRB", "cdr3": "CASSIRSSYEQYF",
     "v.segm": "TRBV19*01", "j.segm": "TRBJ2-7*01"},
    {"record_id": "VDJDB0000000002", "gene": "TRB", "cdr3": "CASSIRSTYEQYF",
     "v.segm": "TRBV19*02", "j.segm": "TRBJ2-7*01"},
    {"record_id": "VDJDB0000000003", "gene": "TRB", "cdr3": "CASSLGGANTEVFF",
     "v.segm": "TRBV1*01", "j.segm": "TRBJ1-1*01"},
])

RESTRICTION = pl.DataFrame([
    {"antigen.epitope": "GILGFVFTL", "antigen.species": "InfluenzaA", "mhc.a": "HLA-A*02:01",
     "mhc.b": "B2M", "mhc.class": "MHCI"},
    {"antigen.epitope": "ASNENMETM", "antigen.species": "InfluenzaA", "mhc.a": "H2-Db",
     "mhc.b": "B2M", "mhc.class": "MHCI"},
])


@pytest.fixture(scope="module")
def corpus() -> dict[str, pl.DataFrame]:
    return build.build(RECORDS, CHAINS, RESTRICTION)


# --------------------------------------------------------------------------------------------
# The MHC dictionary
# --------------------------------------------------------------------------------------------

def test_an_allele_is_truncated_to_the_two_fields_the_database_curates() -> None:
    assert mhc.two_field("HLA-A*02:01:48:02") == "HLA-A*02:01"
    assert mhc.two_field("HLA-A*02:01") == "HLA-A*02:01"
    assert mhc.two_field("HLA-B*07") == "HLA-B*07"
    assert mhc.two_field("H2-Kb") == "H2-Kb"


def test_a_locus_is_the_gene_and_a_murine_haplotype_resolves_to_its_series() -> None:
    assert mhc.locus("HLA-DRB1*04:01") == "HLA-DRB1"
    assert mhc.locus("H2-Kb") == "H2-K"
    assert mhc.locus("H2-IAg7") == "H2-IA"


def test_an_imgt_murine_gene_name_is_its_own_locus() -> None:
    """`H2-Ab1` is a gene, not an allele of `H2-A`, so stripping a trailing letter would be wrong."""
    for gene in ("H2-Aa", "H2-Ab1", "H2-Eb1", "B2M"):
        assert mhc.locus(gene) == gene


def test_the_dictionary_covers_both_chains_and_carries_the_imgt_verdict() -> None:
    d = mhc.dictionary(RESTRICTION)
    assert set(d["chain"]) == {"a", "b"}
    assert d.filter(pl.col("mhc") == "HLA-A*02:01")["status"][0] == "known"
    # Murine and the light chain are checked against proofreading/mhc_nonhuman.tsv, not the HLA
    # database, so `declared` rather than `unknown`.
    assert set(d.filter(pl.col("mhc").is_in(["H2-Db", "B2M"]))["status"]) == {"declared"}
    assert mhc.unrecognised(d).is_empty()


def test_a_call_imgt_has_at_no_depth_is_a_finding() -> None:
    """`HLA-A*08:01` is in the corpus on 74 records and no `HLA-A*08` exists at any resolution."""
    bad = RESTRICTION.with_columns(pl.lit("HLA-A*08:01").alias("mhc.a"))
    found = mhc.unrecognised(mhc.dictionary(bad))
    assert found.height == 1
    assert found["mhc"][0] == "HLA-A*08:01"


# --------------------------------------------------------------------------------------------
# Tokens
# --------------------------------------------------------------------------------------------

def test_every_family_has_a_distinct_prefix() -> None:
    prefixes = [f.prefix for f in tokens.FAMILIES.values()]
    assert len(set(prefixes)) == len(prefixes)


def test_a_longer_prefix_wins_so_kv_never_reads_as_k() -> None:
    assert tokens.family_of("kv:CAS@TRBV9") == "cdr3_kmer_in_v"
    assert tokens.family_of("k:CAS") == "cdr3_kmer"
    assert tokens.family_of("m:HLA-A*02:01") == "mhc_allele"
    assert tokens.family_of("ml:HLA-A") == "mhc_locus"
    assert tokens.family_of("mc:MHCI") == "mhc_class"
    assert tokens.family_of("zz:nothing") is None


def test_kmers_are_every_window_and_nothing_shorter() -> None:
    got = tokens.kmers(pl.Series(["CASSL", "AB"]), 3)
    assert got["kmer"].to_list() == ["ASS", "CAS", "SSL"]
    assert (got["kmer"].str.len_chars() == 3).all()
    assert "AB" not in got["sequence"].to_list()


def test_kmers_of_nothing_is_an_empty_frame_with_the_right_shape() -> None:
    empty = tokens.kmers(pl.Series([], dtype=pl.Utf8), 3)
    assert empty.is_empty()
    assert empty.columns == ["sequence", "kmer"]


def test_a_v_call_loses_its_allele_so_one_gene_is_one_token() -> None:
    out = tokens.receptor(CHAINS, RECORDS)
    v = {t for t in out["term"] if t.startswith("v:")}
    assert v == {"v:TRBV19", "v:TRBV1"}, "TRBV19*01 and *02 are one gene"


def test_a_kmer_is_emitted_both_bare_and_scoped_to_its_v_gene() -> None:
    """The pair is what makes "is the motif specific, or is its V gene?" askable."""
    out = tokens.receptor(CHAINS, RECORDS)
    terms = set(out["term"])
    assert "k:IRS" in terms
    assert "kv:IRS@TRBV19" in terms


def test_the_epitope_is_emitted_whole_and_as_kmers() -> None:
    out = tokens.antigen(RECORDS)
    terms = set(out["term"])
    assert "e:GILGFVFTL" in terms
    assert {"ek:GIL", "ek:ILG", "ek:FTL"} <= terms


def test_the_host_and_the_antigen_organism_are_different_families() -> None:
    """`a:HomoSapiens` is a self-antigen and `s:HomoSapiens` is a human donor."""
    out = tokens.antigen(RECORDS)
    assert "s:HomoSapiens" in set(out["term"])
    assert "a:InfluenzaA" in set(out["term"])
    assert "a:HomoSapiens" not in set(out["term"])


def test_the_three_mhc_granularities_come_from_the_dictionary() -> None:
    out = tokens.restriction(RECORDS, mhc.dictionary(RESTRICTION))
    terms = set(out["term"])
    assert {"m:HLA-A*02:01", "ml:HLA-A", "mc:MHCI"} <= terms
    assert {"m:H2-Db", "ml:H2-D"} <= terms


def test_a_blank_field_emits_no_token() -> None:
    blank = RECORDS.with_columns(pl.lit("").alias("antigen.gene"))
    assert not any(t.startswith("g:") for t in tokens.antigen(blank)["term"])


# --------------------------------------------------------------------------------------------
# The weighting, against a reference implementation
# --------------------------------------------------------------------------------------------

def test_the_weights_match_sklearns_tfidfvectorizer() -> None:
    """Sublinear tf, smoothed idf, L2 per document: exactly `TfidfVectorizer`'s defaults plus
    `sublinear_tf=True`.

    Checked against the library rather than against a frozen table of our own output, because the
    latter only proves the code has not changed and this proves it is right.
    """
    sklearn_text = pytest.importorskip("sklearn.feature_extraction.text")
    documents = {"PMID:1": "alpha beta beta gamma", "PMID:2": "beta gamma gamma gamma",
                 "PMID:3": "alpha alpha delta"}
    counted = pl.DataFrame(
        [{"reference.id": ref, "term": term, "tf": text.split().count(term)}
         for ref, text in documents.items() for term in sorted(set(text.split()))],
        schema={"reference.id": pl.Utf8, "term": pl.Utf8, "tf": pl.UInt32})

    terms, postings = build.weigh(counted, len(documents))
    ours = {(r["reference.id"], r["term"]): r["weight"]
            for r in postings.select("reference.id", "term", "weight").iter_rows(named=True)}

    vec = sklearn_text.TfidfVectorizer(sublinear_tf=True, norm="l2", smooth_idf=True,
                                       token_pattern=r"\S+")
    matrix = vec.fit_transform(list(documents.values()))
    vocabulary = vec.get_feature_names_out()
    for row, ref in enumerate(documents):
        for col, term in enumerate(vocabulary):
            want = matrix[row, col]
            got = ours.get((ref, term), 0.0)
            assert got == pytest.approx(want, abs=1e-9), (ref, term)

    for row in terms.iter_rows(named=True):
        want = math.log((len(documents) + 1) / (row["df"] + 1)) + 1
        assert row["idf"] == pytest.approx(want, abs=1e-12)


def test_every_document_is_l2_normalised(corpus) -> None:
    norms = (corpus["postings"].group_by("document_id")
                               .agg((pl.col("weight") ** 2).sum().sqrt().alias("norm")))
    assert norms["norm"].to_list() == pytest.approx([1.0] * norms.height, abs=1e-12)


def test_ids_are_assigned_from_a_total_order_so_two_builds_agree(corpus) -> None:
    again = build.build(RECORDS, CHAINS, RESTRICTION)
    for name in build.TABLES:
        assert corpus[name].equals(again[name]), name
    shuffled = build.build(RECORDS.reverse(), CHAINS.reverse(), RESTRICTION.reverse())
    assert corpus["documents"].equals(shuffled["documents"])
    assert corpus["terms"].equals(shuffled["terms"])
    assert corpus["postings"].equals(shuffled["postings"])


def test_postings_are_sorted_so_a_term_lookup_is_one_slice(corpus) -> None:
    p = corpus["postings"]
    assert p.equals(p.sort("term_id", "document_id"))


def test_a_document_is_a_reference_and_its_kind_is_read_off_the_identifier() -> None:
    assert build.kind_of("PMID:28629751") == "pubmed"
    assert build.kind_of("https://www.rcsb.org/structure/1AO7") == "pdb"
    assert build.kind_of("https://doi.org/10.1101/2020.05.04.20085779") == "preprint"
    assert build.kind_of("https://github.com/antigenomics/vdjdb-db/issues/193") == "issue"
    assert build.kind_of("http://mediatum.ub.tum.de/doc/1136748") == "other"


def test_a_record_with_no_reference_is_not_a_document() -> None:
    """854 records carry no `reference.id`. A document with no identifier cannot be cited."""
    anonymous = RECORDS.with_columns(
        pl.when(pl.col("record_id") == "VDJDB0000000003").then(pl.lit(""))
          .otherwise(pl.col("reference.id")).alias("reference.id"))
    docs = build.documents(anonymous)
    assert docs.height == 1
    assert docs["reference.id"].to_list() == ["PMID:1"]


def test_the_corpus_builds_without_the_text_family(corpus) -> None:
    """A fork and a first build have no `text_terms.tsv`, and are expected to work."""
    assert not any(f == "word" for f in corpus["terms"]["family"])
    assert corpus["documents"].height == 2


def test_the_text_family_is_added_when_the_counts_are_there() -> None:
    counts = pl.DataFrame({"reference.id": ["PMID:1"], "term": ["influenza"], "tf": [3]},
                          schema={"reference.id": pl.Utf8, "term": pl.Utf8, "tf": pl.UInt32})
    c = build.build(RECORDS, CHAINS, RESTRICTION, text_terms=counts)
    assert "w:influenza" in set(c["terms"]["term"])
    assert c["terms"].filter(pl.col("term") == "w:influenza")["family"][0] == "word"


# --------------------------------------------------------------------------------------------
# Querying
# --------------------------------------------------------------------------------------------

def test_a_query_ranks_the_document_carrying_its_tokens(corpus) -> None:
    hits = query.score(corpus, ["e:GILGFVFTL"])
    assert hits["reference.id"].to_list() == ["PMID:1"]
    assert hits["score"][0] > 0


def test_a_query_of_unknown_tokens_returns_nothing_rather_than_everything(corpus) -> None:
    empty = query.score(corpus, ["e:NOTANEPITOPE"])
    assert empty.is_empty()
    assert empty.columns == ["reference.id", "pmid", "score", "matched"]


def test_the_ranking_is_total_so_a_tie_does_not_reorder(corpus) -> None:
    once = query.score(corpus, ["mc:MHCI"], limit=10)
    twice = query.score(corpus, ["mc:MHCI"], limit=10)
    assert once.equals(twice)
    assert once["reference.id"].to_list() == sorted(once["reference.id"].to_list())


def test_a_lift_of_one_means_the_condition_says_nothing(corpus) -> None:
    """`mc:MHCI` is on both documents, so conditioning on it cannot change any rate."""
    assert query.lift(corpus, "e:GILGFVFTL", ["mc:MHCI"]).lift == pytest.approx(1.0)


def test_a_lift_is_null_where_nothing_carries_the_condition(corpus) -> None:
    """Not a lift of zero, and it must never be averaged as one."""
    got = query.lift(corpus, "k:IRS", ["e:NOTANEPITOPE"])
    assert got.lift is None
    assert got.given_units == 0
    assert "no lift" in str(got)


def test_a_lift_carries_every_number_it_rests_on(corpus) -> None:
    """CLAUDE.md section 5a: a bare ratio is not a result."""
    got = query.lift(corpus, "k:IRS", ["e:GILGFVFTL"])
    text = str(got)
    assert f"{got.both}/{got.given_units}" in text
    assert f"{got.term_units}/{got.units}" in text
    assert got.over in text


def test_the_two_lift_modes_count_different_things(corpus) -> None:
    by_document = query.lift(corpus, "k:IRS", ["e:GILGFVFTL"], over="documents")
    by_occurrence = query.lift(corpus, "k:IRS", ["e:GILGFVFTL"], over="occurrences")
    assert by_document.units == corpus["documents"].height
    assert by_occurrence.units > by_document.units
    with pytest.raises(ValueError, match="over must be one of"):
        query.lift(corpus, "k:IRS", over="records")


def test_an_occurrence_lift_is_scoped_to_the_terms_own_family(corpus) -> None:
    """Otherwise a CDR3 k-mer's rate would move with how many epitopes a paper studied."""
    got = query.lift(corpus, "k:IRS", over="occurrences")
    kmer_total = (corpus["postings"]
                  .join(corpus["terms"].filter(pl.col("family") == "cdr3_kmer").select("term_id"),
                        on="term_id", how="inner")["tf"].sum())
    assert got.units == int(kmer_total)
    assert got.lift == pytest.approx(1.0), "an unconditioned lift is 1 by construction"


# --------------------------------------------------------------------------------------------
# The refsearch contract
# --------------------------------------------------------------------------------------------

def test_a_cdr3_motif_becomes_its_kmers_because_a_query_is_a_motif() -> None:
    assert query.refsearch_query(cdr3="CASSI", species_to_search="HomoSapiens") == [
        "k:CAS", "k:ASS", "k:SSI", "s:HomoSapiens"]


def test_a_motif_shorter_than_k_is_dropped_rather_than_matched_loosely() -> None:
    assert query.refsearch_query(cdr3="CA", species_to_search="HomoSapiens") == ["s:HomoSapiens"]


def test_an_epitope_contributes_both_itself_and_its_kmers() -> None:
    got = query.refsearch_query(epitope="GILGFVFTL", species_to_search="HomoSapiens")
    assert "e:GILGFVFTL" in got
    assert "ek:GIL" in got


def test_the_client_default_species_are_the_three_it_sends() -> None:
    """Read from `vdjdb-web`'s own client: the default when the user picks none."""
    got = query.refsearch_query(epitope="GILGFVFTL")
    assert [t for t in got if t.startswith("s:")] == [
        "s:HomoSapiens", "s:MusMusculus", "s:MacacaMulatta"]


def test_filter_stop_words_removes_only_word_tokens() -> None:
    with_stop = query.refsearch_query(epitope="the", extra_parameters="search_by_antigen")
    assert "w:the" in with_stop
    without = query.refsearch_query(epitope="the",
                                    extra_parameters="search_by_antigen filter_stop_words")
    assert "w:the" not in without
    assert "e:THE" in without, "only the word family is filtered"


# --------------------------------------------------------------------------------------------
# PubMed: parsing and tokenising, with no network
# --------------------------------------------------------------------------------------------

ARTICLE = b"""<?xml version="1.0"?>
<PubmedArticleSet><PubmedArticle><MedlineCitation>
  <PMID Version="1">28629751</PMID>
  <Article>
    <Journal><ISOAbbreviation>J Allergy Clin Immunol</ISOAbbreviation>
      <JournalIssue><PubDate><Year>2017</Year></PubDate></JournalIssue></Journal>
    <ArticleTitle>Unique <i>influenza</i> A memory repertoire</ArticleTitle>
    <Abstract>
      <AbstractText Label="BACKGROUND">The HLA-A2-restricted response.</AbstractText>
      <AbstractText Label="RESULTS">CD8 T-cells in 2019.</AbstractText>
    </Abstract>
  </Article>
</MedlineCitation>
<PubmedData><ArticleIdList>
  <ArticleId IdType="pubmed">28629751</ArticleId>
  <ArticleId IdType="doi">10.1016/j.jaci.2017.05.037</ArticleId>
</ArticleIdList></PubmedData></PubmedArticle></PubmedArticleSet>"""


def test_a_record_is_parsed_whole() -> None:
    (record,) = pubmed.parse(ARTICLE)
    assert record["pmid"] == "28629751"
    assert record["year"] == "2017"
    assert record["journal"] == "J Allergy Clin Immunol"
    assert record["doi"] == "10.1016/j.jaci.2017.05.037"


def test_inline_markup_does_not_truncate_a_title() -> None:
    """`node.text` alone stops at the first `<i>`, losing the rest of the title with no error."""
    (record,) = pubmed.parse(ARTICLE)
    assert record["title"] == "Unique influenza A memory repertoire"


def test_an_abstract_keeps_its_section_labels_as_text() -> None:
    """`BACKGROUND` and `RESULTS` are words a query can match."""
    (record,) = pubmed.parse(ARTICLE)
    assert str(record["abstract"]).startswith("BACKGROUND The HLA-A2-restricted")
    assert "RESULTS CD8 T-cells" in str(record["abstract"])


def test_a_hyphenated_term_stays_one_token() -> None:
    counts = pubmed.tokenise("The HLA-A2-restricted T-cell response to SARS-CoV-2.")
    assert "hla-a2-restricted" in counts
    assert "t-cell" in counts
    assert "sars-cov-2" in counts


def test_short_and_purely_numeric_tokens_are_dropped() -> None:
    counts = pubmed.tokenise("A study of 42 donors in 2019 at pH 7")
    assert "42" not in counts and "2019" not in counts and "7" not in counts
    assert "of" not in counts and "in" not in counts and "at" not in counts
    assert counts["study"] == 1


def test_a_pmid_is_read_tolerantly_and_a_non_pubmed_reference_is_not_one() -> None:
    assert pubmed.pmid_of("PMID:28629751") == "28629751"
    assert pubmed.pmid_of("PMID: 34433824") == "34433824", "22 records carry the stray space"
    assert pubmed.pmid_of("https://www.rcsb.org/structure/1AO7") is None


def test_the_stop_list_holds_words_and_nothing_else() -> None:
    assert "the" in pubmed.STOP_WORDS
    assert "influenza" not in pubmed.STOP_WORDS
    assert all(w.isalpha() for w in pubmed.STOP_WORDS)


def test_loading_absent_text_terms_gives_an_empty_frame_of_the_right_shape(tmp_path) -> None:
    got = pubmed.load_terms(tmp_path / "nothing.tsv")
    assert got.is_empty()
    assert got.columns == ["reference.id", "term", "tf"]


def test_the_corpus_survives_a_round_trip_through_files(corpus, tmp_path) -> None:
    build.write(corpus, tmp_path)
    again = build.read(tmp_path)
    for name in build.TABLES:
        assert corpus[name].equals(again[name]), name


def test_two_spellings_of_one_pmid_are_two_documents(monkeypatch) -> None:
    """The corpus spells PMID 34433824 both `PMID:34433824` and `PMID: 34433824`.

    A one-PMID-to-one-reference mapping kept whichever spelling the iteration reached last and dropped
    the other with nothing raised, which is how the first run of `corpus refs` wrote 609 records where
    610 references are PubMed ones. Both are documents, because a document is the `reference.id` the
    records cite and the records cite both.
    """
    monkeypatch.setattr(pubmed, "fetch", lambda ids: pubmed.parse(ARTICLE))
    records, terms, missing = pubmed.build_tables(
        ["PMID:28629751", "PMID: 28629751", "https://www.rcsb.org/structure/1AO7"])
    assert records["reference.id"].to_list() == ["PMID: 28629751", "PMID:28629751"]
    assert records["abstract_sha256"].n_unique() == 1, "one article, so one digest"
    assert terms["reference.id"].n_unique() == 2
    # The PDB reference is not missing: it is a document with no word family, which is expected and
    # is what `documents.kind` records. Reporting 51 of those every run would bury a dead PMID.
    assert missing == []


def test_a_reference_that_resolves_to_nothing_is_reported(monkeypatch) -> None:
    """Not silently dropped: a `reference.id` is curated, so a dead one is a curation finding."""
    monkeypatch.setattr(pubmed, "fetch", lambda ids: [])
    _, _, missing = pubmed.build_tables(["PMID:1", "PMID:2"])
    assert missing == ["PMID:1", "PMID:2"]


def test_lift_family_agrees_with_lift_term_for_term(corpus) -> None:
    """The batched form is a speedup, so it has to be the same number.

    `lift_family` exists because scoring the 2,342 CDR3 3-mers through `lift` in a loop took 50.7 s
    against 0.01 s batched, about 2,900x (`ROADMAP_local.md` section 57.1). A fast path that disagrees
    with the documented one is worse than a slow path, so this pins them together: every term of every
    family, both `over` modes, and every number rather than only the lift.
    """
    for family in sorted(set(corpus["terms"]["family"])):
        terms = corpus["terms"].filter(pl.col("family") == family)["term"].to_list()
        for over in ("documents", "occurrences"):
            batched = {r["term"]: r for r in
                       query.lift_family(corpus, family, over=over).iter_rows(named=True)}
            assert set(batched) == set(terms), f"{family}/{over}: term set differs"
            for term in terms:
                one, got = query.lift(corpus, term, over=over), batched[term]
                for field in ("units", "given_units", "term_units", "both"):
                    assert got[field] == getattr(one, field), \
                        f"{family}/{over}/{term}: {field} {got[field]} != {getattr(one, field)}"
                assert (got["lift"] is None) == (one.lift is None), f"{term}: lift nullity differs"
                if one.lift is not None:
                    assert got["lift"] == pytest.approx(one.lift, rel=1e-12), f"{term}: lift differs"


def test_lift_family_honours_a_condition_and_the_occurrence_floor(corpus) -> None:
    """Conditioning and `min_units` are the two things the loop did around `lift`, so they move in."""
    epitope = corpus["terms"].filter(pl.col("family") == "epitope")["term"][0]
    conditioned = query.lift_family(corpus, "cdr3_kmer", [epitope], over="occurrences")
    assert conditioned.height, "conditioning on a present epitope returned nothing"
    assert (conditioned["given_units"] <= conditioned["units"]).all()
    floor = int(conditioned["both"].max())
    assert query.lift_family(corpus, "cdr3_kmer", [epitope], over="occurrences",
                             min_units=floor + 1).height == 0, "min_units did not filter"
