"""PubMed author lists as a reviewed input to offline provenance checks."""
from __future__ import annotations

import unicodedata
import urllib.parse
from pathlib import Path
from xml.etree import ElementTree

import polars as pl

AUTHOR_TABLE = Path("proofreading/pubmed_authors.tsv")


def author_key(name: str) -> str:
    """Case/accent/punctuation folding; no fuzzy identity claims."""
    return "".join(c for c in unicodedata.normalize("NFKD", name.casefold()) if c.isalnum())


def parse_authors(xml: bytes) -> pl.DataFrame:
    rows = []
    for article in ElementTree.fromstring(xml).findall(".//PubmedArticle"):
        pmid = article.findtext("./MedlineCitation/PMID", "")
        authors = article.find("./MedlineCitation/Article/AuthorList")
        if authors is None:
            continue
        complete = authors.get("CompleteYN", "Y") == "Y"
        for ordinal, author in enumerate(authors.findall("Author"), 1):
            last = author.findtext("LastName", "")
            first = author.findtext("ForeName", author.findtext("Initials", ""))
            name = f"{first} {last}".strip() if last else author.findtext("CollectiveName", "")
            key = author_key(f"{last} {first[:1]}") if last and first else author_key(name)
            rows.append((f"PMID:{pmid}", ordinal, name, key, complete, bool(last and first)))
    return pl.DataFrame(rows, orient="row", schema={
        "reference.id": pl.String, "ordinal": pl.Int64, "name": pl.String,
        "author.key": pl.String, "complete": pl.Boolean, "individual": pl.Boolean,
    }).sort("reference.id", "ordinal")


def fetch_authors(references: list[str]) -> pl.DataFrame:
    """Batched PubMed retrieval; never called by a build."""
    from ..corpus.pubmed import CHUNK, EUTILS, PMID, http_get

    ids = sorted({match[1] for ref in references if (match := PMID.fullmatch(ref))})
    parts = []
    for start in range(0, len(ids), CHUNK):
        query = urllib.parse.urlencode({"db": "pubmed", "retmode": "xml",
                                        "id": ",".join(ids[start:start + CHUNK])})
        parts.append(parse_authors(http_get(f"{EUTILS}/efetch.fcgi?{query}")))
    return pl.concat(parts) if parts else parse_authors(b"<PubmedArticleSet/>")


def author_pairs(pairs: pl.DataFrame, authors: pl.DataFrame) -> pl.DataFrame:
    """Different senior authors and <1/3 overlap of the smaller full author list.

    Also exclude a senior author appearing anywhere in the other list. Missing,
    truncated or consortium-only metadata cannot establish author independence.
    """
    lists = (authors.sort("reference.id", "ordinal").group_by("reference.id")
             .agg(pl.col("author.key").unique(maintain_order=True).alias("authors"),
                  pl.col("author.key").last().alias("last.author"),
                  pl.col("name").last().alias("last.author.name"),
                  (pl.col("complete").all() & pl.col("individual").last()
                   & (pl.col("author.key") != "").all()).alias("authors.known")))
    joined = pairs.join(lists, on="reference.id", how="left").join(
        lists.rename({c: f"{c}.other" for c in lists.columns}), on="reference.id.other", how="left")
    joined = joined.with_columns(
        pl.col("authors").list.set_intersection("authors.other").list.len().alias("authors.shared"),
        pl.min_horizontal(pl.col("authors").list.len(), pl.col("authors.other").list.len())
        .alias("authors.smaller"),
        (pl.col("authors.known").fill_null(False) & pl.col("authors.known.other").fill_null(False))
        .alias("authors.available"))
    return (joined.with_columns(
        (pl.col("authors.shared") / pl.col("authors.smaller")).alias("authors.overlap"),
        pl.when(~pl.col("authors.available")).then(pl.lit("unknown"))
        .when((pl.col("reference.id") != pl.col("reference.id.other"))
              & ~pl.col("authors").list.contains(pl.col("last.author.other"))
              & ~pl.col("authors.other").list.contains(pl.col("last.author"))
              & (3 * pl.col("authors.shared") < pl.col("authors.smaller")))
        .then(pl.lit("independent"))
        .otherwise(pl.lit("same-or-overlapping-laboratory")).alias("author.independence"))
        .drop("authors", "authors.other", "authors.known", "authors.known.other", "authors.available"))


def categories(records: pl.DataFrame, authors: pl.DataFrame,
               pairs: pl.DataFrame) -> pl.DataFrame:
    """Four evidence classes, with author-qualified support and provenance warnings.

    This audit reports evidence from explicit assay/replicate metadata. It does not
    rewrite confidence scores or the historical chain-support flags.
    """
    from ..score.confidence import SCORE_SIGNATURE

    optional = ["method.verification", "meta.subject.id", "meta.replica.id", "chunk.row", "chunk.id"]
    df = records.with_columns(*(pl.lit("").alias(c) for c in optional if c not in records.columns))
    signature = list(SCORE_SIGNATURE)
    repeat = df.group_by(*signature, "reference.id").agg(
        (pl.col("meta.replica.id").filter(pl.col("meta.replica.id") != "").n_unique() > 1)
        .alias("replica.repeat"),
        (pl.col("meta.subject.id").filter(pl.col("meta.subject.id") != "").n_unique() > 1)
        .alias("subject.repeat"))
    df = df.join(repeat, on=[*signature, "reference.id"], validate="m:1")
    supports = []
    candidates = []
    flagged = pairs.filter(pl.col("review.provenance") & ~pl.col("same.reference"))
    suspect = pl.concat([
        flagged.select("reference.id", "reference.id.other"),
        flagged.select(pl.col("reference.id.other").alias("reference.id"),
                       pl.col("reference.id").alias("reference.id.other")),
    ]).unique().with_columns(pl.lit(True).alias("source.warning"))
    for chain in ["alpha", "beta"]:
        key = ["species", f"cdr3.{chain}", f"v.{chain}", f"j.{chain}",
               "antigen.epitope", "mhc.a", "mhc.b", "mhc.class"]
        unique = df.filter(pl.col(f"cdr3.{chain}") != "").select(*key, "reference.id").unique()
        other = unique.join(unique, on=key, suffix=".other").filter(
            pl.col("reference.id") != pl.col("reference.id.other"))
        other = author_pairs(other, authors).join(
            suspect, on=["reference.id", "reference.id.other"], how="left", validate="m:1")
        other = other.with_columns(pl.col("source.warning").fill_null(False))
        candidates.append(other.select(*key, "reference.id").unique()
                          .with_columns(pl.lit(True).alias(f"cross.reference.{chain}")))
        eligible = other.filter((pl.col("author.independence") == "independent")
                                & ~pl.col("source.warning"))
        supports.append(eligible.select(*key, "reference.id").unique()
                        .with_columns(pl.lit(True).alias(f"author.support.{chain}")))
        warning = other.filter(pl.col("source.warning")).select(*key, "reference.id").unique()
        warning = warning.with_columns(pl.lit(True).alias(f"source.warning.{chain}"))
        for frame in [candidates[-1], supports[-1], warning]:
            df = df.join(frame, on=[*key, "reference.id"], how="left", validate="m:1")
    df = df.with_columns(
        (pl.col("source.warning.alpha").fill_null(False)
         | pl.col("source.warning.beta").fill_null(False)).alias("provenance.warning"),
        pl.col("method.verification").str.strip_chars().ne("").alias("assay.validated"),
        (pl.col("replica.repeat") | pl.col("subject.repeat")).alias("within.study.repeat"),
        *(pl.col(c).fill_null(False) for c in ["author.support.alpha", "author.support.beta",
                                              "cross.reference.alpha", "cross.reference.beta"]))
    df = df.with_columns(
        (pl.col("author.support.alpha") | pl.col("author.support.beta"))
        .alias("independent.support"))
    return (df.with_columns(
        pl.when(pl.col("independent.support")).then(pl.lit("4_independent_corroboration"))
        .when(pl.col("assay.validated")).then(pl.lit("2_observed_and_validated"))
        .when(pl.col("within.study.repeat")).then(pl.lit("3_repeated_within_study"))
        .otherwise(pl.lit("1_observed_only")).alias("validation.category"))
        .select("record_id", "chunk.file", "chunk.row", "chunk.id", "reference.id",
                "validation.category", "assay.validated", "within.study.repeat",
                "author.support.alpha", "author.support.beta", "cross.reference.alpha",
                "cross.reference.beta", "independent.support", "provenance.warning")
        .sort("record_id"))
