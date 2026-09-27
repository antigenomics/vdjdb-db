"""Documents, a vocabulary and postings, with tf-idf weights.

Three tables. Long postings rather than a sparse-matrix format, because the consumer is polars or
duckdb and the query is a join, and because a format nobody has to decode is the one a downstream tool
will actually read.

Every id is assigned from a total order, ``document_id`` by sorted ``reference.id`` and ``term_id`` by
sorted term, so the tables are reproducible without a hash and a diff between two builds is readable.
Postings are sorted by ``(term_id, document_id)``, which makes a term lookup one contiguous slice.
"""
from __future__ import annotations

import re
from pathlib import Path

import polars as pl

from . import mhc, pubmed, tokens

#: A document is a reference, and its kind is read off the identifier. Order matters: the first
#: pattern that matches wins, and ``other`` is what a reference nothing recognises gets. The full
#: resolver set, including the year each kind is dated by, is
#: :mod:`vdjdb.summary.references`; these patterns only have to name the kind.
KINDS: tuple[tuple[str, re.Pattern[str]], ...] = (
    ("pubmed", re.compile(r"^PMID:\s*\d+\s*$")),
    ("pdb", re.compile(r"^https?://www\.rcsb\.org/structure/", re.IGNORECASE)),
    ("preprint", re.compile(r"^https?://(arxiv\.org/abs/|doi\.org/10\.1101/)", re.IGNORECASE)),
    ("issue", re.compile(r"^https?://github\.com/[\w-]+/[\w-]+/issues/\d+")),
)

DOCUMENT_COLUMNS: tuple[str, ...] = (
    "document_id", "reference.id", "kind", "pmid", "year", "n_records", "n_terms",
)
TERM_COLUMNS: tuple[str, ...] = ("term_id", "term", "family", "df", "idf")
POSTING_COLUMNS: tuple[str, ...] = ("document_id", "term_id", "tf", "weight")

#: The three tables plus the MHC dictionary, in the order a reader should meet them.
TABLES: tuple[str, ...] = ("documents", "terms", "postings", "mhc")


def kind_of(reference_id: str) -> str:
    for name, pattern in KINDS:
        if pattern.match(reference_id):
            return name
    return "other"


def documents(records: pl.DataFrame, *, pubmed_records: pl.DataFrame | None = None) -> pl.DataFrame:
    """One row per distinct ``reference.id``, numbered in sorted order.

    ``n_terms`` is filled in by :func:`build` once the postings exist, because it is a property of the
    vocabulary rather than of the reference. A blank ``reference.id`` is dropped: 854 records carry
    none, all from one chunk, and a document with no identifier cannot be cited or scored.
    """
    counts = (records.filter(pl.col("reference.id") != "")
                     .group_by("reference.id").len().rename({"len": "n_records"})
                     .sort("reference.id"))
    out = (counts.with_row_index("document_id")
                 .with_columns(
                     pl.col("document_id").cast(pl.UInt32),
                     pl.col("reference.id")
                       .map_elements(kind_of, return_dtype=pl.Utf8).alias("kind"),
                     pl.col("reference.id")
                       .map_elements(lambda r: pubmed.pmid_of(r) or "", return_dtype=pl.Utf8)
                       .alias("pmid"),
                     pl.col("n_records").cast(pl.UInt32)))
    if pubmed_records is not None and not pubmed_records.is_empty():
        out = (out.join(pubmed_records.select("reference.id", "year"), on="reference.id", how="left")
                  .with_columns(pl.col("year").fill_null("")))
    else:
        out = out.with_columns(pl.lit("").alias("year"))
    return out.with_columns(pl.lit(0, pl.UInt32).alias("n_terms")).select(list(DOCUMENT_COLUMNS))


def occurrences(records: pl.DataFrame, chains: pl.DataFrame, restriction: pl.DataFrame,
                *, text_terms: pl.DataFrame | None = None,
                k: int = tokens.K, root: Path | None = None) -> pl.DataFrame:
    """``(reference.id, term, tf)`` for every family, term frequencies already summed.

    The receptor, antigen and MHC families emit one row per occurrence and are counted here; the text
    family arrives already counted, because its source text is not stored (hard rule 5).
    """
    dictionary = mhc.dictionary(restriction, root=root)
    raw = pl.concat([tokens.receptor(chains, records, k=k),
                     tokens.antigen(records, k=k),
                     tokens.restriction(records, dictionary)], how="vertical")
    counted = (raw.group_by("reference.id", "term").len().rename({"len": "tf"})
                  .with_columns(pl.col("tf").cast(pl.UInt32)))
    if text_terms is not None and not text_terms.is_empty():
        counted = pl.concat([counted, tokens.text(text_terms)], how="vertical")
    return counted.filter(pl.col("reference.id") != "").sort("reference.id", "term")


def weigh(counted: pl.DataFrame, n_documents: int) -> tuple[pl.DataFrame, pl.DataFrame]:
    """``(terms, postings)`` from ``(reference.id, term, tf)``.

    Sublinear term frequency ``1 + log tf``, because a paper reporting ten thousand receptors would
    otherwise dominate every receptor token by arithmetic rather than by relevance. Smoothed inverse
    document frequency ``log((N + 1) / (df + 1)) + 1``, so a term in every document still carries a
    weight of 1 rather than 0 and a term in none cannot divide by zero. L2 normalisation per document,
    so a long abstract and a short one are comparable.

    These are scikit-learn's ``TfidfVectorizer`` conventions, which is the reason to use them: the
    implementation is then checkable against a reference rather than only against itself, and
    ``tests/unit/test_corpus.py`` does exactly that.
    """
    terms = (counted.group_by("term").agg(pl.col("reference.id").n_unique().alias("df"))
                    .sort("term").with_row_index("term_id")
                    .with_columns(
                        pl.col("term_id").cast(pl.UInt32),
                        pl.col("df").cast(pl.UInt32),
                        pl.col("term").map_elements(tokens.family_of, return_dtype=pl.Utf8)
                          .alias("family"),
                        (((n_documents + 1) / (pl.col("df") + 1)).log() + 1).alias("idf"))
                    .select(list(TERM_COLUMNS)))
    postings = (counted.join(terms.select("term_id", "term", "idf"), on="term", how="inner")
                       .with_columns(((pl.col("tf").log() + 1) * pl.col("idf")).alias("weight")))
    # L2 per document, after weighting and before the id join, so the norm covers exactly the terms
    # this document carries.
    norms = (postings.group_by("reference.id")
                     .agg((pl.col("weight") ** 2).sum().sqrt().alias("norm")))
    postings = (postings.join(norms, on="reference.id", how="inner")
                        .with_columns(
                            pl.when(pl.col("norm") > 0)
                              .then(pl.col("weight") / pl.col("norm"))
                              .otherwise(0.0).alias("weight")))
    return terms, postings


def build(records: pl.DataFrame, chains: pl.DataFrame, restriction: pl.DataFrame, *,
          text_terms: pl.DataFrame | None = None, pubmed_records: pl.DataFrame | None = None,
          k: int = tokens.K, root: Path | None = None) -> dict[str, pl.DataFrame]:
    """The four corpus tables, keyed by name.

    Everything is recomputed from the tables handed in; nothing is read back from a previous run
    (hard rule 9). ``text_terms`` and ``pubmed_records`` are the committed inputs
    :mod:`vdjdb.corpus.pubmed` writes, and both are optional: with neither, the corpus builds from the
    receptor, antigen and MHC families and says the text family is absent.
    """
    docs = documents(records, pubmed_records=pubmed_records)
    counted = occurrences(records, chains, restriction, text_terms=text_terms, k=k, root=root)
    terms, postings = weigh(counted, docs.height)

    ids = docs.select("document_id", "reference.id")
    postings = (postings.join(ids, on="reference.id", how="inner")
                        .select(list(POSTING_COLUMNS))
                        .sort("term_id", "document_id"))
    per_document = (postings.group_by("document_id").len().rename({"len": "n_terms"})
                            .with_columns(pl.col("n_terms").cast(pl.UInt32)))
    docs = (docs.drop("n_terms")
                .join(per_document, on="document_id", how="left")
                .with_columns(pl.col("n_terms").fill_null(0).cast(pl.UInt32))
                .select(list(DOCUMENT_COLUMNS)).sort("document_id"))
    return {"documents": docs, "terms": terms, "postings": postings,
            "mhc": mhc.dictionary(restriction, root=root)}


def write(corpus: dict[str, pl.DataFrame], out: Path) -> dict[str, Path]:
    """Parquet and TSV side by side, the same pair the definitive tables ship as."""
    out.mkdir(parents=True, exist_ok=True)
    written = {}
    for name in TABLES:
        frame = corpus[name]
        frame.write_parquet(out / f"{name}.parquet")
        frame.write_csv(out / f"{name}.tsv", separator="\t", quote_style="necessary")
        written[name] = out / f"{name}.parquet"
    return written


def read(directory: Path) -> dict[str, pl.DataFrame]:
    """Read a built corpus back. Parquet, which is the form a query should use."""
    return {name: pl.read_parquet(Path(directory) / f"{name}.parquet")
            for name in TABLES if (Path(directory) / f"{name}.parquet").exists()}
