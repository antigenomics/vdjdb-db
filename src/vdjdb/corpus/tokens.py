"""The token families, one function each, and the k-mer expansion they share.

A token is a prefixed string, so a consumer filters by family with a string comparison and needs no
second table. The prefixes are declared once in :data:`FAMILIES` and every function takes its prefix
from there, so a family cannot be spelled two ways.

Each family exists because a question needs it and no other token can stand in for it, which is the
test for adding one. ``ROADMAP.md`` section 12, phase 17, states the question per family.
"""
from __future__ import annotations

from dataclasses import dataclass

import polars as pl

#: The k of every k-mer family. 3 for both CDR3s and epitopes: over 2,118 epitopes of 7 to 25
#: residues and 180,048 distinct CDR3s, a 3-mer has a document frequency worth an inverse document
#: frequency, while a 5-mer is close to an identifier of its own sequence and a 2-mer is in nearly
#: every document.
K = 3

#: The separator between a k-mer and the gene it is scoped to. Not `*`, which appears in allele
#: names, and not `:`, which separates a prefix from its token.
SCOPE = "@"


@dataclass(frozen=True, slots=True)
class Family:
    """One token family: its prefix, and one line on what it is for."""

    prefix: str
    about: str


FAMILIES: dict[str, Family] = {
    "word": Family("w:", "a word of the title or abstract"),
    "v_gene": Family("v:", "V gene, allele dropped"),
    "j_gene": Family("j:", "J gene, allele dropped"),
    "cdr3_kmer": Family("k:", f"a CDR3 {K}-mer"),
    "cdr3_kmer_in_v": Family("kv:", f"a CDR3 {K}-mer scoped to the V gene carrying it"),
    "epitope": Family("e:", "the epitope sequence"),
    "epitope_kmer": Family("ek:", f"an epitope {K}-mer"),
    "antigen_species": Family("a:", "the antigen's source organism"),
    "host_species": Family("s:", "the host the receptor was sequenced from"),
    "antigen_gene": Family("g:", "the antigen's gene"),
    "mhc_allele": Family("m:", "the presenting allele, two fields"),
    "mhc_locus": Family("ml:", "the presenting locus"),
    "mhc_class": Family("mc:", "MHC class"),
}

#: Prefix to family name, for reading a family off a token.
BY_PREFIX: dict[str, str] = {f.prefix: name for name, f in FAMILIES.items()}


def family_of(token: str) -> str | None:
    """The family a token belongs to, longest prefix first so ``kv:`` never reads as ``k:``."""
    for prefix in sorted(BY_PREFIX, key=len, reverse=True):
        if token.startswith(prefix):
            return BY_PREFIX[prefix]
    return None


def gene(call: str) -> pl.Expr:
    """A V or J call with the allele dropped: ``TRBV9*01`` becomes ``TRBV9``.

    The allele is dropped because it is a resolution the corpus cannot support. 1,302 records call
    `TRAJ24` with no allele at all against 147 that name one, so keeping the allele would put the same
    gene in two tokens and give the commoner spelling the weaker statistics.
    """
    return pl.col(call).str.split("*").list.first()


def kmers(sequences: pl.Series, k: int = K) -> pl.DataFrame:
    """``(sequence, kmer)``, one row per position, over the **distinct** sequences given.

    Vectorised rather than a Python loop per sequence: one ``str.slice`` expression per offset,
    concatenated into a list column and exploded. The offsets run to the longest sequence present, and
    an offset past the end of a given sequence yields a null rather than an empty string, which
    ``drop_nulls`` removes. Null and not empty because empty is this repository's missing marker and
    polars is changing how it reads one inside a list (hard rule 6).

    Deduplicating first is hard rule 4 and not a cache: the k-mers of a sequence are a function of the
    sequence, so computing them once per distinct sequence inside one build cannot change an answer.
    Measured on the corpus, 180,048 distinct CDR3s rather than 286,047 chain rows.
    """
    distinct = sequences.drop_nulls().unique().sort()
    distinct = distinct.filter(distinct.str.len_chars() >= k)
    if distinct.is_empty():
        return pl.DataFrame(schema={"sequence": pl.Utf8, "kmer": pl.Utf8})
    longest = int(distinct.str.len_chars().max() or 0)
    frame = pl.DataFrame({"sequence": distinct})
    return (frame.with_columns(
                pl.concat_list([
                    pl.when(pl.col("sequence").str.len_chars() >= i + k)
                      .then(pl.col("sequence").str.slice(i, k))
                    for i in range(longest - k + 1)
                ]).alias("kmer"))
                # `empty_as_null=False` pins polars 2.0's behaviour, which is the future default and
                # what `release/changelog.py` already asks for. Equivalent here either way, since a
                # sequence shorter than k was filtered out above and so no list is empty.
                .explode("kmer", empty_as_null=False)
                .drop_nulls("kmer")
                .unique(subset=["sequence", "kmer"])
                .sort("sequence", "kmer"))


def _emit(frame: pl.DataFrame, family: str, expr: pl.Expr) -> pl.DataFrame:
    """``(reference.id, term)`` for one family, blanks dropped."""
    prefix = FAMILIES[family].prefix
    return (frame.select("reference.id", (pl.lit(prefix) + expr).alias("term"),
                         expr.alias("__raw"))
                 .filter(pl.col("__raw").is_not_null() & (pl.col("__raw") != ""))
                 .drop("__raw"))


def receptor(chains: pl.DataFrame, records: pl.DataFrame, *, k: int = K) -> pl.DataFrame:
    """``(reference.id, term)`` for the four receptor families, one row per occurrence.

    Occurrences rather than distinct pairs, because a paper reporting a V gene on four hundred
    receptors says something a paper reporting it once does not, and the term frequency is where that
    shows. The sublinear term frequency in :mod:`vdjdb.corpus.build` is what stops a large paper
    dominating.
    """
    joined = chains.join(records.select("record_id", "reference.id"), on="record_id", how="inner")
    parts = [
        _emit(joined, "v_gene", gene("v.segm")),
        _emit(joined, "j_gene", gene("j.segm")),
    ]
    cdr3_kmers = kmers(joined["cdr3"], k)
    with_kmers = joined.join(cdr3_kmers, left_on="cdr3", right_on="sequence", how="inner")
    parts.append(_emit(with_kmers, "cdr3_kmer", pl.col("kmer")))
    parts.append(_emit(with_kmers, "cdr3_kmer_in_v",
                       pl.col("kmer") + pl.lit(SCOPE) + gene("v.segm")))
    return pl.concat(parts, how="vertical")


def antigen(records: pl.DataFrame, *, k: int = K) -> pl.DataFrame:
    """``(reference.id, term)`` for the epitope, its k-mers, the antigen's origin and the host."""
    parts = [
        _emit(records, "epitope", pl.col("antigen.epitope")),
        _emit(records, "antigen_species", pl.col("antigen.species")),
        _emit(records, "antigen_gene", pl.col("antigen.gene")),
        # The host, distinct from the antigen's organism and prefixed differently, because
        # `a:HomoSapiens` (a self-antigen) and `s:HomoSapiens` (a human donor) are different claims
        # about a record. The `refsearch` client sends a host-species filter, so the contract needs
        # this family for the endpoint to be reproducible.
        _emit(records, "host_species", pl.col("species")),
    ]
    epitope_kmers = kmers(records["antigen.epitope"], k)
    with_kmers = records.join(epitope_kmers, left_on="antigen.epitope", right_on="sequence",
                              how="inner")
    parts.append(_emit(with_kmers, "epitope_kmer", pl.col("kmer")))
    return pl.concat(parts, how="vertical")


def restriction(records: pl.DataFrame, dictionary: pl.DataFrame) -> pl.DataFrame:
    """``(reference.id, term)`` for the three MHC granularities, from the dictionary.

    Joined through :mod:`vdjdb.corpus.mhc`'s table rather than parsed here, so the allele, the locus
    and the class a token names are the same ones the dictionary publishes. A record's two chains are
    handled separately: a class II molecule contributes its alpha locus and its beta locus, which is
    two restrictions and not one.
    """
    lookup = dictionary.select("mhc", "allele", "locus", "chain").unique()
    parts = []
    for side in ("a", "b"):
        joined = (records.select("reference.id", pl.col(f"mhc.{side}").alias("mhc"), "mhc.class")
                         .join(lookup.filter(pl.col("chain") == side).drop("chain"),
                               on="mhc", how="inner"))
        parts.append(_emit(joined, "mhc_allele", pl.col("allele")))
        parts.append(_emit(joined, "mhc_locus", pl.col("locus")))
    parts.append(_emit(records, "mhc_class", pl.col("mhc.class")))
    return pl.concat(parts, how="vertical")


def text(term_counts: pl.DataFrame) -> pl.DataFrame:
    """``(reference.id, term, tf)`` from the committed word counts.

    The only family that is already aggregated, because the running text it came from is an input to
    the build and never an output of it (hard rule 5): ``corpus/text_terms.tsv`` carries counts and no
    prose. :mod:`vdjdb.corpus.pubmed` is what writes it.
    """
    prefix = FAMILIES["word"].prefix
    return (term_counts.select("reference.id", (pl.lit(prefix) + pl.col("term")).alias("term"),
                               pl.col("tf").cast(pl.UInt32))
                       .filter(pl.col("term") != prefix))
