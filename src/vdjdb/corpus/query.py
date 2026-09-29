"""Two questions you can ask a corpus.

:func:`score` ranks documents by how well they match a set of tokens, which is what
``vdjdb.com/refsearch/`` returns as ``tf_idf``. :func:`lift` asks whether a token goes with a group of
documents because of itself or because of something it travels with, which the endpoint cannot ask and
which is the reason the corpus is an artifact rather than an index inside a service.

Both are joins over ``postings``. Neither builds anything, so a caller can ask a hundred questions of
one corpus without rebuilding it.
"""
from __future__ import annotations

from dataclasses import dataclass

import polars as pl

from . import pubmed, tokens


def _terms_of(corpus: dict[str, pl.DataFrame], query: list[str]) -> pl.DataFrame:
    """``(term_id, term, idf)`` for the query tokens the vocabulary has. Unknown tokens are dropped."""
    wanted = pl.DataFrame({"term": sorted(set(query))}, schema={"term": pl.Utf8})
    return corpus["terms"].join(wanted, on="term", how="inner").select("term_id", "term", "idf")


def score(corpus: dict[str, pl.DataFrame], query: list[str], *, limit: int = 10) -> pl.DataFrame:
    """Documents ranked by the sum of their matched term weights.

    The sum rather than a cosine against a normalised query vector, because every posting is already
    L2-normalised per document and the query has no length to normalise against. That is the quantity
    the ``refsearch`` endpoint calls ``tf_idf``.

    Returns ``reference.id``, ``pmid``, ``score`` and ``matched`` - how many of the query's tokens the
    document carries, which is what separates a document matching one rare token from one matching
    five. Ties break on ``reference.id``, so the ranking is total and reproducible.
    """
    terms = _terms_of(corpus, query)
    if terms.is_empty():
        return pl.DataFrame(schema={"reference.id": pl.Utf8, "pmid": pl.Utf8,
                                    "score": pl.Float64, "matched": pl.UInt32})
    hits = (corpus["postings"].join(terms.select("term_id"), on="term_id", how="inner")
                              .group_by("document_id")
                              .agg(pl.col("weight").sum().alias("score"),
                                   pl.len().alias("matched")))
    return (hits.join(corpus["documents"].select("document_id", "reference.id", "pmid"),
                      on="document_id", how="inner")
                .select("reference.id", "pmid", "score",
                        pl.col("matched").cast(pl.UInt32))
                .sort(["score", "reference.id"], descending=[True, False])
                .head(limit))


#: What a lift counts. ``documents`` asks how many publications report the two together, which is the
#: retrieval question and the one a study-count answer should use. ``occurrences`` weights by how often
#: a token appears **within its own family**, which is the only mode with useful range for a common
#: token: measured on the corpus, ``k:CAS`` is in 614 of 661 documents, so its document-level lift
#: cannot exceed 1.077 however specific it is.
OVER = ("documents", "occurrences")


@dataclass(frozen=True, slots=True)
class Lift:
    """One lift, with every number it rests on.

    Never report ``lift`` alone: the same 3.0 means something different over four units than over four
    hundred, and ``given_units`` is what says which (CLAUDE.md section 5a). ``over`` says whether a
    unit is a document or an occurrence, because the two are different questions.
    """

    term: str
    given: tuple[str, ...]
    over: str
    units: int
    given_units: int
    term_units: int
    both: int
    rate_given: float
    rate_overall: float
    lift: float | None

    def __str__(self) -> str:
        condition = " + ".join(self.given) or "nothing"
        if self.lift is None:
            return (f"{self.term} given {condition}: no lift, "
                    f"{self.given_units} {self.over} carry the condition")
        return (f"{self.term} given {condition}: lift {self.lift:.3f} over {self.over} "
                f"({self.both}/{self.given_units} = {self.rate_given:.4f} conditioned, "
                f"{self.term_units}/{self.units} = {self.rate_overall:.4f} overall)")


def _documents_with(corpus: dict[str, pl.DataFrame], token: str) -> pl.DataFrame:
    """``(document_id, tf)`` for one token. Empty when the vocabulary does not have it."""
    term = corpus["terms"].filter(pl.col("term") == token)
    if term.is_empty():
        return pl.DataFrame(schema={"document_id": pl.UInt32, "tf": pl.UInt32})
    return (corpus["postings"].join(term.select("term_id"), on="term_id", how="inner")
                              .select("document_id", "tf"))


def _group(corpus: dict[str, pl.DataFrame], given: list[str]) -> pl.DataFrame:
    """``(document_id,)`` for the documents carrying every token in ``given``; all when empty.

    Semi-joins rather than set intersection: the frames are already indexed on ``document_id``, and
    ``Expr.is_in`` against a Series of the same dtype is both slower and deprecated in polars.
    """
    ids = corpus["documents"].select("document_id")
    for token in given:
        ids = ids.join(_documents_with(corpus, token).select("document_id"),
                       on="document_id", how="semi")
    return ids.unique()


def lift(corpus: dict[str, pl.DataFrame], term: str, given: list[str] | None = None, *,
         over: str = "documents") -> Lift:
    """How much more often ``term`` appears among documents carrying every token in ``given``.

    ``lift`` is the conditioned rate over the overall rate: 1.0 means the condition tells you nothing,
    above 1.0 that the two go together, below that they avoid each other. ``None`` where nothing
    carries the condition, which is not a lift of zero and must not be averaged as one.

    This is what separates a CDR3 motif from the V gene that templates it: compare
    ``lift(c, "k:CAS", ["e:KRWIILGLNK"])`` against
    ``lift(c, "k:CAS", ["e:KRWIILGLNK", "v:TRBV9"])``. If conditioning on the V gene leaves the lift
    where it was, the k-mer carries the association; if it moves toward 1.0, the V gene did.

    ⚠ **A species condition is a provenance question, not a specificity one.** ``a:HIV-1`` is a real
    and useful axis - "which papers and receptors are about this species" is where most questions
    start - but the group it selects is a *union over pMHCs*: every receptor reported against some
    epitope of that species, under whatever restriction each study used. A lift over it therefore
    describes that group and is not a motif *for* the pathogen, because its members were shown
    different antigens. Condition on ``e:<epitope>``, or on that plus a restriction, when the claim is
    about recognition. ``docs/standards/terminology.md`` has the distinction.

    Measured on the current corpus over occurrences, the answer is neither: `k:CAS` lifts **0.969** on
    HIV-1 documents (28,422 of 739,216 CDR3 3-mer occurrences against 137,751 of 3,473,003 overall),
    so it is very slightly depleted rather than enriched, which is what a germline-encoded motif looks
    like. Holding TRBV9 moves it to 1.006.

    The instrument has range on the same corpus: conditioning on ``e:GILGFVFTL``, ``k:IRS`` is the
    highest-lifting of the 2,342 CDR3 3-mers with 50 or more occurrences at 2.66x, and the RS-bearing
    3-mers as a family sit far above the 1.17x median - the motif that epitope is known for.

    ``over="documents"`` counts publications, which is the right denominator for "who reported this"
    and lets one small study weigh as much as one large one. ``over="occurrences"`` counts token
    instances within the term's own family, which is the right denominator for "how much of this
    antigen's receptor repertoire carries the motif" and is the only mode with range for a token most
    documents contain.
    """
    if over not in OVER:
        raise ValueError(f"over must be one of {OVER}, not {over!r}")
    conditions = list(given or [])
    hits = _documents_with(corpus, term)
    group = _group(corpus, conditions)
    in_group = hits.join(group, on="document_id", how="semi")

    if over == "documents":
        units = corpus["documents"].height
        term_units, both, given_units = hits.height, in_group.height, group.height
    else:
        # Within the term's own family. Summing every family into one denominator would put `k:CAS`
        # over a total that includes epitope k-mers, MHC tokens and V genes, so the rate would move
        # with how many epitopes a paper studied rather than with the motif. Measured, that dilution
        # flattens the whole comparison to within 2 % of 1.0.
        siblings = (corpus["terms"].filter(pl.col("family") == tokens.family_of(term))
                                   .select("term_id"))
        totals = (corpus["postings"].join(siblings, on="term_id", how="inner")
                                    .group_by("document_id").agg(pl.col("tf").sum().alias("total")))
        units = int(totals["total"].sum() or 0)
        term_units = int(hits["tf"].sum() or 0)
        both = int(in_group["tf"].sum() or 0)
        given_units = int(totals.join(group, on="document_id", how="semi")["total"].sum() or 0)

    rate_given = both / given_units if given_units else 0.0
    rate_overall = term_units / units if units else 0.0
    return Lift(
        term=term, given=tuple(conditions), over=over, units=units, given_units=given_units,
        term_units=term_units, both=both, rate_given=rate_given, rate_overall=rate_overall,
        lift=(rate_given / rate_overall) if (given_units and rate_overall > 0) else None,
    )


def lift_family(corpus: dict[str, pl.DataFrame], family: str,
                given: list[str] | None = None, *, over: str = "documents",
                min_units: int = 0) -> pl.DataFrame:
    """Every term of one ``family``, scored in one pass. The batched form of :func:`lift`.

    One row per term with the same numbers :class:`Lift` carries, so nothing has to be reported
    without its ``given_units`` (``CLAUDE.md`` section 5a). ``min_units`` drops terms whose ``both``
    is below it, which is the occurrence floor a comparison across terms needs.

    **Why this exists.** Calling :func:`lift` in a loop recomputes the condition group and the family
    totals once per term, and both are the same every time: scoring the 2,342 CDR3 3-mers that way
    took **50.7 s**, a quarter of the whole test suite, at about 22 ms a term
    (``ROADMAP_local.md`` section 57.1). Reach order rung 1 - one grouped polars expression - and the
    per-term :func:`lift` stays for the one-shot question it is good at.

    ``tests/unit/test_corpus_query.py`` asserts the two agree term for term, so the fast path cannot
    drift from the one the docstring of :func:`lift` documents.
    """
    if over not in OVER:
        raise ValueError(f"over must be one of {OVER}, not {over!r}")
    conditions = list(given or [])
    members = corpus["terms"].filter(pl.col("family") == family).select("term_id", "term")
    group = _group(corpus, conditions)                      # computed ONCE, not once per term
    posts = corpus["postings"].join(members, on="term_id", how="inner")

    if over == "documents":
        units, given_units = corpus["documents"].height, group.height
        per_term = (posts.group_by("term", maintain_order=True)
                    .agg(pl.col("document_id").n_unique().alias("term_units")))
        both = (posts.join(group, on="document_id", how="semi")
                .group_by("term", maintain_order=True)
                .agg(pl.col("document_id").n_unique().alias("both")))
    else:
        # The family's own denominator, for the reason `lift` gives: summing every family would put a
        # CDR3 k-mer over a total including epitope k-mers and V genes.
        totals = posts.group_by("document_id", maintain_order=True).agg(
            pl.col("tf").sum().alias("total"))
        units = int(totals["total"].sum() or 0)
        given_units = int(totals.join(group, on="document_id", how="semi")["total"].sum() or 0)
        per_term = (posts.group_by("term", maintain_order=True)
                    .agg(pl.col("tf").sum().alias("term_units")))
        both = (posts.join(group, on="document_id", how="semi")
                .group_by("term", maintain_order=True).agg(pl.col("tf").sum().alias("both")))

    out = (per_term.join(both, on="term", how="left")
           .with_columns(pl.col("both").fill_null(0))
           .with_columns(
               pl.lit(units).alias("units"), pl.lit(given_units).alias("given_units"),
               (pl.col("both") / given_units if given_units else pl.lit(0.0)).alias("rate_given"),
               (pl.col("term_units") / units if units else pl.lit(0.0)).alias("rate_overall"))
           .with_columns(
               pl.when((pl.lit(given_units) > 0) & (pl.col("rate_overall") > 0))
               .then(pl.col("rate_given") / pl.col("rate_overall"))
               .otherwise(None).alias("lift")))
    # Sorted, with the term as the tiebreak: two terms on the same lift must not swap between runs
    # (hard rule 7).
    return (out.filter(pl.col("both") >= min_units)
            .select("term", "units", "given_units", "term_units", "both", "rate_given",
                    "rate_overall", "lift")
            .sort("lift", "term", descending=[True, False], nulls_last=True))


def refsearch_query(cdr3: str = "", epitope: str = "", *, extra_parameters: str = "",
                    species_to_search: str = "", k: int = tokens.K) -> list[str]:
    """The token list for a ``vdjdb.com/refsearch/`` request, from the fields the client sends.

    The contract, read from ``vdjdb-web``'s own client: every value is a space-joined string,
    ``extra_parameters`` is drawn from ``search_by_antigen`` and ``filter_stop_words``, and
    ``species_to_search`` defaults to the three species the client sends when the user picks none.

    A CDR3 becomes its k-mers rather than one token, because a query is a motif and not a sequence:
    `CAS` is what a user types, and it is a k-mer. An epitope becomes both the whole sequence and its
    k-mers, so an exact epitope outranks a partial match through carrying more matched tokens.
    """
    flags = set(extra_parameters.split())
    out: list[str] = []
    for motif in cdr3.split():
        upper = motif.upper()
        if len(upper) < k:
            continue
        out += [tokens.FAMILIES["cdr3_kmer"].prefix + upper[i:i + k]
                for i in range(len(upper) - k + 1)]
    for epi in epitope.split():
        upper = epi.upper()
        out.append(tokens.FAMILIES["epitope"].prefix + upper)
        out += [tokens.FAMILIES["epitope_kmer"].prefix + upper[i:i + k]
                for i in range(len(upper) - k + 1)]
        if "search_by_antigen" in flags:
            # The client's `search_by_antigen` asks for the words of the antigen as well as its
            # sequence, which is the text family restricted to the epitope's own spelling.
            out.append(tokens.FAMILIES["word"].prefix + upper.lower())
    species = species_to_search.split() or ["HomoSapiens", "MusMusculus", "MacacaMulatta"]
    out += [tokens.FAMILIES["host_species"].prefix + s for s in species]
    if "filter_stop_words" in flags:
        prefix = tokens.FAMILIES["word"].prefix
        out = [t for t in out
               if not (t.startswith(prefix) and t[len(prefix):] in pubmed.STOP_WORDS)]
    return out
