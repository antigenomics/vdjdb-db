"""A V or J call for the chains that have none (#462, #658).

The README lets a chunk leave ``v.alpha`` or ``j.beta`` blank, and 3,126 of 191,103 distinct
``(species, cdr3, v, j)`` keys do: 644 with no V, 2,943 with no J. Something has to name one, because
``arda.cdr3fix`` repairs a junction against a *named* germline and never proposes one - given a blank
V it reports ``FailedBadSegment``, and the record then fails the legacy build's "a CDR3 needs a V and
a J" filter and is dropped from ``vdjdb.txt``.

Two sources answer, in this order, and neither is the k-mer scanner over ``res/segments.txt`` that
used to (#658):

1. **The recombination model** (``vdjtools.model.infer_nt_batch`` with the blank side marginalised).
   It searches nucleotide histories over the whole germline set and returns the V and J of the most
   likely one, so it uses the junction's interior and not only its ends.
2. **arda's germline anchor table**, where the model declines. Each allele's ``templated_aa`` is the
   junction residues its germline encodes, so the call is the allele agreeing with the most of the
   junction's own end. This is the same shape of answer the k-mer scanner gave, against the
   authoritative IMGT reference rather than a 2023 by-product of an IMGT import, and it is also the
   only source for species no model covers.

Measured against the curated calls with the true call hidden, 1,200 distinct keys per species and
locus, agreement at gene level:

=============  ==========  ==========  ==============  ===============
species/locus  k-mer scan  model       arda germline   population
=============  ==========  ==========  ==============  ===============
human TRB J    95.9 %      97.5 %      96.2 %          1,420 blank
human TRA J    79.3 %      95.8 %      93.6 %          1,259 blank
mouse TRB J    -           -           95.8 %            160 blank
mouse TRA J    -           -           97.2 %             30 blank
human TRB V    0.0 %       23.8 %      15.8 %            241 blank
human TRA V    0.1 %       50.1 %      41.0 %            150 blank
mouse TRB V    -           -           35.8 %            231 blank
mouse TRA V    -           -           10.3 %             22 blank
=============  ==========  ==========  ==============  ===============

So the model leads on both coordinates where it exists, which is why it goes first, and the germline
table beats the scanner it replaces on the one locus where the scanner was not simply broken
(human TRA J, 93.6 % against 79.3 %).

The V column was never a comparison: ``Cdr3Fixer.guess_id`` put ``return ""`` inside the five-prime
loop, so it tried one prefix length and gave up - 3 non-empty guesses in 4,000 sequences, where the
J branch has the same statement correctly in a ``for...else``. VDJdb's V guesser never worked.

**That gap is why only the J proposal ships.** ``j.segm`` carries it, as it has in every release;
``v.segm`` stays blank where the curator left it blank, and the V proposal is reported as
``v.inferred`` and nowhere else. :func:`vdjdb.annotate.cdr3fix.markup` holds the reasoning and the
measurement that settled it.

``v.inferred`` and ``j.inferred`` are non-empty exactly where the record named no segment and a source
proposed one, so a proposal never sits beside a curated call. Read the V accuracy above before using
``v.inferred``: recovered from the junction alone it is right about a quarter of the time for human
TRB, because TRBV contributes only a few residues to it. The J side is reliable.
"""
from __future__ import annotations

from functools import lru_cache

import polars as pl

from .junction import MODELS

#: VDJdb's chain word -> the locus its model and anchor table are keyed under.
LOCI: dict[str, str] = {"alpha": "TRA", "beta": "TRB"}

#: Residues of germline agreement before a name means anything. Two residues of a junction's end
#: match most J alleles of a locus, so a shorter run names whichever sorted first rather than
#: whichever fits.
_MIN_RUN = 3


@lru_cache(maxsize=32)
def _anchors(organism: str, locus: str, segment: str) -> tuple[tuple[str, str], ...]:
    """``(allele, templated_aa)`` for one segment of one locus, best candidate first.

    Functional before ORF and pseudogene, then the longest templated region, then the name: a
    non-functional allele cannot be the segment of an expressed receptor, and where two alleles
    agree with the junction equally far the tie has to break the same way on every host (rule 7).
    """
    from arda.cdr3fix import load_anchors

    rows = [(allele, a.templated_aa, a.functionality)
            for (seg, allele), a in load_anchors(organism).items()
            if seg == segment and a.locus == locus and a.templated_aa]
    rows.sort(key=lambda r: (r[2] != "F", -len(r[1]), r[0]))
    return tuple((allele, templated) for allele, templated, _f in rows)


def _germline_call(cdr3: str, anchors: tuple[tuple[str, str], ...], *, five_prime: bool) -> str:
    """The allele whose germline-templated residues agree with the most of ``cdr3``'s own end.

    ``""`` when nothing reaches :data:`_MIN_RUN`. The V anchor is read forward from Cys104 and the J
    anchor backward from Phe/Trp118, which is the direction each germline actually templates.

    ponytail: a Python loop over ~90 alleles per row. It runs only where the model declined - a few
    hundred rows of a build - so it is 0.1 s and not the bottleneck. If it ever runs on the whole
    corpus, index the anchors by their terminal k-mer first.
    """
    best, longest = "", 0
    for allele, templated in anchors:
        n, limit = 0, min(len(templated), len(cdr3))
        if five_prime:
            while n < limit and templated[n] == cdr3[n]:
                n += 1
        else:
            while n < limit and templated[-1 - n] == cdr3[-1 - n]:
                n += 1
        if n > longest:
            best, longest = allele, n
    return best if longest >= _MIN_RUN else ""


def _locus(keys: pl.DataFrame, gene: str | None) -> pl.Expr:
    """Which locus each row's model and anchor table belong to, from whichever source knows.

    The caller's chain word, else a ``gene`` column (``chains.gene`` already holds the locus), else
    whichever of the two calls is present. A record with neither call and no locus reads as beta,
    which is the majority and is what the retired guesser did.
    """
    if gene:
        return pl.lit(LOCI[gene])
    if "gene" in keys.columns:
        return pl.col("gene")
    return (pl.when(pl.concat_str("v", "j").str.starts_with("TRA")).then(pl.lit("TRA"))
              .otherwise(pl.lit("TRB")))


def propose(keys: pl.DataFrame, gene: str | None = None) -> pl.DataFrame:
    """``keys`` (``species``, ``cdr3``, ``v``, ``j``) plus ``__gv`` / ``__gj``.

    ``v`` and ``j`` are left untouched: they are the join key back to the table, and overwriting them
    made the lookup miss and fanned ``vdjdb_full.txt`` out by 3,266 rows. A call the record already
    carries is never replaced, so a blank is the only cell either source can fill.

    Rows with both calls already present pass straight through, so the model is asked about the 3,126
    keys that need it rather than all 191,103. Row order is not preserved; every caller joins or
    sorts on its own key.
    """
    from arda.cdr3fix import VDJDB_SPECIES
    from vdjtools.model import infer_nt_batch, load_bundled

    from .junction import _resolver

    given = keys.with_columns(pl.col("v").alias("__gv"), pl.col("j").alias("__gj"))
    incomplete = (pl.col("__gv") == "") | (pl.col("__gj") == "")
    todo = given.filter(incomplete)
    if todo.is_empty():
        return given

    parts = [given.filter(~incomplete)]
    # Sorted so the concatenation order never depends on group iteration order (rule 7).
    for (species, loc), group in sorted(
            todo.with_columns(_locus(todo, gene).alias("__locus"))
                .group_by("species", "__locus", maintain_order=True),
            key=lambda kv: kv[0]):
        gv, gj = pl.col("v"), pl.col("j")
        if species in MODELS:
            source, organism = MODELS[species]
            model = load_bundled(loc, source, organism=organism)
            vmap = _resolver(model, "genes_v", "v_allele")
            jmap = _resolver(model, "genes_j", "j_allele")
            # One batched call, nothing around it (hard rule 3). The blank side goes in as None so
            # the model chooses it rather than echoing ours back.
            got = infer_nt_batch(model, group["cdr3"].to_list(),
                                 v=[vmap.get(x) if x else None for x in group["v"]],
                                 j=[jmap.get(x) if x else None for x in group["j"]])
            gv = pl.when(gv == "").then(got["v_call"].fill_null("")).otherwise(gv)
            gj = pl.when(gj == "").then(got["j_call"].fill_null("")).otherwise(gj)
        group = group.with_columns(gv.alias("__gv"), gj.alias("__gj"))

        # arda's germline table answers what the model declined, and is the only source for a
        # species no model covers.
        organism = VDJDB_SPECIES.get(species.lower())
        still_blank = (pl.col("__gv") == "") | (pl.col("__gj") == "")
        if organism is not None and group.filter(still_blank).height:
            av, aj = _anchors(organism, loc, "V"), _anchors(organism, loc, "J")
            group = group.with_columns(
                pl.Series("__gv", [g or _germline_call(c, av, five_prime=True)
                                   for c, g in zip(group["cdr3"], group["__gv"], strict=True)],
                          dtype=pl.Utf8),
                pl.Series("__gj", [g or _germline_call(c, aj, five_prime=False)
                                   for c, g in zip(group["cdr3"], group["__gj"], strict=True)],
                          dtype=pl.Utf8))
        parts.append(group.drop("__locus"))
    return pl.concat(parts, how="vertical")
