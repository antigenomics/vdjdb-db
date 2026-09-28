"""Which alleles can present each epitope VDJdb records. Issue #372.

**This does not second-guess the allele a publication reported.** The restriction on a record is the
author's finding and the authority for it is the paper. What a presentation predictor adds is the
*other* alleles: given a donor's HLA type, which of VDJdb's epitopes could that donor present,
including pairings no publication has reported. That is a filter for a query, not a correction to a
record, and it is what #372 asks for - "see if our epitopes, linked to, say A*02, can be presented by
other alleles".

The one judgement the pipeline does make about `mhc.a` is whether the allele exists. A blank or
spurious allele is a transcription error rather than a finding, and that check belongs in `vdjdb qc`
against `proofreading/mhc_alleles.tsv.gz`, not here.

**The presentation axis, not the specificity axis.** :meth:`mhcmatch.Store.restriction` ranks how
strongly one allele prefers a peptide *relative to other alleles*, and its own docstring warns it can
band a canonical, widely-shared ligand as "weak" while that allele is the correct restriction -
measured on `FLRGRAYGL` it put `HLA-A*02:01` ahead of `HLA-B*08:01` on a neighbour vote of 0.73
against 0.03, contradicting the paper and 80 references. So the scorer here is
:func:`mhcmatch.predict.build_scorer` with ``background="proteome"``, the NetMHCpan ``%Rank_EL``
analogue, which asks "is this presented at all".

One scorer build serves every epitope and every allele: it is memoised on the store and depends only
on the panel, so the whole table is one pass of ``model.score`` per (epitope, allele) - never one
:func:`mhcmatch.predict.predict_windows` call per pair, which re-resolves the panel and re-derives the
cut every time.

The output is a **committed, reviewed input**, not a build artifact: `mhcmatch` fetches its reference
data from HuggingFace on first use, so scoring inside a build would put a network call in the critical
path and make the answer depend on the day. Same treatment as ``summary/reference_years.tsv`` -
refreshed by its own pull request through ``vdjdb promiscuity``, never written by ``vdjdb build``
(hard rule 9).
"""

from __future__ import annotations

from typing import NamedTuple

#: Written to `proofreading/epitope_promiscuity.tsv`, in this order.
COLUMNS = ("epitope", "mhc_a", "percent_rank", "p_present", "band", "recorded")

#: Class I binding cores. An epitope outside this range is not a class I ligand, so the class I
#: scorer has nothing to say about it and it is skipped rather than scored badly.
LENGTHS = range(8, 12)


class Presented(NamedTuple):
    """One epitope scored against one allele of the panel."""

    epitope: str
    mhc_a: str
    percent_rank: float
    p_present: float
    band: str
    #: 1 when VDJdb records this epitope under this allele, 0 when only the predictor pairs them.
    #: A recorded pair is kept whatever it scores - the author reported it, and a reader filtering by
    #: donor type needs to see the record, not the prediction's opinion of it.
    recorded: int


def presented(recorded: dict[str, set[str]], *, species: str = "human",
              tier: str = "shortlist") -> tuple[list[Presented], list[str]]:
    """Score every epitope in ``recorded`` against the whole class I panel.

    ``recorded`` maps an epitope to the alleles VDJdb records it under, in VDJdb's spelling. Returns
    ``(rows, unmatched)``: the rows a reader can filter by donor allele, and the VDJdb allele
    spellings the panel has no entry for. An unmatched allele is a gap in the panel, not a claim
    about the allele - whether an allele exists is IPD-IMGT/HLA's answer, not a predictor's.

    Rows are every allele that presents the epitope (``band`` above non-binder) plus every allele
    VDJdb records for it. Sorted by epitope then ascending %rank, ties by allele, so the file is
    byte-identical across runs (hard rule 7).

    Network-bound on first call: `mhcmatch` fetches its reference data from HuggingFace. That is why
    this is not in a build.
    """
    import mhcmatch
    from mhcmatch.predict import band_for, build_scorer

    store = mhcmatch.Store.from_pmhc(tier=tier, species=species)
    model, cal, _ = build_scorer(store, "mhc1", background="proteome")
    panel = sorted(store.panel_alleles("mhc1", "all"))

    # VDJdb's spelling and the panel's are not the same vocabulary; `panel_alleles` is the one place
    # that folds them together, so translate once rather than comparing strings per row.
    spelling: dict[str, str] = {}
    unmatched: list[str] = []
    for allele in sorted({a for alleles in recorded.values() for a in alleles}):
        hit = store.panel_alleles("mhc1", [allele])
        if hit:
            spelling[allele] = hit[0]
        else:
            unmatched.append(allele)

    rows: list[Presented] = []
    for epitope in sorted(recorded):
        if len(epitope) not in LENGTHS:
            continue
        mine = {spelling[a] for a in recorded[epitope] if a in spelling}
        scored: list[Presented] = []
        for allele in panel:
            s = model.score(epitope, allele)
            if s == float("-inf"):
                continue
            rank = cal.percent_rank(allele, s)
            if rank != rank:                      # nan: this allele has no background
                continue
            band = band_for(rank, "mhc1")
            if band == "non-binder" and allele not in mine:
                continue
            scored.append(Presented(epitope, allele, round(rank, 3),
                                    round(cal.p_present(allele, s), 4), band,
                                    1 if allele in mine else 0))
        rows.extend(sorted(scored, key=lambda p: (p.percent_rank, p.mhc_a)))
    return rows, unmatched
