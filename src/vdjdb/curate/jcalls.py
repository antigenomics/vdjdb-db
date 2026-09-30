"""J calls the junction itself contradicts (#681).

A junction carries the J's own evidence: its 3' end is templated by the J germline, so a record's J
call can be checked against the sequence the same record reports. Match the junction's 3' end -
**the anchor excluded**, so a corrupt anchor cannot vote on its own diagnosis - against every J
germline of that species and take the longest run.

This is not the question `curate.anchors` asks. That one reads the anchor residue of the segment a
record *names* and asks whether the sequence agrees with it; this one asks whether some other gene
explains the whole 3' end better. Measured on the 2026-09-29 build, the two overlap on **19 of 718**
chains, so 699 are found here and by nothing else.

**Advisory, and re-calling a J is a curator's decision, not this module's.** #681 says so directly:
`arda.cdr3fix` already names a J for a record that carries none, and overwriting a
submitted call with an inferred one is a judgement about a publication. What is missing is the list,
which is what this writes.

The rule is stated rather than tuned, because a threshold picked to make a number look right is not
evidence:

* the best-matching gene must be **unique** - a tie between two genes says the junction cannot
  distinguish them, which is a fact about the germlines and not about the record;
* its run must be at least :data:`MIN_RUN` residues. On a 20-letter alphabet a chance run of 5 is
  about 3e-7 per germline, so ~90 germlines still leave it implausible;
* the called gene's own run must be **strictly shorter**. Equal runs mean the sequence is consistent
  with the call, whatever else it is also consistent with.

Under it, 718 chains over 96 chunks, against 285,545 checkable ones - 99.75 % of J calls are not
contradicted by their own junction.
"""
from __future__ import annotations

from collections import defaultdict
from functools import lru_cache

import polars as pl

from .anchors import _anchors

#: Shortest germline run this will call evidence. See the module docstring for why 5.
MIN_RUN = 5

#: Species with a germline reference, and arda's name for each.
ORGANISMS: dict[str, str] = {
    "HomoSapiens": "human", "MusMusculus": "mouse", "MacacaMulatta": "rhesus_monkey",
}

REPORT_COLUMNS: tuple[str, ...] = (
    "chunk.file", "record_id", "gene", "species", "cdr3",
    "j.segm", "called.run", "best.gene", "best.run",
)


@lru_cache(maxsize=8)
def _suffix_index(species: str) -> tuple[dict[str, frozenset[str]], dict[str, str]]:
    """``(suffix -> the genes ending in it, gene -> its longest templated residues)``.

    Indexed rather than scanned. Comparing each junction against all ~90 J germlines costs 7.2 s over
    the corpus's 180,537 distinct ``(species, cdr3, j)`` keys; looking its own suffixes up, longest
    first, is ~18 dictionary probes and costs about a second (CLAUDE.md hard rule 8, rung 1 - this
    never leaves Python, but it stops being quadratic in the germline table).

    Anchors are stripped from both sides here, once, so no caller can forget to.
    """
    organism = ORGANISMS.get(species)
    if organism is None:
        return {}, {}
    longest: dict[str, str] = {}
    for (segment, name), templated in _anchors(organism).items():
        if segment != "J" or not name.startswith("TR"):
            continue
        gene, body = name.split("*")[0], templated[:-1]
        if len(body) > len(longest.get(gene, "")):
            longest[gene] = body
    index: dict[str, set[str]] = defaultdict(set)
    for gene, body in longest.items():
        for n in range(1, len(body) + 1):
            index[body[-n:]].add(gene)
    return {k: frozenset(v) for k, v in index.items()}, longest


def _run(tail: str, germline: str) -> int:
    """Length of the longest common suffix of two strings."""
    n = 0
    while n < len(tail) and n < len(germline) and tail[-1 - n] == germline[-1 - n]:
        n += 1
    return n


def contradicted(species: str, cdr3: str, call: str) -> tuple[str, int, int] | None:
    """``(best gene, its run, the called gene's run)`` when the junction names another J.

    ``None`` when the call is consistent, when no single gene wins, when the run is too short to be
    evidence, or when the species has no germline reference.
    """
    index, longest = _suffix_index(species)
    if not index or not cdr3 or not call:
        return None
    tail = cdr3[:-1]
    for n in range(min(len(tail), max(len(b) for b in longest.values())), MIN_RUN - 1, -1):
        genes = index.get(tail[-n:])
        if not genes:
            continue
        if len(genes) > 1:
            return None
        best = next(iter(genes))
        called = call.split("*")[0]
        if best == called:
            return None
        return best, n, _run(tail, longest.get(called, ""))
    return None


def report(chains: pl.DataFrame, records: pl.DataFrame) -> pl.DataFrame:
    """One row per chain whose J call some other gene explains better.

    Deduplicated to the distinct ``(species, cdr3, j.segm)`` before the scan and joined back
    (CLAUDE.md rule 4): 285,545 checkable chains are 180,537 distinct keys, and the function is
    deterministic in its arguments.
    """
    joined = (chains.join(records.select("record_id", "species", "chunk.file"),
                          on="record_id", how="left")
                    .filter((pl.col("cdr3") != "") & (pl.col("j.segm") != "")))
    keys = joined.select("species", "cdr3", "j.segm").unique()
    found = [(sp, cdr3, j, *hit)
             for sp, cdr3, j in keys.iter_rows()
             if (hit := contradicted(sp, cdr3, j)) is not None]
    schema = {"species": pl.String, "cdr3": pl.String, "j.segm": pl.String,
              "best.gene": pl.String, "best.run": pl.Int64, "called.run": pl.Int64}
    flagged = pl.DataFrame(found, schema=schema, orient="row")
    if flagged.is_empty():
        return pl.DataFrame(schema={c: (pl.Int64 if c.endswith(".run") else pl.String)
                                    for c in REPORT_COLUMNS})
    return (joined.join(flagged, on=["species", "cdr3", "j.segm"], how="inner")
                  .select(REPORT_COLUMNS)
                  .sort("best.run", "chunk.file", "record_id", descending=[True, False, False]))
