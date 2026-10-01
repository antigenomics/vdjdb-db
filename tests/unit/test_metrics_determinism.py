"""No metric may depend on the order of the rows it is handed.

CLAUDE.md hard rule 7: same inputs, same bytes, in another process, on another host, at another core
count. A metric that reads a tie off the incoming row order breaks that without raising, and it
breaks it *selectively* -- purity survived because it sums counts, while precision and the F1 built
on it moved, so the failure looked like a real quality change rather than an instrument fault.

Measured before the fix, human TRB cohort, do-nothing partition, shuffling nothing but the row order:
precision 0.9458, 0.9409, 0.9407, 0.9458. Two CI runs of the same build disagreed at 0.9459 against
0.9407, which is how it was found (`ROADMAP_local.md` §53).

These tests are deliberately built on hand-written frames with **engineered ties**: the real corpus
happens to contain them, but a test that relies on that would stop testing the moment a chunk merges.
"""
from __future__ import annotations

import pandas as pd
import pytest

from vdjdb.emit.vdjdb3 import read_table
from vdjdb.validate import motif_bench as mb
from vdjdb.validate.metrics_lib import binominal_test

#: A cluster split exactly evenly between two epitopes. Whichever one wins becomes `label_cluster`
#: and every one of precision, recall, f1 and ami is computed against it.
TIED = pd.DataFrame({
    "antigen.epitope": ["AAA", "AAA", "BBB", "BBB", "CCC", "CCC", "CCC"],
    "cluster": ["c1", "c1", "c1", "c1", "c2", "c2", "c2"],
})


def test_a_tied_cluster_takes_the_same_epitope_whatever_the_row_order() -> None:
    """The defect, at its smallest: `c1` is 2 AAA and 2 BBB, so `fraction_matched` ties at 0.5."""
    winners = set()
    for order in ([0, 1, 2, 3, 4, 5, 6], [2, 3, 0, 1, 4, 5, 6], [6, 5, 4, 3, 2, 1, 0]):
        res = binominal_test(TIED.iloc[order].reset_index(drop=True), "cluster",
                             "antigen.epitope", compute_pvalue=False)
        winners.add(res.set_index("cluster").loc["c1", "antigen.epitope"])
    assert winners == {"AAA"}, (
        f"a tied cluster took {sorted(winners)} depending on row order; the tie-break must be a "
        f"property of the data, not of the order")


def test_the_tie_break_prefers_the_larger_matched_count_before_the_name() -> None:
    """`fraction_matched` first, then `count_matched`, then the name -- so a genuine majority always
    wins and only a cluster tied on both counts falls through to alphabetical."""
    df = pd.DataFrame({
        "antigen.epitope": ["ZZZ", "ZZZ", "ZZZ", "AAA"],
        "cluster": ["c1"] * 4,
    })
    res = binominal_test(df, "cluster", "antigen.epitope", compute_pvalue=False)
    assert res.set_index("cluster").loc["c1", "antigen.epitope"] == "ZZZ"


def test_every_scored_axis_is_invariant_to_row_order_on_the_real_cohort() -> None:
    """The end-to-end version, on the corpus, with the partition that has the most ties.

    **TRB and two shuffles, not both chains and four.** This exists to catch a *cohort* whose metric is
    order-dependent; the engineered-tie tests above prove the logic, and the four-by-two version cost
    27 s of a 186 s suite for nothing the two-by-one does not say (§57.1). TRB is the chain the CI
    failure appeared on and the one with more tied clusters.

    Skipped without built tables; `build.yml` has them, and this is the shape the CI failure took.
    """
    gene = "TRB"
    from pathlib import Path

    t = Path("out/tables")
    if not (t / "chains.parquet").exists():
        pytest.skip("no built tables; run `vdjdb build --out out/`")
    chains = read_table(t, "chains")
    records = read_table(t, "records")
    cohort = mb.cohort(chains, records, gene=gene)
    assigned = mb.assign(cohort, mb.trivial_members(cohort))
    base = mb.score(assigned)
    for seed in (0, 1):
        got = mb.score(assigned.sample(fraction=1.0, shuffle=True, seed=seed))
        moved = {k: (round(base[k], 6), round(got[k], 6)) for k in base
                 if isinstance(base[k], float) and abs(base[k] - got[k]) > 1e-9}
        assert not moved, f"{gene} moved under row shuffle {seed}: {moved} (before, after)"
