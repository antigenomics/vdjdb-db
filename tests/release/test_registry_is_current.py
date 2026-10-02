"""The committed record registry against the corpus it is supposed to describe (#672).

`registry/records.tsv` is an input: `vdjdb build` reads it and only `vdjdb identity update` writes it.
So it can go stale, and a stale registry is worse than none - it reports amendments nobody made and,
where a key change is ambiguous, retires published identifiers and mints new ones (vdjdb-db#638).

`vdjdb identity check` asserts the seven identifier invariants on the *tables*, all of which hold on a
build from a stale registry. Nothing compared the registry against the corpus until this, which is how
`52cb4e2` left 592 records carrying an amendment for a separator: it gave 591 chain-calls a real IMGT
name, `CHUNK_DEDUP_KEY` carries `v.beta` and `j.beta`, and harmonisation runs *before*
`add_record_ids`.

`CLAUDE.md` asks a **chunk** branch to run `vdjdb identity update`. This is the check that makes the
wider rule enforceable rather than remembered: any branch that moves a value inside `CHUNK_DEDUP_KEY`
refreshes the registry, whether or not it edits `chunks/`.

Marked `release` because it reads the whole corpus, which is a ten-second `build_master` rather than a
unit test's fixture.
"""
from __future__ import annotations

import pytest

from vdjdb.assemble.master import REGISTRY
from vdjdb.config import Paths
from vdjdb.identity import IdentityRegistry, reconcile

pytestmark = pytest.mark.release


@pytest.fixture(scope="module")
def reconciled():
    """The corpus reconciled against the committed registry, without writing anything."""
    from vdjdb.assemble.master import harmonise_all
    from vdjdb.io.chunks import merge_repeated_references, read_chunks

    path = Paths.discover().root / REGISTRY
    if not path.exists():
        pytest.skip(f"no committed registry at {path}")
    # **`harmonise_all`, not a hand-written list of the passes.** This fixture used to name the four
    # it knew about, and adding `harmonise_vocabulary` as a fifth (#637) made it reconcile the
    # committed registry against a frame the build does not produce - 170 amendments reported in the
    # wrong direction, on a registry that was correct. One sequence, two callers.
    df, _ = merge_repeated_references(read_chunks(None))
    df, _ = harmonise_all(df)
    _, _, report = reconcile(df, IdentityRegistry.load(path), release="test")
    return report


def test_the_committed_registry_needs_no_amendment(reconciled):
    """An amendment here is a key that moved since the registry was written, so the file is stale.

    Run `uv run vdjdb identity update` and commit `registry/records.tsv` on the branch that moved it,
    with the message naming what moved and how many records.
    """
    assert not reconciled.amended, (
        f"{len(reconciled.amended)} records' natural keys have moved since the registry was written; "
        f"first five: {reconciled.amended[:5]}")


def test_the_committed_registry_retires_and_adds_nothing(reconciled):
    """Retirement and allocation are what a stale registry costs. A record that cannot be matched
    loses its published `record_id`, and identifiers ship to consumers.
    """
    assert not reconciled.retired, (
        f"{len(reconciled.retired)} records would retire against the current corpus: "
        f"{reconciled.retired[:5]}")
    assert not reconciled.added, (
        f"{len(reconciled.added)} records have no id in the registry: {reconciled.added[:5]}")
