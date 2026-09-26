"""Numbers pasted into ``docs/clustering.md`` must still be the numbers that were measured.

Covers every table row that names a configuration the scorecard knows -- sections 4.3, 4.5 and 9.2,
155 value cells at the time of writing. Section 7's summary table is keyed on the algorithm rather
than the configuration and is not covered; that shape would need its own parser.

The scorecard tables in that document are produced by ``docs/tuning/report.py`` and pasted in by
hand, which is a step a human can get wrong -- and did: a retention of 0.7698 was transcribed as
0.7699, and before that five values in section 4.3 were wrong enough to change what the surrounding
prose claimed. Both inputs are committed, so this check is hermetic: it re-reads
``docs/tuning/scorecard.tsv`` and compares every value cell at the precision the document itself
displays.

It deliberately does **not** re-measure anything (CLAUDE.md 0b). A sweep is the only thing that may
write the scorecard; this only asserts the prose agrees with it.
"""
from __future__ import annotations

import re
from pathlib import Path

import polars as pl
import pytest

DOC = Path("docs/clustering.md")
SCORECARD = Path("docs/tuning/scorecard.tsv")

#: Table header text -> scorecard column. Only value columns; the key columns are parsed separately.
COLUMNS = {"lift": "lift", "f1": "f1", "q": "q", "`q`": "q", "p": "p", "purity": "purity",
           "precision": "precision", "retention": "retention", "epitopes": "epitopes",
           "perc_med": "perc_med", "percolation": "perc_med", "cids": "cids"}

#: The reference rows as section 9.2 names them, mapped onto their scorecard ``algo``.
REFS = {"legacy release": "legacy", "shipped tcrnet": "shipped-tcrnet",
        "shipped tcremp": "shipped-tcremp"}


def _cells(line: str) -> list[str]:
    return [c.strip() for c in line.strip().strip("|").split("|")]


def _key(cells: list[str], gene: str | None) -> tuple[str, str, str] | None:
    """``(gene, algo, config)`` for a body row, or ``None`` if the row is not a scorecard row.

    Position-independent: the configuration is column 1 in section 4.3, column 2 in section 9.2 and
    column 4 in section 4.5, so every leading cell is tried rather than a fixed one.
    """
    algo = config = None
    for cell in cells:
        lead = cell.replace("**", "").replace("`", "").strip()
        if lead in ("TRA", "TRB"):
            gene = lead
        elif lead.lower() in REFS:
            algo, config = REFS[lead.lower()], "--"
        elif m := re.fullmatch(r"(eom|leaf) (mcs\d+ ms\d+)", lead):
            algo, config = f"hdbscan-{m[1]}", m[2]
        elif m := re.fullmatch(r"(hybrid|hybrid-recruited|hybrid-len-\w+|lumbermark)"
                               r" (mcs\d+ M\d+)", lead):
            algo, config = m[1], m[2]
        elif m := re.fullmatch(r"dbscan (coef [\d.]+)", lead):
            algo, config = "dbscan", m[1]
        if algo is not None:
            break
    return (gene, algo, config) if gene and algo else None


def _rows():
    """Every ``(key, header, cells)`` triple in the document that names a scorecard configuration.

    Tables are grouped properly -- a run of consecutive ``|`` lines whose second line is the
    separator -- because a body row that names no configuration must not be mistaken for a header.
    """
    gene, table = None, []
    lines = [*DOC.read_text().splitlines(), ""]
    for line in lines:
        if line.startswith("|"):
            table.append(_cells(line))
            continue
        if len(table) > 2 and set("".join(table[1])) <= set("-: "):
            header = [c.lower().strip("*` ") for c in table[0]]
            for cells in table[2:]:
                if (k := _key(cells, gene)) is not None:
                    yield k, header, cells
        table = []
        if line.startswith("**TRA**"):
            gene = "TRA"
        elif line.startswith("**TRB**"):
            gene = "TRB"


def test_scorecard_tables_in_clustering_md_match_the_measurement():
    if not SCORECARD.exists():
        pytest.skip(f"{SCORECARD} not present")
    sc = pl.read_csv(SCORECARD, separator="\t", infer_schema_length=None)
    idx = {(r["gene"], r["algo"], r["config"]): r for r in sc.iter_rows(named=True)}

    checked, wrong = 0, []
    for (gene, algo, config), header, cells in _rows():
        row = idx.get((gene, algo, config))
        assert row is not None, f"{gene} {algo} {config} is in the document but not the scorecard"
        for head, cell in zip(header, cells, strict=False):
            col = COLUMNS.get(head)
            text = cell.replace("**", "").replace(",", "")
            if col is None or not re.fullmatch(r"-?\d+(\.\d+)?", text):
                continue
            dp = len(text.split(".")[1]) if "." in text else 0
            checked += 1
            if round(float(row[col]), dp) != float(text):
                wrong.append(f"{gene} {algo} {config} {col}: doc {text}, measured {row[col]}")
    assert not wrong, "\n".join(wrong)
    assert checked > 100, f"only {checked} values checked -- the table parser stopped matching"
