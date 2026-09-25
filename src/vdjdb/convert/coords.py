"""Coordinate conversions. **The only module allowed to do this.**

Four coordinate spaces meet in this codebase and three of them disagree with AIRR on two axes at
once (origin and closedness), so an ad-hoc `+1` in a call site is indistinguishable from a correct
one until someone plots a motif two residues off:

======================================  ===============================================
VDJdb ``v.end`` / ``j.start``           0-based, **junction** amino-acid space
``vdjtools`` ``Scenario.v_end``         0-based half-open, CDR3-**nucleotide** space
``arda.annotate.dmap``                  1-based **closed**, junction-nucleotide space
AIRR ``*_sequence_start`` / ``_end``    1-based closed, full-**sequence** nucleotide space
======================================  ===============================================

Every function here is total, pure and its own inverse's counterpart, and every pair has a
round-trip test. Nothing composes them for you: a caller states the conversion it wants, in order,
so the reader of that call site can check it.

**Junction is not CDR3.** VDJdb's ``cdr3`` column holds the *junction*: Cys104 through Phe/Trp118,
**both anchors included**. AIRR's ``junction_aa`` is the same thing, so that mapping is the identity;
AIRR's ``cdr3_aa`` and ``arda``'s ``cdr3_aa`` exclude both anchors and are two residues shorter.
Conflating them silently corrupts every coordinate downstream of it (CLAUDE.md, domain conventions).
"""
from __future__ import annotations

import polars as pl

#: Residues the junction carries at each end that IMGT CDR3 does not: Cys104 and Phe/Trp118.
ANCHOR = 1

#: Nucleotides per amino acid. Named because `* 3` in a coordinate expression reads as a typo.
CODON = 3


# -- junction <-> CDR3 (amino acids) -----------------------------------------------------------

def cdr3_from_junction(junction: str) -> str:
    """Strip the two anchor residues. ``CASSIRSSYEQYF`` -> ``ASSIRSSYEQY``.

    A junction shorter than three residues has no CDR3 at all, and returns empty rather than a
    negative slice.
    """
    return junction[ANCHOR:-ANCHOR] if len(junction) > 2 * ANCHOR else ""


def junction_from_cdr3(cdr3: str, first: str = "C", last: str = "F") -> str:
    """Re-add the anchors. Not the inverse of :func:`cdr3_from_junction` unless they are supplied:
    the CDR3 does not record which residues were removed, and ``last`` is Phe *or* Trp."""
    return f"{first}{cdr3}{last}" if cdr3 else ""


def cdr3_aa(col: str = "cdr3") -> pl.Expr:
    """:func:`cdr3_from_junction` as a polars expression, for a whole column at once."""
    return (pl.when(pl.col(col).str.len_chars() > 2 * ANCHOR)
            .then(pl.col(col).str.slice(ANCHOR, pl.col(col).str.len_chars() - 2 * ANCHOR))
            .otherwise(pl.lit(""))
            .alias("cdr3_aa"))


# -- origin and closedness --------------------------------------------------------------------

def to_one_based(i: int) -> int:
    """0-based index -> 1-based. Undefined for a sentinel; check for -1 before calling."""
    return i + 1


def to_zero_based(i: int) -> int:
    return i - 1


def close_end(end_exclusive: int) -> int:
    """Half-open end -> closed end. The last included position, not one past it."""
    return end_exclusive - 1


def open_end(end_inclusive: int) -> int:
    return end_inclusive + 1


# -- amino acids <-> nucleotides ---------------------------------------------------------------

def aa_to_nt(i: int, *, one_based: bool = False) -> int:
    """Amino-acid index -> the **first** nucleotide of that codon, in the same origin.

    0-based: residue 0 starts at nt 0. 1-based: residue 1 starts at nt 1.
    """
    return i * CODON if not one_based else (i - 1) * CODON + 1


def nt_to_aa(i: int, *, one_based: bool = False) -> int:
    """Nucleotide index -> the amino-acid index of its codon, in the same origin."""
    return i // CODON if not one_based else (i - 1) // CODON + 1


# -- junction space <-> sequence space ---------------------------------------------------------

def to_sequence(pos: int, junction_start: int) -> int:
    """A position inside the junction -> the same position in the full sequence.

    ``junction_start`` is where the junction begins in that sequence, in the *same* origin and
    units as both. Phase 8 supplies it; there is no nucleotide sequence to index into before then.
    """
    return junction_start + pos


def from_sequence(pos: int, junction_start: int) -> int:
    return pos - junction_start
