"""Coordinate conversions. The only module that may do this.

Four coordinate spaces meet in this codebase and three of them disagree with AIRR on two axes at
once (origin and closedness), so an ad-hoc `+1` in a call site is indistinguishable from a correct
one until a motif is plotted two residues off:

======================================  ===============================================
VDJdb ``v.end`` / ``j.start``           0-based, **junction** amino-acid space
``vdjtools`` ``Scenario.v_end``         0-based half-open, CDR3-**nucleotide** space
``arda.annotate.dmap``                  1-based **closed**, junction-nucleotide space
AIRR ``*_sequence_start`` / ``_end``    1-based closed, full-**sequence** nucleotide space
======================================  ===============================================

Every function here is total and pure, and every pair has a round-trip test. Nothing composes them
for you: a caller states the conversion it wants, in order, so the reader of that call site can
check it.

Junction is not CDR3. VDJdb's ``cdr3`` column holds the junction: Cys104 through Phe/Trp118, both
anchors included. AIRR's ``junction_aa`` is the same thing, so that mapping is the identity; AIRR's
``cdr3_aa`` and ``arda``'s ``cdr3_aa`` exclude both anchors and are two residues shorter. Conflating
them corrupts every coordinate downstream, with no error (CLAUDE.md, domain conventions).
"""
from __future__ import annotations

import polars as pl

#: Residues the junction has at each end that IMGT CDR3 does not: Cys104 and Phe/Trp118.
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


def nt_to_aa_boundary(i: int) -> int:
    """A V/J boundary in nucleotides -> the same boundary in ``v.end``/``j.start`` space.

    **Ceiling, not floor, and that was fitted rather than reasoned.** ``v.end`` and ``j.start`` are
    documented as "0-based amino acid, junction space" without stating closedness, and the two
    readings differ by one whenever a boundary falls inside a codon - which is most of the time,
    because a V/J boundary is a nucleotide event and nothing aligns it to a codon edge.

    Fitted twice independently, against the external nucleotide truth in
    ``tests/release/test_cdr3fix_accuracy.py`` and against the recombination model: each of these is
    the **count of residues the segment touches**, so a partly-covered codon counts, and the
    conversion is ``(i + 2) // 3``. Under :func:`nt_to_aa`, which floors, the exact-match rates
    against that truth set drop by tens of percent and no head-to-head count moves.

    Here rather than in either caller, because a coordinate conversion in two places is a
    coordinate conversion that will disagree in one of them (``CLAUDE.md``, the four spaces table).
    """
    return (i + CODON - 1) // CODON


def nt_to_aa_boundary_expr(nt: pl.Expr, *, unmapped: int = -1) -> pl.Expr:
    """:func:`nt_to_aa_boundary` over a column, with nulls and negatives as ``unmapped``.

    The same arithmetic as the scalar above and asserted against it in ``tests/unit/test_coords.py``,
    because two spellings of one conversion is the thing that module comment warns about. A column is
    needed because the build converts 187,055 boundaries at once, not one.
    """
    return (pl.when(nt.is_null() | (nt < 0)).then(pl.lit(unmapped, pl.Int64))
            .otherwise(((nt + CODON - 1) // CODON).cast(pl.Int64)))


# -- junction space <-> sequence space ---------------------------------------------------------

def to_sequence(pos: int, junction_start: int) -> int:
    """A position inside the junction -> the same position in the full sequence.

    ``junction_start`` is where the junction begins in that sequence, in the *same* origin and
    units as both. Phase 8 supplies it; there is no nucleotide sequence to index into before then.
    """
    return junction_start + pos


def from_sequence(pos: int, junction_start: int) -> int:
    return pos - junction_start
