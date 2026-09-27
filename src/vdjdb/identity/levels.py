"""Derived identifiers for the levels above and below a curated record.

A record is a receptor against a presented peptide, and both halves recur across records, so each
half needs a name of its own. Four levels do, and all four key on content that has no separate
existence: there is no such thing as a clonotype apart from its species, gene, CDR3, V and J. So
their ids are hashes of their own keys rather than allocated numbers, and three properties follow
that a counter cannot give (``ROADMAP.md`` section 10.3):

* the order ``chunks/`` is read in cannot reach them;
* adding a chunk names that chunk's new clonotypes and changes no existing id;
* removing a chunk retires only the ids nothing else supported.

:mod:`vdjdb.identity.ids` holds the one exception. ``record_id`` is allocated against a registry
because its purpose is to survive a content change: a curator fixing a CDR3 typo has to keep the
record's id, and no hash of the content can do that.
"""
from __future__ import annotations

import hashlib
from collections.abc import Iterable
from dataclasses import dataclass

import polars as pl

#: Unit separator, so a field containing a tab cannot forge a key boundary. The same choice as
#: :mod:`vdjdb.identity.ids`, because the two produce ids for one database.
SEP = "\x1f"

#: Hex digits of the digest an id keeps: 16, which is 64 bits. At the current 187,935 clonotypes the
#: chance of one collision is around 1 in 10^9, and :func:`duplicate_ids` fails a build rather than
#: letting one corrupt a join. A size decision, not a security one.
ID_HEX = 16


@dataclass(frozen=True, slots=True)
class Level:
    """One derived level: its name, its id prefix, and the columns it keys on."""

    name: str
    prefix: str
    key: tuple[str, ...]

    @property
    def column(self) -> str:
        return f"{self.name}_id"


#: One receptor chain. Every record reporting the same chain shares one id: it is the level motif
#: evidence attaches at, and the level the independent-study support count is measured on
#: (:mod:`vdjdb.assemble.evidence`).
CLONOTYPE = Level("clonotype", "CT", ("species", "gene", "cdr3", "v.segm", "j.segm"))

#: One alpha/beta pair. Keys on its two clonotypes rather than on columns of its own, so it is
#: assigned by :func:`clone_ids` and not by :func:`derive`.
CLONE = Level("clone", "CX", ("clonotype_id",))

#: One presented peptide. ``mhc.class`` is deliberately absent: the alleles determine it, and it adds
#: nothing to the key (2,364 distinct either way on the current build).
PMHC = Level("pmhc", "PM", ("antigen.epitope", "mhc.a", "mhc.b"))

#: The peptide alone, across every allele presenting it.
EPITOPE = Level("epitope", "EP", ("antigen.epitope",))

LEVELS: tuple[Level, ...] = (CLONOTYPE, CLONE, PMHC, EPITOPE)

LEVELS_BY_NAME: dict[str, Level] = {lv.name: lv for lv in LEVELS}

#: Every id prefix in use, including the allocated one, so :func:`level_of` covers all five.
PREFIXES: dict[str, str] = {lv.prefix: lv.name for lv in LEVELS} | {"VDJDB": "record"}


def join_key(parts: Iterable[str]) -> str:
    """The key string an id is hashed from. Empty stands in for a missing field (hard rule 6)."""
    return SEP.join("" if p is None else str(p) for p in parts)


def hash_id(prefix: str, key: str) -> str:
    """One id from an already-joined key.

    sha256 rather than :meth:`polars.Expr.hash`: polars does not specify its hash output across
    versions, and this id ships to consumers, so a dependency upgrade must not be able to renumber
    the database with nothing failing anywhere.
    """
    return prefix + hashlib.sha256(key.encode("utf-8")).hexdigest()[:ID_HEX]


def level_of(identifier: str) -> str | None:
    """The level an id belongs to, read from its prefix. ``None`` if the prefix is not one of ours."""
    for prefix, name in PREFIXES.items():
        if identifier.startswith(prefix):
            return name
    return None


def _assign(df: pl.DataFrame, prefix: str, column: str) -> pl.DataFrame:
    """Hash ``__k`` on its distinct values, join back, drop ``__k``.

    Hashing the distinct set rather than every row is hard rule 4 and not a cache: the hash is a
    pure function of its argument, so computing it once per distinct key inside one build cannot
    change an answer. Measured on the current build, 187,935 distinct clonotype keys cost 83 ms to
    hash and 7 ms to join back.
    """
    uniq = df["__k"].unique().sort()
    table = pl.DataFrame({"__k": uniq, column: [hash_id(prefix, k) for k in uniq]},
                         schema={"__k": pl.Utf8, column: pl.Utf8})
    # `maintain_order="left"` because callers hand in a sorted frame and expect it back sorted;
    # a join does not otherwise promise to keep the left order.
    return df.join(table, on="__k", how="left", maintain_order="left").drop("__k")


def derive(df: pl.DataFrame, level: Level, *, column: str | None = None) -> pl.DataFrame:
    """``df`` with ``level``'s id column appended.

    Every key column has to be present. A missing one is a schema error rather than a missing value,
    so it raises here instead of hashing a column of empties into one id for the whole table.
    """
    missing = [c for c in level.key if c not in df.columns]
    if missing:
        raise KeyError(f"{level.name} id needs {missing}, which {sorted(df.columns)[:6]}... lacks")
    keyed = df.with_columns(
        pl.concat_str([pl.col(c).cast(pl.Utf8).fill_null("") for c in level.key],
                      separator=SEP).alias("__k"))
    return _assign(keyed, level.prefix, column or level.column)


def clone_ids(chains: pl.DataFrame) -> pl.DataFrame:
    """``(record_id, clone_id)`` for every record reporting two chains of different ``gene``.

    A clone is a receptor, so it keys on its two clonotypes and on nothing else, and the same
    alpha/beta pair reported by two papers is one clone in two records. The pair is sorted before
    hashing, so the id does not depend on which chain the reader met first.

    A record with one chain gets no row here, and :func:`vdjdb.assemble.tables.build_chains` fills
    the column with an empty string rather than a null (hard rule 6). Having no clone is the curated
    state and not missing data: 99,459 of 192,753 records carry one chain, against 93,294 carrying
    two, which hold 82,266 distinct clones.
    """
    paired = (chains.group_by("record_id")
                    .agg(pl.col("gene").n_unique().alias("__genes"),
                         pl.col("clonotype_id").cast(pl.Utf8).sort().alias("__pair"))
                    .filter(pl.col("__genes") == 2)
                    .with_columns(pl.col("__pair").list.join(SEP).alias("__k"))
                    .select("record_id", "__k"))
    return _assign(paired, CLONE.prefix, CLONE.column).sort("record_id")


def duplicate_ids(df: pl.DataFrame, column: str, key: Iterable[str]) -> pl.DataFrame:
    """Ids in ``column`` that more than one distinct ``key`` maps to. Empty is the passing state.

    This is invariant 1 of ``ROADMAP.md`` section 10.5, and it is a collision check rather than a
    uniqueness check: an id repeating across rows of the same key is how a derived id is meant to
    behave.
    """
    return (df.select(column, *key).unique()
              .group_by(column).len().filter(pl.col("len") > 1).sort(column))


def recomputed_mismatches(df: pl.DataFrame, level: Level, *,
                          column: str | None = None) -> pl.DataFrame:
    """Rows whose stored id is not the hash of the key on that same row.

    Invariant 2. The id is recomputed here rather than trusted, so a join that silently attached the
    wrong id, or a hash implementation that changed under us, fails instead of shipping.
    """
    name = column or level.column
    both = derive(df.select(*level.key, pl.col(name).alias("__have")), level, column="__want")
    return both.filter(pl.col("__want") != pl.col("__have")).unique(maintain_order=True)
