"""Assemble the master table: read, patch, fix, score, hash.

This is ``runBuidDatabase.py``'s main flow. The legacy version spent nearly all of its 344 s in
seven ``master_table.T.apply(...)`` calls -- transpose the frame, then call a Python function per
row -- plus an ``iterrows()`` over ~192k rows in the score factory.

Everything vectorises except the CDR3 fixer, which is per-sequence. That one is called on the
distinct ``(species, cdr3, v, j)`` set and joined back (CLAUDE.md rule 4): the functions are
deterministic in their arguments, so deduplicating cannot change a value, and the key set is
roughly half the row count.
"""
from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import polars as pl

from ..config import Paths
from ..curate.nomenclature import (
    disambiguate_alleles,
    harmonise_mhc,
    harmonise_references,
    harmonise_segments,
)
from ..curate.patch import apply_antigen_patch
from ..io.chunks import read_chunks
from ..schema import ALL_COLUMNS, FULL_COLUMNS
from ..score.confidence import add_score

#: Concatenated to form ``TCR_hash``, the join key to the structure store.
HASH_FIELDS: tuple[str, ...] = (
    "cdr3.alpha", "v.alpha", "j.alpha", "cdr3.beta", "v.beta", "j.beta",
    "mhc.a", "mhc.b", "antigen.epitope",
)

#: A hash is produced only when both chains and the pMHC are present. The other three of
#: :data:`HASH_FIELDS` may be empty and still hash.
#:
#: In the legacy this fell out of ``NaN`` propagating through ``+``. A chain's V or J is set by the
#: fixer, which returns ``""`` when it cannot name one -- an empty string, which concatenates fine --
#: whereas a chain with no CDR3 at all is never fixed and stays ``NaN``, which poisons the
#: concatenation. So "V could not be identified" still hashes and "there is no alpha chain" does not.
#: Treating an empty V as missing dropped 326 hashes that ship today.
HASH_REQUIRED: tuple[str, ...] = ("cdr3.alpha", "cdr3.beta", "mhc.a", "mhc.b", "antigen.epitope")


#: The members of the legacy ``cdr3fix`` JSON blob, flattened into columns. Order matters only
#: when the blob is reassembled for the legacy export; here it is documentation.
FIX_FIELDS: tuple[tuple[str, str, type], ...] = (
    ("cdr3", "__cdr3", str), ("cdr3_old", "__cdr3old", str),
    ("fixNeeded", "__fixneeded", bool), ("good", "__good", bool),
    ("jCanonical", "__jcanon", bool), ("jFixType", "__jfix", str),
    ("jId", "__j", str), ("jStart", "__jstart", int),
    ("vCanonical", "__vcanon", bool), ("vEnd", "__vend", int),
    ("vFixType", "__vfix", str), ("vId", "__v", str),
)

_FIX_DTYPES = {str: pl.Utf8, bool: pl.Boolean, int: pl.Int64}

#: What the markup engine would have called, kept beside what ships. The engine resolves the allele
#: it aligned against, which is evidence about the sequence; the shipped call is what the curator
#: reported, harmonised to IMGT. Keeping both is the only way a reader can tell a curation decision
#: from a markup one, and it is what `chains` reports as `v.segm.arda` / `j.segm.arda`.
CALL_FIELDS: tuple[tuple[str, str], ...] = (("__varda", "v.segm.arda"),
                                            ("__jarda", "j.segm.arda"))



def fix_cdr3(df: pl.DataFrame, engine: str = "arda") -> pl.DataFrame:
    """Repair the CDR3 and locate the V and J germline parts, per chain.

    Rewrites ``cdr3.*``, ``v.*`` and ``j.*``, and adds one column per :data:`FIX_FIELDS` member
    (``__vend.alpha``, ``__jfix.beta``, ...). The legacy ``cdr3fix`` JSON blob is not produced
    here: a blob is not a variable, and reassembling one is the legacy exporter's job.

    Two engines. ``arda`` is the default and what the shipped build uses; ``legacy`` is the
    vendored k-mer scanner, kept so the swap stays measurable and reversible.

    * ``arda`` -- ``arda.cdr3fix``, arda-mapper >= 2.30.1. Agrees with the legacy on 99.91 % of
      repaired alpha sequences and leads on all four coverage measures: it gains 3,781 alpha and
      1,851 beta V-end mappings and 5,182 alpha and 1,904 beta J-start mappings, against 311 /
      3,936 and 7 / 138 lost. The 3,936 beta V-ends it declines are calls that name a family with
      several functional genes, an ambiguity group, or an allele with no shipped anchor -- a
      curation question, not a markup failure.
    * ``legacy`` -- the vendored k-mer scanner, what every release up to 2026-06-03 shipped.

    The coordinates are also more accurate, not merely more numerous, which is the part coverage
    cannot show. Against external nucleotide truth -- the 8,334 VDJdb human TRB records whose
    ``(species, cdr3, v, j)`` key appears exactly once in ``isalgo/airr_control``, where the
    junction was read with its nucleotides present -- ``j.start`` is exact on 7,987 records for the
    legacy scanner and **8,164** for arda 2.30.1, and arda over-extends by two residues or more on
    **1** record against the scanner's 0. arda 2.29.0 scored 8,041 and over-extended on 138, which
    is why the swap waited for `arda#132`.

    Both are called once per distinct ``(species, cdr3, v, j)`` -- 191,447 keys against 192,753
    rows -- and joined back. Both are deterministic in their arguments, so deduplicating cannot
    change a value (CLAUDE.md rule 4).
    """
    run = _markup_arda if engine == "arda" else _markup_legacy
    for gene in ("alpha", "beta"):
        cdr3, v, j = f"cdr3.{gene}", f"v.{gene}", f"j.{gene}"
        keys = (df.filter(pl.col(cdr3) != "")
                  .select("species", pl.col(cdr3).alias("cdr3"),
                          pl.col(v).alias("v"), pl.col(j).alias("j"))
                  .unique(maintain_order=True)
                  .sort("species", "cdr3", "v", "j"))   # sorted: the join order must not vary
        lookup = run(keys, gene).rename({"cdr3": cdr3, "v": v, "j": j})
        df = (
            df.join(lookup, on=["species", cdr3, v, j], how="left")
            .with_columns(
                # A chain with no CDR3 loses its V and J too. The legacy assigns the fixer result
                # unconditionally -- `x.vId if x else None` -- so a record annotated `v.alpha =
                # TRAV2` but with no alpha CDR3 has that annotation wiped. It is a data loss, and
                # the new format keeps those calls; but it must stay here, because the score
                # signature includes v.alpha, so wiping it merges records the released scores are
                # computed over as one group.
                pl.col("__cdr3").fill_null("").alias(cdr3),
                pl.col("__v").fill_null("").alias(v),
                pl.col("__j").fill_null("").alias(j),
                *(pl.col(tmp).alias(f"{tmp}.{gene}") for _, tmp, _ in FIX_FIELDS),
                *(pl.col(tmp).fill_null("").alias(f"{tmp}.{gene}") for tmp, _ in CALL_FIELDS),
            )
            .drop(*(tmp for _, tmp, _ in FIX_FIELDS), *(tmp for tmp, _ in CALL_FIELDS))
        )
    return df


def _markup_arda(keys: pl.DataFrame, gene: str) -> pl.DataFrame:
    from ..annotate.cdr3fix import markup

    return markup(keys, gene)


def _markup_legacy(keys: pl.DataFrame, gene: str) -> pl.DataFrame:
    """The vendored k-mer scanner. Kept only to attribute the swap; deleted once that is frozen."""
    from ..annotate._legacy_fixer import Cdr3Fixer

    res = Paths.discover().res
    fx = Cdr3Fixer(str(res / "segments.txt"), str(res / "segments.aaparts.txt"))
    out: list[dict] = []
    for species, seq, vid, jid in keys.iter_rows():
        vid = vid or fx.guess_id(seq, species, gene, True) or ""
        jid = jid or fx.guess_id(seq, species, gene, False) or ""
        # Always strings: the legacy relied on guess_id having filled the blank, and
        # `"".split(",")` yields `[""]`, which `fix` treats as "no segment given".
        out.append(fx.fix_both(seq, vid, jid, species).results_to_dict())
    return keys.with_columns(
        *(pl.Series(tmp, [r[key] for r in out], dtype=_FIX_DTYPES[ty])
          for key, tmp, ty in FIX_FIELDS),
        # The k-mer scanner names no allele of its own, so it proposes nothing. The columns exist
        # either way: a table's schema must not depend on which engine ran.
        *(pl.lit("").alias(tmp) for tmp, _ in CALL_FIELDS),
    )


def add_tcr_hash(df: pl.DataFrame) -> pl.DataFrame:
    """sha256 over :data:`HASH_FIELDS`, empty unless every :data:`HASH_REQUIRED` field is present."""
    complete = pl.all_horizontal(*[pl.col(c) != "" for c in HASH_REQUIRED])
    return df.with_columns(
        pl.when(complete)
        .then(pl.concat_str(HASH_FIELDS).map_elements(_sha256, return_dtype=pl.Utf8))
        .otherwise(pl.lit(""))
        .alias("TCR_hash")
    )


def _sha256(s: str) -> str:
    import hashlib
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def build_master(paths: Iterable[Path] | None = None,
                 registry: Path | None = None, *, engine: str = "arda") -> pl.DataFrame:
    """The master table: one row per curated record, fixed, scored, hashed and identified.

    Wide (paired alpha/beta columns) because that is the shape the chunks are written in. It is an
    intermediate: :mod:`vdjdb.assemble.tables` pivots it into the tidy tables everything else
    reads.
    """
    df = read_chunks(paths)
    # As submitted, before any harmonisation touches it. `chains` reports these as
    # `v.segm.submitted` / `j.segm.submitted` / `d.segm.submitted`: a reader comparing them with the
    # shipped call sees exactly what the build decided, which is the difference between a curation
    # record and a black box. Costs three string columns per chain.
    df = df.with_columns(
        *(pl.col(c).alias(f"__sub.{c}") for c in
          ("v.alpha", "j.alpha", "v.beta", "j.beta", "d.beta") if c in df.columns))
    df = apply_antigen_patch(df)
    # IMGT spelling before identity: a record is the same record whether the curator wrote
    # `TRAV14` or `TRAV14/DV4`, so harmonising afterwards would mint a new id for a rename.
    df, _ = harmonise_segments(df)
    # Where two alleles differ inside the junction, the sequence is evidence and the submitted call
    # is not (#327). Runs after the spelling pass so it sees IMGT names.
    df, _ = disambiguate_alleles(df)
    # MHC spelling, the allele that does not exist (#467), and the class-II chain order.
    df, _ = harmonise_mhc(df)
    df, _ = harmonise_references(df)
    # Identity is assigned on what the publications reported, before any repair. Afterwards, CDR3
    # fixing would have merged 215 pairs of records the publications reported separately -- two
    # trimmed sequences repaired to the same full one are still two observations.
    df = add_record_ids(df, registry)
    df = fix_cdr3(df, engine)
    df = add_score(df)
    return add_tcr_hash(df)


def add_record_ids(df: pl.DataFrame, registry: Path | None = None) -> pl.DataFrame:
    """Attach a stable ``record_id`` to every row, reconciled against the committed registry.

    An id survives a content change -- a curator fixing a typo amends a record rather than deleting
    one and creating another -- so external references, structure links and accumulated evidence
    outlive curation.
    """
    from ..identity.ids import IdentityRegistry, reconcile

    path = registry or (Paths.discover().root / "registry" / "records.tsv")
    reg = IdentityRegistry.load(path)
    identified, _, _ = reconcile(df, reg, release="dev")
    return identified


def master_for_output(df: pl.DataFrame) -> pl.DataFrame:
    """Only the 35 released columns, in order."""
    assert set(ALL_COLUMNS) <= set(df.columns)
    return df.select(FULL_COLUMNS)
