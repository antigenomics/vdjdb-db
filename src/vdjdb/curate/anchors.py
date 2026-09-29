"""Junctions that contradict the germline anchors of the V and J they name.

VDJdb's ``cdr3`` is junction space: Cys104 through Phe/Trp118, **both anchors included**. A
submission in IMGT CDR3 space, one carrying framework past an anchor, or one with a mis-read anchor
residue is none of those, and it passes every rule in :mod:`vdjdb.qc.rules` - those check the residue
alphabet and a minimum length, not the ends.

``arda.cdr3fix`` already decides something close to this per chain, and the build already writes the
verdict: ``v.canonical`` and ``j.canonical`` in ``chains``. Measured on the 2026-09-28 build, **990 of
285,989 chains (0.346 %) are flagged non-canonical and all 990 have ``fix.good`` false** - and nothing
read either column. They were written, shipped, and never reported.

Those two columns are not the test used here, because they compare against a fixed Cys / Phe-or-Trp
rather than against the segment's own germline. Measured against the germline, **481 of the 990 are
canonical after all** - the junction ends in the residue its J encodes and only the fixed test
disagrees - and the fixed test separately misses 125 rows whose junction ends in Phe where the
germline says Cys. So the check below reads every chain and keeps arda's verdict beside its own, which
is what makes the disagreement visible instead of averaged away.

**Nothing here fails anything.** A non-canonical junction can be the database's best record of what a
paper published, and the repair below is a proposal with the germline behind it, not a verdict. The
chunk is the data; an edit to it is its own branch with its own issue (``CLAUDE.md``).

The anchor residue is **per segment, read from arda's own table, never hardcoded to F/W**. Mouse
``TRAJ47*01`` is ``HYANKMIC`` and human ``TRAJ35*01`` is ``IGFGNVLHC``, so a junction on either ends
in Cys and an "ends with F or W" test calls a correct record broken: measured, that test flags 125
rows whose junction is canonical and whose anchor table disagrees with it.

Three defects, each with its own repair, told apart by aligning against the germline the record
names:

``under-trimmed``
    framework carried past the anchor. ``YLCSSQEGGYGYTFGSG`` with ``TRBV29-1``/``TRBJ1-2`` is the
    junction plus ``YL`` in front and ``GSG`` behind; arda trims one end and leaves the other, which
    is why these are flagged at all. Repair: trim to the anchors.
``corrupt anchor``
    the anchor residue is mis-read and the rest aligns. ``GASSDTMNTKIL`` and ``WAVRDIYTTAKFIL`` are
    a Cys read as Gly and as Trp. No TCR folds without Cys104, so this is a sequencing or
    transcription error rather than a variant. Repair: substitute the germline residue.
``absent anchor``
    the anchor is missing and the rest aligns. ``ASSNEKLF`` on ``TRBJ1-4`` (``TNEKLFF``) is short a
    Cys in front and an F behind. Repair: add them.

and one where the sequence is right and the **call** is wrong:

``allele mismatch``
    the junction does not match the anchor of the allele the record names, and does match a
    functional sibling allele of the same gene. The repair is then the call, not the sequence.
    Measured: 95 mouse chains name ``TRAJ47``, which resolves to ``*01`` - an ORF allele whose
    templated residues are ``HYANKMIC`` - while every one reads ``DYANKMIF``, which is exactly
    ``TRAJ47*02``, the functional allele. **No record in the corpus reads the ``*01`` signature.**
    Same defect as issue #327, where 66 % of explicit ``TRAJ24*01`` calls carry the ``*02`` motif.

An ORF or pseudogene allele has a non-canonical anchor by definition, and arda records that: of 383 J
entries over four organisms only 14 have ``templated_aa`` not ending in Phe or Trp, and 13 of those 14
are marked ``ORF`` or ``P``. So a junction disagreeing with a non-functional allele is evidence about
the *call*, which is why this class exists and why there is no "the reference is wrong" class.
"""
from __future__ import annotations

from functools import lru_cache

import polars as pl

#: Germline residues that must agree before a repair is proposed. Three is enough to distinguish the
#: three defects on real data and short enough to work on a junction of five residues; the corpus has
#: ``CAMRE``.
MATCH = 3

#: What a functional segment carries at the two anchor positions: Cys104, which the disulphide needs,
#: and Phe or Trp at 118. Used to find the functional sibling allele when the one a record names is an
#: ORF or a pseudogene - not to overrule the table, which marks 13 of its 14 non-canonical J entries
#: ``ORF`` or ``P`` and so already agrees about which alleles are functional.
UNIVERSAL: dict[str, str] = {"V": "C", "J": "FW"}

#: Columns :func:`noncanonical` returns, in order.
COLUMNS: tuple[str, ...] = (
    "chunk.file", "chunk.row", "record_id", "gene", "species", "cdr3", "cdr3.original",
    "v.segm", "j.segm", "v.canonical", "j.canonical", "v.anchor", "j.anchor",
    "defect", "repair", "repair.call",
)


@lru_cache(maxsize=8)
def _anchors(organism: str) -> dict[tuple[str, str], str]:
    """``(segment, allele) -> germline-templated junction residues``, from arda.

    Not a cache of a computed result (rule 9): it is arda's reference table, an input, read once per
    organism instead of once per chain. ``load_anchors`` is itself cached upstream for the same
    reason.

    ``ensure_reference`` first, and it is not optional. Without it ``load_anchors`` returns an empty
    dict rather than raising - the failure documented in :func:`vdjdb.annotate.cdr3fix.ensure_reference`
    - so every junction comes back "no germline" and the check reports nothing while looking like it
    ran. Measured: all 990 flagged chains, before this call was added.
    """
    from arda.cdr3fix import load_anchors

    from ..annotate.cdr3fix import ensure_reference

    ensure_reference()
    table = {k: a.templated_aa for k, a in load_anchors(organism).items() if a.templated_aa}
    if not table:
        raise RuntimeError(f"arda has no anchor table for {organism}: the germline reference is "
                           f"missing or unreadable, so no junction can be checked against it")
    return table


def templated(species: str, segment: str, call: str) -> str | None:
    """The germline-templated junction residues of the segment a record names, or ``None``.

    ``arda.cdr3fix.VDJDB_SPECIES`` is the species translation and the authority - it is keyed
    lowercase, which is the one trap here. An allele-less call falls back to ``*01``: the templated
    residues differ between alleles only where the allele differs, and a record that names no allele
    is not claiming one.
    """
    from arda.cdr3fix import VDJDB_SPECIES

    organism = VDJDB_SPECIES.get((species or "").lower())
    if organism is None or not call:
        return None
    table = _anchors(organism)
    base = call.split("*")[0]
    for name in (call, f"{base}*01", f"{base}*02", f"{base}*03"):
        found = table.get((segment, name))
        if found:
            return found
    return None


@lru_cache(maxsize=4096)
def functional_sibling(species: str, segment: str, call: str, tail: str) -> str | None:
    """A functional allele of the same gene whose anchor is ``tail``, or ``None``.

    Consulted before a sequence repair is proposed, because an ORF or pseudogene allele has a
    non-canonical anchor by definition: a junction disagreeing with one is evidence that the *call* is
    wrong, not that the sequence is. ``tail`` is the junction's last residue for a J, its first for a
    V.

    Only one match counts. Two candidates is not an answer, and taking the first would make the
    proposal depend on dictionary order.
    """
    from arda.cdr3fix import VDJDB_SPECIES

    organism = VDJDB_SPECIES.get((species or "").lower())
    if organism is None or not call:
        return None
    gene = call.split("*")[0]
    matching = sorted(name for (seg, name), residues in _anchors(organism).items()
                      if seg == segment and name.split("*")[0] == gene and name != call
                      and (residues[-1] if segment == "J" else residues[0]) == tail)
    return matching[0] if len(matching) == 1 else None


def classify(cdr3: str, species: str, v: str, j: str) -> tuple[str, str | None, str | None]:
    """``(defect, proposed junction, proposed call)``, or ``("ok", None, None)``.

    Both ends are considered, and a sequence can be wrong at both: ``YFCASSYWVGDTDTQYFGPG`` carries
    framework in front and behind. The V end is repaired first and the J end is then read off the
    result, so the two repairs compose instead of fighting over indices.

    A proposed *call* takes precedence over a proposed sequence: where the junction matches a
    functional sibling allele of the gene the record names, the sequence is right and the allele
    resolution is not, so rewriting the sequence would destroy the evidence for the real defect - the
    reasoning ``MAX_REPLACE = 0`` already applies in :mod:`vdjdb.annotate.cdr3fix`.
    """
    seen: list[str] = []
    out = cdr3
    call_fix: str | None = None
    for end, call, segment in (("V", v, "V"), ("J", j, "J")):
        germline = templated(species, segment, call)
        if not germline:
            # No germline to read: the call is absent, or it is a name the reference does not have.
            # The universal anchor is all that is left, and it is exactly what the retired build's
            # `is_qq_seq_biologically_valid` used for every row. Without this the check declines
            # silently, and 118 chains whose junction ends in neither the anchor nor anything a
            # germline could justify were reported by the old build and by nothing in this one.
            # No repair is proposed: there is no germline to propose from.
            tail = out[:1] if end == "V" else out[-1:]
            seen.append(f"{end} unanchored" if tail not in UNIVERSAL[end]
                        else f"no {end} germline")
            continue
        anchor = germline[0] if end == "V" else germline[-1]
        tail = out[:1] if end == "V" else out[-1:]
        if tail != anchor and tail in UNIVERSAL[end] and anchor not in UNIVERSAL[end]:
            sibling = functional_sibling(species, segment, call, tail)
            if sibling is not None:
                seen.append(f"{end} allele mismatch")
                call_fix = sibling
                continue
        defect, out = (_v_end(out, germline) if end == "V" else _j_end(out, germline))
        seen.append(f"{end} {defect}" if defect != "ok" else "")
    # An end whose segment is not in the reference is *unchecked*, not defective, and saying
    # otherwise buries the finding: 3,730 chains name a V the reference does not have - a missing or
    # misspelled call, which is a different defect with its own rule - against 2,787 that are
    # genuinely short a J anchor. `v.anchor` / `j.anchor` are empty for an unchecked end, so the
    # coverage is still in the report.
    named = [s for s in seen if s and not s.startswith("no ")]
    if not named:
        return "ok", None, None
    return ", ".join(named), (out if out != cdr3 else None), call_fix


def _overlap(a: str, b: str, *, suffix: bool) -> int:
    """Length of the longest common prefix, or suffix, of two strings."""
    if suffix:
        a, b = a[::-1], b[::-1]
    n = 0
    while n < len(a) and n < len(b) and a[n] == b[n]:
        n += 1
    return n


def _v_end(cdr3: str, germline: str) -> tuple[str, str]:
    """Classify the Cys104 end against the V's germline-templated residues.

    ``germline[0]`` is Cys104 itself, so the question is which germline position the junction starts
    at: 0 is correct, 1 means the anchor is missing.
    """
    anchor = germline[0]
    at_anchor = _overlap(cdr3, germline, suffix=False)
    one_past = _overlap(cdr3, germline[1:], suffix=False)
    if cdr3[:1] == anchor and at_anchor >= one_past:
        return "ok", cdr3
    # Framework in front of the anchor. The test is that trimming aligns *strictly better* than not
    # trimming, rather than that it aligns deeply: a V contributes three to five residues to the
    # junction before the N region takes over - `TRAV27*01` is `CAG` - so "deeply" is not available.
    # `YLCSSQEGGYGYTFGSG` on `TRBV29-1*01` (`CSVE`) aligns 0 residues as given and 2 from the Cys.
    for i in range(1, min(6, len(cdr3))):
        if cdr3[i] == anchor and _overlap(cdr3[i:], germline, suffix=False) > at_anchor:
            return "under-trimmed", cdr3[i:]
    if one_past > at_anchor and one_past >= MATCH:
        return "absent anchor", anchor + cdr3
    if _overlap(cdr3[1:], germline[1:], suffix=False) >= MATCH:
        return "corrupt anchor", anchor + cdr3[1:]
    return "unexplained", cdr3


def _j_end(cdr3: str, germline: str) -> tuple[str, str]:
    """Classify the Phe/Trp118 end against the J's germline-templated residues.

    The terminal residue on its own decides nothing here, which is the trap the V end does not have:
    ``TRBJ1-4*01`` is ``TNEKLFF``, so a junction ending ``NEKLF`` ends in the anchor residue and is
    still one short of the anchor. What separates them is how deep the alignment runs - ``CASSNEKLF``
    matches six residues of ``TNEKLF`` and one of ``TNEKLFF``, ``CASSNEKLFF`` the other way round.

    A tie is read as correct. This is a non-blocking alert, and a false alarm costs a curator more
    than a missed one: the sequence can be nibbled back past anything the germline could confirm.
    """
    anchor = germline[-1]
    at_anchor = _overlap(cdr3, germline, suffix=True)
    one_short = _overlap(cdr3, germline[:-1], suffix=True)
    if at_anchor >= one_short and cdr3[-1:] == anchor:
        return "ok", cdr3
    for i in range(1, min(6, len(cdr3))):                    # framework behind the anchor
        head = cdr3[:len(cdr3) - i]
        if cdr3[-1 - i] == anchor and _overlap(head, germline, suffix=True) > at_anchor:
            return "under-trimmed", head
    if one_short > at_anchor and one_short >= MATCH:
        return "absent anchor", cdr3 + anchor
    if _overlap(cdr3[:-1], germline[:-1], suffix=True) >= MATCH:
        return "corrupt anchor", cdr3[:-1] + anchor
    return "unexplained", cdr3


#: ``(chunk column suffix, locus)``, as :mod:`vdjdb.assemble.tables` spells it.
_GENES = (("alpha", "TRA"), ("beta", "TRB"))


def chains_of(master: pl.DataFrame) -> pl.DataFrame:
    """The master table's paired alpha/beta columns as one row per chain.

    Taken from ``master`` rather than from the ``chains`` table so that ``vdjdb submission`` can run
    this on a pull request: that command stops after ``build_master`` on purpose, because the
    annotation stages cost ten times as much and change nothing a curator is deciding.
    """
    parts = []
    for gene, tag in _GENES:
        have = [c for c in (f"cdr3.{gene}", f"v.{gene}", f"j.{gene}") if c in master.columns]
        if len(have) < 3:
            continue
        parts.append(master.select(
            pl.col("record_id") if "record_id" in master.columns else pl.lit("").alias("record_id"),
            pl.lit(tag).alias("gene"),
            pl.col("species"),
            *[pl.col(c) for c in ("chunk.file", "chunk.row") if c in master.columns],
            pl.col(f"cdr3.{gene}").alias("cdr3"),
            # The submitted sequence, which is what a chunk edit would change. Absent when nothing
            # needed fixing, in which case the shipped value *is* the submitted one.
            (pl.col(f"__cdr3old.{gene}") if f"__cdr3old.{gene}" in master.columns
             else pl.lit(None, pl.Utf8)).cast(pl.Utf8).alias("cdr3.original"),
            pl.col(f"v.{gene}").alias("v.segm"),
            pl.col(f"j.{gene}").alias("j.segm"),
            *[(pl.col(f"__{k}.{gene}") if f"__{k}.{gene}" in master.columns
               else pl.lit(None, pl.Boolean)).alias(f"{k[0]}.canonical")
              for k in ("vcanon", "jcanon")],
        ))
    return (pl.concat(parts, how="vertical").filter(pl.col("cdr3") != "")
            if parts else master.head(0))


def noncanonical(master: pl.DataFrame) -> pl.DataFrame:
    """Every chain whose junction contradicts the germline anchors of the V or J it names.

    Every chain is read, not only the ones ``v.canonical`` / ``j.canonical`` flag: that pair compares
    against a fixed Cys / Phe-or-Trp and so both over- and under-reports against the germline (481 and
    125 rows respectively, measured). Both columns are carried through so the disagreement is on the
    face of the report.

    The repair is proposed against the submitted sequence rather than the shipped one, because that is
    what an edit to the chunk would change and ``arda.cdr3fix`` may have repaired one end already -
    ``YLCSSQEGGYGYTFGSG`` ships as ``YLCSSQEGGYGYTF``, the framework trimmed behind the anchor and
    kept in front of it.

    1.4 s over 285,989 chains: pure string comparison against a table read once per organism.
    """
    chains = chains_of(master)
    if chains.is_empty():
        return pl.DataFrame({c: [] for c in COLUMNS}).cast(pl.Utf8)
    rows = chains.to_dicts()
    verdict = [classify(r["cdr3.original"] or r["cdr3"], r["species"] or "",
                        r["v.segm"] or "", r["j.segm"] or "") for r in rows]
    return (chains
            .with_columns(
                pl.Series("defect", [d for d, _, _ in verdict]),
                pl.Series("repair", [p for _, p, _ in verdict], dtype=pl.Utf8),
                pl.Series("repair.call", [c for _, _, c in verdict], dtype=pl.Utf8),
                pl.col("cdr3.original").fill_null(pl.col("cdr3")),
                pl.Series("v.anchor", [(templated(r["species"] or "", "V", r["v.segm"] or "")
                                        or "")[:1] for r in rows]),
                pl.Series("j.anchor", [(templated(r["species"] or "", "J", r["j.segm"] or "")
                                        or "")[-1:] for r in rows]))
            .filter(pl.col("defect") != "ok")
            .select([c for c in COLUMNS if c in set(chains.columns)
                     | {"defect", "repair", "repair.call", "v.anchor", "j.anchor"}])
            .sort("defect", "cdr3.original", "record_id", "gene"))


def report(flagged: pl.DataFrame) -> str:
    """The markdown alert. Empty string when nothing is flagged, so a caller can skip the section."""
    if flagged.is_empty():
        return ""
    repairable = flagged.filter(pl.col("repair").is_not_null())
    recall = flagged.filter(pl.col("repair.call").is_not_null())
    out = ["#### Alert: junctions that contradict their own V or J germline", "",
           f"**This blocks nothing.** {flagged.height:,} chain(s) carry a `cdr3` whose first or last "
           "residue is not the anchor the named segment encodes. VDJdb's `cdr3` is junction space - "
           "Cys104 through Phe/Trp118, both included - so a sequence in IMGT CDR3 space, one carrying "
           "framework past an anchor, or one with a mis-read anchor is not one.", "",
           "The anchor is read from the germline of the segment the record names, not assumed to be "
           "F or W: mouse `TRAJ47*01` is `HYANKMIC`, so a junction on it ends in Cys.", "",
           "| Defect | Chains | Repair proposed |", "|---|---:|---:|"]
    for row in (flagged.group_by("defect").agg(pl.len().alias("n"),
                                               pl.col("repair").is_not_null().sum().alias("fixable"))
                .sort("n", descending=True).iter_rows(named=True)):
        out.append(f"| {row['defect']} | {row['n']:,} | {row['fixable']:,} |")
    out += ["", f"{repairable.height:,} of {flagged.height:,} have a repair the germline supports. "
            "Each is a proposal for the chunk cell, against the submitted sequence rather than the "
            "repaired one, because `arda.cdr3fix` may have already fixed one end.", ""]
    if not recall.is_empty():
        out += [f"For {recall.height:,} of them the **call** is what is wrong, not the sequence: the "
                "junction matches a functional sibling allele of the gene the record names, where the "
                "allele it resolved to is an ORF or a pseudogene. Repairing the sequence there would "
                "destroy the evidence for the real defect.", "",
                "| Junction | Called | Should be | Records |", "|---|---|---|---:|"]
        grouped = (recall.group_by("v.segm", "j.segm", "repair.call")
                   .agg(pl.len().alias("n"), pl.col("cdr3").first().alias("example"))
                   .sort("n", descending=True).head(8))
        for row in grouped.iter_rows(named=True):
            gene = row["repair.call"].split("*")[0]
            called = row["j.segm"] if row["j.segm"].startswith(gene) else row["v.segm"]
            out.append(f"| `{row['example']}` | `{called}` | `{row['repair.call']}` | "
                       f"{row['n']:,} |")
        out.append("")
    if not repairable.is_empty():
        out += ["| Submitted | Proposed | V | J | Defect |", "|---|---|---|---|---|"]
        for row in repairable.head(20).iter_rows(named=True):
            out.append(f"| `{row['cdr3.original']}` | `{row['repair']}` | `{row['v.segm']}` | "
                       f"`{row['j.segm']}` | {row['defect']} |")
        if repairable.height > 20:
            out.append(f"\n...and {repairable.height - 20:,} more, in `out/reports/anchors.tsv`.")
        out.append("")
    return "\n".join(out)
