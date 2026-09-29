"""Whether the recorded MHC could present the recorded epitope at all (ROADMAP phase 9e).

`proofreading/mhc_alleles.tsv.gz` answers whether IPD-IMGT/HLA lists a name. That is a weaker
question than it looks: a name can be in the registry and still be unusable, because what every
presentation model reasons over is the 34 residues lining the binding groove, and a name with no
pseudosequence has none. `mhcmatch` carries those pseudosequences - 20,082 class I keys and 11,048
class II - bundled in the package and loaded in 0.01 s, so the check is offline and deterministic and
belongs in the build (CLAUDE.md hard rule 9). The model-based check, which scores the pair rather
than resolving the name, is phase 9e's fourth and stays in `vdjdb promiscuity`, where a HuggingFace
fetch is allowed.

Four checks, none of them a gate:

* **resolution** - does the call reach a pseudosequence key? `exact` when the name is one, `nearest`
  when `mhcmatch` repairs it by prefix, `none` when nothing matches. Only `none` is a finding:
  `nearest` is overwhelmingly an allele *group* (`HLA-A*02`, which the specification allows and
  10,779 records use) being completed to its first member, which is `mhcmatch` guessing rather than
  the record being wrong. It is carried as a column so a reader can see which scores rest on a guess.
* **the allele's class against `mhc.class`** - the key comes out of the class I FASTA or the class II
  one, so a class II molecule filed as `MHCI` is visible without consulting anything else.
* **the epitope's length against `mhc.class`**, in one direction only. A class I groove is closed at
  both ends, so a 14-mer on `HLA-A*02:01` is a question worth asking; a class II record carrying a
  9-mer is not, because a paper reporting the eluted peptide's core rather than the whole 15-to-25-mer
  is doing something normal. The asymmetry is the difference between 23 pairs worth reading and 41
  pairs of which 18 carry 7,735 records of 9-mer cores.
* **self-consistency** - one MHC molecule filed under both `MHCI` and `MHCII` somewhere in the
  corpus. This needs no authority at all, and it is what catches `H2-IAb` on 77 records of
  `QVYSLIRPNENPAH` where the same molecule is `MHCII` on 537.

**Advisory.** A presentation model is evidence about a pair, never authority over a publication, and
that holds for the name checks too: the corpus contains real alleles `mhcmatch` has no groove for.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from functools import lru_cache

import polars as pl

#: Which bundled FASTA a record's class should draw on. `mhcmatch` names them by locus, not by the
#: `MHCI` / `MHCII` vocabulary VDJdb uses.
FASTA: dict[str, str] = {"MHCI": "mhc1", "MHCII": "mhc2"}

#: The columns of ``out/reports/presentation.tsv``.
REPORT_COLUMNS: tuple[str, ...] = (
    "antigen.epitope", "antigen.species", "mhc.a", "mhc.b", "mhc.class",
    "mhc.key", "mhc.resolution", "mhc.class.inferred", "finding", "records", "references",
)


@lru_cache(maxsize=1)
def _pseudo() -> dict[str, dict[str, str]]:
    from mhcmatch import pseudoseq

    return {cls: pseudoseq.load_pseudo(cls) for cls in ("mhc1", "mhc2")}


def resolve(mhc_a: str, mhc_b: str, mhc_class: str) -> tuple[str, str]:
    """``(pseudosequence key, resolution)`` for one MHC call. Resolution is exact/nearest/none.

    Class II is keyed on the chain pair and class I on ``mhc.a`` alone, which is `mhcmatch`'s own
    split: DR is keyed by its beta chain because DRA is monomorphic, and DP/DQ by the pair.

    ⚠ Both chains are trimmed to two fields first. `mhcmatch`'s `resolve_allele` trims and
    `class2_key` does not, so `HLA-DRB1*11:01:02` builds the key `DRB1_110102` while
    `HLA-DRB1*11:01` builds `DRB1_1101`, and only the second is in the FASTA - same molecule, same
    groove. 29 of the corpus's 354 class II pairs resolve only after the trim. Reported as
    `antigenomics/mhcmatch#3`; remove the trim here when that ships.
    """
    from mhcmatch import pseudoseq

    cls = FASTA.get(mhc_class)
    if cls is None:
        return "", "none"
    keys = _pseudo()[cls]
    if mhc_class == "MHCII":
        if not mhc_b:
            return "", "none"
        key = pseudoseq.class2_key(pseudoseq.trim_allele(mhc_a), pseudoseq.trim_allele(mhc_b))
        if key in keys:
            return key, "exact"
    elif not mhc_a:
        return "", "none"
    try:
        key, exact = pseudoseq.resolve_allele(mhc_a if mhc_class == "MHCI" else key, cls)
    except Exception:  # a name mhcmatch cannot parse is unresolved, not a crash
        return "", "none"
    if key is None:
        return "", "none"
    return key, "exact" if exact else "nearest"


def _class_of_key(key: str) -> str:
    """Which FASTA carries this key, as a VDJdb class name. Empty when neither does."""
    for name, cls in (("MHCI", "mhc1"), ("MHCII", "mhc2")):
        if key in _pseudo()[cls]:
            return name
    return ""


#: Longest peptide a class I groove is asked to hold before this reports it. `mhcmatch.store`'s own
#: length rule is 11, and the corpus has 1,966 class I pairs at or below it against 23 above.
MHCI_MAX_LENGTH = 11


def report(restriction: pl.DataFrame) -> pl.DataFrame:
    """One row per `(epitope, MHC)` pair that fails at least one of the four checks.

    ``restriction`` is the built table: one row per distinct pair, with the record and reference
    counts already aggregated, so this runs over ~2,300 rows rather than ~192,000.
    """
    # How many pairs file each molecule under each class, for the self-consistency check. Keyed on
    # the chains as written, because that is what a curator would have to reconcile, and counted
    # rather than collected so the finding can name which side is the outlier - the answer to "which
    # of these 21 rows is wrong" is the whole of what a reader wants from it.
    classes: dict[tuple[str, str], Counter[str]] = defaultdict(Counter)
    for pair in restriction.iter_rows(named=True):
        classes[(str(pair["mhc.a"] or ""), str(pair["mhc.b"] or ""))][str(pair["mhc.class"] or "")] += 1

    rows: list[dict[str, object]] = []
    for pair in restriction.iter_rows(named=True):
        declared = str(pair["mhc.class"] or "")
        mhc_a, mhc_b = str(pair["mhc.a"] or ""), str(pair["mhc.b"] or "")
        key, resolution = resolve(mhc_a, mhc_b, declared)
        allele_class = _class_of_key(key)
        epitope = str(pair["antigen.epitope"] or "")
        findings = []
        if resolution == "none":
            findings.append("no pseudosequence")
        if allele_class and allele_class != declared:
            findings.append(f"allele is {allele_class}")
        if declared == "MHCI" and len(epitope) > MHCI_MAX_LENGTH:
            findings.append(f"epitope is {len(epitope)} residues")
        filed = classes[(mhc_a, mhc_b)]
        if len(filed) > 1:
            total = sum(filed.values())
            other = ", ".join(f"{c} on {n} of {total}" for c, n in sorted(filed.items())
                              if c != declared)
            findings.append(f"molecule is also {other}")
        if not findings:
            continue
        rows.append({
            "antigen.epitope": epitope, "antigen.species": pair.get("antigen.species", ""),
            "mhc.a": mhc_a, "mhc.b": mhc_b, "mhc.class": declared,
            "mhc.key": key, "mhc.resolution": resolution, "mhc.class.inferred": allele_class,
            "finding": "; ".join(findings),
            "records": pair.get("records", 0), "references": pair.get("references", 0),
        })
    schema = {c: (pl.Int64 if c in ("records", "references") else pl.String)
              for c in REPORT_COLUMNS}
    return (pl.DataFrame(rows, schema=schema)
            .sort("records", "antigen.epitope", "mhc.a", descending=[True, False, False]))


def summarise(flagged: pl.DataFrame) -> pl.DataFrame:
    """The report counted per finding, which is what the build echoes and CI pins."""
    return (flagged.group_by("finding")
                   .agg(pl.len().alias("pairs"), pl.col("records").sum().alias("records"))
                   .sort("records", "finding", descending=[True, False]))
