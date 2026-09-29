"""The legacy proofreading, transcribed. This is the oracle the new rules are measured against.

Two tiers ran in the retired pandas build, and both are here:

* **Per chunk** - ``py_src/ChunkQC.py``: ``check_exist``, ``check_header``, the 13-entry
  ``validators`` dict, the three emptiness checks and the within-chunk duplicate scan.
* **Per master table** - ``py_src/runBuidDatabase.py`` lines 126-152: ``gene_match_check``,
  ``alleles_match_check`` and ``is_qq_seq_biologically_valid``, run after CDR3 repair and written to
  ``vdjdb_full_gene_broken.txt``, ``vdjdb_full_allele_broken.txt`` and
  ``vdjdb_full_cdr3aa_broken.txt``.

Recover either from git: ``git show 2026-06-03:py_src/ChunkQC.py``.

**This is a transcription and not a reimplementation.** Every predicate keeps the legacy form,
including the parts that are wrong, because a corrected oracle cannot answer the question it exists
for: does the new build still catch what the old one caught. Each deliberate defect is named below.

One translation, and only one. Legacy read chunks with ``na_values=['']`` so an empty cell arrived
as ``NaN`` and every validator guarded on ``pd.isnull``. Here the missing marker is the empty string
(``CLAUDE.md`` hard rule 6), so ``pd.isnull(x)`` is transcribed as ``x == ""``. The two select the
same cells; what changes is that three legacy validators which would have raised ``TypeError`` or
``AttributeError`` on a ``NaN`` instead return an answer. Those three are marked ``would have
crashed``, and the parity harness records the answer rather than the crash.

Legacy defects preserved on purpose:

* ``is_MHC_valid`` accepts anything whose first three characters are not ``HLA``, so every murine
  spelling passes unchecked. ``vdjdb.qc.rules`` keeps this, for the reason its docstring gives.
* ``is_MHC_valid`` has no missing guard at all, so on a genuinely absent ``mhc.a`` the legacy build
  raised ``TypeError`` from inside ``re.match``. ``no.mhc`` would have reported the same row, but the
  validator ran first, so the build died before reaching it.
* ``species`` and ``mhc.class`` also have no missing guard: ``NaN.lower()`` raises.
* ``reference.id`` matches ``PMID:`` and ``doi:`` case-sensitively while matching ``unpublished``
  case-insensitively.
* ``alleles_match_check`` compares ``int(allele)`` against a per-gene allele count, so ``TRBV20-1*07``
  passes and ``*08`` fails on a count of 7 whether or not IMGT lists ``*07``. It is a range check
  wearing the clothes of a membership check.
* ``gene_match_check`` and ``alleles_match_check`` read
  ``patches/IGM_nomenclature_table.tsv``, which is human only. The driver therefore ORs both masks
  with ``species != 'HomoSapiens'``, and no non-human call was ever checked.
* ``is_qq_seq_biologically_valid`` hardcodes "starts with C, ends with W or F". Measured, that calls
  125 correct rows broken, because ``TRAJ35*01`` templates ``IGFGNVLHC`` and mouse ``TRAJ47*01``
  templates ``HYANKMIC``, so a junction on either legitimately ends in ``C``.
"""
from __future__ import annotations

import csv
import re
from functools import lru_cache
from pathlib import Path

#: ``py_src/ChunkQC.py`` COMPLEX_COLUMNS + METHOD_COLUMNS + META_COLUMNS, in order.
COMPLEX_COLUMNS = ("cdr3.alpha", "v.alpha", "j.alpha", "cdr3.beta", "v.beta", "d.beta", "j.beta",
                   "species", "mhc.a", "mhc.b", "mhc.class", "antigen.epitope", "antigen.gene",
                   "antigen.species", "reference.id")
METHOD_COLUMNS = ("method.identification", "method.frequency", "method.singlecell",
                  "method.sequencing", "method.verification")
META_COLUMNS = ("meta.study.id", "meta.cell.subset", "meta.subject.cohort", "meta.subject.id",
                "meta.replica.id", "meta.clone.id", "meta.epitope.id", "meta.tissue",
                "meta.donor.MHC", "meta.donor.MHC.method", "meta.structure.id")
ALL_COLS = COMPLEX_COLUMNS + METHOD_COLUMNS + META_COLUMNS

#: ``SIGNATURE_COLS``: what legacy deduplicated and reported duplicates on.
SIGNATURE_COLS = (*COMPLEX_COLUMNS, "meta.study.id", "meta.cell.subset", "meta.subject.cohort",
                  "meta.subject.id", "meta.replica.id", "meta.clone.id", "meta.tissue")

#: Lower-cased, as legacy compared it.
SPECIES_LIST = ("homosapiens", "musmusculus", "rattusnorvegicus", "macacamulatta")


def is_aa_seq_valid(aa_seq: str) -> bool:
    if aa_seq == "":                                       # was pd.isnull
        return True
    return len(aa_seq) > 3 and bool(re.match(r"^[ARNDCQEGHILKMFPSTWYV]+$", aa_seq))


def is_mhc_valid(hla_allele: str) -> bool:
    """No missing guard: would have crashed on an absent value. See the module docstring."""
    return (bool(re.match(r"^HLA-[A-Z]+[0-9]?\*\d{2}(:\d{2,3}){0,3}$", hla_allele))
            or hla_allele[0:3] != "HLA")


def _prefix(prefix: str):
    return lambda x: x.startswith(prefix) if x != "" else True


#: ``validators``, key for key and in the legacy order.
VALIDATORS = {
    "cdr3.alpha": is_aa_seq_valid,
    "v.alpha": _prefix("TRAV"),
    "j.alpha": _prefix("TRAJ"),
    "cdr3.beta": is_aa_seq_valid,
    "v.beta": _prefix("TRBV"),
    "d.beta": _prefix("TRBD"),
    "j.beta": _prefix("TRBJ"),
    "species": lambda x: x.lower() in SPECIES_LIST,        # would have crashed on absent
    "mhc.a": is_mhc_valid,                                 # would have crashed on absent
    "mhc.b": is_mhc_valid,                                 # would have crashed on absent
    "mhc.class": lambda x: x == "MHCI" or x == "MHCII",    # would have crashed on absent
    "antigen.epitope": is_aa_seq_valid,
    "antigen.gene": lambda x: x != "",                     # was pd.notnull
    "reference.id": lambda x: (x.startswith("PMID:") or x.startswith("doi:")
                              or x.startswith("http://") or x.startswith("https://")
                              or "unpublished" in x.lower()) if x != "" else True,
}


def header_error(header: list[str], n_rows: int) -> str | None:
    """``check_exist`` then ``check_header``, in the order ``process_chunk`` called them."""
    if not n_rows:
        return "Empty file"
    if len(header) != len(set(header)):
        return f"Duplicate columns found: {header}"
    missing = set(ALL_COLS).difference(header)
    if missing:
        return f"The following columns are missing: {missing}"
    return None


def check_chunk(rows: list[dict[str, str]]) -> set[tuple[int, str]]:
    """``process_chunk``, as ``{(1-based row, finding)}``.

    ``chunk.row`` is 1-based to line up with :func:`vdjdb.qc.rules.check`; legacy keyed on the pandas
    index, which was 0-based.
    """
    out: set[tuple[int, str]] = set()
    seen: set[tuple[str, ...]] = set()
    for i, row in enumerate(rows, start=1):
        signature = tuple(row.get(c, "") for c in SIGNATURE_COLS)
        if signature in seen:
            out.add((i, "duplicate"))
        seen.add(signature)
        for column, valid in VALIDATORS.items():
            if not valid(row.get(column, "")):
                out.add((i, f"bad {column}"))
        if row.get("cdr3.alpha", "") == "" and row.get("cdr3.beta", "") == "":
            out.add((i, "no.cdr3"))
        if row.get("antigen.epitope", "") == "":
            out.add((i, "no.antigen.seq"))
        if row.get("mhc.a", "") == "" or row.get("mhc.b", "") == "":
            out.add((i, "no.mhc"))
    return out


# --------------------------------------------------------------------------------------------
# The master-table tier: runBuidDatabase.py lines 126-152
# --------------------------------------------------------------------------------------------

@lru_cache(maxsize=2)
def _nomenclature(root: Path) -> tuple[frozenset[str], dict[str, int]]:
    """``IMGT/GENE-DB`` names and their allele counts, from the legacy patch table."""
    text = (root / "patches" / "IGM_nomenclature_table.tsv").read_text()
    # The file opens with a blank line. pandas skipped it (``skip_blank_lines`` defaults on), which
    # is the only reason the legacy read found a header at all.
    table = list(csv.DictReader(text.lstrip("\n").splitlines(), delimiter="\t"))
    return (frozenset(r["IMGT/GENE-DB"] for r in table),
            {r["IMGT/GENE-DB"]: int(r["Number of alleles"]) for r in table})


def gene_match_check(gene_name: str, root: Path) -> bool:
    if gene_name == "":
        return True
    names, _ = _nomenclature(root)
    return gene_name.split("*")[0] in names


def alleles_match_check(gene_name: str, root: Path) -> bool:
    if gene_name == "":
        return True
    parts = gene_name.split("*")
    if len(parts) <= 1:
        return True
    names, counts = _nomenclature(root)
    if parts[0] not in names:
        return True
    return int(parts[1]) <= counts[parts[0]]


def is_qq_seq_biologically_valid(aa_seq: str) -> bool:
    """Legacy's guard was ``not isinstance(aa_seq, str)``, which an absent value satisfied."""
    if aa_seq == "":
        return True
    return aa_seq.startswith("C") and (aa_seq.endswith("W") or aa_seq.endswith("F"))


#: The three master-table findings, and the side file each one filled.
MASTER_FINDINGS = {"gene not in IMGT": "vdjdb_full_gene_broken.txt",
                   "allele out of range": "vdjdb_full_allele_broken.txt",
                   "cdr3 not C..[WF]": "vdjdb_full_cdr3aa_broken.txt"}


def check_master_row(row: dict[str, str], root: Path) -> set[str]:
    """The three checks on one repaired master-table row.

    The two nomenclature masks are ORed with ``species != 'HomoSapiens'`` by the driver, so they
    only ever fired on human rows. The junction-anchor mask was not, and fired on every species.
    """
    out: set[str] = set()
    calls = [row.get(f"{s}.{g}", "") for g in ("alpha", "beta") for s in ("v", "j")]
    if row.get("species", "") == "HomoSapiens":
        if not all(gene_match_check(c, root) for c in calls):
            out.add("gene not in IMGT")
        if not all(alleles_match_check(c, root) for c in calls):
            out.add("allele out of range")
    if not all(is_qq_seq_biologically_valid(row.get(f"cdr3.{g}", ""))
               for g in ("alpha", "beta")):
        out.add("cdr3 not C..[WF]")
    return out
