"""Row-level chunk validation, vectorised.

Ported from ``py_src/ChunkQC.py``, which applied each validator with ``.apply`` per cell and built
two of its three emptiness masks with ``chunk_df.T.apply`` -- transpose, then row-wise Python. Here
every rule is one polars expression over the full frame.

Two behaviours are carried over deliberately:

* **``is_MHC_valid`` passes anything not starting with ``HLA``.** The regex only constrains HLA
  spellings; murine ``H2-Kb`` and friends are accepted unchecked. That is why the murine MHC-II
  fragmentation (``I-Ab`` 3,274 vs ``H2-IAb`` 113 vs ``H2-Ab1`` 9) never tripped QC. Phase 9's
  nomenclature rules are where that gets fixed; changing it here would fail 230 chunks at once.
* **A failing chunk must exit non-zero.** The Groovy build did; the Python port replaced it with
  ``warnings.warn`` and carried on. Measured, 230 of 230 chunks pass today, so restoring the
  non-zero exit needs no quarantine list.
"""
from __future__ import annotations

import polars as pl

from ..schema import CHUNK_DEDUP_KEY, SPECIES

#: The 20 proteinogenic amino acids. Measured: all 305,031 non-empty CDR3 cells in ``chunks/`` are
#: already clean, because unconventional-residue records are quarantined outside ``chunks/``.
AA = "ARNDCQEGHILKMFPSTWYV"

_AA_SEQ = rf"^[{AA}]{{4,}}$"
#: ``HLA-<gene><digit?>*NN(:NN){0,3}``. Anything not starting with ``HLA`` is accepted -- see above.
_HLA = r"^HLA-[A-Z]+[0-9]?\*\d{2}(:\d{2,3}){0,3}$"
_REFERENCE = r"(?i)^(PMID:|doi:|https?://)|unpublished"


def _blank(col: str) -> pl.Expr:
    return pl.col(col) == ""


def _seq_ok(col: str) -> pl.Expr:
    """Empty is allowed; a present sequence must be > 3 residues of the standard 20."""
    return _blank(col) | pl.col(col).str.contains(_AA_SEQ)


def _prefix_ok(col: str, prefix: str) -> pl.Expr:
    return _blank(col) | pl.col(col).str.starts_with(prefix)


#: ``rule id -> expression that is True when the row is GOOD``. Expressed positively so a rule
#: reads as the invariant it protects rather than as the failure it catches.
RULES: dict[str, pl.Expr] = {
    "bad cdr3.alpha": _seq_ok("cdr3.alpha"),
    "bad cdr3.beta": _seq_ok("cdr3.beta"),
    "bad antigen.epitope": _seq_ok("antigen.epitope"),
    "bad v.alpha": _prefix_ok("v.alpha", "TRAV"),
    "bad j.alpha": _prefix_ok("j.alpha", "TRAJ"),
    "bad v.beta": _prefix_ok("v.beta", "TRBV"),
    "bad d.beta": _prefix_ok("d.beta", "TRBD"),
    "bad j.beta": _prefix_ok("j.beta", "TRBJ"),
    "bad species": pl.col("species").is_in(list(SPECIES)),
    "bad mhc.a": _blank("mhc.a") | ~pl.col("mhc.a").str.starts_with("HLA")
                 | pl.col("mhc.a").str.contains(_HLA),
    "bad mhc.b": _blank("mhc.b") | ~pl.col("mhc.b").str.starts_with("HLA")
                 | pl.col("mhc.b").str.contains(_HLA),
    "bad mhc.class": pl.col("mhc.class").is_in(["MHCI", "MHCII"]),
    "bad antigen.gene": ~_blank("antigen.gene"),
    "bad reference.id": _blank("reference.id") | pl.col("reference.id").str.contains(_REFERENCE),
    "no.cdr3": ~(_blank("cdr3.alpha") & _blank("cdr3.beta")),
    "no.antigen.seq": ~_blank("antigen.epitope"),
    "no.mhc": ~(_blank("mhc.a") | _blank("mhc.b")),
    # #561. A paired record whose two chains have the same CDR3 is a transcription error: the
    # beta sequence copied into the alpha field, with the V and J calls left correct. Which chain is
    # wrong cannot be known from the row, so this reports and does not repair -- 99 records on the
    # current corpus, 98 of them from two references. Advisory, so the build does not fail on a
    # defect only a curator can fix.
    "alpha and beta cdr3 identical": (_blank("cdr3.alpha") | _blank("cdr3.beta")
                                      | (pl.col("cdr3.alpha") != pl.col("cdr3.beta"))),
    # A V or J call for a chain whose CDR3 is blank. The call is real information and the row is
    # kept, but the chain cannot reach any output: every shipped table is keyed on the CDR3, so the
    # segment is carried in `chunks/` and dropped by the build, which is the kind of silent loss a
    # submitter should hear about while they can still fix it.
    #
    # Two ways a row gets here. A submission that names the genes and omits the sequence, which is
    # the case this flags for the submitter. And #561, where the alpha CDR3 was a copy of the beta
    # and was cleared, leaving the paper's genuine alpha V and J behind on purpose -- that is a
    # deliberate state, recorded in the row's `comment`, and the flag keeps it visible rather than
    # letting it look like a clean record.
    #
    # Advisory: it reports a chain that cannot ship, not a row that is wrong.
    "segment call with no cdr3": ~(
        (_blank("cdr3.alpha") & ~(_blank("v.alpha") & _blank("j.alpha")))
        | (_blank("cdr3.beta") & ~(_blank("v.beta") & _blank("j.beta")))
    ),

    # `meta.structure.id` is specified as "PDB structure ID if one exists, blank otherwise"
    # (docs/standards/chunk-format.md), and `score.confidence` reads it as the strongest evidence
    # there is: a non-empty value awards 3 outright, above every sequencing and specificity term,
    # because a solved TCR:pMHC complex is direct proof of binding.
    #
    # So the field is not free text, and anything in it that is not a PDB id awards the top score
    # for evidence that does not exist. Measured 2026-09-28: 2,765 of the 3,246 rows carrying a
    # value hold a figure or table reference instead - `Fig.2, Fig. 3, Fig.4, ...` 2,352 rows,
    # `Fig 9, Supp Fig 5, Supp Table 5-8` 400, `Fig3b,Fig3c` 12, `56I` 1. Because the score is a
    # maximum over the sample signature, they reach much further than their own rows: **4,817 of
    # the 6,007 records that score 3 owe it to one of those values**, and blanking the field would
    # take 3,752 of them to 0 and 1,065 to 1.
    #
    # The shape is the PDB entry id: a digit then three alphanumerics, case-insensitive. That is the
    # whole check; whether the entry exists and contains the receptor is a network question and
    # belongs in the structure pass of #402, not in a chunk gate that must run offline.
    #
    # Advisory: the corpus fails it on 2,765 rows today, and what to do with those is a curation
    # decision rather than a submission error.
    "structure id is not a PDB id": (
        _blank("meta.structure.id")
        | pl.col("meta.structure.id").str.contains(r"^[0-9][A-Za-z0-9]{3}$")
    ),
}


def check(df: pl.DataFrame) -> pl.DataFrame:
    """One row per (chunk file, chunk row, failing rule).

    ``duplicate`` marks the second and later rows sharing a :data:`CHUNK_DEDUP_KEY` within one
    chunk. Those are removed by the reader's per-chunk deduplication, so they are reported as a
    curation signal rather than an error that blocks the build.
    """
    flagged = df.with_columns(
        *(expr.not_().alias(f"__{rule}") for rule, expr in RULES.items()),
        (pl.int_range(pl.len()).over(["chunk.file", *CHUNK_DEDUP_KEY]) > 0).alias("__duplicate"),
    )
    cols = [f"__{r}" for r in RULES] + ["__duplicate"]
    return (
        flagged.select("chunk.file", "chunk.row", *cols)
        .unpivot(index=["chunk.file", "chunk.row"], on=cols,
                 variable_name="rule", value_name="failed")
        .filter("failed")
        .with_columns(pl.col("rule").str.strip_prefix("__"))
        .drop("failed")
        .sort("chunk.file", "chunk.row", "rule")
    )


def summarise(findings: pl.DataFrame) -> pl.DataFrame:
    """Findings per rule, most frequent first."""
    return (findings.group_by("rule")
            .agg(pl.len().alias("rows"), pl.col("chunk.file").n_unique().alias("chunks"))
            .sort("rows", descending=True))
