"""IMGT segment nomenclature (#389), driven by the authority table rather than a hand-written map.

``proofreading/imgt_alleles.tsv.gz`` is the authority -- 5,342 alleles across the species VDJdb
carries -- and it has been sitting in the repository unread by any build code. This module is what
makes it live.

**Nothing is invented and nothing is guessed between two candidates.** A call that IMGT already
knows is left alone. Otherwise a small set of *mechanical* respellings is generated and the call is
rewritten **only if exactly one of them is an IMGT name for that species**:

=========================  ==========================================  ====================
respelling                 example                                     chains
=========================  ==========================================  ====================
drop spaces                ``TRBJ 2-7`` -> ``TRBJ2-7``                  2
``.`` -> ``-``             ``TRBJ1.2`` -> ``TRBJ1-2``                   21
``TCR`` -> ``TR``          ``TCRBD2*02`` -> ``TRBD2*02``                5
strip zero padding         ``TRAJ04-1`` -> ``TRAJ4-1``                  —
``-DV`` -> ``/DV``         ``TRAV21-DV12`` -> ``TRAV21/DV12``           181
insert the missing slash   ``TRAV29DV5`` -> ``TRAV29/DV5``              10
drop a D gene's ``-1``     ``TRBD2-1*01`` -> ``TRBD2*01``               224
restore the ``/DV`` name   ``TRAV14`` -> ``TRAV14/DV4``                 1,377
=========================  ==========================================  ====================

Multi-calls are split on ``,`` and ``or``, each member normalised, then rejoined **sorted**, so
``TRBD2,TRBD1`` and ``TRBD1,TRBD2`` stop being two spellings of one fact (89 chains) and
``TRBD2*01 or TCRBD2*02`` becomes machine-readable.

Measured on the corpus: **1,928 of 2,478** non-IMGT calls are resolved this way. The 550 left are
not spelling problems and are reported rather than forced:

* **438 macaque V calls.** 1,333 of 1,771 macaque V calls are valid rhesus IMGT, 206 are valid
  *human* IMGT names applied to macaque records, and 232 are neither. Rewriting a human gene name to
  a rhesus one would be asserting an orthology this build has no basis for.
* **`TRBV8` (20 human), `TRBV7` (17 mouse), `TRAV15` (9 mouse)** and similar: gene names with several
  IMGT candidates and no way to choose. Two candidates is a curation question, not a substitution.
* **`TRAJ16.5`, `TRAJ01-1*01`**: names with no IMGT counterpart at all.

``proofreading/arden.tsv`` supplies the 122 genuinely historical Arden-era names, which no mechanical
rule could derive. None of them appear in the current corpus; it is read so that a chunk carrying one
is normalised on arrival rather than on someone noticing.
"""
from __future__ import annotations

import gzip
import re
from collections.abc import Callable
from functools import lru_cache
from pathlib import Path

import polars as pl

from ..config import Paths

#: VDJdb's species vocabulary -> IMGT's. A species absent here is left untouched: there is no
#: authority to check it against, and an unchecked rewrite is worse than an odd spelling.
IMGT_SPECIES: dict[str, str] = {
    "HomoSapiens": "Homo sapiens",
    "MusMusculus": "Mus musculus",
    "MacacaMulatta": "Macaca mulatta",
    "RattusNorvegicus": "Rattus norvegicus",
}

#: The wide chunk columns holding a segment call.
SEGMENT_COLUMNS: tuple[str, ...] = ("v.alpha", "j.alpha", "v.beta", "d.beta", "j.beta")

#: A curator recording two possible segments writes either a comma or the word "or".
_SPLIT = re.compile(r"\s*(?:,|\bor\b)\s*")


@lru_cache(maxsize=8)
def _imgt(root: Path) -> dict[str, tuple[frozenset[str], frozenset[str]]]:
    """``species -> (allele names, gene names)``, read once per process."""
    with gzip.open(root / "proofreading" / "imgt_alleles.tsv.gz") as fh:
        table = pl.read_csv(fh.read(), separator="\t", infer_schema=False)
    out = {}
    for vdjdb, imgt in IMGT_SPECIES.items():
        rows = table.filter(pl.col("species") == imgt)
        out[vdjdb] = (frozenset(rows["imgt_allele_id"]), frozenset(rows["imgt_gene_id"]))
    return out


@lru_cache(maxsize=8)
def _arden(root: Path) -> dict[str, str]:
    """Arden-era gene names that no mechanical rule derives. Ambiguous entries are dropped."""
    path = root / "proofreading" / "arden.tsv"
    if not path.exists():
        return {}
    table = pl.read_csv(path, separator="\t", infer_schema=False, comment_prefix="#")
    return {a: i for a, i in zip(table["arden_name"], table["imgt_name"], strict=True)
            if i and "," not in i}


def _respellings(call: str, genes: frozenset[str]) -> set[str]:
    """Every mechanical variant of ``call``. Purely syntactic; the caller decides which is real."""
    out = {call.replace(" ", "")}
    out |= {x.replace(".", "-") for x in out}
    out |= {x.replace("TCR", "TR") for x in out}
    out |= {re.sub(r"(?<=[A-Z])0+(\d)", r"\1", x) for x in out}
    out |= {x.replace("-DV", "/DV") for x in out}
    out |= {re.sub(r"(?<=\d)(DV\d)", r"/\1", x) for x in out}
    out |= {re.sub(r"^(TR[AB]D\d)-1", r"\1", x) for x in out}
    # `TRAV14` is IMGT's `TRAV14/DV4`: the gene is shared with the delta locus and IMGT names it
    # once, for both. 1,377 chains write the short form.
    for stem in {x.split("*")[0] for x in out}:
        out |= {g for g in genes if g.startswith(stem + "/DV")}
    return out


def normalise_call(call: str, species: str, root: Path | None = None) -> str | None:
    """The IMGT spelling of ``call``, or ``None`` if it is already IMGT or cannot be resolved.

    ``None`` means "do not rewrite", which covers both "already correct" and "not decidable" -- the
    caller distinguishes them by checking membership itself, and the report does.
    """
    root = root or Paths.discover().root
    tables = _imgt(root)
    if species not in tables or not call:
        return None
    alleles, genes = tables[species]
    known = alleles | genes
    if call in known:
        return None

    arden = _arden(root)
    fixed = []
    for part in (p for p in _SPLIT.split(call.strip()) if p):
        if part in known:
            fixed.append(part)
            continue
        if (a := arden.get(part)) and a in known:
            fixed.append(a)
            continue
        hits = sorted(_respellings(part, genes) & known)
        if len(hits) != 1:          # zero is unknown, several is ambiguous; neither is ours to fix
            return None
        fixed.append(hits[0])
    out = ",".join(sorted(set(fixed)))
    return out if out != call else None


def harmonise_segments(df: pl.DataFrame,
                       root: Path | None = None) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Rewrite the segment columns to IMGT spelling. Returns ``(frame, report)``.

    The report is one row per ``(species, column, from, to)`` with its count -- the audit trail, and
    the source of the declared counts in the difference ledger.
    """
    root = root or Paths.discover().root
    rows: list[dict[str, object]] = []
    for column in SEGMENT_COLUMNS:
        if column not in df.columns:
            continue
        pairs = df.select("species", column).unique().sort("species", column)
        mapping: dict[tuple[str, str], str] = {}
        for species, call in pairs.iter_rows():
            fixed = normalise_call(call, species, root)
            if fixed is not None:
                mapping[(species, call)] = fixed
        if not mapping:
            continue
        # One join, not a per-row call: the distinct (species, call) set is a few hundred pairs.
        lookup = pl.DataFrame(
            {"species": [k[0] for k in mapping], "__from": [k[1] for k in mapping],
             "__to": list(mapping.values())}).sort("species", "__from")
        by_pair = df.join(lookup, left_on=["species", column], right_on=["species", "__from"],
                          how="inner").group_by("species", column, "__to").len().sort(
                              "species", column)
        rows.extend({"column": column, "species": s, "from": f, "to": t, "rows": n}
                    for s, f, t, n in by_pair.iter_rows())
        df = (df.join(lookup, left_on=["species", column], right_on=["species", "__from"],
                      how="left")
              .with_columns(pl.coalesce("__to", column).alias(column))
              .drop("__to"))
    report = (pl.DataFrame(rows, schema={"column": pl.String, "species": pl.String,
                                         "from": pl.String, "to": pl.String, "rows": pl.UInt32})
              .sort("column", "species", "from"))
    return df, report


#: Which legacy columns a wide chunk column becomes. A gene name rewritten in ``v.alpha`` shows up in
#: ``vdjdb_full.txt`` under that name and in ``vdjdb.txt`` / ``vdjdb.slim.txt`` as ``v.segm``.
LEGACY_COLUMNS: dict[str, tuple[str, ...]] = {
    "v.alpha": ("v.alpha", "v.segm"), "j.alpha": ("j.alpha", "j.segm"),
    "v.beta": ("v.beta", "v.segm"), "j.beta": ("j.beta", "j.segm"),
    "d.beta": ("d.beta",),
}

_BEGIN = "# BEGIN generated renames -- `vdjdb rules` writes this block; review it, do not edit it"
_END = "# END generated renames"


def render_renames(report: pl.DataFrame,
                   resolve: Callable[[str, str], str] | None = None) -> str:
    """The ledger's ``[[rename]]`` block, from :func:`harmonise_segments`'s report.

    Generated rather than hand-written, and that is the point: the reviewable artifact is this
    block's diff in a curation pull request, which lists every name the build rewrote. The ledger's
    own job is the complementary one -- proving nothing *else* moved.

    ``resolve(species, call)`` must be the **same** resolution the pipeline applies afterwards. It is
    not optional bookkeeping: the CDR3 fixer maps a segment name onto its germline table and writes
    the resolved name back, so harmonising ``TRAV14`` to ``TRAV14/DV4`` does not put ``TRAV14/DV4``
    in the shipped file -- it puts ``TRAV14/DV4*01``, because the name now *resolves* where the short
    form never did. A rename declaring the intermediate value rewrites the reference and rescues no
    row at all, which is worse than declaring nothing. Pairs that resolve to the same name on both
    sides are dropped: nothing about them reaches the file.
    """
    pairs: dict[tuple[str, str], set[str]] = {}
    rows: dict[tuple[str, str], int] = {}
    for column, species, old, new, n in report.iter_rows():
        if resolve is not None:
            resolved_old, resolved_new = resolve(species, old), resolve(species, new)
            # Only when the fixer leaves the old spelling alone is the rename **injective**. When it
            # does not, it has mapped the odd name onto a real allele -- `get_closest_id` simplifies
            # `TRAV6-7-DV9` to `TRAV6` and then tries `TRAV6-1*01`, `TRAV6-2*01`, ... and takes the
            # first hit -- so the reference ships `TRAV6-1*01` and is indistinguishable from the
            # records that genuinely are TRAV6-1. Declaring that rename rewrites both, and the 15
            # real ones become phantom unmatched rows. Those cases are a declared row delta instead.
            # A multi-call never reaches the file as written: `fix_both` splits it, repairs against
            # each member and keeps the best, so what ships is one of them. That is a selection, not
            # a rename, and declaring it produces a rule that matches nothing.
            if "," in old or "," in new:
                continue
            if resolved_old != old or resolved_new == old:
                continue
            old, new = old, resolved_new
        key = (old, new)
        pairs.setdefault(key, set()).update(LEGACY_COLUMNS.get(column, (column,)))
        rows[key] = rows.get(key, 0) + int(n)
    import json

    # The array-of-tables form, not `rename = [ ... ]`: a bare key after a `[[row_delta]]` header
    # belongs to *that* table, so an inline array appended to the end of the file silently becomes
    # a field of the last section. This form is position-independent.
    out = [_BEGIN]
    for (old, new) in sorted(pairs):
        cols = ",".join(sorted(pairs[(old, new)]))
        out += ["[[rename]]",
                f"columns = {json.dumps(cols)}",
                f"from = {json.dumps(old)}",
                f"to = {json.dumps(new)}",
                f"records = {rows[(old, new)]}",
                ""]
    out.append(_END)
    return "\n".join(out) + "\n"


def legacy_resolver(root: Path | None = None) -> Callable[[str, str], str]:
    """``(species, call) -> the name the CDR3 fixer will write back``.

    Couples the declared renames to the engine that actually ships (`assemble.master.fix_cdr3`
    defaults to ``legacy``). When that default changes, this must change with it.
    """
    from ..annotate._legacy_fixer import Cdr3Fixer

    res = (root or Paths.discover().root) / "res"
    fixer = Cdr3Fixer(str(res / "segments.txt"), str(res / "segments.aaparts.txt"))

    @lru_cache(maxsize=4096)
    def resolve(species: str, call: str) -> str:
        if not call or "," in call:
            return call
        return fixer.get_closest_id(species, call)

    return resolve


def write_renames(report: pl.DataFrame, path: Path,
                  resolve: Callable[[str, str], str] | None = None) -> int:
    """Replace the generated block in ``path`` (appending it if absent). Returns the rename count."""
    block = render_renames(report, resolve)
    text = path.read_text() if path.exists() else ""
    if _BEGIN in text and _END in text:
        head, rest = text.split(_BEGIN, 1)
        _, tail = rest.split(_END, 1)
        text = head + block + tail.lstrip("\n")
    else:
        text = text.rstrip("\n") + "\n\n" + block
    path.write_text(text)
    return block.count("[[rename]]")
