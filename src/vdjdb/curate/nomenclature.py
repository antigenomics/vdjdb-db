"""IMGT segment nomenclature (#389), driven by the authority table rather than a hand-written map.

``proofreading/imgt_alleles.tsv.gz`` is the authority: 5,342 alleles across the species VDJdb
covers. This module reads it.

Nothing is invented and nothing is guessed between two candidates. A call that IMGT already knows
is left alone. Otherwise a small set of mechanical respellings is generated, and the call is
rewritten only if exactly one of them is an IMGT name for that species:

=========================  ==========================================  ====================
respelling                 example                                     chains
=========================  ==========================================  ====================
drop spaces                ``TRBJ 2-7`` -> ``TRBJ2-7``                  2
``.`` -> ``-``             ``TRBJ1.2`` -> ``TRBJ1-2``                   21
``TCR`` -> ``TR``          ``TCRBD2*02`` -> ``TRBD2*02``                5
strip zero padding         ``TRAJ04-1`` -> ``TRAJ4-1``                  -
``-DV`` -> ``/DV``         ``TRAV21-DV12`` -> ``TRAV21/DV12``           181
insert the missing slash   ``TRAV29DV5`` -> ``TRAV29/DV5``              10
drop a D gene's ``-1``     ``TRBD2-1*01`` -> ``TRBD2*01``               224
restore the ``/DV`` name   ``TRAV14`` -> ``TRAV14/DV4``                 1,377
=========================  ==========================================  ====================

Multi-calls are split on ``,`` and ``or``, each member normalised, then rejoined in sorted order,
so ``TRBD2,TRBD1`` and ``TRBD1,TRBD2`` stop being two spellings of one fact (89 chains) and
``TRBD2*01 or TCRBD2*02`` becomes machine-readable.

Measured on the corpus: 1,928 of 2,478 non-IMGT calls are resolved this way. The 550 left are not
spelling problems and are reported rather than forced:

* **438 macaque V calls.** 1,333 of 1,771 macaque V calls are valid rhesus IMGT, 206 are valid
  *human* IMGT names applied to macaque records, and 232 are neither. Rewriting a human gene name to
  a rhesus one would be asserting an orthology this build has no basis for.
* **`TRBV8` (20 human), `TRBV7` (17 mouse), `TRAV15` (9 mouse)** and similar: gene names with several
  IMGT candidates and no way to choose. Two candidates is a curation question, not a substitution.
* **`TRAJ16.5`, `TRAJ01-1*01`**: names with no IMGT counterpart at all.

``proofreading/arden.tsv`` supplies the 122 historical Arden-era names, which no mechanical rule
derives. None appear in the current corpus; the table is read so that a chunk using one is
normalised on arrival.
"""
from __future__ import annotations

import gzip
import json
import re
from collections.abc import Callable
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import polars as pl

from ..config import Paths

#: VDJdb's species vocabulary -> IMGT's. A species absent here is left untouched: there is no
#: authority table to check its calls against.
IMGT_SPECIES: dict[str, str] = {
    "HomoSapiens": "Homo sapiens",
    "MusMusculus": "Mus musculus",
    "MacacaMulatta": "Macaca mulatta",
    "RattusNorvegicus": "Rattus norvegicus",
}

#: The wide chunk columns holding a segment call.
SEGMENT_COLUMNS: tuple[str, ...] = ("v.alpha", "j.alpha", "v.beta", "d.beta", "j.beta")

#: A curator recording two possible segments writes a comma, a semicolon, a plus or the word "or".
#: The three punctuation marks are interchangeable here because no IMGT TR name contains any of them,
#: and a call is rewritten only when *every* part resolves to a real name, so widening the set cannot
#: produce a substitution - it either resolves the whole call or leaves it alone. Measured: 591
#: chain-calls over 9 spellings, the largest mouse ``TRBV12-2+TRBV13-2`` (535) and
#: ``TRBV3-1;TRBV3-2`` (30). The retired build rejected every one of them as `gene not in IMGT`.
_SPLIT = re.compile(r"\s*(?:[,;+]|\bor\b)\s*")


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


#: `TRBVIS1` in `PMID_16237109.txt` is `TRBV1S1` with a capital I for the 1 (#136). Roman numerals
#: only ever appear in the Arden `TRBV<n>S<m>` form, so the rewrite is scoped to that shape and
#: cannot touch an IMGT name.
_ROMAN = {"I": "1", "II": "2", "III": "3", "IV": "4", "V": "5",
          "VI": "6", "VII": "7", "VIII": "8", "IX": "9", "X": "10"}


def _respellings(call: str, genes: frozenset[str]) -> set[str]:
    """Every mechanical variant of ``call``. Purely syntactic; the caller decides which to use."""
    out = {call.replace(" ", "")}
    out |= {x.replace(".", "-") for x in out}
    out |= {x.replace("TCR", "TR") for x in out}
    out |= {re.sub(r"(?<=[A-Z])0+(\d)", r"\1", x) for x in out}
    # A one-digit allele: IMGT writes `*01`, and 26 chains write `*1`, `*2` or `*3` (#402). Only
    # accepted if the padded form is a real allele, which is what the caller checks.
    out |= {re.sub(r"\*(\d)$", r"*0\1", x) for x in out}
    # The Arden roman numeral, `TRBVIS1` -> `TRBV1S1`, which the Arden table then maps (#136).
    out |= {m.group(1) + _ROMAN[m.group(2)] + m.group(3)
            for m in (re.match(r"^(TR[ABGD][VDJ])([IVX]+)(S\d+)$", x) for x in out)
            if m and m.group(2) in _ROMAN}
    # An en or em dash standing in for a hyphen: one `j.beta` reads `TRBJ2<en dash>3*01`, copied out of a
    # typeset table. Neither character appears in any IMGT name, so this cannot touch a correct call,
    # and a literal one is indistinguishable from a hyphen in a diff - which is how it arrived.
    out |= {x.replace("\u2013", "-").replace("\u2014", "-") for x in out}
    out |= {x.replace("-DV", "/DV") for x in out}
    out |= {re.sub(r"(?<=\d)(DV\d)", r"/\1", x) for x in out}
    out |= {re.sub(r"^(TR[AB]D\d)-1", r"\1", x) for x in out}
    # A spurious `-1`: IMGT names the gene `TRBV19`, and 552 chains write `TRBV19-1`. Only `-1`,
    # and only where the family has no second member -- the same criterion arda 2.29 applies to the
    # reverse direction (a family call resolves only when it holds one functional gene). Dropping
    # any `-N` instead would rewrite `TRAJ37-2` to `TRAJ37` and `TRBV13-6` to `TRBV13`, discarding
    # a distinction the curator made rather than repairing a spelling.
    out |= {m.group(1) + (m.group(2) or "")
            for m in (re.match(r"^(TR[ABDG][VDJ]\d+)-1(\*\d+)?$", x) for x in out)
            if m and f"{m.group(1)}-2" not in genes}
    # `TRAV14` is IMGT's `TRAV14/DV4`: the gene is shared with the delta locus and IMGT names it
    # once, for both. 1,377 chains write the short form.
    for stem in {x.split("*")[0] for x in out}:
        out |= {g for g in genes if g.startswith(stem + "/DV")}
    return out


def _expand_slash(part: str, known: frozenset[str]) -> list[str]:
    """Read a ``/`` as "or" and return the names it stands for, or ``[]`` if it cannot be read.

    ``/`` means three different things in a segment call and one of them is not a separator at all:
    ``TRAV14/DV4`` is a single IMGT gene shared with the delta locus. That case never reaches here,
    because the caller checks IMGT membership first. What is left is a curator writing alternatives
    (#136), in three shapes, each accepted only when every name it expands to is a real allele:

    * ``TRBV12-3/TRBV12-4`` - full names;
    * ``TRBV19*01/02`` - one gene, an allele list;
    * ``TRBV12-3/4*01`` - one family, a gene-number list sharing an allele.

    ``TRBV11/2`` matches none of them: it could be `TRBV11-2` or two genes, and guessing which is a
    curation decision rather than a spelling one. It is left alone and reported.
    """
    if "/" not in part:
        return []
    whole = part.split("/")
    if all(x in known for x in whole):
        return whole
    if m := re.match(r"^(.+?)\*(\d+(?:/\d+)+)$", part):
        out = [f"{m.group(1)}*{a.zfill(2)}" for a in m.group(2).split("/")]
        if all(x in known for x in out):
            return out
    if m := re.match(r"^(TR[ABGD][VDJ]\d+)-(\d+(?:/\d+)+)(\*\d+)?$", part):
        out = [f"{m.group(1)}-{n}{m.group(3) or ''}" for n in m.group(2).split("/")]
        if all(x in known for x in out):
            return out
    return []


def normalise_call(call: str, species: str, root: Path | None = None) -> str | None:
    """The IMGT spelling of ``call``, or ``None`` if it is already IMGT or cannot be resolved.

    ``None`` means "do not rewrite", which covers both "already correct" and "not decidable"; the
    caller and the report tell the two apart by checking IMGT membership themselves.
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
        if expanded := _expand_slash(part, known):
            fixed.extend(expanded)
            continue
        variants = _respellings(part, genes)
        hits = sorted(variants & known)
        if not hits:
            # A respelling may land on an *Arden* name rather than an IMGT one: `TRBVIS1` respells
            # to `TRBV1S1`, which the Arden table maps to `TRBV9`. Trying Arden only after the
            # direct lookup and the respellings have both missed keeps the order of evidence -- an
            # IMGT name is never reinterpreted as Arden.
            hits = sorted({a for v in variants if (a := arden.get(v)) and a in known})
        if len(hits) != 1:          # zero is unknown, several is ambiguous; neither is ours to fix
            return None
        fixed.append(hits[0])
    out = ",".join(sorted(set(fixed)))
    return out if out != call else None


def harmonise_segments(df: pl.DataFrame,
                       root: Path | None = None) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Rewrite the segment columns to IMGT spelling. Returns ``(frame, report)``.

    The report is one row per ``(species, column, from, to)`` with its count -- the audit trail, and
    the source of the declared counts in ``rules/expected_diffs.toml``.
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


#: Columns :func:`unresolved` returns, in order.
UNRESOLVED_COLUMNS: tuple[str, ...] = ("species", "column", "call", "part", "chains",
                                       "family.members", "candidates")


def unresolved(df: pl.DataFrame, root: Path | None = None) -> pl.DataFrame:
    """Every segment call IMGT has at neither allele nor gene level, after harmonisation.

    This is the report the retired build wrote as ``vdjdb_full_gene_broken.txt`` and
    ``vdjdb_full_allele_broken.txt``: ``runBuidDatabase.py`` applied ``gene_match_check`` and
    ``alleles_match_check`` to the repaired master table and filed the rows that failed. Nothing
    replaced it. :func:`harmonise_segments` reports what it *rewrote*, and ``build_master`` discards
    even that, so every row naming a gene no authority carries reached every shipped table with no
    report anywhere. Measured 2026-09-29 over the built corpus: 88 rows over 69 distinct names
    and 3,457 chain-calls, against the 2,582 human chain-calls the retired check saw.

    Three improvements on the check it replaces, each of which was a legacy defect:

    * **All four species, not only human.** The legacy table was human immunoglobulin-focused, so the
      driver ORed both masks with ``species != 'HomoSapiens'`` and no macaque or mouse call was ever
      checked. 850 of the 3,457 chain-calls here are not human.
    * **Membership, not a range.** ``alleles_match_check`` compared ``int(allele)`` against a per-gene
      allele count, which passes ``*07`` and fails ``*08`` whether or not IMGT lists either. Its one
      finding on the corpus, ``TRBV28*02``, is an allele IMGT does have.
    * **Per part.** A curator recording two candidates writes ``TRBD1,TRBD2``, which is not an allele
      name and which the legacy check therefore rejected whole. Each part is resolved separately here.

    ``family.members`` and ``candidates`` are what separate the two reasons a call lands here.
    A positive count means the name is a *family* whose members IMGT lists individually
    (``TRBV6`` -> nine of them), so the record is under-specified and only a curator can choose. Zero
    means no IMGT name is near it at all, which is a spelling defect or a gene that species does not
    have. Advisory either way: #389 is the issue that works through them.
    """
    root = root or Paths.discover().root
    tables = _imgt(root)
    rows: list[dict[str, object]] = []
    for column in SEGMENT_COLUMNS:
        if column not in df.columns:
            continue
        for species, call, n in df.group_by("species", column).len().sort("species", column).rows():
            if not call or species not in tables:
                continue
            alleles, genes = tables[species]
            for part in (x for x in _SPLIT.split(call.strip()) if x):
                if part in alleles or part.split("*")[0] in genes:
                    continue
                family = sorted(g for g in genes if g.startswith(part.split("*")[0] + "-"))
                rows.append({"species": species, "column": column, "call": call, "part": part,
                             "chains": n, "family.members": len(family),
                             "candidates": ",".join(family[:4])})
    return (pl.DataFrame(rows, schema={"species": pl.String, "column": pl.String, "call": pl.String,
                                       "part": pl.String, "chains": pl.UInt32,
                                       "family.members": pl.UInt32, "candidates": pl.String})
            .sort("chains", "species", "column", "part", descending=[True, False, False, False]))


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
    """The ``[[rename]]`` block of ``rules/expected_diffs.toml``, from :func:`harmonise_segments`.

    Generated rather than hand-written: the reviewable artifact is this block's diff in a curation
    pull request, which lists every name the build rewrote. ``vdjdb diff`` then checks that nothing
    else moved.

    ``resolve(species, call)`` must be the same resolution the pipeline applies afterwards. The CDR3
    fixer maps a segment name onto its germline table and writes the resolved name back, so
    harmonising ``TRAV14`` to ``TRAV14/DV4`` does not put ``TRAV14/DV4`` in the shipped file -- it
    puts ``TRAV14/DV4*01``, because the name resolves where the short form did not. A rename
    declaring the intermediate value rewrites the reference and matches no row. Pairs that resolve
    to the same name on both sides are dropped: nothing about them reaches the file.
    """
    # Keyed on the species too, because a rename is a statement about one organism's nomenclature:
    # `TRAV14-1*01` is one human record's respelling and 79 mouse records' correct IMGT name, and a
    # key that drops the species declares one rewrite for both (#671).
    pairs: dict[tuple[str, str, str], set[str]] = {}
    rows: dict[tuple[str, str, str], int] = {}
    for column, species, old, new, n in report.iter_rows():
        if resolve is not None:
            resolved_old, resolved_new = resolve(species, old), resolve(species, new)
            # Only when the fixer leaves the old spelling alone is the rename injective. When it
            # does not, it has mapped the odd name onto an existing allele -- `get_closest_id`
            # simplifies `TRAV6-7-DV9` to `TRAV6` and then tries `TRAV6-1*01`, `TRAV6-2*01`, ... and
            # takes the first hit -- so the reference ships `TRAV6-1*01` and is indistinguishable
            # from the records that are TRAV6-1. Declaring that rename rewrites both, and the 15
            # TRAV6-1 records become phantom unmatched rows. Those are a declared row delta instead.
            # A multi-call never reaches the file as written: `fix_both` splits it, repairs against
            # each member and keeps the best, so what ships is one of them. That is a selection, not
            # a rename, and declaring it produces a rule that matches nothing.
            if "," in old or "," in new:
                continue
            if resolved_old != old or resolved_new == old:
                continue
            old, new = old, resolved_new
        key = (species, old, new)
        pairs.setdefault(key, set()).update(LEGACY_COLUMNS.get(column, (column,)))
        rows[key] = rows.get(key, 0) + int(n)
    import json

    # The array-of-tables form, not `rename = [ ... ]`: a bare key after a `[[row_delta]]` header
    # belongs to that table, so an inline array appended to the end of the file becomes a field of
    # the last section with no error. This form is position-independent.
    out = [_BEGIN]
    for key in sorted(pairs):
        species, old, new = key
        out += ["[[rename]]",
                f"columns = {json.dumps(','.join(sorted(pairs[key])))}",
                f"from = {json.dumps(old)}",
                f"to = {json.dumps(new)}",
                f"species = {json.dumps(species)}",
                f"records = {rows[key]}",
                ""]
    out.append(_END)
    return "\n".join(out) + "\n"


def allele_resolver() -> Callable[[str, str], str]:
    """``(species, call) -> the name the markup engine writes back``.

    A declared ``[[rename]]`` has to state the value the shipped build produces, and the markup engine
    names the allele it aligned against - a bare ``TRBV12-3`` comes back ``TRBV12-3*01``, so
    harmonising ``TRAV14`` to ``TRAV14/DV4`` does not put ``TRAV14/DV4`` in the file, it puts
    ``TRAV14/DV4*01``. That engine is ``arda.cdr3fix``, and this is its own ladder,
    ``arda.cdr3fix.resolve_allele`` over arda's anchor table for the organism.

    ``call`` unchanged where the ladder resolves nothing, which is what :func:`render_renames` reads as
    "the engine leaves this name alone". A D call always takes that path: arda's anchor table carries V
    and J only, because a D germline is not part of a V·J scaffold, and the legacy ``vdjdb.txt`` has no
    D column for a rename to touch anyway.

    This used to run the retired k-mer scanner's ``get_closest_id`` over ``res/segments.txt``, whose
    docstring said that coupled the renames to "the engine that actually ships". It had not for some
    time: ``fix_cdr3`` has defaulted to ``arda`` since phase 5, so the renames were resolved by one
    ladder and produced by another (#658).

    ⚠ **The committed block in ``rules/expected_diffs.toml`` was generated by the retired ladder, and
    regenerating it here does not reproduce it: the committed block has 353 renames and a
    regeneration produces 444.** arda's ladder is deliberately stricter - it refuses a family with
    several functional genes where VDJdb's took the lowest-numbered one - so 91 names it declines are
    names ``get_closest_id`` mapped, and
    :func:`render_renames` emits a rename for each. Some are right and some would rewrite the wrong
    rows: ``TRBV13-1*01 -> TRBV13*01`` is one human record's respelling, and the release holds
    ``TRBV13-1*01`` on 653 rows because it is also a real *mouse* gene, so an unscoped rename would
    rewrite all 653. Making the renames species-scoped is the fix and it is its own change with its own
    declared counts: **do not regenerate the block on a branch that is not about that.** #671.
    """
    from arda.cdr3fix import load_anchors, resolve_species

    @lru_cache(maxsize=8)
    def anchors(species: str) -> dict:
        # `vdjdb rules` runs no build, so nothing has pointed arda at a germline reference yet, and
        # without this every call resolves to `""` and every rename loses its `to` value. Silently:
        # `load_anchors` returns an empty table rather than raising.
        from ..annotate.cdr3fix import ensure_reference

        ensure_reference()
        try:
            return load_anchors(resolve_species(species))
        except FileNotFoundError:
            return {}

    @lru_cache(maxsize=4096)
    def resolve(species: str, call: str) -> str:
        if not call or "," in call:
            return call
        segment = call[3:4].upper()
        if segment not in ("V", "J"):
            return call
        from arda.cdr3fix import resolve_allele

        return resolve_allele(call, segment, anchors(species)) or call

    return resolve


def render_allele_renames(report: pl.DataFrame,
                          resolve: Callable[[str, str], str] | None = None) -> str:
    """The ``[[rename]]`` lines for :func:`disambiguate_alleles`, with their evidence.

    An allele correction is not injective on value alone -- the 1,047 TRAJ24 records the CDR3
    identifies as ``*02`` ship the same ``TRAJ24*01`` as the 38 it does not, because the fixer
    resolves a bare ``TRAJ24`` to ``*01``. The rename therefore states the same predicate the rule
    used, and ``vdjdb diff`` applies it to the reference under the same evidence.
    """
    seen: dict[tuple[str, str, str, str, str], int] = {}
    cdr3_for = {"j.alpha": "cdr3,cdr3.alpha", "j.beta": "cdr3,cdr3.beta",
                "v.alpha": "cdr3,cdr3.alpha", "v.beta": "cdr3,cdr3.beta"}
    out = []
    for _issue, column, species, old, new, signature, n in report.iter_rows():
        if resolve is not None:
            old, new = resolve(species, old), resolve(species, new)
        if old == new:
            continue
        key = (column, species, old, new, signature)
        seen[key] = seen.get(key, 0) + int(n)
    for (column, species, old, new, signature), n in sorted(seen.items()):
        cols = ",".join(sorted(set(LEGACY_COLUMNS.get(column, (column,)))))
        out += ["[[rename]]",
                f"columns = {json.dumps(cols)}",
                f"from = {json.dumps(old)}",
                f"to = {json.dumps(new)}",
                f"species = {json.dumps(species)}",
                f"when_columns = {json.dumps(cdr3_for.get(column, 'cdr3'))}",
                f"when_contains = {json.dumps(signature)}",
                f"records = {n}",
                ""]
    return "\n".join(out)


def render_mhc_renames(report: pl.DataFrame) -> str:
    """The ``[[rename]]`` lines for :func:`harmonise_mhc`'s value substitutions.

    Unconditional and injective: unlike a segment call, an MHC string is not touched by the CDR3
    fixer, so what the reference ships is what the chunk wrote. The chain-order swap is not here:
    it rewrites two columns at once and is a declared row delta instead.
    """
    out = []
    for _issue, column, old, new, n in report.iter_rows():
        if "," in column:
            continue
        out += ["[[rename]]", f"columns = {json.dumps(column)}", f"from = {json.dumps(old)}",
                f"to = {json.dumps(new)}", f"records = {n}", ""]
    return "\n".join(out)


def write_renames(report: pl.DataFrame, path: Path,
                  resolve: Callable[[str, str], str] | None = None,
                  alleles: pl.DataFrame | None = None,
                  mhc: pl.DataFrame | None = None,
                  extra_blocks: tuple[str, ...] = ()) -> int:
    """Replace the generated block in ``path`` (appending it if absent). Returns the rename count."""
    block = render_renames(report, resolve)
    for extra in (render_allele_renames(alleles, resolve) if alleles is not None
                  and not alleles.is_empty() else "",
                  render_mhc_renames(mhc) if mhc is not None and not mhc.is_empty() else ""):
        if extra:
            block = block.replace(_END, extra + _END)
    for extra in extra_blocks:
        if extra:
            block = block.replace(_END, extra + _END)
    text = path.read_text() if path.exists() else ""
    if _BEGIN in text and _END in text:
        head, rest = text.split(_BEGIN, 1)
        _, tail = rest.split(_END, 1)
        text = head + block + tail.lstrip("\n")
    else:
        text = text.rstrip("\n") + "\n\n" + block
    path.write_text(text)
    return block.count("[[rename]]")


@dataclass(frozen=True, slots=True)
class AlleleSignature:
    """Two alleles of one gene that the CDR3 itself tells apart.

    Where two alleles differ inside the junction, the sequence is evidence and the call is not: a
    record whose CDR3 has the ``*02`` residues is ``*02``, whatever the submitter wrote.
    """

    issue: str
    species: str
    column: str
    cdr3: str
    prefix: str
    #: CDR3 substring -> the allele it proves. Substrings must be mutually exclusive; a CDR3
    #: matching two of them is a contradiction and is left alone.
    signatures: dict[str, str]


#: #327. TRAJ24*01 encodes ``...GGK**FE**F...`` and *02 ``...GGK**LQ**F...``, two residues apart and
#: both inside the junction. Measured on the corpus: ``WGKFEF`` appears zero times and ``WGKLQF``
#: 1,080 times across the TRAJ24 family, including in 73 of the 111 records explicitly called
#: ``*01``. The original report was that about two thirds of explicit ``*01`` calls are probably
#: ``*02``; no record with an explicit ``*01`` call has the ``*01`` signature. The 364 with neither
#: signature have a CDR3 trimmed short of the anchor and are left alone: no evidence, no correction.
#: #647. The mouse TRAJ47 pair is the same shape with the anchor itself as the discriminator:
#: ``*01`` templates ``HYANKMI**C**`` and is an ORF allele, ``*02`` templates ``DYANKMI**F**``.
#: All 95 mouse TRAJ47 chains read ``ANKMIF`` and none reads ``ANKMIC``, so every one of them is
#: ``*02``: 81 were written without an allele and resolved to ``*01`` by default, 14 name ``*01``
#: outright. Human TRAJ47 is untouched, because there both alleles template ``EYGNKLVF`` and the
#: sequence cannot tell them apart. The three human ``TRAJ47-1`` spellings are a name with no IMGT
#: counterpart and belong to #389; the species scope keeps them out of reach here.
ALLELE_SIGNATURES: tuple[AlleleSignature, ...] = (
    AlleleSignature(issue="#327", species="HomoSapiens", column="j.alpha", cdr3="cdr3.alpha",
                    prefix="TRAJ24",
                    signatures={"WGKLQF": "TRAJ24*02", "WGKFEF": "TRAJ24*01"}),
    AlleleSignature(issue="#647", species="MusMusculus", column="j.alpha", cdr3="cdr3.alpha",
                    prefix="TRAJ47",
                    signatures={"ANKMIF": "TRAJ47*02", "ANKMIC": "TRAJ47*01"}),
)


def disambiguate_alleles(df: pl.DataFrame) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Set the allele from the CDR3 where the sequence decides it. Returns ``(frame, report)``."""
    rows: list[dict[str, object]] = []
    for sig in ALLELE_SIGNATURES:
        if sig.column not in df.columns or sig.cdr3 not in df.columns:
            continue
        family = (pl.col("species") == sig.species) & pl.col(sig.column).str.starts_with(sig.prefix)
        # A CDR3 matching two signatures at once contradicts itself; leave it to a curator.
        hits = [pl.col(sig.cdr3).str.contains(s, literal=True) for s in sig.signatures]
        unambiguous = family & (pl.sum_horizontal(*[h.cast(pl.Int8) for h in hits]) == 1)
        expr = pl.col(sig.column)
        for substring, allele in sig.signatures.items():
            mask = unambiguous & pl.col(sig.cdr3).str.contains(substring, literal=True)
            changed = df.filter(mask & (pl.col(sig.column) != allele))
            if not changed.is_empty():
                rows.extend(
                    {"issue": sig.issue, "column": sig.column, "species": sig.species,
                     "from": f, "to": allele, "signature": substring, "rows": n}
                    for f, n in changed.group_by(sig.column).len().sort(sig.column).iter_rows())
            expr = pl.when(mask).then(pl.lit(allele)).otherwise(expr)
        df = df.with_columns(expr.alias(sig.column))
    report = pl.DataFrame(rows, schema={"issue": pl.String, "column": pl.String,
                                        "species": pl.String, "from": pl.String, "to": pl.String,
                                        "signature": pl.String, "rows": pl.UInt32})
    return df, report.sort("issue", "column", "from")


# ---------------------------------------------------------------------------------------------
# MHC
# ---------------------------------------------------------------------------------------------

#: Reference-scoped corrections belong in ``patches/mhc.dict``, not here: they are data a curator
#: reviews as a diff, with the source that justifies each one in the row beside it. The file lists
#: the murine class-II spellings, the alleles that do not exist in IPD-IMGT/HLA, and -- as comments
#: -- the three cases checked and deliberately left alone.
MHC_PATCH = "mhc.dict"

#: Class-II alpha-chain genes. ``mhc.a`` is the first chain, ``mhc.b`` the second one.
MHC_ALPHA: tuple[str, ...] = ("HLA-DRA", "HLA-DQA1", "HLA-DQA2", "HLA-DPA1", "HLA-DPA2")
MHC_BETA: tuple[str, ...] = ("HLA-DRB1", "HLA-DRB3", "HLA-DRB4", "HLA-DRB5", "HLA-DQB1",
                             "HLA-DQB2", "HLA-DPB1", "HLA-DPB2")

MHC_COLUMNS: tuple[str, ...] = ("mhc.a", "mhc.b")


@lru_cache(maxsize=4)
def _mhc_patch(root: Path) -> tuple[tuple[str, str, str], ...]:
    """``(value, replacement, reference.id or "*")`` from ``patches/mhc.dict``."""
    path = root / "patches" / MHC_PATCH
    if not path.exists():
        return ()
    table = pl.read_csv(path, separator="\t", infer_schema=False, comment_prefix="#")
    return tuple(zip(table["mhc"], table["replacement"], table["reference.id"], strict=True))


def harmonise_mhc(df: pl.DataFrame, root: Path | None = None) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Collapse MHC spellings and put the alpha chain in ``mhc.a``. Returns ``(frame, report)``.

    Two corrections, each with its own evidence:

    * **declared corrections** -- ``patches/mhc.dict``: the murine class-II spellings that split one
      molecule several ways, and the alleles IPD-IMGT/HLA does not list at all. Each row names the
      source that justifies it, and a correction may be scoped to one `reference.id` where that is
      the only place it was verified;
    * **alpha and beta swapped** -- 149 records have a beta-chain gene in ``mhc.a`` and an
      alpha-chain gene in ``mhc.b``. The gene symbol says which chain it is, so this needs no
      judgement.

    #564 (``DPA`` vs ``DPA1``, ``DRA1`` vs ``DRA``) is already fixed in the corpus: zero chunk rows
    have a malformed class-II gene symbol. The issue is stale.
    """
    rows: list[dict[str, object]] = []
    patch = _mhc_patch(root or Paths.discover().root)

    for column in MHC_COLUMNS:
        if column not in df.columns:
            continue
        for value, replacement, scope in patch:
            hit = pl.col(column) == value
            if scope != "*" and "reference.id" in df.columns:
                hit = hit & (pl.col("reference.id") == scope)
            n = df.filter(hit).height
            if n:
                rows.append({"issue": "mhc.dict", "column": column, "from": value,
                             "to": replacement, "rows": n})
                df = df.with_columns(
                    pl.when(hit).then(pl.lit(replacement)).otherwise(pl.col(column)).alias(column))

    if all(c in df.columns for c in MHC_COLUMNS):
        gene = {c: pl.col(c).str.split("*").list.first() for c in MHC_COLUMNS}
        swapped = (gene["mhc.a"].is_in(list(MHC_BETA)) & gene["mhc.b"].is_in(list(MHC_ALPHA)))
        n = df.filter(swapped).height
        if n:
            rows.append({"issue": "mhc-chain-order", "column": "mhc.a,mhc.b",
                         "from": "beta,alpha", "to": "alpha,beta", "rows": n})
            df = df.with_columns(
                pl.when(swapped).then(pl.col("mhc.b")).otherwise(pl.col("mhc.a")).alias("mhc.a"),
                pl.when(swapped).then(pl.col("mhc.a")).otherwise(pl.col("mhc.b")).alias("mhc.b"),
            )
    report = pl.DataFrame(rows, schema={"issue": pl.String, "column": pl.String, "from": pl.String,
                                        "to": pl.String, "rows": pl.Int64})
    return df, report.sort("issue", "column", "from")


# ---------------------------------------------------------------------------------------------
# References
# ---------------------------------------------------------------------------------------------

def harmonise_vocabulary(df: pl.DataFrame,
                         root: Path | None = None) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Apply ``proofreading/vocabulary.tsv``: declared value corrections for the free-text and
    controlled-vocabulary columns (#637, #633).

    `out/reports/lookalikes.tsv` reports values differing from another value only in case or in a
    `-`, `_`, `.` or space, so they never join it - 19 groups over 2 columns when this pass landed.
    That report is the instrument and the table is the fix list; **neither is a gate and neither
    should become one**, because a look-alike pair can be two strains, two serotypes, or one gene
    under two species' symbol conventions. `MBP`/`Mbp` and `G6PC2`/`G6pc2` are both correct: HGNC
    capitalises and MGI title-cases.

    A rewrite is **species-scoped** where a value is right under one convention and wrong under
    another. `POL` is a case variant of `Pol` inside HIV-1 and a separate question inside HCV, so the
    HIV-1 rows are corrected and the 36 HCV rows are not; the table records every such refusal with
    its reason, because a refusal nobody wrote down reads as an omission.

    Runs after :func:`harmonise_references` and before identity is assigned, for the same reason the
    segment passes do: a record is the same record whether the curator wrote `IE1` or `IE-1`, so
    correcting the spelling afterwards would mint a new id for a rename.
    """
    path = (root or Paths.discover().root) / "proofreading" / "vocabulary.tsv"
    empty = pl.DataFrame(schema={"column": pl.String, "species": pl.String, "from": pl.String,
                                 "to": pl.String, "rows": pl.Int64})
    if not path.exists():
        return df, empty
    table = pl.read_csv(path, separator="\t", infer_schema=False, comment_prefix="#").fill_null("")
    rows, exprs = [], {}
    for r in table.sort("column", "species", "from").iter_rows(named=True):
        column = r["column"]
        if column not in df.columns:
            continue
        hit = pl.col(column) == r["from"]
        if r["species"]:
            # `antigen.species` scopes itself; every other column scopes on it.
            scope = column if column == "antigen.species" else "antigen.species"
            if scope not in df.columns:
                continue
            hit = hit & (pl.col(scope) == r["species"])
        n = df.select(hit.sum()).item()
        rows.append({"column": column, "species": r["species"], "from": r["from"], "to": r["to"],
                     "rows": int(n or 0)})
        # Chained into one expression per column, then applied once: two corrections to the same
        # column must both see the *original* value, or `A -> B` followed by `B -> C` would move a
        # cell that was already `B`. Same reasoning as `compare.diff._apply_renames`.
        exprs[column] = pl.when(hit).then(pl.lit(r["to"])).otherwise(
            exprs.get(column, pl.col(column)))
    if exprs:
        df = df.with_columns(*[e.alias(c) for c, e in sorted(exprs.items())])
    report = pl.DataFrame(rows, schema=empty.schema) if rows else empty
    return df, report.filter(pl.col("rows") > 0).sort("column", "from")


def harmonise_references(df: pl.DataFrame,
                         root: Path | None = None) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Replace a non-PMID ``reference.id`` with its PubMed id where one exists (#347).

    Driven by ``proofreading/reference_ids.tsv``, a committed, reviewed input: the build is offline
    and deterministic, so no lookup happens here (CLAUDE.md rule 9).

    #347 asks for DOI and GitHub links to become PMIDs, and most of them cannot. Measured over the
    30,977 records whose reference is not a PMID: the 10x application note is 20,358 of them and is a
    vendor note with no PMID; the eight ``github.com/antigenomics/vdjdb-db/issues/*`` references are
    4,366 records of direct submission, where the issue is the reference; 42 are
    ``rcsb.org/structure/*`` PDB entries; one is a computer-science preprint PubMed does not index;
    and one medRxiv preprint was never indexed. What is left is 668 records across 3 references,
    including a bioRxiv preprint that acquired a PMID after it was submitted.
    """
    path = (root or Paths.discover().root) / "proofreading" / "reference_ids.tsv"
    if not path.exists() or "reference.id" not in df.columns:
        return df, pl.DataFrame(schema={"from": pl.String, "to": pl.String, "rows": pl.Int64})
    table = pl.read_csv(path, separator="\t", infer_schema=False, comment_prefix="#")
    mapping = dict(zip(table["reference.id"], table["pmid"], strict=True))
    rows = [{"from": old, "to": new, "rows": df.filter(pl.col("reference.id") == old).height}
            for old, new in sorted(mapping.items())]
    df = df.with_columns(pl.col("reference.id").replace(mapping))
    report = pl.DataFrame(rows, schema={"from": pl.String, "to": pl.String, "rows": pl.Int64})
    return df, report.filter(pl.col("rows") > 0).sort("from")
