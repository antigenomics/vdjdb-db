# Junction anchors and their repair

> **Authority:** `src/vdjdb/curate/anchors.py`, reading arda's germline reference.
> **Read it with:** `vdjdb submission <chunk>`, or `out/reports/anchors.tsv` after a build.
> **Prose:** [`docs/submission.md`](../docs/submission.md), "Junctions that contradict their own
> germline".

## The rule

VDJdb's `cdr3` is **junction space**: Cys104 through Phe/Trp118, **both anchors included**. A junction
is consistent when its first and last residues are the residues the germline of the V and J the record
names actually encodes.

**The anchor is read from that germline. It is never assumed to be Phe or Trp.** Human `TRAJ35*01`
templates `IGFGNVLHC` and mouse `TRAJ7*01` templates `DYSNNRLTL`, so a correct junction on either ends
in Cys or Leu. Measured over the corpus, a fixed "starts with C, ends with F or W" test calls **481
correct chains broken** and separately **misses 125 that are not** consistent. Earlier revisions of this
file stated that fixed rule and derived its anchors from `chunks/` itself, which is the data validating
itself; both are gone.

`vdjdb qc` does not check this. Its sequence rules check the residue alphabet and a minimum length,
which is why a submission exported in IMGT CDR3 space - short one residue at each end - passes them.

## The defects, as reported

One per end, named for the end and the reason:

| Defect | What the sequence has | Repair proposed |
|---|---|---|
| `V absent anchor` / `J absent anchor` | the anchor residue is missing | the germline residue, added |
| `V corrupt anchor` / `J corrupt anchor` | an anchor is present but is not the germline's | the germline residue, substituted |
| `V under-trimmed` / `J under-trimmed` | framework retained past the anchor | trimmed to the anchor |
| `V allele mismatch` / `J allele mismatch` | the sequence matches a **functional sibling allele** of the named gene | the **call**, not the sequence |
| `V unanchored` / `J unanchored` | no germline to compare against - no call, or a call the reference lacks | none; the universal anchor is all that is left |
| `V unexplained` / `J unexplained` | disagrees with the germline and no rule accounts for it | none |

Three germline residues must agree before a sequence repair is proposed, so a repair is an alignment
result rather than a guess at one letter.

## An allele mismatch is evidence about the call

An ORF or pseudogene allele has a non-canonical anchor by definition, and arda records that faithfully:
of 383 J entries over four organisms only 14 have a `templated_aa` not ending in Phe or Trp, and 13 of
those are marked `ORF` or `P`. So a junction disagreeing with a non-functional allele is evidence that
the **call** is wrong, not the sequence.

95 mouse chains name `TRAJ47`, which resolves to the ORF `*01` (`HYANKMIC`), and every one reads
`DYANKMIF` - exactly `TRAJ47*02`, the functional allele. No corpus record reads the `*01` signature.
This is the same defect as [#327](https://github.com/antigenomics/vdjdb-db/issues/327), where 66 % of
explicit `TRAJ24*01` calls carry the `*02` motif, and rewriting the sequence there would destroy the
evidence for it.

A missing Cys104 deserves a second look even where the rest of the record is clean: no TCR folds
without it, so a first residue that is not Cys while the body still aligns to the V is a sequencing or
transcription error rather than a variant.

## Two independent sources for the anchor

The germline of the named segment is one. The corpus is the other, and they disagree in useful places -
a germline that is an ORF allele has lost its anchor, while 285,989 curated chains have not.

**A consensus over already-stable records is not the data validating itself.** The map below requires
**10 records and 50 % agreement** per gene before it will answer, so a single bad record cannot move a
gene's anchor and a gene with thin coverage is left alone. Where the two sources agree, a repair is
well founded. Where they disagree, the disagreement is the finding: mouse `TRAJ47` germline ends in Cys
because `*01` is an ORF, and every one of the 95 corpus chains on it reads the functional `*02`
signature. The corpus is right and the call is wrong.

Use the germline to name the defect and the corpus consensus to sanity-check the residue it proposes.

### 3.2 Build anchor maps from existing chunks

Run once per session to build lookup tables:

```python
import csv, glob, re
from collections import Counter, defaultdict

FF_J_GENES = frozenset({"TRAJ36", "TRBJ1-1", "TRBJ1-4", "TRBJ2-1", "TRBJ2-2"})

def _normalize_gene(g: str) -> str:
    g = g.strip().upper()
    g = re.sub(r"\*\d+$", "", g)
    g = re.sub(r"^TCRA", "TRA", g)
    g = re.sub(r"^TCRB", "TRB", g)
    return g

def build_anchor_maps(chunks_glob="chunks/*.txt"):
    """
    Returns:
        j_map: {normalized_j: (ctx2, terminal)}
               ctx2 = 2-char anchor immediately before the terminal residue(s)
               terminal = "FF" or "F" (or "W" for W-terminal genes)
        v_map: {normalized_v: ctx2}
               ctx2 = 2-char anchor immediately after leading C
    """
    j_ctx: dict[str, Counter] = defaultdict(Counter)
    v_ctx: dict[str, Counter] = defaultdict(Counter)

    for path in glob.glob(chunks_glob):
        with open(path) as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                for cdr3_col, j_col, v_col in [
                    ('cdr3.beta', 'j.beta', 'v.beta'),
                    ('cdr3.alpha', 'j.alpha', 'v.alpha'),
                ]:
                    cdr3 = (row.get(cdr3_col) or '').strip()
                    j    = (row.get(j_col) or '').strip()
                    v    = (row.get(v_col) or '').strip()
                    if len(cdr3) < 5:
                        continue
                    jn = _normalize_gene(j) if j else ''
                    vn = _normalize_gene(v) if v else ''
                    if v and cdr3.startswith('C'):
                        v_ctx[vn][cdr3[1:3]] += 1
                    if j:
                        if jn in FF_J_GENES and cdr3.endswith('FF'):
                            j_ctx[jn][cdr3[-4:-2]] += 1
                        elif jn not in FF_J_GENES and cdr3.endswith(('F', 'W')):
                            j_ctx[jn][cdr3[-3:-1]] += 1

    j_map: dict[str, tuple[str, str]] = {}
    for jn, counts in j_ctx.items():
        total = sum(counts.values())
        top, top_n = counts.most_common(1)[0]
        if total >= 10 and top_n / total >= 0.5:
            terminal = "FF" if jn in FF_J_GENES else "F"
            j_map[jn] = (top, terminal)

    v_map: dict[str, str] = {}
    for vn, counts in v_ctx.items():
        total = sum(counts.values())
        top, top_n = counts.most_common(1)[0]
        if total >= 10 and top_n / total >= 0.5:
            v_map[vn] = top

    return j_map, v_map
```

Require at least 10 records and 50% consensus for a gene to appear in the anchor maps. Genes with insufficient data are not repaired.

## Why the proposal is against the submitted sequence

`arda.cdr3fix` repairs most of these on the way through the build, so the shipped `cdr3` is usually
right and **the chunk keeps the wrong sequence** - the next export of that data is wrong again. The
repair is therefore computed against the value the chunk holds, not the value that ships:
`YLCSSQEGGYGYTFGSG` ships as `YLCSSQEGGYGYTF`, framework trimmed behind the anchor and kept in front
of it.

## Applying one

A repair is a chunk edit: its own branch, its own issue, and a message naming the files, the rows and
the reason. Nothing applies one automatically, and neither should you in bulk -
[`/vdjdb-proofread`](../skills/vdjdb-proofread/SKILL.md) step 3 and
[`/vdjdb-publish`](../skills/vdjdb-publish/SKILL.md) have the procedure.
