# CDR3 Canonical Repair Reference

> **Used by:** `skills/vdjdb-proofread/SKILL.md` Step 5a
> **Patch source:** VDJdb `chunks/` germline anchor statistics

---

## 1. Canonicality Rules

A CDR3 is canonical if:

1. It starts with `C`.
2. It ends with `F` or `W`.

**Exception — double-terminal J genes:** for certain J genes the germline encodes two consecutive terminal residues (`FF`). A CDR3 assigned to one of these J genes must end with `FF` to be canonical.

---

## 2. Double-Terminal J Genes

| Chain | J genes with `FF` terminal |
|---|---|
| Beta | `TRBJ1-1`, `TRBJ1-4`, `TRBJ2-1`, `TRBJ2-2` |
| Alpha | `TRAJ36` |

This set was derived from the distribution of CDR3 terminal residues across VDJdb `chunks/` records. To verify any gene, count CDR3s assigned to it and inspect the terminal residue distribution:

```bash
# Example: confirm TRBJ1-1 ends with FF
grep -h "TRBJ1-1" chunks/*.txt | awk -F'\t' '{print substr($5, length($5)-1)}' | sort | uniq -c | sort -rn | head
```

---

## 3. Repair Algorithm

When a CDR3 fails the canonical check, attempt repair before flagging the row as non-canonical. Repair uses 2-character germline anchor sequences derived from existing VDJdb records for the same V or J gene.

### 3.1 Gene name normalisation (before anchor lookup)

Strip allele suffix and normalise common prefix variants:

```python
import re

def _normalize_gene(g: str) -> str:
    g = g.strip().upper()
    g = re.sub(r"\*\d+$", "", g)   # strip allele
    g = re.sub(r"^TCRA", "TRA", g)  # Adaptive prefix
    g = re.sub(r"^TCRB", "TRB", g)  # Adaptive prefix
    return g
```

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

### 3.3 Apply repair to a single CDR3

```python
def try_fix_cdr3(cdr3: str, v_gene: str, j_gene: str,
                 v_map: dict[str, str],
                 j_map: dict[str, tuple[str, str]]) -> tuple[str, bool]:
    """
    Attempt a one-letter CDR3 boundary fix using V/J germline anchor.
    Returns (fixed_cdr3, was_fixed). At most one repair per call.
    """
    if not cdr3:
        return cdr3, False

    fixed = cdr3
    jn = _normalize_gene(j_gene) if j_gene else ''
    vn = _normalize_gene(v_gene) if v_gene else ''
    is_ff = jn in FF_J_GENES

    # Missing leading C
    if not fixed.startswith('C') and vn in v_map:
        if fixed[:2] == v_map[vn]:
            fixed = 'C' + fixed

    # Missing terminal residue(s)
    j_ctx = j_map.get(jn)
    if j_ctx:
        ctx2, terminal = j_ctx
        if is_ff:
            if not fixed.endswith('FF'):
                if fixed.endswith('F') and len(fixed) >= 3 and fixed[-3:-1] == ctx2:
                    fixed = fixed + 'F'                        # single F → FF
                elif not fixed.endswith(('F', 'W')) and len(fixed) >= 2 and fixed[-2:] == ctx2:
                    fixed = fixed + 'FF'                       # none → FF
        else:
            if not fixed.endswith(('F', 'W')) and len(fixed) >= 2 and fixed[-2:] == ctx2:
                fixed = fixed + terminal                        # none → F (or W)

    return fixed, fixed != cdr3
```

---

## 4. Repair Decision Rules

| Condition | Repair | Log entry |
|---|---|---|
| CDR3 missing leading `C`; first 2 AAs match V anchor | Prepend `C` | `REPAIR leading-C: <old> → <new>` |
| CDR3 missing terminal `F`/`W`; last 2 AAs match J anchor | Append `F` or `W` | `REPAIR terminal-FW: <old> → <new>` |
| J gene is double-terminal; CDR3 ends with single `F`; penultimate 2 AAs match J anchor | Append second `F` | `REPAIR double-F: <old> → <new>` |
| Anchor match fails (< 50% consensus or < 10 records) | Do NOT repair | Flag as non-canonical, keep as-is |
| Repaired CDR3 still fails another QC check | Reject repair | Flag row |

Always log every repair. Report repair counts by type in the proofreading summary.

---

## 5. Batch Repair Script

For chunks with many non-canonical CDR3s:

```python
def repair_chunk(rows: list[dict], j_map, v_map) -> tuple[list[dict], list[str]]:
    """
    Apply try_fix_cdr3 to every row in-place.
    Returns (repaired_rows, repair_log).
    """
    log = []
    for row in rows:
        for cdr3_col, v_col, j_col in [
            ('cdr3.beta',  'v.beta',  'j.beta'),
            ('cdr3.alpha', 'v.alpha', 'j.alpha'),
        ]:
            cdr3 = (row.get(cdr3_col) or '').strip()
            if not cdr3:
                continue
            fixed, was_fixed = try_fix_cdr3(
                cdr3,
                row.get(v_col, ''),
                row.get(j_col, ''),
                v_map, j_map,
            )
            if was_fixed:
                log.append(f"ROW {row.get('chunk.id', '?')} {cdr3_col}: {cdr3!r} → {fixed!r}")
                row[cdr3_col] = fixed
    return rows, log
```
