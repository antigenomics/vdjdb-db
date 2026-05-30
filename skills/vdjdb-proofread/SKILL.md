---
name: vdjdb-proofread
description: Run QC scripts on a VDJdb chunk, report every error with a suggested fix, verify the output of previous /extract and /format steps, estimate confidence scores, and flag gaps in current py_src QC coverage.
---

# /proofread — VDJdb Chunk Proofreading Skill

## Purpose

Validate a VDJdb chunk file against all available QC tools in `py_src/`, report every problem with a specific suggested fix, estimate confidence scores, and identify any error patterns that the current scripts do not yet detect. This is the third and final stage: extract → format → **proofread**.

## Invocation

```
/proofread [path-to-tsv]
```

The file may be from `chunks/`, `chunks_unformatted/`, or the output of `/extract` or `/format`.

---

## Step 1 — Pre-Script Structural Check

Before running any Python, verify:

| Check | Pass condition | Fail action |
|---|---|---|
| Header columns | All 31 columns in `ALL_COLS` (from `py_src/ChunkQC.py`) are present | Report missing columns; halt |
| No duplicate column names | Column names are unique | Report duplicates; halt |
| `chunk.id` present | First column named `chunk.id` with sequential integers from 1 | Report and fix |
| `meta.subset.frequency` present | Column 24 exists (not in `ChunkQC.py` `ALL_COLS` — see Gap #8) | Note as absent if missing |
| Encoding | File is valid UTF-8 | Report encoding errors |
| Line endings | Unix LF (not Windows CRLF) | Report; convert with `sed -i 's/\r//'` |
| Delimiter | Tab-separated (TSV, not CSV) | Report mixed delimiters |
| Extra columns | Only columns in canonical order plus optional `comment` | Flag extra columns |

**Canonical column order** (from actual `chunks/` files):
```
chunk.id | cdr3.alpha | v.alpha | j.alpha | cdr3.beta | v.beta | d.beta | j.beta |
species | mhc.a | mhc.b | mhc.class | antigen.epitope | antigen.gene | antigen.species |
reference.id | method.identification | method.frequency | method.singlecell |
method.sequencing | method.verification | meta.study.id | meta.cell.subset |
meta.subset.frequency | meta.subject.cohort | meta.subject.id | meta.replica.id |
meta.clone.id | meta.epitope.id | meta.tissue | meta.donor.MHC |
meta.donor.MHC.method | meta.structure.id | [comment]
```

> Note: `meta.subset.frequency` (column 24) is present in all validated `chunks/` files but is **missing from `py_src/ChunkQC.py`'s `META_COLUMNS`**. This is a known gap (Gap #8 below).

---

## Step 2 — Run ChunkQC

Navigate to `py_src/` before running (ChunkQC.py loads `../patches/IGM_nomenclature_table.tsv` with a relative path):

```python
import sys
sys.path.insert(0, 'py_src/')
import pandas as pd
from ChunkQC import ChunkQC, gene_match_check, alleles_match_check, is_qq_seq_biologically_valid

df = pd.read_csv('<chunk_file>', sep='\t')
qc = ChunkQC(df)
errors = qc.process_chunk()
```

### Validators applied by `ChunkQC.process_chunk()`

**Per-row field validators:**
| Column | Rule | Error code |
|---|---|---|
| `cdr3.alpha` | Standard AA only (`ARNDCQEGHILKMFPSTWYV`), length > 3, or null | `bad cdr3.alpha` |
| `v.alpha` | Starts with `TRAV`, or null | `bad v.alpha` |
| `j.alpha` | Starts with `TRAJ`, or null | `bad j.alpha` |
| `cdr3.beta` | Standard AA only, length > 3, or null | `bad cdr3.beta` |
| `v.beta` | Starts with `TRBV`, or null | `bad v.beta` |
| `d.beta` | Starts with `TRBD`, or null | `bad d.beta` |
| `j.beta` | Starts with `TRBJ`, or null | `bad j.beta` |
| `species` | In `['homosapiens', 'musmusculus', 'rattusnorvegicus', 'macacamulatta']` (case-insensitive) | `bad species` |
| `mhc.a` | Matches `HLA-[A-Z]+[0-9]?\*\d{2}(:\d{2,3}){0,3}` OR does not start with `HLA` | `bad mhc.a` |
| `mhc.b` | Same regex as `mhc.a` | `bad mhc.b` |
| `mhc.class` | Exactly `MHCI` or `MHCII` | `bad mhc.class` |
| `antigen.epitope` | Standard AA, length > 3, or null | `bad antigen.epitope` |
| `antigen.gene` | Not null | `bad antigen.gene` |
| `reference.id` | Starts with `PMID:`, `doi:`, `http://`, `https://`, or contains `unpublished` (or null) | `bad reference.id` |

**Cross-row validators:**
| Check | Error code |
|---|---|
| At least one of `cdr3.alpha`, `cdr3.beta` is non-null | `no.cdr3` |
| `antigen.epitope` is non-null | `no.antigen.seq` |
| Both `mhc.a` AND `mhc.b` are non-null | `no.mhc` |
| Row is not a duplicate of another row (on SIGNATURE_COLS) | `duplicate` |

**SIGNATURE_COLS** (used for duplicate detection):
`cdr3.alpha, v.alpha, j.alpha, cdr3.beta, v.beta, d.beta, j.beta, species, mhc.a, mhc.b, mhc.class, antigen.epitope, antigen.gene, antigen.species, reference.id, meta.study.id, meta.cell.subset, meta.subject.cohort, meta.subject.id, meta.replica.id, meta.clone.id, meta.tissue`

**Extended validators** (call separately):
- `gene_match_check(gene_name)`: gene name (without allele) must exist in `patches/IGM_nomenclature_table.tsv`
- `alleles_match_check(gene_name)`: allele number must not exceed known allele count for that gene
- `is_qq_seq_biologically_valid(aa_seq)`: CDR3 must start with `C` and end with `F` or `W`

---

## Step 3 — Report All Problems

**Report format for every error:**
```
ROW <chunk.id> (row <N>) | COLUMN <field> | VALUE "<value>" | ERROR <error_code>
SUGGESTED FIX: <specific action>
```

**Common fixes:**
| Error | Suggested fix |
|---|---|
| `bad v.alpha` (doesn't start with TRAV) | Check `patches/nomenclature.conversions`; look up in `proofreading/imgt_alleles.tsv.gz`; confirm gene is TCR alpha V-gene |
| `bad species` | Normalise to exact VDJdb value (see `proofreading/mhc.md` section 7); check for typos |
| `bad mhc.a` | Check `proofreading/mhc_alleles.tsv.gz`; ensure `HLA-` prefix and `*` separator; see `proofreading/mhc.md` |
| `bad mhc.class` | Must be exactly `MHCI` or `MHCII` — check capitalisation and spelling |
| `bad reference.id` | Convert to `PMID:`, `doi:`, or full preprint URL format |
| `no.cdr3` | Both CDR3 fields are null — at least one is required; check if data was extracted correctly |
| `no.antigen.seq` | `antigen.epitope` is null — this is required; fill from source or flag row for exclusion |
| `no.mhc` | `mhc.a` or `mhc.b` is null — both required if MHC data is known |
| `duplicate` | Exact duplicate of another row on SIGNATURE_COLS — check if intentional (different samples); if so, fill a differentiating meta field |

**Gene not found in IMGT — diagnosis flowchart:**
1. Does it contain spaces? → strip all spaces
2. Does it start with `TCRB`/`TCRA`/`TCRG`/`TCRD`? → Adaptive prefix; replace with `TRB`/`TRA`/`TRG`/`TRD`
3. Does the subgroup number have a leading zero (e.g., `TRBV06-`)? → strip leading zero → `TRBV6-`
4. Does the cluster number have a leading zero (e.g., `-06`, `-02`, `-01`)? → strip leading zero
5. Does the result still not exist in IMGT? → try dropping the cluster suffix entirely (Adaptive adds `-01` to all genes, but many IMGT genes have no cluster suffix)
6. Does it match `TRxVnSn` (Arden pattern with TR prefix)? → look up in `patches/nomenclature.conversions` / `proofreading/arden.tsv`
7. Does it match `BVnSn` (Arden without TR prefix)? → prepend `TR` and check conversions
8. Still unresolved? → flag, report to curator, ask user

---

## Step 4 — Enhanced Gene Validation Against IMGT

For each non-null V/D/J field, run the following checks using `proofreading/imgt_alleles.tsv.gz`
(more authoritative than `patches/IGM_nomenclature_table.tsv`).

Columns: `species | imgt_gene_id | imgt_allele_id | functionality | region_type | accession`

```bash
# Check gene exists in IMGT (gene-level, strip allele suffix first)
gzip -dc proofreading/imgt_alleles.tsv.gz | awk -F'\t' '$2=="TRBV12-3" && $1=="Homo sapiens" {print "FOUND:", $3, $4; exit}'

# Check a specific allele exists
gzip -dc proofreading/imgt_alleles.tsv.gz | awk -F'\t' '$3=="TRBV12-3*02" {print "FOUND:", $4; exit}'

# List all alleles for a gene (human)
gzip -dc proofreading/imgt_alleles.tsv.gz | awk -F'\t' '$2=="TRBV12-3" && $1=="Homo sapiens" {print $3, $4}'

# Check functionality of a specific allele
gzip -dc proofreading/imgt_alleles.tsv.gz | awk -F'\t' '$3=="TRBV7-9*08" {print $4}'
```

Report:
- Gene name not found in `imgt_gene_id` column (even if prefix is correct) — flag; then apply diagnosis flowchart above
- **Adaptive ImmunoSEQ pattern detected** (starts with `TCRB`/`TCRA`, or has zero-padded subgroup like `TRBV06-`, or has zero-padded cluster like `TRBV7-06`/`TRBV4-01`) — convert per `proofreading/imgt.md` §9.2; see format skill §2 for algorithm
- Gene found but allele has `functionality = P` — flag as biologically suspicious
- Specific allele string not found in `imgt_allele_id` — flag as invalid allele; check `patches/nomenclature.conversions` for Arden names (pattern: `TRxVnSn`, e.g., `TRBV1S1`)
- Arden-style names (containing `S` digit after gene type letter) — look up in `patches/nomenclature.conversions` / `proofreading/arden.tsv` and convert

---

## Step 5 — Canonical CDR3 Biology Check

For every non-null CDR3 (`cdr3.alpha`, `cdr3.beta`), verify:
- Starts with `C`
- Ends with `F` or `W`
- Length is biologically reasonable: 8–20 AA for alpha chain, 10–20 AA for beta chain (flag outside these ranges)

**Critical rule — when to use `chunks_with_unconventional_aa/`:**
> `chunks_with_unconventional_aa/` is reserved **exclusively** for CDR3 sequences containing **non-standard amino acids** — i.e., residues outside the canonical 20 (e.g., `X`, `B`, `Z`, `U`, modified residues, L-amino acid designators). It is **NOT** for CDR3s that start without `C`, end with something other than `F/W`, or have unusual lengths — those are canonical composition failures, not non-standard amino acid failures.

| CDR3 issue | Action |
|---|---|
| Contains non-20-AA character (`X`, `B`, `#`, etc.) | Exclude the row entirely (invalid data — likely sequencing noise or data entry artefact) |
| Does not start with `C` | Attempt germline repair (Step 5a); if repair fails, keep in `chunks/` — flag as non-canonical |
| Ends with residue other than `F`/`W` (e.g., `C`, `L`, `P`) | Attempt germline repair (Step 5a); if repair fails, keep in `chunks/` — flag as non-canonical |
| Unusual length (< 8 or > 20 AA) but valid composition | Keep in `chunks/` — flag in extraction log |
| **Contains modified residues or non-standard chemistry** | Move to `chunks_with_unconventional_aa/` — confirm with user first |

**Report** all non-canonical CDR3s in the proofread log with a count per category; do **not** remove them or move them to `chunks_with_unconventional_aa/` unless they contain non-20-AA characters.

---

## Step 5a — CDR3 Canonical Repair Using V/J Germline Context

When a CDR3 fails the canonical check (missing leading `C`, missing terminal `F`/`W`, or single `F`/`W` where double is expected), attempt repair using V and J germline anchor sequences from the VDJdb chunk database before flagging the row as non-canonical.

**This repair is applicable to data sources that systematically trim CDR3 boundaries** (e.g., IEDB dumps, Adaptive ImmunoSEQ exports, some older publications).

### 5a.1 — J genes with double-terminal residues

Certain J genes encode `FF` or `WW` at the end of the CDR3. If a CDR3 ends with a single `F` or `W` and the assigned J gene is in the double-terminal list, the trimmed terminal residue must be restored.

**Beta chain double-terminal J genes** (end with `FF`):
`TRBJ1-1`, `TRBJ1-4`, `TRBJ2-1`, `TRBJ2-2`

**Alpha chain double-terminal J genes** (end with `FF`):
`TRAJ36`

To confirm: look up the J gene's germline CDR3 ending in the VDJdb chunks:

```bash
# Find all CDR3s assigned to this J gene in existing chunks
grep -h "TRBJ1-1" chunks/*.txt | awk -F'\t' '{print $5}' | sort | uniq -c | sort -rn | head -20
```

**Repair rule for double-terminal:** if `CDR3[-1] == 'F'` (or `'W'`) and `CDR3[-3:-1]` matches the expected penultimate 2 AAs from the germline J (derived from existing VDJdb records for that J gene), append one `F` (or `W`).

### 5a.2 — Missing terminal F/W (single instead of double, or completely absent)

For any CDR3 that does not end with `F`/`W`:

1. Extract the last 2–3 AAs of the CDR3: call them `tail`.
2. Look up all VDJdb records with the same `j.beta` (or `j.alpha`) from existing `chunks/` to determine the expected J-gene anchor sequence (last 2–3 AAs before the terminal F/W).
3. If `tail` matches the expected anchor (case-insensitive), the terminal `F` or `W` was trimmed — **append it**.
4. If the J gene is double-terminal and the CDR3 ends with one `F`, apply 5a.1 to add the second `F`.

```python
import csv, glob, collections

def get_j_anchors(j_gene, chain='beta', n=3):
    """Return Counter of CDR3 endings (last n AAs before terminal F/W) for a given J gene."""
    col = 'cdr3.beta' if chain == 'beta' else 'cdr3.alpha'
    j_col = 'j.beta' if chain == 'beta' else 'j.alpha'
    anchors = collections.Counter()
    for f in glob.glob('chunks/*.txt'):
        with open(f) as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                cdr3 = row.get(col, '').strip()
                j = row.get(j_col, '').strip()
                if j == j_gene and cdr3.endswith(('F', 'W')) and len(cdr3) >= n + 1:
                    anchors[cdr3[-(n+1):-1]] += 1
    return anchors

# Example: check if CDR3 tail matches TRBJ1-4 anchor
anchors = get_j_anchors('TRBJ1-4')
most_common_anchor = anchors.most_common(1)[0][0]  # e.g. 'EKL'
```

If the CDR3's last `n` AAs match the most common anchor for that J gene → append `F` (or `W`, whichever the J gene uses).

### 5a.3 — Missing leading C (V-gene conserved cysteine)

For any CDR3 that does not start with `C`:

1. Extract the first 2 AAs of the CDR3: call them `head`.
2. Look up all VDJdb records with the same `v.beta` (or `v.alpha`) from existing `chunks/`.
3. Determine the expected 2 AAs immediately following the conserved `C` (positions 2–3 of the germline CDR3 start).
4. If `head` matches positions 2–3 of the expected V-gene anchor, the leading `C` was trimmed — **prepend it**.

```python
def get_v_anchors(v_gene, chain='beta', n=2):
    """Return Counter of CDR3 starts (first n AAs after leading C) for a given V gene."""
    col = 'cdr3.beta' if chain == 'beta' else 'cdr3.alpha'
    v_col = 'v.beta' if chain == 'beta' else 'v.alpha'
    anchors = collections.Counter()
    for f in glob.glob('chunks/*.txt'):
        with open(f) as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                cdr3 = row.get(col, '').strip()
                v = row.get(v_col, '').strip()
                if v == v_gene and cdr3.startswith('C') and len(cdr3) >= n + 1:
                    anchors[cdr3[1:n+1]] += 1
    return anchors

# If CDR3 starts with 'AS' and TRBV6-2 anchor is 'AS' → prepend 'C'
```

### 5a.4 — Repair decision rules

| Condition | Repair | Log entry |
|---|---|---|
| CDR3 missing leading `C`, head matches V anchor | Prepend `C` | `REPAIR leading-C: <old> → <new> (V-anchor match)` |
| CDR3 missing terminal `F`/`W`, tail matches J anchor | Append `F` or `W` | `REPAIR terminal-FW: <old> → <new> (J-anchor match)` |
| J gene is double-terminal, ends with single `F`, penultimate matches J anchor | Append second `F` | `REPAIR double-F: <old> → <new> (double-terminal J)` |
| Head/tail does NOT match any anchor (< 10% frequency) | Do NOT repair | Flag as non-canonical, keep as-is |
| CDR3 is already canonical after repair | Accept repaired form | |
| Repaired CDR3 fails any other QC check | Reject repair, flag row | |

**Always log every repair** — report count of repaired rows by repair type in Step 10 summary.

### 5a.5 — Batch repair script

For large chunks (>100 non-canonical rows), run as a batch rather than row-by-row:

```python
DOUBLE_TERMINAL_J_BETA = {'TRBJ1-1', 'TRBJ1-4', 'TRBJ2-1', 'TRBJ2-2'}
DOUBLE_TERMINAL_J_ALPHA = {'TRAJ36'}

def repair_cdr3(row, v_anchors_cache, j_anchors_cache):
    repaired = dict(row)
    for chain, cdr3_col, v_col, j_col, dbl_j in [
        ('beta', 'cdr3.beta', 'v.beta', 'j.beta', DOUBLE_TERMINAL_J_BETA),
        ('alpha', 'cdr3.alpha', 'v.alpha', 'j.alpha', DOUBLE_TERMINAL_J_ALPHA),
    ]:
        cdr3 = repaired.get(cdr3_col, '').strip()
        v = repaired.get(v_col, '').strip()
        j = repaired.get(j_col, '').strip()
        if not cdr3:
            continue
        # Missing leading C
        if not cdr3.startswith('C') and v:
            anchors = v_anchors_cache.get(v) or get_v_anchors(v, chain)
            v_anchors_cache[v] = anchors
            if anchors and cdr3[:2] == anchors.most_common(1)[0][0]:
                repaired[cdr3_col] = 'C' + cdr3
                cdr3 = repaired[cdr3_col]
        # Missing terminal F/W
        if cdr3 and cdr3[-1] not in 'FW' and j:
            anchors = j_anchors_cache.get(j) or get_j_anchors(j, chain)
            j_anchors_cache[j] = anchors
            terminal = 'F'  # default; inspect J gene name for W-terminal genes
            if anchors and cdr3[-3:] == anchors.most_common(1)[0][0] + terminal:
                pass  # already ends correctly with the tail
            elif anchors and cdr3[-3:-1] in [a for a, _ in anchors.most_common(5)]:
                repaired[cdr3_col] = cdr3 + terminal
                cdr3 = repaired[cdr3_col]
        # Double-terminal J: single F → FF
        if cdr3 and cdr3.endswith('F') and not cdr3.endswith('FF') and j in dbl_j:
            anchors = j_anchors_cache.get(j) or get_j_anchors(j, chain, n=3)
            j_anchors_cache[j] = anchors
            if anchors and cdr3[-4:-1] in [a for a, _ in anchors.most_common(5)]:
                repaired[cdr3_col] = cdr3 + 'F'
    return repaired
```

---

## Step 6 — MHC Consistency Checks (Beyond ChunkQC)

These checks are **not** performed by `ChunkQC.py` — apply them manually:

| Check | Rule | Fix |
|---|---|---|
| MHC-I / B2M | If `mhc.class = MHCI`, then `mhc.b` must be `B2M` | Set `mhc.b = B2M` |
| MHC-II / B2M mismatch | If `mhc.class = MHCII`, then `mhc.b` must NOT be `B2M` | Fill correct β-chain allele |
| mhc.a → class inference | If `mhc.a` starts with `HLA-A/B/C/E/F/G`, class must be `MHCI` | Correct `mhc.class` |
| mhc.a → class inference | If `mhc.a` starts with `HLA-DR/DQ/DP/DO/DM`, class must be `MHCII` | Correct `mhc.class` |
| HLA allele in mhc_alleles.tsv.gz | Human `mhc.a/mhc.b` starting with `HLA-` should exist in `proofreading/mhc_alleles.tsv.gz` | Note if not found |

---

## Step 7 — Cross-Check Previous Pipeline Steps

If extraction/format logs exist for this chunk:
1. Read `<basename>_extraction_log.txt` — verify logged verifications match current file values
2. Read `<basename>_format_log.txt` — verify normalised values were applied correctly
3. **Spot-check 5 random rows**: grep for their CDR3 sequences and epitopes in source files (if available in the session)

---

## Step 8 — Score Estimation

Based on method fields, estimate the confidence score (0–3) each row would receive from `py_src/ScoreFactory.py`. Flag score-0 rows and suggest which fields to fill.

**Score rubric summary:**
| Score | Meaning | Typical criteria |
|---|---|---|
| 3 | Very high | `meta.structure.id` filled (PDB), OR direct binding assay + single-cell sequencing |
| 2 | High | Verified by re-staining/co-culture; single-cell or Sanger with frequency data |
| 1 | Moderate | No verification OR poor sequencing confidence |
| 0 | Low | Missing method information |

**Score shortcut**: any row with a non-null `meta.structure.id` gets score 3 automatically — verify PDB ID is 4 alphanumeric characters.

---

## Step 9 — Identify Gaps in Current QC Coverage

Document any data quality problem that `ChunkQC.py` does NOT currently detect. Use this pre-seeded list as a starting point, then add any new findings:

### Known gaps (pre-seeded from codebase analysis)

| # | Check | Description | Suggested Python validator |
|---|---|---|---|
| 1 | MHC-I/B2M | `MHCI` row where `mhc.b != 'B2M'` | `lambda r: r['mhc.b'] == 'B2M' if r['mhc.class'] == 'MHCI' else True` |
| 2 | MHC-II/B2M mismatch | `MHCII` row where `mhc.b == 'B2M'` | `lambda r: r['mhc.b'] != 'B2M' if r['mhc.class'] == 'MHCII' else True` |
| 3 | PDB ID format | `meta.structure.id` is not exactly 4 alphanumeric characters | `lambda x: bool(re.match(r'^[A-Za-z0-9]{4}$', x)) if pd.notnull(x) else True` |
| 4 | Comment length | `comment` field exceeds 140 characters | `lambda x: len(x) <= 140 if pd.notnull(x) else True` |
| 5 | Epitope case | `antigen.epitope` contains lowercase letters | `lambda x: x == x.upper() if pd.notnull(x) else True` |
| 6 | DOI URL prefix | `reference.id` starts with `https://doi.org/` instead of `doi:` | `lambda x: not x.startswith('https://doi.org/') if pd.notnull(x) else True` |
| 7 | Single-cell + Sanger | `method.singlecell == 'yes'` and `method.sequencing == 'sanger'` is biologically inconsistent | Cross-field check |
| 8 | `meta.subset.frequency` missing from `META_COLUMNS` | Present in all actual chunk files (column 24) but absent from `ChunkQC.py`'s `META_COLUMNS` list | Add `"meta.subset.frequency"` to `META_COLUMNS` after `"meta.cell.subset"` |
| 9 | Pseudogene V/J assignment | Gene name resolves to `functionality = P` in `proofreading/imgt_alleles.tsv.gz` | Query imgt_alleles.tsv.gz |
| 10 | Invalid allele number | Allele number exceeds `allele_count` in `proofreading/imgt_alleles.tsv.gz` | Query imgt_alleles.tsv.gz |
| 11 | HLA allele not in IPD-IMGT/HLA | Human `mhc.a/mhc.b` (starting with `HLA-`) not found in `proofreading/mhc_alleles.tsv.gz` | Query mhc_alleles.tsv.gz |
| 12 | Unconfirmed HLA allele | Allele in `proofreading/mhc_alleles.tsv.gz` with `confirmed = Unconfirmed` | Query mhc_alleles.tsv.gz |

For each new gap found during this session:
- Document: check description, example failing row, suggested Python code
- Note severity: **blocking** (must fix before merging) vs **warning** (note but may accept)
- Record in `skills/memory.md` under "Known Issues in py_src"

---

## Step 10 — Summary Report

Conclude with a structured summary:

```
=== PROOFREADING SUMMARY ===
File: <filename>
Total rows: N
Rows with errors: N
Clean rows: N

Error breakdown:
  bad cdr3.beta:      N rows
  bad mhc.a:          N rows
  no.antigen.seq:     N rows
  [etc.]

Beyond-ChunkQC findings:
  MHC-I/B2M mismatch: N rows
  [etc.]

Score distribution (estimated):
  Score 3: N rows
  Score 2: N rows
  Score 1: N rows
  Score 0: N rows

Proposed QC additions: N new checks identified
  [list titles]

RECOMMENDATION: [ready for chunks/ | fix N issues first | move N rows to chunks_with_unconventional_aa/]
```

Optionally write this to `<input_basename>_proofread_report.txt` if the user requests it.

---

## Reference Files

| File | Role |
|---|---|
| `py_src/ChunkQC.py` | Primary QC implementation; run this first |
| `py_src/ScoreFactory.py` | Confidence score computation |
| `proofreading/imgt_alleles.tsv.gz` | IMGT V/D/J gene authority (beyond ChunkQC) |
| `proofreading/imgt.md` | IMGT nomenclature rules |
| `proofreading/mhc_alleles.tsv.gz` | HLA allele authority (beyond ChunkQC) |
| `proofreading/mhc.md` | MHC/HLA naming rules |
| `patches/IGM_nomenclature_table.tsv` | Secondary IMGT fallback (used by ChunkQC internally) |
| `chunks/` | 175+ validated reference chunks |
| `chunks_with_unconventional_aa/` | Where non-canonical CDR3s go |
| `chunks_negative/` | Chunks that failed QC and were excluded |
| `skills/memory.md` | Running log; append new py_src gaps found here |

---

## Runtime Note

`ChunkQC.py` loads `../patches/IGM_nomenclature_table.tsv` with a **relative path** at import time. The script must be run from inside the `py_src/` directory, or you must patch the path before importing:

```python
import os
os.chdir('py_src/')
from ChunkQC import ChunkQC
```
