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

### Column-shift detection

A column shift occurs when the file's header and data rows have different column counts, or when the header is missing the leading `chunk.id` column. This causes every field to map to the wrong column name — so `antigen.species` might contain a submitter name, `antigen.gene` might contain a reference ID, etc.

**Detect before running ChunkQC:**

```python
with open(chunk_file) as f:
    header = f.readline().rstrip('\n').split('\t')

    data_cols = []
    for _ in range(5):
        line = f.readline()
        if not line:
            break
        data_cols.append(len(line.rstrip('\n').split('\t')))
# Check 1: header vs data column count
for n in data_cols:
    if n != len(header):
        print(f"COLUMN SHIFT: header={len(header)} cols, data row={n} cols — delta={len(header)-n}")

# Check 2: first header column should be 'chunk.id'
if header[0] != 'chunk.id':
    print(f"MISSING chunk.id: first column is {header[0]!r}")
```

**Confirm a shift with content-based sanity checks** — even when column counts match, a shift may be present if:

```python
import csv, re

VALID_SPECIES  = {'HomoSapiens','MusMusculus','RattusNorvegicus','MacacaMulatta','GallusGallus'}
VALID_ANTIGENS = re.compile(r'^[ARNDCQEGHILKMFPSTWYV]{4,}$')
AA_ONLY        = re.compile(r'^[ARNDCQEGHILKMFPSTWYV]+$')

with open(chunk_file) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for i, row in enumerate(reader):
        sp  = row.get('antigen.species', '')
        epi = row.get('antigen.epitope', '')
        ref = row.get('reference.id', '')
        # red flags for column shift:
        if sp and sp not in VALID_SPECIES and 'synthetic' not in sp.lower():
            if not any(sp.startswith(p) for p in ('EBV','CMV','HIV','DENV','HCV','HSV','HBV',
                                                    'HTLV','InfluenzaA','YFV','HPV','VSV','M.')):
                print(f"ROW {i+2}: suspicious antigen.species={sp!r} — possible column shift")
        if epi and not VALID_ANTIGENS.match(epi):
            print(f"ROW {i+2}: antigen.epitope={epi!r} contains non-AA chars — possible shift")
        if ref and not (ref.startswith('PMID:') or ref.startswith('doi:') or
                        ref.startswith('http') or 'unpublished' in ref):
            print(f"ROW {i+2}: reference.id={ref!r} — not a valid reference format")
```

**`reference.id` valid formats:**

| Format | Example | Notes |
|---|---|---|
| `PMID:` prefix | `PMID:8906788` | PubMed ID — preferred |
| `doi:` prefix | `doi:10.1038/ncomms3623` | DOI without https |
| `https://` URL | `https://biorxiv.org/content/...` | arXiv, bioRxiv, other preprints |
| `http(s)://` URL | `https://www.rcsb.org/structure/1AO7` | PDB URL for unpublished structures |
| `unpublished` | `unpublished` | Acceptable only if no DOI/PDB exists |

**For PDB structures without a PMID**, use the canonical PDB URL: `https://www.rcsb.org/structure/{PDB_ID_UPPER}`. This is preferable to leaving blank or using `unpublished`.

Optionally stop after the first 10 rows by adding `if i >= 9: break` inside the loop above.

If a shift is confirmed:
1. Report the delta (header has N more columns than data, or vice versa).
2. Identify which extra header columns are spurious (e.g. `submitter`, `optional columns...`).
3. Fix by either: (a) removing ghost header columns so counts match, or (b) prepending a missing `chunk.id` column to data rows.
4. Re-run the full column validation after repair.

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

**Method field validators** (apply after ChunkQC):

| Field | Valid values | Common mistakes to fix |
|---|---|---|
| `method.identification` | tokens from README set + `structural`, comma-separated | Blank is an error (score penalty); `crystal structure` → `structural`; `tetramer sort` → `tetramer-sort`; see rules below |
| `method.singlecell` | `yes` or blank | Any other value (e.g. `no`, `true`, `single-cell`) → blank or `yes` |
| `method.sequencing` | `sanger`, `amplicon-seq`, `rna-seq`, or blank | `illumina`, `nextseq`, `miseq` → `amplicon-seq` (those are platforms, not methods); `Single cell` → `rna-seq` + set `method.singlecell=yes`; `RNA-seq` → `rna-seq`; `amplicon` → `amplicon-seq` |
| `method.verification` | tokens from the README set, comma-separated | Software names (`mixcr`, `cellranger`) → blank; sort methods misplaced here (`tetramer-sort`, `multimer-sort`) → replace with stain form (`tetramer-stain`, `multimer-stain`); `antigen-coated-targets` → `antigen-loaded-targets`; plain text descriptions → blank |
| `method.frequency` | `N/M` (count/total, e.g. `1/56`) or blank | **Percentage format (`%`) is not valid**. A repeated identical `%` value across all clones for one epitope (e.g., `36.1%` on all NP-reactive rows) means it is a group-level figure (% of tetramer+ cells reactive to that epitope) — move to `meta.subset.frequency` and blank `method.frequency`. A varied `%` value per clone from a sequencing experiment (e.g., `2.40%`, `2.48%`) may be a per-clone repertoire fraction — note in extraction log; ideally convert to `N/M` using the paper's denominators, or blank if denominator is unknown. |

**`method.identification` — restoration rules (apply when blank):**

```
VALID_ID_TOKENS = {
    'tetramer-sort', 'dextramer-sort', 'pelimer-sort', 'pentamer-sort',
    'multimer-sort', 'cd8null-tetramer', 'tetramer-umi',
    'antigen-loaded-targets', 'antigen-expressing-targets', 'beads',
    'cultured-T-cells', 'limiting-dilution-cloning', 'structural',
}
BAD_ID = {
    'crystal structure': 'structural',  # non-standard freetext → canonical
    'tetramer sort':     'tetramer-sort',  # missing hyphen
}
```

**Rule 1 — Swap detection:** If `method.verification` is non-blank AND `method.identification` is blank, the values are likely swapped. Move the verification value to identification and blank verification. Verify by checking the paper: identification = how the antigen-specific T cells were found; verification = how cloned TCRs were re-tested.

**Rule 2 — PDB/structural entries:** If `meta.structure.id` is non-null and `method.identification` is blank (and the file is `PDB_Database.txt` or the chunk has no other method information), set both `method.identification = structural` and `method.verification = structural`.

**Rule 3 — Pre-tetramer era (papers published before ~2000):** Tetramers became available in 1996 and were not widespread until ~2000. For blank identification in entries from papers with PMID < ~10,500,000 that report only 1–5 records (T cell clones), use `antigen-loaded-targets,limiting-dilution-cloning`. Use just `limiting-dilution-cloning` if the paper exclusively describes limiting-dilution steps without explicit antigen-stimulation detail.

**Rule 4 — Infer from file context:** If all other rows in the same chunk use a single identification method and the blank row is isolated (no special MHC or epitope anomaly), fill with that method. Verify against the paper PMID if unsure.

**Rule 5 — PubMed abstract lookup:** For ambiguous cases, fetch `https://pubmed.ncbi.nlm.nih.gov/<PMID>/` and look for: "tetramer", "sorted", "FACS" → `tetramer-sort`; "stimulated", "functional", "ELISpot", "IFN-γ", "killing assay" → `antigen-loaded-targets`; "expressing", "transfected", "transformed" → `antigen-expressing-targets`; "clone", "limiting dilution" + pre-2000 paper → `limiting-dilution-cloning`.

```python
VALID_SEQUENCING = {'sanger', 'amplicon-seq', 'rna-seq', ''}
VALID_SINGLECELL = {'yes', ''}
VALID_ID_TOKENS = {
    'tetramer-sort','dextramer-sort','pelimer-sort','pentamer-sort','multimer-sort',
    'cd8null-tetramer','tetramer-umi','antigen-loaded-targets',
    'antigen-expressing-targets','beads','cultured-T-cells',
    'limiting-dilution-cloning','structural',
}
VALID_VERIF_TOKENS = {
    'tetramer-stain','dextramer-stain','pelimer-stain','pentamer-stain',
    'multimer-stain','beads','restimulation','co-culture',
    'antigen-loaded-targets','antigen-expressing-targets','direct','structural',
}
# Canonical software/platform → blank or correct method
BAD_SEQ = {
    'illumina': 'amplicon-seq', 'nextseq': 'amplicon-seq', 'miseq': 'amplicon-seq',
    'amplicon': 'amplicon-seq', 'RNA-seq': 'rna-seq', 'Single cell': 'rna-seq',
}
BAD_ID = {
    'crystal structure': 'structural', 'tetramer sort': 'tetramer-sort',
}
BAD_VERIF = {
    'mixcr': '', 'cellranger': '', 'CTL clone': '',
    'tetramer-sort': 'tetramer-stain', 'multimer-sort': 'multimer-stain',
    'pentamer-sort': 'pentamer-stain', 'dextramer-sort': 'dextramer-stain',
    'antigen-coated-targets': 'antigen-loaded-targets',
}

for row in rows:
    seq = (row.get('method.sequencing','') or '').strip()
    sc  = (row.get('method.singlecell','') or '').strip()
    ver = (row.get('method.verification','') or '').strip()
    if seq not in VALID_SEQUENCING:
        print(f"BAD method.sequencing={seq!r}")
    if sc not in VALID_SINGLECELL:
        print(f"BAD method.singlecell={sc!r}")
    if ver:
        for token in [t.strip() for t in ver.split(',')]:
            if token and token not in VALID_VERIF_TOKENS:
                print(f"BAD method.verification token={token!r}")
```

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

When a CDR3 fails the canonical check (missing leading `C`, missing terminal `F`/`W`, or single `F`/`W` where double is expected), attempt repair using V and J germline anchor sequences derived from existing VDJdb chunks.

Full algorithm, anchor map construction, and batch repair script: **`proofreading/cdr3_repair.md`**.

**Summary of repair rules:**

| Condition | Repair | Log entry |
|---|---|---|
| CDR3 missing leading `C`; first 2 AAs match V anchor | Prepend `C` | `REPAIR leading-C: <old> → <new>` |
| CDR3 missing terminal `F`/`W`; last 2 AAs match J anchor | Append `F` or `W` | `REPAIR terminal-FW: <old> → <new>` |
| J gene is double-terminal; CDR3 ends with single `F`; penultimate 2 AAs match J anchor | Append second `F` | `REPAIR double-F: <old> → <new>` |
| Anchor match fails | Do NOT repair | Flag as non-canonical, keep as-is |
| Repaired CDR3 fails another QC check | Reject repair | Flag row |

**Double-terminal J genes** (CDR3 must end with `FF`): `TRBJ1-1`, `TRBJ1-4`, `TRBJ2-1`, `TRBJ2-2`, `TRAJ36`.

Always log every repair. Report repair counts by type in Step 10 summary.

---

## Step 6a — Antigen Harmonization Trigger

Before MHC checks, scan `antigen.gene` and `antigen.species` for spurious values using the detectors from `/harmonize`:

```python
spurious_gene_rows = [r for r in rows if is_spurious_gene(str(r.get('antigen.gene', '')))]
spurious_species_rows = [r for r in rows if is_spurious_species(str(r.get('antigen.species', '')))]
consistency_issues = check_consistency(rows)  # same epitope → multiple gene/species values
```

Where `is_spurious_gene`, `is_spurious_species`, and `check_consistency` are defined in `/harmonize`.

If **any** of these return non-empty results:
1. Report count and up to 5 examples of each category.
2. Ask: "Run `/harmonize` to fix antigen.gene/species automatically? [y/n]"
3. If yes: run `/harmonize [path]`, then re-run ChunkQC to verify no regressions.

### 6a-i — Blank antigen.species / antigen.gene detection

**Scan for blanks first** (separate from the spurious-value scan above):

```python
blank_species = [r for r in rows if not str(r.get('antigen.species', '')).strip()]
blank_gene    = [r for r in rows if not str(r.get('antigen.gene', '')).strip()]
# gene blank is acceptable only when antigen.species is 'Synthetic'
real_blank_gene = [r for r in blank_gene if r.get('antigen.species', '') != 'Synthetic']
```

Report distinct `(antigen.epitope, reference.id)` pairs for each category.

**Resolution procedure — IEDB lookup:**

1. Collect all unique blank-species/gene epitope sequences from the chunk.

2. **Choose lookup method** — ask the user which is available:
   - **Local dump** (fast, requires the file): ask `"Do you have the IEDB epitope dump? If so, provide the path to epitope_full_v3.tsv.gz"`. Use the path provided.
   - **IEDB API** (no local file needed): query `https://query-api.iedb.org/epitope_search` with `linear_sequence=<EPITOPE>` and parse JSON results.
   - **PubMed MCP** (available in-session): use `mcp__claude_ai_PubMed__get_article_metadata` with the chunk's `reference.id` PMIDs; read the title/abstract for antigen context.

3. **Local dump lookup** (preferred when available):

```python
import gzip, csv
from collections import Counter

def iedb_lookup(dump_path, epitopes):
    target = {e.upper() for e in epitopes}
    results = {}
    with gzip.open(dump_path, 'rt') as f:
        reader = csv.reader(f, delimiter='\t')
        next(reader); next(reader)  # skip 2 header rows
        for row in reader:
            if len(row) < 3: continue
            epi = row[2].strip().upper()
            if epi not in target: continue
            # natural source: cols 9,13; analog/mimotope: cols 24,28
            gene = row[9].strip() or (row[24].strip() if len(row) > 24 else '')
            org  = row[13].strip() or (row[28].strip() if len(row) > 28 else '')
            results.setdefault(epi, []).append((gene, org))
    return {
        e: (Counter(h[0] for h in hits if h[0]).most_common(1)[0][0] if any(h[0] for h in hits) else '',
            Counter(h[1] for h in hits if h[1]).most_common(1)[0][0] if any(h[1] for h in hits) else '')
        for e, hits in results.items()
    }
```

4. **IEDB API lookup** (when no local dump):

```python
import urllib.request, json

def iedb_api_lookup(epitope):
    url = f"https://query-api.iedb.org/epitope_search?linear_sequence={epitope}&limit=5"
    with urllib.request.urlopen(url) as resp:
        data = json.loads(resp.read())
    hits = data.get('results', [])
    if not hits: return ('', '')
    h = hits[0]
    gene = h.get('source_molecule_name', '') or h.get('molecule_parent_name', '')
    org = h.get('source_organism_name', '')
    return gene, org
```

5. Map IEDB organism names to VDJdb `antigen.species` CamelCase:

| IEDB organism | VDJdb antigen.species |
|---|---|
| Homo sapiens | HomoSapiens |
| Mus musculus | MusMusculus |
| Gallus gallus | GallusGallus |
| Human herpesvirus 4 / Epstein-Barr virus | EBV |
| Human herpesvirus 5 / cytomegalovirus | CMV |
| Dengue virus | DENV |
| HIV / Human immunodeficiency virus | HIV |
| Influenza A virus | InfluenzaA |
| Columba livia | ColumbaLivia |
| Manduca sexta | ManducaSexta |
| Synthetic / mimotope (no natural source) | Synthetic |

6. Map IEDB protein names to gene symbols (HGNC for human, MGI for mouse, standard virus gene names). The IEDB `Source Molecule` field (col 9) often contains synonyms — use the canonical gene symbol, not the full protein name.

7. When IEDB has no match: use PubMed MCP to fetch the abstract for the chunk's `reference.id` PMID and infer from title/abstract context (cancer neoantigen study → HomoSapiens; cross-reactive T cell study → source species varies per epitope; autoantigen study → typically HomoSapiens or MusMusculus).

8. For intentionally synthetic peptides (mimotopes, modified sequences, designed peptides): set `antigen.species = Synthetic`; leave `antigen.gene` blank — this is the only valid case for a blank gene.

### 6a-ii — Synthetic casing normalization

VDJdb uses CamelCase for all species values. Normalize any lowercase variant:

```python
for row in rows:
    if row.get('antigen.species', '').lower() == 'synthetic':
        row['antigen.species'] = 'Synthetic'
```

---

## Step 6 — MHC Validation (Beyond ChunkQC)

Run all sub-steps below on every chunk, regardless of species. The order matters: fix structural/naming errors first, then validate allele identity, then cross-check class consistency.

---

### 6.0 — Quick scan (run first, before any manual checks)

```bash
# 1. Any blank mhc.a or mhc.b
awk -F'\t' 'NR>1 && ($10=="" || $11=="") {print NR, $9, $12, $10, $11}' <chunk_file>

# 2. Combined HLA α/β collapsed into mhc.a (slash present, mhc.b blank)
awk -F'\t' 'NR>1 && $10~/\// && $11=="" {print NR, $10}' <chunk_file>

# 3. Missing digit in DP/DQ gene names (HLA-DPA* → HLA-DPA1*, etc.)
awk -F'\t' 'NR>1 && ($10~/^HLA-DP[AB]\*/ || $11~/^HLA-DP[AB]\*/ || $10~/^HLA-DQ[AB]\*/ || $11~/^HLA-DQ[AB]\*/) {print NR,$10,$11}' <chunk_file>

# 4. Spurious digit in DRA gene name (HLA-DRA1* → HLA-DRA*)
awk -F'\t' 'NR>1 && ($10~/^HLA-DRA1\*/ || $11~/^HLA-DRA1\*/) {print NR,$10,$11}' <chunk_file>

# 5. Missing HLA- prefix in human MHCII alleles
awk -F'\t' 'NR>1 && ($10~/^D[PQR][ABMNO]/ || $11~/^D[PQR][ABMNO]/) {print NR,$10,$11}' <chunk_file>

# 6. mhc.b is B2M but mhc.class is MHCII
awk -F'\t' 'NR>1 && $12=="MHCII" && $11=="B2M" {print NR,$10,$11,$12}' <chunk_file>

# 7. mhc.b is not B2M but mhc.class is MHCI
awk -F'\t' 'NR>1 && $12=="MHCI" && $11!="B2M" && $11!="" {print NR,$10,$11,$12}' <chunk_file>
```

Fix all findings from this scan before proceeding.

---

### 6.1 — mhc.class ↔ mhc.a/mhc.b Correspondence

Apply deterministically from `mhc.a` gene prefix. No lookup needed.

**Human (HomoSapiens):**

| `mhc.a` prefix | `mhc.class` | `mhc.b` must be |
|---|---|---|
| `HLA-A`, `HLA-B`, `HLA-C` | `MHCI` | `B2M` |
| `HLA-E`, `HLA-F`, `HLA-G` | `MHCI` | `B2M` |
| `HLA-DRA` | `MHCII` | `HLA-DRB1*xx:xx` (the paired β-chain allele) |
| `HLA-DRB1`–`HLA-DRB5` | `MHCII` | `HLA-DRA*01:01` (monomorphic α-chain) |
| `HLA-DQA1` | `MHCII` | `HLA-DQB1*xx:xx` |
| `HLA-DQB1` | `MHCII` | `HLA-DQA1*xx:xx` |
| `HLA-DPA1` | `MHCII` | `HLA-DPB1*xx:xx` |
| `HLA-DPB1` | `MHCII` | `HLA-DPA1*xx:xx` |

If `mhc.class` is inconsistent with the `mhc.a` prefix, correct `mhc.class` to match the gene.

**Mouse (MusMusculus):**

| `mhc.a` pattern | `mhc.class` | `mhc.b` must be |
|---|---|---|
| `H-2Db`, `H-2Kb`, `H-2Ld`, `H-2Dd`, etc. | `MHCI` | `B2M` |
| `H2-IAb`, `H2-IAd`, `H2-IEd`, etc. | `MHCII` | same as `mhc.a` (VDJdb canonical) |

**Other species:** see §6.4.

---

### 6.2 — Human HLA Allele Validation Against `mhc_alleles.tsv.gz`

For every `mhc.a` and `mhc.b` value that starts with `HLA-` (human entries), validate against `proofreading/mhc_alleles.tsv.gz` (IPD-IMGT/HLA 3.64.0, 46,005 alleles).

**Validation procedure:**

```bash
# 2-field allele (most common in VDJdb): use prefix match
# Example: validate HLA-DPB1*04:01
gzip -dc proofreading/mhc_alleles.tsv.gz | awk -F'\t' '$2 ~ /^HLA-DPB1\*04:01:/' | head -3

# 1-field allele (low resolution): prefix match on antigen group
# Example: validate HLA-DRB1*15
gzip -dc proofreading/mhc_alleles.tsv.gz | awk -F'\t' '$2 ~ /^HLA-DRB1\*15:/' | head -3

# 4-field allele: exact match
gzip -dc proofreading/mhc_alleles.tsv.gz | awk -F'\t' '$2 == "HLA-A*02:01:01:01"'
```

**Interpretation:**
- **Rows returned**: allele (group) exists — check whether any row has `confirmed = Confirmed`
- **No rows returned**: allele does not exist in IPD-IMGT/HLA → flag; apply §11 fixes from `proofreading/mhc.md`; if still absent, escalate (see "Ambiguous allele" below)
- **All rows `Unconfirmed`**: note in proofreading log; do not reject, but flag

**When a human allele is not found or ambiguous — cross-check against existing VDJdb data:**

1. Grep for the same `antigen.epitope` value in all existing `chunks/` files:
   ```bash
   grep -r "<EPITOPE>" /path/to/chunks/ | cut -f10,11,12 | sort | uniq -c | sort -rn
   ```
2. If other chunks use the same epitope with a well-validated allele, adopt that value (same epitope → same MHC restriction is a strong prior).
3. If the allele differs from what the paper reports, note the discrepancy; do NOT silently overwrite — flag and ask the user.
4. If the epitope is novel (no existing VDJdb rows), validate the allele against the paper text and `mhc_alleles.tsv.gz`. If still ambiguous, ask the user.

**Common naming errors to fix before re-validating** (full list in `proofreading/mhc.md` §11):

| Wrong | Correct |
|---|---|
| `HLA-DPA*01:03` | `HLA-DPA1*01:03` |
| `HLA-DPB*04:01` | `HLA-DPB1*04:01` |
| `HLA-DQA*01` | `HLA-DQA1*01` |
| `HLA-DRA1*01` | `HLA-DRA*01` |
| `DPA1*02:02` (no prefix) | `HLA-DPA1*02:02` |

---

### 6.3 — Mouse MHC Validation

Mouse entries do not use the IPD-IMGT/HLA database. Apply these rules instead.

**Rule 1 — B2M for mouse MHCI alleles:**
Any row with `species = MusMusculus`, `mhc.class = MHCI`, and `mhc.a` starting with `H-2` or `H2-` must have `mhc.b = B2M`.
Set it if blank; flag if set to anything else.

**Rule 2 — Self-fill for mouse MHCII alleles:**
Any row with `species = MusMusculus`, `mhc.class = MHCII`, and `mhc.a` starting with `H2-I` or `I-` must have `mhc.b = mhc.a`.
This is the canonical VDJdb convention for mouse class II (the same allele string fills both fields).
Set `mhc.b = mhc.a` if `mhc.b` is blank or inconsistent.

**Rule 3 — Cross-check new allele values against existing `chunks/`:**
For a newly added mouse allele (e.g., `H2-Kb`, `H2-IAd`, `H2-IEb`), verify the exact string matches what is already in VDJdb:

```bash
# List all mouse MHC-I alleles used in existing chunks
cat chunks/*.txt | awk -F'\t' 'NR>1 && $9=="MusMusculus" && $12=="MHCI" {print $10}' | sort | uniq -c | sort -rn | head -20

# List all mouse MHC-II alleles used in existing chunks
cat chunks/*.txt | awk -F'\t' 'NR>1 && $9=="MusMusculus" && $12=="MHCII" {print $10, $11}' | sort | uniq -c | sort -rn | head -20
```

If the new allele string (e.g., `H-2Kb` vs `H2-Kb`) differs from what is already in the database, normalise to the existing form. If the allele itself is novel (new haplotype, new strain), check the paper for the exact designation and note it in the proofreading log.

**Mouse allele normalisation:**

| Wrong form | Correct form | Rule |
|---|---|---|
| `H2-Db` | `H-2Db` | hyphen between H and 2 |
| `IAb` | `I-Ab` | hyphen after I |
| `H-2D^b` | `H-2Db` | no superscript notation |
| `IEb/d` | check paper | ambiguous — ask user |

---

### 6.4 — Other Species (Mamu, Rat, Novel)

For non-human, non-mouse species, apply in order:

1. **Check `proofreading/mhc.md` §7** for the species-specific naming conventions (Mamu, RT1, etc.)
2. **Search existing chunks** for the same species and epitope to find validated allele strings:
   ```bash
   cat chunks/*.txt | awk -F'\t' 'NR>1 && $9=="<Species>" {print $10, $11, $12, $13}' | sort | uniq -c | sort -rn | head -20
   ```
3. **Cross-check mhc.class**: apply the same logic as §6.1 — if the gene name implies class I or II, correct `mhc.class` to match.
4. **If allele is novel or ambiguous**: do a literature search (PubMed via MCP tool) for the epitope + species + MHC combination, then ask the user to confirm before writing the value.
5. **Never guess**: if the species-specific convention is unclear and no VDJdb precedent exists, ask the user explicitly — provide the paper text and proposed value for confirmation.

---

### 6.5 — Post-Fix Verification

After all MHC corrections, re-run scan 6.0 and confirm zero output for each check. Then verify:

```bash
# Confirm mhc.class distribution is sane
awk -F'\t' 'NR>1 {print $12}' <chunk_file> | sort | uniq -c

# Confirm mhc.b = B2M for all MHCI rows
awk -F'\t' 'NR>1 && $12=="MHCI" && $11!="B2M" {print NR,$10,$11,$12}' <chunk_file>

# Confirm no mhc.b = B2M for MHCII rows
awk -F'\t' 'NR>1 && $12=="MHCII" && $11=="B2M" {print NR,$10,$11,$12}' <chunk_file>
```

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
| 13 | Combined MHC-II α/β in `mhc.a` | `mhc.a` contains a `/` (e.g., `HLA-DQA1*01:02/DQB1*06:02`) with `mhc.b` blank — the α-chain and β-chain are collapsed into one field | Split on `/`: `mhc.a` = part before slash (including `HLA-` prefix); `mhc.b` = `HLA-` + part after slash. Check: `re.search(r'^(HLA-\S+?)/([A-Z]\S+)$', mhc_a)` where `mhc_b == ''` |
| 14 | Percentage in `method.frequency` | `method.frequency` contains `%` instead of count/total (e.g., `36.1%` instead of `13/36`). Percentages are not a valid VDJdb frequency format. Particularly suspicious when the same percentage repeats across all clones for a given epitope (indicating it is a group-level statistic, not a per-clone frequency) | **Do not blindly convert to N/M** — the denominator is often unknown from the paper. Check: if the same value repeats for all clones of one epitope, it likely represents the frequency of that epitope-reactive fraction (e.g., % of tetramer-positive cells) and should be moved to `meta.subset.frequency`. If it is truly a per-clone repertoire frequency (e.g., from high-throughput sequencing), retain as a note in the extraction log and leave blank or convert if the count/total can be determined from the paper. |
| 15 | Blank MHC fields not blocked early | Rows with blank `mhc.a` or `mhc.b` can persist unless explicitly scanned pre/post-proofread | Add explicit audit: `((mhc.a == '') or (mhc.b == ''))` and fail proofreading unless a deterministic repair rule is applied |
| 16 | Mouse MHCII missing `mhc.b` | For `MusMusculus` + `MHCII`, rows often have `mhc.a` filled (e.g., `H2-IEd`) and blank `mhc.b`, despite canonical VDJdb representation using the same allele string in both fields in this dataset | Auto-repair validator: `if species == 'MusMusculus' and mhc.class == 'MHCII' and mhc.a and not mhc.b: mhc.b = mhc.a` → **✅ RESOLVED June 2026** (150 rows filled) |
| 17 | MHC-II gene name digit errors | Three related issues: (a) `HLA-DPA*`/`HLA-DPB*`/`HLA-DQA*` missing trailing `1` (correct: `HLA-DPA1*`, `HLA-DPB1*`, `HLA-DQA1*`); (b) `HLA-DRA1*` with spurious `1` (correct: `HLA-DRA*` — DRA has no digit suffix); (c) `DPA1*`/`DPB1*` without `HLA-` prefix. See `proofreading/mhc.md` §11 for scan commands. → **✅ RESOLVED June 2026** (1417 rows across 5 files) |
| 18 | Blank `antigen.species` | `antigen.species` is empty in non-synthetic records. `ChunkQC` currently flags blank `antigen.gene` (code `bad antigen.gene`) but does **not** flag blank `antigen.species`. Detected in 5 chunk files (1841 rows total): PMID_35667687.txt (1359 rows, MusMusculus/G6pc2), PMID_30418433.txt (346 rows, HomoSapiens neoantigens), PMID_31685621.txt (41 rows, HomoSapiens neoantigens), small_datasets_2026-05-29.txt (94 rows, mixed species). → **✅ RESOLVED June 2026** (IEDB lookup + per-epitope mapping applied; see `fix_antigen_fields.py` in session) |
| 19 | Blank `antigen.gene` (non-synthetic) | `antigen.gene` is empty and `antigen.species` is not `Synthetic`. `ChunkQC` flags this as `bad antigen.gene` but provides no repair guidance. Fix via IEDB lookup (see Step 6a-i). Blanks are acceptable only when `antigen.species == 'Synthetic'`. → **✅ RESOLVED June 2026** (same batch as Gap #18) |
| 20 | `antigen.species` casing: `synthetic` vs `Synthetic` | VDJdb uses CamelCase for all species values, so the canonical form is `Synthetic` (capital S). Found 47 records with lowercase `synthetic` across 3 files (PDB_Database.txt, PMID_29275860.txt, PMID_39286976.txt). Validator: `antigen.species.lower() == 'synthetic' and antigen.species != 'Synthetic'`. → **✅ RESOLVED June 2026** (47 rows normalized) |
| 21 | Species non-canonical variants | Several species names used inconsistently: `HIV` (should be `HIV-1`), `HPV-16` (should be `HPV16`), `HPV-18` (should be `HPV18`), `TriticumAestivum` (should be `Wheat`), `MycobacteriumTuberculosis` (should be `M.tuberculosis`). All are now in `proofreading/species_aliases.tsv`. Scan: `antigen.species not in CANONICAL_SPECIES_SET` — use `is_spurious_species()` from `/harmonize`. → **✅ RESOLVED June 2026** (67 rows across 6 files) |
| 22 | Blank antigen fields in PDB_Database | `antigen.species` or `antigen.gene` blank for PDB structural entries. These genuinely unknown values should be filled with `"Unknown"` to signal curated-but-unknown status vs unchecked blanks. → **✅ RESOLVED June 2026** (8 rows in PDB_Database.txt) |
| 23 | Blank/unpublished `reference.id` in PDB_Database | 43 rows had blank or `"unpublished"` `reference.id`. For PDB structures with no publication, use the canonical PDB URL: `https://www.rcsb.org/structure/{PDB_ID_UPPER}`. → **✅ RESOLVED June 2026** (43 rows fixed) |
| 24 | Gene name for viral nucleocapsid cross-species | `antigen.gene = "Nucleocapsid"` used for both SARS-CoV-2 (canonical VDJdb: `Nucleocapsid`) and InfluenzaA (canonical VDJdb: `NP`). The two should not be conflated. For InfluenzaA rows, `Nucleocapsid` → `NP`. → **✅ RESOLVED June 2026** (104 rows in PMID_31811120.txt) |
| 25 | SARS-CoV mislabeled as SARS-CoV-2 | PMID_34793243.txt contained 2 SARS-CoV-2 epitopes (Spike + ORF1ab) labeled `SARS-CoV`. Confirmed via PubMed abstract (paper is about SARS-CoV-2 CD8 T cells). PMID_38866784.txt retains `SARS-CoV` deliberately — that paper explicitly studies cross-reactive T cells targeting SARS-CoV from COVID-19 convalescents. Always verify SARS-CoV vs SARS-CoV-2 via PubMed abstract before relabeling. → **✅ RESOLVED June 2026** (2 rows) |

### Resolved gaps

| # | Resolution date | Details |
|---|---|---|
| 13 | June 2026 | Combined HLA α/β collapses (e.g., `HLA-DQA1*01:02/DQB1*06:02`) detected in 628 rows across 4 files (PMID_30541895, PMID_33837283, PMID_35675811, small_datasets_2026-05-29). Split using regex `r'^(HLA-\S+?)/([A-Z]\S+)$'`: mhc.a = part before slash, mhc.b = `HLA-` + part after slash. **✅ FULLY RESOLVED** |
| 14 | In progress | Percentage-format frequencies (e.g., `36.1%`) identified in 36 files; flag them as likely group-level statistics when identical across all clones of one epitope. Do not auto-convert to `N/M` unless the denominator is known. If the same percentage repeats for all clones of an epitope, move it to `meta.subset.frequency` and blank `method.frequency`; if it varies per clone from a sequencing experiment, keep the curation note and resolve from the source paper when possible. |
| 15 | June 2026 | Blank mhc.a/mhc.b audit: implemented pre/post-proofread checks. Mouse MHCII blanks resolved (150 rows). Human MHCII blanks resolved (1048 rows total: 153 DRB→HLA-DRA*01:01 fill, 467 DP pairs from PubMed-validated canonical pairings, 428 DQ pairs from narcolepsy literature). **✅ FULLY RESOLVED** |
| 16 | June 2026 | Mouse MHCII self-fill rule applied to 150 rows; validator confirmed all resolved. **✅ FULLY RESOLVED** |
| 17 | June 2026 | MHC-II gene name digit errors: (a) `HLA-DPA*` / `HLA-DPB*` / `HLA-DQA*` missing trailing `1` — fixed 526+526+16 rows in PMID_35750048.txt, nguyen-etal-2023.txt, goncharov-taa-2022-01-27.txt; (b) `HLA-DRA1*01` with spurious `1` — fixed 341 rows in drlcook-etal-2020-02-01.txt → `HLA-DRA*01` (gene is `DRA`, not `DRA1`); (c) missing `HLA-` prefix (`DPA1*02:02`, `DPB1*05:01`) — fixed 4+4 rows in PMID_37418020.txt. See `proofreading/mhc.md` §11. **✅ FULLY RESOLVED** |

### Recommended fills for future submissions

When encountering blank `mhc.b` in human MHCII rows:

**DRB-only rows** (mhc.a matches `HLA-DRB[1-5]*`, mhc.b blank):
→ Set `mhc.b = HLA-DRA*01:01` (canonical, population wild-type ~97%)

**DPA-only rows** (mhc.a matches `HLA-DPA1*`, mhc.b blank):
→ Do NOT auto-fill; requires paper-specific allele information

**DPB-only rows** (mhc.a matches `HLA-DPB1*04:01`, mhc.b blank):
→ Set `mhc.b = HLA-DPA1*01:03` if paper is COVID vaccine/post-COVID study (1034+ confirmed instances in VDJdb)
→ Otherwise, require author statement or leave blank for manual curation

**DQB-only rows** (mhc.a matches `HLA-DQB1*`, mhc.b blank):
→ Check paper context:
  - If narcolepsy study with DQB1*06:02 → Set `mhc.a = HLA-DQA1*05:01`, mhc.b = `HLA-DQB1*06:02` (DQ0602 haplotype, 98% of NT1)
  - Otherwise, require author statement or manual curation

**Mouse MHCII rows** (species=MusMusculus, mhc.class=MHCII, mhc.a=H2-IA*/H2-IE*, mhc.b blank):
→ Set `mhc.b = mhc.a` (canonical VDJdb representation)

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
| `proofreading/cdr3_repair.md` | CDR3 canonical repair algorithm and batch script |
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
