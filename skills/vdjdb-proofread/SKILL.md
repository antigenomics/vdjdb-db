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
| **Combined α/β in `mhc.a`** (Gap #13) | `mhc.a` matches `HLA-DXA*/DYB*...` or similar (contains `/`), `mhc.b` is blank — both chains are collapsed into `mhc.a`. Applies to HLA-DQ, HLA-DP, HLA-DR heterodimers. | Split on `/`: `mhc.a` = prefix before slash (e.g., `HLA-DQA1*01:02`); `mhc.b` = `HLA-` + suffix after slash (e.g., `HLA-DQB1*06:02`). Detection: `bool(re.match(r'^HLA-\S+/\S+$', mhc_a)) and mhc_b == ''` |

**Quick scan command** (run before proofreading any MHC-II chunk):
```bash
awk -F'\t' 'NR>1 && $10~/\// && ($11=="" || $11~/^\s*$/) {print NR, $10}' <chunk_file>
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
