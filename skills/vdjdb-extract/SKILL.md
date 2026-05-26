---
name: vdjdb-extract
description: Extract TCR:pMHC specificity data from raw source files (papers, supplementary tables, XLS, PDF, 10X output, AIRR-format) and produce a VDJdb-formatted TSV chunk ready for /format and /proofread.
---

# /extract — VDJdb Data Extraction Skill

## Purpose

Extract T-cell receptor antigen-specificity records from arbitrary source files and write a single TSV in VDJdb chunk format. Every value written to the output **must be verified back against the original source** to prevent hallucinations. This skill feeds directly into `/format` and `/proofread`.

## Invocation

```
/extract [path-to-folder-or-file]
```

The input may be a folder or individual file(s) containing any mix of:
- Supplementary Excel/CSV/TSV tables
- PDF manuscripts or supplementary PDFs
- Plain-text files (FASTA, exported tables)
- 10X Genomics `filtered_contig_annotations.csv` / clonotype files
- AIRR-format TSV files (`productive`, `v_call`, `j_call`, `cdr3_aa` columns)
- Jupyter notebooks or scripts used by the authors

If the user limits scope (e.g., "only beta chains", "skip MHC data"), respect that limit and note it in the extraction log.

---

## ⚠️ Absolute Requirements (Non-Negotiable)

### CDR3 sequences
- Must contain **only standard amino acids**: `ARNDCQEGHILKMFPSTWYV`
- Must be **canonical**: starts with `C`, ends with `F` or `W`
- Minimum length: **4 residues**
- If a CDR3 fails the canonical check: **flag it explicitly** in the log and ask the user whether to include it in `chunks_with_unconventional_aa/` instead of `chunks/`. Never silently drop or silently accept.

### Epitope sequences
- Must be **standard amino acids only**
- If authors describe a chemical modification, a non-peptide antigen, or a long peptide pool: **flag prominently** in the log and ask before including

### References — strict enforcement
Only these formats are acceptable in `reference.id`:
| Format | Example | Notes |
|---|---|---|
| `PMID:XXXXXXX` | `PMID:28975614` | Strongly preferred; numeric only after colon |
| `doi:10.XXXX/...` | `doi:10.1016/j.immuni.2023.01.001` | Lowercase `doi:`, no URL prefix |
| Preprint URL | `https://www.biorxiv.org/content/10.1101/2024.01.01.123456` | Full URL |
| Unpublished | `unpublished: Submitter Name YYYY-MM-DD` | For submissions without a publication |

**NEVER invent or guess a PMID.** If uncertain, leave blank and ask the user. DOI and preprint URLs must be quoted exactly from the source.

### Hallucination prevention (mandatory for every extracted value)
After extracting any amino acid sequence, gene name, species name, MHC allele, or reference ID:
1. Run a grep or direct text search in the **original source file** to confirm the exact string is present
2. Log the verification result (found / not found / found with minor variant)
3. If the value **cannot be confirmed** in the source: mark as `[UNVERIFIED]` and do **not** include it without explicit user approval

---

## Step-by-Step Workflow

### Step 1 — Survey the source folder

1. List all files and identify types (PDF, XLS, TSV, CSV, FASTQ, etc.)
2. Note which files are likely to contain: TCR sequences, antigen/epitope data, MHC/HLA data, methods, references
3. Share the inventory with the user before proceeding if the folder contains more than 3 files or has an unclear structure

### Step 2 — Build a cross-reference graph

Source data is often spread across multiple files. Build an explicit graph:
1. Identify all **ID columns** in each file: barcode, clone ID, sample ID, donor ID, clonotype ID, barcode, well ID, etc.
2. Determine which IDs appear in multiple files and can be used to join records
3. Perform the join; log the join keys used

**Ambiguities** (one-to-many links, missing join keys, contradictory values between files) must be logged **immediately**. Do not silently pick one option.

### Step 3 — Extract TCR complex fields

For each record, extract the following. Leave blank if absent — **never use a placeholder string** (`NA`, `N/A`, `null`, `nan`, `-`, `.`).

| VDJdb field | What to look for | Verification |
|---|---|---|
| `cdr3.alpha` | Alpha chain CDR3 amino acid sequence | grep in source |
| `v.alpha` | TRAV gene (IMGT style preferred; note if not IMGT) | grep in source |
| `j.alpha` | TRAJ gene | grep in source |
| `cdr3.beta` | Beta chain CDR3 amino acid sequence | grep in source |
| `v.beta` | TRBV gene | grep in source |
| `d.beta` | TRBD gene (often missing; leave blank) | grep if present |
| `j.beta` | TRBJ gene | grep in source |
| `species` | Organism (`HomoSapiens`, `MusMusculus`, `RattusNorvegicus`, `MacacaMulatta`) | confirm from text |
| `mhc.a` | First MHC chain (e.g., `HLA-A*02:01`, `H-2Db`) | grep in source |
| `mhc.b` | Second MHC chain (`B2M` for MHC-I; β-chain allele for MHC-II) | grep in source |
| `mhc.class` | `MHCI` or `MHCII` | confirm from context |
| `antigen.epitope` | Epitope amino acid sequence | grep in source |
| `antigen.gene` | Antigen gene name (e.g., `pp65`, `NP`, `MART-1`) | grep in source |
| `antigen.species` | Antigen origin (e.g., `CMV`, `InfluenzaA`, `HomoSapiens`) | confirm from text |
| `reference.id` | PMID, DOI, or preprint URL | strict format check + grep |

**At least one of `cdr3.alpha` or `cdr3.beta` must be non-blank per row.**
**Both `mhc.a` and `mhc.b` must be filled if MHC data exists.**

Cross-reference epitopes against `patches/antigen_epitope_species_gene.dict` — if the epitope is already known, use the dict's gene and species values.

### Step 4 — Extract method fields

These fields directly affect the VDJdb confidence score (0–3) computed by `py_src/ScoreFactory.py`. Extract them carefully from the methods section.

| Field | Recognised values | Notes |
|---|---|---|
| `method.identification` | `tetramer-sort`, `dextramer-sort`, `pelimer-sort`, `pentamer-sort`, `antigen-loaded-targets`, `antigen-expressing-targets`, `beads`, `cultured-T-cells`, `limiting-dilution-cloning`, `tetramer-umi`, `cd8null-tetramer` | Multiple values comma-separated, no spaces |
| `method.frequency` | `X/X` (e.g., `7/30`), `X%`, or decimal fraction | Fraction format preferred |
| `method.singlecell` | `yes` if single-cell sequencing was used; blank otherwise | |
| `method.sequencing` | `sanger`, `rna-seq`, `amplicon-seq` | |
| `method.verification` | `tetramer-stain`, `dextramer-stain`, `direct`, `restimulation`, `co-culture`, `antigen-loaded-targets`, `antigen-expressing-targets`, `beads` | |

If the paper uses a method not in the above lists, **do not force it into an existing category**. Log it as a "Novel method" candidate for extending the VDJdb specification.

> **Score shortcut:** If a PDB structure ID is available, record it in `meta.structure.id` — this grants score 3 automatically, bypassing all other scoring logic.

### Step 5 — Extract metadata fields

Fill as many of the 12 meta columns as the source supports. Leave others blank.

| Field | Description |
|---|---|
| `meta.study.id` | Internal study identifier used in the paper |
| `meta.cell.subset` | T cell subset (`CD8+`, `CD4+CD25+`, etc.) |
| `meta.subset.frequency` | Clone frequency within the cell subset |
| `meta.subject.cohort` | Donor cohort (`healthy`, `HIV+`, `CMV-seroneg`, etc.) |
| `meta.subject.id` | Donor/patient identifier |
| `meta.replica.id` | Replicate or timepoint identifier |
| `meta.clone.id` | T cell clone identifier or barcode |
| `meta.epitope.id` | Short epitope label from the paper (e.g., `FL10`) |
| `meta.tissue` | `PBMC`, `spleen`, `TIL`, `TCL`, etc. |
| `meta.donor.MHC` | Donor HLA typing if reported |
| `meta.donor.MHC.method` | HLA typing method if reported |
| `meta.structure.id` | PDB ID if a structure was solved for this complex |
| `comment` | Any important note not captured elsewhere (max 140 characters) |

### Step 6 — Assemble the output TSV

**Canonical column order** (must match exactly):
```
chunk.id
cdr3.alpha  v.alpha  j.alpha  cdr3.beta  v.beta  d.beta  j.beta
species  mhc.a  mhc.b  mhc.class  antigen.epitope  antigen.gene  antigen.species
reference.id
method.identification  method.frequency  method.singlecell  method.sequencing  method.verification
meta.study.id  meta.cell.subset  meta.subset.frequency  meta.subject.cohort  meta.subject.id
meta.replica.id  meta.clone.id  meta.epitope.id  meta.tissue
meta.donor.MHC  meta.donor.MHC.method  meta.structure.id
comment
```

**Formatting rules:**
- Tab-separated, UTF-8, Unix line endings (LF)
- `chunk.id`: sequential integers starting from 1
- Blank fields: truly empty (no quotes, no `NA`, no `-`)
- No trailing whitespace; no quoted fields (TSV, not CSV)
- `comment` is optional — include column only if at least one row has a comment

### Step 7 — Write the extraction log

Write `<output_basename>_extraction_log.txt` containing:

1. **Source inventory**: all files found, their types, and assigned roles
2. **Cross-reference graph**: which ID columns were used to link which files, and join cardinality
3. **Ambiguities**: every unclear/ambiguous case and the decision made (or flagged for user)
4. **Verification results**: for each field, confirm or report failure to verify in source
5. **Unverified values**: any value included at user direction despite not being directly verifiable
6. **Novel methods**: identification/verification methods not in the current VDJdb vocabulary
7. **Excluded records**: rows dropped and the reason (failed canonical check, missing required fields, etc.)
8. **Scope restrictions**: any user-imposed limits on what was extracted

---

## Validation Reference Files

| File | Purpose |
|---|---|
| `proofreading/imgt_alleles.tsv` | **Primary** V/D/J gene ID authority — check all extracted gene names here first |
| `proofreading/imgt.md` | IMGT nomenclature rules and gene structure explanation |
| `proofreading/mhc_alleles.tsv` | **Primary** HLA allele authority — check all extracted human MHC alleles here |
| `proofreading/mhc.md` | MHC/HLA naming rules, class I vs II, non-human conventions |
| `patches/IGM_nomenclature_table.tsv` | Secondary V/D/J fallback (existing repo file) |
| `patches/nomenclature.conversions` | Old-style → IMGT gene name conversions |
| `patches/antigen_epitope_species_gene.dict` | Known epitope → antigen.gene/antigen.species mappings |
| `py_src/ScoreFactory.py` | Confidence score logic and method vocabulary |
| `README.md` | VDJdb column specification |

---

## Output

- **Primary output**: `<PMID_xxxxxxx>_unformatted.txt` — raw VDJdb TSV (gene names not yet IMGT-normalised)
- **Log**: `<PMID_xxxxxxx>_extraction_log.txt` — full provenance and ambiguity record

Suggest the output filename based on PMID if available (e.g., `PMID_40713946_unformatted.txt`), otherwise use the folder name or ask the user.

---

## Next Steps

1. Run `/format` on the output to normalise gene names, species, and MHC alleles to IMGT standard
2. Run `/proofread` to validate against QC scripts in `py_src/`
