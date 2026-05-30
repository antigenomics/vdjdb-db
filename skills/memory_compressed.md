# VDJdb Skills — Memory (Compressed)

**Last compressed:** 2026-05-26 (rev 3)
**Covers sessions through:** 2026-05-27
**Full log:** `skills/memory.md`

---

## Current State

### Chunks In Progress

| Chunk file | Source | Stage | Blocker |
|---|---|---|---|

*None.*

### Recently Completed Chunks

| Chunk file | PMID | Completed date | Notes |
|---|---|---|---|
| `PMID_41842944.txt` | 41842944 | 2026-05-27 | NF9/HLA-A*24:02; 285 rows; vaccinated donors; tetramer-sort; Sanger |
| `PMID_42125653.txt` | 42125653 | 2026-05-27 | QI9/HLA-A*24:02; 168 rows; mixed vaccinated/convalescent |
| `PMID_40877317.txt` | 40877317 | 2026-05-27 | KF9/HLA-C*12:02; 166 rows (52 paired + 114 beta-only); convalescent |

---

## Active Decisions

| Decision | Date | Summary |
|---|---|---|
| `imgt_alleles.tsv.gz` is primary gene authority | 2026-05-26 | Allele-level (one row per allele); supersedes `patches/IGM_nomenclature_table.tsv` |
| `mhc_alleles.tsv.gz` is primary HLA authority (human only) | 2026-05-26 | Non-human MHC rules documented in `proofreading/mhc.md` |
| Both TSV files gzip-compressed; use `gzip -dc` in queries | 2026-05-26 | User requested; 286K+2.2M → 29K+322K |
| `meta.subset.frequency` gap left open | 2026-05-26 | Missing from `ChunkQC.py` META_COLUMNS but present in all real chunks; Gap #8 |
| `TRAV23S1→TRAV27` retained despite IMGT table discrepancy | 2026-05-27 | IMGT TRAV table says Arden 23S1→TRAV21; existing VDJdb conversion (→TRAV27) from issues #298/#299 takes precedence |
| Adaptive ImmunoSEQ naming documented + fixed in 4 chunks | 2026-05-27 | Adaptive uses `TCRB`/`TCRA` prefix and zero-padded subgroup/cluster; 17 fixes across 4 files; see `proofreading/imgt.md` §9.2 |
| `PDB_Database.txt` rows 258/259 left as-is despite QC duplicate flag | 2026-05-27 | Differ only in `meta.structure.id` (7sg1 vs 7sg2); root cause is Gap #3 (meta.structure.id not in SIGNATURE_COLS) |

---

## Open Questions

| Question | Date | Context |
|---|---|---|
| Add `meta.subset.frequency` to `ChunkQC.py` `META_COLUMNS`? | 2026-05-26 | Present in all real chunks; looks like a bug |
| Non-human MHC reference file (`mhc_nonhuman.tsv`) needed? | 2026-05-26 | IMGT-MHC database at ebi.ac.uk/ipd/mhc could be source |
| How to handle 10X paired alpha/beta in `/extract`? | 2026-05-26 | Join on barcode; no test case yet |

---

## Key Facts (Memorise These)

### VDJdb column order (from actual `chunks/` files)
```
chunk.id | cdr3.alpha | v.alpha | j.alpha | cdr3.beta | v.beta | d.beta | j.beta |
species | mhc.a | mhc.b | mhc.class | antigen.epitope | antigen.gene | antigen.species |
reference.id | method.identification | method.frequency | method.singlecell |
method.sequencing | method.verification | meta.study.id | meta.cell.subset |
meta.subset.frequency | meta.subject.cohort | meta.subject.id | meta.replica.id |
meta.clone.id | meta.epitope.id | meta.tissue | meta.donor.MHC |
meta.donor.MHC.method | meta.structure.id | [comment]
```
Total: 34 columns (with `comment`); 33 without. `meta.subset.frequency` at position 24.

### Species (exact values, case-sensitive)
`HomoSapiens` | `MusMusculus` | `RattusNorvegicus` | `MacacaMulatta`

### CDR3 canonical form
Standard AA only: `ARNDCQEGHILKMFPSTWYV` | Starts `C` | Ends `F` or `W` | Length ≥ 4

### MHC class rules
- `MHCI`: `mhc.b = B2M` (always literal); `mhc.a` = HLA-A/B/C/E/F/G allele
- `MHCII`: `mhc.a` = α-chain allele, `mhc.b` = β-chain allele (never B2M)

### Reference ID formats
- `PMID:XXXXXXX` (preferred)
- `doi:10.XXXX/...` (no URL prefix)
- `https://www.biorxiv.org/...` (full URL)
- `unpublished: Name YYYY-MM-DD`

### Confidence scores
- Score 3: `meta.structure.id` filled (PDB ID, exactly 4 alphanumeric chars) OR direct binding
- Score 2: verification + good sequencing confidence
- Score 1: no verification OR poor sequencing
- Score 0: no method info

### Reference files
| File | Purpose |
|---|---|
| `proofreading/imgt_alleles.tsv.gz` | IMGT allele authority — **allele-level** (5,342 TR alleles; 492 human; release 2026-05-25_GENEDB_202621-7) |
| `proofreading/imgt.md` | IMGT nomenclature + Arden 1995 nomenclature documentation |
| `proofreading/mhc_alleles.tsv.gz` | HLA allele authority (46,005 rows; IPD-IMGT/HLA 3.64.0) |
| `proofreading/mhc.md` | MHC/HLA naming, non-conventional HLA (serological, old format, high-res, expression suffixes) |
| `proofreading/arden.tsv` | Comprehensive Arden 1995 → IMGT conversion table (~123 rows; TRAV, TRAJ, TRBV, TRBJ; includes CDR3-verified TRAJ entries and TRBJ2 dual-numbering notes) |
| `patches/nomenclature.conversions` | Machine-readable Arden / old-style → IMGT gene name conversions (used by fix scripts) |
| `patches/antigen_epitope_species_gene.dict` | Known epitope → gene/species mappings |
| `py_src/ChunkQC.py` | QC implementation (run from `py_src/` directory) |
| `py_src/ScoreFactory.py` | Confidence score computation |

---

## Known py_src Gaps

| # | Check | Description | Status |
|---|---|---|---|
| 1 | MHC-I/B2M | MHCI rows must have `mhc.b = B2M` | open |
| 2 | MHC-II/B2M mismatch | MHCII rows must NOT have `mhc.b = B2M` | open |
| 3 | PDB ID format | `meta.structure.id` must be exactly 4 alphanumeric chars | open |
| 4 | Comment length | `comment` must be ≤ 140 characters | open |
| 5 | Epitope case | `antigen.epitope` must be uppercase | open |
| 6 | DOI URL prefix | `reference.id` starting with `https://doi.org/` should use `doi:` prefix | open |
| 7 | Single-cell+Sanger | `method.singlecell=yes` + `method.sequencing=sanger` is biologically inconsistent | open |
| 8 | `meta.subset.frequency` | Present in all real chunk files but missing from `META_COLUMNS` in `ChunkQC.py` | open |
| 9 | Pseudogene V/J | Gene `functionality = P` in `imgt_alleles.tsv.gz` | open |
| 10 | Invalid allele | Allele number > `allele_count` in `imgt_alleles.tsv.gz` | open |
| 11 | HLA not in IMGTHLA | Human MHC allele not found in `mhc_alleles.tsv.gz` | open |
| 12 | Unconfirmed HLA | Allele has `confirmed = Unconfirmed` in `mhc_alleles.tsv.gz` | open |

---

## Vocabulary Gaps Pending Resolution

| Field | Proposed term | Paper | Status |
|---|---|---|---|

*None yet.*

---

## Reference Data Update Schedule

| File | Source | Update trigger | Re-extract command hint |
|---|---|---|---|
| `proofreading/imgt_alleles.tsv.gz` | genedb-releases (weekly) | Allele name fails lookup in `imgt_allele_id` | Download latest FASTA from genedb-releases → filter TR headers → gzip |
| `proofreading/mhc_alleles.tsv.gz` | ANHIG/IMGTHLA (quarterly: Jan/Apr/Jul/Oct) | HLA allele name fails lookup | Fetch Allelelist.txt + Allele_status.txt from IMGTHLA Latest → join → gzip |
