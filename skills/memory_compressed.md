# VDJdb Skills — Memory (Compressed)

**Last compressed:** 2026-05-26
**Covers sessions through:** 2026-05-26
**Full log:** `skills/memory.md`

---

## Current State

### Chunks In Progress

| Chunk file | Source | Stage | Blocker |
|---|---|---|---|

*None — pipeline just set up.*

### Recently Completed Chunks

| Chunk file | PMID | Completed date | Notes |
|---|---|---|---|

*None yet.*

---

## Active Decisions

| Decision | Date | Summary |
|---|---|---|
| `imgt_alleles.tsv` is primary gene authority | 2026-05-26 | Supersedes `patches/IGM_nomenclature_table.tsv` for validation; use IGM file as fallback |
| `mhc_alleles.tsv` is primary HLA authority (human only) | 2026-05-26 | Non-human MHC rules documented in `proofreading/mhc.md` |
| `meta.subset.frequency` gap left open | 2026-05-26 | Missing from `ChunkQC.py` META_COLUMNS but present in all real chunks; documented as Gap #8 |

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
| `proofreading/imgt_alleles.tsv` | IMGT gene authority (4,693 TR gene rows; release 2026-05-25) |
| `proofreading/imgt.md` | IMGT nomenclature documentation |
| `proofreading/mhc_alleles.tsv` | HLA allele authority (46,005 rows; IPD-IMGT/HLA 3.64.0) |
| `proofreading/mhc.md` | MHC/HLA nomenclature documentation |
| `patches/nomenclature.conversions` | Old gene name → IMGT conversions |
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
| 9 | Pseudogene V/J | Gene `functionality = P` in `imgt_alleles.tsv` | open |
| 10 | Invalid allele | Allele number > `allele_count` in `imgt_alleles.tsv` | open |
| 11 | HLA not in IMGTHLA | Human MHC allele not found in `mhc_alleles.tsv` | open |
| 12 | Unconfirmed HLA | Allele has `confirmed = Unconfirmed` in `mhc_alleles.tsv` | open |

---

## Vocabulary Gaps Pending Resolution

| Field | Proposed term | Paper | Status |
|---|---|---|---|

*None yet.*

---

## Reference Data Update Schedule

| File | Source | Update trigger |
|---|---|---|
| `proofreading/imgt_alleles.tsv` | genedb-releases (weekly) | Re-extract when a gene name fails lookup |
| `proofreading/mhc_alleles.tsv` | ANHIG/IMGTHLA (quarterly: Jan/Apr/Jul/Oct) | Re-extract after each major release |
