# VDJdb Skills — Session Memory

This file is the **full running log** for all Claude Code sessions working on the VDJdb curation pipeline using the `/extract`, `/format`, and `/proofread` skills.

**Rules:**
- Append entries; never delete
- Compress into `memory_compressed.md` when this file exceeds ~300 lines
- Always update the "Chunks in progress" table in `memory_compressed.md` at the end of each session

---

## Session Log

<!-- Format:
### [YYYY-MM-DD] <Session title>
**Skills used:** /extract, /format, /proofread
**Source:** [paper title or PMID or submitter]
**Input files:** [list]
**Output files:** [list]
**Summary:** [1–3 sentences describing what was done and outcome]
-->

### [2026-05-27] Extract/format/proofread: chihiro26 (3 SARS-CoV-2 papers)

**Skills used:** /extract, /format, /proofread
**Source:** Submitter Chihiro; XLSX `Summary of TCR clonotypes 260522 Chihiro.xlsx` + `readme.md` at `/Users/aekoneva/data/vdjdb-papers/chihiro26`
**Input files:** `Summary of TCR clonotypes 260522 Chihiro.xlsx` (3 sheets), `readme.md`
**Output files:**
- `chunks/PMID_41842944.txt` (285 rows — NF9/HLA-A*24:02, vaccinated donors, PMID:41842944)
- `chunks/PMID_42125653.txt` (168 rows — QI9/HLA-A*24:02, mixed vaccinated/convalescent, PMID:42125653)
- `chunks/PMID_40877317.txt` (166 rows — KF9/HLA-C*12:02, convalescent, 52 paired + 114 beta-only, PMID:40877317)
**Summary:** Three SARS-CoV-2 CD8+ tetramer-sort datasets curated from a single submitter XLSX. Sheet 3 had a dual-section structure (paired alpha+beta on the left; beta-only from a second experiment on the right). All three sheets had embedded sub-header rows repeated mid-table. TRBJ/TRBD columns were swapped in sheets 2 and 3 right-section. Two gene cells contained allele+functionality code (`TRAV16*01 F`). One cell had an Excel formula artifact (`TRAJ3+D107:D1082`). 36 rows excluded total; 25 CDR3 sequences flagged as non-canonical (kept in chunks — handled by VDJdb build process).

**Key decisions:**
- `tetramer-sort` used (not `multimer-sort`) — papers explicitly used tetramers; tetramer-sort is most common in VDJdb (~50K rows)
- Non-canonical CDR3s (no leading C, ends with C/L/P) kept in chunks/ — downstream build marks them; NOT moved to chunks_with_unconventional_aa/ (that directory is only for non-20-AA characters)
- Rows with `#`/`X` in CDR3 excluded (non-20-AA → data artefacts)
- TRBJ/TRBD swap corrected by gene name prefix detection
- Sheet 3 right section (beta-only): alpha columns left empty; valid VDJdb rows with only beta chain

**Pitfalls documented:** Excel sub-header row filtering, `*01 F` allele stripping, formula artifacts, TRBJ/TRBD swap, `tetramer-sort` vs `multimer-sort` distinction, CDR3 non-canonical handling.

### [2026-05-27] Adaptive Biotech ImmunoSEQ naming — fix zero-padded gene names

**Skills used:** none (bulk fix + documentation session)
**Source:** GitHub issue #495 follow-up; `chunks/PMID_40577481.txt` user-flagged; tcrdist3 docs https://tcrdist3.readthedocs.io/en/latest/adaptive.html
**Input files:** `chunks/PMID_40577481.txt`, `chunks/goncharov-taa-2020-11-02.txt`, `chunks/goncharov-gluten-2023-05-06.txt`, `chunks/goncharov-taa-2020-11-24.txt`; `proofreading/imgt_alleles.tsv.gz` (verification)
**Output files:**
- 4 `chunks/*.txt` files modified (17 total field changes)
- `proofreading/imgt.md` updated (new Section 9.2: Adaptive ImmunoSEQ naming; Section 10 table expanded)
- `skills/vdjdb-format/SKILL.md` updated (Adaptive detection as step 2 in gene normalisation)
- `skills/vdjdb-proofread/SKILL.md` updated (gene diagnosis flowchart; Adaptive pattern in enhanced validation)
- `skills/vdjdb-extract/SKILL.md` updated (Adaptive listed as input source type)
**Summary:** Fixed Adaptive Biotech ImmunoSEQ-style gene names across 4 chunk files: (1) zero-padded subgroup (`TRBV06-5`→`TRBV6-5`, `TRBV06-8`→`TRBV6-8`); (2) zero-padded cluster (`TRBV7-06`→`TRBV7-6`, `TRBV4-01`→`TRBV4-1`, `TRBV5-04`→`TRBV5-4`, `TRBV11-02`→`TRBV11-2`, `TRBV25-01`→`TRBV25-1`); (3) spurious cluster on single-cluster genes (`TRBV19-01`→`TRBV19`, `TRBV9-01`→`TRBV9`, `TRBV15-02`→`TRBV15`, `TRAV21-01`→`TRAV21`). All conversions verified in `imgt_alleles.tsv.gz`. Documented Adaptive convention in `imgt.md` §9.2 and added conversion algorithm to all three skill files.

### [2026-05-27] Issue #495 — fix Arden/space/dual-usage gene names in chunks

**Skills used:** none (bulk fix + documentation session)
**Source:** GitHub issue #495 (reporter: JamieHeather), `antigenomics/vdjdb-db`
**Input files:** All `chunks/*.txt` files; `patches/nomenclature.conversions`; IMGT TRAV/TRBV nomenclature tables (online); `proofreading/imgt_alleles.tsv.gz`
**Output files:**
- 26 `chunks/*.txt` files modified (530 total field changes)
- `chunks/PMID_24906112.txt` additionally deduplicated (49→46 rows; 3 trailing-space duplicates removed)
- `patches/nomenclature.conversions` updated (29 new Arden→IMGT conversions appended)
- `proofreading/arden.tsv` created (comprehensive Arden→IMGT conversion table, ~123 rows)
- `proofreading/imgt.md` updated (Section 9.1 expanded; TRBJ2 dual-numbering documented; CDR3-verified TRAJ entries; arden.tsv reference added)
**Summary:** Fixed three categories of gene name errors reported in issue #495: (1) trailing/internal whitespace stripped from gene names across all chunks; (2) Arden 1995 nomenclature names with spurious TR-prefix (e.g., `TRAV32S1`→`TRAV25`, `TRAV18S1`→`TRAV24`) converted to IMGT using IMGT TRAV/TRBV nomenclature tables; (3) incomplete dual-usage names corrected (`TRAV14`→`TRAV14/DV4`, `TRAV23`→`TRAV23/DV6`, `TRBV17-1`→`TRBV17`). Three TRAJ Arden names without IMGT table coverage were resolved via CDR3 motif cross-referencing within VDJdb (`TRAJ4S1`→`TRAJ58`, `TRAJ17S8`→`TRAJ48`, `TRAJ1S8`→`TRAJ33`). TRBJ2 cluster-relative numbering (PMID_16237109: `TRBJ2S1`–`TRBJ2S6`) documented separately from global Arden numbering (`TRBJ2S7`–`TRBJ2S13`). TRAV23S1→TRAV27 retained despite IMGT table suggesting TRAV21 (established in previous issues #298/#299).

**Key decisions:**
- `TRAV23S1→TRAV27` retained (existing VDJdb conversion from issues #298/#299 takes precedence over IMGT table discrepancy)
- `PDB_Database.txt` false-positive duplicate (rows 258/259, differ only in `meta.structure.id`) left as-is; root cause is Gap #3 (meta.structure.id not in SIGNATURE_COLS)
- 3 genuine duplicates in `PMID_24906112.txt` removed (were masked by trailing spaces in `j.beta`)

### [2026-05-26] IMGT allele expansion, Arden nomenclature, HLA resolution, gzip
**Skills used:** none (skills update session)
**Source:** N/A
**Input files:** IMGT FASTA `IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+inframeP` (release 2026-05-25_GENEDB_202621-7); ANHIG/IMGTHLA `Allelelist.txt` + `Allele_status.txt` (release 3.64.0); PubMed PMID:8550092, PMID:8550093; `patches/HLA2026.pdf`; `patches/nomenclature.conversions`
**Output files:**
- `proofreading/imgt_alleles.tsv.gz` (re-extracted at allele level; 5,342 unique TR alleles; 492 human)
- `proofreading/mhc_alleles.tsv.gz` (compressed)
- `proofreading/imgt.md` (updated: allele-level format, Arden nomenclature section added)
- `proofreading/mhc.md` (updated: non-conventional HLA naming section, zcat queries, numbered section 9 split)
- All 3 SKILL.md files (updated: .tsv → .tsv.gz, allele-level queries, HLA resolution notes in format skill)
**Summary:** Rebuilt `imgt_alleles.tsv.gz` at allele granularity (one row per allele, e.g., TRBV7-9*01/02/.../08 as separate rows) using the IMGT FASTA file. Added Arden 1995 nomenclature documentation (PMID:8550092/8550093) with conversion table reference. Added non-conventional HLA naming section to mhc.md (serological names A2/B7, old 4-digit format A*0201, high-res 4-field alleles, expression suffixes). Compressed both TSV files with gzip; updated all skill files to use zcat queries.

### [2026-05-26] Initial skills setup
**Skills used:** none (setup session)
**Source:** N/A
**Input files:** N/A
**Output files:**
- `skills/vdjdb-extract/SKILL.md`
- `skills/vdjdb-format/SKILL.md`
- `skills/vdjdb-proofread/SKILL.md`
- `proofreading/imgt.md`
- `proofreading/imgt_alleles.tsv.gz` (IMGT/GENE-DB release 2026-05-25_GENEDB_202621-7; 4,693 TR gene rows)
- `proofreading/mhc.md`
- `proofreading/mhc_alleles.tsv.gz` (IPD-IMGT/HLA release 3.64.0, 2026-04-16; 46,005 allele rows)
- `skills/memory.md`
- `skills/memory_compressed.md`
**Summary:** Created the three VDJdb curation skill files and their reference data. The `proofreading/` directory contains authoritative IMGT and IMGT/HLA nomenclature data extracted from official sources. Identified a pre-existing gap in `ChunkQC.py`: `meta.subset.frequency` is present in all validated chunk files (column 24) but is missing from `META_COLUMNS` in the QC script.

---

## Decisions Made

| Date | Decision | Rationale |
|---|---|---|
| 2026-05-26 | Use `proofreading/imgt_alleles.tsv.gz` as primary gene authority, with `patches/IGM_nomenclature_table.tsv` as secondary fallback | IMGT/GENE-DB via genedb-releases is more complete and up-to-date; covers all TR genes across all species at allele resolution |
| 2026-05-26 | Use `proofreading/mhc_alleles.tsv.gz` (ANHIG/IMGTHLA) for human HLA validation; use `proofreading/mhc.md` for non-human MHC rules | IMGTHLA is the authoritative source for HLA alleles; non-human MHC is not in IMGTHLA so documented as rules |
| 2026-05-26 | `meta.subset.frequency` gap in ChunkQC.py left unfixed for now; documented in proofread skill as Gap #8 | Fixing the script requires testing; flagged for maintainer attention |
| 2026-05-26 | Skill files placed in `skills/<skill-name>/SKILL.md` (project-local) | User specified this path structure; consistent with Claude Code project skill conventions |
| 2026-05-26 | `imgt_alleles.tsv.gz` is allele-level (one row per allele), not gene-level | User requested that all alleles (TRBV7-9*01, *02, etc.) be recorded; sourced from FASTA not GeneList |
| 2026-05-26 | Both reference TSVs compressed with gzip; all queries use `zcat` | User requested compression; reduces storage from 2.5 MB to 351 KB total |

---

## Open Questions

| Date opened | Topic | Status | Notes |
|---|---|---|---|
| 2026-05-26 | Should `meta.subset.frequency` be added to `ChunkQC.py`'s `META_COLUMNS`? | open | Present in all validated chunks; omission appears to be a bug |
| 2026-05-26 | Should non-human MHC alleles be added to a separate reference file (e.g., `proofreading/mhc_nonhuman.tsv`)? | open | Currently only rules in `proofreading/mhc.md`; IMGT-MHC has a database at ebi.ac.uk/ipd/mhc |
| 2026-05-26 | How should `/extract` handle 10X Genomics data where both alpha and beta are in the same row (`filtered_contig_annotations.csv` pairs them by barcode)? | open | Mention in skill but no test case yet |

---

## Data Provenance Tracking

<!-- Format per chunk:
#### Chunk: <filename>
- Source paper: <PMID or title>
- Source files: <list>
- Extraction: <date>, log: <path>
- Format: <date>, log: <path>
- Proofread: <date>, report: <path>
- Status: in-progress | ready-for-review | merged-to-chunks | withheld
- Issues/blockers: <free text>
-->

*No chunks in progress yet.*

---

## Reference ID Registry

Track all reference IDs seen across sessions to prevent duplicates and catch format errors.

| Reference ID | Type | Short title | Status |
|---|---|---|---|
| | | | |

*Registry is empty — populate as chunks are processed.*

---

## Vocabulary Extension Proposals

Novel methods or field values encountered during extraction/formatting that are not yet in the VDJdb spec.

| Date | Field | Proposed value | Example paper | Rationale | Status |
|---|---|---|---|---|---|
| | | | | | proposed / accepted / rejected |

*No proposals yet.*

---

## Known Issues in py_src

Bugs or gaps in the QC scripts identified during proofreading sessions.

| Date | Script | Issue | Suggested fix | Severity | Status |
|---|---|---|---|---|---|
| 2026-05-26 | `ChunkQC.py` | `meta.subset.frequency` (column 24 in all validated chunks) is missing from `META_COLUMNS`; `check_header()` does not require it | Add `"meta.subset.frequency"` to `META_COLUMNS` after `"meta.cell.subset"` | Warning | open |
| 2026-05-26 | `ChunkQC.py` | No check that MHCI rows have `mhc.b = B2M` | Add cross-field validator | Warning | open |
| 2026-05-26 | `ChunkQC.py` | No check that MHCII rows do not have `mhc.b = B2M` | Add cross-field validator | Warning | open |
| 2026-05-26 | `ChunkQC.py` | No check that `meta.structure.id` is exactly 4 alphanumeric characters (valid PDB ID format) | Add `lambda x: bool(re.match(r'^[A-Za-z0-9]{4}$', x)) if pd.notnull(x) else True` | Warning | open |
| 2026-05-26 | `ChunkQC.py` | No check that `comment` field ≤ 140 characters | Add length check | Warning | open |
| 2026-05-26 | `ChunkQC.py` | No check that `antigen.epitope` is uppercase-only | Add `lambda x: x == x.upper() if pd.notnull(x) else True` | Warning | open |
| 2026-05-26 | `ChunkQC.py` | `reference.id` starting with `https://doi.org/` passes validation but should use `doi:` prefix | Add normalisation check | Warning | open |
| 2026-05-26 | `ChunkQC.py` | `method.singlecell = yes` + `method.sequencing = sanger` combination is biologically inconsistent | Add cross-field validator | Warning | open |
