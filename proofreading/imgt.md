# IMGT TCR Gene Nomenclature Reference

> **Authority file:** `proofreading/imgt_alleles.tsv.gz`
> **Format:** One row per allele (e.g., TRBV7-9*01, TRBV7-9*02 … are separate rows)
> **Source:** [IMGT/GENE-DB](https://www.imgt.org/genedb/) via [JamieHeather/genedb-releases](https://github.com/JamieHeather/genedb-releases), FASTA `IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+inframeP`
> **Release used:** `2026-05-25_GENEDB_202621-7` (IMGT/GENE-DB 202621, week 21 of 2026)
> **Contents:** 5,342 unique TR allele entries; 492 human TR allele entries
> **Update cadence:** genedb-releases is updated weekly. Re-extract when a gene or allele name is not found in the current file.

---

## 1. What is IMGT?

IMGT® (international ImMunoGeneTics information system®, https://www.imgt.org) is the global reference for immunogenetics and immunoinformatics, maintained by CNRS, Université de Montpellier. It provides standardised nomenclature, sequences, and 3D structures for immunoglobulins (IG), T-cell receptors (TR), major histocompatibility complex (MHC), and related proteins across species.

**IMGT/GENE-DB** contains, as of the 2026-05-25 release:
- 12,824+ genes across 42 species
- 18,088+ alleles
- Covers all IG and TR loci

---

## 2. TCR Loci and Gene Types

| Locus | Human chromosome | Gene segment types |
|---|---|---|
| **TRA** | 14q11.2 | TRAV, TRAJ, TRAC |
| **TRB** | 7q34 | TRBV, TRBD, TRBJ, TRBC |
| **TRD** | 14q11.2 (embedded within TRA locus) | TRDV, TRDD, TRDJ, TRDC |
| **TRG** | 7p14 | TRGV, TRGJ, TRGC |

> Note: The TRD locus is physically located inside the TRA locus. Some V genes are shared between TRA and TRD rearrangements (dual-usage TRAV/TRDV genes, e.g., TRAV29/DV5).

---

## 3. Gene and Allele Naming Structure

### Format
```
<LOCUS><TYPE><SUBGROUP>[-<CLUSTER>][*<ALLELE_NUMBER>]
```

### Components
| Component | Meaning | Example |
|---|---|---|
| Locus prefix | TRA, TRB, TRD, TRG | TRB |
| Type letter | V=Variable, D=Diversity, J=Joining, C=Constant | V |
| Subgroup number | Gene family (1-based integer) | 12 |
| `-<Cluster>` | Distinct gene within the family (optional) | -3 |
| `*<Allele>` | Allele number (2-digit, zero-padded); `*01` = reference | *02 |

### Examples
| Name | Meaning |
|---|---|
| `TRAV1-1` | TRA-locus, Variable gene, subgroup 1, cluster 1 — any allele accepted |
| `TRAV1-1*01` | Same gene, reference allele |
| `TRAV1-1*02` | Same gene, second allele |
| `TRBV12-3*02` | TRB V gene subgroup 12, cluster 3, allele 02 |
| `TRAJ45*01` | TRA J gene number 45, reference allele |
| `TRBD1` | TRB D gene number 1 (D genes often have no subgroups) |

### VDJdb allele recording convention
- Chunks may record just the **gene name** (`TRBV7-9`) when the exact allele is unknown, or a **specific allele** (`TRBV7-9*01`) when typing data are available
- Both forms are valid — the proofreading check strips `*NN` before gene-level lookup
- If an allele is given, verify it exists in `imgt_alleles.tsv.gz` as a full `imgt_allele_id` entry

### Special designations
- **Unmapped genes**: marked with `S` (sequential) — e.g., `TRBVS1` (position in genome not yet determined)
- **Orphon genes**: located outside the main locus, marked with `/OR`
- **Dual-usage genes**: TRAV genes that can rearrange with both TRAJ and TRDD/TRDJ — written as `TRAV<n>/DV<n>`

---

## 4. Gene Functionality Codes

| Code | Full name | Meaning |
|---|---|---|
| `F` | Functional | Protein-coding, no known defect |
| `(F)` | Functional with note | Functional but mapped to two loci or with ambiguity |
| `[F]` | Inferred functional | Rearranged/expressed; functionality inferred |
| `ORF` | Open Reading Frame | Intact coding sequence but non-functional due to regulatory/structural issue |
| `P` | Pseudogene | Contains stop codons, frameshifts, or other defects |

> **VDJdb usage note:** VDJdb records should use Functional (`F`, `(F)`, `[F]`) or ORF genes. Pseudogene (`P`) assignments should be flagged during proofreading — they likely indicate an error in gene assignment.

---

## 5. Species in IMGT/GENE-DB Relevant to VDJdb

| VDJdb species name | IMGT species name | Common name |
|---|---|---|
| `HomoSapiens` | Homo sapiens | Human |
| `MusMusculus` | Mus musculus | Mouse |
| `RattusNorvegicus` | Rattus norvegicus | Rat |
| `MacacaMulatta` | Macaca mulatta | Rhesus macaque |

The `imgt_alleles.tsv.gz` file contains TR allele entries for all four species plus many others. Filter by the `species` column when validating.

---

## 6. Human TCR Allele Counts (from `imgt_alleles.tsv.gz`, release 2026-05-25)

| Gene type | Allele count (Homo sapiens) |
|---|---|
| TRAV | ~100 alleles across ~70 genes |
| TRAJ | ~61 alleles (mostly 1 per gene) |
| TRAC | 1 allele |
| TRBV | ~250+ alleles across ~65 genes |
| TRBD | 2 alleles |
| TRBJ | ~14 alleles (mostly 1 per gene) |
| TRBC | 2 alleles |
| TRDV | ~20 alleles across ~14 genes (many shared with TRAV) |
| TRDD | 3 alleles |
| TRDJ | 4 alleles |
| TRDC | 1 allele |

---

## 7. Allele Numbering Rules

1. Reference allele is always `*01`
2. Subsequent alleles numbered `*02`, `*03`, etc. in order of discovery
3. Alleles for a given gene can be found by filtering `imgt_allele_id` in `imgt_alleles.tsv.gz` where `imgt_gene_id` matches
4. If an allele string (e.g., `TRBV7-9*09`) is not in `imgt_allele_id`: either the allele doesn't exist or the gene name is wrong

---

## 8. How to Use `imgt_alleles.tsv.gz` for Validation

### Column reference
```
species  imgt_gene_id  imgt_allele_id  functionality  region_type  accession
```

- `species`: IMGT species name (e.g., `Homo sapiens`)
- `imgt_gene_id`: bare gene name without allele (e.g., `TRBV7-9`)
- `imgt_allele_id`: full allele name (e.g., `TRBV7-9*01`)
- `functionality`: `F`, `(F)`, `[F]`, `ORF`, `P`
- `region_type`: `V-REGION`, `J-REGION`, `D-REGION`, `C-REGION` (or exon for constant genes)
- `accession`: IMGT/LIGM-DB accession number

### Validation queries (decompress on-the-fly with zcat)

**Check if a specific allele exists (human):**
```bash
gzip -dc proofreading/imgt_alleles.tsv.gz | awk -F'\t' '$1=="Homo sapiens" && $3=="TRBV7-9*02"'
```

**Find all alleles for a gene (human):**
```bash
gzip -dc proofreading/imgt_alleles.tsv.gz | awk -F'\t' '$1=="Homo sapiens" && $2=="TRBV7-9" {print $3, $4}'
```

**Check if gene name is valid (any species):**
```bash
gzip -dc proofreading/imgt_alleles.tsv.gz | awk -F'\t' '$2=="TRBV7-9"' | head -3
```

**Check functionality of a specific allele:**
```bash
gzip -dc proofreading/imgt_alleles.tsv.gz | awk -F'\t' '$3=="TRBV7-9*08" {print $4}'
```

### Validation rules for VDJdb
1. **Gene-level check** (strip `*NN` suffix): gene must be in `imgt_gene_id` column for the correct species
2. **Allele-level check** (if allele given): full `gene*allele` must be in `imgt_allele_id` column
3. **Functionality check**: flag if `functionality` is `P`; flag with note if `(F)` or `[F]`
4. For human TCR genes: name must start with `TRAV`, `TRAJ`, `TRBV`, `TRBD`, `TRBJ`, `TRDV`, `TRDD`, `TRDJ`, `TRGV`, `TRGJ`

---

## 9. Older Nomenclatures and Conversion

### 9.1 Arden Nomenclature (1995)

**Reference papers:**
- Human: Arden B, Clark SP, Kabelitz D, Mak TW. *"Human T-cell receptor variable gene segment families."* Immunogenetics. 1995;42(6):455-500. **PMID:8550092** doi:10.1007/BF00172176
- Mouse: Arden B, Clark SP, Kabelitz D, Mak TW. *"Mouse T-cell receptor variable gene segment families."* Immunogenetics. 1995;42(6):501-30. **PMID:8550093** doi:10.1007/BF00172177

**How Arden names work:**
The Arden papers named V gene segments by subfamily and discovery order: `TR<LOCUS>V<SUBFAMILY>S<SEQUENTIAL>`. The same pattern applies to J genes. Sequential numbering is within each subfamily:

| Arden name | Current IMGT name | Notes |
|---|---|---|
| `TRBV1S1` | `TRBV9` | — |
| `TRBV2S1` | `TRBV20-1` | — |
| `TRBV4S1` | `TRBV29-1` | — |
| `TRBV6S4` | `TRBV7-9` | — |
| `TRBV6S7` | `TRBV7-1` | — |
| `TRBV7S1` | `TRBV4-1` | — |
| `TRAV23S1` | `TRAV27` | ⚠ Discrepancy: IMGT TRAV table maps Arden 23S1→TRAV21; existing VDJdb conversion (→TRAV27) retained from issues #298/#299 |
| `TRAV15S1` | `TRAV5,TRAV13-2` | Ambiguous primer — comma-separated list |
| `TRAV6S1` | `TRAV14/DV4` | Dual-usage gene |
| `TRAV18S1` | `TRAV24` | Sequence-verified (issue #495) |
| `TRAV32S1` | `TRAV25` | Sequence-verified (issue #495) |

**Comprehensive conversion table:** `proofreading/arden.tsv`
All known human Arden→IMGT conversions in TSV format with columns: `arden_name | imgt_name | locus | segment | notes | source`. Use this as the primary lookup; `patches/nomenclature.conversions` is the machine-readable form used by fix scripts.

**TRAJ conversion — CDR3-verified entries:**
The IMGT TRAJ nomenclature table was not accessible during issue #495 processing. Three TRAJ Arden names were resolved via CDR3 sequence cross-referencing within VDJdb:

| Arden name | IMGT name | CDR3 evidence |
|---|---|---|
| `TRAJ4S1` | `TRAJ58` | CDR3 `CAVSKANGSRLTF` — RLTF motif matches TRAJ58 in all other entries |
| `TRAJ17S8` | `TRAJ48` | CDR3 `CAVLFGNEKLTF` — cross-entry with TRAV13-2/TRAJ48 |
| `TRAJ1S8` | `TRAJ33` | CDR3 `CILRVLGSNYQLIW` — QLIW motif matches TRAJ33 in all other entries |

**TRBJ2 dual-numbering issue:**
Two incompatible numbering schemes exist for TRBJ2 cluster genes:

| Scheme | Range | Used in |
|---|---|---|
| Arden 1995 **global** numbering | `TRBJ2S7`–`TRBJ2S13` | Most 1990s papers; the S-numbers are globally sequential across J1+J2 (J1 has 6 members, so J2 starts at S7) |
| **Cluster-relative** numbering | `TRBJ2S1`–`TRBJ2S6` | PMID:16237109 and some later papers; S1 = first member of J2 cluster = TRBJ2-1 |

Both map to the same IMGT genes:

| Global Arden | Cluster-relative Arden | IMGT name |
|---|---|---|
| `TRBJ2S7` | `TRBJ2S1` | `TRBJ2-1` |
| `TRBJ2S8` | `TRBJ2S2` | `TRBJ2-2` |
| `TRBJ2S9` | `TRBJ2S3` | `TRBJ2-3` |
| `TRBJ2S10` | `TRBJ2S4` | `TRBJ2-4` |
| `TRBJ2S11` | `TRBJ2S5` | `TRBJ2-5` |
| — | `TRBJ2S6` | `TRBJ2-6` |
| `TRBJ2S13` | — | `TRBJ2-7` |

When a paper uses `TRBJ2S1`–`TRBJ2S6`, confirm whether it also uses `TRBJ1S1`–`TRBJ1S5` for J1 (cluster-relative) or `TRBJ1S1`–`TRBJ1S5` and `TRBJ2S7`+ for J2 (global). Check `patches/nomenclature.conversions` — both row sets are present.

**Conversion:** See `proofreading/arden.tsv` for the comprehensive table and `patches/nomenclature.conversions` for the machine-readable mapping. Note that some Arden names map to **multiple** IMGT genes (ambiguous primer sets); when this occurs, document the ambiguity in the extraction log rather than picking one.

### 9.2 Author-specific labels

Papers from the 1990s–2000s sometimes use non-IMGT labels (e.g., `BV20S1`, `AV2S1`, `hVβ2`, `Vβ17`). These usually follow the Arden pattern but may also be author-specific. Always check `patches/nomenclature.conversions` first, then search the literature if not found.

### 9.3 Primer-derived ambiguities

Some older multiplex RT-PCR-based experiments use primers that detect multiple V genes. In this case the `v.beta` field may list multiple IMGT names comma-separated (e.g., `TRBV6-2,TRBV6-3`). This is acceptable in VDJdb — do not reduce to a single gene.

---

## 10. Common Naming Errors and Corrections

| Observed | Likely issue | Action |
|---|---|---|
| `TRBV 7-9` | Space in gene name | Strip space → `TRBV7-9` |
| `TRBV7-9*09` | Allele not in IMGT | Check `imgt_alleles.tsv.gz`; may be a typo of `*07` or `*08` |
| `BV20S1` | Arden nomenclature without TR prefix | Check `nomenclature.conversions` after prepending `TRBV` → `TRBV20S1` is not in conversions; try `TRBV2S1` → `TRBV20-1` |
| `Vβ17` | Informal label | Not mappable without species context; ask user |
| `TRAV14/DV4` | Dual-usage gene, looks wrong | Correct — this is a legitimate IMGT name for a dual-usage gene |
| `TRBVS1` | Unmapped gene (old IMGT notation) | May exist in `imgt_alleles.tsv.gz`; check |

---

## 11. External Resources

- IMGT/GENE-DB web interface: https://www.imgt.org/genedb/
- IMGT Repertoire (human TR loci and genes): https://www.imgt.org/IMGTrepertoire/LocusGenes/listIG_TR/TR/human/Hu_TRgroup.html
- IMGT nomenclature rules: https://www.imgt.org/IMGTScientificChart/Nomenclature/IMGTnomenclature.php
- IMGT TRAV nomenclature table (Arden→IMGT column): https://www.imgt.org/IMGTrepertoire/index.php?section=LocusGenes&repertoire=nomenclatures&species=human&group=TRAV
- IMGT TRBV nomenclature table (Arden→IMGT column): https://www.imgt.org/IMGTrepertoire/index.php?section=LocusGenes&repertoire=nomenclatures&species=human&group=TRBV
- genedb-releases (data source for `imgt_alleles.tsv.gz`): https://github.com/JamieHeather/genedb-releases
- Arden human paper: https://doi.org/10.1007/BF00172176 (PMID:8550092)
- Arden mouse paper: https://doi.org/10.1007/BF00172177 (PMID:8550093)
- Comprehensive Arden conversion table: `proofreading/arden.tsv`
- Machine-readable conversion table: `patches/nomenclature.conversions`
