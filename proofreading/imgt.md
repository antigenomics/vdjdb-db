# IMGT TCR Gene Nomenclature Reference

> **Authority file:** `proofreading/imgt_alleles.tsv`
> **Source:** [IMGT/GENE-DB](https://www.imgt.org/genedb/) via [JamieHeather/genedb-releases](https://github.com/JamieHeather/genedb-releases)
> **Release used:** `2026-05-25_GENEDB_202621-7` (IMGT/GENE-DB 202621, week 21 of 2026)
> **Update cadence:** genedb-releases is updated weekly. Re-extract `imgt_alleles.tsv` when onboarding data with gene names not found in the current file.

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

## 3. Gene Naming Structure

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
| `TRAV1-1` | TRA-locus, Variable gene, subgroup 1, cluster 1 (no allele = any allele) |
| `TRAV1-1*01` | Same gene, reference allele |
| `TRAV1-1*02` | Same gene, second allele |
| `TRBV12-3*02` | TRB V gene subgroup 12, cluster 3, allele 02 |
| `TRAJ45*01` | TRA J gene number 45, reference allele |
| `TRBD1` | TRB D gene number 1 (D genes often have no subgroups) |

### Special designations
- **Unmapped genes**: marked with `S` (sequential) — e.g., `TRBVS1` (position in genome not yet determined)
- **Orphon genes**: located outside the main locus, marked with `/OR` — e.g., `TRAV14/DV4` (dual-usage gene)
- **Dual-usage genes**: TRAV genes that can rearrange with both TRAJ and TRDD/TRDJ — written as `TRAV<n>/DV<n>`

---

## 4. Gene Functionality Codes

| Code | Full name | Meaning |
|---|---|---|
| `F` | Functional | Protein-coding, no known defect |
| `(F)` | Functional (mapped to two loci or ambiguous) | Functional but with an annotation note |
| `ORF` | Open Reading Frame | Intact coding sequence but non-functional due to regulatory or structural issue |
| `P` | Pseudogene | Contains stop codons, frameshifts, or other defects that prevent protein production |
| `[F]` | Functional (in square brackets) | Rearranged and expressed sequence where functionality is inferred |

> **VDJdb usage note:** VDJdb records should typically use Functional (`F`) or ORF genes. Pseudogene (`P`) gene assignments should be flagged during proofreading — they may indicate an error in gene assignment.

---

## 5. Species in IMGT/GENE-DB Relevant to VDJdb

| VDJdb species name | IMGT species name | Common name |
|---|---|---|
| `HomoSapiens` | Homo sapiens | Human |
| `MusMusculus` | Mus musculus | Mouse |
| `RattusNorvegicus` | Rattus norvegicus | Rat |
| `MacacaMulatta` | Macaca mulatta | Rhesus macaque |

The `imgt_alleles.tsv` file contains TR gene entries for all four species plus many others. Filter by `species` column when validating.

---

## 6. Human TCR Gene Counts (IMGT/GENE-DB 2026-05-25 release)

| Gene type | Count (Homo sapiens) |
|---|---|
| TRAV | ~70 genes, ~100 alleles |
| TRAJ | ~61 genes |
| TRAC | 1 gene |
| TRBV | ~65 genes, ~200+ alleles |
| TRBD | 2 genes |
| TRBJ | ~14 genes |
| TRBC | 2 genes |
| TRDV | ~14 genes (many shared with TRAV) |
| TRDD | 3 genes |
| TRDJ | 4 genes |
| TRDC | 1 gene |

---

## 7. Allele Numbering Rules

1. Reference allele is always `*01`
2. Subsequent alleles are numbered `*02`, `*03`, etc. in order of discovery
3. Allele numbers within a gene go up to the value in the `allele_count` column of `imgt_alleles.tsv`
4. If an allele number exceeds `allele_count` for its gene, it is invalid — flag during proofreading

---

## 8. How to Use `imgt_alleles.tsv` for Validation

### Column reference
```
species  imgt_gene_id  functionality  gene_definition  allele_count  chromosome  chromosomal_localization  ligm_db_accessions  ncbi_gene_id
```

### Validation queries

**Check if a gene name is valid (any species):**
```bash
grep -P "^.*\t<GENE_NAME>\t" proofreading/imgt_alleles.tsv
```

**Check if a gene name is valid for Homo sapiens:**
```bash
grep -P "^Homo sapiens\t<GENE_NAME>\t" proofreading/imgt_alleles.tsv
```

**Check allele count for a gene:**
```bash
awk -F'\t' '$2 == "TRBV12-3" {print $5}' proofreading/imgt_alleles.tsv
```

**Check functionality of a gene:**
```bash
awk -F'\t' '$2 == "TRBV12-3" {print $3}' proofreading/imgt_alleles.tsv
```

### Validation rules for VDJdb
1. Strip the `*NN` allele suffix before looking up in `imgt_gene_id` column
2. If the gene name (without allele) is NOT in `imgt_gene_id`: invalid gene name — flag
3. If the gene IS in the table but `functionality` is `P`: pseudogene — flag as suspicious
4. If the gene has an allele suffix and allele number > `allele_count`: invalid allele — flag
5. For human TCR genes: name must start with `TRAV`, `TRAJ`, `TRBV`, `TRBD`, `TRBJ`, `TRDV`, `TRDD`, `TRDJ`, `TRGV`, `TRGJ`

---

## 9. Common Naming Errors and Corrections

See `patches/nomenclature.conversions` for a conversion table of old-style gene names to current IMGT names. Common patterns:

| Old-style name | Current IMGT name |
|---|---|
| `TRAV23S1` | `TRAV27` |
| `TRBV1S1` | `TRBV9` |
| S-numbered genes (e.g., `TRBV1S3`) | Look up in conversions file |

If a gene name is not found in `imgt_alleles.tsv` AND not in `nomenclature.conversions`, it may be:
- A typographical error
- An author-specific label (e.g., `BV20S1`) that needs manual mapping
- A gene from a species not yet in IMGT

---

## 10. External Resources

- IMGT/GENE-DB web interface: https://www.imgt.org/genedb/
- IMGT Repertoire (human TR loci and genes): https://www.imgt.org/IMGTrepertoire/LocusGenes/listIG_TR/TR/human/Hu_TRgroup.html
- IMGT nomenclature rules: https://www.imgt.org/IMGTScientificChart/Nomenclature/IMGTnomenclature.php
- genedb-releases (data source for `imgt_alleles.tsv`): https://github.com/JamieHeather/genedb-releases
