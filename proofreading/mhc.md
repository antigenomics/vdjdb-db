# MHC / HLA Nomenclature Reference

> **Authority file:** `proofreading/mhc_alleles.tsv.gz`
> **Sources:**
> - [IPD-IMGT/HLA Database](https://www.ebi.ac.uk/ipd/imgt/hla/) via [ANHIG/IMGTHLA](https://github.com/ANHIG/IMGTHLA)
> - WHO Nomenclature Committee report: *Nomenclature for Factors of the HLA System, 2026* (Marsh et al., HLA, 2026; 107:e70595; doi:10.1111/tan.70595)
> **IMGTHLA release used:** 3.64.0 (2026-04-16); `mhc_alleles.tsv.gz` contains **46,005 allele entries**
> **Update cadence:** IPD-IMGT/HLA releases quarterly (January, April, July, October). Re-extract `mhc_alleles.tsv.gz` after each major release.

---

## 1. Overview of IPD-IMGT/HLA

The IPD-IMGT/HLA database is the authoritative repository for human HLA sequences and nomenclature, maintained by the Anthony Nolan Research Institute (UCL) and the European Bioinformatics Institute (EMBL-EBI). As of Release 3.64.0 (April 2026), the database contains **46,005 named HLA alleles** across all loci.

The WHO Nomenclature Committee for Factors of the HLA System oversees official naming. New alleles require submission to the IPD-IMGT/HLA database for official designation.

---

## 2. HLA Allele Naming Structure

### Format
```
HLA-<GENE>*<FIELD1>:<FIELD2>[:<FIELD3>[:<FIELD4>]]
```

### Fields
| Field | Meaning | Resolution level | Example |
|---|---|---|---|
| `HLA-<GENE>` | Gene name with HLA prefix | Gene | `HLA-A` |
| `*<FIELD1>` | Serological antigen group (2+ digits) | Serological | `*02` |
| `:<FIELD2>` | Protein-level differences in coding region (2+ digits) | Protein | `:01` |
| `:<FIELD3>` | Synonymous coding changes (optional) | Genomic synonymous | `:01` |
| `:<FIELD4>` | Non-coding region differences (optional) | Full genomic | `:02` |

### Examples
| Allele name | Resolution | Notes |
|---|---|---|
| `HLA-A*02` | Low (antigen group only) | Acceptable in VDJdb; note higher resolution preferred |
| `HLA-A*02:01` | Standard (protein level) | **Preferred for VDJdb** |
| `HLA-A*02:01:01` | High (genomic synonymous) | Accepted |
| `HLA-A*02:01:01:01` | Highest (full genomic + non-coding) | Accepted; store full string |
| `HLA-DRB1*03:01` | Standard class II | Preferred |

### Special allele suffixes
| Suffix | Meaning |
|---|---|
| `N` | Null allele — not expressed (e.g., `HLA-A*01:01:01:02N`) |
| `L` | Low surface expression |
| `S` | Secreted form only |
| `C` | Expressed in cytoplasm only |
| `A` | Aberrant expression |
| `Q` | Questionable expression |

---

## 3. HLA Genes — Complete List (from HLA2026.pdf, Table 1)

### Class I genes (α-chain; paired with B2M for surface expression)
| Gene | Previous names | Molecular characteristics |
|---|---|---|
| `HLA-A` | — | Class I α-chain (classical) |
| `HLA-B` | — | Class I α-chain (classical) |
| `HLA-C` | — | Class I α-chain (classical) |
| `HLA-E` | E, '6.2' | Associated with class I, 6.2-kB HindIII fragment |
| `HLA-F` | F, '5.4' | Associated with class I, 5.4-kB HindIII fragment |
| `HLA-G` | G, '6.0' | Associated with class I, 6.0-kB HindIII fragment |
| `HLA-H` | H, AR, '12.4', HLA-54 | Class I pseudogene |
| `HLA-J` | cda12, HLA-59 | Class I pseudogene |
| `HLA-K` | HLA-70 | Class I pseudogene |
| `HLA-L` | HLA-92 | Class I pseudogene |
| `HLA-N` | HLA-30 | Class I gene fragment |
| `HLA-P` | HLA-90 | Class I gene fragment |
| `HLA-R` | HLA-OLI | Class I pseudogene (newly named 2026) |
| `HLA-S` | HLA-17 | Class I gene fragment |
| `HLA-T` | HLA-16 | Class I gene fragment |
| `HLA-U` | HLA-21 | Class I gene fragment |
| `HLA-V` | HLA-75 | Class I gene fragment |
| `HLA-W` | HLA-80 | Class I gene fragment |
| `HLA-X` | HLA-X | Class I gene fragment |
| `HLA-Y` | HLA-BEL/COQ/DEL | Class I gene fragment |
| `HLA-Z` | HLA-Z1 | Class I gene fragment (in Class II region) |

### Class II genes
| Gene | Previous names | Molecular characteristics |
|---|---|---|
| `HLA-DRA` | DRα | DR α chain (essentially monomorphic) |
| `HLA-DRB1` | DRβI, DR1B | DR β1 chain (most polymorphic DR gene) |
| `HLA-DRB2` | DRβII | Pseudogene with DR β-like sequences |
| `HLA-DRB3` | DRβIII, DR3B | DR β3 chain (DR52, Dw24/25/26) |
| `HLA-DRB4` | DRβIV, DR4B | DR β4 chain (DR53) |
| `HLA-DRB5` | DRβIII | DR β5 chain (DR51) |
| `HLA-DRB6` | DRBX, DRBσ | DRB pseudogene (DR1, DR2, DR10 haplotypes) |
| `HLA-DRB7` | DRBψ1 | DRB pseudogene (DR4, DR7, DR9 haplotypes) |
| `HLA-DRB8` | DRBψ2 | DRB pseudogene (DR4, DR7, DR9 haplotypes) |
| `HLA-DRB9` | M4.2 βexon | DRB pseudogene, isolated fragment |
| `HLA-DQA1` | DQα1, DQ1A | DQ α chain |
| `HLA-DQB1` | DQβ1, DQ1B | DQ β chain |
| `HLA-DQA2` | DXα, DQ2A | DQ α-chain-related (not known to be expressed) |
| `HLA-DQB2` | DXβ, DQ2B | DQ β-chain-related (not known to be expressed) |
| `HLA-DQB3` | DVβ, DQB3 | DQ β-chain-related (not known to be expressed) |
| `HLA-DOA` | DNA, DZα, DOα | DO α chain |
| `HLA-DOB` | DOβ | DO β chain |
| `HLA-DMA` | RING6 | DM α chain |
| `HLA-DMB` | RING7 | DM β chain |
| `HLA-DPA1` | DPα1, DP1A | DP α chain |
| `HLA-DPB1` | DPβ1, DP1B | DP β chain |
| `HLA-DPA2` | DPα2, DP2A | DP α-chain-related pseudogene |
| `HLA-DPA3` | DPA3 | DP α-chain-related pseudogene |
| `HLA-DPB2` | DPβ2, DP2B | DP β-chain-related pseudogene |

### Other HLA-region genes
| Gene | Previous names | Molecular characteristics |
|---|---|---|
| `TAP1` | ABCB2, RING4, Y3, PSF1 | ABC transporter |
| `TAP2` | ABCB3, RING11, Y1, PSF2 | ABC transporter |
| `PSMB9` | LMP2, RING12 | Proteasome-related |
| `PSMB8` | LMP7, RING10 | Proteasome-related |

---

## 4. MHC Class I vs Class II: VDJdb Conventions

### Class I entries
```
mhc.a   = HLA-A*02:01          ← classical α-chain allele
mhc.b   = B2M                   ← always literal "B2M" for class I
mhc.class = MHCI
```

### Class II entries
```
mhc.a   = HLA-DQA1*01:02       ← α-chain allele
mhc.b   = HLA-DQB1*03:02       ← β-chain allele
mhc.class = MHCII
```

### Class II DR entries (special case: DRA is monomorphic)
```
mhc.a   = HLA-DRA*01:01        ← or just HLA-DRA (non-polymorphic)
mhc.b   = HLA-DRB1*03:01       ← the polymorphic β-chain
mhc.class = MHCII
```

---

## 5. B2M (Beta-2 Microglobulin)

B2M is the invariant light chain of all MHC class I molecules. It is **not an HLA gene** — it is encoded on chromosome 15 (not in the MHC locus on chromosome 6). In VDJdb, `mhc.b` is always the literal string `B2M` for any MHC class I entry. Never substitute `beta-2-microglobulin`, `β2m`, `B2M*`, `β2-microglobulin`, or any allele-qualified version.

---

## 6. Cross-Checking mhc.class Against mhc.a

| If `mhc.a` starts with... | Then `mhc.class` should be... | And `mhc.b` should be... |
|---|---|---|
| `HLA-A`, `HLA-B`, `HLA-C` | `MHCI` | `B2M` |
| `HLA-E`, `HLA-F`, `HLA-G` | `MHCI` | `B2M` |
| `HLA-DR`, `HLA-DQ`, `HLA-DP` | `MHCII` | The corresponding β-chain allele (NOT `B2M`) |
| `HLA-DO`, `HLA-DM` | `MHCII` | The corresponding β-chain allele |

---

## 7. Non-Human MHC Naming

### Mouse (MusMusculus) — H-2 system
| Class | Example names | Notes |
|---|---|---|
| Class I | `H-2Db`, `H-2Kb`, `H-2Ld`, `H-2Dd` | Format: `H-2<locus><haplotype>` |
| Class II | `I-Ab`, `I-Ad`, `I-Ak`, `I-As`, `I-Eb` | Format: `I-<locus><haplotype>` |

**Common normalizations needed:**
- `H2-Db` → `H-2Db` (add hyphen)
- `IAb` → `I-Ab` (add hyphen after I)
- `H-2D^b` → `H-2Db`

For mouse class I: `mhc.b = B2M`, `mhc.class = MHCI`
For mouse class II: `mhc.b` = β-chain allele (e.g., `I-Ab`), `mhc.class = MHCII`

### Macaque (MacacaMulatta) — Mamu system
| Class | Example names |
|---|---|
| Class I | `Mamu-A*01`, `Mamu-A*02`, `Mamu-B*01`, `Mamu-B*17` |
| Class II | `Mamu-DRA*01`, `Mamu-DRB1*03`, `Mamu-DPB1*01` |

### Rat (RattusNorvegicus) — RT1 system
| Class | Example names |
|---|---|
| Class I | `RT1-A`, `RT1-CE16`, `RT1-Aw2` |
| Class II | `RT1-B`, `RT1-D`, `RT1-Da`, `RT1-Db1` |

---

## 8. HLA Allele Database Statistics (HLA2026.pdf, December 2025)

| Gene | Allele count | Notes |
|---|---|---|
| HLA-A | 9,022 | Most studied class I locus |
| HLA-B | 10,876 | Largest number of alleles |
| HLA-C | 9,031 | |
| HLA-DRB1 | 3,924 | Most polymorphic class II locus |
| HLA-DQB1 | 3,022 | |
| HLA-DPB1 | 3,075 | |
| HLA-DQA1 | ~500 | Less polymorphic α-chain |
| HLA-DPA1 | ~100 | Nearly monomorphic α-chain |

**Total as of IPD-IMGT/HLA 3.63.0 (January 2026): 43,758 alleles**
**Total as of IPD-IMGT/HLA 3.64.0 (April 2026): 46,005 alleles in `mhc_alleles.tsv.gz`**

---

## 9. Non-Conventional and Historical HLA Naming

Papers (especially pre-2000) may use non-standard HLA nomenclature. These must be normalised to the canonical `HLA-<GENE>*<FIELD1>:<FIELD2>` format.

### 9.1 Old serological antigen names

Before molecular typing, alleles were named by the antigen they defined. These appear without the `HLA-` prefix, without `*`, and without `:`.

| Serological name | Modern canonical name | Notes |
|---|---|---|
| `A2`, `HLA-A2` | `HLA-A*02` (low res) or `HLA-A*02:01` (if context known) | Most common HLA-A allele globally |
| `A3` | `HLA-A*03:01` | |
| `B7` | `HLA-B*07:02` | |
| `B27` | `HLA-B*27:05` (most common) | Many B27 subtypes exist — do NOT assume 05 |
| `DR4` | `HLA-DRB1*04:XX` | Very ambiguous — multiple DRB1*04 alleles |
| `DQ2` | Usually `HLA-DQB1*02:01` or `*02:02` | Ambiguous without typing |
| `B44` | `HLA-B*44:02` or `*44:03` | Two common subtypes |

> **Rule:** When a serological antigen name is given without molecular typing, record at low resolution (e.g., `HLA-A*02`) and note in the extraction log that higher resolution was not available.

### 9.2 Old molecular format (no colon separator, no field 2)

Older publications (mid-1990s to mid-2000s) sometimes use:
- `A*0201` → normalise to `HLA-A*02:01` (insert `HLA-` prefix and colon after field 1)
- `DRB1*0401` → `HLA-DRB1*04:01`
- `B*0702` → `HLA-B*07:02`

Pattern: `<GENE>*<4-digit>` → `HLA-<GENE>*<2-digit>:<2-digit>`

### 9.3 WHO serology notation (no `*`)

Some databases use `A0201`, `B0702` — no `*`, no `HLA-`:
- `A0201` → `HLA-A*02:01`
- `B0702` → `HLA-B*07:02`
- `DRB10401` → `HLA-DRB1*04:01`

### 9.4 High-resolution (3- and 4-field) alleles

Modern NGS typing produces up to 4-field resolution. VDJdb accepts and stores the full string:
| Format | Example | Resolution level |
|---|---|---|
| 1-field | `HLA-A*02` | Antigen group only — acceptable, note in log |
| 2-field | `HLA-A*02:01` | Protein level — **preferred** |
| 3-field | `HLA-A*02:01:01` | Synonymous coding changes |
| 4-field | `HLA-A*02:01:01:01` | Non-coding differences — accepted, store full string |

When looking up 3- or 4-field alleles in `mhc_alleles.tsv.gz`, match on the full string in `allele_name`. When looking up low-resolution (1-field), use prefix matching.

### 9.5 Null, low-expression, and aberrant alleles

Alleles with expression suffixes:
| Suffix | Meaning | Action |
|---|---|---|
| `N` | Null — not expressed (e.g., `HLA-A*01:01:01:02N`) | Include suffix; note in log |
| `L` | Low surface expression | Include suffix; note |
| `S`, `C`, `A`, `Q` | Secreted / cytoplasm / aberrant / questionable | Include suffix; note |

### 9.6 Non-HLA naming in context of human class I

Some older literature refers to non-classical class I molecules without full allele typing:
- `HLA-E` (no allele) → leave as-is; note single known common allele is `HLA-E*01:01` or `HLA-E*01:03`
- `MICA`, `MICB` → these are MHC-class-I-related chains, NOT HLA; should appear in `mhc.a` for NKT or innate studies; flag for review

---

## 10. How to Use `mhc_alleles.tsv.gz` for Validation

### Column reference
```
allele_id  allele_name  confirmed  sequence_type  partial  cell_count  group_count
```
- `allele_id`: IPD-IMGT/HLA internal ID (format `HLAnnnnn`)
- `allele_name`: Full canonical name with `HLA-` prefix (e.g., `HLA-A*02:01`)
- `confirmed`: `Confirmed` or `Unconfirmed`
- `sequence_type`: `gDNA` or `cDNA`
- `partial`: `Full` or `Partial`

> **Important:** All entries in `mhc_alleles.tsv.gz` are stored at full (4-field) resolution (e.g., `HLA-A*02:01:01:01`). A VDJdb entry like `HLA-A*02:01` (2-field) will **not** exact-match any row — use **prefix matching** (`$2 ~ /^HLA-A\*02:01:/`) to validate it. There are 504 rows starting with `HLA-A*02:01:` confirming this protein group exists.

### Validation queries (decompress on-the-fly with gzip -dc)

**Look up a specific 2-field allele (prefix match — all entries in file are at 4-field resolution):**
```bash
gzip -dc proofreading/mhc_alleles.tsv.gz | awk -F'\t' '$2 ~ /^HLA-A\*02:01:/'
```

**Look up an exact 4-field allele:**
```bash
gzip -dc proofreading/mhc_alleles.tsv.gz | awk -F'\t' '$2 == "HLA-A*02:01:01:01"'
```

**Check if an allele is confirmed (use prefix match for 2-field):**
```bash
gzip -dc proofreading/mhc_alleles.tsv.gz | awk -F'\t' '$2 ~ /^HLA-A\*02:01:/ {print $2, $3}' | head -3
```

**Find all alleles at a locus (prefix match):**
```bash
gzip -dc proofreading/mhc_alleles.tsv.gz | awk -F'\t' '$2 ~ /^HLA-DRB1/ {print $2, $3}' | head -20
```

**Low-resolution lookup (1-field, e.g., HLA-A*02):**
```bash
gzip -dc proofreading/mhc_alleles.tsv.gz | awk -F'\t' '$2 ~ /^HLA-A\*02:/' | head -5
# If any rows are returned, the antigen group exists
```

### Validation rules for VDJdb
1. For human entries where `mhc.a` starts with `HLA-`: look up in `mhc_alleles.tsv.gz`
   - **4-field allele** (e.g., `HLA-A*02:01:01:01`): exact match on `allele_name`
   - **2- or 3-field allele** (e.g., `HLA-A*02:01`): use **prefix match** (`$2 ~ /^HLA-A\*02:01:/`) — file stores only 4-field entries
   - **Low-resolution** (1-field, e.g., `HLA-A*02`): prefix match acceptable (`$2 ~ /^HLA-A\*02:/`); flag and note in log
   - **Not found**: check for common naming variants (sections 9.1–9.3); flag if truly absent
   - **Unconfirmed**: note in log
2. For non-human entries: use species-specific naming rules above (not validated against this file)
3. `mhc.b = B2M` is never in the allele list (it is not an HLA allele) — this is correct

---

## 11. Common MHC-II Gene Naming Errors (Data Entry)

These errors have been observed in VDJdb chunks and must be corrected during proofreading.

### 11.1 Missing digit in DP/DQ gene names

Several HLA class II genes include a digit suffix (e.g., `DPA1`, `DPB1`, `DQA1`, `DQB1`). It is a common data-entry mistake to omit this digit:

| Wrong form | Correct form | Gene | Notes |
|---|---|---|---|
| `HLA-DPA*01:03` | `HLA-DPA1*01:03` | DPA1 | DP α chain — missing `1` |
| `HLA-DPB*04:01` | `HLA-DPB1*04:01` | DPB1 | DP β chain — missing `1` |
| `HLA-DQA*01` | `HLA-DQA1*01` | DQA1 | DQ α chain — missing `1` |
| `HLA-DQB*06:02` | `HLA-DQB1*06:02` | DQB1 | DQ β chain — missing `1` |

**Note:** `HLA-DRA` is correct — DRA has no digit suffix (the gene is simply `DRA`, not `DRA1`). Do NOT write `HLA-DRA1*`.

### 11.2 `HLA-DRA1*` — incorrect extra digit

`HLA-DRA1*` is wrong. The gene name is `DRA` (no trailing digit). Always use `HLA-DRA*01:01` (or `HLA-DRA*01` for low-resolution).

| Wrong form | Correct form |
|---|---|
| `HLA-DRA1*01` | `HLA-DRA*01` |
| `HLA-DRA1*01:01` | `HLA-DRA*01:01` |

### 11.3 Missing `HLA-` prefix in mhc.a / mhc.b

Human MHC-II alleles in `mhc.a` and `mhc.b` must always include the `HLA-` prefix:

| Wrong form | Correct form |
|---|---|
| `DPA1*02:02` | `HLA-DPA1*02:02` |
| `DPB1*05:01` | `HLA-DPB1*05:01` |
| `DRB1*15:01` | `HLA-DRB1*15:01` |
| `DQB1*06:02` | `HLA-DQB1*06:02` |

**Quick scan for all three error types:**

```bash
# Missing digit in DP/DQ gene names
awk -F'\t' 'NR>1 && ($10 ~ /^HLA-DP[AB]\*/ || $11 ~ /^HLA-DP[AB]\*/ || $10 ~ /^HLA-DQ[AB]\*/ || $11 ~ /^HLA-DQ[AB]\*/) {print NR, $10, $11}' <chunk_file>

# HLA-DRA1* (extra digit)
awk -F'\t' 'NR>1 && ($10 ~ /^HLA-DRA1\*/ || $11 ~ /^HLA-DRA1\*/) {print NR, $10, $11}' <chunk_file>

# Missing HLA- prefix (starts with DP/DQ/DR but not HLA-)
awk -F'\t' 'NR>1 && ($10 ~ /^D[PQR][ABMNO]/ || $11 ~ /^D[PQR][ABMNO]/) {print NR, $10, $11}' <chunk_file>
```

---

## 12. Reference Information

- **WHO Nomenclature Report 2026:** Marsh et al., *HLA* 2026; 107:e70595. doi:10.1111/tan.70595
- **IPD-IMGT/HLA Database:** https://www.ebi.ac.uk/ipd/imgt/hla/
- **IMGTHLA GitHub:** https://github.com/ANHIG/IMGTHLA
- **HLA nomenclature website:** https://hla.alleles.org
- **IMGT-MHC (non-human):** https://www.ebi.ac.uk/ipd/mhc/
