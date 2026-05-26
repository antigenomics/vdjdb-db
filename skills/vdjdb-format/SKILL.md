---
name: vdjdb-format
description: Standardise a raw or partially-formatted VDJdb TSV chunk — normalising species names, IMGT V/D/J gene IDs, IMGT-HLA MHC alleles, and method vocabulary — and produce a properly-named chunk file ready for proofreading.
---

# /format — VDJdb Chunk Formatting Skill

## Purpose

Take the output of `/extract` (or any TSV resembling a chunk file) and standardise all controlled-vocabulary fields to the VDJdb / IMGT specification. Cross-check naming conventions against existing `chunks/` files to ensure internal consistency. This is the second stage: extract → **format** → proofread.

## Invocation

```
/format [path-to-tsv]
```

The TSV must have a VDJdb-compatible header (see canonical column order in `/extract`). If the file has structural problems (missing columns, wrong separator), halt and report — this is a job for `/proofread` Step 1, not format.

---

## Standardisation Rules

### 1. Species Names

Normalise to the exact VDJdb-accepted values (case-sensitive):

| Normalise FROM | Normalise TO |
|---|---|
| `Homo sapiens`, `human`, `H. sapiens`, `hs`, `Human` | `HomoSapiens` |
| `Mus musculus`, `mouse`, `M. musculus`, `mm`, `Mouse` | `MusMusculus` |
| `Rattus norvegicus`, `rat`, `Rattus` | `RattusNorvegicus` |
| `Macaca mulatta`, `rhesus`, `macaque`, `NHP` | `MacacaMulatta` |

If a species is not in the above list:
1. Log it as a candidate for extending `speciesList` in `py_src/ChunkQC.py`
2. Ask the user whether to include or exclude those rows

---

### 2. IMGT V/D/J Gene IDs

**Primary authority:** `proofreading/imgt_alleles.tsv` (column `imgt_gene_id`)
**Conversion table:** `patches/nomenclature.conversions`
**Secondary fallback:** `patches/IGM_nomenclature_table.tsv`

**Rules (apply in order):**

1. **Strip whitespace**: remove all spaces within the gene name (`TRBV 7` → `TRBV7`, `TRAV 12-2` → `TRAV12-2`)

2. **Look up in `imgt_alleles.tsv`** (strip allele suffix `*NN` before lookup):
   - If found → keep (or correct capitalisation to match)
   - If not found → check `patches/nomenclature.conversions` for a mapping
   - If found in conversions → apply the conversion and log `old_name → new_name`
   - If not found in either → flag as unresolvable; ask user

3. **Validate allele number** (if present, e.g., `TRBV12-3*02`):
   - Look up `allele_count` for that gene in `proofreading/imgt_alleles.tsv`
   - If allele number > `allele_count`: flag as invalid allele; suggest dropping the allele suffix or contacting IMGT

4. **Check functionality**: if `functionality` is `P` (pseudogene) in `imgt_alleles.tsv`: flag as biologically suspicious

5. **Consistency check against existing chunks**:
   ```bash
   grep -h "" chunks/*.txt | cut -f3 | sort -u | grep "^TRAV"  # check v.alpha values
   ```
   If the same gene appears with different notation in existing chunks (e.g., `TRAV13-1` vs `TRAV13`), standardise to the IMGT-canonical form.

6. **Multiple gene possibilities** (comma-separated ambiguous assignments): check each against `imgt_alleles.tsv`, keep all valid candidates comma-separated without spaces (e.g., `TRBV7-2,TRBV7-3`)

---

### 3. MHC Alleles

**Primary authority:** `proofreading/mhc_alleles.tsv` (column `allele_name`)
**Reference:** `proofreading/mhc.md`

#### Human MHC (HLA)

Target format: `HLA-<GENE>*<FIELD1>:<FIELD2>` (e.g., `HLA-A*02:01`)

| Problem | Fix |
|---|---|
| `A02`, `A0201` (old serological) | → `HLA-A*02:01` if unambiguous; flag if ambiguous |
| `HLA-A*02` (low resolution, no field 2) | Keep as-is; note in log that higher resolution preferred |
| `https://doi.org/...` prefix on allele | Remove prefix |
| `HLA-A 02:01` (space) | → `HLA-A*02:01` |
| Confirmed in `mhc_alleles.tsv` | Log confirmation status from `confirmed` column |

**MHC-I second chain**: always normalise to literal `B2M` — never `beta-2-microglobulin`, `β2m`, `b2m`, `B2M*01`, etc.

**mhc.class cross-check:**
| If `mhc.a` starts with... | `mhc.class` must be | `mhc.b` must be |
|---|---|---|
| `HLA-A`, `HLA-B`, `HLA-C`, `HLA-E`, `HLA-F`, `HLA-G` | `MHCI` | `B2M` |
| `HLA-DR`, `HLA-DQ`, `HLA-DP`, `HLA-DO` | `MHCII` | HLA β-chain allele |

#### Mouse MHC (H-2)

| Normalise FROM | Normalise TO |
|---|---|
| `H2-Db`, `H2Db` | `H-2Db` |
| `IAb`, `I-Ab`, `IA-b` | `I-Ab` |
| `H-2D^b` | `H-2Db` |

For mouse class I: `mhc.b = B2M`
For mouse class II: `mhc.b` = the β-chain name (e.g., `I-Ab`)

---

### 4. Method Vocabulary

Normalise `method.identification` and `method.verification` to VDJdb-recognised terms.

**Common normalizations:**
| Author writes | Normalise to | Notes |
|---|---|---|
| `FACS sort with pMHC tetramer` | `tetramer-sort` | |
| `pMHC multimer sort` | `tetramer-sort` | |
| `dextramer sort` | `dextramer-sort` | |
| `pentamer sort` | `pentamer-sort` or `pelimer-sort` | Confirm which reagent |
| `ELISpot` | Do NOT map automatically | Log and ask user |
| `51Cr release assay` | Do NOT map automatically | Log and ask user |

**Rule:** If an identification or verification method has no close equivalent in the current VDJdb vocabulary, do NOT force it. Instead:
1. Leave a descriptive string in the field (for reference)
2. Document it under "Vocabulary gaps" in the format log
3. Suggest adding it as a new term via a note to the database maintainers

---

### 5. Reference IDs

Enforce correct format:

| Problem | Fix |
|---|---|
| `https://doi.org/10.1016/...` | → `doi:10.1016/...` (remove URL prefix) |
| `http://dx.doi.org/10.1016/...` | → `doi:10.1016/...` |
| `doi: 10.1016/...` (space after colon) | → `doi:10.1016/...` |
| `pubmed:12345678` or `PubMed:12345678` | → `PMID:12345678` |
| Bare number `12345678` | Ask if it is a PMID; if confirmed → `PMID:12345678` |

---

### 6. Antigen Fields

Cross-reference `patches/antigen_epitope_species_gene.dict`:
- If `antigen.epitope` exists in the dict → use the dict's `antigen.gene` and `antigen.species` (this ensures consistency with the full database)
- If the epitope is new → keep author-provided gene/species values, note in log

---

### 7. Chunk ID

After all formatting changes, renumber `chunk.id` sequentially from 1 (integer, no leading zeros).

---

## Output Filename

Prefer `PMID_<pubmed_id>.txt` (e.g., `PMID_28975614.txt`).

If no PMID is available:
1. Ask the user for the preferred name
2. Alternatives: `doi_<mangled_doi>.txt`, submitter-date format
3. Check that the chosen name does not duplicate an existing file in `chunks/`

Suggested placement: `chunks_unformatted/` if uncertain about QC status; `chunks/` only after `/proofread` passes.

---

## Format Log

Write `<output_basename>_format_log.txt` containing:

1. **Changes made**: for each change — field name, old value, new value, source of normalisation (imgt_alleles.tsv / mhc_alleles.tsv / nomenclature.conversions / manual)
2. **Unresolvable fields**: fields that could not be normalised and why
3. **Vocabulary gaps**: novel method/verification terms encountered
4. **Allele resolution notes**: alleles that exist in `mhc_alleles.tsv` at low resolution only
5. **Consistency discrepancies**: naming differences found vs existing `chunks/` files
6. **Pseudogene warnings**: gene names whose `functionality = P` in `imgt_alleles.tsv`
7. **Unconfirmed HLA alleles**: alleles present in `mhc_alleles.tsv` with `confirmed = Unconfirmed`

---

## Reference Files

| File | Role |
|---|---|
| `proofreading/imgt_alleles.tsv` | **Primary** IMGT V/D/J gene authority |
| `proofreading/imgt.md` | IMGT nomenclature rules |
| `proofreading/mhc_alleles.tsv` | **Primary** HLA allele authority |
| `proofreading/mhc.md` | MHC/HLA naming rules, class I vs II, non-human |
| `patches/nomenclature.conversions` | Old → current IMGT gene name mappings |
| `patches/IGM_nomenclature_table.tsv` | Secondary IMGT fallback (existing repo file) |
| `patches/antigen_epitope_species_gene.dict` | Known epitope → gene/species mappings |
| `py_src/ScoreFactory.py` | Method vocabulary and scoring logic |
| `py_src/ChunkQC.py` | ALL_COLS definition (canonical column list) |
| `chunks/*.txt` | Reference for consistency checks |
| `README.md` | Full VDJdb specification |

---

## Next Step

After formatting, run `/proofread [path-to-formatted-tsv]` to validate against `py_src/ChunkQC.py` and all other QC checks.
