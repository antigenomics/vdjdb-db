---
name: vdjdb-harmonize
description: Harmonize antigen.gene and antigen.species fields in a VDJdb chunk to canonical VDJdb naming. Detects spurious gene/species names, resolves inconsistencies (same epitope → multiple names), and warns about epitopes that are exact substrings of longer epitopes. Invoked standalone or from /proofread when spurious values are detected.
---

# /harmonize — VDJdb Antigen Gene/Species Harmonization Skill

## Purpose

Normalize `antigen.gene` and `antigen.species` values in a chunk to the VDJdb canonical vocabulary. This skill uses a unified, layered lookup system drawing from two patch sources:

- **`patches/antigen_epitope_species_gene.dict`** — epitope-keyed authority: maps known epitopes directly to their canonical species + gene. This is the highest-confidence source.
- **`proofreading/gene_aliases.tsv`** — free-text gene name aliases → VDJdb gene name.
- **`proofreading/species_aliases.tsv`** — species substring fragments → VDJdb species name. Checked in order; first match wins.

## Invocation

```
/harmonize [path-to-tsv]
```

Or called automatically from `/proofread` (Step 6a) when spurious antigen fields are detected.

---

## Step 1 — Load Lookup Tables

Load all three sources at the start of the session:

```python
import csv, re

# 1a. Epitope → (species, gene) from patches/
epitope_dict = {}  # {epitope: (species, gene)}
with open('patches/antigen_epitope_species_gene.dict') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        ep = row['antigen.epitope'].strip()
        epitope_dict[ep] = (row['antigen.species'].strip(), row['antigen.gene'].strip())

# 1b. Gene alias table
gene_map = {}  # {raw_name: canonical}
with open('proofreading/gene_aliases.tsv') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'): continue
        parts = line.split('\t')
        if len(parts) >= 2 and parts[0] not in ('source_name', 'vdjdb_name'):
            gene_map[parts[0].strip()] = parts[1].strip()

# 1c. Species fragment list (order-sensitive)
species_fragments = []  # [(fragment_lower, canonical), ...]
with open('proofreading/species_aliases.tsv') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'): continue
        parts = line.split('\t')
        if len(parts) >= 2 and parts[0] != 'fragment':
            species_fragments.append((parts[0].strip().lower(), parts[1].strip()))
```

---

## Step 2 — Harmonization Functions

Use this unified pipeline for every row:

```python
GENE_STRIP_SUFFIXES = (' protein', ' glycoprotein', ' polyprotein', ' precursor')

def harmonize_gene(raw: str) -> str:
    """Map raw antigen.gene to VDJdb canonical name."""
    if not raw:
        return raw
    # Strip embedded [species] annotation (e.g. "pp65 [CMV]")
    s = re.sub(r'\s*\[[^\]]+\]\s*$', '', raw).strip()
    # Exact match
    if s in gene_map:
        return gene_map[s]
    # Strip common trailing suffixes and retry
    for suffix in GENE_STRIP_SUFFIXES:
        if s.lower().endswith(suffix) and len(s) > len(suffix):
            stripped = s[:-len(suffix)].strip()
            if stripped in gene_map:
                return gene_map[stripped]
            if len(stripped) >= 2:
                return stripped   # cleaned fallback
    return s

def harmonize_species(raw: str) -> str:
    """Map raw antigen.species to VDJdb canonical name via substring fragments."""
    if not raw:
        return raw
    low = raw.lower()
    for fragment, canonical in species_fragments:
        if fragment in low:
            return canonical
    return raw  # unchanged if no match

def harmonize_row(row: dict) -> dict:
    """Apply full harmonization pipeline to a single row."""
    ep = row.get('antigen.epitope', '').strip()

    # Priority 1: epitope dict (most authoritative)
    if ep in epitope_dict:
        canon_species, canon_gene = epitope_dict[ep]
        row['antigen.species'] = canon_species
        row['antigen.gene'] = canon_gene
        return row

    # Priority 2: alias tables
    raw_gene = row.get('antigen.gene', '').strip()
    raw_species = row.get('antigen.species', '').strip()
    new_gene = harmonize_gene(raw_gene)
    new_species = harmonize_species(raw_species)

    if new_gene != raw_gene:
        row['antigen.gene'] = new_gene
    if new_species != raw_species:
        row['antigen.species'] = new_species

    return row
```

---

## Step 3 — Detect Spurious Values

Before harmonizing, scan the chunk and flag any rows where `antigen.gene` or `antigen.species` looks suspicious. These are the patterns that trigger automatic `/harmonize` invocation from `/proofread`:

### 3a. Spurious `antigen.gene` indicators

| Pattern | Example | Action |
|---|---|---|
| Contains `[...]` species annotation | `pp65 [CMV]` | Strip annotation, re-lookup |
| Ends with ` protein`, ` glycoprotein`, ` polyprotein`, ` precursor` | `Spike glycoprotein` | Strip suffix, canonical lookup |
| Length > 30 characters | `RNA-directed RNA polymerase catalytic subunit` | Alias lookup; flag if no match |
| Contains spaces AND is not a known multi-word canonical name | `Nucleoprotein M1` | Flag for manual review |
| Matches `Polyprotein` or `polyprotein` exactly | | Flag — epitopes from polyprotein should have specific gene assigned; look up by epitope in dict |
| Null / empty | | Flag `no.antigen.gene` |

```python
KNOWN_MULTIWORD_GENES = {'PB1-F2', 'non-structural', 'Large T antigen'}  # extend as needed

def is_spurious_gene(gene: str) -> bool:
    if not gene: return True
    if re.search(r'\[.+\]', gene): return True
    if any(gene.lower().endswith(s) for s in GENE_STRIP_SUFFIXES): return True
    if len(gene) > 30 and gene not in KNOWN_MULTIWORD_GENES: return True
    if gene.lower() in ('polyprotein', 'unknown', 'na', 'n/a'): return True
    return False
```

### 3b. Spurious `antigen.species` indicators

| Pattern | Example | Action |
|---|---|---|
| Contains full scientific name with spaces | `Human herpesvirus 4` | Fragment match lookup |
| Known alias variants | `HIV`, `HTLV`, `Influenza`, `IAV` | Map to canonical |
| Capitalization mismatch | `Epstein barr virus`, `influenzaA` | Normalise |
| Not in known canonical set | anything not in the list below | Flag for review |

**Known canonical `antigen.species` values** (derive from existing chunks):
```
EBV, CMV, MCMV, HSV-1, HSV-2, VZV, InfluenzaA, InfluenzaB, SARS-CoV-2, SARS-CoV,
HCoV-OC43, HCoV-HKU1, HIV-1, HCV, HBV, YFV, DENV, DENV1, DENV2, DENV3, DENV3/4,
LCMV, HTLV-1, MLV, RSV, MCPyV, HomoSapiens, MusMusculus, RattusNorvegicus,
MacacaMulatta, GallusGallus, M.tuberculosis, PlasmodiumFalciparum, PlasmodiumBerghei,
Trypanosoma cruzi, SIV, CoxsackievirusB, Wheat, ManducaSexta
```

```python
CANONICAL_SPECIES = {
    'EBV', 'CMV', 'MCMV', 'HSV-1', 'HSV-2', 'VZV', 'InfluenzaA', 'InfluenzaB',
    'SARS-CoV-2', 'SARS-CoV', 'HCoV-OC43', 'HCoV-HKU1', 'HIV-1', 'HCV', 'HBV',
    'YFV', 'DENV', 'DENV1', 'DENV2', 'DENV3', 'DENV3/4', 'LCMV', 'HTLV-1',
    'MLV', 'RSV', 'MCPyV', 'HomoSapiens', 'MusMusculus', 'RattusNorvegicus',
    'MacacaMulatta', 'GallusGallus', 'M.tuberculosis', 'PlasmodiumFalciparum',
    'PlasmodiumBerghei', 'Trypanosoma cruzi', 'SIV', 'CoxsackievirusB',
    'Wheat', 'ManducaSexta',
}

def is_spurious_species(species: str) -> bool:
    if not species: return True
    if species in CANONICAL_SPECIES: return False
    if ' ' in species: return True   # multi-word → likely not normalised
    return True  # unknown single token → flag
```

---

## Step 4 — Cross-Chunk Consistency Check

After harmonizing individual rows, check for same-epitope inconsistencies across the chunk (and optionally across all of `chunks/`):

```python
from collections import defaultdict

def check_consistency(rows: list[dict]) -> list[str]:
    warnings = []
    ep_genes = defaultdict(set)
    ep_species = defaultdict(set)
    for r in rows:
        ep = r.get('antigen.epitope', '').strip()
        g = r.get('antigen.gene', '').strip()
        s = r.get('antigen.species', '').strip()
        if ep:
            if g: ep_genes[ep].add(g)
            if s: ep_species[ep].add(s)
    for ep, genes in ep_genes.items():
        if len(genes) > 1:
            warnings.append(f'INCONSISTENT antigen.gene for epitope {ep!r}: {genes} — check patches/antigen_epitope_species_gene.dict')
    for ep, sps in ep_species.items():
        if len(sps) > 1:
            warnings.append(f'INCONSISTENT antigen.species for epitope {ep!r}: {sps} — check patches/antigen_epitope_species_gene.dict')
    return warnings
```

If the epitope dict has a definitive entry, the inconsistency will be resolved by Step 2. Report remaining inconsistencies (those the dict does not cover) for manual curation.

---

## Step 5 — Epitope Substring Warning

Check whether any epitope in the chunk is an **exact substring** of a longer epitope in the same chunk, or in the global VDJdb database. This detects truncation artefacts (e.g., `EPLPQGQLTAY` being a substring of `GPEPLPQGQLTAY`).

```python
def check_epitope_substrings(rows: list[dict], global_epitopes: set[str] = None) -> list[str]:
    """
    Warn when epitope A is an exact substring of longer epitope B.
    Checks within the chunk and against global_epitopes if provided.
    """
    warnings = []
    chunk_epitopes = {r.get('antigen.epitope', '').strip() for r in rows if r.get('antigen.epitope')}
    all_epitopes = chunk_epitopes | (global_epitopes or set())

    for ep in sorted(chunk_epitopes):
        if len(ep) < 4: continue
        for longer in all_epitopes:
            if longer != ep and ep in longer:
                src = 'chunk' if longer in chunk_epitopes else 'global VDJdb'
                warnings.append(
                    f'SUBSTRING WARNING: {ep!r} is contained in longer epitope {longer!r} ({src}) — '
                    f'possible truncation artefact; verify correct epitope boundaries'
                )
    return warnings
```

To load global VDJdb epitopes for the cross-chunk check:

```python
import glob

def load_global_epitopes() -> set[str]:
    eps = set()
    for f in glob.glob('chunks/*.txt'):
        with open(f) as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                ep = row.get('antigen.epitope', '').strip()
                if ep: eps.add(ep)
    return eps
```

---

## Step 6 — Apply and Report

Run the full harmonization on the chunk, report all changes and warnings:

```python
def harmonize_chunk(path: str, check_global: bool = True) -> None:
    with open(path) as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    fieldnames = list(rows[0].keys()) if rows else []

    changes = []
    for i, row in enumerate(rows):
        orig_gene = row.get('antigen.gene', '')
        orig_species = row.get('antigen.species', '')
        row = harmonize_row(row)
        if row['antigen.gene'] != orig_gene:
            changes.append(f'  ROW {row["chunk.id"]}: antigen.gene {orig_gene!r} → {row["antigen.gene"]!r}')
        if row['antigen.species'] != orig_species:
            changes.append(f'  ROW {row["chunk.id"]}: antigen.species {orig_species!r} → {row["antigen.species"]!r}')
        rows[i] = row

    # Consistency check
    consistency_warnings = check_consistency(rows)

    # Substring check
    global_eps = load_global_epitopes() if check_global else None
    substring_warnings = check_epitope_substrings(rows, global_eps)

    # Write back
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

    # Report
    print(f'=== HARMONIZATION REPORT: {path} ===')
    print(f'Rows processed: {len(rows)}')
    print(f'Fields changed: {len(changes)}')
    for c in changes: print(c)
    if consistency_warnings:
        print(f'\nConsistency warnings ({len(consistency_warnings)}):')
        for w in consistency_warnings: print(f'  {w}')
    if substring_warnings:
        print(f'\nSubstring warnings ({len(substring_warnings)}):')
        for w in substring_warnings: print(f'  {w}')
    if not changes and not consistency_warnings and not substring_warnings:
        print('No issues found.')
```

---

## Integration with /proofread

When `/proofread` is running **Step 6** (MHC Consistency) or scanning `antigen.gene`/`antigen.species` fields, it should invoke harmonization automatically if any of these are detected:

- `is_spurious_gene(row['antigen.gene'])` returns True for ≥1 row
- `is_spurious_species(row['antigen.species'])` returns True for ≥1 row
- `check_consistency(rows)` returns any warnings

**From /proofread Step 6, add:**

> **Step 6a — Antigen Harmonization Trigger**
>
> After running ChunkQC, scan all `antigen.gene` and `antigen.species` values using the spurious-value detectors from `/harmonize`. If any are flagged:
> 1. Report the count and examples to the user.
> 2. Ask: "Run `/harmonize` to fix these automatically? [y/n]"
> 3. If yes: run the full harmonization pipeline, then re-run ChunkQC to verify no regressions.

---

## Patch File Maintenance

When harmonization fails for a value (no alias match, no epitope dict entry), **add the mapping** to the appropriate patch file rather than leaving it unfixed:

- New **epitope → species/gene** mapping → append to `patches/antigen_epitope_species_gene.dict`
- New **gene alias** (free-text → canonical) → append to `proofreading/gene_aliases.tsv`
- New **species fragment** (substring → canonical) → append to `proofreading/species_aliases.tsv` (put more-specific fragments before less-specific ones)

---

## Reference Files

| File | Role |
|---|---|
| `patches/antigen_epitope_species_gene.dict` | Epitope-keyed authority: epitope → (species, gene) |
| `proofreading/gene_aliases.tsv` | Free-text gene alias table → VDJdb gene name |
| `proofreading/species_aliases.tsv` | Species fragment → VDJdb canonical (order-sensitive) |
| `py_src/ChunkQC.py` | Run after harmonization to verify no regressions |
