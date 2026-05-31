---
name: vdjdb-duplicates
description: Identify and classify duplicate TCR records across all VDJdb chunks at three resolution levels (beta-only, paired, and same-epitope multi-MHC), categorise by publication source, author overlap, and flag spurious high-frequency records.
---

# /vdjdb-duplicates — VDJdb Duplicate & Consistency Audit Skill

## Purpose

Scan all `chunks/` files to:
1. Find duplicate TCRs at two CDR3 resolution levels.
2. Classify duplicates as within-publication, cross-publication same-lab, or genuinely independent.
3. Check author overlap between publications sharing many TCR sequences.
4. Flag extremely frequent records that suggest read-count inflation rather than unique T cell clones.
5. Report epitopes presented by multiple distinct MHC molecules or with inconsistent allele resolution.

## Invocation

```
/vdjdb-duplicates
```

No arguments. Run from the repo root.

---

## Step 1 — Load All Chunks

```python
import csv, glob, re
from collections import defaultdict, Counter

rows_all = []
for path in sorted(glob.glob('chunks/*.txt')):
    fname = path.split('/')[-1]
    with open(path) as f:
        for row in csv.DictReader(f, delimiter='\t'):
            row['_file'] = fname
            rows_all.append(row)

def v(row, col): return (row.get(col) or '').strip()
```

---

## Step 2 — Beta-Only Duplicates

Key: `(cdr3.beta, v.beta, antigen.epitope)` — records sharing the same beta chain and epitope regardless of alpha chain or donor metadata.

```python
key_beta = defaultdict(list)
for row in rows_all:
    cb = v(row,'cdr3.beta'); vb = v(row,'v.beta'); ep = v(row,'antigen.epitope')
    if cb and ep:
        key_beta[(cb, vb, ep)].append(row)

dups_beta = {k: vs for k, vs in key_beta.items() if len(vs) > 1}
```

**Classify each duplicate group:**

| Category | Criterion |
|---|---|
| Within-file | All rows in one chunk file |
| Cross-file, same reference | Multiple files, all share the same `reference.id` |
| Cross-file, same lab | Multiple PMIDs but overlapping authors (check Step 4) |
| Cross-file, independent | Multiple PMIDs with no shared authors |

Count and report each category. List the top 30 groups by frequency.

---

## Step 3 — Paired Duplicates

Key: `(cdr3.alpha, v.alpha, j.alpha, cdr3.beta, v.beta, j.beta, antigen.epitope)` — exact paired chain duplicates. Only applied to rows where both CDR3 chains are non-empty.

```python
key_pair = defaultdict(list)
for row in rows_all:
    ca = v(row,'cdr3.alpha'); cb = v(row,'cdr3.beta'); ep = v(row,'antigen.epitope')
    if ca and cb and ep:
        k = (ca, v(row,'v.alpha'), v(row,'j.alpha'),
             cb, v(row,'v.beta'),  v(row,'j.beta'), ep)
        key_pair[k].append(row)

dups_pair = {k: vs for k, vs in key_pair.items() if len(vs) > 1}
```

Report same categories as Step 2.

---

## Step 4 — Author Overlap Check

For each cross-file duplicate group involving ≥2 distinct PMIDs, retrieve author lists and compute overlap:

```python
import urllib.request, json, time

def get_authors(pmid):
    url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={pmid}&retmode=json'
    with urllib.request.urlopen(url, timeout=10) as r:
        data = json.load(r)
    return [a['name'] for a in data['result'][pmid].get('authors', [])]

# Build PMID-pair → shared author count
pmid_pairs = Counter()
for k, vs in dups_beta.items():
    refs = sorted({v(r,'reference.id') for r in vs if v(r,'reference.id').startswith('PMID:')})
    for i in range(len(refs)):
        for j in range(i+1, len(refs)):
            pmid_pairs[(refs[i], refs[j])] += 1

# For top pairs (> threshold shared duplicates), check authors
OVERLAP_THRESHOLD = 20
for (r1, r2), cnt in pmid_pairs.most_common():
    if cnt < OVERLAP_THRESHOLD: break
    p1, p2 = r1.replace('PMID:',''), r2.replace('PMID:','')
    try:
        a1 = get_authors(p1); time.sleep(0.4)
        a2 = get_authors(p2); time.sleep(0.4)
        shared = set(a1) & set(a2)
        print(f"{r1} × {r2}: {cnt} shared groups, {len(shared)} common authors")
        if shared:
            print(f"  Authors: {', '.join(sorted(shared)[:6])}")
    except Exception as e:
        print(f"  ERROR: {e}")
```

**Classification rule:**
- ≥3 shared authors → "same lab, follow-up publication" (expected overlap, not spurious)
- 0–2 shared authors → "independent replication" (genuine public clonotype)

---

## Step 5 — High-Frequency Record Detection

Within-file groups with frequency ≥50 from ≤3 distinct donors are suspicious. Distinguish two legitimate assay causes before flagging as data errors:

**Pattern A — Single TCR tested against many epitopes (combinatorial assay):**
The same CDR3 appears with many different epitopes in one study. Each row is a genuine specificity claim. High n within one epitope group from few donors still indicates read-depth, not clonal abundance.

**Pattern B — Pool of TCRs tested against one epitope/pattern:**
Many different CDR3s all assigned the same epitope from a single donor. Here n reflects different clones in the repertoire, not read inflation. This is biologically expected in large repertoire studies.

Distinguishing them: if few donors but many distinct CDR3s → Pattern B (normal). If few donors with the SAME CDR3 repeated n times → Pattern A / read inflation.

```python
FREQ_THRESHOLD = 50
DONOR_THRESHOLD = 3

for k, vs in sorted(dups_beta.items(), key=lambda x: -len(x[1])):
    files = {r['_file'] for r in vs}
    if len(files) > 1: continue   # only single-file
    donors = {v(r,'meta.subject.id') for r in vs}
    clones = {v(r,'meta.clone.id') for r in vs if v(r,'meta.clone.id')}
    if len(vs) >= FREQ_THRESHOLD and len(donors) <= DONOR_THRESHOLD:
        cb, vb, ep = k
        refs = {v(r,'reference.id') for r in vs}
        pattern = 'READ-INFLATION' if len(clones) <= 1 else 'DEEP-REPERTOIRE'
        print(f"{pattern} n={len(vs):4d} CDR3b={cb:22} ep={ep:15} donors={len(donors)} ref={list(refs)[0]}")
```

**Known confirmed cases from vdjdb-db audit (2026-05-31):**

| File | Epitope | CDR3b | n | Donors | Pattern | Cause |
|---|---|---|---|---|---|---|
| PMID_41315082 | VEALYLVCG | CASSEAGTGGYEQYF | 530 | 2 | A | scTCR-seq read depth; same clone seen many times across cells |
| PMID_39746936 | VISNDVCAQV | multiple | 160–271 | 1 | A | Single donor deep-seq; each clone's frequency encoded as row count |
| PMID_34811538 | RAKFKQLL/CLGGLLTMV | multiple | 107–241 | 2–3 | B | Bulk repertoire depth; many clones, legitimate |
| 10xgenomics-2019-07-09 | IVTDFSVIK/RAKFKQLL | multiple | 100–133 | 2 | B | 10x Genomics multiplexed assay; cell barcodes give multiplicity |

**Interpretation:** Pattern A records represent sequencing read depth not unique T cells — flag in release notes but do not remove. Pattern B records are biologically valid and expected.

---

## Step 6 — Multi-MHC Epitope Report

```python
key_mhc = defaultdict(set)
for row in rows_all:
    ep = v(row,'antigen.epitope'); mhca = v(row,'mhc.a')
    if ep and mhca: key_mhc[ep].add(mhca)

# Flag epitopes with distinct HLA genes (not just allele sub-typing)
for ep, alleles in sorted(key_mhc.items(), key=lambda x: -len(x[1])):
    genes = set()
    for a in alleles:
        m = re.match(r'(HLA-[A-Z0-9]+|H2-\w+)', a)
        if m: genes.add(m.group(1))
    if len(genes) > 1:
        print(f"{ep:20} {len(alleles)} alleles, {len(genes)} genes: {sorted(genes)}")

# Flag allele resolution inconsistencies (same gene, coarse + fine)
for ep, alleles in key_mhc.items():
    coarse = {a for a in alleles if '*' in a and ':' not in a}
    fine   = {a for a in alleles if ':' in a}
    if coarse and fine:
        coarse_g = {a.split('*')[0] for a in coarse}
        fine_g   = {a.split('*')[0] for a in fine}
        if coarse_g & fine_g:
            print(f"RESOLUTION MIX {ep}: coarse={sorted(coarse)[:2]} fine={sorted(fine)[:2]}")
```

---

## Step 7 — Summary Report

```
=== VDJDB DUPLICATE AUDIT SUMMARY ===
Total rows: N
Total chunks: N

BETA-ONLY DUPLICATES
  Unique duplicate groups: N
  Total redundant rows:    N
  Within-file:             N groups
  Cross-file:              N groups
    Same-lab (≥3 shared authors): N
    Independent replication:       N

TOP CROSS-PUBLICATION PAIRS (by shared CDR3b+Vb+epitope groups):
  [N]  PMIDX × PMIDY — [shared authors count] common authors — [lab relationship]
  ...

HIGH-FREQUENCY SUSPECTS (≥50 copies, ≤3 donors, single file):
  [N]  CDR3b / epitope — file — likely cause

MULTI-MHC EPITOPES (distinct HLA genes):
  [N epitopes] — list top cases

MHC FORMAT ISSUES (allele resolution inconsistency):
  [N epitopes] — coarse vs fine allele mix

RECOMMENDATIONS:
  - Flag high-frequency within-file duplicates in the DB release notes
  - Resolve coarse/fine allele mix per paper (check original publication)
  - Same-lab cross-publication duplicates: expected, keep all records
  - Independent cross-publication public clonotypes: expected, keep all records
```

---

## Known Findings from 2026-05-31 Audit

### Cross-publication duplicate patterns

| PMID pair | Shared groups | Shared authors | Relationship |
|---|---|---|---|
| PMID:37749325 × PMID:40694338 | 921 | 17 (Kedzierska K et al.) | Same lab, follow-up study |
| PMID:28423320 × PMID:37749325 | 171 | 0 | Public clonotypes (GIL, CMV) |
| PMID:23267020 × PMID:24512815 | 126 | 3 (Koning D, van Baarle D) | Same lab, method comparison |
| PMID:35589842 × PMID:40713946 | 95 | 0 | Neoantigen public clonotypes |
| PMID:28250417 × PMID:28629751 | 38 | 4 (Selin LK et al.) | Same lab, influenza repertoire |
| PMID:18802118 × PMID:21562156 | 35 | 7 (Kalams SA et al.) | Same lab, HIV studies |
| PMID:19017975 × PMID:21135165 | 31 | 7 (Price DA, Douek DC) | Same lab, longitudinal HIV/CMV |

**Interpretation:** The dominant source of cross-publication duplicates is same-lab follow-up publications. True independent replication (0 shared authors) reflects genuinely public clonotypes for immunodominant epitopes (GIL, GLCTLVAML, NLVPMVATV).

### Most replicated epitopes (cross-publication)

| Epitope | Cross-pub duplicate rows | Notes |
|---|---|---|
| GILGFVFTL | 7,436 | Influenza GIL — the most public CD8 epitope in humans |
| GLCTLVAML | 1,032 | EBV BMLF1 — highly public, 10+ studies |
| FLRGRAYGL | 546 | EBV EBNA3A — restricted by HLA-B*08 AND HLA-A*02:01 (genuine bi-restriction) |
| NLVPMVATV | 377 | CMV pp65 — dominant CMV epitope |
| YLQPRTFLL | 298 | SARS-CoV-2 Spike — multiple COVID-19 cohort studies |

### Spurious high-frequency records

Extremely frequent within-file records (≥100 copies, ≤3 donors) reflect sequencing read depth, not individual T cell clones. Known cases:

- **PMID_41315082**: CASSEAGTGGYEQYF/VEALYLVCG n=530, 2 donors — high-throughput scTCR-seq
- **PMID_39746936**: 5 CDR3b/VISNDVCAQV combinations n=160–271, 1 donor
- **PMID_34811538**: Multiple CDR3b n=107–241, 2–3 donors — EBV/beta-cell antigen study
- **10xgenomics-2019-07-09**: n=100–133, 2 donors — multiplexed 10x Genomics data

These are **not data errors** but should be noted in release documentation.

### Genuine multi-MHC restriction

| Epitope | MHC genes | Notes |
|---|---|---|
| FLRGRAYGL | HLA-A + HLA-B | Published cross-restriction A*02 and B*08:01 |
| RAKFKQLL (EBV BZLF1) | HLA-A + HLA-B | Atypical; primary restriction is B*08; A*02 entries warrant review |
| RPPIFIRRL | HLA-A + HLA-B | Cross-restriction A*02 and B*07 — published |

### MHC notation fixes applied (2026-05-31)

- `H-2Db`/`H-2Kb`/`H-2Kd`/`H-2Ld`/`H-2Dd` → `H2-Db`/`H2-Kb`/`H2-Kd`/`H2-Ld`/`H2-Dd` (3,001 rows: strip hyphen between H and 2)
- `H-2KB` → `H2-Kb` (capitalization fix)
- `H2 class I` + SIINFEKL → `H2-Kb` (OVA/C57BL6 context)
- `H2-b class I` + SIINFEKL → `H2-Kb`
- `H2-b class II` + SIINFEKL → **removed** (MHC-I epitope mislabeled as class II)
- `H2 class II` + Ins2/SHLVEALYLVCGERG → `H2-IAg7` (NOD mouse T1D context)
- `H2-d class II` + SFERFEIFPKE → `H2-IEd` (BALB/c HA restriction)
