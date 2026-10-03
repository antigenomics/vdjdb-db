---
name: vdjdb-harmonize
description: Canonicalise antigen.gene and antigen.species in a VDJdb chunk against the epitope dictionary and the gene and species alias tables, detect spurious values (UniProt descriptions, "[species]" annotations, "Probable"/"Chain A," prefixes, multi-word organism names), resolve blanks from the publication, flag one epitope that is a substring of another, and extend the alias tables where they have no entry. Use when antigen.gene or antigen.species carry free text rather than symbols, or when vdjdb-proofread reports antigen fields that disagree across a chunk.
---

# vdjdb-harmonize

Bring `antigen.gene` and `antigen.species` to the names VDJdb records, and extend the tables where
they have no answer. Runs standalone or from [proofread](../vdjdb-proofread/SKILL.md) step 5.

Read [`skills/AUTHORITIES.md`](../AUTHORITIES.md) first.

## Invocation

```
/vdjdb-harmonize [path-to-tsv]
```

## Why this is a skill and not a build stage

`patches/antigen_epitope_species_gene.dict` is applied by the build (`vdjdb.curate.patch`) and covers
515 epitopes, which is 166,153 of the 203,348 chunk rows - 81.7 %, measured 2026-09-29. Nothing in
the build reads the two alias tables, and that is deliberate: mapping free text onto a gene symbol is
a curation
decision, made once, recorded in a table, and applied to the chunk before it lands. The build applies
declared patches; it does not guess at prose.

So the deliverable here is two things, and the second matters more: **the chunk, and the table rows
that made the chunk's values derivable.** A value normalised without a table row behind it comes back
on the next submission of the same data.

## The three sources, in priority order

| Priority | Source | Keyed on | Rows |
|---|---|---|---:|
| 1 | `patches/antigen_epitope_species_gene.dict` | the epitope, exactly | 515 |
| 2 | `proofreading/gene_aliases.tsv` | free-text gene name, exact after stripping | 202 |
| 3 | `proofreading/species_aliases.tsv` | lowercase substring of the organism, **first match wins** | 84 |

An epitope in source 1 settles both fields; skip 2 and 3 for that row. The species table is ordered,
so a more specific fragment must sit above a less specific one (`human herpesvirus 4` before
`herpesvirus`) - when you add a row, add it in the right place, not at the end.

## Step 1 - load and apply

```python
import csv, re

def _table(path, header):
    """Two-column TSV, comments skipped, header row skipped by its own first field."""
    with open(path) as f:
        rows = [r for r in csv.reader(f, delimiter='\t')
                if r and not r[0].startswith('#') and r[0] != header]
    return [(r[0].strip(), r[1].strip()) for r in rows]

epitopes = {}                     # epitope -> (species, gene)
with open('patches/antigen_epitope_species_gene.dict') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        epitopes[r['antigen.epitope'].strip()] = (r['antigen.species'].strip(),
                                                  r['antigen.gene'].strip())
genes = dict(_table('proofreading/gene_aliases.tsv', 'source_name'))
species = _table('proofreading/species_aliases.tsv', 'fragment')   # a list: the order is the rule

SUFFIXES = (' protein', ' glycoprotein', ' polyprotein', ' precursor')

def gene(raw):
    s = re.sub(r'\s*\[[^\]]+\]\s*$', '', raw).strip()     # drop a trailing "[CMV]"
    if s in genes:
        return genes[s]
    for suf in SUFFIXES:
        if s.lower().endswith(suf) and len(s) > len(suf):
            stripped = s[:-len(suf)].strip()
            return genes.get(stripped, stripped)
    return s

def organism(raw):
    low = raw.lower()
    return next((c for frag, c in species if frag in low), raw)

def harmonise(row):
    ep = row.get('antigen.epitope', '').strip()
    if ep in epitopes:
        row['antigen.species'], row['antigen.gene'] = epitopes[ep]
    else:
        row['antigen.gene'] = gene(row.get('antigen.gene', '').strip())
        row['antigen.species'] = organism(row.get('antigen.species', '').strip())
    return row
```

## Step 2 - what counts as spurious

These are the patterns that say a value is a description rather than a name. They are what triggers
this skill from `/vdjdb-proofread`, and what to fix here.

**`antigen.gene`:**

| Pattern | Example | Fix |
|---|---|---|
| a `[species]` annotation | `pp65 [CMV]` | strip, re-look-up |
| a trailing ` protein` / ` glycoprotein` / ` polyprotein` / ` precursor` | `Spike glycoprotein` | strip, re-look-up |
| a `Probable ` / `Putative ` / `Chain A, ` / `MULTISPECIES:` prefix | `Chain A, Nucleoprotein` | strip, re-look-up |
| a full UniProt description, or a comma | `Sterol-4-alpha-carboxylate 3-dehydrogenase, decarboxylating` | alias table; add the row |
| exactly `Polyprotein` | | resolve by epitope - a polyprotein epitope has a specific gene |
| blank, and `antigen.species` is not `Synthetic` | | step 4 |

**`antigen.species`:**

| Pattern | Example | Fix |
|---|---|---|
| a multi-word scientific name | `Human herpesvirus 4` | fragment table → `EBV` |
| a parenthetical common name | `Columba livia (carrier pigeon)` | strip, CamelCase → `ColumbaLivia` |
| casing | `Epstein barr virus`, `influenzaA`, `synthetic` | normalise; `Synthetic` is capitalised |
| blank | | step 4 |

**Conventions:** two-word binomials become CamelCase with no space (`BacillusSubtilis`); well-known
pathogens keep their established abbreviation (`EBV`, `CMV`, `HIV-1`, `SARS-CoV-2`, `InfluenzaA`); a
genus with unknown species stays one CamelCase word (`Bacillus`); human and mouse genes take the HGNC
or MGI symbol (`PABPC1`, `G6pc2`); viral genes take the short name the literature uses (`pp65`,
`BMLF1`, `Gag`, `Tax`).

Do not hold a hardcoded list of accepted species in this file - it drifts the moment a chunk lands.
Derive it:

```bash
cut -f15 chunks/*.tsv | sort | uniq -c | sort -rn
```

A value absent from that output is new to VDJdb, which `vdjdb submission` also reports, and new is not
the same as wrong.

## Step 3 - consistency, within the chunk and against the corpus

Two checks, both advisory, both needing a curator:

**One epitope with two gene or species values.** After step 1, group the chunk's rows by
`antigen.epitope` and report any epitope with more than one `antigen.gene` or `antigen.species`. What
remains after the dict has been applied is a gap in the dict: resolve it and add the row.

**Spellings that differ only in case or a separator.** `IE1` and `IE-1` are two `antigen.gene` values
for one CMV gene, and no query filtering on one finds the other. The build computes this over the
whole corpus:

```bash
uv run vdjdb build --out out/        # writes out/reports/lookalikes.tsv
```

Sort that report on `same.species` - `true` means the two spellings sit on one `antigen.species`, so
one of them is wrong. `false` can be correct: HGNC capitalises human symbols and MGI title-cases
mouse ones, so `MBP` and `Mbp` are two conventions for one gene and both stay.

## Step 4 - blank antigen fields

A blank `antigen.gene` is advisory. A mimotope or designed peptide may have no source gene;
other papers may not identify it. Keep unreported provenance blank when the paper establishes
a junction and epitope. Review a blank `antigen.species` by the same rule.

Resolve in this order:

1. **The epitope dictionary** - it may already carry the epitope.
2. **The rest of the corpus** - `grep -h '<EPITOPE>' chunks/*.tsv | cut -f13,14,15 | sort -u`. Another
   paper's curated answer for the same peptide is the strongest available prior.
3. **The publication.** Fetch the abstract for the row's `reference.id` and read the antigen context.
   A neoantigen study is `HomoSapiens`; a cross-reactivity study varies per epitope.
4. **Genuinely unknown, and curated as such** - for structural entries with no reported antigen, write
   `Unknown` rather than leaving blank, so a checked-and-unknown value is distinguishable from an
   unchecked one.

Confirm species assignments that a related pathogen could explain. `SARS-CoV` and `SARS-CoV-2` share
epitopes, and the abstract is what settles which the paper studied.

## Step 5 - epitope substrings

Report any epitope in the chunk that is an exact substring of a longer epitope, in the chunk or in the
corpus: `EPLPQGQLTAY` inside `GPEPLPQGQLTAY`. Skip anything under 4 residues.

This is a warning and never a fix. A shorter peptide can be the genuine minimal epitope, or it can be
an export that lost its flanks. The paper decides. State both candidates and ask.

## Step 6 - report

State: rows processed; every changed cell as `ROW <chunk.id>: <field> <old> -> <new>` with the source
that decided it; consistency findings the dict does not cover; substring warnings; and **the table
rows added**, by file. Then re-run `vdjdb qc` on the result.

If the chunk is already in `chunks/`, changing it is a chunk edit: its own branch, its own issue, and
a message naming the files, the row counts and the authority. See invariant 1.
