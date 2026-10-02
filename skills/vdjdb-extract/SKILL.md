---
name: vdjdb-extract
description: Extract TCR:pMHC specificity records from raw submission sources - supplementary XLS/CSV tables, PDF manuscripts, 10x Genomics contig and clonotype files, AIRR Rearrangement TSVs, Adaptive ImmunoSEQ exports - into a VDJdb chunk TSV, verifying every extracted sequence, gene, allele and reference id back against the source. Use when a paper, dataset or submitter folder has to become a chunk file, before vdjdb-format and vdjdb-proofread.
---

# vdjdb-extract

Turn raw source files into one VDJdb chunk TSV plus an extraction log. First of three stages:
**extract** → [format](../vdjdb-format/SKILL.md) → [proofread](../vdjdb-proofread/SKILL.md).

Read [`skills/AUTHORITIES.md`](../AUTHORITIES.md) first. Invariant 5 there governs this skill: no
value reaches the output that was not found in the source.

## Invocation

```
/vdjdb-extract [path-to-folder-or-file]
```

Respect any scope the user sets ("beta chains only", "skip the MHC typing") and record the limit in
the log.

An existing reference is not a completeness check. Reconcile every supplied sequence table against
its publication's observations, including negative tables and supplementary validation experiments.
Record an explicit disposition for every remainder. Apply a curator-approved historical cutoff
only to the chunks it covers; never replace this review with a filter for newly seen publications.

## Step 1 - inventory the source

List every file, name its type, and say what each is likely to hold: TCR sequences, epitope and MHC
data, methods, references. Show the inventory to the user before extracting if there are more than
three files or the structure is unclear.

## Step 2 - build the join graph explicitly

Submission data is almost always spread over several files. Identify every id column in each file
(barcode, clonotype id, clone id, well, sample, donor), determine which ids appear in more than one
file, perform the join, and log the keys and the cardinality.

Log every ambiguity **as you find it**: one-to-many links, a missing join key, two files disagreeing
on the same cell. Never pick one option silently.

## Step 3 - extract the complex fields

One row per record, reporting **both chains of one clone in that row** - `cdr3.alpha` and `cdr3.beta`
are columns of the same row, not two rows. At least one of the two must be non-blank.

`cdr3.alpha` `v.alpha` `j.alpha` `cdr3.beta` `v.beta` `d.beta` `j.beta` `species` `mhc.a` `mhc.b`
`mhc.class` `antigen.epitope` `antigen.gene` `antigen.species` `reference.id`

For every sequence, gene name, allele, species and reference id: search the original file for the
exact string, and log found / not found / found as a variant. A value that cannot be confirmed is
`[UNVERIFIED]` and is withheld pending user approval.

Leave a field blank when the source does not report it. Never a placeholder (invariant 2).

Cross-reference each epitope against `patches/antigen_epitope_species_gene.dict`; where the epitope
is already known, take the dict's gene and species so the chunk agrees with the rest of the database.

### Sequence admission

A TCR junction **starts with `C`, ends with `F` or `W`, and carries no other cysteine**. That is the
definition of the region, not a statistical tendency. Import is the cheapest place to catch a breach of
it - before the file has a name, an issue or a commit - and it is where the submitter still has their
own source open.

**Nothing is dropped for failing the definition.** The record is kept, flagged, and warned about. A
non-canonical junction can be the best record of what a publication reported, and people filter on the
flag: the build ships `v.canonical`, `j.canonical` and `cdr3.one.cysteine` on `chains`, plus the same
three asked of the sequence as submitted, and `vdjdb qc` reports `internal cysteine in cdr3.alpha` /
`.beta` on the raw chunk.

| Case | Action |
|---|---|
| Any character outside the 20 canonical amino acids (`X`, `B`, `*`, `#`) | drop the row, log it - this is an invalid alphabet, not a non-canonical junction |
| Fewer than 4 residues | drop the row, log it |
| Does not start with `C`, or does not end with `F` or `W` | **keep it, and work out which of the four causes below it is.** Fix it at the source where you can, flag it where you cannot, and log either way |
| Carries a cysteine after the first residue | **keep it and warn.** Rare rather than impossible - the Jurkat receptor has one, and a disulphide-bonded loop is a real thing to report - so check it against the source and say so in the log |
| Genuinely modified or non-natural residues | ask the user, then `chunks_with_unconventional_aa/` |

Four causes of a missing anchor, in the order they are worth checking:

1. **The source is in IMGT CDR3 space.** Both anchors are absent because that coordinate system
   excludes them. The column is `cdr3_aa`, not `junction_aa` - invariant 3. The whole file is affected,
   not one row, so check one sequence and you have checked all of them. This is the one to fix rather
   than flag: re-read the right column.
2. **Framework was left in.** `YLCSSQEGGYGYTFGSG` carries residues in front of the Cys and behind the
   Phe. Trim to the anchors.
3. **The anchor was mis-read.** `GASSDTMNTKIL` for `CASSDTMNTKIL`. One residue, and the V germline says
   which it should be.
4. **The V or J call is wrong, and the sequence is right.** Mouse `TRAJ47` resolves to an ORF `*01`
   whose anchor is genuinely lost, and every corpus chain on it reads the functional `*02` signature.
   Repair the call, not the sequence - `proofreading/cdr3_repair.md`.

Cause 1 is a whole-file mistake and worth stopping for. Causes 2 to 4 are per-row: fix what the source
supports and flag the rest. Log the count per cause.

Epitopes are canonical amino acids only. A chemical modification, a non-peptide antigen or a peptide
pool is flagged and asked about before it is included.

### `reference.id`

| Form | Example |
|---|---|
| `PMID:<digits>` - preferred | `PMID:28975614` |
| `doi:<doi>` - lowercase prefix, no URL | `doi:10.1016/j.immuni.2023.01.001` |
| preprint URL | `https://www.biorxiv.org/content/10.1101/2024.01.01.123456` |
| PDB entry with no publication | `https://www.rcsb.org/structure/1AO7` |
| no publication at all | `unpublished: <Submitter Name> <YYYY-MM-DD>` |

Never infer a PMID. If it is not in the source, leave the field blank and ask.

## Step 4 - extract the method fields

`method.identification` `method.frequency` `method.singlecell` `method.sequencing`
`method.verification`

These determine `vdjdb.score`. Take the vocabulary from
[`docs/standards/chunk-format.md`](../../docs/standards/chunk-format.md) and the score rules from
[`docs/standards/confidence-score.md`](../../docs/standards/confidence-score.md). Two rules:

- **Record what the source says.** If the authors write only "multimer" or only "sorted", the term is
  `multimer-sort`. Do not promote it to `tetramer-sort` because tetramers are commoner in VDJdb.
- **If no term fits, do not force one.** Leave the author's wording in the field, log it under
  "vocabulary gaps", and propose the new term to the maintainers.

`method.frequency` is a count over a total (`7/30`) or a fraction. A percentage is not that form - see
[format](../vdjdb-format/SKILL.md) for where a group-level percentage belongs instead.

## Step 5 - extract the meta fields

`meta.study.id` `meta.cell.subset` `meta.subset.frequency` `meta.subject.cohort` `meta.subject.id`
`meta.replica.id` `meta.clone.id` `meta.epitope.id` `meta.tissue` `meta.donor.MHC`
`meta.donor.MHC.method` `meta.structure.id`

Fill as many as the source supports; leave the rest blank. `meta.*` and `method.*` describe the
record, not the act of curating it. Only `submitter`, `comment` and `chunk.id` are curation
properties. `comment` is under 140 characters.

A `meta.structure.id` is a four-character PDB entry id. Anything else there - a figure or table
reference - is a curation decision, not a fill: `vdjdb qc` reports it as
`structure id is not a PDB id` and does not fail.

## Step 6 - write the TSV

**Take the header from a shipping chunk rather than any list, including the ones above:**

```bash
head -1 chunks/PMID_28423320.tsv
```

That is the 33-column canonical header, `chunk.id` first. 173 of 230 chunks carry exactly it.

- tab-separated, UTF-8, LF line endings, no quoting, no trailing whitespace
- `chunk.id`: integers from 1
- blank means empty, per invariant 2
- add `comment` as a 34th column only if at least one row has one

## Step 7 - write the extraction log

`<basename>_extraction_log.txt`, with: the source inventory and the role assigned to each file; the
join graph, keys and cardinality; every ambiguity and what was decided or escalated; the verification
result per field; any `[UNVERIFIED]` value the user approved; vocabulary gaps; dropped rows and why;
and the scope limits the user set.

## Output and next step

`<PMID_xxxxxxx>_unformatted.txt` and `<PMID_xxxxxxx>_extraction_log.txt`, named from the PMID where
there is one. Gene names are not yet normalised - that is the next stage.

Then run `/vdjdb-format` on the output.

## Source-specific traps

**Excel** (always load with `data_only=True`, or formula cells arrive as their formula text):

1. **Repeated header rows mid-table** marking a new donor or group. They show as rows with `TCRα`,
   `TCRβ`, `TRAV`, `TRBV`, `CDR3α`, `CDR3β` as literal cell values. Filter on the gene columns as
   well as the CDR3 columns - some have a blank CDR3 and the marker only in a gene column.
2. **Allele plus functionality code** in one cell: `TRAV16*01 F`. `\*\d+\s*$` misses it because `F`
   follows the space. Strip from the `*`: `re.sub(r'\*.*$', '', v).strip()`.
3. **Formula artefacts**: `TRAJ3+D107:D1082` is a gene name plus a cell reference. Take
   `val.split('+')[0].strip()`.
4. **J and D columns swapped** relative to their own headers. Decide by the gene prefix, never by the
   column: `TRBJ2-7*01` is a J call wherever it sits.
5. **Copy-paste characters** in CDR3 cells (`#`, `X`, `*`). Drop those rows and log them.

**Adaptive Biotech ImmunoSEQ**: gene columns use a `TCRB`/`TCRA` prefix and zero-padded numbers
(`TCRBV06-05*01`). Do not convert them here - flag them all for `/vdjdb-format`, which applies the
rules in `proofreading/imgt.md` §9.2.

**AIRR Rearrangement TSV**: `junction_aa` is VDJdb's `cdr3`; `cdr3_aa` is **two residues shorter** and
is not (invariant 3). Confirm which column you are reading before writing a single row.

**10x Genomics**: `filtered_contig_annotations.csv` pairs chains by `barcode`; the clonotype file
aggregates them. One clone becomes one row with both chains, not two rows.
