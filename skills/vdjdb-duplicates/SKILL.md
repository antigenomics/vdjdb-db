---
name: vdjdb-duplicates
description: Audit repeated TCR records across the whole VDJdb corpus - which clonotypes and pMHC pairs recur, whether a recurrence is a same-lab follow-up publication or independent replication (by PubMed author overlap), which within-chunk multiplicities are sequencing read depth rather than distinct clones, and which epitopes are recorded against inconsistent or genuinely multiple MHC restrictions. Use for a corpus-wide consistency audit or before landing a large chunk that may already be in the database.
---

# vdjdb-duplicates

Measure how records repeat across the corpus, and separate the three reasons they do.

Read [`skills/AUTHORITIES.md`](../AUTHORITIES.md) first.

## Invocation

```
/vdjdb-duplicates
```

No arguments. Needs a build: `uv run vdjdb build --out out/`.

## What a repeat is, and is not

**Two rows in two different chunks are independent reports, never duplicates**, even when every
field matches. A chunk is one publication, so a matching row in a second publication is a second
laboratory finding the same receptor against the same peptide. That is signal: it is what raises
`vdjdb.score` and what the motif clustering is tuned against. Deduplication in this database is
**within a chunk only**, on `schema.CHUNK_DEDUP_KEY`, and the build already does it - the count it
removes is the gap between 203,308 raw rows and 192,753 released ones, reported as the advisory
`duplicate` rule.

So this audit is not looking for errors to delete. It is separating three things that look alike in a
row count:

| Reason a record repeats | What it means | Action |
|---|---|---|
| Same lab, follow-up publication | expected; the same cohort re-sequenced or re-analysed | keep both; note the relationship |
| Independent replication | a public clonotype, the strongest evidence the database holds | keep both; this is the finding |
| Within-chunk multiplicity from read depth | one clone counted once per cell or per read | keep; state it in the release notes so nobody reads the count as clonal abundance |

## Step 1 - load the built tables

The build assigns the ids this audit needs, so do not rebuild the keys by hand:

- `clonotype_id` on `chains.tsv` keys `(species, gene, cdr3, v.segm, j.segm)`
- `pmhc_id` on `records.tsv` keys `(antigen.epitope, mhc.a, mhc.b)`
- `epitope_id` keys the epitope alone

```python
import polars as pl
ch = pl.read_csv('out/tables/chains.tsv', separator='\t', infer_schema_length=0)
rec = pl.read_csv('out/tables/records.tsv', separator='\t', infer_schema_length=0)
d = ch.join(rec.select('record_id', 'pmhc_id', 'epitope_id', 'reference.id', 'chunk.file',
                       'meta.subject.id', 'meta.clone.id'), on='record_id', how='left')
```

A receptor against a peptide is `(clonotype_id, pmhc_id)`. Group on that pair, not on a
hand-assembled tuple of `cdr3.beta` and `antigen.epitope` - the pair follows the harmonised call and
the repaired sequence, which is what the database actually ships.

## Step 2 - classify each recurring pair

For every `(clonotype_id, pmhc_id)` with more than one record:

| Class | Test |
|---|---|
| within one chunk | one distinct `chunk.file` |
| across chunks, one reference | several files, one `reference.id` |
| across chunks, same lab | several PMIDs with overlapping author lists (step 3) |
| across chunks, independent | several PMIDs with no shared authors |

Report the counts per class and the largest groups. `vdjdb submission <chunk>` gives the same
relationship for one chunk against the corpus, without a build of your own, and is the right tool when
the question is about one submission rather than the whole database.

## Step 3 - author overlap

For pairs of PMIDs that share more than about 20 recurring groups, fetch both author lists and
intersect them:

```python
import urllib.request, json, time

def authors(pmid):
    url = ('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
           f'?db=pubmed&id={pmid}&retmode=json')
    with urllib.request.urlopen(url, timeout=10) as r:
        return [a['name'] for a in json.load(r)['result'][pmid].get('authors', [])]
```

Rate-limit to under three requests a second, or NCBI will throttle without saying so.

Three or more shared authors reads as one group publishing twice. Zero to two reads as independent
replication. State the count, not the verdict alone - the threshold is a convention, and a large
consortium paper breaks it.

`summary/reference_years.tsv` already holds every reference's year and is offline, so use it for
chronology rather than a second round of requests. Publication order is what distinguishes a
follow-up from a re-analysis.

## Step 4 - within-chunk multiplicity

A clonotype appearing many times inside one chunk, from one or two donors, is usually the assay's
resolution rather than a repertoire measurement. Two patterns, and they need opposite readings:

**Read depth.** Few donors, and the *same* clonotype repeated - one clone observed once per cell in
single-cell sequencing, or once per read. The count is instrument output, not clonal abundance.

**Deep repertoire.** Few donors, and *many distinct* clonotypes against one epitope. That is what a
bulk repertoire study of an immunodominant epitope looks like, and it is the data doing its job.

Distinguish on distinct `clonotype_id` per group, and on whether `meta.clone.id` varies. Neither is a
defect and neither is removed. The first is stated in the release notes; the second is left alone.

Where the chunk's own `method.frequency` or `meta.subset.frequency` already records abundance, say so:
a row count standing in for a frequency that the chunk reports properly is a submission that could be
collapsed, and that is a question for the submitter.

## Step 5 - one epitope, several MHC restrictions

Group `records.tsv` by `epitope_id` and list the distinct `mhc.a`. Then split the result three ways,
because they need three different responses:

**Different HLA genes.** Genuine cross-restriction exists and is published - an epitope presented by
both an A and a B allele. Report it; do not correct it. Where the corpus carries one such pairing on
very few records against many on another, that minority is worth checking against its paper.

**Different alleles of one gene.** Normal. Populations differ.

**The same allele at two resolutions** - `HLA-A*02` beside `HLA-A*02:01`. That is a spelling problem,
not biology: the two do not join, so a query on either misses the other. This is the case to fix, and
it is fixed by finding the resolution the paper reports, per chunk.

For "which other alleles could present this epitope", do not derive it here. That is
`uv run vdjdb promiscuity`, which answers it against `mhcmatch` predictions and writes
`proofreading/epitope_promiscuity.tsv`, marking which pairings the database already records. It makes
no claim about any record's `mhc.a`.

## Step 6 - spellings that differ only in case or a separator

`out/reports/lookalikes.tsv`, written by the build, is this check over every value column. Sort on
`same.species`: `true` means two spellings sit on one organism, so one is wrong. `false` can be
correct - HGNC capitalises human gene symbols and MGI title-cases mouse ones.

## Step 7 - report

```
=== VDJDB RECURRENCE AUDIT ===
records: N        chunks: N        distinct clonotypes: N        distinct pMHC: N

recurring (clonotype, pMHC) pairs: N
  within one chunk:            N groups
  across chunks, one reference: N groups
  across chunks, same lab:      N groups   (>=3 shared authors)
  across chunks, independent:   N groups

largest cross-publication pairs:
  N groups   PMID:X x PMID:Y   S shared authors   <relationship>

within-chunk multiplicity, >=50 records from <=3 donors:
  N groups   <read depth | deep repertoire>   <chunk>   <epitope>

epitopes with more than one MHC gene:        N      (cross-restriction, keep)
epitopes with one allele at two resolutions: N      (spelling, fix per chunk)

look-alike values with same.species=true:    N

RECOMMENDATIONS:
  - <allele resolution to settle, per chunk and paper>
  - <read-depth chunks to name in the release notes>
```

Findings that are corpus-wide and mechanical belong on the tracker, not in this file: a numbered issue
that a chunk branch can close, with the row count it fires on today. Findings that need a paper read
belong on that paper's PMID issue.
