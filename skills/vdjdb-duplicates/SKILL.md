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

A matching receptor-pMHC is a candidate for review, not permission to delete a row.
Different references, donors, methods, subsets, tissues, clone IDs or other reported metadata
identify distinct observations and must be preserved. The same reference may occur in an aggregate
chunk and a paper chunk, so file boundaries alone do not establish independence.

## Step 1 - compare paired receptors

Read [completeness and observation identity](../../docs/standards/chunk-format.md#completeness-and-observation-identity).
Read chunks with `read_chunks(deduplicate=False)` so the audit can detect what deduplication
would remove. Group within species on `v.alpha`, `j.alpha`, `cdr3.alpha`, `v.beta`, `j.beta`,
`cdr3.beta`, `antigen.epitope`, `mhc.a`, `mhc.b` and `mhc.class`. Keep both chains in the same
key: sharing a beta chain while reporting different alpha chains is not a duplicate receptor.
Repeat the comparison on harmonised values and retain the submitted-value comparison.

## Step 2 - classify each matching group

Compare every `method.*` and `meta.*` field and the reference, including optional fields.
Retain differing observations. A blank and a populated field need review; do not discard the
populated value by keeping the first row. Merge only confirmed duplicate reports, preserving
complementary metadata. Curation serials, submitter and comment do not by themselves establish
independent biological observations.

Report within-file groups, cross-file groups sharing one reference, and cross-reference groups
separately. Record the decision and evidence on the paper issue. `vdjdb submission <chunk>` is a
useful recurrence summary, but a shared single-chain clonotype is not this paired-receptor audit.

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
