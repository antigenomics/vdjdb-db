# Submit and review records

A chunk records observations from one publication. Start with the
[TSV template](https://raw.githubusercontent.com/antigenomics/vdjdb-db/dev/template.tsv) or
[Excel template](https://raw.githubusercontent.com/antigenomics/vdjdb-db/dev/template.xlsx).
The [chunk reference](standards/chunk-format.md) defines the columns; leave unknown values blank.
The Excel template formats cells as text to preserve allele names and fractions.

## Find the publication issue

Search open and closed issues by PMID, DOI and publication URL, including issue bodies. Check
existing chunks and `proofreading/reference_ids.tsv` before creating a duplicate submission.
A preprint and its published version may describe the same work; use a verified PubMed version
link to reconcile them while preserving existing observations and record identifiers.

The [PMID–DOI table](https://github.com/antigenomics/vdjdb-db/blob/master/proofreading/pubmed_doi.tsv)
lists verified publication identifiers. A blank DOI means the retrieved PubMed record supplied none.

New paper issues use the title `PMID:<id>` and label `vdjdb-records-paper-pending`. Their body is
the citation retrieved from PubMed, linked to `https://pubmed.ncbi.nlm.nih.gov/<id>/`. When a DOI
or URL issue resolves to a PMID, keep the earlier identifier in its body. Public provenance text
uses the paper, its abstract or PubMed record; patent and PDB records are also supported.

## Prepare the chunk

Check every supplied sequence table, including negative results, against the imported observations.
An existing publication or chunk does not establish that its later tables are complete. Record each
table's imported observations and any unresolved remainder before closing the publication issue.

Use `chunks/PMID_<id>.tsv` for a paper. Each row needs an epitope and at least one alpha or beta
junction. Neither can be invented or inferred. VDJdb stores junctions with both anchors included,
not the shorter IMGT CDR3 region. See [sequence repair](standards/cdr3-fixing.md).

Use the supplied definitive tables first. Routine boundary trimming, supported nomenclature
conversion and allele resolution happen during database assembly. Keep the submitted cells when
the build already resolves them, and record the outcome as **resolved during database build**.
Do not retrieve a paper repeatedly to settle a case the existing machinery handles. If a paper is
inaccessible, retain supported values from the supplied table and record any remaining question.

Check the assembled output: not every missing call can be recovered. A proposed V is reported as
inferred and is not silently promoted to a paper-reported V call. Class I uses `B2M` as its partner;
class II partner inference must follow the supported MHC rules, not a guessed allele pairing.

Keep experimental details in their structured columns: `method.*` describes the assay,
`meta.*` describes subjects, donors, tissues, subsets and clones. Reserve `comment` for information
that has no suitable field. Do not manufacture donor identities or assay evidence.

## Review completeness and repeated observations

```bash
uv run vdjdb qc chunks/PMID_<id>.tsv
uv run vdjdb submission chunks/PMID_<id>.tsv
```

Read the score distribution and fatal/advisory findings. Compare matches using both chains'
V/J calls and junctions, epitope and MHC, then inspect all observation metadata. Different methods can validate the same
sample within one publication. Independent validation across studies requires evidence that the
underlying samples or experiments are independent; different publication IDs alone do not establish
this. Repeated curation of the same assay is reused data. Preserve existing PDB and mixed-paper
records unless a documented correction is necessary.

Internal cysteine, an unusual terminal residue or a non-functional segment is an advisory finding,
not sufficient reason to delete a source-supported observation. Inspect assembly repair and flags.
Missing epitope, missing both junctions, an unresolved pairing or a conflicting positive/negative
outcome requires action before the affected observation can enter the positive build.

## Separate blockers from follow-up questions

| Question | Action |
|---|---|
| Missing receptor or epitope, ambiguous pairing, unresolved restriction, conflicting assay outcome | Keep affected observations pending; ask a specific source-based question |
| Unreported donor genotype, optional frequency, extra method detail that does not change identity or the supported score | Leave the field blank; keep a nonblocking follow-up issue |
| Terminal flank or allele spelling handled by assembly | Retain the source cell; note the build resolution |
| More source material requested after supplied observations are fully imported | Track the request separately from the completed import |

Noncritical metadata must not hold up a supported record when its assay and verification evidence
already determine the [confidence score](standards/confidence-score.md). Record the exact missing
field and its impact. Do not change a score or infer an experimental result to avoid a question.

## Choose the destination

- `chunks/`: supported positive observations used by the build.
- `chunks_negative/`: explicitly negative assay observations; excluded from the positive build.
  A receptor can have a positive binding observation and a negative functional-assay observation.
- `pending/`: readable material whose import is blocked. Preserve unresolved source rows and state
  which evidence or build capability would unblock them.
- `withheld/`: material that cannot yet be assessed because its format needs re-export.

Do not silently discard unresolved rows. Keep the issue open while source observations remain
pending, withheld or in an unmerged PR. An excluded original file preserved for unresolved rows is
not evidence that every row in it is still missing; list the unresolved subset explicitly.

## Publish through dev

Start a data branch from current `dev`. Normally use one chunk per branch and PR; a maintainer may
explicitly group a batch, with a per-file manifest and publication links. Keep already-open PRs
separate unless consolidation was requested.

```bash
git switch -c codex/chunk-PMID_<id> origin/dev
uv run vdjdb identity update
```

Fetch Git LFS inputs with `git lfs pull` before building. The identity registry and large compressed
negative chunks use LFS. Commit the chunk and `registry/records.tsv.gz` together. Review additions, amendments and retirements;
existing record IDs must survive. The commit names the changed files, per-file row counts, source
issue, reason, build result and curator decision where applicable.

Open the PR against `dev`. Read `chunk-check` and fix fatal findings. Put code, test expectations
and release-comparison declarations in a **separate validation PR**. Review each changed count:
corpus growth can change release projections, motif baselines and reference counts. A green test
must follow an explained measurement, not a blanket tolerance increase.

Run full CI on the combined data and validation candidate before integration. When merging separate
PRs sequentially, refresh the validation branch onto the new dev, verify its complete tree matches
the green candidate, and wait for its required check. Chunks reach `master` only through a later
`dev` release merge.

## Close out the issue

After the data and validation are merged into `dev`, link the PRs and record the imported scope.
Close the paper issue only when no source rows, excluded observations, unmerged changes or material
questions remain. GitHub may not auto-close an issue for a PR targeting `dev`; verify its state.

If only optional metadata remains, keep a focused follow-up open and mark it nonblocking. Remove
the paper-pending label once the import is complete, so it is not counted as an unprocessed paper.
Umbrella issues stay open until all their component submissions are accounted for. Deferred
10X reconciliation is tracked separately and is not part of a routine paper closeout.

## Curation skills

The repository's [curation skills](https://github.com/antigenomics/vdjdb-db/tree/master/skills) guide
extraction, formatting, proofreading, publication, harmonisation and duplicate review. Use the
[structural collection skill](https://github.com/antigenomics/vdjdb-db/tree/master/skills/vdjdb-structures)
for deposited TCR-pMHC complexes and the local tcren pipeline. They use the
same CLI and authority tables as the build. Existing authorization for a batch applies throughout;
ask for new scientific decisions, not repeated permission for already-authorized routine actions.

Publication issues use `PMID:<id>` as the title and a linked citation from the PubMed API's NBIB
record: authors, article title, journal, date and pages. Keep source attachments and one concise
open question or disposition. A preprint without a PMID uses its title and DOI; do not guess an
identifier. Closed imports link their current chunk and imported scope. A duplicate names its
active intake; an excluded or unsuitable submission states why it does not enter the build.
Quarantined submissions remain open.

Metadata identifiers must be reported in the cited publication or its supplementary tables.
Do not copy generated export identifiers, joined identifier lists or reference-derived labels into
`meta.clone.id`, donor, study or epitope fields. If the paper does not supply an identifier, leave
the field blank. Use the PMID in `reference.id`; do not encode it as a clone identifier.

## Check experimental provenance through repertoire overlap

Run `vdjdb overlap --by-sample` for the complete report, including pairs with no matches.
Add `--submitted` to compare original junctions before repair. `vdjdb submission` includes a
compact overlap screen; CI attaches the complete corpus report to each pull request.

Within each species and `(epitope, mhc.a, mhc.b, mhc.class)`, compare distinct alpha junctions,
beta junctions and paired junctions separately. For two sets of sizes `n1` and `n2`, report the
intersection `k`, size product `n1*n2`, containment `k/min(n1,n2)` and smoothed ratio
`(k+1)/(n1*n2+1)`. Repeated cells, peptide variants and technical replicates must not inflate the
set sizes. Missing chains do not match. Donor identifiers are scoped to their publication;
blank donor identifiers mean unknown, not a shared sample. Chunk-level sizes can pool donors.

The operational review threshold is at least five shared junctions and 20% containment of the
smaller set. It requests source tracing; it is not a significance test or a deletion rule.
Inspect smaller datasets even when they fall below the threshold. Use YLQPRTFLL and NLVPMVATV
as empirical comparison strata, with matching chain mode and restrictions, and include zero-match
pairs. These database comparisons include follow-up studies and reused data, so they do not by
themselves estimate independent-donor collision probabilities.

A sequence generation probability is not its probability in a selected tetramer-positive
repertoire. Selection, expansion, public clonotypes, chain pairing and ascertainment affect
sharing. Do not turn a typical `pgen` range into a universal overlap cutoff or multiply it by
chunk sizes to claim a p-value. See the [OLGA methods paper](https://doi.org/10.1093/bioinformatics/btz035).

For every flagged group, trace the manuscript, table, donor/cohort, receptor construct and assay.
Compare both chains' V/J calls and all `method.*` and `meta.*` fields before deciding:

- **Within-study validation:** the same sample or construct measured by different experiments.
  Preserve the reported observations and describe validation within that study.
- **Independent validation:** separately obtained samples or independently performed experiments
  confirmed by source evidence. Preserve both reports; state whether they share a construct.
- **Reused measurement:** the same source table or assay reprinted under another identifier.
  Correct provenance or repeated curation on a data-issue branch; do not count reuse as independent validation.
- **Conflicting attribution:** unmatched samples, assay labels or references. Resolve each matched
  observation against its original source, splitting mixed-source records when the evidence supports it.

An overlap flag alone never rewrites references, merges records or changes confidence scores.
Record per-observation decisions and source locations on the data issue before changing a chunk.

`paired.pmhc` additionally compares distinct paired-receptor/pMHC observations across peptides
within species. This catches reprinted mutational scans even when each peptide has fewer than
five receptors. Its denominator counts receptor/pMHC observations, not distinct receptors.

Measured on 2026-10-07, using repaired junctions, different reference IDs, matched species/pMHC
and both set sizes at least 20 (including unknown or pooled donors):

| Peptide | Mode | Group pairs | Shared junctions summed | Size products summed | Pooled ratio +1 |
|---|---|---:|---:|---:|---:|
| NLVPMVATV | alpha | 27 | 63 | 4,114,595 | 1.56e-5 |
| NLVPMVATV | beta | 98 | 222 | 11,063,981 | 2.02e-5 |
| NLVPMVATV | paired | 34 | 0 | 4,619,954 | 2.16e-7 |
| YLQPRTFLL | alpha | 31 | 241 | 609,001 | 3.97e-4 |
| YLQPRTFLL | beta | 45 | 248 | 446,070 | 5.58e-4 |
| YLQPRTFLL | paired | 18 | 16 | 353,274 | 4.81e-5 |

The pooled ratio is `(sum(shared)+1)/(sum(n1*n2)+1)`, with one pseudocount for the aggregate.
Junctions shared by several group pairs are counted once per comparison. These are corpus
screens, not verified independent-donor benchmarks. Regenerate the full TSV for current counts;
singleton sets have pseudocount-dominated ratios and should not calibrate a large-repertoire screen.

## Classify validation evidence

Keep these evidence classes separate; a record can have several kinds of support. When a single
label is needed, use the order 4, 2, 3, 1:

| Class | Evidence required | Browser selection |
|---|---|---|
| 1. Observed only | Initial capture/identification, without reported subsequent verification or replication | No additional validation requirement |
| 2. Observed and validated | A reported receptor verification assay, such as cloning followed by pMHC staining, target stimulation or direct binding measurement | Inspect `method.verification` and assay confidence |
| 3. Repeated within a study | Matching receptor/pMHC with distinct reported experiments, sampling time points, replicas or donors within the same reference | Same study validation |
| 4. Independently corroborated | Matching chain/receptor and pMHC in author-independent references, with reused-source observations excluded | Independent validation |

Replication within one study is stronger than one observation but does not replace receptor
verification. Different assays on the same sample are valid within-study evidence. A method name
change or a curation serial alone is not proof of a second experiment.

`vdjdb overlap` also writes a per-record `*-categories.tsv` audit. It uses nonempty
`method.verification` for reported assay validation and distinct nonempty replica/subject identifiers
for documented within-study repetition. Replication that the source reports only in prose requires
curation; blank identifiers do not establish it. The audit includes independent support separately
for alpha and beta. One supported chain does not establish independent capture of the paired receptor.

Author independence uses ordered PubMed author lists in `proofreading/pubmed_authors.tsv`:

- Senior authors must differ, and neither senior author may appear in the other paper's author list.
- Shared authors must be fewer than one third of the smaller list, using the strict inequality
  `3*shared < min(n_authors_1,n_authors_2)`.
- Missing, truncated or consortium-only senior-author metadata yields unknown. Non-PMID references
  require source-based review; do not treat missing metadata as disjoint author lists.

Names are matched conservatively by surname and first initial, with case, accents and punctuation
folded. Homonyms can require manual disambiguation. The author screen is a laboratory-independence
rule, not proof that an assay was repeated: a collection can reprint another laboratory's data.
High repertoire overlap requests original-table tracing even when authors qualify as independent.

Refresh this input explicitly with `vdjdb refs-authors`; review and commit the retrieved table.
Builds and CI read it offline. The audit excludes cross-reference support between high-overlap
reference pairs until source tracing establishes independence, while preserving direct verification
and within-study evidence. Other qualified supporting references remain available.

The browser currently reads `evidence.validation.same.study` and
`evidence.validation.independent`. The historical build's independent flag counts references per
chain/epitope, while its same-study flag has no producer. The new categories audit checks these
claims using assay, replicate and author metadata; it does not silently change the published
confidence scores, motif tuning objective or historical evidence flags.

Applying the author rule to the size-at-least-20 comparison above leaves the following paired
junction results. Unknown/pooled donor groups remain included; laboratory independence does not
by itself prove independent samples or assays.

| Peptide | Author-qualified group pairs | Shared pairs summed | Size products summed | Pooled ratio +1 |
|---|---:|---:|---:|---:|
| NLVPMVATV | 20 | 0 | 1,384,473 | 7.22e-7 |
| YLQPRTFLL | 11 | 13 | 186,835 | 7.49e-5 |
