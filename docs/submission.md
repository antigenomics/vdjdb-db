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
V/J calls and junctions, epitope and MHC, then inspect all observation metadata. Different donors,
methods, subsets or publications are independent observations. Repeated curation of the same
publication is different from independent replication. Preserve existing PDB and mixed-paper
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

The repository's [curation skills](https://github.com/antigenomics/vdjdb-db/tree/dev/skills) guide
extraction, formatting, proofreading, publication, harmonisation and duplicate review. They use the
same CLI and authority tables as the build. Existing authorization for a batch applies throughout;
ask for new scientific decisions, not repeated permission for already-authorized routine actions.

Metadata identifiers must be reported in the cited publication or its supplementary tables.
Do not copy generated export identifiers, joined identifier lists or reference-derived labels into
`meta.clone.id`, donor, study or epitope fields. If the paper does not supply an identifier, leave
the field blank. Use the PMID in `reference.id`; do not encode it as a clone identifier.
