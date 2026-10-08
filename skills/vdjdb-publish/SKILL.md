---
name: vdjdb-publish
description: Publish proofread VDJdb chunks through dev, reconcile publication issues, preserve record identity, keep validation changes separate, and close completed imports while tracking unresolved observations and optional metadata follow-ups. Use when a chunk is ready to publish or an imported paper needs closeout.
---

# vdjdb-publish

Read [AUTHORITIES.md](../AUTHORITIES.md) and the
[submission guide](../../docs/submission.md). Use existing batch authorization throughout;
ask only for a new decision outside that scope.

## 1. Reconcile the submission

Search open and closed issues by PMID, DOI and URL, including bodies. Check existing chunks and
`proofreading/reference_ids.tsv`. Verify preprint-to-publication mappings with PubMed before
merging their provenance; matching sequences alone do not establish that two papers are one work.

Create a missing issue as `PMID:<id>`, label `vdjdb-records-paper-pending`, with the retrieved
PubMed NBIB citation and record link: author names, title, journal, date and pages in one linked
paragraph. Preserve source attachments and one concise blocker or disposition. Use a concise
title and DOI for an unindexed preprint. A multi-PDB paper has one issue; use
[structures](../vdjdb-structures/SKILL.md) and place supported rows in `chunks/PDB_Database.tsv`
with their publication reference and individual structure identifiers.
Public provenance cites only the paper, abstract, PubMed,
patents or PDB. Patent and mixed-paper inputs need their own manifest of references and issues.

Preserve existing chunk content. Audit both-chain V/J/junction, epitope and MHC matches together
with donor, method, subset, clone and other observation metadata. Never append or deduplicate
using beta junction plus epitope alone.

For an extension, apply the curator's scope to existing records. If they instruct keeping those
records, preserve their source cells and add supported new observations. Change an existing field
only when requested. Do not make re-proofreading or optional metadata completion a prerequisite
for the extension. Reconcile identity, run the combined build and finish the authorized merge.

## 2. Check readiness

```bash
uv run vdjdb qc chunks/PMID_<id>.tsv
uv run vdjdb submission chunks/PMID_<id>.tsv
```

Read the reports, including the confidence-score histogram. Retain source cells for routine
repairs that assembly performs. Say **resolved during database build**, and inspect the built
result instead of manually repeating the repair at import.

A noncritical metadata question is not an import blocker when the supported assay evidence already
determines the score and observation identity. Leave unknown fields blank and keep a focused
follow-up issue open. Pairing, epitope/restriction and assay-outcome ambiguities remain blockers for
the affected observations. Preserve unresolved material in `pending/` or `withheld/` as described
in the submission guide; explicit negative observations go to `chunks_negative/`.

## 3. Commit data separately

Start from current `dev`, never `master`, and name the branch for the chunk or data issue.
Normally use one PR per chunk. Group remaining chunks only when the user authorizes a batch;
retain existing PRs unless they explicitly request consolidation.

```bash
git switch -c codex/chunk-PMID_<id> origin/dev
uv run vdjdb identity update
git add chunks/PMID_<id>.tsv registry/records.tsv.gz
```

The registry diff must account for this submission only. Preserve old IDs; investigate unexpected
amendments or retirements before committing. Never stage unrelated files or reset another task's
staging area.

The commit describes changed files, per-file added/removed/amended row counts, the source issue,
why they changed, what the build shows, and who made any curation judgement. Use `Refs #N` while
material work remains, or `Closes #N` when the import is complete. A batch lists every chunk and
its issue. Never mix chunk edits with code or validation changes.

## 4. Validate the combined candidate

Keep release-comparison declarations and tests in a separate validation PR. Rebuild from source
and inspect `vdjdb diff --report`: each changed row bucket needs an explained count and note.
Do not refreeze all expectations automatically. Check tables, AIRR projection counts, reference
years and corpus/motif measurements when the import affects them.

A fixed released clustering scored on a larger corpus can change its measured statistics.
Distinguish that from a current-versus-release regression on the same cohort. Keep tolerances and
acceptance criteria unchanged unless the curator explicitly approves a different criterion.
Keep one reviewed baseline instead of duplicating corpus-dependent constants in tests.

Run required PR checks and full CI on the combined data-plus-validation candidate. Fix failures;
never merge a red candidate. If legacy compatibility was requested, test the new chunk with the
specified legacy checkout and retain the result in the local execution record.

## 5. Merge and close out

Target `dev`. After the combined full build and required checks pass, merge data, refresh the
separate validation branch on the resulting dev, and compare its complete tree with the green
candidate. Wait for the refreshed required check, then merge validation and verify the final tree.
Do not promote to `master` unless requested.

Before every PR merge, check both CI results and code review. Inspect the final diff and all
available review summaries, inline comments and unresolved threads; resolve substantive findings.
Record the review outcome and CI evidence in the local execution record. An empty review list
does not replace reviewing the diff. Checks and reviews must cover the final head or a
proven-identical tree; after a material revision, review the changed diff again and obtain fresh
applicable checks. A stale review or check does not cover subsequent changes.

Read issue comments and the pending/withheld/negative manifests before closing. Record merged PRs,
imported scope and any leftover rows. Close only when no material work remains; GitHub may not
auto-close issues when the target is `dev`. If optional metadata remains, keep a nonblocking
follow-up and remove the pending-paper label from the completed import. Leave umbrella issues
open while any component is still pending.

Verify the actual `reference.id` rows, including verified preprint/publication aliases and aggregate
chunks; a PMID-shaped filename or a matching receptor alone does not prove coverage. Closed imports
link the current path and scope. A duplicate links its active intake. An unsuitable submission has
a brief source-supported exclusion reason; a quarantined or incomplete submission remains open.
Build and CI issues need their implementation and validation evidence, not a publication chunk.

Report prepared, merged and deferred counts separately. Update the local roadmap and per-file
manifest. Delete only branches whose commits are integrated, with no open PR or attached working
tree; preserve unresolved source tables and diagnostic evidence.
