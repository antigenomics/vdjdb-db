---
name: vdjdb-build-integrity
description: Validate VDJdb build code, CI, output contracts and release promotions against a released corpus. Use for build changes, dependency upgrades, integration audits and release orchestration; route manual source extraction and curation decisions to the curation skills.
---

# Verify the build and release lifecycle

Read repository `AGENTS.md`, [source authorities](../AUTHORITIES.md),
[build and release](../../docs/builds.md),
[outputs](../../docs/outputs.md) and [build integrity](../../docs/build-integrity.md).
Read `docs/denoising.md` before motif changes. Inspect the current local roadmap and worktree
before starting; preserve existing source work and unresolved follow-ups.

## Define the candidate

Record the base commit, candidate tree, dependency lock and reference release. Use the reference
configured by CI and the comparison declarations; changing the baseline is a reviewed input change.
Keep chunk edits on their own data branches. A code/CI audit does not silently normalize sources.

The build reads committed/reviewed inputs and recomputes its outputs. Downloads are inputs;
computed indexes, tables and model results must not persist between database builds. Exercise
success and failure cleanup when introducing temporary storage. Use published dependencies in
the build, with cross-repository fixes released before their version bounds move here.

Inspect workflow event and path filters: changes to any measured input, schema, dependency or
build implementation must reach full validation. A fast chunk check or green documentation job
is not proof of a full database build. Check artifact provenance by exact commit/run ID.

## Test the affected contracts

| Change | Required evidence |
|---|---|
| Schema or projection | Generated schema, dtype/null contract, joins, exact legacy column order and metadata agreement |
| Annotation or harmonisation | Literal source-to-build original CDR3/V/J equality, repair cases and unresolved flags, identity amendments |
| Deduplication or identity | Pre-repair identity, retained record IDs, within-publication decisions and preserved independent reports |
| Sampling, clustering or parallelism | Repeated digests in separate processes, seed/hash/worker-count checks, relevant motif metrics |
| CI, packaging or publication | Event/selector cases, all three bundle layouts, exact tested artifact and live consumer/publication checks |

Use focused tests for the changed behavior, then run release-aware checks for any measured corpus
value. A plain test run can skip release markers. Build fresh tables and projections, run strict
QC, identity and corpus invariants, the complete release comparison, motif metrics and the offline
summary checks. Validate documentation with Sphinx warnings treated as errors.

Re-measure declared differences and explain each affected bucket. Keep source/behavior corrections
separate from baseline counts; do not regenerate expectations automatically or widen tolerances
merely to obtain a pass. For a build-only change, compare input subtree hashes as well as results.
When assembled `content_hash` changes, run the documented registry update even if no chunk changed.

Legacy metadata must agree with table columns in order. Check the actual vdjdb-web motif parser's
19/27-column positional contracts when changing its inputs. Distinguish parser/schema validation
from starting the Play application and from deploying it. Check primary and AIRR projection
invariants independently of legacy compatibility.

Measure stage time, input size, machine/core count and peak memory where relevant. Vectorize,
then batch existing native calls before parallelizing. File a measured bottleneck in the owning
repository with its profile and the limit of the chosen approach.

## Promote the tested tree

Follow feature/data branch -> dev -> master; master accepts dev or an authorized hotfix.
Required green checks apply to the final tree, not an earlier head. Inspect the full CI stages,
test totals, skipped tests and warnings; preserve unresolved scientific or numerical defects on
focused issues. A waived gate is a limitation, not a successful test.

Before every PR merge, check both CI results and code review. Inspect the final diff and all
available review summaries, inline comments and unresolved threads; resolve substantive findings.
Record the review outcome together with CI evidence in the local execution record. An empty review
list does not replace reviewing the diff. Reviews must cover the final head or a proven-identical
tree; after a material revision, review the changed diff again and obtain fresh applicable checks.
A stale review or check does not cover subsequent changes.

Promote to master only within the user's authorization. Verify the master build after merging.
Documentation triggered by a successful build selects that exact run; other documentation updates
select the latest successful master build. Do not combine mutually exclusive run/branch selectors.
A missing required artifact must fail publication. Check the deployed summary bytes against the
selected artifact after the deployment succeeds.

Prepare all three bundles with `vdjdb release --dry-run`, inspect manifests and generated changelog,
and publish a tag only when requested. Keep the live-master dashboard distinct from dated release
snapshots. Update local execution/merge records, reconcile completed issues, back-sync hotfixes to
dev and delete only ancestry-verified integrated branches. Retain source proofs and unmerged work.

For a final corpus audit, follow [metadata and reference checks](../../docs/build-integrity.md#check-metadata-and-references).
Record source duplicate counts separately from assembled observations and retained experimental
replications. Reopened curation issues do not invalidate a build; retain their agreed scope and
proofreading/maintenance labels without changing supported source observations.

## Choose assessment scope

Use [check scopes](../../docs/builds.md#choose-the-scope-of-a-check) to select the work needed
for the current decision. Start chunk triage with `vdjdb qc <chunk>` and `vdjdb submission <chunk>`.
For a peptide/MHC question, run `vdjdb assess-epitopes <chunk>` with the separately fetched pinned
reference, or request opt-in `assessment.yml` CI on those tracked paths when CI is authorized.
For an existing built corpus, `--tables` recomputes assessment from records; previous predictions
are never reused. Selected-input support counts describe that input, not the full database.

Read `reported`, `assessment_status`, `allele_resolution`, predicted presenter and core fields
with their [output contract](../../docs/outputs.md).
Keep receptor, MHC and parent species separate. Preserve literal assay peptides and source-supported
parent genes/species. Predicted cores and TCR-facing sequences are advisory comparison representations,
not measured minimal recognition epitopes or contact maps. Class-II registers depend on the molecule;
class-I footprint cores can omit insertions. Do not trim source peptides, infer an unreported DP/DQ
partner, merge records or assign motif groups from these predictions. Structure inputs require the
full appropriate ligand and source evidence.

Partial checks accelerate review. Full integration still recomputes junction annotation, corpus,
motifs, dashboard and release measurements on the combined candidate; a partial artifact does not
establish those contracts. Record the scope, selected inputs and revision with each result.
