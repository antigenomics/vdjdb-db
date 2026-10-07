---
name: vdjdb-build-integrity
description: Validate VDJdb build code, CI, output contracts and release promotions against a released corpus. Use for build changes, dependency upgrades, integration audits and release orchestration; route manual source extraction and curation decisions to the curation skills.
---

# Verify the build and release lifecycle

Read repository `AGENTS.md`, [build and release](../../docs/builds.md),
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

Promote to master only within the user's authorization. Verify the master build after merging.
Documentation triggered by a successful build selects that exact run; other documentation updates
select the latest successful master build. Do not combine mutually exclusive run/branch selectors.
A missing required artifact must fail publication. Check the deployed summary bytes against the
selected artifact after the deployment succeeds.

Prepare all three bundles with `vdjdb release --dry-run`, inspect manifests and generated changelog,
and publish a tag only when requested. Keep the live-master dashboard distinct from dated release
snapshots. Update local execution/merge records, reconcile completed issues, back-sync hotfixes to
dev and delete only ancestry-verified integrated branches. Retain source proofs and unmerged work.
