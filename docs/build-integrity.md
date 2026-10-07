# Build integrity

A release must preserve the source observations, satisfy the output contracts and explain its
changes from the previous release. This guide covers code and CI. Use the
[submission guide](submission.md) for source extraction and curation decisions.

## Keep inputs and outputs separate

Build from a recorded commit and dependency lock. Preserve literal submitted junctions and V/J
calls alongside their repaired values. Assign record identity before repair. Keep independent
publication reports and distinct experiments separate.

Recompute every derived output on each build. Downloaded references are inputs; generated indexes,
intermediate tables and model results are not retained between builds. Use the package seed and
explicit sorting so hash order, scheduling and worker count cannot change the answer.

A build change that changes record content requires registry reconciliation, even without a chunk
edit. Compare record IDs, scores, original calls and amended content separately. Report unresolved
annotations alongside successful repairs.

## Validate a candidate

The [build commands](builds.md#local-build) produce fresh tables, projections, motifs and the offline
summary. Validate against the release configured by CI. Set `VDJDB_REFERENCE_ZIP` when running the
suite; otherwise release tests are skipped.

Check these independently:

- Strict input QC and record-identity invariants.
- Primary, legacy and AIRR schemas, joins and bundle members.
- Complete release comparison, with every changed bucket measured and explained.
- Motif quality, positional consumer contracts and repeated-build determinism.
- Summary structure, embedded plots and visual checks.
- Documentation generated from the same schema and built with warnings treated as errors.

Legacy table columns and metadata must agree in order. vdjdb-web also reads motif files positionally;
validate its parser before changing them. A parser check does not establish that the application
starts or that a deployment has been updated.

Do not refresh expected counts or loosen tolerances automatically. A source addition, a correction
and a code regression need different explanations. Keep known numerical or scientific limitations
on focused issues and state which validation remains unavailable.

Measure slow stages with input size, wall time, core count and memory. File bottlenecks where the
slow code is maintained. Optimize with vectorized expressions and batch calls before adding workers.

## Integrate and publish

Merge the validated candidate through dev, then promote the final tested tree to master when
requested. Required checks must pass on that tree. Verify the master build and inspect its artifacts
after merging; a green fast check alone does not establish release integrity.

The dashboard follows successful master builds. A build-triggered documentation run must download
that exact run's artifact and fail if it is unavailable. Other documentation updates select the
latest successful master build. Verify the published fragment against its selected artifact.
Dated release downloads remain separate snapshots.

Inspect all three bundles with `vdjdb release --dry-run` before publishing an authorized tag.
Release notes come from `vdjdb changelog`. Back-sync hotfixes to dev, reconcile completed issues and
remove only integrated branches and obsolete generated output. Preserve source inputs, evidence
and unmerged work.

The repository's [build integrity skill](https://github.com/antigenomics/vdjdb-db/tree/master/skills/vdjdb-build-integrity)
provides the operational checks for this lifecycle.
