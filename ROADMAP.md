# VDJdb build roadmap

## Current work

The Python build migration is complete. Publication curation, metadata review and release
preparation use the [submission guide](docs/submission.md), [build commands](docs/builds.md)
and [integration checks](docs/build-integrity.md).

Use the [issue queue](https://github.com/antigenomics/vdjdb-db/issues) for current submissions and
unresolved defects. Completed migration checklists and dated issue-count snapshots have been
removed. Their history remains in Git. Local execution state is in the untracked
`ROADMAP_local.md`; `CHANGELOG_local.md` records completed integrations.

Sections 0-12 retain their numbers so existing code and documentation references still resolve.
Method definitions belong to `docs/`, not to a second specification in this roadmap.

## 0. Hard rules

1. Legacy motif files are positional contracts: 27 columns in `motif_pwms.txt`, 19 in
   `cluster_members.txt`. Keep their order and types compatible with `vdjdb-web`.
2. `vdjdb.meta.txt` must match `vdjdb.txt` column-for-column and in order. Generate metadata from
   the schema registry.
3. Use batch APIs once per batch, including `vdjtools.model.annotate_junctions`, `TCREmp.embed`
   and `vdjtools.overlap.tcrnet`. Do not wrap internally threaded tools in a per-record pool.
4. Deduplicate expensive deterministic annotation inputs within each build, then join results
   back to every observation. Preserve independent publications and distinct experiments.
5. Backgrounds are build inputs. Ship derived statistics, never background sequences or samples.
6. Empty string is the only missing marker for string fields in the pipeline.
7. Reproduce the same bytes across processes, hosts and core counts. Use `vdjdb.config.SEED`,
   sort unordered results, give ties an explicit order, and test every sampling or parallel stage.
8. Optimize in order: one Polars/NumPy expression, one batched call into an existing native
   library, then new native code. Profile before adding concurrency.
9. Recompute every output from source inputs on every build. Do not retain computed
   intermediates for reuse. Downloaded raw inputs and committed, reviewed reference tables are
   inputs; refresh reference tables separately from builds.

## 1. Build architecture

`chunks/`, `patches/`, `proofreading/` and the identity registry are inputs. Assembly produces
linked definitive tables; legacy and AIRR files are projections of those tables.
See [build and release](docs/builds.md), [outputs](docs/outputs.md) and the
[completed library migration](docs/builds/library-pipeline-migration.md).

The retired Docker/Groovy/GitLab pipeline is historical reference. Its defect list is no longer
an implementation backlog. `attic/` remains because schema tests read its original declarations.

## 2. Release assets

The release produces primary, legacy and AIRR bundles, a manifest and checksums.
[Build outputs](docs/outputs.md) specifies members, column orders and compatibility.
Legacy member basenames remain `vdjdb.txt`, `vdjdb.slim.txt` and `vdjdb_full.txt`.

## 3. Remaining work and release gates

### 3.1 Consumer compatibility

The old vdjmatch manifest-selection gate is resolved: published vdjmatch 0.4.0 selects the
requested bundle role and validates its manifest. Before publishing a release, verify the
intended consumers against its actual bundles using [build integrity](docs/build-integrity.md).
Do not infer a published database release from a merged build change.

### 3.2 Optional runner deployment

`build.yml` supports choosing a runner. No self-hosted runner is registered as of the
2026-10-17 inventory. Registering one and verifying identical canonical outputs on hosted and
self-hosted runners remain optional deployment work, requiring the author's request.
Hosted CI remains the verified release route.

### 3.3 Deferred consumer features

Reported-peptide assessment is implemented. Comparing flank variants in consumer queries,
representing alternative plausible class-II registers, and choosing structure-ready ligands
with their required flanks remain separate follow-ups. A predicted core does not authorize
trimming a submitted peptide or merging recognition assignments.
Optional corpus services and alternative weighting are described in
[the corpus specification](docs/standards/corpus.md), rather than repeated as migration tasks.

## 4. Integration workflow

Use `master` -> `dev` -> feature or chunk branch -> `dev` -> `master`.
Chunk changes and code changes are separate commits. A source change updates the record registry
before CDR3 repair. Commits resolving a tracker issue include `Closes #N`.

Run focused checks first, then the complete combined-corpus checks before promotion.
See [submission](docs/submission.md) and [build integrity](docs/build-integrity.md).
A source addition or correction must explain its measured release differences.

## 4a. Issue tracker and curation queue

The [live tracker](https://github.com/antigenomics/vdjdb-db/issues) owns the queue. Read each
issue's labels and disposition: a paper submission is not automatically a software defect.
Do not maintain another dated issue-count table here.

`pending/` and `withheld/` preserve excluded submissions. Closing an issue or completing a
migration does not make those observations positive database records. Their admission requires
source-level evidence and the [submission checks](docs/submission.md).

## 5. Release comparison

`vdjdb diff` compares member sets, raw and canonical digests, and keyed row/cell differences.
Every changed cell and every changed row bucket needs a declared reason and measured count.
A changed predicate, stale count or undeclared difference fails the gate.

[Build and release](docs/builds.md) defines canonical equality, identity-column corrections,
`[[row_delta]]` and members measured by other instruments. Motif quality and dashboard checks
remain required even when those members are excluded from cell attribution.

## 6. Legacy deprecation path

Legacy stops shipping only when all of the following hold:

1. Every known consumer supports manifest role selection.
2. Consumers use `latest-version.txt` only for historical URLs.
3. `vdjdb-web` reads the new format's generated metadata.
4. Two consecutive releases ship both formats with comparable use of the new format.

Until then, `latest-version.txt` points at the legacy zip. Legacy is emitted from the definitive
tables, never assembled through a second pipeline.

## 7. Measurements and coordinate conventions

Current measurements come from the selected build's reports and release tests. Historical
benchmark values remain in Git and the reviewed [clustering measurements](docs/clustering.md).
Do not use an old migration corpus size or timing as a description of the current database.

### 7.1 Coordinate spaces

VDJdb `cdr3` is junction space, including Cys104 and Phe/Trp118. AIRR `cdr3_aa` excludes both
anchors. All coordinate conversions belong to `vdjdb.convert.coords` and need round-trip tests.

| Space | Convention |
|---|---|
| VDJdb `v.end` / `j.start` | 0-based amino-acid indices in junction space |
| vdjtools `Scenario.v_end` / `j_start` | 0-based half-open CDR3 nucleotide space |
| arda `annotate.dmap` | 1-based closed junction nucleotide space |
| AIRR | 1-based closed sequence space |

See [AIRR mapping](docs/standards/airr-mapping.md) and [CDR3 repair](docs/standards/cdr3-fixing.md).

## 8. Motif inference contracts

[Clustering](docs/clustering.md) and [denoising](docs/denoising.md) define algorithms, formulas,
parameters, observation units and acceptance. The pointers below preserve existing section references.

### 8.1 TCRNET neighbour statistics

Neighbour counts, finite-background p-values and degree-zero handling are specified in
[denoising](docs/denoising.md).

### 8.2 Comparators and metric definitions

Compare the same cohort under the same metric definitions. [Clustering](docs/clustering.md)
distinguishes published reference partitions, measured alternatives and production algorithms.

### 8.3 Radius estimation

[Clustering](docs/clustering.md) defines the production radius and diagnostic Kneedle curve.
Do not describe a diagnostic knee as the production estimator.

### 8.4 Clustering topology

Fit chain-wide geometry on deduplicated receptors using stable sampling without replacement,
then cluster within each epitope. [Clustering](docs/clustering.md) specifies these two scopes.

### 8.5 Motif emission

[Outputs](docs/outputs.md) and [denoising](docs/denoising.md) specify PWM normalization,
background fallback, imputation flags and membership. Every emitted member needs its motif.

### 8.6 Length strata and cluster ids

[Clustering](docs/clustering.md) and [outputs](docs/outputs.md) specify length-stratified legacy
projection and cluster identifiers. Keep the web parser's positional contract.

### 8.7 Background choice

Use amino-acid controls for amino-acid neighbourhoods. [Denoising](docs/denoising.md) defines
control populations, sampling and denominators. Fetched controls are inputs, never outputs.

### 8.8 Memory and throughput

Fit and project embeddings in bounded batches as specified in [clustering](docs/clustering.md).
Rebuild every fit from the current distinct receptor cohort. Store timing reports, not reusable fits.

### 8.9 Evaluation helpers

Use [clustering](docs/clustering.md) for purity, retention, coverage and Mayer metrics.
Metrics with different unclustered-record conventions are not interchangeable.

## 9. Source and method authority

Publications establish experimental observations. [Column definitions](docs/standards/columns.md)
explain representative peptide parents, source recovery and marked modal MHC completion.
Predictions and inferred annotations remain distinguishable from reported measurements.

[Outputs](docs/outputs.md), [confidence scoring](docs/standards/confidence-score.md),
[terminology](docs/standards/terminology.md) and the schema registry own the remaining contracts.

## 10. Record identity and evidence

### 10.1 Record identity

`record_id` is allocated and survives an unambiguous amendment; `content_hash` detects changed
content. [Identifiers](docs/standards/identity.md) defines reconciliation and ambiguity handling.

### 10.2 Evidence model

Keep observation-level evidence and independent-publication support separate. See
[outputs](docs/outputs.md) and [confidence scoring](docs/standards/confidence-score.md).

### 10.3 Identity levels

Clonotype, clone, pMHC and epitope ids are derived from their own keys. Record ids use the
committed registry. [Identifiers](docs/standards/identity.md) defines all five levels.

### 10.4 Lifecycle

[Identifiers](docs/standards/identity.md) defines amendment, retirement, replacement and lookup.
Never allocate from an empty registry during an ordinary build.

### 10.5 Consistency checks

`vdjdb identity check` and the permuted-input tests enforce the
[identifier invariants](docs/standards/identity.md#invariants). Run them before integration.

### 10.6 Promiscuity annotations

MHC presentation predictions are annotations, never identity keys. See
[identifiers](docs/standards/identity.md) and [outputs](docs/outputs.md).

### 10.6a Reported peptides and binding cores

Preserve reported assay peptides. Allele-conditioned cores, presentation scores and uncertainty
belong to `epitope_assessment`; they do not rewrite identity or prove equivalent recognition.
See [outputs](docs/outputs.md) and [assessment commands](docs/builds.md).

### 10.7 Legacy structure ids

Preserve the legacy `TCR_hash` recipe and structure links. See
[identifiers](docs/standards/identity.md#legacy-structure-ids).

## 11. Tuning and validation

### 11.1 Independent-study support

Keep tuning and held-out validation separate. [Denoising](docs/denoising.md) defines the
independent-study objective, display exclusions and two-stage acceptance. [Clustering](docs/clustering.md)
defines matched-cohort comparisons and the current quality gates. A corpus change does not
justify silently resetting a historical baseline or relaxing a tolerance.

### 11.2 Held-out validation

Held-out records stay outside the repository and release bundles. `vdjdb.validate.guard` checks
both names and content shape. Read the set only through `$VDJDB_TCRVDB` and report aggregate
metrics; do not tune on it or emit record-level verdicts. See [denoising](docs/denoising.md)
and [outputs](docs/outputs.md).

## 12. Completed migration reference

The phase 0-17 implementation checklists are complete, apart from the optional runner deployment
in section 3.2. The [pre-cleanup plan](https://github.com/antigenomics/vdjdb-db/blob/c73919d108a8dc6fd5d7d27f7d49efbbd99d78ea/ROADMAP.md)
retains their history and acceptance evidence. They are not an active task list.

### Phase 14 release policy

The former phase 14 item 4 remains a release constraint: do not emit `-ZENODO` twin tags or
automate deposition. A DOI supplied by the author is a separate input.

## Record scoring v3.2 release

V3.2a exports assay-validation probability and v3.2b exports independent-study
validation probability, with separate frozen artifacts and output identities.
Offline fitting is in `processing/fit_record_score_v3_2.py`; private labels use
`TCRVDB_PATH` or `--tcrvdb-path`. See the
[scoring specification](docs/standards/record-scoring-v3.md).
