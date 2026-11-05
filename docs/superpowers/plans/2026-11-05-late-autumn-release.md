# Late Autumn Release Implementation Plan

> Execute in this chat using superpowers:executing-plans. Track completion here; do not publish before all acceptance gates pass.

**Goal:** Publish “Late Autumn Release” on GitHub with the new database format and legacy/web bundles, verified on Aldan-3 and manually proofread.

**Architecture:** Reuse the existing build, release packager, manifests, consumer loaders and web installer. Pin one source commit and one candidate cohort; validate the exact bytes that will be published. Keep scientific generation separate from packaging and interoperability checks.

**Tech stack:** uv/Python, existing R/pandoc summary, GitHub Releases/Actions, Aldan-3 client/SLURM, Python VDJmatch and Scala/Angular VDJdb-web.

**Spec:** User request in this chat: new release, new-format and legacy ZIPs, date tag, title “Late Autumn Release”, Aldan-3 interoperability, manual visual proofreading; review date tags with optional NOTE suffix before deciding.

## Constraints and review focus

- Preserve accepted scientific results and all existing source edits. No new controls, model fitting, thresholds or unrelated analyses.
- No public private provenance, filesystem paths, credentials or identifying stakeholder labels.
- Generated TSV/TXT above 1,000,000 bytes must be deterministic gzip, including ZIP members.
- Tags/assets are immutable. Do not create a Zenodo twin or mutate old releases.
- Use scheduled compute nodes for Aldan-3 execution; frontend only for orchestration/I/O. Use project storage for bulk inputs, preview transfers and verify digests.
- Confirm exact tag selection, multi-ZIP role selection, gzip handling, full-pair identities and prepared structure bindings. Never treat earlier Stratos checks as Aldan-3 evidence.
- Preserve complete full-length epitope identities even where visual abbreviations collide; do not silently hide labels to cure overlaps.
- Publication, source integration and download-pointer finalization require independent verification. Do not infer one from another.

## 1. Release/tag and source inventory

- [x] Inspect current checkout, repository instructions, release workflow and existing candidate evidence. Working branch is clean; it contains release tooling not yet on master.
- [x] List published GitHub releases: public history includes plain dates, date-NOTE suffixes and historical pyvdjdb-date prefixes; no published vYYYY.MM.PATCH release was listed.
- [ ] Inspect full tag parsing/selection paths in producer, Python vdjmatch, web installer, legacy JVM consumer and workflow/finalizer tests. Python acquisition uses an exact URL-quoted tag, manifest role and digest; the web installer permits safe date/suffix tags and enforces exact identity. Legacy web local loading uses resolved plain/gzip paths, not tag parsing; automatic JVM download still needs inspection.
- [ ] Compare choices: YYYY-MM-DD for an ordinary snapshot; YYYY-MM-DD-NOTE for a distinct explicitly named variant; vYYYY.MM.PATCH for a deliberate policy migration. Recommend the least disruptive interoperable convention and record the decision before editing. Candidate date is 2026-11-05; title is independent of tag.
- [ ] Inspect current remote master/dev/feature divergence, required checks and available CI artifacts. Select the exact source commit and integration route without discarding local work.

**Acceptance:** One documented tag, unoccupied remotely, accepted end-to-end by producers and consumers; truthful source identity and explicit distinction between assembly, scoring and packaging revisions.

## 2. Candidate provenance and build scope

- [ ] Verify existing BUILD/SOURCES receipts, actual input hashes, structure metadata pin and build artifact identities. Existing recorded candidate has 315,921 records, 516,101 chains and 1,625,234 evidence rows; these are recorded counts, not a new-build claim.
- [ ] Determine whether the existing accepted build matches intended current inputs. Package completed outputs without rerunning analyses if it does. If a fresh database build is required, use the existing complete pipeline from raw inputs, with pinned environment and no intermediate reuse; record resource/disk budget and monitor completion.
- [ ] Check fresh-build consequences against the accepted previous release using existing diff and identity checks, explaining each material source difference. Stop before publishing weakened standing results.
- [ ] Check frozen confidence/motif/structure outputs and current IDs remain bound to the selected cohort. Do not regenerate accepted scientific outputs solely for packaging.

**Acceptance:** Cohort/source provenance is reconciled; required stage reports and identity invariants pass; no unexplained record loss, duplicate identity or mixed build stages.

## 3. Manual visual and table proofreading

- [x] Run existing summary structural/style check: 9 headings, 5 tables, 8 images match baseline.
- [x] Extract and visually inspect all eight current static panels on white backgrounds. Initial findings: crowded/overlapping flow labels, misspelled “utoimmune” in source prose, COVID source predicate includes SARS-CoV-1 as well as SARS-CoV-2. These findings are not yet fixed.
- [ ] Audit tables for duplicate keys, accidental abbreviations, typos, missing labels, inconsistent units/denominators and agreement with selected source counts. Check COVID specificity and aggregation explicitly; do not change scientific thresholds.
- [ ] Inspect splash light/dark variants and interactive summary, legends, tooltips, ordering, labels, clipping, all table sections and any release infographic. Distinguish repeated legitimate observations from duplicate curation.
- [ ] Trace each defect to owning code/callers. Correct only release-facing defects; use existing plotting backend and installed tools. For nontrivial logic add a focused failing regression check before the fix.
- [ ] Re-render affected summary outputs from the same selected inputs, rerun automated checks and visually inspect all changed panels/tables at readable scale. Keep the review record with findings, fixes and verification evidence.

**Acceptance:** No unresolved release-blocking visual/table defect; every changed display is manually inspected; final images and tables belong to the shipped candidate.

## 4. Build immutable candidate bundles

- [x] Make the minimum required tag-policy compatibility changes in release/bundle.py and owning workflow/help/docs. Date and date-NOTE tags now validate real calendar dates and safe suffix characters; CalVer remains accepted. Added exact tag/URL and immutability regression checks. Observed three new date-tag cases fail before the change; all 28 release/finalization tests pass afterward (2.34 seconds).
- [ ] Run focused release/manifest/finalization/summary tests with bounded timeouts. Preserve existing supported consumer behavior.
- [ ] Package primary, legacy, AIRR, web and paired web-structures when valid accepted preparation is available, using existing vdjdb release APIs. Include manifest.json, SHA256SUMS, lifecycle asset and concise release notes.
- [ ] Verify every bundle/member checksum, member uniqueness, safe names, CRCs, gzip thresholds, embedded download pointer and manifest/source/tag consistency. Confirm web table, metrics and structures remain sealed to one cohort.
- [ ] Measure disk footprint before copying full structure ZIPs; avoid duplicate archives and retain only requested final outputs and compact receipts.

**Acceptance:** All packaged bytes are verified and frozen for the Aldan-3 checks; required roles are present and cannot be confused by consumers.

## 5. Aldan-3 interoperability against actual bundles

- [x] Verify live access through aldan3 doctor; SSH and existing authentication work outside the local network sandbox.
- [ ] Inventory existing project paths, runtime versions, consumer source revisions, Java/Python environments and scheduler capacity; select isolated test destinations and resources.
- [ ] Preview then transfer only required source/environment recipes, candidate archives and public test inputs; compare destination hashes. Keep runtime/auth configuration separate.
- [ ] Submit bounded SLURM checks with discoverable logs and exit status. Exercise Python vdjmatch exact-tag/default/primary/legacy/AIRR acquisition and actual receptor annotation using existing public fixtures; validate row counts, full junctions, pairing and score/output semantics. Check historical and date-NOTE discovery compatibility separately where affected.
- [ ] Install candidate web roles through tools/install_release.py; verify loader counts, metadata, motif tables/charts and structure bindings. Run actual Scala loader/backend and annotation paths in the supported legacy runtime; include Browse/overview/motif/structure/export journeys against an isolated preview when feasible. Do not replace production configuration.
- [ ] Exercise legacy JVM vdjmatch dependency as used by vdjdb-web, not only the Python namesake. Use one batch of query receptors rather than per-record engine calls.
- [ ] Retrieve compact receipts, versions, exact archive digests and failure logs. Investigate any stalled/failed job and rerun only the affected check after correction. Record actual execution boundaries.

**Acceptance:** Both consumers actually run on scheduled Aldan-3 nodes using the final candidate bytes; expected counts, annotation fixture results and web journeys pass. No substitution of local or Stratos evidence.

## 6. Integrate and publish

- [ ] Review final changes for correctness/confidentiality, integrate through the repository's supported branch/check route and verify the visible integration checkout. Do not publish a tag that falsely claims to represent a different source build.
- [ ] Confirm tag/title/date, previous-release comparison, role descriptions, counts and known caveats in release notes. Include concise usage instructions for new-format and legacy/web consumers.
- [ ] Publish the immutable GitHub release titled exactly “Late Autumn Release”, following existing release protections and required checks. If a human GitHub environment gate appears, report the exact pending action without bypassing it.
- [ ] Read the published release back, verify tag points to pinned SHA, asset names/sizes and downloaded manifest/checksums; smoke-test real published acquisition on Aldan-3.
- [ ] Finalize latest-version.txt only after publication succeeds, preserve concurrent source commits and verify line 1 returns the intended legacy asset successfully.
- [ ] Deliver release URL, exact tag, asset list, compact Aldan-3 and visual-review evidence, and any remaining boundary. Mark the goal complete only after all required work is done.

**Acceptance:** Public release and final pointer resolve to verified assets; title/date/source match; final consumer checks use published assets.

## Current decision and evidence ledger

- Release convention: retain YYYY-MM-DD and safe optional NOTE suffixes; CalVer remains accepted for compatibility. Public history lists 43 releases and no published CalVer tag. Never create a Zenodo twin. Preparation crossed midnight, so the intended ordinary-release tag is now **2026-11-06**, title **Late Autumn Release**; check uniqueness again immediately before publication.
- Execute inline as the author requested; the goal remains active. Visible source branch is codex/legacy-web-interop. Commits 1f790ef and a56238d are pushed; subsequent joint-context summary fixes are under focused verification. No release, date tag or public download-pointer change has been published.
- Full build [42103107759](https://github.com/antigenomics/vdjdb-db/actions/runs/42103107759), source 1f790efdbc8d4cfd8d7d3d5d1eabe1e0c8fa1907, completed successfully. It retained reports/scores/summary, not core database tables, because standalone dispatch had no asset-tag input. Expose the existing candidate tag/web/previous inputs and run final source once to retain the requested actual ZIPs. This run also needs the later summary corrections; do not treat old local tables as that fresh deliverable.
- Static review found overlapping alluvial labels, three-residue epitope identities merging distinct epitopes, SARS-CoV rows in the COVID section, the autoimmune typo, inaccurate per-combination/log2 captions, and paired HLA stars disappearing as Markdown emphasis. Corrections preserve thresholds, complete epitope strings and safe HTML escaping. Self-antigen panel height increases from 20 to 40 inches. Automatic repelled labels were rejected after visual inspection; no R dependency was added.
- The working cohort contains 307,789 slim rows and 516,101 long chain observations. COVID selection retains 26,192 human SARS-CoV-2 slim rows and excludes 378 labelled SARS-CoV. These are working-build observations, pending final ZIP identity verification.
- All eight original panels and both splash themes were visually inspected. Splash: 105 rows, 105 distinct junctions, zero duplicate epitope/motif/junction keys. Offline R table proofs covered all 607 original summary rows (7/87/226/203/84), with zero exact duplicate rows or malformed widths. These proofs are content review, not a claim of browser CSS inspection. Changed local panels were re-inspected after the joint-context correction, including all self-antigen labels at native resolution, full MHC-II labels, heatmap and chord. Extra sector gaps and smaller chord labels remove remaining overlap. Final artifact equivalence remains to verify. All 51 changed/new table rows were visually read against the previously reviewed rows; final table sizes are 7/87/226/203/82, with no duplicate rows.
- Manual HLA-table review found two artificial cross-locus pairs. Source long annotations are correct: one slim receptor has separate DP and DR observations from separate references. Independent expansion of collapsed MHC/V lists fabricated associations and attributed both references to each. Static context panels and interactive V/HLA/HLA-table paths now join long observations to the existing six-field slim key, retain IDs, count distinct slim IDs and validate complete cohort coverage. A joint fixture failed before the fix and passes afterward, checking MHC/V coupling and per-context study counts through real R and Python paths.
- Matched working-cohort V/HLA audit: all 307,789 slim IDs retained; 516,101 long rows; heatmap cells 2,880 → 2,867; 13 unsupported cells removed; no added cells; 15 shared-cell counts corrected. Marginal cohort data and scientific model outputs are unchanged. R/Python full-cohort heatmap parity is exact across all 2,867 cells. Keep exact cell deltas privately and report final-artifact quantities after verification.
- Source nomenclature variations are retained, rather than silently changing curated organism/gene identities. One labelled note in the summary explains that aliases remain separate groups. HLA display retains full chain names and first-field resolution; chord labels retain TRBV names; its caption correctly says at least 50 records without changing the threshold.
- Reproduced and fixed gzip changelog reading for directories and ZIPs. Date-tag, immutable-tag, plain/gzip, prepared-asset and real R/Pandoc HLA checks pass. Final focused suite: 86 passed in 17.49 seconds. The workflow lifecycle selector now retains both physical plain and gzip forms; a focused glob test verifies both.
- Accepted structure ZIP: 1,320,957,573 bytes, SHA256 50000dce647be2c1e7298a377ba276e32d2c85f86b09c9be5d3d58881c9eb393. Previewed retrieval/transfer succeeded; the scheduled Aldan-3 source/hash receipt matches. No structure or motif analysis was rerun for this transfer.
- Aldan-3 project storage is tracked and isolated. Setup job 1509654 ran on a compute node with 8 CPUs/16 GiB: Python VDJmatch 1.2.2 release-loader suite 42 passed in 3.85 seconds, with resolved dependencies recorded privately. These are unit checks, not final-candidate interoperability.
- Runtime setup: overlay whiteouts are unsupported even on node-local ext4. VFS construction succeeded but copies approached the 20 GB task budget; cancelled job 1509667 was observed completing, and no scientific computation was cancelled. Native runtime extraction job 1509671 reuses the exact stored runroot, copies and checks the pinned Linux Java 8/Node 11/Yarn binaries, then removes only the isolated task-created image store before backend tests. Last observed project free space recovered to about 34 GiB. Monitor job and footprint; never infer cancellation/completion from transport timeouts.
- Historical JVM auto-download source opens plain tables and unzips into a fixed location, so gzip/versioned members need consumer support. Updated web local loading handles gzip. Author explicitly chose updating the JVM consumer. The isolated JVM maintenance branch adds shared gzip loading and safe versioned installation before pointer replacement. Aldan-3 job 1509681: transport/recovery test and four existing loader tests pass; JAR/install build passes. The Python checkout is preserved.

- Final fresh core build: CI run 42113249883 at 3c978e7, dated candidate inputs enabled. Later changes only adjust chord label spacing and correct interactive publication counting (2,156 collapsed reference sets → 845 distinct nonempty references on the working cohort). Re-render the summary from the actual final core tables, then package using the final source; retain separate core-build and summary/package source receipts. No assembly, motif or frozen-score rerun is needed for these presentation-only changes.
- Aldan-3 web retry job 1509679 passes 339 backend tests across 43 suites; one pre-existing test is cancelled. The alternate sbt downloader resolves the old repository stall; native pinned Java 8 runtime is retained in project storage. Actual final-candidate install and annotation checks remain open.

- Final fresh CI run 42113249883 succeeded and retained all four unprepared candidate roles. Artifact comparison found the accepted rich snapshot predates the documented shared-locus dependency updates (12c547c): 3 long junctions, 698 V-call cells, 2 legacy score cells, 1,000 hash cells differ. The accepted source is recorded in out/BUILD.json as assembly 5ffa45d and frozen scoring 69b326a; its 32 immutable deliverable inputs are hash-inventoried. The author has been asked to select preserved rich/structure deliverables or fresh current-source data plus fresh structure preparation; do not publish either choice before the answer. All 605 summary table rows agree across these cohorts, but that does not imply raw table identity.
- Required source lint failed on formatting introduced in the prepared gzip branch. Mechanical import cleanup and line wrapping repair all 64 diagnostics; Ruff and strict mypy pass. Focused release/compression/summary suite: 94 passed in 18.97 seconds. No scientific producing logic changed in this cleanup.
- JVM maintenance PR62/63 integrated; immutable tags 1.3.2 and corrected 1.3.3 published. Final 1.3.3 retains safe ancillary legacy assets and passes transport/recovery plus four loader checks on Aldan-3. Web source409b1855 resolves the selected JVM directory and published1.3.3; Aldan job1509685 passes339 tests across43 suites (one pre-existing cancellation). Web PR creation is denied through both current API integrations; investigate the current repository integration policy before altering master. Python1.2.3 PR64 to dev is merged with green CI; release PR65 to master is running final checks.
