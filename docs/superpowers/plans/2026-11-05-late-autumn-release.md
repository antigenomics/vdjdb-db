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

## Decision ledger

- Ruling: retain date tags and optional NOTE suffix support, with 2026-11-05 proposed for the ordinary Late Autumn snapshot. Public history has 43 releases with several naming variants and no published CalVer release; role-based consumers do not require a naming migration. Do not create a ZENODO twin. Final tag uniqueness will be checked again at publication.
- Ruling: execute the plan inline as requested, without a separate plan-approval ceremony. The user explicitly asked to set a goal and follow the plan.
- Release tooling changes remain uncommitted while the complete change set is prepared. No release, tag or download-pointer update has been published.
- Live Aldan-3 project storage has about 35 GiB free; select a bounded candidate/test footprint before transfer. Existing home-directory CI outputs are old and contain plain tables; they are not the selected new release candidate.
- Verified visual defects and corrections: SARS-CoV-2 summary now excludes 378 human slim rows labelled SARS-CoV while retaining the 26,192 SARS-CoV-2 rows; complete epitope sequences replace three-residue graph identities, preventing unintended merged strata. All source record counts and thresholds are unchanged. Self-antigen panel height increases from 20 to 40 inches to retain readable labels; repelled labels were trialled and rejected after visual inspection, with no new R dependency retained. Corrected the autoimmune typo and explained the actual per-combination thresholds and log2 ribbon widths.
- Current checks: 53 focused release/finalization/prepared-asset/summary tests pass (3.25 seconds); full R render passes the 9-heading/5-table/8-image fingerprint with the intentional self-panel height change. Final caption edits still need their final render. All eight initial static panels and both splash themes were visually inspected; label-column strips from the taller panel are readable. Splash: 105 rows, 105 distinct junctions, no duplicate epitope/motif/junction keys.
- Extracted summary tables have 7, 87, 226, 203 and 84 data rows respectively, with zero exact duplicate rows and zero malformed row widths. Individual content/layout review remains pending.
- Compatibility risk found by source inspection: the historical JVM automatic downloader unzips the pointer asset and opens plain vdjdb.txt directly; its code does not resolve gzip or a versioned archive subdirectory. The updated web application's configured local loader handles gzip. This must be resolved or explicitly accommodated before finalizing the public pointer; local web loading alone does not prove automatic legacy updates.
- Reviewed sealed structure ZIP retrieval was previewed and completed locally; its digest still needs verification before packaging/transfer. Scientific structures and motifs are not being recomputed.
