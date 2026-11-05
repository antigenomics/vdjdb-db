# Late Autumn Release Implementation Plan

> Completed 2026-11-06. Final publication and public-consumer evidence are recorded below.

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
- [x] Inspect full tag parsing/selection paths in producer, Python vdjmatch, web installer, legacy JVM consumer and workflow/finalizer tests. Python acquisition uses an exact URL-quoted tag, manifest role and digest; the web installer permits safe date/suffix tags and enforces exact identity. Legacy web local loading uses resolved plain/gzip paths, not tag parsing; JVM 1.3.3 now supports versioned gzip bundles.
- [x] Compare choices: YYYY-MM-DD for an ordinary snapshot; YYYY-MM-DD-NOTE for a distinct explicitly named variant; vYYYY.MM.PATCH for a deliberate policy migration. Recommend the least disruptive interoperable convention and record the decision before editing. Selected date is 2026-11-06 after preparation crossed midnight; title is independent of tag.
- [x] Inspect current remote master/dev/feature divergence, required checks and available CI artifacts. Select the exact source commit and integration route without discarding local work.

**Acceptance:** One documented tag, unoccupied remotely, accepted end-to-end by producers and consumers; truthful source identity and explicit distinction between assembly, scoring and packaging revisions.

## 2. Candidate provenance and build scope

- [x] Verify existing BUILD/SOURCES receipts, actual input hashes, structure metadata pin and build artifact identities. Existing recorded candidate has 315,921 records, 516,101 chains and 1,625,234 evidence rows; these are recorded counts, not a new-build claim.
- [x] Determine whether the existing accepted build matches intended current inputs. Package completed outputs without rerunning analyses if it does. If a fresh database build is required, use the existing complete pipeline from raw inputs, with pinned environment and no intermediate reuse; record resource/disk budget and monitor completion.
- [x] Check fresh-build consequences against the accepted previous release using existing diff and identity checks, explaining each material source difference. Stop before publishing weakened standing results.
- [x] Check frozen confidence/motif/structure outputs and current IDs remain bound to the selected cohort. Do not regenerate accepted scientific outputs solely for packaging.

**Acceptance:** Cohort/source provenance is reconciled; required stage reports and identity invariants pass; no unexplained record loss, duplicate identity or mixed build stages.

## 3. Manual visual and table proofreading

- [x] Run existing summary structural/style check: 9 headings, 5 tables, 8 images match baseline.
- [x] Extract and visually inspect all eight current static panels on white backgrounds. Initial findings: crowded/overlapping flow labels, misspelled “utoimmune” in source prose, COVID source predicate includes SARS-CoV-1 as well as SARS-CoV-2. These findings were corrected and the changed outputs re-inspected.
- [x] Audit tables for duplicate keys, accidental abbreviations, typos, missing labels, inconsistent units/denominators and agreement with selected source counts. Check COVID specificity and aggregation explicitly; do not change scientific thresholds.
- [x] Inspect splash light/dark variants and interactive summary, legends, tooltips, ordering, labels, clipping, all table sections and any release infographic. Distinguish repeated legitimate observations from duplicate curation.
- [x] Trace each defect to owning code/callers. Correct only release-facing defects; use existing plotting backend and installed tools. For nontrivial logic add a focused failing regression check before the fix.
- [x] Re-render affected summary outputs from the same selected inputs, rerun automated checks and visually inspect all changed panels/tables at readable scale. Keep the review record with findings, fixes and verification evidence.

**Acceptance:** No unresolved release-blocking visual/table defect; every changed display is manually inspected; final images and tables belong to the shipped candidate.

## 4. Build immutable candidate bundles

- [x] Make the minimum required tag-policy compatibility changes in release/bundle.py and owning workflow/help/docs. Date and date-NOTE tags now validate real calendar dates and safe suffix characters; CalVer remains accepted. Added exact tag/URL and immutability regression checks. Observed three new date-tag cases fail before the change; all 28 release/finalization tests pass afterward (2.34 seconds).
- [x] Run focused release/manifest/finalization/summary tests with bounded timeouts. Preserve existing supported consumer behavior.
- [x] Package primary, legacy, AIRR, web and paired web-structures when valid accepted preparation is available, using existing vdjdb release APIs. Include manifest.json, SHA256SUMS, lifecycle asset and concise release notes.
- [x] Verify every bundle/member checksum, member uniqueness, safe names, CRCs, gzip thresholds, embedded download pointer and manifest/source/tag consistency. Confirm web table, metrics and structures remain sealed to one cohort.
- [x] Measure disk footprint before copying full structure ZIPs; avoid duplicate archives and retain only requested final outputs and compact receipts.

**Acceptance:** All packaged bytes are verified and frozen for the Aldan-3 checks; required roles are present and cannot be confused by consumers.

## 5. Aldan-3 interoperability against actual bundles

- [x] Verify live access through aldan3 doctor; SSH and existing authentication work outside the local network sandbox.
- [x] Inventory existing project paths, runtime versions, consumer source revisions, Java/Python environments and scheduler capacity; select isolated test destinations and resources.
- [x] Preview then transfer only required source/environment recipes, candidate archives and public test inputs; compare destination hashes. Keep runtime/auth configuration separate.
- [x] Submit bounded SLURM checks with discoverable logs and exit status. Exercise Python vdjmatch exact-tag/default/primary/legacy/AIRR acquisition and actual receptor annotation using existing public fixtures; validate row counts, full junctions, pairing and score/output semantics. Check historical and date-NOTE discovery compatibility separately where affected.
- [x] Install candidate web roles through tools/install_release.py; verify loader counts, metadata, motif tables/charts and structure bindings. Run actual Scala loader/backend and annotation paths in the supported legacy runtime; include Browse/overview/motif/structure/export journeys against an isolated preview when feasible. Do not replace production configuration.
- [x] Exercise legacy JVM vdjmatch dependency as used by vdjdb-web, not only the Python namesake. Use one batch of query receptors rather than per-record engine calls.
- [x] Retrieve compact receipts, versions, exact archive digests and failure logs. Investigate any stalled/failed job and rerun only the affected check after correction. Record actual execution boundaries.

**Acceptance:** Both consumers actually run on scheduled Aldan-3 nodes using the final candidate bytes; expected counts, annotation fixture results and web journeys pass. No substitution of local or Stratos evidence.

## 6. Integrate and publish

- [x] Review final changes for correctness/confidentiality, integrate through the repository's supported branch/check route and verify the visible integration checkout. Do not publish a tag that falsely claims to represent a different source build.
- [x] Confirm tag/title/date, previous-release comparison, role descriptions, counts and known caveats in release notes. Include concise usage instructions for new-format and legacy/web consumers.
- [x] Publish the immutable GitHub release titled exactly “Late Autumn Release”, following existing release protections and required checks. Respect GitHub environment gates; publication approval must be explicit.
- [x] Read the published release back, verify tag points to pinned SHA, asset names/sizes and downloaded manifest/checksums; smoke-test real published acquisition on Aldan-3.
- [x] Finalize latest-version.txt only after publication succeeds, preserve concurrent source commits and verify line 1 returns the intended legacy asset successfully.
- [x] Deliver release URL, exact tag, asset list, compact Aldan-3 and visual-review evidence, and any remaining boundary. Mark the goal complete only after all required work is done.

**Acceptance:** Public release and final pointer resolve to verified assets; title/date/source match; final consumer checks use published assets.

## Current decision and evidence ledger

- Published [Late Autumn Release](https://github.com/antigenomics/vdjdb-db/releases/tag/2026-11-06), tag **2026-11-06**. Retain ordinary YYYY-MM-DD tags and safe optional NOTE suffixes; no policy migration or Zenodo twin. Exact packaging source: **8c4cdba9a3bcc0f6e286f2e58cb88ebd09a73bcc**, verified through the public annotated tag and manifest.
- The author explicitly selected the accepted rich snapshot and its matching prepared structures. Assembly source **5ffa45de203e1660dbfa797e4dd5e73d231362e1** and frozen confidence-v3.2 scoring source **69b326a10f5fab1c7b265dc7e97f790752c28aeb** are recorded in public BUILD_PROVENANCE.json. Thirty scientific inputs remain byte-identical; summary/date metadata were updated. No accepted motif, score or structure analysis was regenerated for packaging.
- The current-source CI cohort was compared and deliberately excluded: documented shared-locus dependency changes affect 3 junctions, 698 V-call cells, 2 legacy score cells and 1,000 hash cells. Cohort selection is the author's explicit decision. Both have the same 605 summary rows; this was never treated as raw table equivalence.
- Published assets: primary, legacy, AIRR, web and web-structures ZIPs; manifest, SHA256SUMS, BUILD_PROVENANCE.json, RELEASE_NOTES.md and compressed identity lifecycle. Every GitHub-computed asset digest equals the verified final local bytes. ZIP boundaries/CRCs, safe unique members, deterministic gzip and the 1,000,000-byte table limit passed. All structure downloads were checked against sealed members.
- Counts: **315,921 records**, **516,101 chain observations**, **1,625,234 evidence rows**, **307,789 slim records**, **200,180 paired receptor groups**. Web has 27 columns; legacy has 22. Prepared structures preserve **8,346 hashes / 50,086 sealed files**; web metrics contain 11,234 hashes and cover every served hash.
- Manual review covered all eight original panels, both splash themes and all 607 original table rows; all 51 changed/new rows were read after correction. Final tables contain 605 rows (7/87/226/203/82), with no exact duplicate rows or malformed widths. Changed plots were re-inspected at readable/native scale. Offline content proofs and direct image inspection are the evidence; no browser-CSS inspection is claimed.
- Corrected full epitope flow labels/overlap, SARS-CoV-2 selection, autoimmune typo, HLA stars/full beta-chain names, chart captions, chord spacing and distinct publication count (845). Source aliases remain explicit. Joint MHC/V/reference panels now count observed contexts rather than Cartesian products of independently collapsed lists: 13 unsupported heatmap cells removed, 15 shared counts corrected; all 307,789 slim identities retained, with exact R/Python parity across 2,867 final cells. Scientific thresholds and model outputs were preserved.
- Focused source checks: 94 passed; Ruff and strict mypy pass. Full initial CI passed 1,629 tests with 22 skips. Final required PR build [42122148230](https://github.com/antigenomics/vdjdb-db/actions/runs/42122148230) succeeded. Database PR1416/1417 integrated; master release merge **173e914674789908b78583db243ea44d841e5b78** contains the exact packaging source.
- JVM gzip/versioned transport and ancillary asset recovery were integrated as immutable **1.3.3** with JitPack availability verified. Python gzip release-loader work was integrated and **1.2.3** published to GitHub and PyPI; the PyPI wheel is byte-identical to the actual Aldan-tested wheel. Web PR212 merged as **561e7a0fb3d6ab13e790b164f0c2973392fe0adf**. PR creation remains enabled for **collaborators only**. Integration and protected-environment publication approvals are complete.
- Aldan-3 uses isolated tracked project storage, pinned Java8/Node11 runtime, scheduled eight-CPU nodes and bounded jobs. Web backend suite: 339 passed across 43 suites, one pre-existing cancellation, using the actual published JVM dependency. Final archive checks: job1509698 verifies exact default/primary/legacy identities, shared AIRR identities, stable record IDs, 200,180 pair groups and six/six batch annotations. Job1509699 verifies JVM default/full/slim readers and six/six batch queries, then actual web readers, metadata, suggestions, summary, both motif families and every served structure.
- AIRR Rearrangement/Reactivity omit receptor species on all 516,101 rows. Every exported junction, V/J, pMHC, antigen and reference field matches by stable cell/record+locus identity; pair membership is preserved. Release notes direct species-filtered consumers to primary/legacy. No species was guessed or manufactured.
- Public Python acquisition passed on Aldan (job1509700): real GitHub discovery/downloads and digests for default/primary/legacy/AIRR, all identity/pair checks and six/six batch annotations. Public web acquisition/installation passed as job1509709 (26.46 seconds). Public JVM default-pointer acquisition passed as job1509708 (54.05 seconds), loading all 516,101 full / 307,789 slim rows and matching six/six batch queries. Both used the public GitHub release; no discovery stub was used.
- latest-version.txt finalized only after publication and public legacy asset HTTP206 verification, through the existing concurrent-master-safe helper. Pointer commit **cd451bc** is integrated in the visible checkout; public line1 names the dated legacy ZIP. The immutable release tag still points to 8c4cdba. No production service configuration was changed.
