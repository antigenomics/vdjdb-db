# Legacy web interoperability audit

## Scope and acceptance plan

Deliver a separate legacy-web projection from the rich build, with surgical web
adapters. Preserve ordinary legacy exports, rich tables, curation, scoring and
motif analysis. A web rewrite belongs to a separate project.

1. Trace file/schema consumers → test plain/gzip loading and the original 22 fields.
2. Join shipped evidence → verify flags and selected-hash metric provenance.
3. Review ambiguous links → retain paired junction, peptide and cited-PDB identities.
4. Prepare maps through TCRen → check canonical coordinates, shared overlay scale,
   residue indexing, metadata interface ipTM and served download filenames.
5. Host a frozen source package on Stratos → manually exercise Browse, evidence,
   both motif methods, structures, overlays, hover and downloads.
6. Compare large-epitope responsiveness with the live site → report measured scope,
   cohort sizes and host differences; commit only verified compatibility changes.

## Comprehensive UI benchmark extension

Goal: explain the summer/production count discrepancy and exercise all available
legacy UI workflows on Stratos with both summer and current web co-deliverables.
Use the writing-plans, test-execution, systematic-debugging and executing-plans
skills, plus the browser's local-web-development guidance for manual inspection.
The author has explicitly authorized execution and functional fixes.

Execution remains in the visible `codex/legacy-*` checkouts. Preserve other edits,
rich-build outputs, release locks, accepted motif results and the published TCRen
potential. Production is a read-only baseline; reversible uploads and exports
belong on Stratos. Keep raw public exports and identifying mappings under ignored
`out/interop/`; publish aggregate findings and neutral reproductions only.

### Review focus

- Total versus filtered versus paired counts must remain distinct and trace to source rows.
- Filters must behave identically after paging, sorting, reset, refresh and deep links.
- Summer lacks optional TCREMP/evidence/chart inputs; absence must not break other pages.
- Empty/invalid searches, missing structures and stale selected IDs must recover visibly.
- Large NLV card lists, overlays and mobile layouts must remain usable without errors.

### Task 1: Trace database counts

Owning code: web `Database.scala`, `DatabaseMetadata.scala`, `SearchTable.scala`,
`filters_meta/meta-filters.ts`, and the deployed/released table provenance.

- [x] Verify current live build identity and raw metadata count; independently query zero-filter and default-filter counts.
- [x] Verify summer ZIP checksum, source row count and loaded count without changing its tables.
- [x] Compare live and summer row multisets using available deployed inputs or a public export; group every added/removed row by source and transformation.
- [x] Test the specific suspected cause against exact row counts; distinguish data selection, deduplication, parsing and display filtering.
- [x] Record a reproducible arithmetic reconciliation, evidence and any unverified boundary.

Count reconciliation: [legacy-web-record-counts.md](legacy-web-record-counts.md).

### Task 2: Navigation, overview and visual baseline

Owning files: existing navbar, home, overview and shared responsive components.

- [ ] Exercise every navbar/menu destination, internal cross-page link, browser back/forward and reload.
- [ ] Inspect Home/About/documentation links and Overview figures for loading, clipping and missing assets.
- [ ] Inspect desktop and 375px mobile layouts, keyboard focus, dropdowns, tooltips and cookie controls.
- [ ] Record console errors, broken resources, clipped controls and misleading loading/empty states.

### Task 3: Browse and evidence

Owning files: existing search components, filters, `SearchTable`, `DatabaseAPI` and websocket actor.

- [x] Verify totals, default filters and named epitope counts for both databases against source predicates.
- [ ] Manually exercise CDR3 exact/substring/pattern and chain; antigen sequence/list/origin/gene; MHC class/allele; metadata/confidence/evidence filters, including combinations and negatives.
- [x] Verify reset/refresh, every page-size option, first/previous/next/last, ascending/descending sort and paired-row expansion.
- [ ] Follow M/S evidence badges and external references; verify sequence, gene and epitope identity at destinations.
- [x] Download TSV with paired export off/on; verify header, row count, sequence contents and pairing behavior against the selected records.
- [x] Exercise empty results and invalid patterns; verify visible recovery and a subsequent successful query.

### Task 4: Paratopes

Owning files: existing motif service/components, SVG logos, scatter bridge and method adapters.

- [x] Inspect specific alpha/beta GIL/NLV/GLC motifs and their sequences against source members/PWMs in both releases.
- [ ] Exercise all supplied methods, epitope tree, multiple selection, collapse/clear and background subtraction.
- [ ] Exercise CDR3 exact/substring search, chain selection, rankings and no-hit recovery.
- [ ] Exercise scatter/legend selection and deep-link/reload synchronization with the correct logo.
- [ ] Download motif/member exports and inspect sequences, cluster IDs, counts and filenames.
- [x] Verify summer's missing optional method/charts do not create errors or mislabel another method's results.

### Task 5: Structures

Owning files: existing structure service/cards, SVG index, hover controller and download adapter.

- [ ] Verify MHC/epitope selectors and CDR3 chain/substring search, including no-hit recovery.
- [x] Inspect predicted/native paired sequences, metadata ipTM, native 100%, positive peptide contacts and exclusion of zero-contact maps.
- [ ] Exercise selection 0→5, deselection, sixth-selection handling, front promotion, zoom/pan/reset, hover and selected-ID reloads.
- [ ] Follow alpha/beta motif and Show table links; inspect destination sequence and evidence filters.
- [x] Download PDB, contact table, CA coordinates and bundle; verify names, types, CRC, row identities and contact-map agreement.

### Task 6: Annotation and auxiliary workflows

Owning files: existing guest annotation/upload/sample/summary/export components.

- [ ] On Stratos, exercise existing public synthetic/demo inputs, supported format selection, upload and resulting sample navigation.
- [ ] Exercise available sample settings, annotation launch/status/results, summaries and exports with explicit bounded runtimes.
- [ ] Test user-visible invalid input recovery with a small neutral fixture; do not create accounts, send mail or change production state.
- [ ] Record any authentication-only controls that cannot be exercised with existing authorized sessions.

### Task 7: Measured regression checks, fixes and delivery

- [x] Benchmark fixed Browse, motif and NLV structure journeys with exact cohorts, elapsed times, response sizes, DOM counts and server memory; use three bounded requests per comparable API case.
- [x] For every reproducible defect: record steps/expected/actual; add the smallest failing check; fix the owning shared cause; rerun that check and its manual journey.
- [x] Run affected Scala/frontend/Python checks first, then required full suites/build once for the final code state; cap ordinary tests at 5–10 minutes and diagnose stalls before rerunning.
- [x] Conduct one independent final source review, reconcile important findings, and recheck corrected paths.
- [x] Update aggregate audit/results and limitations, commit changes, preserve the preview, and checksum-sync final sources to Stratos.

### Task 8: Versioned producer/consumer handoff (author extension)

Ownership: `vdjdb-db/src/vdjdb/release`, release CLI/workflows and tests;
`vdjmatch/src/vdjmatch/db/vdjdb.py` and fetch tests; legacy web preparation tools,
new stdlib release installer, and focused offline adapter tests.

- [x] Strict canonical monthly tags `vYYYY.MM.PATCH`, starting at patch 0; historical exact date tags remain readable.
- [x] Pin publication to a source SHA, verify immutable tags, checksum the manifest, resolve legacy assets by role and resolve `--previous` to a real downloaded ZIP.
- [x] Preserve compressed physical table names in every release role; teach shared VDJmatch ingestion to accept gzip while retaining primary as its default role.
- [x] Bind prepared maps and metadata to the exact web table and source metadata using schema-versioned digest records. Existing geometry and accepted outputs stay unchanged.
- [x] Package optional complete structure assets under role `web-structures` alongside `web`; preserve default rich/legacy/AIRR release roles and latest-version's legacy URL.
- [x] Implement a stdlib web installer: exact-tag/manifest-role selection, checksums, safe fresh staging, full validation, and immutable release-directory installation. Never overlay old plain/gzip files or silently replace existing assets in place.
- [x] Require complete structure assets for normal installation; make table-only mode explicit for releases intentionally lacking maps. Keep charts under the same selected database root; restart/recreate after selecting a release.
- [x] Validate prepared-asset consistency and failed-update preservation with offline fixtures, then exercise the actual prepared Stratos co-deliverable without recomputing motif/structure analyses.
- [x] Document exact fetch/update commands and source/version identities, review cross-repo changes, commit and final-sync changed sources.

Shared prepared-asset contract: schema version 1. `structure/asset-manifest.json`
records `web_table_sha256`, `source_metadata_sha256`, structure hashes and served
file paths/bytes/SHA256. `prepared-assets.json` beside the web table records those
two input digests plus `structures_metadata_sha256` and
`structure_manifest_sha256`. The web archive includes this record; the separate
structure archive includes `structure/asset-manifest.json`. Validate all digest
bindings before packaging or installation. Hash actual physical gzip bytes;
never infer asset identity from a directory name, tag spelling, or existing hash.

Report each journey as pass/fixed/blocked/not applicable with exact observations.
Previously verified evidence is retained; this extension does not authorize motif
reanalysis, production deployment or a web rewrite.

## Design

`vdjdb make web` emits 27 columns; the ordinary legacy projection keeps 22.
`vdjdb release --with-web` opts into an additional ZIP without changing existing
bundle roles. Generated large tables remain deterministic gzip, including in the
web ZIP. The web resolves logical table names as plain/gzip streams.

The web-only same-study flag means nonempty `method.verification`, matching the
live site's effective predicate. Independent and native evidence come from the
rich view. Model contacts/quality use evidence for the selected source hash.
Native links prefer their cited PDB. The rich evidence producer is unchanged.

The reviewed model table resolves all 54 ambiguous records (108 associations,
18 models). All preserve both full junctions and GILGFVFTL. Current allele
agreement resolves 52; original source V calls resolve two ties. None has a PDB
ID. The separate native table reviews 19 hash groups with 47 candidate PDB IDs,
using reported experimental resolution and a lexical tie-break. Unknown or
changed ambiguities fail instead of silently choosing a model.

Motif charts use shipped embedding coordinates and PWMs: no re-clustering or
re-embedding. There are 456 HTML charts with one shared Plotly asset. Trace names
remain cluster IDs for scatter/PWM synchronization.

TCRen owns structure parsing, batched receptor/MHC annotation, canonical projection,
contact geometry and SVG rendering. A small SVG adapter preserves legacy hover
and overlay behavior, including PDB insertion codes. All maps share coordinate
bounds. Exact metadata junctions disambiguate legacy display chain roles before
complex selection. Predicted interface ipTM is read from `tcr_pmhc_iptm`. By the requested
experimental-display convention, native ipTM is 1 (100%); native model confidence
and model percentiles are absent. Native and predicted structures both require
positive CDR3-peptide contacts for display. The companion contact count is derived
from the actual TCRen map; predicted evidence flags and model scores retain their
source meaning. Contact downloads use the same TCRen tables as rendering. This
corrects 253 native download/map discrepancies in the old supplied contact tables.
Observed zero values remain observations.

## Verified acceptance (2026-11-05)

- Web projection: 516,101 chain rows; 17,866 same-study flags before default
  noncanonical/unmapped filtering. Hosted Browse: 508,168 default matching rows;
  same-study filter: 17,292. The live release has a different cohort.
- All 8,346 referenced metadata hashes have raw-coordinate assets and full/simplified
  maps: 16,692 HTML maps and 8,346 CRC-checked download bundles. Full/simplified
  residue positions agree for every structure. Metadata has 11,234 unique hashes.
- The hosted API exposes exactly 7,855 positive-contact structures: 7,551 predicted
  and 304 native. The 481 zero-contact predicted and 10 zero-contact native
  structures are excluded. Every displayed predicted ipTM equals the companion
  metadata; all displayed native ipTM values equal 1, with no model confidence.
- Database schema/new-format/release checks: 71 passed, 3 skipped.
- Final Stratos backend: 338 passed, one canceled; production package built.
  The focused native-policy suites passed 27 tests.
- Captured-browser frontend suite: 115 passed. Python adapter checks: 10 passed.
- TCRen shared-bounds/project-table checks: 13 passed.
- Real-coordinate rendering check, 64 structures: serial 2.368 s; 12 workers
  1.620 s. Annotation ran once for that check, outside the rendering timings.
- Live manual baseline: Browse 223,344 matching / 228,214 total; motif deep-link
  GILGFVFTL opened cluster H.B.GILGFVFTL.1 with scatter and logo; NLVPMVATV
  structure view loaded 3,355 cards, five selected, metrics and overlay controls.
- Hosted manual checks: Browse, same-study evidence filtering, Browse-to-motif
  link, TCRNet scatter/logo and TCREMP scatter/logo load without browser errors.

## Manual audit and performance

Hosted native MHC-II PKYVKQNTLKLAT loads two cards and overlaid maps, positive
contact counts, native labels and 100% ipTM. Residue hover reports author numbering
(`PRO 306`); contact hover reports both endpoints (`LYS 307 : GLU 94`). Zoom,
front promotion and ZIP download work; all 8,346 download bundles have CRC-valid members whose bytes match the served files. Show table opens Browse with eight matching chain rows.

The hosted NLVPMVATV view renders 3,340 cards and five selected overlay maps,
with no captured browser errors. Its Alpha motif link opens TCRNet cluster
H.A.NLVPMVATV.30, size 24, with scatter and logo. Browse-to-motif deep links and
TCREMP scatter/legend-to-logo selection were also exercised.

Identical bounded API requests were sent three times to each host. Median elapsed
seconds (including transport):

| Request | Live vdjdb.com | Stratos preview | Observations |
|---|---:|---:|---|
| NLV Browse, first 25 rows | 0.8213 | 0.3353 | 34,364 live / 35,947 new matching chain rows |
| NLV structures | 2.2068 | 0.8648 | 3,355 structures on each host |

The structure response was 2,944,668 bytes live and 2,941,773 bytes on Stratos.
Those timings precede the final geometry-based exclusion of 15 zero-contact
NLV maps. The final raw-coordinate build took 360.93 seconds with one batched search per supported
reference and the existing 12-worker rendering API. One process-tree rendering
sample measured 26,578.6 MiB RSS and 51 threads; the host has 128 GB RAM.

## Summer and new-release loading matrix

The latest published asset remains `2026-06-03-ZENODO/vdjdb-2026-06-03.zip`,
41,017,409 bytes; SHA256
`9ae7724ab2ce9fdb09b1fb82a90c7462f75269f14f669238f8c154a80239cfc1`.
The identical production web binary is tested with separately mounted summer and
new-release databases and isolated runtime state. Summer tables remain unchanged;
prepared maps and their optional companion metadata are separate additions.

All counts below are chain records. Defaults exclude `web.cdr3fix.nc=yes` and
`web.cdr3fix.unmp=yes`; epitope filters add exact sequence equality. Both the HTTP
API and visible Browse counts were checked against the original table predicates.

| Browse quantity | Summer release | New legacy web co-deliverable |
|---|---:|---:|
| Total records | 284,546 | 516,101 |
| Default filtered records | 275,804 | 508,168 |
| NLVPMVATV, default filters | 34,547 | 35,636 |
| GILGFVFTL, default filters | 21,574 | 33,224 |
| GLCTLVAML, default filters | 13,555 | 14,893 |
| PKYVKQNTLKLAT, default filters | 549 | 600 |

Visually inspected influenza paratope H.B.GILGFVFTL.1 in both databases:
1,054 summer / 1,118 new members, 13 positions, TRBV19*01/TRBJ2-7*01,
with the CASS…EQYF pattern and matching residue ordering. The new alpha NLV
H.A.NLVPMVATV.30 has 24 members, 12 positions and TRAV12-2*01/TRAJ15*01;
its scatter and logo load together. Its paired structure card shows
CAENQAGTALIF/CASSAQGVVTEAFF, 15 contacts and predicted ipTM 0.74.
Selection removal/restoration changes five overlays to four and back, with
162 → 131 → 162 visible residue nodes and no browser errors.

For the four named epitopes, every source PWM count/frequency was compared to
the served API: the new release passes 18,267 cells across 691 TCRNet clusters
and 34,132 cells across 1,220 TCREMP clusters. No motif analysis was rerun. The summer release passes 15,317 PWM cells across
645 TCRNet clusters. It supplies no TCREMP tables; selecting that optional method
shows an empty tree without a browser error.
The summer structure-to-alpha link correctly targets H.A.NLVPMVATV.42 (25 members),
where the new build targets H.A.NLVPMVATV.30 (24 members); both preserve the
12-position motif sequence and TRAV12-2*01/TRAJ15*01 identity. Summer NLV displays
3,343 model cards with five selected maps, positive contacts and source ipTM.

Summer structure preparation independently covers 11,092 referenced raw structures:
22,184 HTML maps and 11,092 checked ZIPs. Every contact download equals its map's
edges and every ZIP member equals its served file. The runtime exposes exactly
10,425 positive-contact structures, including 203 native; 667 zero-contact maps
are excluded. All receptor/peptide identities pass, with 38 metadata-confirmed
chain role corrections. This fresh build took 462.99 seconds. Neither database's
maps, annotations or projections were reused as computed intermediates.

Native PKY cards show CAVSESPFGNEKLTF/CASSSTGLPYGYTF (14 contacts) and
CAVSEQDDKIIF/CATSDESYGYTF (15 contacts), both ipTM 100%. Native contact hover
reports author-numbered endpoints (for example GLU 102 : THR 97).

## Reviewed deposition differences

All 17 native identity discrepancies were inspected. Eight are exact peptide
segments in fused constructs, recovered from raw observed coordinates. Six have
alpha-junction differences, two lack a metadata terminal peptide residue, and one
contains modified residue F2F (displayed as X). The nine latter cases are gated by
`vdjdb-web/tools/native_structure_variants.tsv`; their maps and download bundles
carry explicit expected/observed labels. No residue is invented or coordinate
changed. Unknown mismatches fail together for review. Eight chain roles are
restored from exact metadata junctions before selecting the primary complex.

## Current tool compatibility and delivery

An isolated Stratos environment passed 84 motif/export tests with seqtree 1.3.3,
arda-mapper 2.38.3, vdjtools 5.9.1, vdjmatch 1.2.1 and mirpy-lib 4.10.1. Map
preparation used TCRen 4.4.2 plus the tested shared-bounds addition, ARDA 2.38.3
and seqtree 1.3.3. The existing rich-build lock, scored tables, motif coordinates
and PWMs remain unchanged. All 22 ordinary legacy fields match in 516,101 rows
except 2,320 reviewed web-only structure-hash cells.

The web co-deliverable is `out/web-release/vdjdb-web-interop-audit.zip`:
57,001,726 bytes, 467 members. Maps/downloads are separate assets on Stratos.
No generated plain TSV/TXT exceeds 1,000,000 bytes. The runbook is
`vdjdb-web/tools/README.interop.md`. Private input hashes, recipes, API timing and
source-sync receipts are retained under `out/interop/`; source syncing uses
checksums, backups and preservation of destination-only files.

## Limits

The timing comparison traverses different hosts/networks and the Browse cohorts
differ, so it characterizes this audited setup rather than a controlled same-host
benchmark. The test host uses isolated runtime state and loopback access. No
production deployment, release or push is part of this delivery.

## Versioned handoff validation (2026-11-05)

The complete local `v2026.11.0` candidate contains all five roles: primary,
legacy, AIRR, web and web-structures. All six SHA256SUMS entries (including the
manifest) pass. Packaging preserves all 32 accepted source files (254,795,745
bytes) exactly; rich tables retain 315,921 records, 516,101 chains and 1,625,234
evidence rows. Legacy exports retain 516,101 main rows, 307,789 slim rows and
315,921 full rows. No scientific analysis was rerun.

The actual complete web installer passes with 516,101 rows, 11,234 metric rows,
8,346 structure hashes and 50,086 sealed files. The running isolated preview
selects this immutable installed root, with charts and structures under that same
root. Original archives and the prior container remain available for recovery.
The historical summer installer passes with 284,546 rows and the exact original
decoded table SHA256. Its historical metadata omits a hidden hash row, puts the
score row in another position and drops a trailing comment; historical-only
staging normalizes those descriptors without changing table bytes. Modern
manifest releases retain strict schema and byte validation.

Final synchronization verified all 4,215 selected source files across the eight
software repositories with zero differences or conflicts. Remote Git heads and
destination-only files were preserved. Publication and production service
changes remain outside this validation run.

Final UI suite: 124 passed; final backend suite: 339 passed, with the existing
external-input control benchmark explicitly skipped. The final binary again
loads summer's 284,546 total / 275,804 default-filtered rows, expands the reviewed
four-row alpha query to eight paired rows, and shows native sequence-search maps
with 14 contacts and ipTM 1.00 without console errors. Its absent TCREMP method
remains an empty optional tree. The current full installed root is restored after
this sequential check.

The detailed findings, measured counts, defect fixes and verification boundaries
are in [legacy-web-ui-audit.md](legacy-web-ui-audit.md). Unchecked broad checklist
items must not be interpreted as fully verified: the demo account forbids uploads,
mobile viewport resizing did not take effect, and some exhaustive control/route
combinations remain outside the recorded manual journeys. The audit reports
these boundaries instead of claiming every possible UI state is bug-free.

The popup lifecycle correction also passes its packaged manual journey: Quick
view is absent after hiding and after navigation to the structure page. Exact-tag
fetching of all actual Python database roles passes against the full candidate
artifact bytes; only GitHub discovery was replaced by the local test fixture.

A final shortcut check found 381 SARS-CoV-1 rows entering the COVID-19 preset
through its broad SARS-CoV substring. The shared preset now uses SARS-CoV-2.
Real-service collection tests preserve Influenza/self-antigen behavior; all 124
frontend tests pass. No database rows or scientific outputs were changed.
