# Summer archive versus production record counts

The difference is in the record contents, not the Browse default filters.
The public production export was fetched with no filters and paired export
explicitly disabled. Its 228,214 data rows equal both the HTTP search count and
the Browse websocket count. The summer archive contains 284,546 data rows and
loads exactly that many in the audited Stratos web instance.

| Chain-record count | Production | Summer archive |
|---|---:|---:|
| Loaded total / no filters | 228,214 | 284,546 |
| Default filters | 223,344 | 275,804 |
| Excluded by default filters | 4,870 | 8,742 |

Defaults exclude `web.cdr3fix.nc=yes` and `web.cdr3fix.unmp=yes`.
The total-count gap is 284,546 − 228,214 = **56,332**; filtering is a separate step.

## Exact content reconciliation

Matching multisets by chain, full junction, V/J, species, both MHC chains/class,
epitope, antigen gene/species and reference preserves repeated observations.
Complex IDs, regenerated structure hashes and added evidence columns are not
biological record identities and are excluded from this comparison.

Every one of the 226,494 identity rows in the `2025-12-29` reference is present
in production, including its multiplicity. Production adds exactly 1,720 rows:

- 1,540 rows from six studies: PMID:40892914 (817), PMID:40250246 (585),
  PMID:41610844 (120), PMID:41644530 (10), PMID:41424382 (4), and PMID:42129507 (4).
- 180 rows citing 90 RCSB structures.

Thus 226,494 + 1,540 + 180 = **228,214**. This identifies the deployed content as
the older reference plus additions; it does not establish when a server file was
copied or which preprocessing command was run.

Comparing the summer archive directly with production gives 220,413 shared
identity rows, 64,133 summer-only rows and 7,801 production-only rows:
64,133 − 7,801 = **56,332**. The counts by reference independently reconcile to
56,479 additional summer rows minus 147 additional production rows. This second
comparison groups reannotated receptor identities under their original reference.

Examples of the source differences:

| Reference | Summer rows | Production rows |
|---|---:|---:|
| PMID:32184241 | 23,475 | 0 |
| PMID:28636589 | 6,320 | 1,843 |
| PMID:30737354 | 2,761 | 0 |
| PMID:37133356 | 2,484 | 0 |
| PMID:36134660 | 2,414 | 0 |

The summer raw table has exactly 22 fields in all 284,546 rows. Its checksum is
`9ae7724ab2ce9fdb09b1fb82a90c7462f75269f14f669238f8c154a80239cfc1`.
The latest GitHub asset is named `vdjdb-2026-06-03.zip` under tag
`2026-06-03-ZENODO`; the observed asset was uploaded on 2026-08-26. A release
name alone therefore does not identify the deployed input.

## Loader and filter verification

`DatabaseMetadata.createFromInstance` reports the loaded row collection size.
Browse's `Total` uses that metadata; `Found` uses the search result size.
The live no-filter HTTP request, websocket request and full public export all
agree at 228,214. The unchanged summer input independently agrees with its
Stratos metadata and no-filter query at 284,546; no summer rows are dropped.
Default-filter queries produce the separate counts above. These checks explain
the observed gap without a loader or filtering change.

The production management SSH route was unavailable during this audit, so the
mounted production file was not inspected directly. The complete public export
and code trace establish the record-level difference; they do not prove the
server file's physical provenance. Large live metadata responses also stalled
on both local and Stratos network paths; the small count responses and full
static export succeeded. That transport issue is recorded separately from counts.

Private raw exports, source checksums, complete reference deltas and reproducible
checks remain under ignored `out/interop/`. See the comprehensive UI benchmark
in [legacy-web-interop-plan.md](legacy-web-interop-plan.md).
