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
positive CDR3-peptide contacts for display. Observed zero values remain observations.

## Verified acceptance (2026-11-05)

- Web projection: 516,101 chain rows; 17,866 same-study flags before default
  noncanonical/unmapped filtering. Hosted Browse: 508,168 default matching rows;
  same-study filter: 17,292. The live release has a different cohort.
- All 8,346 referenced metadata hashes have raw-coordinate assets and full/simplified
  maps: 16,692 HTML maps and 8,346 CRC-checked download bundles. Full/simplified
  residue positions agree for every structure. Metadata has 11,234 unique hashes.
- The hosted API exposes exactly 7,871 positive-contact structures: 7,570 predicted
  and 301 native. The 462 zero-contact predicted and 13 zero/missing-contact native
  structures are excluded. Every displayed predicted ipTM equals the companion
  metadata; all displayed native ipTM values equal 1, with no model confidence.
- Database schema/new-format/release checks: 71 passed, 3 skipped.
- Final Stratos backend: 338 passed, one canceled; production package built.
  The focused native-policy suites passed 27 tests.
- Captured-browser frontend suite: 115 passed. Python adapter checks: 7 passed.
- TCRen shared-bounds/project-table checks: 13 passed.
- Real-coordinate rendering check, 64 structures: serial 2.368 s; 12 workers
  1.620 s. Annotation ran once for that check, outside the rendering timings.
- Live manual baseline: Browse 223,344 matching / 228,214 total; motif deep-link
  GILGFVFTL opened cluster H.B.GILGFVFTL.1 with scatter and logo; NLVPMVATV
  structure view loaded 3,355 cards, five selected, metrics and overlay controls.
- Hosted manual checks: Browse, same-study evidence filtering, Browse-to-motif
  link, TCRNet scatter/logo and TCREMP scatter/logo load without browser errors.

## Manual audit and performance

Hosted native MHC-II PKYVKQNTLKLAT loads three cards and overlaid maps, positive
contact counts, native labels and 100% ipTM. Residue hover reports author numbering
(`PRO 306`); contact hover reports both endpoints (`LYS 307 : GLU 94`). Zoom,
front promotion and ZIP download work; the downloaded 343,493-byte bundle has
three CRC-valid members. Show table opens Browse with eight matching chain rows.

The hosted NLVPMVATV view renders 3,355 cards and five selected overlay maps,
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
The final package again returned 3,355 NLV structures in 0.8739 seconds. The
full raw-coordinate build took 361.23 seconds with one batched search per supported
reference and the existing 12-worker rendering API. One process-tree rendering
sample measured 26,578.6 MiB RSS and 51 threads; the host has 128 GB RAM.

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
57,002,006 bytes, 467 members. Maps/downloads are separate assets on Stratos.
No generated plain TSV/TXT exceeds 1,000,000 bytes. The runbook is
`vdjdb-web/tools/README.interop.md`. Private input hashes, recipes, API timing and
source-sync receipts are retained under `out/interop/`; source syncing uses
checksums, backups and preservation of destination-only files.

## Limits

The timing comparison traverses different hosts/networks and the Browse cohorts
differ, so it characterizes this audited setup rather than a controlled same-host
benchmark. The test host uses isolated runtime state and loopback access. No
production deployment, release or push is part of this delivery.
