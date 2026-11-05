# Legacy web audit and release handoff

This audit keeps the legacy application functional with summer and current
co-deliverables. It does not replace the application or change the accepted rich
scientific outputs. The execution checklist is in
[legacy-web-interop-plan.md](legacy-web-interop-plan.md).

## Counts and content identity

The summer release contains **284,546 rows**; the public production export
contains **228,214 rows**. These are different database contents, not a loader
or default-filter discrepancy. Their exact net difference is **56,332 rows**:
64,133 summer-only identities minus 7,801 production-only identities, with
220,413 shared identities (multiplicity retained). Production contains all
226,494 December-2025 identities plus 1,720 additions. The complete trace and
source checks are in [legacy-web-record-counts.md](legacy-web-record-counts.md).

The current web co-deliverable loads **516,101 rows**. Its default query retains
508,168; MHC II retains 28,646; MHC II plus native evidence retains 176. These
three filtered counts were checked against the source predicates. Summer loads
284,546 and its default query retains 275,804. Totals and filtered counts remain
separate throughout.

## Functional evidence

| Journey | Observed result |
| --- | --- |
| Navigation and overview | Main destinations load; eight overview images render without desktop horizontal overflow. Broken field/submission documentation links were replaced with verified targets. |
| Browse input and recovery | Exact public sequence queries, invalid patterns, no-hit queries and subsequent recovery work. CAVSESPFGNEKLTF gives four alpha rows; appending partners gives eight alpha/beta rows. |
| Browse filters | Antigen source and peptide autocomplete, MHC class, native evidence and species controls work. Cattle plus unmapped rows exposes all 52 source cattle rows; deselecting every species gives no results. |
| Paging and sorting | All 25/50/100 page sizes work. The 35,636-row default-filtered NLV query ends at page 357 with 36 rows at size 100. Boundary navigation and both CDR3 sort directions work; pagination supports Enter. |
| TSV export | The four-row paired fixture exports four rows with pairing disabled and eight with it enabled. Header, sequences, row counts and response digests were inspected. Completed delivery removes temporary files by design. |
| Paratopes | Specific alpha/beta motifs retain source sequences, V/J calls, cluster sizes and method identity. For GIL/NLV/GLC/PKY, served values match 18,267 current TCRNet PWM cells, 34,132 current TCREMP cells and 15,317 summer TCRNet cells. Summer supplies no TCREMP tables and degrades to an empty optional-method tree. |
| Motif exploration | Tree/scatter selection, sequence search and background subtraction were exercised. Current member TSV exports match all 17 source columns in 47 TCRNet rows and five TCREMP rows. Summer/current NLV alpha IDs differ correctly with their underlying clusters (25 versus 24 members). |
| Structures | Current NLV shows 3,340 positive-contact cards; summer shows 3,343. Five-map overlays, removal/restoration, front promotion, zoom/reset, hover, linked motifs and Browse destinations were inspected. |
| Native structures | Zero-contact native maps are excluded. Retained native structures show ipTM 1.00 (100%). PKY maps have the expected paired sequences and 14/15 peptide contacts; predicted ipTM comes from metadata. |
| Structure downloads | Structure, contacts, C-alpha and bundle options were exercised. All served per-hash bundle members, CRCs and contact-map/download agreement were checked across the accepted asset sets. |
| Annotation | The existing demo account runs a public three-clonotype synthetic fixture (175 synthetic reads), returns all three expected NLV/GIL/GLC clonotypes, exposes match details and renders summary/pie options. Distance labels now select their actual modes. |

## Defects corrected

- Native labels were not associated with checkbox/radio controls. One shared
  association restores label activation and preserves native disabled behavior.
- Cattle records could not be selected explicitly in the species form.
- Pagination anchors lacked keyboard activation; native buttons retain the
  existing page-selection logic and visual state.
- Busy search actions could be triggered while another search ran. Their guard
  now follows actual request activity and permits retry after an initial error.
- A displaced `@Singleton` created multiple user providers and seeded the same
  demo file four times at startup. The provider is again a singleton; a fresh
  isolated runtime seeds exactly one sample. Earlier test state was preserved.
- Quick-view popup cleanup could leave an attached body element after hiding or
  navigation. Explicit hide, disabling and destruction share one cleanup path.

## Versioned producer/consumer contract

New dataset versions use immutable `vYYYY.MM.PATCH` tags, starting each month at
patch zero. Publication is pinned to a source SHA; a concurrent source update
cannot silently change the release or prevent safe pointer finalization.

The complete local candidate contains primary, legacy, AIRR, web and
web-structures roles. All manifest and archive checksums pass. All 32 accepted
input files (254,795,745 bytes) remain byte-identical. Rich tables retain 315,921
records, 516,101 chains and 1,625,234 evidence rows. No motifs, scores or maps
were recomputed for this packaging check.

VDJmatch exact-tag fetching was exercised with the actual candidate artifacts
through a loopback fixture, mocking only GitHub release discovery. Default and
explicit primary, legacy and AIRR each load 516,101 rows. The legacy result keeps
its `.txt.gz` filename; default and primary share an identity; a wrong exact tag
fails before download. Web roles are never substituted for Python database roles.

The web installer installs an immutable version directory after checking the
manifest, member sets, checksums, table schema and prepared-asset bindings.
The actual full current install contains 516,101 table rows, 11,234 metric rows,
8,346 structure hashes and 50,086 sealed files. The running preview selects that
installed root for database, charts and structures together. The original summer
ZIP also installs successfully in explicit table-only mode, preserving decoded
table bytes; historical-only staging repairs its old metadata descriptors.

Use the database release instructions in [builds.md](builds.md), the VDJmatch
`update --pin TAG --asset ROLE` command, and the legacy web
`tools/install_release.py --tag TAG --destination DIRECTORY` runbook. Normal web
installation requires both web roles; table-only mode is explicit. Service
restart/recreation follows selection of a validated installed root.

## Measured response times

Final installed-release HTTP measurements use three requests per journey on the
same Stratos preview; counts and response sizes were identical across trials.

| NLV journey | Returned items | Response bytes | Median seconds |
| --- | ---: | ---: | ---: |
| Browse, no default filters | 35,947 records | 28,255 | 0.3071 |
| Structures | 3,340 cards | 2,928,502 | 0.8509 |
| TCRNet motifs | 240 clusters | 1,000,152 | 0.1739 |
| TCREMP motifs | 584 clusters | 2,454,736 | 0.3570 |

The earlier production Browse median was 0.8213 seconds for 34,364 records.
Different hosts, networks and contents prevent interpreting this as a controlled
speedup. The final fixed preview shows no material responsiveness regression in
the exercised journeys.

## Verification boundaries

The full backend suite passes 339 tests; one pre-existing external-input control
benchmark is explicitly skipped. The frontend suite passes 122 tests, including actual templates and
native label, paging, export and popup interactions. Release tooling and consumer checks remain offline except
for the explicitly identified original release/export inputs.

The demo account intentionally has no upload allowance. Registered/token upload,
account management and email workflows were not manually exercised without an
existing authorized session; format conversion has dedicated backend coverage.
The browser viewport override returned a 1280px viewport after requesting 375px,
so mobile visual conformance is unverified. Production filesystem provenance is
unverified because management access was unavailable; the complete public export
establishes its content identity. No production deployment or release publication
was performed.
