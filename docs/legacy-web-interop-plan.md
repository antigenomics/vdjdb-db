# Legacy web interoperability audit and plan

Scope: preserve the rich database and its build stages. Adapt the existing legacy web's readers and supply compatible visualization assets; no frontend rewrite or scientific reanalysis.

1. Trace legacy table schemas, ZIP packaging and local loaders. Verify positional schemas and pairing against current build artifacts.
2. Add streaming gzip support at the web file boundary, shared by Browse, motif and structure readers. Check plain/gzip equivalence and existing web specs.
3. Trace evidence joins, motif badges/tree navigation and structure identifiers against emitted metadata. Preserve missing-value semantics and metadata-derived interface ipTM.
4. Use TCRen's existing structure/visualization APIs for new complementarity maps; retain old maps until verified replacements exist. Batch structure annotation; no reconstructed scoring formula.
5. Manually exercise current vdjdb.com Browse, evidence, Motif and Structure; repeat locally with the legacy co-deliverable, including the large NLVPMVATV case. Record actual checks and limitations below.
6. Run focused database/web checks, inspect rich-output regression coverage and generated table sizes. Commit only reviewed changes in the visible checkouts.

## Findings and verification

- Legacy table emission retains the positional 22-column database, 19-column membership and 27-column PWM contracts.
- Release packaging streams gzip sources into declared legacy ZIP member names. Local build directories instead retain `.gz`, which the current web readers do not discover. This affects Browse startup, motifs, structure rows and metrics.
- Structure pages require separately prepared `structure/<hash>.html` assets; absence prunes records from the structure tree. The legacy database ZIP alone does not contain those assets.
- Existing web metadata preparation explicitly takes interface ipTM from `tcr_pmhc_iptm`, not global `iptm`; TCRen does not calculate ipTM.
- Existing structure sync copies historical SVG maps from the dataset mirror. This does not establish that new maps are produced through current TCRen.

Status: audit and implementation in progress; no deployment or performance equivalence claimed.

## Implemented boundary

`vdjdb make web --tables out/tables --out out/web` emits a separate 27-column
web table plus matching metadata, preserving the original 22-column legacy files.
Evidence flags use shipped long evidence; model contact/quality flags preserve the
existing production predicates (>0 contacts; ranking confidence >=0.85). Same-study
validation follows the rich view's existing producer (currently false); this is an
open evidence gap, not evidence of absent biological validation.

The web readers now discover plain or gzip tables and stream decompression for
Browse, motif PWMs/members, structure table/header and per-model metrics. Plain
release ZIP members retain their existing names. No rich build stage was changed.

New maps can be generated with `vdjdb-web/tools/build_structure_maps.py` in the
TCRen environment using an explicit path/hash TSV, then served under `structure/`.
The renderer uses TCRen projection, contact tables and SVG APIs. TCR and MHC
annotation each use their existing batch API. No ipTM is computed from coordinates;
`build_structures_metadata.py` reads `tcr_pmhc_iptm` from model metadata.

## Checks so far

- 108 affected Scala web tests passed, including a new logical-name gzip load check.
- 103 database release/motif/structure unit tests passed.
- 53 schema/new-format tests passed, 3 skipped; web projection preserves all 22
  legacy columns exactly in the synthetic check, including zero-contact semantics.
- Live Browse displayed 223,344 matching records / 228,214 total and evidence controls.
- A live Browse motif link opened GILGFVFTL, selected cluster H.B.GILGFVFTL.1
  (896 sequences, junction length 13), with plot and logo present.
- Stratos Docker is reachable; native Java/sbt are absent. The author is syncing
  `~/hf` and `~/vcs` in another chat. Do not race that transfer or overwrite it.

## Remaining acceptance checks

Finish full-data web projection inspection, metadata and recovered structure-ID
joins, same-study evidence parity, gzip motif equivalence, TCRen renderer smoke,
Stratos hosting with fresh source identity, and manual Browse/evidence/motif/structure
checks there. Measure large-epitope responsiveness before claiming nonregression.
No production deployment or runtime performance equivalence has been established.

## Current audit update

All 54 ambiguous records were reviewed (108 candidate associations, 18 unique
models). Every selected model preserves both junctions and GILGFVFTL. In 52
records current V/J allele agreement resolves the choice. For VDJDB0000319328
and VDJDB0000319378, both candidates disagree with the repaired TRAV16 call;
original ISSUE_289 source rows 13247 and 13297 (one-based data rows) specify
TRAV12-2 and TRAV12-3 respectively. Those calls distinguish the representatives.
None of these records has a PDB ID. Choices and reasons are versioned in the
reviewed rules table; unknown ambiguous cases still fail.

Full-data projection: 516,101 chain rows; 2,006 structure-hash remaps, all other
legacy fields equal. Evidence counts are 145,788 independent-study chain rows,
22,898 model-contact chain rows, 11,692 model-quality chain rows and 758 native
chain rows (metrics restricted to the selected representative).

The live NLVPMVATV structure page loaded 3,355 cards, five selected structures,
contact counts, interface ipTM and model confidence, motif buttons, downloads,
and the beginning of the complementarity overlay. A full DOM snapshot timed out
on this large page; screenshot inspection succeeded. No timing equivalence claim.

The metadata-preparation tool incorrectly treated observed zero ipTM/confidence
as missing for percentile calculations (`value or NaN`); fixed to preserve zero.
The focused Python adapter checks pass. TCRen SVG rendering was smoke-checked
with synthetic geometry; real structure annotation/export remains to be checked
on Stratos. The Stratos sync monitor reported roughly two hours remaining at
20:15 Moscow time; its owning chat manages that transfer.


## Delivery checkpoint

The opt-in ZIP was built and inspected: 49,572,834 bytes, eight members; large
TSV/TXT members retain gzip compression. Final focused checks: 69 database tests
passed / 3 skipped; 108 affected Scala tests passed; two standalone Python adapter
checks passed. Rich/ordinary legacy release roles remain unchanged. Vectorized
Polars joins/pivots suffice for export; existing batched native annotation APIs
own TCRen map preparation. No new parallel search pool was introduced.

Stratos toolchain image `vdjdb-legacy-interop:local` built successfully from pinned
Java 8 and Node 11 images using `vdjdb-web/tools/Dockerfile.interop`. The source
sync is still active; runtime hosting and measured nonregression are pending.
Keep the goal active until the remaining acceptance checks above are resolved.
