---
name: vdjdb-use
description: Query VDJdb with vdjmatch, inspect specificity-record evidence and confidence, retrieve released TCRNET or TCREmp motifs, and download linked native or predicted structures. Use for consuming the database, not curating submissions or rebuilding it.
---

# Use VDJdb observations and evidence

Read the relevant section of [the usage guide](../../../docs/using-vdjdb.md).
Use existing released data or the user's completed build. Do not launch a build,
motif recomputation, model fitting or new scientific evaluation to answer a lookup.

## Select inputs and tools

- Record the database snapshot, checksum and tool versions. Primary Parquet tables
  preserve record IDs and linked evidence; the legacy projection serves older readers.
  Check actual file availability in the selected snapshot before promising a score or motif.
- Use the installed `vdjmatch annotate --help` and its published workflow for receptor
  annotation. Pass `--vdjdb` for a pinned ZIP or a user-specified build's `tables/`
  directory. Do not search a parent folder containing several builds or choose one by mtime.
- Use `curated-benchmark` for the recommended annotation population; use explicit
  `all-valid` when the user requests the full eligible reference or historical reproduction.
  State which population was searched. Batch queries once with an explicit CPU budget.
- Require anchor-inclusive junctions. Preserve chain, species, available V/J calls,
  peptide and MHC context, and reported pairing. Do not fabricate anchors or pairs.

## Follow the evidence

Join `records`, `chains`, `evidence` and `record_scores` by `record_id`. A target-level
vdjmatch result is not a record-ID mapping. Recover exact observations using chain and
pMHC keys, retaining multiple publications; distinguish exact matches from neighbours.
Keep report, chain, clonotype, publication and structure counts separate.

Read experimental identification and verification methods alongside the legacy 0-3
score. For v3.2, report `validation_score` and `independent_validation_score` separately
with their versions and model digests. The latter predicts another publication's report
at the same pMHC, not an independent assay outcome. Missing decisions are not negatives;
`paired_lower` and `paired_upper` are sensitivity bounds, not confidence intervals.

Retrieve motif membership and PWM tables from the same release. Keep TCRNET and TCREmp
separate. Filter species, chain and peptide before selecting `cid`, and select one
junction length before making a PWM or logo. Motif groups pool MHC restrictions; do
not describe an exported MHC annotation as evidence of an allele-specific motif test.

## Retrieve coordinates

Use the Hugging Face dataset's current metadata and archive inventory at one pinned
revision. Join `record_links.source_row` through `asset_sources` to `asset_inventory`;
select the intended original/repaired sequence state and `exact_tcr_pmhc_match`.
Compare the metadata's database revision and sequences with the user's snapshot.
Download only required archives and extract only selected members, preserving existing
files. The supplemental catalogue has a separate `structure_id` join and sequence-state
field; consult its README instead of applying the historical schema to it.

Coordinate availability is not experimental validation. Preserve native/predicted
provenance and distinguish generator confidence from assay evidence. Source positive,
negative or pairing flags do not establish outcomes. Never admit a withdrawn assignment
merely because its coordinates remain in an archive. An empty lookup is limited to the
chosen catalogue and snapshot.

Use [the structure-curation skill](../../curator/vdjdb-structures/SKILL.md) only when the task is to
extract, proofread or submit a structural observation. Return the requested records,
motifs or files with their identifiers and provenance; keep interpretation limits together.
