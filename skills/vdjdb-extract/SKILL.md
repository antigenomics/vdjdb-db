---
name: vdjdb-extract
description: Extract TCR:pMHC specificity records from raw submission sources - supplementary XLS/CSV tables, PDF manuscripts, 10x Genomics contig and clonotype files, AIRR Rearrangement TSVs, Adaptive ImmunoSEQ exports - into a VDJdb chunk TSV, verifying every extracted sequence, gene, allele and reference id back against the source. Use when a paper, dataset or submitter folder has to become a chunk file, before vdjdb-format and vdjdb-proofread.
---

# Extract publication observations

Read [AUTHORITIES.md](../AUTHORITIES.md) and the
[chunk specification](../../docs/standards/chunk-format.md). Respect the requested source and
technology scope. A reference already present in VDJdb does not establish complete coverage.

For deposited complexes, use [structures](../vdjdb-structures/SKILL.md) and the local tcren
pipeline. New observations go to the PDB aggregate with their publication reference; inventory
reused templates separately. Inspect complete polymer sequences as well as resolved coordinates.

## Reconcile the sources

Inventory the supplied article, supplements and processed assignments. Record which tables
contain sequences, pairing, peptides, restrictions and assays. Join on reported clone, donor,
well or barcode identifiers and record join cardinality. Do not infer alpha/beta pairing from
row order or oligoclonal coexpression. Preserve ambiguous extra chains for review.

One row reports one observation, with both chains of a reported pair on that row. Verify literal
sequences, calls and identifiers against the source. A positive specificity needs at least one
reported junction, a literal peptide and supported restriction. A gene or protein label alone
is not a peptide assignment. Never derive the tested peptide from a canonical protein sequence.

Separate positive observations, explicit assay negatives, unassigned sequences, reused literature
comparisons and unavailable assignments. Give every supplied observation a disposition. Do not
turn an assay failure into a negative if a separate binding identification remains supported.

## Preserve the submitted biology

The canonical junction starts with `C` and ends with `F` or `W`.
Use AIRR `junction_aa`, which includes both anchors; `cdr3_aa` excludes them. If only an IMGT core
is printed, retain that literal core and inspect the build repair. Preserve terminal flanks and
reported V/J values when assembly resolves them. Internal cysteines and noncanonical anchors
are review findings, not permission to discard an observation.

Invalid amino-acid symbols, incomplete sequences and chemical modifications require a recorded
source disposition outside the positive build. Keep the original source; do not silently delete
it or substitute standard residues. Ask only for a curation decision that is still unresolved.

For spreadsheets, inspect repeated headers, formula cells, functionality suffixes and swapped
J/D fields. Do not strip allele suffixes merely because a cell also contains a functionality code.
For 10x, barcode pairing and clonotype aggregation have different meanings; respect the published
assignment and the user's choice of vendor versus reanalysis data.

## Record methods and provenance

Use source-supported vocabulary for identification, sequencing, pairing and verification.
`cultured-T-cells` requires culture before receptor sequencing. Later culture of a reconstructed
receptor does not establish it. Do not infer tetramers from an unspecified multimer or sorting
step. Verification describes the retested receptor, not the initial identification.

Keep clonotype frequency separate from subset frequency. Reported counts, totals, fractions and
percentages retain their meanings. Optional donor and method detail stays blank when unknown.
Identifiers in metadata must come from the publication, not generated export row labels.

Use a verified PMID, DOI, patent or PDB reference. A submission without a publication uses its
VDJdb issue URL. Never guess a PMID. Keep extraction evidence in the local execution record,
with source table/row, literal values, joins, exclusions and unresolved questions.

Write the current schema header with UTF-8, tabs, LF endings and empty cells for missing values.
Continue with [format](../vdjdb-format/SKILL.md) and [proofread](../vdjdb-proofread/SKILL.md).
