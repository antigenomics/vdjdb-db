---
name: vdjdb-structures
description: Extract and proofread deposited TCR-pMHC observations with the local tcren pipeline, then update the PDB aggregate chunk. Use for a structural publication or missing PDB entries; exclude predicted docking models and complexes without supported receptor-pMHC assignments.
---

# Collect structural observations

Read [AUTHORITIES.md](../AUTHORITIES.md), [extract](../vdjdb-extract/SKILL.md) and
[proofread](../vdjdb-proofread/SKILL.md). The structural aggregate is
`chunks/PDB_Database.tsv`; retain the publication reference on each row and the accession in
`meta.structure.id`. A multi-PDB paper has one intake issue, not a chunk for each structure.

## Establish the source set

Read the paper, supplements and PDB primary citations. Separate newly deposited structures,
reused templates, new assays and predicted models. Inventory every accession with its reported
receptor, peptide, restriction and assay outcome. A cited template is not a new experimental
observation. Distinct deposited structures can remain distinct observations under one publication.

Fetch mmCIF coordinates and entry/polymer metadata. Check each response independently; a failed
first download must not hide the other accession results. Unavailable coordinates are a source
blocker: keep the issue open with the accessions and admission condition. Do not make empty or
sequence-free chunks. Preserve downloadable source attachments when shortening issue text.

## Use tcren locally

Use `~/vcs/code/tcren` and read its `AGENTS.md` and
`~/vcs/code/tcren/skills/tcren/SKILL.md`. Install it into an
isolated extraction environment when needed. This is a curation tool, not a source-checkout
replacement for this repository's published build dependencies.

Parse the complete set with `tcren.structure.parse_structure`, annotate with
`tcren.annotation.batch.annotate_batch`, classify with the resulting precomputed records, and
run `tcren.mhc.annotate_mhc_batch` once over the set. Annotation is batched; never put mmseqs
calls in a per-structure loop or a process pool. Resolve reference roots explicitly and use a
separate output/data root for the extraction.

Keep chain identifiers, literal polymer sequences, coordinate sequences, annotations and allele
candidates in the extraction proof. Inspect chain roles and the peptide in the groove; short
protein fragments are not automatically epitopes. Single-chain pMHC constructs need their reported
peptide separated from linkers, beta-2-microglobulin and the MHC chain.

## Check what the coordinates omit

Compare coordinate residues with the deposited complete polymer sequence. Unresolved residues
can shorten a junction without changing the submitted construct. Recover missing residues only
from the complete deposited sequence or an explicit source sequence; do not copy a template's
CDR3 into an engineered receptor.

VDJdb uses `junction_aa`, including Cys104 and Phe/Trp118. tcren CDR3 region markup excludes the
anchors. Check the FR3/FR4 sequence and the J motif; adding `C` and `F` by concatenation is not a
valid conversion. Keep inferred gene/allele calls distinguishable from published assignments.

An MHC groove sequence can match several alleles equally. A top hit, including a null-expression
allele, is not proof that the experiment used it. Use the reported construct and primary PDB
annotation to resolve the allele; retain reported resolution or flag a material contradiction.
Inspect engineered MHC substitutions and nonstandard peptide residues explicitly.

## Check specificity and publish

A covalent clamp, forced docking model or common orientation does not establish specificity.
Review uncoupled binding and functional assays. Failure to activate does not erase separately
measured binding; preserve the distinct assay outcomes. Engineered receptors require their own
literal sequences and evidence, not the parental receptor's assignments. Culture of a reconstructed
receptor does not establish culture before sequencing.

Format supported observations in the current chunk schema. Keep unavailable, unassigned,
explicitly negative and chemically modified observations in their appropriate source dispositions.
Use [duplicates](../vdjdb-duplicates/SKILL.md) to compare complete observations and
[publish](../vdjdb-publish/SKILL.md) for a PDB/data-issue branch, registry reconciliation,
separate validation and release-aware checks. Close the intake only when its material scope is
accounted for; a PDB citation match alone is not complete curation.

## Reported peptides and predicted cores

Preserve the complete assay peptide and source-supported MHC/provenance. Use standalone
`vdjdb assess-epitopes` during review when predicted presenters, binding registers or TCR-facing
representations help resolve a comparison. Follow the
[assessment guidance](../vdjdb-build-integrity/SKILL.md#choose-assessment-scope). Predictions
do not supply a missing experimental epitope or establish equivalent recognition. A predicted
core can omit required flanks or insertions and does not replace a structure's full ligand.

## Audit hosted coordinate files against current records

Use the independent audit after current tables have been assembled:

```bash
uv run python processing/proofread_structure_files.py --structures out/inputs/structures \
  --tables out/tables --out out/reports/structures --fetch
```

`proofreading/structure_inputs.tsv` pins the hosted raw files by revision, size and SHA256.
Check its current download sizes and available disk first.
`--fetch` stages those inputs separately; omit it when checking an existing local mirror.
The audit verifies every digest and streams all declared archives without extraction. It checks
current links against reviewed coordinate-sequence metadata and verifies coverage for every native
accession in current records. CI does not independently re-extract chain sequences.
The current scope is 367 native accessions across all chunks; the older 360-accession count
covers only the PDB aggregate chunk. No chunk or coordinate is edited.

Inspect `structure_files.parquet`, `native_structure_files.parquet`,
`supplemental_structure_files.parquet` and `structure_files.json`
in the report directory. Keep original/repaired sequence matches, exact receptor:pMHC matches
and historical lineage links distinct. A historical association is not sequence proof.
The [supplemental catalogue](https://huggingface.co/datasets/isalgo/vdjdb_structure_models/tree/main/metadata/supplemental)
contains 11,911 distinct coordinate contents, with 11,388 new and 523 existing contents. Its
16,118 sequence-state links connect 4,777 structures to 7,292 current records; they do not count
independent experiments. Use its README and verification report for dispositions and input identity.
Audit supplemental raw assets and their current links separately. Ignore input positive/negative,
pairing and validation flags for linkage. These inputs identify receptor/pMHC sequences and files,
not binding outcomes. Outcome labels come only from curated negative chunks or reviewed 10x
removals where applicable. Require literal coordinate/sequence/context agreement regardless of flags;
keep native/predicted coordinate provenance separate from experimental validation.
Optional `vdjdb build --structure-supplemental <metadata.tsv.gz>` merges their raw metadata within
one build. Leave it off in default hosted assembly/scoring until matched-cohort extraction is validated.
CI uploads these reports separately from the release bundles. See
[structure cross-links](../../docs/outputs.md#structure-files-and-current-record-cross-links).

For record confidence v3, complete motifs and run `vdjdb score-records --tables out/tables`
with the `confidence` extra. The frozen model uses accepted predicted-structure metrics;
native coordinates remain separate evidence. Parameters stay frozen as available evidence changes;
held-out assay metrics use the complete original structural inputs. Default membership is hard;
exploratory `--soft-membership` and its batched search/count fields are documented in the
specification. See
[record confidence v3](../../docs/standards/record-scoring-v3.md) for its precise inputs and outputs.
