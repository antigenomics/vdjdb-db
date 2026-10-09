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
