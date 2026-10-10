---
name: vdjdb-proofread
description: Validate a VDJdb chunk with vdjdb qc and vdjdb submission, explain every finding with a specific suggested fix and the authority behind it, resolve method-field and MHC problems the rules cannot decide, read the junction-versus-germline and score reports, and decide whether the chunk lands in chunks/, pending/ or withheld/. Use as the last gate before a chunk is committed, or on any chunk suspected of having a column shift, bad gene or allele names, blank MHC fields or missing method information.
---

# Proofread a submission

Read [AUTHORITIES.md](../AUTHORITIES.md), the
[chunk specification](../../docs/standards/chunk-format.md), and
[submission guide](../../docs/submission.md). Review the extraction evidence if present.

For structural sources, follow [structures](../vdjdb-structures/SKILL.md): check complete polymer
sequences, unresolved coordinates, allele ties and evidence from forced or clamped complexes.

## Run the checks

```bash
uv run vdjdb qc <file> --strict --report out/reports/qc.tsv
uv run vdjdb submission <file>
```

Read per-row findings and `qc-summary.tsv`; their advisory flag determines whether a finding is
fatal. Do not maintain a second list of rules or copy historical finding counts into this skill.
Inspect header shape and content together to catch a field shift that preserves row width.

Every positive observation needs a reported junction and literal epitope. Class I requires its
valid first chain and `B2M`. A partially reported class-II restriction may retain an unknown
partner blank. Both MHC chains absent, contradictory class, invalid names and unsupported
positive assignments block the affected records. Missing V/J or optional peptide provenance does
not justify discarding a sequence-bearing observation.

## Inspect sequence repair

Build from the supplied source values. Compare `cdr3.original`, `v.segm.submitted` and
`j.segm.submitted` with the repaired sequence, shipped calls and fix types. Inspect at least one
sequence-bearing observation per chunk and include repair/call-change cases when present. Count
unresolved cases separately from successful repairs.

The canonical junction starts with `C` (Cys104) and ends with `F` or `W` (Phe/Trp118).
Missing anchors, internal cysteines, partial germlines
and conflicting calls require interpretation with the source and the germline reports. They are
not all transcription errors. Preserve literal cores/flanks when the build resolves them. Do not
edit a source sequence merely to make a canonical flag true. A proposed V is inferred, not reported.

## Review metadata and repeats

Use the score histogram and assay description to check identification versus verification,
sequencing, pairing and culture before sequencing. A low score alone does not prove a method-field
error. Optional detail stays blank when the paper does not provide it.

For multiple-chain clonotypes, verify the
[candidate-pair convention](../../docs/standards/chunk-format.md#multiple-chains-in-one-clonotype):
all chains come from the same source clonotype, every candidate retains the original
`meta.clone.id` and experiment metadata, and `method.pairing=ambiguous` plus the comment/issue
describe unresolved pairing. Check that expansion preserves source counts without treating the
candidate rows as independent clonotypes. Do not quarantine solely because more than one chain
was reported.

Count distinct values and blanks in every field. Inspect frequent junctions and singleton metadata
with source locations. Use [harmonize](../vdjdb-harmonize/SKILL.md) for provenance aliases and
[duplicates](../vdjdb-duplicates/SKILL.md) for paired receptor-pMHC recurrence. Compare complete
experiment metadata before merging; preserve different experiments and publication reports.

Separate closed vocabularies from descriptive metadata. Cohort, subset, tissue and source identifiers
allow free text. Check aliases in the assembled output before editing source spellings. Prefer a
verified PMID, but retain vendor submission and unpublished structure references when no publication
link is established. A reanalysis does not supply independent measurements of the original dataset.

## Resolve and report

Fix supported mechanical defects under the existing authorization. Request a curator decision
only for an unresolved material contradiction, with the evidence and affected rows. Avoid repeated
retrieval attempts for inaccessible papers; use supplied definitive source tables within scope.

Report raw/built counts, score histogram, fatal/advisory findings, successful/unresolved repairs,
repeated-observation decisions and remaining questions. Keep source proofs in the local execution
record, not in chunk comments.

A current-format submission blocked by a missing build reference goes unchanged to `pending/`.
An unreadable older-format submission needing re-export goes to `withheld/`. Preserve explicit
negatives in `chunks_negative/`; they do not enter the positive build. Name the path, actionable
blocker and condition for admission on the issue, and leave it open. Follow
[publish](../vdjdb-publish/SKILL.md) for separate data/validation commits and release-aware gates.

## Choose assessment scope

Start triage with `qc` and corpus-relative `submission`. For peptide/MHC review, use
`vdjdb assess-epitopes <chunk>` with the pinned reference, or opt-in `assessment.yml` CI when
authorized. Selected support counts describe only the chosen input. Follow the
[assessment guidance](../vdjdb-build-integrity/SKILL.md#choose-assessment-scope) for coverage,
core/TCR-facing fields and partial versus full validation. Preserve reported assay peptides,
restrictions and parent provenance; shared predicted cores do not justify trimming or merging.
Full combined-corpus integration checks remain required before merging.
