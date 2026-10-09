---
name: vdjdb-harmonize
description: Canonicalise antigen.gene and antigen.species in a VDJdb chunk against the epitope dictionary and the gene and species alias tables, detect spurious values (UniProt descriptions, "[species]" annotations, "Probable"/"Chain A," prefixes, multi-word organism names), resolve blanks from the publication, flag one epitope that is a substring of another, and extend the alias tables where they have no entry. Use when antigen.gene or antigen.species carry free text rather than symbols, or when vdjdb-proofread reports antigen fields that disagree across a chunk.
---

# Review peptide provenance and aliases

Read [AUTHORITIES.md](../AUTHORITIES.md) and
[terminology](../../docs/standards/terminology.md). `antigen.gene` and `antigen.species` describe
peptide provenance. Recognition is assigned to the peptide-MHC complex.

## Inspect existing harmonisation

Use a fresh `vdjdb build --out out/`. Read `harmonisation.tsv`, `lookalikes.tsv` and
`epitope-sources.tsv` in `out/reports/`. The build applies declared epitope corrections and the
reviewed gene/species alias tables; preserve source cells already resolved by those rules.

The authorities are `patches/antigen_epitope_species_gene.dict`,
`proofreading/gene_aliases.tsv` and `proofreading/species_aliases.tsv`. Species substring rules
are ordered: more specific matches precede broader ones. Review a new alias in its governing
organism and retain evidence for the addition.

## Resolve findings from the publication

Inspect free-text protein descriptions, bracketed organism names, prefixes, case/separator
collisions and blanks. Human and mouse gene symbols can differ correctly. A same-species spelling
collision is a review candidate, not proof of a typo.

One peptide can occur in several organisms or proteins. Retain source-supported provenance;
do not force a single global assignment from sequence equality. Source variation, synthetic
peptides and unknown parent genes also require different readings. Unknown optional provenance
stays blank and is advisory. Do not write an `Unknown` placeholder or infer the organism from
another paper reporting the same peptide.

Report peptide substrings as candidates for source review. A minimal epitope and a longer tested
peptide can both be valid. Chemical modifications cannot be represented by silently substituting
standard amino acids; preserve the source and link the representation question to its issue.

## Record the result

For each change, record file/row, old/new value and source or authority. Add reviewed alias/rule
entries when needed, then inspect the resulting build. Chunk edits require a chunk/data-issue
branch, identity update and separate validation. Leave optional provenance follow-ups open without
blocking otherwise supported observations. Re-run `vdjdb qc` and the submission report.

## Choose assessment scope

Start triage with `qc` and corpus-relative `submission`. For peptide/MHC review, use
`vdjdb assess-epitopes <chunk>` with the pinned reference, or opt-in `assessment.yml` CI when
authorized. Selected support counts describe only the chosen input. Follow the
[assessment guidance](../vdjdb-build-integrity/SKILL.md#choose-assessment-scope) for coverage,
core/TCR-facing fields and partial versus full validation. Preserve reported assay peptides,
restrictions and parent provenance; shared predicted cores do not justify trimming or merging.
Full combined-corpus integration checks remain required before merging.
