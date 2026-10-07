---
name: vdjdb-duplicates
description: Review repeated paired receptor-pMHC observations and metadata synonyms within and across VDJdb publications, preserving distinct experiments and reports. Use for a corpus audit or a submission that may overlap existing records.
---

# Review repeated observations

Read [AUTHORITIES.md](../AUTHORITIES.md) and
[observation identity](../../docs/standards/chunk-format.md#completeness-and-observation-identity).
Use a fresh build and read raw chunks with `read_chunks(deduplicate=False)`.

## Compare receptors and experiments

Before importing a collection, measure epitope plus CDR3 overlap against both the last release
and current chunks, within species and for each chain. Compare submitted and repaired junctions
separately. Report matched and unmatched rows and distinct pairs per original publication, with
the source files that matched. This broad screen detects earlier imports; it does not prove
that paired receptors, restrictions or experiments are identical.

For substantial overlap, reconcile the matched observations with their primary papers before
creating a chunk. Verify paired chains, V/J calls, restriction and experiment metadata using
the comparison below. Keep independently reported experiments; do not import the same source
observation again under a collection's reference. Review unmatched measurements against the
original assay before classifying them as positive or negative.

Within species, group both chains together using `cdr3.alpha`, `v.alpha`, `j.alpha`,
`cdr3.beta`, `v.beta`, `j.beta`, `antigen.epitope`, `mhc.a`, `mhc.b` and `mhc.class`.
Compare submitted values first, then harmonised values. Sharing a beta chain does not establish
that paired receptors match.

For each repeated group, compare the reference and every `method.*` and `meta.*` field,
including donor, clone, replica, tissue, subset, pairing and frequency fields. Distinct reported
experiments remain separate. Submitter, comment and curation serial alone do not distinguish
experiments. A blank versus a populated value requires review before merging complementary data.

Different publications remain separate reports. Author overlap can suggest a follow-up or reused
cohort, but an author-count threshold cannot prove independent experiments or justify deletion.
Reconcile a preprint and publication only with a verified publication link and matching source
observations. Preserve separately reported PDB structure observations.

Report exact within-file duplicates, same-publication cross-file groups, distinct experiments
within a publication, and cross-publication recurrence separately. `vdjdb submission <chunk>`
provides a recurrence summary; it does not replace the paired comparison.

## Inspect terms and missing values

Count every field's distinct values and blanks. Inspect the most frequent junctions and singleton
metadata values with their file/row locations. A singleton clone identifier or a common public
junction is expected and is not an error by itself.

Read `out/reports/lookalikes.tsv` and `epitope-sources.tsv`. Case or separator differences are
review candidates. Human and mouse gene spelling can differ; antigen aliases need the governing
organism. Compare method token sets rather than treating token order as a different experiment.

Several restrictions for one peptide can be supported observations. A broad allele call and a
more precise call are different reported resolutions; do not upgrade the broad call without
source evidence. A predicted presentation is not evidence that the experiment used that allele.

## Record decisions

Keep per-group evidence and counts in the local execution record. Report denominators as
observations, clonotypes or publications. Link unresolved source questions to the publication issue;
create a data issue for a mechanical corpus correction. Apply chunk changes only on the branch
for that chunk or data issue, with identity reconciliation and release-aware validation.
