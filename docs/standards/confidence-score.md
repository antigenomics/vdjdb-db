# The confidence score

The final stage of database processing assigns a confidence score to each TCR:peptide:MHC complex, computed from the reported ``method`` entries.

The score combines sequencing confidence, identification confidence and verification
confidence. The rule summaries below are generated from the scoring functions.
Sequencing confidence uses the sequencing method, single-cell flag,
frequency and available support counts. Identification confidence uses the assignment
method and frequency; verification confidence uses the reported verification method.

The following thresholds describe `vdjdb.score.confidence`. Frequency is a fraction
from 0 to 1. Within each component, the first matching row takes precedence;
method matches use substrings of the lowercased field.

| Component | Criterion | Points |
|---|---|---|
| Sequencing | Nonempty single-cell flag other than `no` | 3 |
| Sequencing | `sanger` | 3 with support count at least 2; otherwise 2 |
| Sequencing | `amplicon-seq` | 3 with frequency at least 0.01 and support count at least 2 when known; otherwise 1 |
| Sequencing | Other or missing method | 1 |
| Identification | `culture` or `cloning` | 1 with frequency at least 0.5; otherwise 0 |
| Identification | `sort`, `beads`, `separation` or `stain` | 1 with frequency at least 0.05; otherwise 0 |
| Identification | `targets` | 1 with frequency at least 0.25; otherwise 0 |
| Verification | `direct` | 3 |
| Verification | `targets` | 2 |
| Verification | `sort`, `beads`, `separation` or `stain` | 1 |

Unmatched identification or verification methods contribute 0 points. The Sanger
support count uses `method.frequency.count` where available, otherwise the numerator
of a ratio frequency. Unknown counts do not disqualify the amplicon frequency rule.

For an observation without a structure, the score is the smaller of sequencing
confidence and the sum of identification and verification confidence. A verification
method recognised by the scorer sets sequencing confidence to 3. A nonempty
`meta.structure.id` sets the score to 3; QC requires a PDB entry identifier in that field.

Among records that point to the same unique complex entry, that is the same set of unique **complex** fields - independent submissions, replicas, and so on - the maximal score is taken.

This maximum groups the alpha and beta junctions and V/J calls, receptor species,
both MHC chains, MHC class and peptide. It does not group by donor or publication.
A shared score therefore does not mean every observation received the strongest assay.

score | description
------|----------------------
0     | Low confidence/no information - a critical aspect of sequencing/specificity validation is missing
1     | Moderate confidence - no verification / poor TCR sequence confidence
2     | High confidence - has some specificity verification, good TCR sequence confidence
3     | Very high confidence - has extensive verification or structural data



## Interpreting the score

The score ranks support for the sequence-to-pMHC annotation. It is not a probability,
an affinity estimate, a measure of killing or cytokine production, or a clinical-response
prediction. Multimer binding and functional responses establish different properties;
retain the methods alongside the score. Single-cell sequencing improves sequence
confidence but still requires evidence linking that sequence to the reported pMHC.
See {doc}`chunk-format` for assay attribution and unresolved peptide pools.

Method vocabulary and score matching are distinct. A source-reported method can pass
vocabulary review without matching a scoring rule. Do not rename a readout to obtain
more points or infer an unreported validation experiment from a high shared score.

## Score rules

The table below is rendered from `vdjdb.score.confidence` while this page builds.

```{vdjdb-score-rules}
```
