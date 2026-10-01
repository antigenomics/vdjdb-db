# The confidence score

The final stage of database processing assigns a confidence score to each TCR:peptide:MHC complex, computed from the reported ``method`` entries.

The score evaluates TCR sequence confidence, identification confidence and verification confidence, by the following criteria:

1. Ensuring TCR sequence is correctly identified according to ``method.sequencing`` and ``method.singlecell`` (1-3 points)
    * sanger - several cells sequenced (2+ cells sequenced according to ``method.frequency``) - 2 points, otherwise 1
    * amplicon-seq - frequency is at least ``0.01`` **and** at least 2 reads (``method.frequency.count``) - 3 points, otherwise 1
    * single-cell - 3 points if performed
2. Initial identification of TCR:pMHC is correct according to ``method.identification`` (0-1 point)
    * sort-based - frequency is higher than ``0.1`` according to ``method.frequency``)
    * culture-based - frequency is higher than ``0.5``
    * limiting dilution/culture prior to sequencing - ``method.frequency`` is ambiguous here, check whether it is higher than ``0.5``
3. Verification T-cell specificity (0-3 points)
    * direct method - 3 points, e.g. has PDB id (``meta.structure.id`` is not empty) or some other method that directly evaluates TCR:pMHC binding
    * target stimulation-based - 2 points
    * staining-based - 1 points
    * If verification is performed, the TCR sequence is assumed to be known, so the score from part ``1.`` is set to 3

The final score is the minimum of the part ``1.`` score and the sum of the part ``2.`` and part ``3.`` scores.

Among records that point to the same unique complex entry, that is the same set of unique **complex** fields - independent submissions, replicas, and so on - the maximal score is taken.

score | description
------|----------------------
0     | Low confidence/no information - a critical aspect of sequencing/specificity validation is missing
1     | Moderate confidence - no verification / poor TCR sequence confidence
2     | High confidence - has some specificity verification, good TCR sequence confidence
3     | Very high confidence - has extensive verification or structural data



## Score rules

The table below is rendered from `vdjdb.score.confidence` while this page builds.

```{vdjdb-score-rules}
```
