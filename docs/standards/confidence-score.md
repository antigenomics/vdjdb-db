# The confidence score


At the final stage of database processing, TCR:peptide:MHC complexes are assigned with confidence scores. Scores are computed according to reported **method** entries.

VDJdb scoring is performed by evaluating TCR sequence, identification and verification confidence based on the following criteria:

1. Ensuring TCR sequence is correctly identified according to ``method.sequencing`` and ``method.singlecell`` (1-3 points)
    * sanger - several cells sequenced (2+ cells sequenced according to ``method.frequency``) - 2 points, otherwise 1
    * amplicon-seq - frequency is higher than ``0.01`` - 2 points, otherwise 0
    * single-cell - 3 points if performed
2. Initial identification of TCR:pMHC is correct according to ``method.identification`` (0-1 point)
    * sort-based - frequency is higher than ``0.1`` according to ``method.frequency``)
    * culture-based - frequency is higher than ``0.5``
    * limiting dilution/culture prior to sequencing - the ``method.frequency`` becomes somewhat ambigous, check if it is higher than ``0.5``
3. Verification T-cell specificity (0-3 points)
    * direct method - 3 points, e.g. has PDB id (``meta.structure.id`` is not empty) or some other method that directly evaluates TCR:pMHC binding
    * target stimulation-based - 2 points
    * staining-based - 1 points
    * If verification is performed, then the TCR sequence is assumed to be known, so score from ``1.`` is set to 3

The final score is then calculated as minimal between score from part ``1.`` and sum of scores from part ``2.`` and part ``3.``.

Maximal score is then selected among different records (independent submissions, replicas, etc) pointing to the same unique complex entry (i.e. set of unique **complex** fields).

score | description
------|----------------------
0     | Low confidence/no information - a critical aspect of sequencing/specificity validation is missing
1     | Moderate confidence - no verification / poor TCR sequence confidence
2     | High confidence - has some specificity verification, good TCR sequence confidence
3     | Very high confidence - has extensive verification or structural data



## The rules, from the source

Rendered from `vdjdb.score.confidence` while this page builds, so a change to the score
cannot leave the prose behind.

```{vdjdb-score-rules}
```
