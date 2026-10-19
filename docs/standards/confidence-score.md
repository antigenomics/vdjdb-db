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

## Optional evidence classifier

The `confidence` extra provides an experimental classifier separate from the build's
existing score. It uses nine features from completed build tables: the existing legacy
score, TRA/TRB publication totals, the integrated donor total, four TRA/TRB motif flags
(TCRNET and TCREmp), and mean predicted ipTM. It does not reconstruct assays from chunks
or model allele-specific effects. Fitted models and outcomes are caller-supplied inputs.

```python
from vdjdb.score.evidence import FEATURES, categorical, integrate
from vdjdb.score.bayes import fit, predict_marginals

features = integrate(records, chains, evidence, queries)
x = categorical(features.select(FEATURES))
model = fit(x.to_dicts(), outcomes)  # independently supplied binary outcomes
probabilities = predict_marginals(model, x.to_dicts())
```

Frames use the internal names returned by `vdjdb.emit.vdjdb3.read_table`.
Queries supply `species`, `antigen.epitope`, `mhc.a`, `mhc.b`, `mhc.class`,
`a_junction`, `a_v`, `a_j`, `b_junction`, `b_v`, `b_j`. Junctions include both
conserved anchors; an absent chain is empty. V/J calls are grouped at gene resolution.

The legacy score is used as built, taking the maximum across the integrated set.
Publication totals count distinct references per chain/peptide/V-gene key. Donors
count distinct nonempty publication/donor contexts for the integrated receptor:pMHC.
Motif membership means a matching built entry belongs to a motif; explicit nonmembers
are zero, while missing annotations remain unknown. The existing structure links supply
ipTM, averaged over distinct linked models.

Counts use bins `0`, `1`, `2`, `3+`; ipTM uses below 0.5, 0.5 to below 0.8, and
at least 0.8. The legacy score keeps 0-3. The Bayesian model learns a tree over the
categorical evidence conditional on the outcome, with Dirichlet smoothing. Missing
fields and states absent from training are marginalized; all missing evidence returns
the fitted prior. The returned model serializes to JSON. The caller owns grouped
validation and the interpretation of the outcome; assay validation is not an affinity.

`predict_marginals` returns `paired`, `TRA` and `TRB` probabilities from the same fitted
network. The chain marginals condition only on that chain's publication and motif evidence;
the other chain and paired legacy/donor/structure evidence are integrated out. Each predicts
validation of the full interaction. It does not reinterpret a paired assay as a single-chain
binding experiment. No extra model or imputed partner is needed.

For a simpler comparator, `fit(..., structure="naive")` retains the class-to-feature
edges and omits the feature-dependence tree. `categorical(..., iptm_bounds=(low, high))`
accepts two ordered cutoffs for the three ipTM grades; the defaults remain `(0.5, 0.8)`.
Store any selected cutoffs with the fitted model and use the same encoding at inference.
Choose smoothing, structure and cutoffs inside training groups; assess the selected
procedure on held-out groups, including inputs with missing structural or motif evidence.

`equivalent_sample_size` also accepts a map with every feature name and `validated` as
keys. This permits stronger smoothing of one sparse conditional table while retaining
the other priors. If ipTM is a leaf, changing only its conditional table leaves every
prediction without observed ipTM unchanged: marginalizing that table sums it to one.
The supplied model must retain its grade boundaries and all fitted tables together;
paired and chain-only predictions use that same network.

`structure="tan_leaf"` learns the nonstructural tree first and attaches ipTM to the
feature with greatest conditional mutual information given validation. Lexical ties
are deterministic. This uses at most seventeen edges while preventing structural grading
from changing the remaining graph. ipTM has no children, so missing structural evidence
is marginalized without altering the nonstructural network.

`structure="tan_independent_iptm"` removes the additional feature parent of ipTM.
It retains the nonstructural dependence tree, with at most sixteen edges, and models
structural evidence as conditionally independent of those features given validation.

For deployment, `predict_evidence(model, numeric_frame)` uses the saved `iptm_bounds`
and returns paired, TRA and TRB probabilities from that one fitted network. Omitted
feature columns are unobserved; extra identity columns are not predictors. The model
does not take an epitope name as an input. Keep the input evidence definitions identical
to training; application to an unmeasured epitope does not establish its validation accuracy.

`score_records(model, numeric_frame)` retains all three forecasts and adds
`record_score` (equal to `paired`), `paired_lower`, `paired_upper`, `pair_supported`,
`direct_pair_evidence` and nullable `positive`. Missing pair evidence never forces the score to zero.
The existing BN already combines the chain features jointly; TRA and TRB are
correlated views of the same assay endpoint, so do not multiply their probabilities.

The range varies missing chain features and missing ipTM over their observed training
states, including every joint completion. Unreported metadata stays marginalized and
the nominal forecast is included. Without direct pair evidence, the range also includes
the TRA-only and TRB-only forecasts, requiring their agreement for a confident call.
States absent from training supply no fitted evidence
and are treated as missing, as in ordinary prediction. This is a sensitivity envelope,
not a posterior credible interval or calibrated coverage guarantee. It identifies
forecasts that depend on unavailable chain or structure evidence. It can be broad
without an ipTM or a supported partner, while keeping the paired point forecast usable.

One threshold, 0.5, applies throughout. `positive` is true if the lower limit is at least
0.5, false if the upper limit is below 0.5, and unknown otherwise. Report unknown calls
and their count; do not turn them into negatives or drop them from evaluation.
`pair_supported` describes positive support for both chains or direct pair evidence
(an exact-reference annotation, donor support or matching predicted structure). It no
longer gates the forecast or forces the decision. The fitted network and TRA/TRB probabilities
are unchanged. Direct evidence for a single-chain record cannot be reused as exact-pair
evidence for that chain with a new partner; an unrelated model cannot supply ipTM.

Motif inputs to `predict`, `predict_marginals` and numeric `predict_evidence` may also be
probabilities on [0,1]. A supplied membership probability `q` mixes the predictions
conditional on nonmembership and membership, with weights `1-q` and `q`. Values 0 and 1
recover the hard predictions; missing remains marginalized. With several uncertain motif
inputs, the mixing distribution assumes their uncertainty is independent. This assumption
concerns the supplied memberships; the fitted BN still models its learned dependencies.
A supplied soft membership is retained as a mixture while unavailable evidence is varied.

For example, an uncertain TRA membership can be evaluated without refitting:

```python
from vdjdb.score.bayes import predict_marginals, score_records

soft = predict_marginals(model, [{"a_tcrnet": 0.7, "b_studies": "1"}])
scored = score_records(model, features.select(FEATURES))
```

A Hamming-1 neighbourhood can supply an exploratory local membership fraction. Keep
species, chain and pMHC context fixed, retain the producer's V/J matching convention,
count distinct annotated neighbours once, and report both member and total counts.
Use one batched sequence search with substitutions only; include the exact match and
all one-substitution neighbours. An empty annotated neighbourhood is unknown. This
fraction is not a calibrated binding probability. Preserve exact record identity and
publication counts: a neighbouring junction is sequence support, not an independent
report of the same receptor. Validate any such feature change on held-out groups before
replacing the hard-membership inputs of an accepted model.
