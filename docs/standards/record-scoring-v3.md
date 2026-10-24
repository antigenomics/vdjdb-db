# Record confidence v3.2a and v3.2b

The default scorer exports two separate probabilities for each record:

- **v3.2a, `validation_score`**: TCRvdb assay-validation probability from integrated receptor:pMHC evidence. `record_score` and `paired` remain aliases of this forecast.
- **v3.2b, `independent_validation_score`**: probability that either reported chain is also reported by another publication at the same pMHC, from focal-record evidence.

The legacy `vdjdb.score` remains the original 0-3 annotation score. These probabilities
are separate columns in the record-score tables; they do not change legacy positional files.

## Run v3.2

```bash
uv sync --extra motifs --extra confidence
uv run vdjdb score-records --tables out/tables
```

A new build uses its completed motif evidence and the existing separated producer
indices on its full cohort. No training occurs during scoring. To use declared fixed
continuous motif evidence unchanged, select it explicitly and write to another directory:

```bash
uv run vdjdb score-records --tables /path/to/fixed/tables \
  --motif-scores /path/to/fixed/record_scores.parquet --out out/v32-scores
```

The fixed input contains `record_id` and `a_tcrnet`, `b_tcrnet`, `a_tcremp`,
`b_tcremp`. Conflicting values within an interaction fail. Unmatched interactions
retain the available hard evidence. This input is an explicitly declared reference,
not a feature cache. The command refuses to overwrite that input. It does not rebuild
geometry or motifs when `--motif-scores` is supplied.

Use `--legacy-v31` for the original v3.1 artifact or `--legacy-v3` for original v3.
Both preserve their previous encodings and inference. `--hard-membership` and
`--soft-membership` remain legacy-model options; v3.2a requires the evidence encoding
on which its parameters were fitted. `--threshold` changes the v3.2a envelope decision,
not either probability. Its default is 0.5, with no new threshold optimization.
The independent-study probability has no binary decision column.

## Output and model identities

`record_scores.parquet` and deterministic `record_scores.tsv.gz` contain one row per
current `record_id`. All records are scored; the fitting subsample does not remove
database records from exports.

| Column | Meaning |
|---|---|
| `validation_score` | v3.2a assay-validation probability, [0,1] |
| `record_score`, `paired` | Aliases of `validation_score` |
| `validation_score_version` | `3.2a` |
| `model_sha256` | Exact assay-model bytes |
| `independent_validation_score` | v3.2b additional-publication probability, [0,1] |
| `independent_validation_score_version` | `3.2b` |
| `independent_model_sha256` | Exact independent-study-model bytes |
| `score_version` | `3.2`, identifying the paired score-table release |
| `TRA`, `TRB` | v3.2a conditioned only on that chain's study/motif evidence |
| `paired_lower`, `paired_upper` | v3.2a sensitivity range over legal completions of missing predictors |
| `positive` | True when lower bound >= threshold, false when upper bound < threshold; otherwise missing |

The envelope includes the nominal probability and, without direct pair evidence,
the two chain forecasts. It is a sensitivity range, not a statistical confidence
interval. `direct_pair_evidence` and `pair_supported` remain descriptive.
Both score columns use 12 significant figures on export; decisions use full precision.

## Fitted models

The frozen artifacts are `src/vdjdb/score/record_v3_2.json` and
`src/vdjdb/score/record_v3_2b.json`. The production implementation uses NumPy/Polars
inference; graph learning and EM run only in `processing/`.

For v3.2a, integrated legacy confidence is encoded as 0, 1, 2+.
An experimentally known structure sets that confidence to 2+ and ipTM to 1.0.
A predicted structure alone does not trigger this override. Count bins remain
0/1/2/3+. Exact ipTM boundaries and continuous-motif bins are frozen in the artifact.
Four motif states divide each side of the producer membership boundary at 0.5.
No bin fitting occurs in a build. Missing predictors are marginalized exactly.

Graph candidates were HC-BIC, HC-AIC, conditional-MI ARACNE-style graphs and refitted
legacy structures. The accepted model uses AIC graph proposals and explicit epitope
contexts: GLC, YLQ and TTD use local parameters; GIL uses the labeled-only global fit;
other epitopes use the selected-unlabeled global fit. Unrecognized epitopes use that
global fit. Epitope identity selects a parameter context and is required at inference.

The latent outcome is binary. Observed outcomes are clamped during Dirichlet-smoothed
MAP EM; unknown outcomes remain latent. Unknown records are never assigned manual
positive labels. The smoothing mass is 10 per node. All known labels are used for
fitting and fine-tuning, without folds. AIC/BIC diagnostics use observed-data
likelihood, marginalizing unknown outcomes. Search ties are resolved lexically.

Selected ambiguous outcomes and selected records from three additional epitopes enter
the accepted global fit with unknown outcomes. Adding all retained database records
at unit weight failed the regression gates and is not the accepted fitting policy.
The fixed motif evidence was used unchanged. The permanent 30-per-group study
manifest affects fitting/prior cohorts only; it does not alter source chunks or
rebuild the pooled motif geometry.

V3.2b uses focal-record confidence and structural evidence plus fixed hard motif
memberships. It excludes study counts, donor counts, publication IDs and
cross-publication confidence maxima. Publication IDs construct the outcome only:
at least one nonempty chain clonotype must occur in another publication at the same
`pmhc_id`. Two records from one publication do not supply replication.
“No additional report” is the other class and does not mean assay nonvalidation.
Its native-structure override is also 2+ and ipTM=1.0.

## Standalone fitting with private inputs

The standalone program follows the processing-script pattern used for
[10X assays](../tenx-dextramer.md):

```bash
uv sync --extra score-training
TCRVDB_PATH=/path/to/private/validation.csv \
uv run python processing/fit_record_score_v3_2.py \
  --tables /path/to/fixed/tables \
  --selected-records /path/to/private/selected.csv \
  --selection-manifest /path/to/private/selection.tsv \
  --motif-scores /path/to/fixed/record_scores.parquet \
  --task both --out /path/to/new-fitting-results
```

`--tcrvdb-path` overrides `TCRVDB_PATH`. There is no developer-specific fallback.
The validation CSV requires unique `tcr_id`, binary `pos` and `ambiguous`.
Ambiguous labels remain unknown. The selected CSV requires `tcr_id`, `epitope_aa`,
`hla_long`, `cdr3_alpha_aa`, `cdr3_beta_aa`, `TRAV_IMGT`, `TRAJ_IMGT`,
`TRBV_IMGT`, `TRBJ_IMGT` and `vdjdb_score`. A ZIP containing `test.csv` is also
accepted. Junctions include both conserved anchors; gene aliases are normalized
without modifying submitted junctions.

The selection TSV requires `record_id` and exactly 30 retained records from each
of the three engineered phage-display groups. It is supplied privately; the runner
does not redraw it. Other records, including target-cell assays, remain included.

Default `--methods frozen --arms selected` refits the accepted graph/bin specification.
This reproduces parameters without claiming to repeat earlier exploratory search
trajectories. `--methods hc-bic hc-aic aracne` explicitly requests new graph searches.
`--arms labeled selected all-vdjdb` requests the recorded prior alternatives.
The accepted GIL labeled-only context is explicit in the runner.
`--task independent` fits v3.2b without assay labels or selected-record input.

Private raw data is not distributed. Run outputs include model parameters,
predictions, metrics, convergence and a private input-hash manifest. Output directories
must be new. Distributed parameter artifacts contain no source paths or receptor IDs.

## Evaluation contract

Report ROC AUC, average precision (AP) and trapezoidal precision-recall area
(PR-AUC) separately on the same labeled rows. Unknown outcomes receive predictions
and are excluded from metrics. An epitope without both classes has no AUC/AP/PR-AUC
estimate. The global endpoint is the unweighted mean ROC AUC over the five labeled
epitopes, not pooled AUC.

Differences use 10,000 paired within-epitope/class-stratified bootstrap draws with
`vdjdb.config.SEED`. Recompute macro differences within each draw. A regression is
a candidate-minus-comparator 95% interval wholly below zero. Check every evaluable
epitope, historical GLC/YLQ labels separately, and the five-epitope mean. Preserve
point losses even when the interval crosses zero.

Audit aligned forecasts with the separate script:

```bash
uv run python processing/audit_record_scores.py --cases /path/to/private/cases.tsv \
  --candidate v32a --comparators v3 v31 incumbent --out /path/to/new-audit
```

The flat TSV requires `endpoint`, `epitope`, `case`, binary `label`, `predictor`
and finite `score`. Each predictor must contain identical case/label sets.
Use separate endpoint values for current and historical assay labels. Omit unknown
labels explicitly. The script writes AP, PR-AUC, ROC AUC and paired difference
intervals, including per-endpoint macro summaries; it fails if any CI upper bound
is below zero. This input and its results remain private.

These are fitted-cohort comparisons after fine-tuning on the labels. Their intervals
condition on the fitted predictions and do not establish generalization after
selection. V3.2b is evaluated against corpus-derived additional-report labels;
the fixed motif parameters had previously been tuned using replication evidence.
Its evaluation does not constitute independent assay validation.

## Historical v3.1 and original v3 specification

The remainder describes the preserved v3.1/v3 models and their original evaluations.
Its default and training statements apply to those explicit legacy modes.


Record confidence v3.1 is the default. Its frozen model is fitted on the released
positive/negative outcomes for five epitopes, with 1,782 labelled interactions.
The 252 ambiguous outcomes receive scores and are excluded from binary metrics.
The original 614-interaction GLC/YLQ model remains available as `--legacy-v3`.

The model combines publication, donor, motif and predicted-structure evidence
in one frozen Bayesian network. It returns a paired forecast and two chain-only views
of validation of the same receptor:pMHC interaction. The existing `vdjdb.score` remains
the legacy 0-3 score; v3 is a separate table linked by `record_id`.

![Record scoring pipeline](../_static/record-score-v3/workflow.svg)

**Figure 1. Record scoring pipeline.** Current records and completed chain-level motif
and record-level structural evidence are integrated at gene-level receptor:pMHC identity.
One fitted network supplies all three forecasts. Missing evidence is marginalized;
missing chain and structural evidence also defines a sensitivity envelope. The result
is joined back to every observational record ID, including independent reports of the
same interaction. Model parameters are inputs; no fitting occurs during a build.

## Use the score

Choose `paired` for the available combined evidence, or `TRA` / `TRB` for a chain-only
view of the same paired assay endpoint. Use the sensitivity envelope and `positive`
together when a definite call is needed. See [missing evidence and decisions](#missing-evidence-and-decisions)
before filtering records.

### Run the pipeline

```bash
uv sync --extra motifs --extra confidence --extra test
uv run vdjdb build --out out/ --structure-metadata input-metadata.tsv.gz
uv run vdjdb motifs --tables out/tables --out out/motifs
uv run vdjdb score-records --tables out/tables --legacy-v31
```

The metadata argument supplies the reviewed structural input described in
[build instructions](../builds.md). A build without it can still be scored; unavailable
ipTM is marginalized. Scoring requires the completed motif stage and refuses a table
with no motif evidence. CI runs scoring after motifs with packaged v3.1. `--legacy-v3` selects the unchanged
original model and hard motif evidence. `--hard-membership` selects hard evidence
with the chosen model; `--threshold` overrides its decision cutoff without changing
probabilities. For example:

```bash
uv run vdjdb score-records --tables out/tables --legacy-v3 --out out/legacy-v3
uv run vdjdb score-records --tables out/tables --legacy-v31 --hard-membership --out out/v31-hard
```

`--model` selects an explicitly supplied frozen JSON artifact; `--out` selects a separate output
directory. Default output is the input table directory.

`record_scores.parquet` and `record_scores.tsv.gz` contain exactly one row per current
`record_id`, sorted by ID. TSV has a deterministic gzip header; unknown values are empty
cells. Parquet retains typed nulls. No column is appended to a legacy positional file.
The six probability/envelope columns use 12 significant figures in both formats.
Inference preserves the fitted network's node and factor order after pruning, so
Python hash seeds cannot change tensor contraction order.
Decisions use full-precision inference;
input evidence and fitted parameters are unchanged.
The primary bundle includes both score files when that stage has run; the hosted build
requires their production and verifies their cardinality. The model SHA256 is recorded
on every row. No stored feature or prediction rows are reused between builds.
Normal hosted builds retain both files in the `record-confidence-v3` artifact for 90 days.

| Columns | Type / interpretation |
|---|---|
| `record_id` | Current observation identifier; unique primary key |
| `a_studies`, `b_studies` | Distinct publication counts for the TRA/TRB chain, peptide and V gene |
| `donors` | Distinct nonempty publication/donor contexts for the exact integrated receptor:pMHC |
| `a_tcrnet`, `b_tcrnet`, `a_tcremp`, `b_tcremp` | Separated producer indices in v3.1; hard membership in legacy mode; missing if unannotated |
| `*_hard`, `*_support` | Original binary assignments and bounded producer support in continuous mode |
| `iptm` | Mean ipTM on [0,1] over distinct accepted predicted models of that interaction |
| `legacy_score` | Maximum built legacy score over the integrated interaction |
| `paired`, `record_score` | Identical probabilities from all available evidence, on [0,1] |
| `TRA`, `TRB` | Probabilities conditional only on that chain's publication and motif evidence |
| `paired_lower`, `paired_upper` | Sensitivity envelope on [0,1], including the nominal paired forecast |
| `direct_pair_evidence` | An exact-pair reference annotation, donor support or accepted predicted structure |
| `pair_supported` | Direct pair evidence or positive support for both chains; descriptive |
| `positive` | True if lower limit ≥ the decision cutoff, false if upper limit is below it, missing otherwise |
| `decision_threshold` | Global cutoff used by v3.1; 0.26 in the packaged model |
| `score_version`, `model_sha256`, `soft_membership` | Version `3.1` by default or `3` for legacy, exact frozen model identity and the inference mode |

Publication counts use the established chain/peptide/V-gene grouping, not a count of
independent experiments on the paired construct. Donor labels are scoped to their
publication. Repeated annotations of one structural model contribute once to mean ipTM.
Native experimental structures do not contribute predicted ipTM and were excluded from
model fitting and evaluation. The source feature definitions match the accepted model.

## Model and parameter defaults

V3.1 uses `vdjdb.score/record_v3_1.json`: the same independent-ipTM-leaf network
architecture, fitted on the five-epitope labelled cohort. Equivalent sample size is
10. Training-only three-grade ipTM boundaries are 0.83 and 0.87. Original producer memberships define the fitted binary states;
separated indices provide narrow interpolation at inference. No epitope identity
is a predictor. Missing evidence remains marginalized.

The global decision cutoff is **0.26**, rounded to two significant digits from
the F1 optimum on all grouped out-of-fold forecasts. The nested evaluation below learns its global
cutoff separately within each outer training fold. Per-epitope F1 policies are
research outputs; production uses one scalar cutoff. `--threshold 0.5` is an
explicit override. The cutoff changes decisions and leaves probabilities, AP
and ROC-AUC unchanged.

### Legacy model and defaults

The legacy deployment fit uses 614 labelled interactions, nine predictors and sixteen directed
edges. Validation is a parent of every predictor. The nonstructural predictors also
form a conditional dependence tree; ipTM has only validation as its parent and no
children. Removing unavailable ipTM therefore preserves the fitted chain-evidence
network. No epitope name, MHC allele, sequence embedding or outcome is a predictor.
The endpoint is the measured assay validation criterion `padj < 1e-5`.

| Setting | Accepted deployment value | Other implemented or assessed settings |
|---|---|---|
| Model | Packaged `vdjdb.score/record_v3.json`; one all-labelled fit | `--model` explicitly selects another artifact; no automatic model search |
| Dependence graph | `tan_independent_iptm` | Generic `fit()` supports `tan`, `tan_leaf`, `naive`; these do not replace the packaged fit |
| Dirichlet equivalent sample size | 10 for every node | `fit()` accepts a positive scalar or complete node map; stronger smoothing was not promoted |
| Count encoding | 0, 1, 2, 3+ | Fixed encoding; no build-time bin fitting |
| ipTM encoding | Low below ≈0.74; mid below ≈0.84; high otherwise | Exact boundaries are saved in the [packaged model](https://github.com/antigenomics/vdjdb-db/blob/master/src/vdjdb/score/record_v3.json); generic `categorical()` defaults to 0.5/0.8 unless supplied |
| Motif membership | Hard, producer-owned membership | `soft_membership=False` throughout the inference API; fractional memberships require explicit `True` |
| Missing evidence | Exact marginalization, including legal states absent from training | No zero imputation or complete-case exclusion |
| Decision threshold | 0.5, fixed | No per-epitope threshold tuning during deployment |
| Prior / temperature | Saved fitted prior; no post-fit temperature correction | Prior 0.5/0.7/0.9 and temperature comparisons were diagnostics, not deployment options |
| Random seed | `vdjdb.config.SEED`, stored as 20260925 | Inference does not sample |

`fit()` is a generic training API whose default graph is `tan`; calling it is not the
production pipeline. Production reads the accepted parameters without fitting. Any
new fitted settings require their own grouped validation before selection.

Explicit soft membership mixes the nonmember/member posteriors with supplied weights `1-q` and
`q`. Several uncertain memberships use an independent mixing distribution, while the
BN retains its fitted feature dependencies. Hard endpoints are unchanged. Direct
Hamming-1 membership-fraction substitution reduced paired AP from 0.911 to 0.909
for GLC and from 0.869 to 0.868 for YLQ, and lost two correct GLC TRB calls.
It remains off. `score-records --soft-membership` explicitly selects this exploratory
variant. It performs one substitutions-only batched native search, includes exact
neighbours, retains each producer's species/chain/peptide and V/J scope, and adds
`*_members` / `*_neighbours` integer columns for each motif fraction. Publication
counts and record identity stay unchanged. Legacy hard scoring runs no neighbourhood search.

### Separated producer indices

V3.1 uses separated producer indices by default. `--continuous-membership`
explicitly selects the same mode for another model, including legacy v3.
It cannot be combined with `--soft-membership`. For hard membership $b\in\{0,1\}$,
bounded support $u\in[0,1]$ and width $w=1/12$, the index is

$$s=(1-w)b+wu.$$

Nonmembers score in $[0,1/12]$ and members in $[11/12,1]$, with a minimum
between-class gap of $5/6$. Thresholding a chain index at 0.5 exactly recovers
its original binary assignment. Unannotated inputs stay missing. Original
`*_hard` assignments determine pair-support flags.

TCRNET support is $u=d/(d+1)$, where $d$ counts other distinct annotated
Hamming-1 junctions in the producer's species/chain/peptide scope. The native
search includes the exact reference hit, which is removed from this degree.
TCREmp support is $u=e/(e+d_2)$, where $e>0$ is its pooled reference DBSCAN
radius and $d_2$ is distance to the second retained core. This supplies density
variation for an exact self-hit. Geometry is rebuilt from the reference;
query outcomes never alter it. Eligible groups without retained cores have
support zero. The original hard-assignment eligibility remains authoritative.

Continuous output records `continuous_membership=True`,
`motif_score_version=separated-producer-v2`, `motif_score_width`, `*_support`,
`*_hard` and neighbourhood counts. `soft_membership=True` identifies the
posterior-mixture inference used by both fractional modes. These indices are
interpolation weights, with no claim of calibrated membership probabilities.
The chain-index cutoff of 0.5 is separate from the final forecast's decision cutoff.

`score.evidence.producer_motifs` implements this integration. The earlier
Hamming-1 member-fraction index remains available through `separated_motifs`.
`score.membership.separated(binary, support, width=...)` accepts any finite
$0\le w<0.5$; its default width remains 1/3 for earlier research callers.
Zero width recovers hard evidence wherever both inputs are present.
`motifs.tcremp.ClusterReference.support` supplies bounded core geometry for
standalone research. Numerical reproduction follows the
[clustering runtime contract](../clustering.md).

`validate.membership` provides training-only F1 cutoff selection, a monotone paired
formula, and ranking metrics without reading an assay. Its paired formula is
$\operatorname{sigmoid}(c+w_a s_a+w_b s_b+w_{ab}s_a s_b)$, with nonnegative
weights, an optional interaction and a fixed L2 penalty of 0.001 on the weights.
Here $s_a,s_b$ are the TRA/TRB indices, $c$ is the fitted intercept, and
$\operatorname{sigmoid}(z)=1/(1+\exp(-z))$. Coefficients and cutoffs require
grouped training data; they are not packaged production parameters.

AP uses step integration and PR-AUC uses trapezoids, reported separately.
Normalized partial AP at recall $R$ is
$R^{-1}\int_0^R P(r)\,dr$, using AP's constant precision over each tied score
group's recall increment. Precision at recall $R$ uses the first attainable
threshold reaching that recall, keeping ties together. ROC-AUC and standardized
partial ROC-AUC at FPR 0.1 are separate metrics. Macro averages give each peptide
equal weight. A global or peptide-specific decision threshold changes F1 and
confusion counts but does not change these ranking metrics. Cutoffs and candidate
selection must use training folds, with shared receptor groups kept together.
`fit_thresholds` fits peptide cutoffs only where training rows include both
outcomes. `apply_thresholds` uses the learned global cutoff for an unseen peptide
or one with a single training outcome. The research runner records fallback use.

## Missing evidence and decisions

TRA and TRB condition on the respective chain's three predictors. They predict the
paired assay endpoint, not a single-chain binding experiment. Do not multiply them
with each other or the paired probability: these forecasts share a model and evidence.

The sensitivity envelope completes missing chain variables and ipTM with states
observed during training, considering their joint combinations. Metadata stays
marginalized and the nominal forecast is included. Without direct pair evidence,
the envelope also includes both chain-only forecasts. Unreported or unsupported
states do not create positive evidence. The envelope is not a credible interval.

Missing pair evidence never sets `record_score` to zero. Neither support flag gates
the forecast. `positive` is a separate conservative decision: an envelope crossing
the selected global cutoff yields unknown. Report the unknown count alongside definite calls. An absent
partner or an unrelated structural model cannot supply direct evidence for a new pair.

## Five-epitope assessment

The new endpoint uses the released positive/negative labels. Transitive groups
keep every shared TRA or TRB junction in one outer fold. Five outer and three
inner grouped folds separate fitting, candidate selection and threshold selection.
Reference geometry is fixed and transductive. The single-chain views predict the
paired assay endpoint. Scores are retained for all 2,034 rows, including 252
ambiguous outcomes; binary metrics use 1,782 labelled rows.

| Epitope | Labelled interactions | AP | ROC-AUC | F1 with training-selected global cutoff |
|---|---:|---:|---:|---:|
| GILGFVFTL | 575 | 0.937955 | 0.895412 | 0.873536 |
| GLCTLVAML | 166 | 0.884662 | 0.787500 | 0.788136 |
| NLVPMVATV | 317 | 0.759791 | 0.898036 | 0.717949 |
| TTDPSFLGRY | 342 | 0.837451 | 0.807832 | 0.637993 |
| YLQPRTFLL | 382 | 0.857628 | 0.903735 | 0.802198 |
| Equal-epitope macro mean | 1,782 across 5 epitopes | 0.855498 | 0.858503 | 0.763962 |

Producer interpolation was selected in each outer training fold. Matched hard
refits give macro AP 0.851150 and macro ROC-AUC 0.856594. Macro F1 is 0.739737
at a fixed 0.5 cutoff, 0.763962 with training-selected global cutoffs and 0.793708
with training-selected per-epitope cutoffs. Threshold policies leave ranking
metrics unchanged. The final deployment cutoff uses all out-of-fold forecasts;
the table reports the nested training-policy assessment.

The original measured cohort and endpoint remain a separate legacy comparison.
On those 193 GLC and 421 YLQ cases, the new grouped forecasts have AP 0.880402
and 0.871931, respectively, versus legacy 0.911197 and 0.869028. The GLC reduction
is retained explicitly; the original model and its historical results remain
available. Released labels and the original adjusted-P endpoint use different
outcome definitions and denominators. Full predictions, source hashes, coverage,
cutoff policies and TRA/TRB/paired assessments remain private.

## Legacy performance with continuous and binary scores

The accepted legacy comparison uses five shared-junction grouped held-out folds, with all
193 GLC and 421 YLQ interactions retained for each predictor. AP is step-integrated
average precision, including tied scores. Accuracy uses the recorded classification
rules; ipTM decisions use the original training-fold cutoff, and motif decisions use
membership. V3 uses its fixed 0.5 threshold. No cutoff is reselected here.

Balanced accuracy is suitable for binary decisions and unequal class counts:

$$\mathrm{BA}=\tfrac12\left(\frac{TP}{TP+FN}+\frac{TN}{TN+FP}\right).$$

TP and FN count validated interactions; TN and FP count nonvalidated interactions.
Thus BA averages sensitivity and specificity. AP also accepts discrete scores; ties
are evaluated together. The [scikit-learn definition](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.balanced_accuracy_score.html)
provides the metric contract. The exported comparison additionally reports Matthews
correlation coefficient, FP/(FP+TN) and FP/(TP+FP).

| Epitope | Interactions (validated / nonvalidated) | Predictor | AP | Correct calls / interactions | Balanced accuracy |
|---|---:|---|---:|---:|---:|
| GLC | 193 (134 / 59) | ipTM | 0.882 | 146 / 193 | 0.720 |
| GLC | 193 (134 / 59) | TRA TCREmp | 0.824 | 137 / 193 | 0.739 |
| GLC | 193 (134 / 59) | V3 paired | **0.911** | **149 / 193** | **0.746** |
| YLQ | 421 (172 / 249) | ipTM | 0.683 | 300 / 421 | 0.691 |
| YLQ | 421 (172 / 249) | TRA TCREmp | 0.763 | 362 / 421 | 0.848 |
| YLQ | 421 (172 / 249) | V3 paired | **0.869** | **367 / 421** | **0.873** |

TRA TCREmp is the strongest individual motif baseline on these three metrics in each
epitope. [The complete 28-row comparison](../_static/record-score-v3/metrics.tsv) includes
TRA/TRB/both/any TCRNET and TCREmp membership, publication counts and ipTM. V3 exceeds
every supplied individual baseline on AP, accuracy and balanced accuracy in both epitopes.

![Grouped held-out score comparison](../_static/record-score-v3/performance.svg)

**Figure 2. Grouped held-out validation performance.** (a) GLC, 193 interactions;
(b) YLQ, 421 interactions. Orange circles show ipTM, green squares TRA TCREmp membership,
and blue filled circles v3. AP, accuracy and balanced accuracy share a 0-1 scale.
The points summarize the same cases in each epitope; there are no error bars or
additional significance tests. Full confusion counts are in the linked table.

`vdjdb score-metrics predictions.parquet --out metrics.tsv` checks exact case, fold and
outcome agreement before comparing predictors. Its input columns are `case`, `epitope`,
`predictor`, `fold`, `y`, `score`, `predicted`. It fails if the `v3` predictor regresses a
supplied baseline on AP, accuracy or balanced accuracy. This command evaluates saved
held-out forecasts without refitting. Assay-level identities and inputs remain private.
CI checks the reviewed aggregate comparisons, packaged model identity and fixed synthetic
forecast probes; it does not claim to rerun the private assay validation.

## False positives and negative-source calls

At the fixed nominal decision threshold, the assay confusion counts are:

| Epitope | TP | FP | FN | TN | False-positive fraction FP / (FP + TN) | False-discovery fraction FP / (TP + FP) |
|---|---:|---:|---:|---:|---:|---:|
| GLC | 109 | 19 | 25 | 40 | 19 / 59 (32.20%) | 19 / 128 (14.84%) |
| YLQ | 151 | 33 | 21 | 216 | 33 / 249 (13.25%) | 33 / 184 (17.93%) |
| Both | 260 | 52 | 46 | 256 | 52 / 308 (16.88%) | 52 / 312 (16.67%) |

The separate source-negative assessment uses all 37 files in `chunks_negative/`, with
repeated source observations retained and no fitting on these rows:

| Scope | Source observations | Nominal positive | Confident positive | Confident negative | Unknown |
|---|---:|---:|---:|---:|---:|
| All files | 6,062,284 | 2,155,249 (35.6%) | 3,747 (0.0618%) | 16 | 6,058,521 |
| Human class I with specified context | 5,015,831 | 1,116,387 (22.3%) | 3,740 (0.0746%) | 6 | 5,012,085 |

Percentages divide by the observation count in the same row. Confident calls use the
envelope, while nominal calls threshold `record_score`. The source-negative positive
fraction is not an estimate of the false-record prevalence in positive VDJdb records.
The model avoids confident acceptance for most of this collection by returning unknown:
definite-call coverage is only 3,763 / 6,062,284 (0.0621%). It does not establish those
unknown observations as negatives. The assumed 10-50% false-record prevalence explored
previously is a prior scenario, not a measured corpus error fraction.

## Python API and implementation

```python
from vdjdb.emit.vdjdb3 import read_table
from vdjdb.score.evidence import integrate
from vdjdb.score.pipeline import load_model
from vdjdb.score.bayes import score_records

records = read_table(tables, "records")
chains = read_table(tables, "chains")
evidence = read_table(tables, "evidence")
features = integrate(records, chains, evidence)
model, model_sha256 = load_model(version='3.1')
forecasts = score_records(model, features)
```

`integrate()` owns feature definitions, identity and aggregation. `interactions()`
creates canonical record-to-interaction keys. `bayes` owns fitting, exact inference,
chain views, optional soft inference and missing-evidence ranges. `pipeline.record_scores()`
scores distinct interactions once, then restores every record ID through a validated
many-to-one join. `pipeline.run()` reads completed tables and writes both compressed
formats. Existing Polars joins and batched inference suffice; there is no worker pool.

For explicit soft inference use `score_records(model, numeric_frame,
soft_membership=True)` or the same keyword on `predict`, `predict_marginals` and
`predict_evidence`. Fractions are rejected by default. Omitted evidence columns are
unobserved, and extra identity columns are ignored by numeric inference. Count and
ipTM validation reject invalid ranges and noninteger counts.

For the separated-neighbour option, the pipeline preserves the original hard
assignments as well as the experimental continuous features:

```python
from vdjdb.score.pipeline import record_scores

optional = record_scores(
    records, chains, evidence, model,
    model_sha256=model_sha256, continuous_membership=True,
)
```

Standalone research cutoffs use a separate training-only policy:

```python
from vdjdb.validate.membership import fit_thresholds, apply_thresholds

policy = fit_thresholds(training_labels, training_scores, training_peptides)
cutoffs = apply_thresholds(policy, test_peptides)
predicted = test_scores >= cutoffs
```

The training scores must be out-of-fold predictions if their model was fitted
on these labels. Keep related receptor observations in the same fold. These
cutoffs do not change the production model's envelope decision at its cutoff (v3.1: 0.26; legacy: 0.5).

## Limits of the evidence

These two-epitope results are grouped development evidence after model refinement,
not an untouched external estimate of the complete selection procedure. They establish
assay validation discrimination, not affinity, function or clinical response. Native
experimental structures are excluded. CIN/NLV have scored examples without measured
validation endpoints, so accuracy there is unmeasured. Cross-epitope calibration is
not established. The chosen fit preserves accepted TRA/TRB/no-ipTM marginals; its YLQ
paired AP is 0.00198 below an earlier paired-only candidate whose marginal regressions
prevented selection. Confidence envelopes describe sensitivity to unavailable evidence
and do not guarantee frequentist coverage. Source-negative labels and withdrawn 10X
assignments are distinct; a discarded assignment is not an experimentally negative TCR.
