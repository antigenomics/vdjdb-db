# Motif denoising and tuning

**Status: normative.** This is the specification the motif stage is built and judged against. When
this document and a measurement disagree, the measurement is reported and this document is revised.

## Sources

Two papers define the method.

1. **Shugay M, Luppov DV, Vlasova EK, Chudakov DM.** *Towards high-quality large-scale T cell
   receptor antigen specificity data: challenges and promises.* **Nature Methods**, 3 September 2026.
   [10.1038/s41592-026-03227-2](https://doi.org/10.1038/s41592-026-03227-2) · PMID 42693187.
   *(Record verified against PubMed.)*
2. **Luppov DV, Koneva AE, Bagaev DV, Alexandrova AV, Vlasova EK, Chudakov DM, Motozono C,
   Sewell AK, Shugay M.** *VDJdb in 2026: boosting T-cell receptor recognition evidence using
   paratope embeddings and AI-based structure prediction.* **Nucleic Acids Research**, 2026,
   Database Issue. *In press - DOI assigned at production.*

Two further documents are cited where they are used: the IMMREP25 audit
(`~/vcs/manuscripts/2026-immrep25-audit`) for the calibrated cohort ladder, and the validation
manuscript (`~/vcs/manuscripts/2026-vdjdb-db-validation`) for the held-out check.

---

## 1. False positives in assay-derived records

The motif stage starts from the assumption that a substantial fraction of assay-derived records are
false positives. Multimer and screen-derived records include bystanders and mis-assigned
specificities, so a record's presence in the database is not evidence that the receptor binds the
epitope it is filed under.

VDJdb therefore reports, per record, how much of the database's own structure supports it.

## 2. Motif membership as the noise filter

> Sequences not associated with any statistically significant motif are flagged as likely false
> positives.

That definition was established for TCRNET and carried forward to the embedding-based annotation. A
motif is several independent receptors arriving at the same paratope solution for the same pMHC. A
record that no other record resembles has nothing supporting it beyond the single assay that
produced it.

Three consequences for the code:

1. `in.motif` is a per-record evidence column, not a display feature. The validation manuscript
   scores it on the same axis as the confidence score and the independent-study count. Emit it for
   every record, including the records that get `false`.
2. A cluster is judged on being right, not on being large. Coverage is a reported quantity, never a
   maximised one (§5).
3. Absence of a motif is a weak signal, not a verdict. Low-generation-probability receptors
   (autoimmunity-associated variants are the named case) and under-sampled epitopes are motif-less
   for reasons unrelated to being wrong. The flag is advisory and never a deletion.

## 3. The noise signature

The computable read on which records are spurious:

| CDR3 is associated with… | …and sits in | reading |
|---|---|---|
| one or a few epitopes | a motif | specific, supported |
| several epitopes | a few motif contexts | cross-reactivity |
| an unexpectedly high number of epitopes | no motif at all | false positive |

The third row is the actionable one: the noise signature is promiscuity together with
motif-lessness, and neither alone. A CDR3 in many epitopes that still sits in coherent motifs is
cross-reactive; a CDR3 in many epitopes with no motif anywhere is a sequence that keeps turning up
in assays.

Do not use promiscuity as a filter on its own, and report the epitope count and the motif assignment
together.

## 4. TCRNET sensitivity limits

TCRNET is highly stringent - a Hamming-1 CDR3 graph, with limited handling of CDR1 and CDR2 - and
that stringency restricts sensitivity. Two named failure modes:

- **Fragmentation.** One motif is split across clusters that differ only in CDR3 length. The
  documented case is HLA-A\*02:YLQPRTFLL TRB, where TCRNET clusters C21 and C24 have highly similar
  length-16 logos and stay separate; the embedding-based annotation merges them into one cluster
  (C84) spanning several raw CDR3 lengths, on the shared CSAR/GELFF-16 motif.
- **CDR1/CDR2 blindness.** V-gene identity accounts for two thirds of the paratope. A CDR3-only
  distance cannot see it, so V-gene bias (TRBV20-1 for HLA-A\*02:GLCTLVAML, TRAV12-2 for
  HLA-A\*02:LLWNGPMAV and HLA-A\*02:ELAGIGILTV) is invisible to the graph that defines TCRNET's
  balls.

The embedding aggregates V segments with variable CDR3 lengths, which gives larger, less fragmented
clusters and higher coverage on both chains.

An epitope may have several disjoint motifs; HLA-A\*02:NLVPMVATV is the named case. Do not assume
one motif per epitope, and do not score a clustering by how close it comes to one cluster per
epitope.

## 5. Coverage and neighbourhood width

Fragmentation is a failure. Percolation is a failure only where it fuses *different* motifs: an
epitope with one prominent motif (A\*02 GIL) should have one dominant cluster, and a featureless one
(A\*02 NLV) should have many, so the largest cluster's share is a property of the epitope and is
recorded, not gated (the note on percolation at the top of `docs/clustering.md`). Fixing fragmentation by fusing motifs is not progress,
and purity is what shows it.

Widening the similarity ball raises the fraction of clustered records monotonically, by recruiting
the bystanders and mis-assigned records of §1. Once those sit inside a named motif with a logo they
read as epitope-specific signal, so every aggregate coverage metric improves while the database gets
worse.

Measured on this corpus: TRA connected components at scope `1,0,0,1` → `2,0,0,2` → `3,0,0,3` →
`4,0,0,4` raise the clustered-record fraction 0.2105 → 0.4149 → 0.4692 → 0.4807, while the
enrichment of clustered clonotypes for independent replication falls from 1.96× to 1.30× at the
second step alone.

> **Rule.** Never rank motif configurations on coverage, retention, or the clustered fraction.
> Report them; do not optimise them.

### 5.1 The inflation mechanism

The combinatorial bound and the measured inflation differ by an order of magnitude.

The set of peptides within $s$ substitutions of a CDR3 of length $L$ has size

$$V_s(L)=\sum_{j=0}^{s}\binom{L}{j}\,19^{j},$$

so for $L=14$: $V_1=267$, $V_2=33{,}118$, $V_3=2{,}529{,}794$. Scope 2 is a 124-fold larger ball
than scope 1. If sequence space were uniformly occupied, the chance neighbour rate would inflate by
that factor.

The measured inflation is 12-fold. Measured on human TRB, 124,819 scored clonotypes, background
$M=10^{6}$: the median $\hat p_\theta=(n_{\mathrm{control}}+1)/(M+1)$ is $1.00\times10^{-6}$ at scope 1 and
$1.20\times10^{-5}$ at scope 2 (means $3.48\times10^{-6}$ and $4.28\times10^{-5}$), a 12-fold
inflation against a 124-fold bound. Repertoire space is concentrated rather than uniform, so ball
volume overstates the null by an order of magnitude. Never substitute $V_s$ for a measured
$\hat p_\theta$.

Chance recruitment through the degree test is not the mechanism. With
$D\sim\mathrm{Bin}(n_e-1,\hat p_\theta)$ and recruitment at $D\ge d_{\min}$, the chance-recruitment
rate is $\alpha_\theta(n_e)=\Pr[D\ge d_{\min}]$. Computed on all 87 human TRB epitopes with $\ge100$
clonotypes: $\alpha\le0.0011$ at scope 1 and $\alpha\le0.094$ at scope 2 (mean $0.0027$), the
maximum falling on VEALYLVCG. The degree test does not become vacuous at scope 2. ($\alpha$ still
depends on $n_e$, so a fixed $p$ is not a fixed noise level, though the dependence is small at these
values.)

Untested neighbour recruitment is not the mechanism either: Stage I admits a clonotype for being
within scope of an enriched one, without its passing any test, so a denser ball would be expected to
admit more untested members, and instead untested members fall from 15.9 % of the clustered set at
scope 1 to 12.7 % at scope 2.

The mechanism is that the test's own operating point moves. The same $p<0.05$ passes 35,775
clonotypes at scope 1 and 54,339 at scope 2, a 52 % increase, because within-sample degree grows
faster with the ball than the background expectation does. Those extra passes are passes of the
test, and they are enriched for bystanders, so the independent-replication lift falls from 1.96× to
1.30× across the same step.

Tightening the threshold does not recover it: at scope 2, $p<0.01$ gives lift 1.35 against 1.30 for
$p<0.05$. A threshold cannot compensate for the scope, because the $p$-value is already computed
against a scope-matched background; widening the ball changes what the test is about, not how strict
it is.

`vdjdb.validate.noise.ball_volume` and `chance_recruitment` compute the first two quantities.
$\phi_e$, the share of epitope $e$'s records that bind it, is unknown by construction, so
$\mathbb{E}[F_e]=(1-\phi_e)\,n_e\,\alpha_\theta(n_e)$ is used as a sensitivity curve over plausible
$\phi$ and never as a point estimate.

### 5.2 Interpretation of the lift

Write $T$ for "record is a true binder" with prevalence $\phi=\Pr(T)$, $M$ for "in a motif", $R$ for
"independently replicated", and $r_1=\Pr(R\mid T)$, $r_0=\Pr(R\mid\neg T)$. If $M$ and $R$ are
conditionally independent given $T$, both following from a receptor's being a convergent solution,
then

$$\Pr(R\mid M)=r_1\Pr(T\mid M)+r_0\bigl(1-\Pr(T\mid M)\bigr),\qquad
\Pr(R)=r_1\phi+r_0(1-\phi),$$

and when $r_0\ll r_1$, which states that a mis-assigned record is not independently re-reported
except by coincidence, the lift reduces to

$$\mathrm{lift}=\frac{\Pr(R\mid M)}{\Pr(R)}\;\approx\;\frac{\Pr(T\mid M)}{\phi}.$$

The independent-study lift therefore estimates how much motif membership enriches for true binders,
which is the quantity the denoising claim needs, and it inverts to
$\Pr(T\mid M)\approx \mathrm{lift}\times\phi$.

⚠ **The conditional-independence assumption can fail through publicity, so check it rather than
assume it.** A high-generation-probability clonotype is more likely both to be re-observed by a second
laboratory and to have sequence neighbours, so part of a raw lift can be $P_{\mathrm{gen}}$ rather
than specificity; the IMMREP25 audit finds this confound dominant in a weak cohort.

The control is the audit's own: permute $R$ within strata of $\log_{10}P_{\mathrm{gen}}$, which
preserves each clonotype's publicity and destroys only its association with the clustering, then
report

$$\mathrm{lift}_{\mathrm{ctrl}}=\frac{\mathrm{lift}}{\mathbb{E}[\mathrm{lift}\mid\text{permuted}]}.$$

A ratio near 1 means the raw lift was the covariate.

On the shipped 2026-06-03 annotation, the within-$P_{\mathrm{gen}}$ null sits at 1.043 (TRA) and
0.907 (TRB) against raw lifts of 1.721 and 1.355, giving controlled ratios of 1.650 and 1.494
($p<0.005$, 200 permutations). Publicity explains almost none of the enrichment here, so on this
corpus the raw lift can be read as specificity enrichment. That is a measurement on one corpus and
not a general licence: recompute it whenever the clustering or the corpus changes.

Because the confound is small here, $\Pr(T\mid M)\approx\mathrm{lift}\times\phi$ is usable. Do not
quote $\phi\le1/\mathrm{lift}$ as a bound on how much of VDJdb is wrong: it inherits every remaining
assumption ($r_0\ll r_1$ in particular), and the scored cohort is not the database.

Both are computed by `vdjdb.validate.noise.controlled_lift`, with `pgen_stratum` supplying the
strata from the `cdr3nt.pgen` column the build already generates. Report the raw lift and the
controlled lift together, never the raw one alone.

## 6. The three tuning instruments

| | question it answers | when it is used |
|---|---|---|
| §11.1 independent-study lift | does motif membership predict a record an independent laboratory also reported? | picks the shipped configuration |
| Q = 2hp/(h+p) | is the motif structure sound, neither fragmented nor percolated? | compares algorithms and releases |
| MATCHMAKERS / TCRvdb | does the denoising recover functionally validated records? | held out, touched once, aggregate only |

### 6.1 The tuning objective

The signal is how many distinct publications report the same clonotype against the same epitope. A
clustering that finds convergent selection recovers those clonotypes.

It is scored as F1 of `clustered` predicting `replicated`, with the lift over the base rate reported
beside it (`vdjdb.motifs.tcremp._objective`). F1 and not recall, for the §5 reason: clustering
everything wins recall and loses precision. Precision here is the fraction of the clonotypes put in
motifs that a second laboratory independently reported, and bystanders are by construction not
independently replicated.

Quote lift with the base rate, which is ~2.2 % of clonotype-epitope pairs, so that a precision of
0.15 is a 6.8× enrichment and means nothing stated alone.

#### The non-display denominator

A display-selected record - `method.identification` containing `display` - is not an independent
natural observation: a library panned against one pMHC yields thousands of receptors one
substitution apart by construction, and the panning experiment is a single `reference.id`. Every
display clonotype therefore contributes zero independently-replicated pairs while occupying a
denominator slot.

Measured on human TRB, 2026-09-26:

| Cohort | Clonotype-epitope pairs | Replicated | Base rate |
|---|---:|---:|---:|
| all | 116,053 | 2,771 | 2.3877 % |
| excluding display | 86,361 | 2,771 | 3.2086 % |

The replicated count is identical in both rows: all 2,771 replicated pairs sit outside the display
set, and display contributes 29,692 pairs that can only ever be false positives for this objective.
Scoring lift on the full cohort deflates precision by an amount that depends on how aggressively a
given parameter clusters the display block rather than on whether it finds convergent selection.

Lift is therefore scored on the non-display subset, as `vdjdb.motifs.tcremp.fit_coef` documents. Two
consequences:

- Clustering still runs on everything. The hold-out is the *scoring* denominator, not an exclusion
  from the data. Whether display records should also be held out of the clustering is a separate,
  open question (ROADMAP §30.6).
- A lift figure is not comparable across denominators. The same clustering reads 5.31× on the
  non-display cohort and a much lower number on the full one. State which cohort a lift is on
  wherever it is quoted in this repository.

The tuning sweeps report both conventions side by side.

### 6.2 The structural instrument Q

`Q = 2hp/(h+p)` is the homogeneity–parsimony trade-off (Tiffeau-Mayer, arXiv:2607.20799), taken from
the authors' `clustereval`. Homogeneity `h` asks whether clusters predict the epitope; parsimony `p`
asks whether they do it without shattering each epitope into singletons. Both are normalised by the
maximum attainable for the label partition, so epitope count and size divide out.

It scores §4's failure and §5's failure with one number:

- shatter an epitope into singletons → `h = 1`, `p = 0`, `Q = 0`
- percolate it into a single cluster → `p = 1`, `h = 0`, `Q = 0`

It replaces the purity/retention pair, two numbers that trade against each other and can be moved in
opposite directions. Computed by `vdjdb.validate.qscore`.

⚠ **Q is keyed on the clonotype, never on `(epitope, clonotype)`.** Production runs motif detection
inside one epitope's record set, so no other epitope's records are ever in the input and `h` is 1.000
by construction there. Keyed on the clonotype, it measures whether a CDR3 filed under two epitopes
ends up in one motif.

### 6.3 The held-out check

MATCHMAKERS / TCRvdb is read only via `$VDJDB_TCRVDB`, enforced by `vdjdb.validate.guard`. Never tune
on it, never ship, commit or redistribute it, and report aggregate metrics only, never per-record
verdicts. Two epitopes is too narrow a basis for a global hyperparameter, and tuning on it would
remove the only independent check this project has.

## 7. The acceptance bar and the baseline

The bar, from the published validation of the current annotation: precision maintained, recall and F1
improved. A configuration that trades precision for recall has not met the bar; it is a different
operating point and is reported as one.

The baseline is measured on the shipped 2026-06-03 `cluster_members.txt`; these are the numbers a
rebuild has to beat, and they were not computed during the purity-era tuning.

| chain | lift | controlled | F1 | precision | recall | clustered | base rate |
|---|---:|---:|---:|---:|---:|---:|---:|
| TRA | 1.721 | 1.650 | 0.0849 | 0.0469 | 0.4458 | 16,662 of 64,309 | 2.73 % |
| TRB | 1.355 | 1.494 | 0.0565 | 0.0303 | 0.4188 | 38,587 of 124,826 | 2.23 % |

Structural floor, same files, `vdjdb.validate.qscore`: TRA $Q=0.1691$ ($h=0.962$, $p=0.093$),
TRB $Q=0.4433$ ($h=0.993$, $p=0.285$).

### 7.1 The decision rule

Two stages, applied in order:

1. **Admissibility.** Four quantities must each be at or above the annotation being replaced:
   $Q$, purity, precision, and the number of epitopes that get at least one cluster. Nothing that
   regresses any of them ships, whatever it gains elsewhere.
2. **Selection.** Among admissible configurations, maximise the independent-study lift.

Never the reverse order, and never either stage alone.

#### The do-nothing partition

The four axes are not four independent guarantees. The do-nothing partition, which puts every
clonotype of an epitope in one cluster and excludes nothing, is an instrument reading rather than an
algorithm; it returns $Q=0.7940$, purity $0.9204$, precision $0.9237$ and coverage $118/118$ on human TRA, and
$Q=0.8632$, purity $0.9343$, precision $0.9455$, coverage $178/178$ on TRB. Its lift is $1.000$ by
construction and its $F_1$ sits at the floor.

On TRA that partition clears all four admissibility axes today, and stage 2 alone excludes it. $Q$
and coverage rise toward it on both chains; purity and precision fall only on TRB, where legacy's
$0.9790$ is above its $0.9343$.

Two rules follow, and they are the operative form of stage 1:

- $Q$ and coverage are guards against the shattering corner, not evidence of quality. Read them as
  floors, never as a ranking, and never across configurations at different retentions.
- An absolute purity floor must exceed the measured purity of that partition, per chain and per
  build, rather than being a round number fixed in advance. On human TRB that purity is $0.9343$, so
  a floor of $0.93$ admits a clustering that has clustered nothing. The floor is $0.94$, the lowest
  round number above the do-nothing purity of both chains, and it is usable only on TRB: the window
  there is $[0.9343, 0.9829]$, while on TRA it is empty by $0.0030$, because no configuration that
  clears $Q$ and coverage reaches a purity above the $0.9204$ that doing nothing scores. TRA is
  therefore guarded by the legacy-relative bar and by stage 2, not by an absolute floor.
  `docs/clustering.md` §8 has the full audit and both windows.

#### Epitope coverage as an admissibility axis

$Q$, purity and precision were the original three axes, and they exclude the *shattering* corner: the
configurations that maximise lift reach $Q=0.023$ at parsimony $0.012$ (§6.2). They do not exclude a
second degenerate corner, which the per-epitope breakdown makes visible.

Stage 2 maximises a **pooled** lift. Lift falls monotonically as a clustering's radius widens (§5),
so stage 2 always prefers the narrowest admissible radius, and a narrow radius finds clusters only
where the data is densest, which means in fewer epitopes. Measured on human TRB, all cells admissible
on the original three axes:

| configuration | pooled lift | epitopes covered, of 178 | epitopes with local lift > 1 |
|---|---:|---:|---:|
| the shipped 2026-06-03 annotation | 2.855 | 103 | 42 |
| TCREMP `coef` 1.15 | **4.490** | 88 | 36 |
| TCREMP `coef` 1.55 | 3.585 | 105 | 42 |
| TCREMP `coef` 1.7 | 3.079 | **109** | **44** |

`coef` 1.15 wins stage 2 by a wide margin and abandons fifteen epitopes the file it replaces covered,
dropping from 42 to 36 epitopes where the clustering beats local chance. Its pooled lift is higher
because it is computed over the clonotypes it kept: it clusters nearly the same *total* (37,094
against 37,210) concentrated into fifteen fewer epitopes and 746 clusters instead of 1,074.

That is a worse database. An epitope with no motif receives no denoising at all (§2), so trading
epitope breadth for a pooled average optimises the metric against the purpose of the stage. Coverage
therefore sits in stage 1, where a regression is disqualifying, rather than in stage 2 where it could
be traded away.

Coverage is counted on epitopes, not clonotypes, because that is the unit the curation acts on. The
per-epitope breakdown is `vdjdb.validate.motif_bench.per_epitope`, and the pooled scorecard is never
reported without it: a pooled retention of 0.33 is a different database depending on whether every
epitope is a third clustered or a third of them are fully clustered and the rest not at all.

## 8. Orthogonal evidence

Structural model confidence is an independent discriminator of true binders from noise: in IMMREP25
only structure-based methods gained meaningfully over random on unseen epitopes, and TCRbridge-style
use of model confidence separates binders from noise. VDJdb reports CDR3–peptide contact counts,
interface ipTM with percentile ranks, and an outlier-binding-mode flag.

These are a second axis of evidence, computed per record and reported alongside the motif flag. They
are not a way to recover coverage lost at the motif stage, and combining them into a single score is
a separate decision with its own validation.

## 9. Checklist for any change to the motif stage

0. Does it beat the do-nothing partition - one cluster per epitope, nothing excluded - on lift and
   $F_1$? That partition clears every other axis on at least one chain (§7.1), so ask this first.
1. Does it change the §11.1 lift and F1, per chain? Report both, with the base rate and n, and the
   publicity-controlled lift beside the raw one (§5.2).
2. Does it change $Q$, with `h` and `p` separately? A gain in one that is a loss in the other is not
   a gain.
3. Does it change the fraction of clustered records? Report it; do not rank on it.
4. Does it fragment a known motif or percolate an epitope? The A\*02:YLQPRTFLL CSAR/GELFF-16 merge
   and the A\*02:NLVPMVATV multi-motif structure are the two standing regression cases.
5. Does it push the chance-recruitment rate $\alpha$ up (§5.1)? At $\alpha\to1$ the enrichment test
   is not discriminating and coverage is meaningless, whatever it reads.
6. Is the result stable to the control draw, the control size, and the seed?
7. TCRvdb is not consulted. It is touched once, at the end, in aggregate.

---

## Notation

| symbol | meaning |
|---|---|
| $n_e$ | records for epitope $e$ in the scored sample |
| $\phi_e$ | fraction of them that bind $e$; unknown by construction |
| $s$, $\theta$ | substitutions allowed; the similarity ball they define |
| $V_s(L)$ | number of peptides within $s$ substitutions of a length-$L$ CDR3 |
| $\hat p_\theta$ | background per-pair probability, $(n_{\mathrm{control}}+1)/(M+1)$ |
| $M$ | control size (1,000,000); also, in §5.2, the event "in a motif" |
| $d_{\min}$ | degree floor (`min_degree`) |
| $\alpha_\theta(n_e)$ | chance-recruitment rate of a non-convergent record |
| $T,R$ | record is a true binder; record is independently replicated |
| $h,p,Q$ | homogeneity, parsimony, and $Q=2hp/(h+p)$ |

Rendered with MathJax on both renderers: Sphinx enables it by default (`sphinx.ext.mathjax`) and
`conf.py` turns on `myst_enable_extensions = ["dollarmath"]` for the `$...$`/`$$...$$` delimiters,
while GitHub renders them natively as-is.
