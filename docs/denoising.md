# Denoising VDJdb — what motifs are for, and how to tune them

**Status: normative.** This is the specification the motif stage is built and judged against. When
this document and a measurement disagree, the measurement is reported and this document is revised
deliberately — not silently worked around.

## Sources

Two papers define the method. Neither is restated from memory; every claim below is traceable to one
of them.

1. **Shugay M, Luppov DV, Vlasova EK, Chudakov DM.** *Towards high-quality large-scale T cell
   receptor antigen specificity data: challenges and promises.* **Nature Methods**, 3 September 2026.
   [10.1038/s41592-026-03227-2](https://doi.org/10.1038/s41592-026-03227-2) · PMID 42693187.
   *(Record verified against PubMed.)*
2. **Luppov DV, Koneva AE, Bagaev DV, Alexandrova AV, Vlasova EK, Chudakov DM, Motozono C,
   Sewell AK, Shugay M.** *VDJdb in 2026: boosting T-cell receptor recognition evidence using
   paratope embeddings and AI-based structure prediction.* **Nucleic Acids Research**, 2026,
   Database Issue. *In press — DOI assigned at production.*

Two further documents are the local instruments rather than the doctrine, and are cited where used:
the IMMREP25 audit (`~/vcs/manuscripts/2026-immrep25-audit`) for the calibrated cohort ladder, and
the validation manuscript (`~/vcs/manuscripts/2026-vdjdb-db-validation`) for the held-out check.

---

## 1. The premise: a substantial fraction of assay-derived records are false positives

This is the starting assumption, not a hedge. Multimer and screen-derived records carry bystanders
and mis-assigned specificities, so **a record's presence in the database is not evidence that the
receptor binds the epitope it is filed under**.

Everything downstream follows from that. VDJdb's job is not only to accumulate records but to carry,
per record, how much of the database's own structure supports it.

## 2. The rule: motif membership is the noise filter

> Sequences not associated with any statistically significant motif are flagged as likely false
> positives.

That is the operative definition, established for TCRNET and carried forward to the embedding-based
annotation. A motif is convergent selection made visible: several independent receptors arriving at
the same paratope solution for the same pMHC. A record that no other record resembles has nothing
supporting it beyond the single assay that produced it.

**Three consequences the code must honour.**

1. **`in.motif` is a per-record evidence column, not a display feature.** The validation manuscript
   scores it on the same axis as the confidence score and the independent-study count. It has to be
   emitted for every record, including the records that get `false`.
2. **A cluster's job is to be *right*, not to be *large*.** Coverage is a reported quantity, never a
   maximised one — see §5.
3. **Absence of a motif is a weak signal, not a verdict.** Low-generation-probability receptors
   (autoimmunity-associated variants are the named case) and under-sampled epitopes will be
   motif-less for reasons that have nothing to do with being wrong. The flag is advisory and is
   never a deletion.

## 3. The noise signature: many epitopes, no motif

The concrete, computable read on which records are spurious:

| CDR3 is associated with… | …and sits in | reading |
|---|---|---|
| one or a few epitopes | a motif | specific, supported |
| several epitopes | a few motif contexts | **genuine cross-reactivity** |
| an unexpectedly high number of epitopes | **no motif at all** | **false positive** |

The third row is the actionable one: promiscuity *plus* motif-lessness is the noise signature, and
either alone is not. A CDR3 in many epitopes that still sits in coherent motifs is cross-reactive,
which is real biology; a CDR3 in many epitopes with no motif anywhere is a sequence that keeps
turning up in assays.

**So promiscuity must never be used as a filter on its own**, and the epitope count and the motif
assignment have to be reported together.

## 4. Sensitivity: why TCRNET alone is not enough

TCRNET is **highly stringent** — a Hamming-1 CDR3 graph, with limited handling of CDR1 and CDR2 —
and that stringency restricts sensitivity. Two named failure modes:

- **Fragmentation.** One motif is split across clusters that differ only in CDR3 length. The
  documented case is HLA-A\*02:YLQPRTFLL TRB, where TCRNET clusters C21 and C24 carry highly similar
  length-16 logos and stay separate; the embedding-based annotation merges them into one cluster
  (C84) spanning several raw CDR3 lengths, on the shared CSAR/GELFF-16 motif.
- **CDR1/CDR2 blindness.** V-gene identity carries two thirds of the paratope. A CDR3-only distance
  cannot see it, which is why V-gene bias (TRBV20-1 for HLA-A\*02:GLCTLVAML, TRAV12-2 for
  HLA-A\*02:LLWNGPMAV and HLA-A\*02:ELAGIGILTV) is invisible to the graph that defines TCRNET's balls.

The embedding aggregates V segments with variable CDR3 lengths, which is what buys larger, less
fragmented clusters and higher coverage on both chains.

**An epitope may carry several disjoint motifs.** HLA-A\*02:NLVPMVATV is the named case. Any
procedure that assumes one motif per epitope, or that scores a clustering by how close it comes to
one cluster per epitope, is measuring the wrong thing.

## 5. The trap: coverage bought by widening the neighbourhood

**Fragmentation and percolation are both failures, and fixing one by causing the other is not
progress.**

Widening the similarity ball raises the fraction of clustered records monotonically. It does so by
recruiting exactly the bystanders and mis-assigned records §1 warns about — and once they sit inside
a named motif with a logo, the clustering has **laundered them into apparent epitope-specific
signal**. Every aggregate coverage metric improves while the database gets worse.

This is not hypothetical; it is measured on this corpus. TRA connected components at
scope `1,0,0,1` → `2,0,0,2` → `3,0,0,3` → `4,0,0,4` raise clustered-record fraction
0.2105 → 0.4149 → 0.4692 → 0.4807, while the enrichment of clustered clonotypes for independent
replication falls from **1.96×** to **1.30×** at the second step alone. More coverage, less signal.

> **Rule.** Never rank motif configurations on coverage, retention, or the clustered fraction.
> Report them; do not optimise them.

### 5.1 The mechanism, measured — and two plausible explanations the data rejects

The combinatorial argument is the obvious starting point and it is **not** what happens. It is worth
writing down both, because the gap between them is itself the finding.

**The bound.** The set of peptides within $s$ substitutions of a CDR3 of length $L$ has size

$$V_s(L)=\sum_{j=0}^{s}\binom{L}{j}\,19^{j},$$

so for $L=14$: $V_1=267$, $V_2=33{,}118$, $V_3=2{,}529{,}794$. Scope 2 is a **124-fold** larger ball
than scope 1. If sequence space were uniformly occupied, the chance neighbour rate would inflate by
that factor.

**What it actually inflates by: 12×.** Measured on human TRB, 124,819 scored clonotypes, background
$M=10^{6}$: the median $\hat p_\theta=(n_{\mathrm{control}}+1)/(M+1)$ is
$1.00\times10^{-6}$ at scope 1 and $1.20\times10^{-5}$ at scope 2 (means $3.48\times10^{-6}$ and
$4.28\times10^{-5}$). That is a **12-fold** inflation against a 124-fold bound. Real repertoire space
is concentrated, not uniform, so ball volume overstates the null by an order of magnitude. **Never
substitute $V_s$ for a measured $\hat p_\theta$.**

**Rejected explanation 1 — chance recruitment through the degree test.** With
$D\sim\mathrm{Bin}(n_e-1,\hat p_\theta)$ and recruitment at $D\ge d_{\min}$, the chance-recruitment
rate is $\alpha_\theta(n_e)=\Pr[D\ge d_{\min}]$. Computed on all 87 human TRB epitopes with $\ge100$
clonotypes: $\alpha\le0.0011$ at scope 1 and $\alpha\le0.094$ at scope 2 (mean $0.0027$), the maximum
falling on VEALYLVCG. **The degree test does not become vacuous at scope 2**, and chance recruitment
through it is not the mechanism. ($\alpha$ still depends on $n_e$, so a fixed $p$ is still not a
fixed noise level — just not by enough to matter here.)

**Rejected explanation 2 — untested neighbour recruitment.** Stage I admits a clonotype for being
within scope of an enriched one, without its passing any test, so a denser ball should drag in more
free riders. It does the opposite: untested members are **15.9 %** of the clustered set at scope 1
and **12.7 %** at scope 2.

**What it is: the test's own operating point moves.** The same $p<0.05$ passes **35,775** clonotypes
at scope 1 and **54,339** at scope 2, a **52 % increase**, because within-sample degree grows faster
with the ball than the background expectation does. The extra passes are real passes of a real test —
and they are enriched for bystanders, which is why the independent-replication lift falls from
**1.96×** to **1.30×** across the same step.

Tightening the threshold does not recover it: at scope 2, $p<0.01$ gives lift 1.35 against $p<0.05$'s
1.30. The scope is not a knob that a threshold can compensate for, because the $p$-value is already
computed against a scope-matched background — widening the ball changes what the test is *about*,
not merely how strict it is.

`vdjdb.validate.noise.ball_volume` and `chance_recruitment` compute the first two; $\phi_e$, the
share of epitope $e$'s records that really bind it, is unknown by construction, so
$\mathbb{E}[F_e]=(1-\phi_e)\,n_e\,\alpha_\theta(n_e)$ is used as a sensitivity curve over plausible
$\phi$ and never as a point estimate.

### 5.2 What the lift actually estimates

Write $T$ for "record is a true binder" with prevalence $\phi=\Pr(T)$, $M$ for "in a motif", $R$ for
"independently replicated", and $r_1=\Pr(R\mid T)$, $r_0=\Pr(R\mid\neg T)$. If $M$ and $R$ are
conditionally independent given $T$ — both being consequences of a receptor's being a real,
convergent solution — then

$$\Pr(R\mid M)=r_1\Pr(T\mid M)+r_0\bigl(1-\Pr(T\mid M)\bigr),\qquad
\Pr(R)=r_1\phi+r_0(1-\phi),$$

and when $r_0\ll r_1$, which is the statement that a mis-assigned record is not independently
re-reported except by coincidence, the lift collapses to

$$\mathrm{lift}=\frac{\Pr(R\mid M)}{\Pr(R)}\;\approx\;\frac{\Pr(T\mid M)}{\phi}.$$

So **the independent-study lift is an estimate of how much motif membership enriches for true
binders**, which is exactly the quantity the denoising claim needs, and it inverts to
$\Pr(T\mid M)\approx \mathrm{lift}\times\phi$.

⚠ **The conditional-independence assumption can fail through publicity**, and it has to be checked
rather than assumed. A high-generation-probability clonotype is more likely *both* to be re-observed
by a second laboratory *and* to have sequence neighbours, so part of a raw lift could be
$P_{\mathrm{gen}}$ rather than specificity — the IMMREP25 audit finds exactly this confound dominant
in a weak cohort.

The control is the audit's own: permute $R$ **within strata of $\log_{10}P_{\mathrm{gen}}$**, which
preserves each clonotype's publicity and destroys only its association with the clustering, then
report

$$\mathrm{lift}_{\mathrm{ctrl}}=\frac{\mathrm{lift}}{\mathbb{E}[\mathrm{lift}\mid\text{permuted}]}.$$

A ratio near 1 says the raw lift *was* the covariate.

**Measured on VDJdb, it is not.** On the shipped 2026-06-03 annotation, the within-$P_{\mathrm{gen}}$
null sits at **1.043** (TRA) and **0.907** (TRB) against raw lifts of 1.721 and 1.355, giving
controlled ratios of **1.650** and **1.494** ($p<0.005$, 200 permutations). Publicity explains
essentially none of the enrichment here, so on this corpus the raw lift can be read as specificity
enrichment. That is a measurement on one corpus and not a general licence: recompute it whenever the
clustering or the corpus changes.

Because the confound is small here, $\Pr(T\mid M)\approx\mathrm{lift}\times\phi$ is usable — but
$\phi\le1/\mathrm{lift}$ still should not be quoted as a bound on how much of VDJdb is wrong, since
it inherits every remaining assumption ($r_0\ll r_1$ in particular) and the scored cohort is not the
database.

`vdjdb.validate.noise.controlled_lift`, with `pgen_stratum` supplying the strata from the
`cdr3nt.pgen` column the build already generates. **Report the raw lift and the controlled lift
together, never the raw one alone.**

## 6. What to optimise instead — three instruments, three jobs

| | question it answers | when it is used |
|---|---|---|
| **§11.1 independent-study lift** | does motif membership predict a record an independent laboratory also reported? | **picks the shipped configuration** |
| **Q = 2hp/(h+p)** | is the motif structure sound — neither fragmented nor percolated? | compares algorithms and releases |
| **MATCHMAKERS / TCRvdb** | does the denoising recover functionally validated records? | **held out, touched once, aggregate only** |

### 6.1 The tuning objective — independent replication

The signal is how many distinct publications report the same clonotype against the same epitope. A
clustering that is finding real convergent selection recovers exactly those clonotypes.

Scored as **F1 of `clustered` predicting `replicated`, with the lift over the base rate reported
beside it** (`vdjdb.motifs.tcremp._objective`). **F1 and not recall**, for the §5 reason: clustering
everything wins recall and loses precision, and precision here *is* "of the clonotypes we put in
motifs, how many did a second laboratory independently report". Bystanders are, by construction, not
independently replicated.

Lift is quoted with the base rate, always — the base rate is ~2.2 % of clonotype-epitope pairs, so a
precision of 0.15 is a 6.8× enrichment and means nothing stated alone.

#### The denominator is the non-display cohort, and it has to be stated

A **display-selected** record — `method.identification` containing `display` — is not an independent
natural observation: a library panned against one pMHC yields thousands of receptors one
substitution apart *by construction*, and the whole panning experiment is a single `reference.id`.
So every display clonotype contributes **zero** independently-replicated pairs while occupying a
denominator slot.

Measured on human TRB, 2026-09-26:

| Cohort | Clonotype-epitope pairs | Replicated | Base rate |
|---|---:|---:|---:|
| all | 116,053 | 2,771 | 2.3877 % |
| excluding display | 86,361 | 2,771 | 3.2086 % |

The replicated count is **identical** in both rows — which is the measurement, not an assumption:
all 2,771 replicated pairs sit outside the display set, and display contributes 29,692 pairs that
can only ever be false positives for this objective. Scoring lift on the full cohort therefore
deflates precision mechanically, by an amount that depends on how aggressively a given parameter
clusters the display block rather than on whether it is finding convergent selection.

**So lift is scored on the non-display subset**, which is what `vdjdb.motifs.tcremp.fit_coef` has
always documented. Two consequences worth stating plainly, because getting this wrong is silent:

- **Clustering still runs on everything.** The hold-out is the *scoring* denominator, not an
  exclusion from the data. Whether display records should also be held out of the clustering is a
  separate, open question (ROADMAP §30.6).
- **A lift figure is not comparable across denominators.** The same clustering reads 5.31× on the
  non-display cohort and a much lower number on the full one. Any lift quoted anywhere in this
  repository states which cohort it is on, or it is not a number.

Both conventions are reported side by side in the tuning sweeps for exactly this reason.

### 6.2 The structural instrument — Q

`Q = 2hp/(h+p)`, the homogeneity–parsimony trade-off (Tiffeau-Mayer, arXiv:2607.20799), from the
authors' `clustereval`. Homogeneity `h` asks whether clusters predict the epitope; parsimony `p` asks
whether they do it without shattering each epitope into singletons. Both are normalised by the
maximum attainable for the label partition, so epitope count and size divide out.

It is the right instrument here because **it scores §4's failure and §5's failure with one number**:

- shatter an epitope into singletons → `h = 1`, `p = 0`, **Q = 0**
- percolate it into a single cluster → `p = 1`, `h = 0`, **Q = 0**

That is why it replaces the purity/retention pair, which is two numbers that trade against each
other and can be moved in opposite directions to flatter effect. `vdjdb.validate.qscore`.

⚠ **Q is keyed on the clonotype, never on `(epitope, clonotype)`.** Production runs motif detection
inside one epitope's record set, so no other epitope's records are ever in the input and `h` is 1.000
by construction there. Keyed on the clonotype it measures something real: whether a CDR3 filed under
two epitopes ends up in one motif.

### 6.3 The held-out check

MATCHMAKERS / TCRvdb, read only via `$VDJDB_TCRVDB`, enforced by `vdjdb.validate.guard`. **Never
tuned on. Never shipped, committed or redistributed. Aggregate metrics only — never per-record
verdicts.** Two epitopes is too narrow a basis for a global hyperparameter, and tuning on it would
destroy the only independent read this project has.

## 7. The acceptance bar, and the baseline it is measured against

From the published validation of the current annotation: **precision maintained, recall and F1
improved**. A configuration that trades precision for recall has not met the bar; it has found a
different operating point and must be reported as one.

**The baseline, measured on the shipped 2026-06-03 `cluster_members.txt`.** These are the numbers a
rebuild has to beat, and they were never computed during the purity-era tuning because that era
never asked this question:

| chain | lift | controlled | F1 | precision | recall | clustered | base rate |
|---|---:|---:|---:|---:|---:|---:|---:|
| TRA | 1.721 | 1.650 | 0.0849 | 0.0469 | 0.4458 | 16,662 of 64,309 | 2.73 % |
| TRB | 1.355 | 1.494 | 0.0565 | 0.0303 | 0.4188 | 38,587 of 124,826 | 2.23 % |

Structural floor, same files, `vdjdb.validate.qscore`: **TRA $Q=0.1691$** ($h=0.962$, $p=0.093$),
**TRB $Q=0.4433$** ($h=0.993$, $p=0.285$).

### 7.1 The decision rule

Two stages, applied in order — neither alone is sufficient, and the sweeps show why:

1. **Admissibility.** Four quantities must each be at or above the annotation being replaced:
   $Q$, purity, precision, and **the number of epitopes that get at least one cluster**. Nothing
   that regresses any of them ships, whatever it gains elsewhere.
2. **Selection.** Among admissible configurations, **maximise the independent-study lift**.

Never the reverse order, and never either stage alone.

#### What stage 1 can and cannot do — measured

The four axes are not four independent guarantees. Score the partition that puts **every clonotype of
an epitope in one cluster and excludes nothing** — an instrument reading, not an algorithm — and it
returns $Q=0.7940$, purity $0.9204$, precision $0.9237$ and coverage $118/118$ on human TRA, and
$Q=0.8632$, purity $0.9343$, precision $0.9455$, coverage $178/178$ on TRB. Its lift is $1.000$ by
construction and its $F_1$ sits at the floor.

So on TRA **that partition clears all four admissibility axes today**, and what excludes it is stage
2 alone. $Q$ and coverage rise toward it on both chains; only purity and precision fall, and only on
TRB, where legacy's $0.9790$ is above its $0.9343$.

Two rules follow, and they are the operative form of stage 1:

- **$Q$ and coverage are guards against the shattering corner, not evidence of quality.** Read them
  as floors, never as a ranking, and never across configurations at different retentions.
- **An absolute purity floor must exceed the measured purity of that partition**, per chain, per
  build — not a round number fixed in advance. On human TRB that is $0.9343$, so a floor of $0.93$
  admits a clustering that has clustered nothing. **The floor is $0.94$**, the lowest round number
  above the do-nothing purity of both chains, and it is usable only on TRB: the window there is
  $[0.9343, 0.9829]$, while on TRA it is **empty by $0.0030$** — no configuration that clears $Q$ and
  coverage reaches a purity above the $0.9204$ that doing nothing scores. TRA is therefore guarded by
  the legacy-relative bar and by stage 2, not by an absolute floor. `docs/clustering.md` §8 carries
  the full audit and both windows.

#### Why epitope coverage is an admissibility axis and not a tiebreak

$Q$, purity and precision were the original three, and they exclude the *shattering* corner: the
configurations that maximise lift reach $Q=0.023$ at parsimony $0.012$ (§6.2). They do **not** exclude
a second degenerate corner, and it took the per-epitope breakdown to see it.

Stage 2 maximises a **pooled** lift. Lift falls monotonically as a clustering's radius widens (§5), so
stage 2 always prefers the narrowest admissible radius — and a narrow radius finds clusters only where
the data is densest, which means **in fewer epitopes**. Measured on human TRB, all cells admissible on
the original three axes:

| configuration | pooled lift | epitopes covered, of 178 | epitopes with local lift > 1 |
|---|---:|---:|---:|
| the shipped 2026-06-03 annotation | 2.855 | 103 | 42 |
| TCREMP `coef` 1.15 | **4.490** | **88** | **36** |
| TCREMP `coef` 1.55 | 3.585 | 105 | 42 |
| TCREMP `coef` 1.7 | 3.079 | 109 | 44 |

`coef` 1.15 wins stage 2 by a wide margin and **abandons fifteen epitopes the file it replaces
covered**, dropping from 42 to 36 epitopes where the clustering beats local chance. Its pooled lift is
higher because it is computed over the clonotypes it kept; it clusters essentially the same *total*
(37,094 against 37,210) concentrated into fifteen fewer epitopes and 746 clusters instead of 1,074.

That is a worse database. **An epitope with no motif receives no denoising at all** — the whole point
of the stage (§2) — so trading epitope breadth for a pooled average optimises the metric against the
purpose. Coverage therefore sits in stage 1, where a regression is simply disqualifying, rather than
in stage 2 where it could be traded away.

Counted on **epitopes, not clonotypes**, because that is the unit the curation acts on. The
per-epitope breakdown is `vdjdb.validate.motif_bench.per_epitope`, and the pooled scorecard must never
be reported without it: a pooled retention of 0.33 is a different database depending on whether every
epitope is a third clustered or a third of them are fully clustered and the rest not at all.

## 8. Orthogonal evidence, not a substitute

Structural model confidence is an **independent** discriminator of true binders from noise — the
IMMREP25 result was that only structure-based methods gained meaningfully over random on unseen
epitopes, and TCRbridge-style use of model confidence separates binders from noise. VDJdb carries
CDR3–peptide contact counts, interface ipTM with percentile ranks, and an outlier-binding-mode flag.

These are a **second axis of evidence**, computed per record and reported alongside the motif flag.
They are not a way to recover coverage lost at the motif stage, and combining them into a single
score is a separate decision with its own validation.

## 9. Checklist for any change to the motif stage

0. Does it beat the **do-nothing partition** — one cluster per epitope, nothing excluded — on lift
   and $F_1$? Every other axis is cleared by that partition on at least one chain (§7.1), so this is
   the first question and not the last.
1. Does it change the **§11.1 lift and F1**, per chain? Report both, with the base rate and n —
   and the **publicity-controlled** lift beside the raw one (§5.2).
2. Does it change **Q**, with `h` and `p` separately? A gain in one that is a loss in the other is
   not a gain.
3. Does it change the **fraction of clustered records**? Report it. It is not an argument.
4. Does it **fragment** a known motif or **percolate** an epitope? The A\*02:YLQPRTFLL CSAR/GELFF-16
   merge and the A\*02:NLVPMVATV multi-motif structure are the two standing regression cases.
5. Does it push the **chance-recruitment rate** $\alpha$ up (§5.1)? At $\alpha\to1$ the enrichment
   test is not discriminating and coverage is meaningless, whatever it reads.
6. Is the result **stable** to the control draw, the control size, and the seed?
7. TCRvdb is **not** consulted. It is touched once, at the end, in aggregate.

---

## Notation

| symbol | meaning |
|---|---|
| $n_e$ | records for epitope $e$ in the scored sample |
| $\phi_e$ | fraction of them that really bind $e$ — unknown by construction |
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
