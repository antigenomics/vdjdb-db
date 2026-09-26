# Clustering VDJdb motifs — the four partitions, their parameters, and their measured scorecards

Which algorithm turns a set of candidate clonotypes into motif clusters, what each one's parameters
do, and what each one measures at. **`docs/denoising.md` decides *what* to optimise; this document
catalogues *what is available to optimise over*.** Read that one first — a scorecard here is
meaningless without §7.1's two-stage rule.

Becomes `docs/standards/clustering.rst` when the Sphinx site lands (ROADMAP phase 13); Markdown until
then so it is useful now. Carries MathJax, so `conf.py` needs
`myst_enable_extensions = ["dollarmath"]`.

Status key: **shipped** — the default in `TUNED` · **wired** — implemented, tested, measured, not
default · **rejected** — measured and ruled out, kept only so the ruling stays checkable.

---

## 0. The shared contract

Every algorithm here produces the same thing: a `labels` array, one entry per clonotype, `-1` for
noise and otherwise an integer **unique across epitopes**. `vdjdb.motifs.tcremp.clusters(...,
labels=...)` takes it from there, so nothing downstream knows which algorithm ran.

That matters for a reason worth stating: the machinery *after* labelling is not neutral.

- **Clusters are split by CDR3 length** before the `min_cluster` filter, because a PWM cannot span
  lengths and `vdjdb-web` splits by `len` anyway (ROADMAP §8.6). So a cluster of 20 spanning four
  lengths becomes four clusters of ~5, and at `min_cluster = 5` some of them disappear. **An algorithm
  that produces length-heterogeneous clusters is penalised by the emit path, not by the metric.**
- **`min_cluster` is applied after the split**, so it is not the same knob as HDBSCAN's
  `min_cluster_size` and the two must not be conflated. Every measurement below fixes
  `min_cluster = 5`, the shipped value, so the algorithm is the only thing varying.
- **Cluster ids follow content, not order** — size descending, then the alphabetically first member —
  so a cluster whose membership is unchanged keeps its id across builds without anything being stored
  (hard rule 9).

All four are seeded from `vdjdb.config.SEED` and verified byte-reproducible across processes
(hard rule 7).

---

## 1. Connected components — **shipped (TCRNET)**

`vdjdb.motifs.cluster._components`. The transitive closure of the similarity graph: two clonotypes
are in one cluster if a path of within-`scope` edges joins them.

**Parameters: none.** That is the point — it has nothing to tune, so it cannot be overfitted, and it
is what the legacy Rmd did.

| Knob | Where | Default | Effect |
|---|---|---|---|
| `scope` | `tcrnet.TUNED` | `1,0,0,1` | edge definition: `subs,ins,dels,total`. The dangerous one — see `denoising.md` §5 |
| `p` | `tcrnet.TUNED` | `0.01` | enrichment threshold that decides which clonotypes enter the graph |
| `min_cluster` | `cluster.MIN_CLUSTER` | `5` | smallest surviving cluster, post length-split |

**The failure mode is percolation.** One spurious edge merges two motifs; at scope 2 a chain of such
edges can swallow an epitope. Measured: at scope `1,0,0,1` the largest component holds ≥ 90 % of a
TRA epitope's clustered clonotypes in **22 of 118** epitopes; median percolation 0.3958.

---

## 2. CPM Leiden — **rejected**

`vdjdb.motifs.cluster._leiden`, reached by setting `tcrnet.TUNED[gene]["resolution"]` to a positive
number. Runs *inside* each connected component, so it can only ever subdivide what §1 found — it
never merges across components and never rescues a clonotype components missed.

### 2.1 Parameters, and why each is what it is

```python
igraph.Graph(n=n, edges=edges).community_leiden(
    objective_function="CPM", resolution=resolution, n_iterations=-1)
```

| Parameter | Value | Why |
|---|---|---|
| `objective_function` | `"CPM"` | **Not modularity.** Modularity's resolution is scaled by the graph's total edge count, so the same number means something different for every epitope — a global setting would be a different filter on a 30-clonotype epitope than on a 30,000 one. CPM's resolution is an absolute internal-density threshold and transfers. |
| `resolution` | `None` (components) | A community is kept while its internal edge density exceeds $\gamma$. $\gamma = 0$ reproduces connected components **exactly** — asserted in `tests/unit/test_motifs.py`, not assumed. |
| `n_iterations` | `-1` | Iterate to convergence rather than a fixed count. A fixed count makes the answer depend on how long it ran. |
| seed | `config.SEED` via `random.seed` | Leiden randomises node visit order and the refinement step, and **python-igraph draws from Python's own global RNG** rather than accepting a seed argument — so the seed must be set on the `random` module immediately before the call, or the partition is not reproducible (hard rule 7). |

CPM maximises

$$\mathcal{H} = \sum_{c} \left( e_c - \gamma \binom{n_c}{2} \right)$$

over partitions, where $e_c$ is the edge count inside community $c$, $n_c$ its size, and $\gamma$ the
resolution. A community survives while its internal density $e_c / \binom{n_c}{2}$ exceeds $\gamma$,
so $\gamma$ reads directly as **the minimum density of a motif**.

### 2.2 The resolution ladder, measured

Human, scope `1,0,0,1`, `p` 0.01, `min_cluster` 5. `res = None` is connected components — the shipped
default and the first row of each table.

**TRA** (legacy bar: `Q` 0.1691, purity 0.8658, precision 0.8567):

| resolution | lift | Q | h | parsimony | purity | precision | retention | cids | clonotypes clustered |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **None** | 2.051 | **0.1798** | 0.9600 | **0.0992** | 0.8754 | 0.8679 | **0.2323** | 859 | **15,608** |
| 0.05 | 2.136 | 0.1218 | 0.9653 | 0.0650 | 0.8943 | 0.8841 | 0.2211 | 1,376 | 14,554 |
| 0.1 | 2.289 | 0.0929 | 0.9696 | 0.0488 | **0.8961** | **0.8853** | 0.2030 | 1,573 | 12,798 |
| 0.2 | 2.949 | 0.0538 | 0.9759 | 0.0276 | 0.8871 | 0.8747 | 0.1581 | 1,380 | 8,591 |
| 0.3 | 4.869 | 0.0259 | 0.9838 | 0.0131 | 0.8824 | 0.8560 | 0.1001 | 526 | 3,803 |
| 0.5 | **5.349** | 0.0233 | **0.9849** | 0.0118 | 0.8852 | 0.8576 | 0.0951 | 452 | 3,352 |
| 0.8 | **5.349** | 0.0233 | **0.9849** | 0.0118 | 0.8852 | 0.8576 | 0.0951 | 452 | 3,352 |

**TRB** (legacy bar: `Q` 0.4433, purity 0.9790, precision 0.9756):

| resolution | lift | Q | h | parsimony | purity | precision | retention | cids | clonotypes clustered |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **None** | 1.501 | **0.4428** | 0.9920 | **0.2850** | **0.9790** | 0.9761 | **0.3337** | 657 | **37,898** |
| 0.05 | **1.518** | 0.2243 | 0.9937 | 0.1264 | 0.9785 | 0.9758 | 0.3181 | 1,949 | 35,822 |
| 0.1 | 1.488 | 0.1429 | 0.9944 | 0.0770 | 0.9782 | 0.9751 | 0.3011 | 3,374 | 33,570 |
| 0.2 | 1.465 | 0.1067 | 0.9958 | 0.0564 | 0.9799 | **0.9774** | 0.2471 | 3,263 | 26,572 |
| 0.3 | 1.293 | 0.0879 | 0.9976 | 0.0460 | 0.9856 | 0.9769 | 0.1868 | 2,103 | 20,011 |
| 0.5 | 1.198 | 0.0854 | **0.9977** | 0.0446 | **0.9856** | 0.9767 | 0.1782 | 1,967 | 19,168 |
| 0.8 | 1.198 | 0.0854 | **0.9977** | 0.0446 | **0.9856** | 0.9767 | 0.1782 | 1,967 | 19,168 |

Lift is on the full cohort here, so these rows are comparable **to each other** but not to
`denoising.md`'s non-display figures.

### 2.3 What the ladder says

**Homogeneity rises while parsimony collapses.** TRA: `h` 0.9600 → 0.9849 as parsimony goes
0.0992 → 0.0118, an 8.4-fold fall. TRB: `h` 0.9920 → 0.9977 against parsimony 0.2850 → 0.0446, 6.4-fold.
That *is* shattering, stated in the two terms: every cluster gets purer because every cluster gets
smaller, and past some point they are small enough to be uninformative. **`Q` is the only one of the
seven columns that notices** — purity and precision barely move (TRA precision actually *falls*,
0.8679 → 0.8576) and lift rises, so any criterion built on purity, precision or lift alone ranks the
degenerate end top.

**On TRA the lift gain is real but bought with coverage.** Lift 2.051 → 5.349 is a 2.6× improvement
in how often a clustered clonotype is independently replicated. It is not an artefact. But the
clustered set shrinks from **15,608 clonotypes to 3,352** — Leiden keeps the densest cores and
discards the rest, so it is not finding better motifs, it is finding *fewer, safer* ones. For a
database whose purpose is denoising the whole corpus, dropping 78 % of what components clustered is
not a trade §7.1 can accept.

**On TRB it does not even buy lift.** The maximum is 1.518 at resolution 0.05, a 1.1 % gain over
components' 1.501, and it falls to 1.198 by 0.5 — while `Q` drops to a fifth. There is no resolution
at which TRB Leiden is worth running.

**The ladder saturates between 0.5 and 0.8** — the two rows are bit-identical on both chains. Once
every surviving community is dense enough to clear $\gamma$, raising $\gamma$ further changes nothing,
so the ladder has a top and `resolution = 0.8` is not "more aggressive" than `0.5`. Anyone extending
the sweep upward is measuring the same partition twice.

**Verdict: no Leiden cell is admissible on either chain.** `resolution` stays in `TUNED`, set to
`None`, and the implementation and its tests stay — the finding is that the knob does not help on this
corpus, and that finding has to remain falsifiable when the corpus changes.

---

## 3. DBSCAN — **shipped (TCREMP)**

`vdjdb.motifs.tcremp.cluster_labels`. Density clustering in the prototype-distance embedding, run
per epitope at one radius estimated on the chain's pooled geometry.

| Parameter | Where | Default | Effect |
|---|---|---|---|
| `eps` | `coef × mean(1st-NN)` | TRA 1.8, TRB 1.15 × | the radius. **The** knob |
| `min_samples` | `tcremp.MIN_SAMPLES` | `2` | the TCREMP paper's Table-1 value, and the smallest number that can be a cluster |
| `n_components` | `tcremp.TUNED` | `50` | PCA dimensions before clustering |
| `min_cluster` | `tcremp.TUNED` | `5` | post length-split filter |

**`eps` is pooled, DBSCAN is per-epitope.** Re-estimating the radius on a per-epitope $n$ of 30–300
is exactly the degenerate regime a knee fails in (§5), so the geometry is estimated once per chain and
the clustering runs per epitope (ROADMAP §8.4).

**The frontier is monotone**, which is what makes tuning a bisection rather than a search: lift falls
as `eps` widens and `Q` rises, so the best admissible cell is always the *smallest* `eps` clearing the
`Q` bar. Both chains therefore sit on the boundary by construction — TRB clears by 0.0004, TRA's next
step down misses by 0.0010.

**`min_cluster` 5 dominates 3** on lift, purity **and** precision at every `coef` measured, so it is
not a trade: clusters of 3–4 contribute low-confidence members, not signal.

---

## 4. HDBSCAN — **rejected**

`vdjdb.motifs.tcremp.hdbscan_labels`, `sklearn.cluster.HDBSCAN` — **no dependency beyond the sklearn
already used for the PCA and DBSCAN.**

### 4.1 Why it was worth measuring

DBSCAN commits to **one radius for every epitope**, and VDJdb's epitopes do not share a density: a
display-selected library yields thousands of receptors one substitution apart by construction (29,688
of 192,753 records) while a 30-record epitope is sparse. HDBSCAN condenses a cluster hierarchy over
mutual-reachability distance and selects by stability, so there is no global radius to be wrong. It
also folds `min_samples` and the `min_cluster` post-filter into one `min_cluster_size`.

### 4.2 Parameters

| Parameter | Swept | Effect |
|---|---|---|
| `min_cluster_size` | 3, 5, 10, 20 | smallest cluster the hierarchy will emit. HDBSCAN's one real knob |
| `cluster_selection_method` | `eom`, `leaf` | excess-of-mass takes the most stable cut (parsimonious); `leaf` takes every condensed-tree leaf (shatters deliberately) |
| `min_samples` | 2, 5 | how conservative the core-distance estimate is; higher declares more points noise |
| `cluster_selection_epsilon` | **not swept** | the anti-shatter knob, and **it is not usable**: see below |
| `allow_single_cluster` | `False` (default) | left alone; `True` would let an epitope legitimately percolate |

**`cluster_selection_epsilon` is broken in sklearn 1.9.1.** On a 12-cell probe at 0.5 and 1.0 × the
chain's mean 1-NN distance, **5 cells raised
`TypeError: only 0-dimensional arrays can be converted to Python scalars`** from
`sklearn/cluster/_hdbscan/_tree.pyx:588` (`traverse_upwards`), and the cells that did run returned
**bit-identical** results to `epsilon = 0.0`. So the one parameter that would directly counter
shattering is unavailable, which is worth recording: it is the reason the `leaf` rows below cannot be
repaired.

### 4.3 The scorecard — 32 cells, human, `min_cluster` fixed at 5

**TRA** (bar: `Q` 0.1691, purity 0.8658, precision 0.8567), ranked by lift:

| `min_cluster_size` | method | `min_samples` | lift | Q | h | parsimony | purity | precision | retention | cids | admissible |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|:-:|
| 3 | leaf | 5 | **1.507** | 0.2023 | 0.948 | 0.113 | 0.8991 | 0.8954 | 0.3000 | 2,384 | — |
| 5 | leaf | 5 | 1.397 | 0.2069 | 0.947 | 0.116 | 0.9012 | 0.8968 | 0.3344 | 2,619 | **yes** |
| 3 | eom | 2 | 1.389 | 0.2120 | 0.946 | 0.119 | **0.9061** | **0.9032** | 0.3244 | 2,567 | **yes** |
| 5 | leaf | 2 | 1.364 | 0.1878 | 0.949 | 0.104 | 0.8982 | 0.8946 | 0.3523 | 3,199 | **yes** |
| 10 | leaf | 5 | 1.265 | 0.2492 | 0.943 | 0.144 | 0.8999 | 0.8949 | 0.4016 | 2,669 | **yes** |
| 5 | eom | 2 | 1.257 | 0.2616 | 0.940 | 0.152 | 0.9017 | 0.8978 | 0.4366 | 3,236 | **yes** |
| 20 | eom | 2 | 1.078 | 0.3979 | 0.922 | 0.254 | 0.8773 | 0.8741 | 0.5276 | 2,639 | **yes** |
| 20 | eom | 5 | 1.009 | **0.5344** | 0.913 | **0.378** | 0.8825 | 0.8788 | **0.5454** | 1,630 | **yes** |

**TRB** (bar: `Q` 0.4433, purity 0.9790, precision 0.9756), ranked by lift:

| `min_cluster_size` | method | `min_samples` | lift | Q | h | parsimony | purity | precision | retention | cids | admissible |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|:-:|
| 3 | leaf | 5 | **1.571** | 0.1707 | 0.980 | 0.093 | 0.9360 | 0.9317 | 0.3905 | 4,379 | — |
| 5 | leaf | 5 | 1.470 | 0.2114 | 0.976 | 0.119 | 0.9357 | 0.9321 | 0.4610 | 4,852 | — |
| 5 | leaf | 2 | 1.404 | 0.2355 | 0.974 | 0.134 | **0.9403** | 0.9372 | 0.5548 | 7,152 | — |
| 5 | eom | 5 | 1.168 | 0.4939 | 0.962 | 0.332 | 0.9396 | 0.9393 | 0.7427 | 4,351 | — |
| 10 | eom | 5 | 1.122 | 0.5154 | 0.962 | 0.352 | 0.9401 | 0.9399 | 0.7758 | 3,960 | — |
| 20 | eom | 2 | 1.116 | 0.5375 | 0.961 | 0.373 | **0.9412** | 0.9392 | **0.8042** | 3,827 | — |
| 20 | eom | 5 | 1.092 | **0.5525** | 0.961 | **0.388** | **0.9412** | **0.9402** | 0.7963 | 3,376 | — |

Full 32 rows in `/tmp/hdbscan_sweep.csv`; the eight per chain above are the ranked extremes and every
admissible TRA cell of interest.

### 4.4 What the scorecard says

**12 of 32 cells are admissible, all of them on TRA, none on TRB.** On TRB purity tops out at 0.9412
against a bar of 0.9790 — every single cell regresses it, by 4 points.

**The `eom` / `leaf` axis behaves exactly as the theory predicts, which is the useful part.** `leaf`
buys lift and destroys parsimony (TRB `leaf` parsimony 0.093–0.134 against `eom`'s 0.332–0.388);
`eom` buys parsimony and retention and gives up lift. The two failure modes bracket the truth and
neither end is admissible on TRB.

**HDBSCAN's real strength is retention and structure, not selectivity.** TRB `eom mcs20 ms2` clusters
**80.4 %** of clonotypes at `Q` 0.5375 — against DBSCAN's 32.3 % at `Q` 0.4437. That is a
*structurally better* partition by the instrument `Q`, and it covers far more of the database. It
fails only because purity falls below legacy's, and it fails on lift badly: 1.116 against DBSCAN's
4.490.

**Verdict: rejected under §7.1 on both chains** — no admissible cell on TRB, and on TRA the best
admissible lift is 1.397 against DBSCAN's 1.855 and TCRNET's 1.977. It stays wired and tested because
the retention result is worth revisiting if the acceptance bar ever moves off legacy purity — see §6.

---

## 5. Choosing `eps` by a knee — **wired, not in the operating path**

`vdjdb.motifs.tcremp.knee`. The published TCREMP protocol picks `eps` from the knee of the sorted
k-distance curve; this repo does not, and the reason is measured.

**The reference implementation is degenerate at production scale.** Kneedle with a degree-10
polynomial fit on pooled human TRB — 112,983 clonotypes — returns **knee index 1, fraction 0.000**. A
degree-10 fit over 113,000 points oscillates, and the reported knee is the first oscillation rather
than a feature of the data. Downstream that puts `eps` below the data entirely and retention
collapses to a few per cent.

`knee` states Kneedle's difference curve directly, with no polynomial:

$$x_i = \frac{i}{n-1}, \qquad y_i = \frac{c_i - \min c}{\max c - \min c}, \qquad
\text{knee} = \arg\max_i (y_i - x_i), \qquad \text{strength} = \max_i (y_i - x_i)$$

for a sorted curve $c$ (signs reversed for a convex curve). Three guards, each against an observed
failure rather than an imagined one:

| Guard | Default | What it catches |
|---|---|---|
| resample onto a fixed grid | 1,000 points | the degree-10 oscillation. Makes the knee a property of the curve's **shape**: measured invariant at 0.799 / 0.800 / 0.800 for one elbow at $n$ = 300 / 30,000 / 300,000 |
| `min_strength` | 0.05 | a curve with no knee. $\max(y-x)$ is **exactly 0** for a straight line, so the statistic that locates a knee is the same one that reports there is none — not a separate test |
| `floor_frac` / `ceil_frac` | 0.05 / 0.95 | a knee pinned to an end. At the bottom `eps` falls below the data; at the top the epitope percolates |

`strength` is scale-invariant by construction, so it is comparable across chains with different
embedding scales (TRA mean 1-NN 9.66, TRB 9.19).

**What runs instead is `eps = coef × mean(1st-NN distance)`**, with `coef` fitted against independent
replication. That is not a knee method and must not be described as one (ROADMAP §8.3).

---

## 6. The per-epitope view — what every pooled number hides

`vdjdb.validate.motif_bench.per_epitope`, counted on **clonotypes** so a clonotype reported by forty
papers does not weigh forty times in its own epitope's retention, and reported over **epitopes** because
that is the unit curation acts on.

A pooled retention of 0.33 is a different database depending on whether every epitope is a third
clustered or a third of them are fully clustered and the rest not at all. On TRB it is much closer to
the second, and no pooled column says so.

### 6.1 The shipped clusterings, per epitope

**TRA — 118 epitopes** with ≥ 30 records:

| method | epitopes covered | none at all | retention median | retention q75 | ≥90 % percolated | percolation median | clonotypes clustered | cids | epitopes with lift > 1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| legacy | 103 | 15 | 0.0647 | 0.1603 | 24 | 0.4717 | 15,063 | 2,214 | 29 |
| TCRNET | 103 | 15 | 0.0706 | 0.1624 | **22** | **0.3958** | 15,480 | 2,173 | 33 |
| TCREMP | **104** | **14** | **0.0850** | **0.2115** | 23 | 0.4834 | **17,380** | **2,926** | **34** |

**TRB — 178 epitopes**:

| method | epitopes covered | none at all | retention median | retention q75 | ≥90 % percolated | percolation median | clonotypes clustered | cids | epitopes with lift > 1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| legacy | 103 | 75 | 0.0089 | 0.0941 | 42 | 0.7692 | 37,210 | 1,074 | 42 |
| TCRNET | **109** | **69** | **0.0142** | 0.0941 | 41 | 0.7059 | 37,109 | 1,051 | **43** |
| TCREMP | 105 | 73 | 0.0128 | **0.1139** | **40** | **0.6667** | **42,802** | **1,156** | 42 |

### 6.2 What it says that the pooled scorecard does not

**TRB is a two-population database and TRA is not.** 75 of 178 TRB epitopes get no motif from the
shipped annotation and the median epitope has **0.9 %** of its clonotypes clustered, while the pooled
retention is 32 %. The pooled figure is carried by a minority of large, deeply-sampled epitopes —
GILGFVFTL, NLVPMVATV, YLQPRTFLL. On TRA the median epitope has 6.5 % clustered against a pooled 21 %,
still skewed but far less bimodal. **Any statement about "VDJdb's motif coverage" that quotes one
number for TRB is describing about sixty epitopes and implying 178.**

**Percolation is the dominant structural defect on TRB, not shattering.** The median TRB epitope has
**77 %** of its clustered clonotypes in a single cluster, and 42 of 178 are ≥ 90 % collapsed. Both
rewrites reduce it (TCRNET to 0.7059, TCREMP to 0.6667) but neither solves it. This is the direction a
future motif change should attack, and it is why §2's Leiden was worth trying at all — subdividing
percolated components is the right idea; the resolutions that do it enough also shatter everything
else.

**Coverage and pooled lift pull in opposite directions**, which is the finding that changed the
shipped parameters. See `denoising.md` §7.1: TCREMP TRB at `coef` 1.15 maximised pooled lift at 4.490
and covered 88 of 178 epitopes, fifteen fewer than the annotation it replaces. At 1.55 it covers 105
and still leads legacy by 25.6 % on lift.

**`lift > 1` per epitope is the honest coverage-weighted read.** Of the epitopes where a local base
rate exists at all, legacy beats chance in 42 of them on TRB and 29 on TRA; the rewrites reach 43/42
and 33/34. That is a modest gain, and stating it beside a pooled "+25.6 % lift" is the difference
between a result and a headline.

### 6.3 The report

`out/reports/motifs_per_epitope.tsv`, 888 rows — one per (gene, method, epitope) — written by the
build and uploaded as a CI artifact, never shipped in a bundle (`outputs.md` §5). Columns:
`clonotypes`, `clustered`, `retention`, `clusters`, `largest_cluster`, `mean_cluster_size`,
`singleton_clusters`, `percolation`, `replicated`, `tp`, `precision`, `lift`.

`lift` is **null, not zero**, where an epitope has no replicated clonotype: there is no base rate to
lift against, and averaging a fabricated zero is how a per-epitope mean ends up below the pooled figure
for no reason. 47 of 118 TRA epitopes and 60 of 178 TRB epitopes have a defined lift.

---

## 7. Summary, and what would change a verdict

| Algorithm | Status | TRA lift | TRB lift | TRB epitopes covered | Distinguishing property |
|---|---|---:|---:|---:|---|
| the shipped 2026-06-03 annotation | baseline | 1.703 | 2.855 | 103 / 178 | what every bar is set from |
| connected components | **shipped**, TCRNET | **1.977** | 3.194 | **109** | no parameters to overfit |
| CPM Leiden | rejected | — | — | — | lift 2.6× on TRA at parsimony ÷ 8.4 |
| DBSCAN | **shipped**, TCREMP | 1.855 | **3.586** | 105 | monotone frontier; tuning is a bisection |
| HDBSCAN | rejected | 1.397 | — | 119–146, none admissible | retention 0.80 at `Q` 0.5375; purity −4 pts |

Leiden has no admissible cell on either chain; HDBSCAN has twelve on TRA and none on TRB. The two
rejected rows carry a dash where no admissible configuration exists, rather than their best
inadmissible number — an inadmissible lift is not a lift this project can spend.

Both rejections are **conditional on the acceptance bar being legacy's purity**, and that is the
assumption most worth attacking. Legacy purity is not a law of nature — it is the 2026-06-03 file's
number, and `denoising.md` §7 adopts it so that no release ever regresses. Three things would reopen
these verdicts, in descending order of how likely they are to matter:

1. **A bar that is absolute rather than relative.** If purity ≥ 0.93 were acceptable on TRB, HDBSCAN
   `eom mcs20 ms2` ships instead — 80 % retention at `Q` 0.5375 — and the database gains motif
   coverage over most of its epitopes rather than a third of them. That is a curation decision, not a
   measurement.
2. ~~Per-epitope coverage entering the criterion.~~ **Settled**: it is now a fourth admissibility
   axis (`denoising.md` §7.1). Stage 2 maximises pooled lift, lift prefers the narrowest admissible
   radius, and a narrow radius covers fewer epitopes -- measured, TCREMP TRB at `coef` 1.15 covered
   88 of 178 against legacy's 103. The shipped TRB radius moved to 1.55 as a result.
3. **`cluster_selection_epsilon` becoming usable** in a later sklearn, which is the only knob that
   would let HDBSCAN's `leaf` lift be had without its shattering.

A rejected algorithm stays implemented and tested. The measurements above are the reason it is not
default, and a measurement that cannot be re-run is not a reason.
