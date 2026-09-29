# Motif clustering algorithms, parameters and scorecards

Which algorithm turns a set of candidate clonotypes into motif clusters, what each one's parameters
do, and what each one measures at. `docs/denoising.md` defines what to optimise and this document
catalogues what is available to optimise over, so read `docs/denoising.md` first: a scorecard here is
read against its §7.1 two-stage rule.

Rendered on the documentation site from this file by `myst-parser`. It stays Markdown rather than becoming `.rst`, because it is cited by section number from `ROADMAP.md`, from the package docstrings and from `docs/tuning/`, and converting it would break those references.

Status key: **shipped** - the default in `TUNED` · **wired** - implemented, tested, measured, not
default · **rejected** - measured and ruled out, kept so the ruling stays checkable ·
**measured** - scored on the same instruments from a script under `docs/tuning/`, not in the package.

Every number below comes from `docs/tuning/scorecard.tsv`, 252 configurations scored through one
harness on one cohort, and the markdown is generated from it by `docs/tuning/report.py`, so no figure
in this document is transcribed by hand. §11 has the commands.

**The cohort is `chunks/` as of 2026-09-26, and the legacy bar below is that cohort's.** The released
`cluster_members.txt` never changes, but the corpus it is scored against does, so the bar moves when a
chunk merges: three merges on 2026-09-28 took TRA's `Q` from 0.1691 to 0.1668, its purity from 0.8658
to 0.8641 and its precision from 0.8567 to 0.8548, and TRB's `Q` from 0.4433 to 0.4391. **The numbers
here are not restated for that**, because every configuration in the scorecard was scored against the
same bar in the same run and re-stating one column would break the comparison the document exists to
make. The current values of every axis for the released files, the latest release and the build are
measured on every build and recorded in `rules/motif_metrics.tsv`, which is what a regression is gated
against (`vdjdb motif-metrics`).

---

## 0. Shared contract

Every algorithm here produces a `labels` array, one entry per clonotype, `-1` for noise and
otherwise an integer unique across epitopes. `vdjdb.motifs.tcremp.clusters(..., labels=...)`
consumes it, so nothing downstream knows which algorithm ran.

The machinery after labelling is not neutral:

- Clusters are split by CDR3 length before the `min_cluster` filter, because a PWM cannot span
  lengths and `vdjdb-web` splits by `len` anyway (ROADMAP §8.6), so a cluster of 20 spanning four
  lengths becomes four clusters of ~5 and at `min_cluster = 5` some of them disappear. An algorithm
  that produces length-heterogeneous clusters is penalised by the emit path, not by the metric.
- `min_cluster` is applied after the split, so it is not the same knob as HDBSCAN's
  `min_cluster_size` and the two must not be conflated. Every measurement below fixes
  `min_cluster = 5`, the shipped value, so the algorithm is the only thing varying.
- Cluster ids follow content, not order, size descending and then the alphabetically first member,
  so a cluster whose membership is unchanged keeps its id across builds without anything being
  stored (hard rule 9).

All four are seeded from `vdjdb.config.SEED` and verified byte-reproducible across processes
(hard rule 7).

---

## 1. Connected components - shipped (TCRNET)

`vdjdb.motifs.cluster._components`. The transitive closure of the similarity graph: two clonotypes
are in one cluster if a path of within-`scope` edges joins them. Parameters: none, so there is
nothing to tune and nothing to overfit; it is what the legacy Rmd did.

| Knob | Where | Default | Effect |
|---|---|---|---|
| `scope` | `tcrnet.TUNED` | `1,0,0,1` | edge definition: `subs,ins,dels,total` - see `denoising.md` §5 |
| `p` | `tcrnet.TUNED` | `0.01` | enrichment threshold that decides which clonotypes enter the graph |
| `min_cluster` | `cluster.MIN_CLUSTER` | `5` | smallest surviving cluster, post length-split |

The failure mode is percolation: one spurious edge merges two motifs, and at scope 2 a chain of such
edges can span an epitope. At scope `1,0,0,1` the largest component holds ≥ 90 % of a TRA epitope's
clustered clonotypes in 22 of 118 epitopes; median percolation 0.3958.

---

## 2. CPM Leiden - rejected

`vdjdb.motifs.cluster._leiden`, reached by setting `tcrnet.TUNED[gene]["resolution"]` to a positive
number. It runs inside each connected component, so it subdivides what §1 found: it never merges
across components and never recovers a clonotype components missed.

### 2.1 Parameters

```python
igraph.Graph(n=n, edges=edges).community_leiden(
    objective_function="CPM", resolution=resolution, n_iterations=-1)
```

| Parameter | Value | Notes |
|---|---|---|
| `objective_function` | `"CPM"` | Not modularity: modularity's resolution is scaled by the graph's total edge count, so one number is a different filter on a 30-clonotype epitope than on a 30,000 one. CPM's resolution is an absolute internal-density threshold and transfers across epitopes. |
| `resolution` | `None` (components) | A community is kept while its internal edge density exceeds $\gamma$. $\gamma = 0$ reproduces connected components exactly, asserted in `tests/unit/test_motifs.py`. |
| `n_iterations` | `-1` | Iterate to convergence rather than a fixed count, which would make the answer depend on how long it ran. |
| seed | `config.SEED` via `random.seed` | Leiden randomises node visit order and the refinement step, and python-igraph draws from Python's own global RNG rather than accepting a seed argument, so the seed must be set on the `random` module immediately before the call, or the partition is not reproducible (hard rule 7). |

CPM maximises

$$\mathcal{H} = \sum_{c} \left( e_c - \gamma \binom{n_c}{2} \right)$$

over partitions, where $e_c$ is the edge count inside community $c$, $n_c$ its size, and $\gamma$ the
resolution. A community survives while its internal density $e_c / \binom{n_c}{2}$ exceeds $\gamma$,
so $\gamma$ reads directly as the minimum density of a motif.

### 2.2 The resolution ladder

Human, scope `1,0,0,1`, `p` 0.01, `min_cluster` 5. `res = None` is connected components, the shipped
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

Lift is on the full cohort here, so these rows are comparable to each other but not to
`denoising.md`'s non-display figures.

### 2.3 Ladder results

Homogeneity rises while parsimony falls: on TRA `h` goes 0.9600 → 0.9849 as parsimony goes
0.0992 → 0.0118, an 8.4-fold fall, and on TRB `h` goes 0.9920 → 0.9977 against parsimony
0.2850 → 0.0446, a 6.4-fold fall. That is shattering: every cluster gets purer because every
cluster gets smaller, and past some point they are small enough to be uninformative. `Q` is the only
one of the seven columns that records it: purity and precision barely move (TRA precision falls,
0.8679 → 0.8576) and lift rises, so a criterion built on purity, precision or lift alone ranks the
degenerate end top.

On TRA lift rises at the cost of coverage. Lift 2.051 → 5.349 is a 2.6× improvement in how often a
clustered clonotype is independently replicated, but the clustered set shrinks from 15,608
clonotypes to 3,352, because Leiden keeps the densest cores and discards the rest. Dropping 78 % of
what components clustered is not a trade §7.1 can accept.

On TRB there is no lift gain: the maximum is 1.518 at resolution 0.05, a 1.1 % gain over components'
1.501, and it falls to 1.198 by 0.5 while `Q` drops to a fifth of its value. No resolution improves
TRB over connected components.

The ladder saturates between 0.5 and 0.8, where the two rows are bit-identical on both chains. Once
every surviving community is dense enough to clear $\gamma$, raising $\gamma$ further changes
nothing, so the ladder has a top, `resolution = 0.8` is not more aggressive than `0.5`, and
extending the sweep upward measures the same partition twice.

No Leiden cell is admissible on either chain. `resolution` stays in `TUNED`, set to `None`, and the
implementation and its tests stay, so the measurement can be re-run when the corpus changes.

---

## 3. DBSCAN - shipped (TCREMP)

`vdjdb.motifs.tcremp.cluster_labels`. Density clustering in the prototype-distance embedding, run
per epitope at one radius estimated on the chain's pooled geometry.

| Parameter | Where | Default | Effect |
|---|---|---|---|
| `eps` | `coef × mean(1st-NN)` | TRA 1.8, TRB 1.15 × | the radius, and the primary knob |
| `min_samples` | `tcremp.MIN_SAMPLES` | `2` | the TCREMP paper's Table-1 value, and the smallest number that can be a cluster |
| `n_components` | `tcremp.TUNED` | `50` | PCA dimensions before clustering |
| `min_cluster` | `tcremp.TUNED` | `5` | post length-split filter |

`eps` is pooled and DBSCAN is per-epitope. Re-estimating the radius on a per-epitope $n$ of 30–300
is the degenerate regime a knee fails in (§5), so the geometry is estimated once per chain and the
clustering runs per epitope (ROADMAP §8.4).

The frontier is monotone, so tuning is a bisection rather than a search: lift falls as `eps` widens
and `Q` rises, so the best admissible cell is the smallest `eps` clearing the `Q` bar, and both
chains therefore sit on the boundary by construction, TRB clearing by 0.0004 and TRA's next step
down missing by 0.0010. `min_cluster` 5 dominates 3 on lift, purity and precision at every `coef`
measured, so there is no trade-off between them: clusters of 3–4 contribute low-confidence members,
not signal.

---

## 4. HDBSCAN - rejected

`vdjdb.motifs.tcremp.hdbscan_labels`, `sklearn.cluster.HDBSCAN`, with no dependency beyond the
sklearn already used for the PCA and DBSCAN.

### 4.1 Motivation

DBSCAN commits to one radius for every epitope, and VDJdb's epitopes do not share a density: a
display-selected library yields thousands of receptors one substitution apart by construction (29,688
of 192,753 records) while a 30-record epitope is sparse. HDBSCAN condenses a cluster hierarchy over
mutual-reachability distance and selects by stability, so there is no global radius to be wrong, and
it folds `min_samples` and the `min_cluster` post-filter into one `min_cluster_size`.

### 4.2 Parameters

| Parameter | Swept | Effect |
|---|---|---|
| `min_cluster_size` | 3, 5, 10, 20 | smallest cluster the hierarchy will emit, and HDBSCAN's primary knob |
| `cluster_selection_method` | `eom`, `leaf` | excess-of-mass takes the most stable cut (parsimonious); `leaf` takes every condensed-tree leaf, which shatters by design |
| `min_samples` | 2, 5 | how conservative the core-distance estimate is; higher declares more points noise |
| `cluster_selection_epsilon` | not swept | the anti-shatter knob, and not usable; see below |
| `allow_single_cluster` | `False` (default) | left alone; `True` would let an epitope legitimately percolate |

`cluster_selection_epsilon` is broken in sklearn 1.9.1. On a 12-cell probe at 0.5 and 1.0 × the
chain's mean 1-NN distance, 5 cells raised
`TypeError: only 0-dimensional arrays can be converted to Python scalars` from
`sklearn/cluster/_hdbscan/_tree.pyx:588` (`traverse_upwards`), and the cells that did run returned
bit-identical results to `epsilon = 0.0`, so the one parameter that would directly counter
shattering is unavailable and the `leaf` rows below cannot be repaired.

### 4.3 Scorecard - 56 cells per chain, human, `min_cluster` fixed at 5

`min_cluster_size` × {3, 5, 8, 10, 15, 20, 30, 50}, `cluster_selection_method` × {`eom`, `leaf`},
`min_samples` × {2, 3, 5, 10}, cells with `min_samples > min_cluster_size` skipped: 56 per chain.
The full grid with every instrument is `docs/tuning/scorecard.tsv`; the extremes are below.

**TRA** - bar: `Q` 0.1691, purity 0.8658, precision 0.8567, coverage 103 epitopes.
Each row is the grid's extreme on one column; `p` is `Q`'s parsimony term.

| config | lift | f1 | q | p | purity | retention | epitopes | perc_med |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| leaf mcs5 ms5 | **1.397** | **0.0799** | 0.1876 | 0.104 | 0.9002 | 0.3344 | 114 | 0.377 |
| eom mcs3 ms2 | 1.389 | 0.0789 | 0.1700 | 0.093 | **0.9061** | 0.3244 | 114 | 0.355 |
| eom mcs5 ms2 | 1.257 | 0.0734 | 0.2616 | 0.152 | 0.9017 | 0.4366 | **116** | 0.306 |
| leaf mcs50 ms5 | 1.143 | 0.0669 | 0.3467 | 0.213 | 0.8869 | 0.4304 | 108 | **0.241** |
| eom mcs30 ms5 | 0.921 | 0.0548 | **0.6489** | 0.507 | 0.8833 | **0.5623** | 111 | 0.286 |

**TRB** - bar: `Q` 0.4433, purity 0.9790, precision 0.9756, coverage 103 epitopes.

| config | lift | f1 | q | p | purity | retention | epitopes | perc_med |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| leaf mcs3 ms3 | **1.538** | **0.0884** | 0.1346 | 0.072 | 0.9337 | 0.3478 | 147 | 0.455 |
| leaf mcs15 ms3 | 1.222 | 0.0743 | 0.3331 | 0.201 | 0.9326 | 0.6270 | 130 | **0.215** |
| eom mcs5 ms2 | 1.218 | 0.0744 | 0.4427 | 0.287 | 0.9389 | 0.7240 | **165** | 0.357 |
| eom mcs20 ms3 | 1.119 | 0.0689 | 0.6254 | 0.464 | **0.9433** | 0.8326 | 125 | 0.244 |
| eom mcs30 ms2 | 1.055 | 0.0649 | **0.6397** | 0.479 | 0.9424 | **0.8336** | 115 | 0.247 |

### 4.4 Scorecard results

Under the legacy-relative bar, 51 of 56 TRA cells are admissible and 0 of 56 TRB cells are: on TRB
purity tops out at 0.9433 against a bar of 0.9790, so every cell regresses it by 3.6 points. The 51
admissible TRA cells are not on their own a result, since §8 shows the do-nothing partition is
admissible on TRA too.

Under the absolute floor of 0.94 the two chains swap places, 0 of 56 on TRA and 11 of 56 on TRB.
TRA's highest HDBSCAN purity is 0.9061, so the floor excludes the grid there, and §8.1 shows it
excludes every other algorithm's TRA cells with it. Both verdicts are columns in `scorecard.tsv`
(`admissible`, `admissible_legacy_bar`); neither replaces the other.

The `eom` / `leaf` axis behaves as theory predicts: `leaf` raises lift and lowers parsimony, TRB
`leaf` parsimony spanning 0.068–0.286 against `eom`'s 0.148–0.500, while `eom` raises parsimony,
retention and coverage and gives up lift, and neither end is admissible on TRB.

HDBSCAN's strength is coverage and structure rather than selectivity. TRB `eom mcs5 ms2` reaches 165
of 178 epitopes at retention 0.7240 and median percolation 0.357, against legacy's 103 epitopes,
0.3218 and 0.769, and it leads on every axis except the one being optimised, its lift being 1.218
against shipped TCREMP's 3.586.

### 4.5 Under the absolute purity floor of 0.94

The floor §7 lists as an open question was measured; §8.1 gives the measurement behind choosing 0.94
on both chains rather than 0.93.

| gene | purity floor | HDBSCAN cells admissible of 56 | best admissible cell | lift | F1 | retention | epitopes |
|---|---|---:|---|---:|---:|---:|---:|
| TRA | legacy, 0.8658 | 51 | `leaf mcs5 ms5` | 1.397 | 0.0799 | 0.3344 | 114 |
| TRA | **0.9400** | **0** | -- | -- | -- | -- | -- |
| TRB | legacy, 0.9790 | 0 | -- | -- | -- | -- | -- |
| TRB | 0.9300 | 25 | `eom mcs5 ms3` | 1.194 | 0.0731 | 0.7367 | 160 |
| TRB | **0.9400** | **11** | `eom mcs8 ms3` | 1.155 | 0.0709 | 0.7698 | 155 |

At a floor of 0.94 HDBSCAN becomes admissible on TRB and loses stage 2 by 3.1×, lift 1.155 against
the shipped TCREMP's 3.586 and F1 0.0709 against 0.1888, so the floor can be set without changing
what ships. Rewriting stage 2 as well would give a different database: 155 of 178 epitopes carrying
motifs instead of 105, at 77.0 % retention instead of 37.5 %, with a third of legacy's percolation,
and a clustered clonotype 16 % more likely than chance to have been seen by a second laboratory
rather than 259 % more likely.

On TRA the same floor admits nothing, HDBSCAN included, since its best purity is 0.9061 against the
floor's 0.94. TRA's admissible reading stays the legacy-relative one, and §8.1 measures why no
absolute floor can replace it there. HDBSCAN is rejected on the objective rather than on the bar; it
stays wired and tested, and §8 records what the bar can and cannot do.

---

## 5. Choosing `eps` by a knee - wired, not in the operating path

`vdjdb.motifs.tcremp.knee`. The published TCREMP protocol picks `eps` from the knee of the sorted
k-distance curve; this repo does not.

The reference implementation is degenerate at production scale: Kneedle with a degree-10 polynomial
fit on pooled human TRB, 112,983 clonotypes, returns knee index 1, fraction 0.000. A degree-10 fit
over 113,000 points oscillates, and the reported knee is the first oscillation rather than a feature
of the data, which puts `eps` below the data entirely and retention falls to a few per cent.

`knee` states Kneedle's difference curve directly, with no polynomial:

$$x_i = \frac{i}{n-1}, \qquad y_i = \frac{c_i - \min c}{\max c - \min c}, \qquad
\text{knee} = \arg\max_i (y_i - x_i), \qquad \text{strength} = \max_i (y_i - x_i)$$

for a sorted curve $c$ (signs reversed for a convex curve). Three guards, each against an observed
failure:

| Guard | Default | What it catches |
|---|---|---|
| resample onto a fixed grid | 1,000 points | the degree-10 oscillation; makes the knee a property of the curve's shape, measured invariant at 0.799 / 0.800 / 0.800 for one elbow at $n$ = 300 / 30,000 / 300,000 |
| `min_strength` | 0.05 | a curve with no knee. $\max(y-x)$ is exactly 0 for a straight line, so the statistic that locates a knee is the same one that reports there is none, and no separate test is needed |
| `floor_frac` / `ceil_frac` | 0.05 / 0.95 | a knee pinned to an end. At the bottom `eps` falls below the data; at the top the epitope percolates |

`strength` is scale-invariant by construction, so it is comparable across chains with different
embedding scales (TRA mean 1-NN 9.66, TRB 9.19). What runs instead is
`eps = coef × mean(1st-NN distance)`, with `coef` fitted against independent replication, which is
not a knee method and must not be described as one (ROADMAP §8.3).

---

## 6. The per-epitope view

`vdjdb.validate.motif_bench.per_epitope`, counted on clonotypes so a clonotype reported by forty
papers does not weigh forty times in its own epitope's retention, and reported over epitopes, the
unit curation acts on. A pooled retention of 0.33 is a different database depending on whether every
epitope is a third clustered or a third of them are fully clustered and the rest not at all; on TRB
it is closer to the second, and no pooled column reports that.

### 6.1 The shipped clusterings, per epitope

**TRA - 118 epitopes** with ≥ 30 records:

| method | epitopes covered | none at all | retention median | retention q75 | ≥90 % percolated | percolation median | clonotypes clustered | cids | epitopes with lift > 1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| legacy | 103 | 15 | 0.0647 | 0.1603 | 24 | 0.4717 | 15,063 | 2,214 | 29 |
| TCRNET | 103 | 15 | 0.0706 | 0.1624 | **22** | **0.3958** | 15,480 | 2,173 | 33 |
| TCREMP | **104** | **14** | **0.0850** | **0.2115** | 23 | 0.4834 | **17,380** | **2,926** | **34** |

**TRB - 178 epitopes**:

| method | epitopes covered | none at all | retention median | retention q75 | ≥90 % percolated | percolation median | clonotypes clustered | cids | epitopes with lift > 1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| legacy | 103 | 75 | 0.0089 | 0.0941 | 42 | 0.7692 | 37,210 | 1,074 | 42 |
| TCRNET | **109** | **69** | **0.0142** | 0.0941 | 41 | 0.7059 | 37,109 | 1,051 | **43** |
| TCREMP | 105 | 73 | 0.0128 | **0.1139** | **40** | **0.6667** | **42,802** | **1,156** | 42 |

### 6.2 Per-epitope results

TRB is a two-population database and TRA is not. 75 of 178 TRB epitopes get no motif from the
shipped annotation, and the median epitope has 0.9 % of its clonotypes clustered against a pooled
retention of 32 %, the pooled figure coming from a minority of large, deeply-sampled epitopes such as
GILGFVFTL, NLVPMVATV and YLQPRTFLL. On TRA the median epitope has 6.5 % clustered against a pooled
21 %, still skewed but far less bimodal. A single coverage number for TRB describes about sixty
epitopes and implies 178.

Percolation is the dominant structural defect on TRB, not shattering. The median TRB epitope has
77 % of its clustered clonotypes in a single cluster, 42 of 178 are ≥ 90 % collapsed, and both
rewrites reduce it (TCRNET to 0.7059, TCREMP to 0.6667) without removing it. It is the target for a
future motif change, and the reason §2 measured Leiden: subdividing percolated components addresses
it, but the resolutions that subdivide enough also shatter the rest.

Coverage and pooled lift move in opposite directions, and that is the finding that changed the
shipped parameters. See `denoising.md` §7.1, where TCREMP TRB at `coef` 1.15 maximised pooled lift at
4.490 and covered 88 of 178 epitopes, fifteen fewer than the annotation it replaces; at 1.55 it
covers 105 and leads legacy by 25.6 % on lift.

`lift > 1` per epitope is the coverage-weighted read. Of the epitopes where a local base rate exists
at all, legacy beats chance in 42 of them on TRB and 29 on TRA, while the rewrites reach 43 (TCRNET)
and 42 (TCREMP) on TRB and 33 and 34 on TRA.

### 6.3 The report

`out/reports/motifs_per_epitope.tsv`, 888 rows, one per (gene, method, epitope), written by the
build and uploaded as a CI artifact, never shipped in a bundle (`outputs.md` §5). Columns:
`clonotypes`, `clustered`, `retention`, `clusters`, `largest_cluster`, `mean_cluster_size`,
`singleton_clusters`, `percolation`, `replicated`, `tp`, `precision`, `lift`.

`lift` is null, not zero, where an epitope has no replicated clonotype: there is no base rate to
lift against, and averaging a zero there would put the per-epitope mean below the pooled figure. 47
of 118 TRA epitopes and 60 of 178 TRB epitopes have a defined lift.

---

## 7. Summary

| Algorithm | Status | TRA lift | TRA F1 | TRB lift | TRB F1 | TRB epitopes | Distinguishing property |
|---|---|---:|---:|---:|---:|---:|---|
| the shipped 2026-06-03 annotation | baseline | 1.703 | 0.0932 | 2.855 | 0.1489 | 103 / 178 | the source of every bar |
| connected components | **shipped**, TCRNET | 1.977 | 0.1086 | 3.194 | 0.1669 | **109** | no parameters to overfit |
| CPM Leiden | rejected | -- | -- | -- | -- | -- | lift 2.6× on TRA at parsimony ÷ 8.4 |
| DBSCAN | **shipped**, TCREMP | 1.855 | 0.1032 | **3.586** | **0.1888** | 105 | monotone frontier; tuning is a bisection |
| HDBSCAN | rejected | 1.397 | 0.0799 | 1.155 | 0.0709 | 155 at floor 0.94 | coverage 160/178 at retention 0.74 |
| Lumbermark | measured | -- | -- | 1.148 | 0.0706 | 175 | percolation 0.224 vs legacy's 0.769 |
| TCRNET gate + Lumbermark | measured | **3.618** | **0.1702** | 4.089 | 0.1943 | 82, below the floor | best objective measured; fails coverage |

A dash is where no admissible configuration exists, so the rejected rows show the best admissible
number or nothing. The two measured rows show their distinguishing configuration instead; for
Lumbermark that is not one of the six TRB cells that clear the 0.94 floor, which sit at lift
1.009–1.069 (§9.3). Across 252 configurations and six algorithms, under both the legacy-relative bar
and the absolute 0.94 floor, nothing admissible beats TCREMP at `coef` 1.55 on TRB or TCRNET on TRA,
and §§8–9 record the supporting measurements.

### What would change a verdict

1. ~~A bar that is absolute rather than relative.~~ Settled at 0.94 on both chains, and it changes
   nothing that ships; measured in §4.5 and §8.1. On TRB the floor takes stage 1 from 1 admissible
   configuration of 122 to 23, namely 11 `hdbscan-eom`, 6 Lumbermark, 3 wider DBSCAN radii and 2
   `hybrid-len-all`, and stage 2 rejects all 22 newcomers, because the shipped `coef` 1.55 still has
   the highest admissible lift at 3.585 and the selected radius does not move on either chain. On
   TRA the floor admits nothing at all, including both shipped methods; §8.1 measures the window as
   empty by 0.0030 of purity and states what guards TRA instead.
2. ~~Per-epitope coverage entering the criterion.~~ Settled as a fourth admissibility axis
   (`denoising.md` §7.1). Measured, TCREMP TRB at `coef` 1.15 covered 88 of 178 against legacy's 103;
   the shipped TRB radius moved to 1.55 as a result.
3. A gate with the recruited set's coverage and the enriched set's selectivity. §9 shows the gate
   and the partition are separable, and that the gate is what moves the objective, TRA F1 0.1086 →
   0.1702 on the enriched gate alone; nothing measured keeps both, and that is the open question
   this document leaves.
4. `cluster_selection_epsilon` becoming usable in a later sklearn, the one knob that would give
   HDBSCAN's `leaf` lift without its shattering.

A rejected algorithm stays implemented and tested, so the measurement that made it non-default can
be re-run.

---

## 8. Instrument audit - the do-nothing partition

Every sweep above includes one extra row, `trivial`: one cluster per epitope, every clonotype in it,
nothing excluded. It is not an algorithm but the reading an instrument gives when no clustering has
happened at all, and §7.1's admissibility axes do not survive it.

| gene | axis | do-nothing partition | legacy bar | does the bar exclude it? |
|---|---|---:|---:|---|
| TRA | `Q` | 0.7940 | 0.1691 | **no** |
| TRA | purity | 0.9204 | 0.8658 | **no** |
| TRA | precision | 0.9237 | 0.8567 | **no** |
| TRA | epitope coverage | 118 | 103 | **no** |
| TRA | lift | 1.0000 | 1.7027 | yes |
| TRA | F1 | 0.0600 | 0.0932 | yes |
| TRB | `Q` | 0.8632 | 0.4433 | **no** |
| TRB | purity | 0.9343 | 0.9790 | yes |
| TRB | precision | 0.9455 | 0.9756 | yes |
| TRB | epitope coverage | 178 | 103 | **no** |
| TRB | lift | 1.0000 | 2.8552 | yes |
| TRB | F1 | 0.0622 | 0.1489 | yes |

`vdjdb.validate.motif_bench.trivial_members` builds it, so the audit re-runs from the package and
not from a sweep script. Two forms of it exist, and every shipped clustering is split by CDR3 length
before it reaches a release (§0), so the split form is the one a bar has to clear; the unsplit form,
one cluster per epitope, scores higher on `Q` (0.9208 on TRA, 0.9654 on TRB) and lower on purity
(0.9144, 0.9253). A floor of 0.93 would exclude the unsplit partition and admit the split one; the
table above is the split form.

Three consequences:

`Q` is not an admissibility axis, since doing nothing scores `Q` 0.7940 on TRA and 0.8632 on TRB,
higher than any clustering measured here on either chain across 252 configurations. Its parsimony
term charges for shattering, and the partition that shatters least is the one with one cluster per
class, so `Q` discriminates between clusterings of comparable retention and says nothing across
retentions. It stays in the scorecard as a shatter guard and comes out of the admissibility rule.

Epitope coverage behaves the same way. Coverage is maximised by claiming every clonotype, so a
coverage floor is a floor against narrowness rather than a quality bar; it still serves the purpose
§7.1 added it for, stopping stage 2 from raising lift by abandoning epitopes, but it is not evidence
that a clustering is good.

On TRA the four-axis rule already admits the do-nothing partition. Legacy TRA purity is 0.8658,
below the 0.9204 that doing nothing scores, so `trivial` clears all four axes today, and what
excludes it is stage 2: its lift is exactly 1.000 by construction and its F1 sits at the floor.

### 8.0 It also beats the shipped clustering at reproducing the shipped clustering

The tables a release ships were never compared against the previous release's, because `vdjdb diff`
keys `cluster_members.txt` on `cid` and a cid is `<species>.<chain>.<epitope>.<n>` with `n` a
position in a sorted list, so one renumbered cluster reads as the entire file replaced. Since
2026-09-29 `vdjdb motif-metrics` measures the partition instead
(`vdjdb.compare.clustering.agreement`): of every clonotype the 2026-06-03 release clustered, what
fraction of its cluster-mates does a source still give it?

| gene | source | clonotypes the release clustered | also placed here | cluster-mates kept | co-clustered pairs kept | ARI |
|---|---|---:|---:|---:|---:|---:|
| TRA | TCRNET (shipped) | 15,016 | 12,011 | 0.7024 | 0.7194 | 0.9945 |
| TRA | TCREMP | 15,016 | 8,530 | 0.3430 | 0.4039 | 0.8161 |
| TRA | the do-nothing partition | 15,016 | 13,108 | **0.7490** | 0.7345 | 0.1872 |
| TRB | TCRNET (shipped) | 36,906 | 35,095 | 0.9224 | 0.9952 | 0.9999 |
| TRB | TCREMP | 36,906 | 32,047 | 0.7921 | 0.9966 | 0.9977 |
| TRB | the do-nothing partition | 36,906 | 35,649 | **0.9394** | 0.9991 | 0.9939 |

So doing nothing reproduces the released clustering's neighbourhoods better than either shipped
method does, on both chains. That is a statement about the released clustering, not about ours: it
puts **19,971 of its 36,906 TRB clonotypes in one cluster**, `H.B.SLLMWITQV.1`, which is 54 % of the
chain's clonotypes and **199,410,435 of the file's 210,634,576 co-clustered pairs, 94.7 %**. A
partition that groups every clonotype of an epitope reproduces that blob exactly.

Two consequences, the same shape as the three above.

**Agreement with the last release is not an admissibility axis.** It is gated against its own
committed baseline, not against `latest` and not against `trivial`, for the reason every other row of
§8 gives: the do-nothing partition clears it.

**A pair-weighted agreement is a measurement of one cluster.** `partition_pairs_preserved` and the
adjusted Rand index read 0.9991 and 0.9939 for doing nothing on TRB, while the per-clonotype figure
separates it from TCRNET by 1.7 points. Both are recorded, because the blob is worth watching, and
only the per-clonotype one is gated.

### 8.1 The purity floor - 0.94 on both chains

On TRB the legacy purity bar of 0.9790 is the one axis the do-nothing partition fails, and relaxing
it to 0.93 removes that, because doing nothing scores 0.9343, 0.0043 above the proposed floor. An
absolute floor has to clear the measured do-nothing purity per chain rather than be a round number
chosen in advance.

The floor is 0.94 on both chains (`docs/tuning/sweeps.py`, `BAR`), the lowest round number above the
do-nothing purity of either chain, 0.9204 on TRA and 0.9343 on TRB. Every cell in `scorecard.tsv`
has two verdicts, `admissible` under that floor and `admissible_legacy_bar` under the
legacy-relative one, so neither reading has to be recomputed to be checked.

| TRB purity floor | HDBSCAN cells admissible of 56 | do-nothing partition admitted? |
|---|---:|---|
| legacy, 0.9790 | 0 | no |
| 0.9300 | 25 | **yes** |
| 0.9343 | 23 | **yes** (the floor equals its purity) |
| **0.9400** | **11** | **no** |

#### Floor sweep and admissible window, 122 configurations per chain

Both tables are generated by `docs/tuning/report.py` from `scorecard.tsv`; a cell is admissible at an
absolute floor when purity and precision clear it and `Q` and epitope coverage clear legacy's.

| gene | purity floor | cells admissible | of | do-nothing admitted? | best by lift | lift | F1 | retention | epitopes |
|---|---|---:|---:|---|---|---:|---:|---:|---:|
| TRA | 0.8658 (legacy) | 70 | 122 | **yes** | `dbscan coef 1.8` | 1.854 | 0.1032 | 0.2509 | 104 |
| TRA | 0.9204 | 0 | 122 | **yes** | -- | -- | -- | -- | -- |
| TRA | 0.9300 | 0 | 122 | no | -- | -- | -- | -- | -- |
| TRA | **0.9400** | **0** | 122 | no | -- | -- | -- | -- | -- |
| TRB | 0.9300 | 37 | 122 | **yes** | `dbscan coef 1.55` | 3.585 | 0.1887 | 0.3747 | 105 |
| TRB | 0.9343 | 35 | 122 | **yes** | `dbscan coef 1.55` | 3.585 | 0.1887 | 0.3747 | 105 |
| TRB | **0.9400** | **23** | 122 | no | `dbscan coef 1.55` | 3.585 | 0.1887 | 0.3747 | 105 |
| TRB | 0.9790 (legacy) | 1 | 122 | no | `dbscan coef 1.55` | 3.585 | 0.1887 | 0.3747 | 105 |

The floor takes TRB from one admissible configuration to twenty-three, and stage 2 selects the same
one: under the legacy bar exactly one candidate clears stage 1, `dbscan coef 1.55`, which is what
ships, and at 0.94 it is joined by 11 `hdbscan-eom` cells, 6 Lumbermark, 3 wider DBSCAN radii and 2
`hybrid-len-all`, each of which loses on lift. The relaxation adds 22 candidates and changes nothing.

A floor is usable only if it excludes the do-nothing partition and admits something, so its ceiling
is the highest purity reached by any configuration that already clears the other two axes, `Q` and
epitope coverage, both legacy-relative.

| gene | do-nothing purity | cells clearing `Q` and coverage | highest purity among them | window | width |
|---|---:|---|---:|---|---:|
| TRA | 0.9204 | 74 of 122 | 0.9174 | **empty** | −0.0030 |
| TRB | 0.9343 | 37 of 122 | 0.9829 | **[0.9343, 0.9829]** | 0.0486 |

On TRB the window is 0.0486 of purity wide and 0.94 sits inside it with margin at both ends, 5.4×
the [0.9343, 0.9433] window HDBSCAN alone offers, because the incumbent DBSCAN radius sits at 0.9829
and the ceiling is set by the best method rather than by the method being tested against it.

On TRA the window is empty by 0.0030: across six algorithms and 122 configurations, the highest
purity any cell reaches while clearing `Q` and coverage is 0.9174 (`hybrid-len-all mcs5 M1`), below
the 0.9204 that doing nothing scores. Every absolute purity floor that excludes the do-nothing
partition on TRA therefore also excludes every configuration that clears the other two axes,
including both shipped methods, whose TRA purity is 0.8754 (TCRNET) and 0.9104 (TCREMP), and at 0.94
specifically 0 of 122 TRA cells are admissible. Purity is therefore not the axis that separates
signal from nothing on TRA, and §8's third consequence identifies stage 2 as the axis that is; TRA
keeps `admissible_legacy_bar` as its usable reading, and the 0.94 column records what an absolute
floor would cost there.

---

## 9. Lumbermark and the TCRNET-gated hybrid - measured, not shipped

### 9.1 Motivation

§6.2 measured the defect: TRB's dominant failure is percolation, not shattering, with the median
epitope putting 0.706 of its clustered clonotypes in one cluster under the shipped TCRNET and 41 of
178 epitopes ≥90 % collapsed. DBSCAN percolates because one radius chains through density bridges;
Leiden (§2) and HDBSCAN `leaf` (§4) break the bridges by shattering, which is the other failure.

Lumbermark (Gagolewski 2026, [arXiv:2604.07143](https://arxiv.org/abs/2604.07143), `lumbermark` on
PyPI, AGPL-3) addresses that gap. It cuts the *k*−1 longest edges of the M-mutual-reachability
minimum spanning tree, HDBSCAN's own internal structure, subject to every resulting component being
at least `min_cluster_size`, after pruning the tree's leaves. Mutual reachability pulls low-density
points apart so a bridge becomes a long edge, the leaf pruning puts the cut on a protruding edge
rather than an outlier's stalk, and the size floor stops the cut sequence shattering. It beats Genie
and HDBSCAN on 61 reference datasets by adjusted Rand index. Requesting *k* = *n*−1 makes
`min_cluster_size` the only granularity knob: the algorithm saturates and warns rather than failing,
so there is no cluster count to choose.

It has no noise label: it partitions every point, and the paper lists outlier detection as future
work. VDJdb already supplies that half, since TCRNET's enrichment test against a matched background
repertoire is a statistical noise model rather than a density heuristic. The variants measured each
differ from their comparator in exactly one factor:

| variant | vertex set | partition | isolates |
|---|---|---|---|
| `lumbermark` | everything | Lumbermark | the partition's ceiling with no gate |
| `hybrid` | TCRNET **enriched** | Lumbermark | the gate, against `lumbermark` |
| `hybrid-recruited` | TCRNET **enriched + neighbours** | Lumbermark | the partition, against shipped TCRNET |
| `hybrid-len` | as above, within CDR3-length strata | Lumbermark | the emit path's length split |

`hybrid-recruited` and the shipped TCRNET partition the same vertices, so the only difference between
them is the partition step.

### 9.2 Measurements

Lumbermark gives the lowest percolation measured in this document. TRB `mcs10 M1`: median
per-epitope percolation 0.224 against legacy's 0.769 and shipped TCREMP's 0.667, at retention 0.8345
over 175 of 178 epitopes, purity 0.9445, which is the highest purity of any configuration that
retains four-fifths of the cohort (23 of the 115 TRB cells clear retention 0.80, and it leads all of
them), above HDBSCAN's best of 0.9433 and above the do-nothing partition's 0.9343. That applies to
the inclusive regime only, since DBSCAN reaches purity 0.9920 on TRB at retention 0.1189, the trade
this document is about, and Lumbermark sits at the other end with lift 1.148.

The gated hybrid gives the highest selectivity measured in this document. On TRA, `hybrid mcs3 M5`
reaches F1 0.1702 and lift 3.618, against the shipped TCRNET's F1 0.1086 and lift 1.977 and legacy's
F1 0.0932 and lift 1.703, which is +57 % F1 over the best shipped TRA configuration on the objective
§11.1 tunes against, and the largest single improvement measured on either chain.

| gene | configuration | lift | F1 | `Q` | purity | retention | epitopes | percolation |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| TRA | legacy release | 1.703 | 0.0932 | 0.1691 | 0.8658 | 0.2105 | 103 | 0.472 |
| TRA | shipped TCRNET | 1.977 | 0.1086 | 0.1798 | 0.8754 | 0.2323 | 103 | 0.396 |
| TRA | `hybrid-recruited mcs3 M1` | 2.368 | 0.1244 | 0.0759 | 0.9025 | 0.1816 | 97 | 0.375 |
| TRA | **`hybrid mcs3 M5`** | **3.618** | **0.1702** | 0.0452 | 0.9040 | 0.1255 | 86 | 0.500 |
| TRB | legacy release | 2.855 | 0.1489 | 0.4433 | 0.9790 | 0.3218 | 103 | 0.769 |
| TRB | shipped TCRNET | 3.194 | 0.1669 | 0.4428 | 0.9790 | 0.3337 | 109 | 0.706 |
| TRB | shipped TCREMP | 3.586 | 0.1888 | 0.4947 | 0.9829 | 0.3745 | 105 | 0.667 |
| TRB | `hybrid-recruited mcs3 M1` | 3.234 | 0.1649 | 0.1110 | 0.9786 | 0.2864 | 102 | 0.500 |
| TRB | **`hybrid mcs5 M10`** | **4.089** | **0.1943** | 0.1388 | 0.9816 | 0.2847 | 82 | 0.527 |

Three readings:

The gate moves the objective and the partition moves percolation. `hybrid-recruited` against the
shipped TCRNET - same vertices, Lumbermark instead of connected components - moves TRB percolation
0.706 → 0.500 and leaves F1 flat (0.1669 → 0.1649) while `Q` falls 0.4428 → 0.1110 on 4,702 cids
against 657, so swapping the partition alone de-percolates and costs parsimony, and tightening the
gate from recruited to enriched is what moves lift and F1.

Both hybrids fail the coverage axis, by 17 epitopes on TRA and 23 on TRB, which is the narrow corner
§7.1's coverage floor exists to exclude.

Clustering inside CDR3-length strata does not recover the coverage. `hybrid-len` was run because the
emit path splits every cluster by length before applying `min_cluster` (§0), so a cluster of 8 across
three lengths is removed entirely at the floor of 5 and clustering within a stratum should have
recovered it. It recovers some coverage and gives back the same amount of objective, TRA 86 → 93
epitopes at F1 0.1702 → 0.1436 and TRB 82 → 91 epitopes at F1 0.1943 → 0.1832.

### 9.3 Status

Neither is shipped, and the reason is coverage rather than quality. The gated hybrid is inadmissible
under both rules on both chains, 0 of 15 cells each for `hybrid` and `hybrid-recruited` under the
legacy bar and under the 0.94 floor alike, because the enriched gate abandons epitopes: the two
variants that clear the objective cover 82–97, and `hybrid-recruited` tops out at 102 on both chains,
one epitope short of the 103 floor.

Lumbermark is the exception, on TRB and under the new floor only. Six of its 15 TRB cells, `mcs20`
and `mcs50` at every `M`, clear all four axes at 0.94, with purity 0.9405–0.9439, `Q` 0.4800–0.5681
against legacy's 0.4433, and 124–156 epitopes against a floor of 103. Stage 2 then rejects all six at
lift 1.009–1.069 against shipped TCREMP's 3.586, the same shape as HDBSCAN, where the floor admits
the inclusive regime and the objective rejects it.

The gate and the partition are separable and were measured separately. The shipped methods each
couple one representation to one partition, TCRNET to connected components on the Hamming-1 graph and
TCREMP to DBSCAN in the embedding, while the hybrid shows the two choices are independent, that
VDJdb's best noise model is its enrichment test rather than any density heuristic, and that the
embedding's MST is where percolation can be fixed. A method that kept the recruited gate's coverage
and the enriched gate's selectivity would beat everything in this document, and nothing measured here
does both.

One thing deliberately not tried: admitting clusters individually by the within-stratum
sum-of-hypergeometrics null in `vdjdb.validate.noise`. Admitting a cluster because its members are
enriched for independent replication, and then scoring the result by independent-study lift, fits on
the test statistic, so the apparent improvement would not be informative.

---

## 10. Related methods in the literature

The TCR clustering literature is large and mostly concerns similarity, which is the half of this
problem VDJdb does not have trouble with.

| tool | representation | partition | noise model | relation to this document |
|---|---|---|---|---|
| GLIPH (Glanville et al. 2017, *Nature* 547:94–98, [doi](https://doi.org/10.1038/nature22976)) | CDR3 motif enrichment against a reference set, plus global similarity | specificity groups | enrichment *p* | the idea TCRNET implements, with a matched background rather than a naive reference |
| GLIPH2 (Huang et al. 2020, *Nat Biotechnol* 38:1194–1202, [doi](https://doi.org/10.1038/s41587-020-0505-4)) | as above, scaled to millions | greedy grouping | enrichment *p* | a scaling result; VDJdb's *n* is 10⁵ |
| GIANA (Zhang et al. 2021, *Nat Commun* 12:4699, [doi](https://doi.org/10.1038/s41467-021-25006-7)) | isometric embedding of the CDR3 | nearest-neighbour hashing | none | 600× TCRdist's speed at equal specificity; a scaling result |
| clusTCR (Valkiers et al. 2021, *Bioinformatics* 37:4865–4867, [doi](https://doi.org/10.1093/bioinformatics/btab446)) | Hamming-1 graph inside *k*-means superclusters | Markov clustering (MCL) | none | the one partition on the graph side not measured here; see below |
| Vujovic et al. 2020, *Comput Struct Biotechnol J* 18:2166–2173, [doi](https://doi.org/10.1016/j.csbj.2020.06.041) | -- | -- | -- | the survey of the above; confirms the pattern |

*(Bibliographic records retrieved from PubMed.)*

Each of them clusters and none of them selects: they answer which receptors are similar, so their
published evaluations are purity-like and their parameters are similarity thresholds, while VDJdb's
motif stage answers which similarity is evidence, so the objective here is independent-study
replication (§11.1) and what §9 measures is a gate plus partition composition rather than a new
distance.

MCL is deliberately left untested: it is the clusTCR partition and it resists percolation by flow
simulation, so it is the natural candidate on the graph side, but §2 already measured CPM Leiden on
that same enriched graph and found the failure is not the partition's granularity, homogeneity rising
while parsimony falls 8.4-fold on TRA and TRB gaining 1.1 % at best. MCL's inflation parameter is a
smoother knob than CPM's resolution rather than a different mechanism, and it is recorded here as an
open axis with a stated prior.

The clustering-methods literature supplied the one idea that did change a measurement,
mutual-reachability MSTs with a size-floored longest-edge cut (§9), from the line of work running
Campello's HDBSCAN through Genie (Gagolewski et al. 2016, *Inf Sci* 363:8–23) to Lumbermark. `Q`
comes from the same line of work, Tiffeau-Mayer,
[arXiv:2607.20799](https://arxiv.org/abs/2607.20799).

---

## 11. Figures

Rendered from the committed tables by `gnuplot`; regenerate with the commands in each script header.

| figure | script | shows |
|---|---|---|
| `out/reports/tuning/tuning_TRA.svg`, `tuning_TRB.svg` | `docs/tuning/tuning.gp` | all six instruments against retention, every configuration, both chains |
| `out/reports/tuning/per_epitope.svg` | `docs/tuning/per_epitope.gp` | the per-epitope retention and percolation ECDFs behind §6 |

The first figure shows §8's measurements: lift and F1 fall monotonically as retention
rises, `Q` and coverage rise toward the grey do-nothing marker at retention 1.0, and purity is the
only panel where that marker sits below the measured clusterings.

### Reproducing the tables

```bash
uv run python docs/tuning/sweeps.py all      # -> out/reports/tuning/*.csv   (~50 min)
uv run python docs/tuning/report.py          # -> docs/tuning/*.tsv, tables.md
gnuplot -e "chain='TRB'" docs/tuning/tuning.gp
```

`docs/tuning/scorecard.tsv` is the full 252-row measurement, one row per configuration and every
instrument, a committed and reviewed input to this document, refreshed by its own pull request and
never written by a build (hard rule 9).
