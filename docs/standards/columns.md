# Column reference

Every table below is built from `vdjdb.schema` while this page renders. There is no
generated file to refresh, so no copy of these tables can disagree with the build.


## Representative parents and completed MHC pairs

`antigen.gene` and `antigen.species` give one representative parent label per literal
peptide. A verified source identification resolves an incorrect or missing label first.
For remaining alternatives, the reviewed epitope dictionary selects the most frequent
complete `(parent species, parent gene)` pair across distinct observation records.
Known synthetic-library ligands use `Synthetic`; a synthesized natural peptide retains
its reported biological parent. Gene and species are selected together. Equal counts
use alphabetical species, then gene, as the deterministic tie-break.

For epitope $e$, let $n_e(s,g)$ count records with nonmissing parent species $s$ and
gene $g$. The representative pair is

$$
(s_e,g_e)=\operatorname*{arg\,max}_{(s,g)} n_e(s,g).
$$

The source chunks retain reported alternatives; the dictionary supplies the canonical
representative in assembled outputs. This label does not establish a unique biological
origin for a peptide shared between proteins or organisms.

MHC completion preserves a reported chain and first recovers a partner from the
publication or construct description. Otherwise, the most commonly used compatible
partner is imputed. Let $C(b)$ contain compatible alpha alleles for known beta allele
$b$, restricted to experimentally supported alternatives when the source supplies them.
Let $N(a,b)$ count curated observation records with that pair, combining additional
HLA allele fields at two-field protein resolution. Choose

$$
\hat a(b)=\operatorname*{arg\,max}_{a\in C(b)}N(a,b).
$$

The same rule completes a missing beta chain conditional on its known alpha chain.
When no conditional pair has support, use the most frequent partner within the
compatible DR, DQ, DP, I-A or I-E family; equal counts use alphabetical allele order.
Murine I-A/I-E shorthand follows the molecule labels used in the corpus. Known
source errors are corrected before counting. For example, PMID35914171 specifies
DQA1*01:02/DQB1*06:02 constructs; its former DQA1*05:01 annotation must not determine
the mode. See the [primary construct acknowledgment](https://pmc.ncbi.nlm.nih.gov/articles/PMC9371724/).

Completed assignments are reviewed inputs committed in chunks. Each imputed partner
is marked in the curation `comment`, together with the submitted pair. It is an
annotation estimate, not an additional restriction assay or donor genotype. Reported
methods and the confidence-score algorithm remain unchanged. Source-recovered partners
are distinguished from imputed ones in the comments.

## `vdjdb.txt` - the release table

```{vdjdb-schema} vdjdb
:columns: name, title, comment
```

## `vdjdb.txt` - the vdjdb-web variant

The variant served by `vdjdb-web` adds five `evidence.*` columns and a `TCR_hash` row that
the release table does not have.

```{vdjdb-schema} vdjdb-web
:columns: name, title, comment
```

## `vdjdb.slim.txt`

```{vdjdb-schema} slim
:columns: name, title, comment
```

## `vdjdb_full.txt`

```{vdjdb-schema} full
:columns: name, title
```

## The definitive tables

`records` and `chains` are the database. Every shipped file is a join and a pivot away
from them.

### `records`

```{vdjdb-schema} records
:columns: ships_as, title, comment
```

### `chains`

```{vdjdb-schema} chains
:columns: ships_as, title, comment
```

### `epitopes`

```{vdjdb-schema} epitopes
:columns: ships_as, title
```

### `restriction`

```{vdjdb-schema} restriction
:columns: ships_as, title
```

### `evidence`

```{vdjdb-schema} evidence
:columns: ships_as, title
```

### `epitope_assessment`

```{vdjdb-schema} epitope_assessment
:columns: ships_as, title
```

## The motif files

⚠ **These two files are parsed positionally by `vdjdb-web`**, which hands Tablesaw a fixed
column-type array with no header check. An inserted, removed or reordered column mistypes
or shifts every column after it, without an error. Column order is a contract.

### `cluster_members.txt`

```{vdjdb-schema} cluster_members
:columns: name, title
```

### `motif_pwms.txt`

```{vdjdb-schema} motif_pwms
:columns: name, title
```
