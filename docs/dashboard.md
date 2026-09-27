# The summary dashboard

The dashboard summarises the current contents of `chunks/`: records by species and chain, growth by
publication year, and breakdowns by antigen, antigen origin and MHC allele. It is served in two
places:

- `vdjdb-web` injects a fragment of it into the `/overview` page at [vdjdb.com](https://vdjdb.com).
- This page shows the same fragment below, rebuilt by CI on every push to `master`, so it
  reflects records added since the last release.

## What is shown, and how each number is computed

Eleven blocks, read top to bottom on the page: five tables and eight figures. Every one is computed
in `summary/vdjdb_summary.Rmd` from the legacy projection of the current build, and nothing else.
This section is the panel-by-panel reading of that script, so a number on the page can be traced to
the rows it came from.

### Conventions that apply to more than one panel

**Two source files, and they count different things.** Nine of the eleven blocks read
`vdjdb.slim.txt`; the by-year figure reads `vdjdb_full.txt`.

| | `vdjdb.slim.txt` | `vdjdb_full.txt` |
|---|---|---|
| A row is | one clonotype against one epitope | one curated record, both chains in one row |
| The same CDR3 in five donors | one row | five rows |
| V/J | one representative call per clonotype | as curated |
| `vdjdb.score` | the best across the merged records | per record |
| Panels using it | all but the by-year figure | the by-year figure only |

So the `Records` column of the first table is a clonotype count, not a record count, and it is
smaller than the record count in the release notes. The script's own caption says this; the
distinction is repeated here because it is the single most common misreading of the page.

**"Studies" counts every reference on the row.** Four blocks report a `Studies` column, computed as
`length(unique(unlist(strsplit(reference.id, ","))))`. A slim row merged from several publications
carries them comma-joined in `reference.id`, and 5,630 of 197,670 slim rows (2.8 %) do.

Until 2026-09-27 the column kept only field 1 of that split, so a merged row reported its first
publication and dropped the rest: **527 of the 636 distinct non-blank references** in the slim table,
missing 109. Rendering the same build both ways gives the size of the correction per row of the first
table:

| Species | Chain | Field 1 only | Every field |
|---|---|---|---|
| HomoSapiens | TRA | 348 | 401 |
| HomoSapiens | TRB | 431 | 483 |
| MacacaMulatta | TRB | 2 | 2 |
| MusMusculus | TRA | 71 | 84 |
| MusMusculus | TRB | 88 | 151 |

`vdjdb.assemble.evidence.support_counts` counts distinct `reference.id` on `records`, where the column
holds one value per row, and the R side now reads the same quantity.
`tests/unit/test_summary.py` fails if the truncating form comes back.

**Allele truncation.** Panels that group by MHC allele cut it at the first colon, so
`HLA-A*02:01:48` and `HLA-A*02:01` are one group. Panels that group by V gene cut at the `*`, so
`TRBV9*01` and `TRBV9*02` are one group.

**Comma-separated cells are split into rows, then counted by record.** `mhc.a`, `mhc.b` and `v.segm`
can each hold several calls. The heatmap and the HLA table expand them with `separate_rows`, which
multiplies a row with two alleles and two V calls into four, then count `length(unique(id))` against
a row number taken before the expansion. A record with three alleles is therefore counted once in
each of the three alleles' cells and three times in the page total, which is what "records supporting
this allele" means and not what "records" means elsewhere.

### 1. Render stamp

`Last updated on <date>`, from `Sys.time()` at render time.

This is the one part of the fragment that is not a function of the data: the same build rendered on
two days produces two different bytes. It is a declared exception to hard rule 7 rather than an
oversight, because the page's purpose includes saying how fresh it is. `vdjdb diff` never sees it:
the fragment is written to the build root and not to the legacy directory, so it is outside the
compared file set.

### 2. Record statistics by species and TCR chain (table)

Source `vdjdb.slim.txt`, grouped by `(species, gene)`.

| Column | Computed as |
|---|---|
| Records | rows in the group |
| Paired records | rows with `complex.id != "0"`, that is, a chain whose partner chain is also curated |
| Unique epitopes | distinct `antigen.epitope` |
| Studies | distinct first-listed `reference.id`, see above |

No filter, so every species in the build appears. The current build has three: `HomoSapiens`,
`MusMusculus` and `MacacaMulatta`. `RattusNorvegicus` is in the vocabulary and has no records.

### 3. Record statistics by year (figure, four panels)

The only panel reading `vdjdb_full.txt`, because it counts what a publication added and a slim row
has already merged publications together.

1. Drop `MacacaMulatta`.
2. Build two composite keys per record: `tcr_key` is `v.alpha j.alpha cdr3.alpha v.beta j.beta
   cdr3.beta` joined by spaces, `mhc_key` is `mhc.a mhc.b`. `chains` is `paired` when both CDR3s are
   present, otherwise `TRA` or `TRB`.
3. Reduce to distinct `(reference.id, tcr_key, mhc_key, chains, antigen.epitope, species)`.
4. Join the publication year from `summary/reference_years.tsv`. A reference with no year stops the
   render rather than dropping out of the plot.
5. For each of the four metrics - `tcr_key`, `antigen.epitope`, `reference.id`, `mhc_key` - and each
   `chains` group, take each key's **earliest** year, count keys per year, and cumulate.

So a curve is "how many distinct X had been reported by the end of year Y", and a TCR reported again
in a later paper does not lift the curve twice. The four panels are unique TCRs, unique epitopes,
number of studies, and number of MHC alleles; colour is the chain group.

Callouts come from `summary/annotations.tsv`, which carries `panel, year, label, hjust, vjust` and no
coordinates. The segment is drawn from zero to the series value in that year and the label floats
5 % of the panel maximum above it, so the annotation stays correct as the database grows. That is
issue #460's second half; the hardcoded y positions it replaced were right for a 2022 database.

The cumulation is a first-year tabulation rather than the cross join it replaced. The earlier form
paired every year with every other year, about 7 million rows over roughly 35 distinct years, and
was the only part of the render that could plausibly exhaust a 16 GB runner.

### 4. Summary by antigen and antigen origin (table)

Source `vdjdb.slim.txt`, `species == "HomoSapiens"`, grouped by `antigen.species` alone, sorted by
records descending. Columns: records, distinct epitopes, Studies. `antigen.gene` is deliberately not
in the grouping, so one row covers a whole source organism.

### 5. COVID-19 (figure and table)

Source `vdjdb.slim.txt`, `species == "HomoSapiens"` and `antigen.species` beginning `SARS-CoV`.

The MHC label is built before grouping: `mhc.a` and `mhc.b` are cut at the first comma or colon, and
for class II the two are joined as `<mhc.a>/<characters 7 to 15 of mhc.b>`, which drops the locus
prefix from the beta chain. Class I uses `mhc.a` alone.

The alluvial figure has four axes read left to right: `antigen.gene`, the HLA label without its
`HLA-` prefix, the **first three residues** of the epitope, and the TCR chain. Ribbon height is
`log2(records)`, so a stratum twice as tall holds four times the records. Ribbon colour indexes the
three-residue epitope stub and carries no other meaning.

The figure keeps groups with 30 or more records. The caption says "epitopes with less than 30 records
in total were not counted", but the filter is applied per `(gene, HLA, epitope, chain, studies)`
group rather than per epitope, so an epitope whose records are split across two HLAs can be dropped
from both while its total exceeds 30. The caption overstates what the filter does.

The table below the figure pivots the same frame so TRA and TRB become columns, and keeps rows where
`TRA + TRB >= 10`.

### 6. Self-antigen (figure and table)

Identical to the COVID-19 block in every step, with `antigen.species` beginning `HomoSapiens`, a
threshold of 10 rather than 30, and the `Accent` palette rather than `Set3`. The caption's "at least
10 records" carries the same per-group reading as above.

### 7. Distribution of VDJdb confidence scores (figure)

Source `vdjdb.slim.txt`, `species == "HomoSapiens"`, grouped by `(mhc.class, gene, vdjdb.score)`.
Bars are dodged within each `mhc.class gene` combination and the y axis is `log10` of the row count,
so the visible difference between adjacent bars is a ratio.

The score in slim is the best across the records merged into that clonotype, so this is a
distribution over clonotypes at their highest-scoring report.

### 8. Spectratype (two figures)

Both from `vdjdb.slim.txt`, `species == "HomoSapiens"`.

**CDR3 length.** A histogram of `nchar(cdr3)` at binwidth 1, faceted by `gene`. `cdr3` is junction
space, both anchors included, so these lengths are two residues longer than an AIRR `cdr3_aa`. Fill
indexes the epitope ordered by epitope length, which reads as a rough class I to class II gradient
rather than as an epitope legend. The dotted line is a density estimate with `adjust = 3.0`, rescaled
to counts. The x axis is limited to 5 to 25 residues, and lengths outside that range are dropped from
the plot without a warning: measured on the current build, 19 human slim rows, over a full range of 4
to 126 residues.

**Epitope length.** A histogram of `nchar(antigen.epitope)` at binwidth 1, faceted by `mhc.class`
with free scales, fill indexing the CDR3 ordered by CDR3 length. The two facets do not share a y
axis, so class I and class II heights are not comparable by eye.

### 9. V gene and MHC allele usage (figure)

Source `vdjdb.slim.txt`, `species == "HomoSapiens"`. `mhc.a`, `mhc.b` and `v.segm` are expanded to
one row per call, alleles are cut at the first colon and V genes at the `*`, and each cell counts the
distinct pre-expansion records.

Two filters, both on marginal totals rather than on the cell: an MHC combination is kept when it has
10 or more records summed over all V genes, and a V gene is kept when it has 10 or more summed over
all MHC combinations. Faceting is `gene` by `mhc.class` with free scales and space, so the four
panels have different extents.

Fill is `pmin(records, 1000)` on a log scale, so every cell at or above 1,000 records is the same
colour and the scale does not distinguish the largest cells from each other.

### 10. TRBV and HLA class I correspondence (figure)

The chord diagram, from the same frame as the heatmap, restricted to `gene == "TRB"`,
`mhc.class == "MHCI"` and cells with 50 or more records. Labels drop the `HLA-` prefix and rewrite
`TRBV9` as `Vb9`.

The matrix is not plotted as counts. It is divided by its row sums, then by its column sums, then
multiplied by the grand total, which makes each band the ratio of observed records to the number
expected if V gene and allele were independent at those margins. A band wider than its neighbours is
an enrichment, not a larger absolute count. Band colour runs over the range of that same matrix,
lowest as light yellow and highest as red.

### 11. Detailed summary for HLA (table)

Source `vdjdb.slim.txt`, `species == "HomoSapiens"`, `mhc.a` and `mhc.b` expanded to one row per
call and cut at the first colon, grouped by the resulting pair and sorted by records descending.
Columns: distinct pre-expansion records, distinct epitopes, Studies. Class I rows carry `B2M` as the
second chain, which is how the data spells the class I light chain rather than a curation artifact.

## Rendering it

```bash
uv run vdjdb summary --legacy out/legacy
```

The command renders `summary/vdjdb_summary.Rmd` with knitr against the legacy projection of the
current build, converts knitr's intermediate to the fragment with a single pandoc pass using
`summary/embed.tpl` and `summary/embed.lua`, then checks the result. R runs once; pandoc runs twice.

Output is `summary/vdjdb_summary_embed.html`.

| Option | Default | Effect |
|---|---|---|
| `--legacy <dir>` | `out/legacy` | the legacy projection to read |
| `--reference <file>` | none | a previous fragment, enabling the SSIM comparison |
| `--assets <dir>` | none | write figures as files instead of inlining them |
| `--verbose` | off | show knitr's chunk-by-chunk progress |

Requirements: R 4.5.3 with the packages listed in `summary/install-packages.R`, pandoc 3.10, and the
`summary` Python extra (`uv sync --extra summary`).

The render makes no network call. Publication years come from `summary/reference_years.tsv`, which
`vdjdb refs` rebuilds as a separate, committed input. If any reference in the database has no year
in that table the render stops rather than plotting an incomplete year axis.

## Fragment requirements

`vdjdb-web` locates content in the fragment by pattern, so three properties are required. Breaking
any of them produces an empty `/overview` with no error. `summary/check_summary.py` asserts all
three on every render.

1. Each image's base64 payload is a single unbroken line. The Scala side matches it with a regular
   expression that does not cross newlines.
2. There are no `<div>` wrappers.
3. Tables carry Semantic UI's classes. Pandoc does not emit them; the Lua filter adds them.

The filter works on pandoc's document tree rather than on rendered HTML text: it keeps the blocks
between the two markers in the Rmd, drops the R source, unwraps printed output, and adds the table
classes and the responsive image style. The pass does not use `--section-divs`, so no `<div>`
wrappers are produced in the first place.

## Checking a render

`summary/check_summary.py` runs three layers.

| Layer | Requires | Detects |
|---|---|---|
| structural | nothing beyond the fragment | a heading, table or image appearing or vanishing; a `dpi`, `fig.retina` or `fig.width` change, read from each PNG's IHDR |
| style | Pillow | a swapped palette, checked against the declared ColorBrewer anchors with an anti-assertion against viridis; a blank or collapsed facet, from ink fraction |
| perceptual | a previous fragment via `--reference` | a panel that has become a solid block, by SSIM per panel |

SSIM thresholds are asymmetric because the database grows: below 0.55 fails, below 0.80 warns.
They test for gross corruption, not pixel equality.

## Figure encoding

Figures are inlined as base64 by default, because that is how `vdjdb-web` finds them. Measured on
the current dashboard:

| Form | HTML per page load | Figures |
|---|---:|---|
| inlined (default) | 5.14 MB | inside the HTML, re-sent on every load |
| external (`--assets`) | 94.8 KB | 8 PNGs totalling 3.78 MB, cacheable separately |

Base64 accounts for about 1.3 MB of the inlined size as encoding overhead.

`--assets <dir>` produces the external form. It is not the default and cannot become one until
`vdjdb-web` changes: its Scala side locates images by matching `data:image/png;base64`, and a
fragment with external `src` attributes renders eight broken images.

## Previewing locally

Opening the fragment in a browser shows unstyled tables, because the Semantic UI classes come from
`vdjdb-web`'s stylesheet rather than from the fragment. `summary/preview/index.html` loads that
stylesheet from a CDN and fetches the fragment into a container sized to `vdjdb-web`'s content
column:

```bash
uv run python -m http.server -d summary 8000   # then open /preview/
```

The page below uses the same harness, with the fragment from the last successful build on `master`.

```{raw} html
<iframe src="_static/dashboard.html" title="VDJdb summary dashboard"
        style="width:100%;height:70vh;border:1px solid #d4d4d5;border-radius:4px"></iframe>
```

## Files

| Path | Role |
|---|---|
| `summary/vdjdb_summary.Rmd` | the release dashboard |
| `summary/vdjdb_paper_figures.Rmd` | figures for the manuscript, not part of a release |
| `summary/embed.tpl`, `summary/embed.lua` | the pandoc template and filter that produce the fragment |
| `summary/panels.py` | the Python port of the portable panels |
| `summary/check_summary.py` | the three-layer check |
| `summary/fingerprint.json` | committed structural fingerprint |
| `summary/reference_years.tsv` | publication years, rebuilt by `vdjdb refs` |
| `summary/annotations.tsv` | callout labels for the by-year panels |
| `summary/install-packages.R` | the R package list and its snapshot date |
| `summary/preview/index.html` | the local preview harness |
