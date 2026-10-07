# The summary dashboard

The dashboard shows database counts, publication-year growth, peptide provenance, confidence
scores, junction lengths, and V-gene/MHC usage. The embedded view below comes from the latest
successful database build on `master`. Downloads on the releases page may contain an older snapshot.

## Reading the counts

Most panels read `vdjdb.slim.txt`, where repeated observations are combined into a clonotype-pMHC
row. The growth panels read `vdjdb_full.txt`, where each row is one publication's observation and
may contain both chains. These counts answer different questions.

| Label | Meaning |
|---|---|
| Records in the species/chain table | Rows in the slim projection |
| Paired records | Slim rows whose partner chain is reported |
| Studies | Every distinct reference, including all references in a combined slim row |
| Growth by year | Observations, peptides, genes and references grouped by publication year |
| Confidence score | Evidence reported for the observation; slim rows use the highest contributing score |
| Junction length | Cys104 through Phe/Trp118, including both anchors |

MHC usage panels combine allele resolutions at the first colon; V-gene panels combine alleles
at `*`. Rows with several calls contribute once to each applicable group. Group totals therefore
need not sum to the database total. Peptide-source groups describe provenance; recognition is
assigned to a peptide-MHC complex. See [terminology](standards/terminology.md).

Unknown optional metadata remains blank. Species with unsupported germline models retain their
reported sequences and calls without borrowing another species' annotation.

## Current master summary

```{raw} html
<iframe src="_static/dashboard.html" title="VDJdb summary dashboard"
        style="width:100%;height:70vh;border:1px solid #d4d4d5;border-radius:4px"></iframe>
```

## Rendering it

```bash
uv sync --extra summary
uv run vdjdb summary --legacy out/legacy
```

The checked fragment is `summary/vdjdb_summary_embed.html`. The static render requires R 4.5.3,
the packages in `summary/install-packages.R`, and pandoc 3.10. Publication years come from the
reviewed `summary/reference_years.tsv` input. Rendering is offline and fails if a reference lacks
a year. Refresh that input separately with `vdjdb refs`.

Open `summary/preview/index.html` to preview the fragment with the website's table styles.
Use `--reference <fragment>` to compare the figures with an earlier render, or `--verbose` to
show rendering progress.

## Fragment requirements

The website expects unbroken base64 image payloads, no `div` wrappers, and Semantic UI table
classes. Keep the default inline images for the website. `--assets <dir>` writes external figures
for other uses and does not produce the fragment the website expects.

`summary/check_summary.py` checks heading/table/image structure, PNG sizes, palettes and blank
panels. With `--reference`, it also checks image similarity for gross rendering corruption.
Changing database counts is expected; these checks do not require identical pixels.

## The interactive version

The interactive dashboard adds hover values and filters for growth and V-gene/MHC panels.
It uses the same calculations as the static summary.

```bash
uv run vdjdb summary --no-static --interactive --legacy out/legacy --out out/summary
```

This command needs the `summary` Python extra and no R installation. The HTML loads Plotly from
a CDN, so viewing it requires network access. CI publishes the static master summary above;
a separately generated interactive file describes the build used to produce that file.

## Files

| Path | Role |
|---|---|
| `summary/vdjdb_summary.Rmd` | Static dashboard |
| `summary/panels.py` | Shared calculations and interactive panels |
| `summary/embed.tpl`, `summary/embed.lua` | Website fragment conversion |
| `summary/check_summary.py`, `summary/fingerprint.json` | Render checks |
| `summary/reference_years.tsv`, `summary/annotations.tsv` | Reviewed year and figure-label inputs |
| `summary/install-packages.R` | R package requirements |
