# The summary dashboard

The dashboard summarises the current contents of `chunks/`: records by species and chain, growth by
publication year, and breakdowns by antigen, antigen origin and MHC allele. It is served in two
places:

- `vdjdb-web` injects a fragment of it into the `/overview` page at [vdjdb.com](https://vdjdb.com).
- This page shows the same fragment below, rebuilt by CI on every push to `master`, so it
  reflects records added since the last release.

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
