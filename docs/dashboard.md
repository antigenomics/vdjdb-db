# The summary dashboard

`vdjdb-web` injects a fragment of this dashboard into its `/overview` page. The fragment is produced
by `uv run vdjdb summary`, which renders `summary/vdjdb_summary.Rmd` against the legacy projection
of the current build, then runs **one more pandoc pass** over knitr's intermediate with
`summary/embed.html` and `summary/embed.lua` to emit the fragment directly, and verifies the result
against a committed fingerprint. The R never runs twice; only pandoc does.

That pass replaced a script that line-scanned pandoc's finished output for `<div`,
`<pre class="r">` and `<table>` — three guesses about markup, each of which silently blanks
`/overview` when it stops matching. Stating the transforms on the document tree removes the guesses:
the filter keeps the blocks between the two markers, drops the R source, unwraps printed output, and
adds the table classes and the responsive image style. The `<div>` stripping is gone entirely,
because this pass never passes `--section-divs` and so emits none.

**The render makes no network call.** Publication years come from `summary/reference_years.tsv`,
resolved once by `vdjdb refs` and committed; the document stops rather than plotting an incomplete
year axis if any reference is missing from it.

## Three properties are load-bearing

Any one of them failing silently blanks `/overview`, so `summary/check_summary.py` asserts all three
on every render:

1. each image's base64 payload is a **single unbroken line** — the Scala side matches it with a
   regex that does not cross newlines;
2. there are **no `<div>` wrappers**;
3. tables carry Semantic UI's classes, which pandoc does not emit and the filter injects.

## What the checker compares

| layer | what it needs | what it catches |
|---|---|---|
| structural | nothing beyond the fragment | a heading, table or image appearing or vanishing; a `dpi`/`fig.retina`/`fig.width` regression, from the PNG IHDR |
| style | Pillow | a swapped palette (an anti-assertion against viridis anchors), a blank or collapsed facet |
| perceptual | a previous release's fragment | "panel 7 is a solid grey block now", via SSIM per panel |

The perceptual thresholds are loose and asymmetric because the database grows: fail below 0.55, warn
below 0.80. The point is not pixel equality.

## Seeing it as the site will

Opening the fragment directly shows **unstyled** tables — the Semantic UI classes come from
vdjdb-web's own bundle, not from the fragment. `summary/preview/index.html` loads that stylesheet
from a CDN and fetches the fragment into a container sized to vdjdb-web's content column:

```bash
uv run python -m http.server -d summary 8000   # then open /preview/
```

```{raw} html
<iframe src="_static/dashboard.html" title="VDJdb summary dashboard"
        style="width:100%;height:70vh;border:1px solid #d4d4d5;border-radius:4px"></iframe>
```
