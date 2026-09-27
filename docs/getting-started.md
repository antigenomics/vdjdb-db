# Getting started

VDJdb is a curated database of T-cell receptor sequences of known antigen specificity, served at
<https://vdjdb.com>. `chunks/` is the data, one file per publication; everything else in the
repository validates it, assembles it and publishes it as a release.

## Getting the data

Download the latest release zip from
[the releases page](https://github.com/antigenomics/vdjdb-db/releases). It contains
`vdjdb.txt` (the full table), `vdjdb.slim.txt` (one row per CDR3-antigen pair, easy to
parse with R or pandas), `vdjdb_full.txt` (paired-chain records), the two motif tables and
their metadata files.

`vdjmatch` resolves and downloads it for you; `vdjdb-web` serves it at vdjdb.com.

## Building it yourself

```bash
uv sync
uv run vdjdb qc                               # chunk validation, fail-fast
uv run vdjdb build --out out/                 # the definitive tables, then every projection
uv run vdjdb make legacy --tables out/tables  # the legacy files
uv run vdjdb convert airr --tables out/tables # AIRR Rearrangement + Reactivity
uv run vdjdb motifs --tables out/tables       # TCRNET and TCREMP
uv run vdjdb summary --legacy out/legacy      # the dashboard, offline
```

## Citing

Cite the most recent paper: Daniil V. Luppov, Anna E. Koneva, Dmitry V. Bagaev, Anastasiia V. Alexandrova, Elizaveta K. Vlasova, Dmitry M. Chudakov, Chihiro Motozono, Andrew K. Sewell & Mikhail Shugay. VDJdb in 2026: boosting T-cell receptor recognition evidence using paratope embeddings and AI-based structure prediction. *Nucleic Acids Research*, 2026. [doi:10.1093/nar/gkag904](https://doi.org/10.1093/nar/gkag904)
