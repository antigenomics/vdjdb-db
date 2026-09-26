# Getting started

VDJdb is a curated database of T-cell receptor sequences of known antigen specificity,
served at <https://vdjdb.com>. `chunks/` **is the data** -- one file per publication --
and everything else in the repository is machinery for validating it, assembling it and
publishing it. The release is the product.

## Getting the data

Download the latest release zip from
[the releases page](https://github.com/antigenomics/vdjdb-db/releases). It carries
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

Please cite the database using the **most recent** paper ``Mikhail Goncharov, Dmitry Bagaev, Dmitrii Shcherbinin, Ivan Zvyagin, Dmitry Bolotin, Paul G. Thomas, Anastasia A. Minervina, Mikhail V. Pogorelyy, Kristin Ladell, James E. McLaren, David A. Price, Thi H. O. Nguyen, Louise C. Rowntree, E. Bridie Clemens, Katherine Kedzierska, Garry Dolton, Cristina Rafael Rius, Andrew Sewell, Jerome Samir, Fabio Luciani, Ksenia V. Zornikova, Alexandra A. Khmelevskaya, Saveliy A. Sheetikov, Grigory A. Efimov, Dmitry Chudakov & Mikhail Shugay. VDJdb in the pandemic era: a compendium of T cell receptors specific for SARS-CoV-2. Nature Methods 2022.`` [doi:10.1038/s41592-022-01578-0](https://doi.org/10.1038/s41592-022-01578-0).

