# VDJDB: A curated database of T-cell receptor sequences of known antigen specificity

[![Docs](https://img.shields.io/badge/docs-docs.isalgo.dev-blue)](https://docs.isalgo.dev/vdjdb-db/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22104776.svg)](https://zenodo.org/records/22104776)
[![Build](https://github.com/antigenomics/vdjdb-db/actions/workflows/build.yml/badge.svg)](https://github.com/antigenomics/vdjdb-db/actions/workflows/build.yml)

![Splash](images/vdjdb-splash.png)

The primary goal of VDJdb is to facilitate access to existing information on T-cell receptor antigen specificities, i.e. the ability to recognize certain epitopes in certain MHC contexts.

Our mission is to both aggregate the scarce TCR specificity information available so far and to create a curated repository to store such data.

In addition to routine database updates providing the most up-to-date information, we make our best to ensure data consistency and fight irregularities in TCR specificity reporting with a complex database validation scheme:

* We take into account all available information on experimental setup used to identify antigen-specific TCR sequences and assign a single confidence score to highlight the most reliable records at the database generation stage.
* Each database record is also automatically checked against a database of V/J segment germline sequences to ensure standardized and consistent reporting of V-J junctions and CDR3 sequences that define T-cell clones.

This repository hosts the submissions to the database and the build that validates, assembles and publishes it. **`chunks/` is the data** - one file per publication - and everything else is machinery.

## Documentation

**<https://docs.isalgo.dev/vdjdb-db/>** is the full specification. Every column table, vocabulary
and score rule on that site is rendered from the build's own field registry while the page builds,
so it cannot disagree with the code.

Two parts of it answer most questions:

- **[Specification](https://docs.isalgo.dev/vdjdb-db/standards/chunk-format.html)** - what a
  submission may contain, and
  **[every shipped column](https://docs.isalgo.dev/vdjdb-db/standards/columns.html)**, table by
  table, generated from the registry.
- **[Dashboard](https://docs.isalgo.dev/vdjdb-db/dashboard.html)** - the summary panels for the
  current state of `chunks/`, rebuilt by CI on every push to `master`. It is not tied to a release,
  so records added since the last zip show up there as they land.

Readable in the tree as well:

| | |
|---|---|
| [Getting started](docs/getting-started.md) | what is in a release, and how to build one |
| [The chunk format](docs/standards/chunk-format.md) | every complex, method and meta column a submission may carry |
| [Column reference](docs/standards/columns.md) | every shipped table, generated from the registry |
| [The confidence score](docs/standards/confidence-score.md) | 0–3, and what each level asserts |
| [CDR3 fixing](docs/standards/cdr3-fixing.md) | how V/J anchors are repaired, and what `cdr3fix` records |
| [AIRR mapping](docs/standards/airr-mapping.md) | which VDJdb column is which AIRR field, and which deliberately is not |
| [Build outputs](docs/outputs.md) | every file the build produces and its contract |
| [Denoising](docs/denoising.md) | what the motif stage is for, and the rule that tunes it |
| [Clustering](docs/clustering.md) | six algorithms, 252 configurations, one harness |
| [Submitting and curating](docs/submission.md) | the submission guide and the curation skills |
| [Building and releasing](docs/builds.md) | commands, how a build is compared against the last release, reproducibility |
| [The dashboard](docs/dashboard.md) | what vdjdb-web's `/overview` is, and how it is checked |

## Using the data

Download the latest release zip from
[the releases page](https://github.com/antigenomics/vdjdb-db/releases). A web GUI is at
[vdjdb.com](https://vdjdb.com), served by [VDJdb-web](https://github.com/antigenomics/vdjdb-web).

[vdjmatch](https://github.com/antigenomics/vdjmatch) can resolve a release for you when annotating
repertoires, but that path is **work in progress**: it takes the first `.zip` asset it finds, so it
picks an arbitrary one as soon as a release carries more than one. Selecting by role from
`manifest.json` is the fix, and it ships in vdjmatch before VDJdb publishes a multi-zip release.

## Building it

```bash
uv sync
uv run vdjdb qc                               # chunk validation, fail-fast
uv run vdjdb build --out out/                 # the definitive tables, then every projection
uv run vdjdb make legacy --tables out/tables  # the legacy files
uv run vdjdb motifs --tables out/tables       # TCRNET + TCREMP
uv run vdjdb summary --legacy out/legacy      # the dashboard, offline
uv run vdjdb diff <reference.zip> out/legacy  # the difference ledger
uv run pytest -q
```

`ROADMAP.md` is the migration plan and the record of what each phase measured.

## Contributing

New records are submitted as chunks - see [the submission guide](docs/submission.md). A chunk pull
request is checked in under three minutes by `chunk-check`; the specification it is checked against
is the documentation above.

## Citing

Please cite the **most recent** paper:

> Daniil V. Luppov, Anna E. Koneva, Dmitry V. Bagaev, Anastasiia V. Alexandrova, Elizaveta K.
> Vlasova, Dmitry M. Chudakov, Chihiro Motozono, Andrew K. Sewell & Mikhail Shugay. VDJdb in 2026:
> boosting T-cell receptor recognition evidence using paratope embeddings and AI-based structure
> prediction. *Nucleic Acids Research*, 2026.
> [doi:10.1093/nar/gkag904](https://doi.org/10.1093/nar/gkag904)

A release of the database itself is archived at
[doi:10.5281/zenodo.22104776](https://doi.org/10.5281/zenodo.22104776).
