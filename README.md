# VDJDB: A curated database of T-cell receptor sequences of known antigen specificity

[![Docs](https://img.shields.io/badge/docs-docs.isalgo.dev-blue)](https://docs.isalgo.dev/vdjdb-db/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22104776.svg)](https://zenodo.org/records/22104776)
[![Build](https://github.com/antigenomics/vdjdb-db/actions/workflows/build.yml/badge.svg)](https://github.com/antigenomics/vdjdb-db/actions/workflows/build.yml)

![Splash](images/vdjdb-splash.png)

VDJdb aggregates published information on T-cell receptor antigen specificity - the ability to
recognize certain epitopes in certain MHC contexts - and curates it into a single repository.

Each observation links a receptor to a peptide-MHC complex and its experimental evidence.
See [terminology](https://docs.isalgo.dev/vdjdb-db/standards/terminology.html) for the distinction
between the epitope, MHC restriction and source organism.

Routine updates keep the database current, and a validation scheme standardizes how specificity is
reported:

* All available information on the experimental setup used to identify an epitope-specific TCR
  sequence is taken into account and reduced to a single confidence score, assigned at the database
  generation stage, which highlights the most reliable records.
* Each record is also checked automatically against a database of V/J segment germline sequences,
  which standardizes reporting of the V-J junctions and CDR3 sequences that define a T-cell clone.

This repository holds the submissions to the database and the build that validates, assembles and
publishes it. `chunks/` is the data, one file per publication; everything else is machinery.

## Documentation

Start at [VDJdb documentation](https://docs.isalgo.dev/vdjdb-db/):

- **Tutorial:** [explore your first records](docs/getting-started.md).
- **How-to guides:** [submit records](docs/submission.md), [build a release](docs/builds.md),
  [check experimental provenance](docs/submission.md#check-experimental-provenance-through-repertoire-overlap)
  and [validation evidence](docs/submission.md#classify-validation-evidence),
  or [verify build integrity](docs/build-integrity.md).
- **Reference:** [chunk format](docs/standards/chunk-format.md),
  [columns](docs/standards/columns.md), [scoring](docs/standards/confidence-score.md) and
  [output files](docs/outputs.md).
- **Explanation:** [sequence repair](docs/standards/cdr3-fixing.md),
  [motifs and denoising](docs/denoising.md), and [the dashboard](docs/dashboard.md).

Column tables, vocabularies and score rules are generated from the package declarations when the
site builds. The [summary dashboard](https://docs.isalgo.dev/vdjdb-db/dashboard.html) follows the
latest successful build on `master`; released downloads are dated snapshots.

For corpus checks, distinguish repeated source rows from assembled observations. The build removes
duplicate curation within a publication and preserves distinct experiments. See
[metadata and reference checks](docs/build-integrity.md#check-metadata-and-references).

## Using the data

Download the latest release zip from
[the releases page](https://github.com/antigenomics/vdjdb-db/releases). A web GUI is at
[vdjdb.com](https://vdjdb.com), served by [VDJdb-web](https://github.com/antigenomics/vdjdb-web).

[vdjmatch](https://github.com/antigenomics/vdjmatch) can resolve a release for you when annotating
repertoires. That path is work in progress; see the vdjmatch repository. Standalone annotation from a
downloaded release is planned.

## Building it

```bash
uv sync --extra motifs --extra summary --extra test
uv run vdjdb qc                               # chunk validation, fail-fast
uv run vdjdb build --out out/                 # the definitive tables, then every projection
uv run vdjdb make legacy --tables out/tables  # the legacy files
uv run vdjdb convert airr --tables out/tables # AIRR tables
uv run vdjdb motifs --tables out/tables       # TCRNET + TCREMP
uv run vdjdb summary --legacy out/legacy      # the dashboard, offline
uv run vdjdb identity check --tables out/tables # the identifier invariants
uv run vdjdb corpus build --tables out/tables # the reference corpus: tf-idf over 12 token families
uv run vdjdb diff <reference.zip> out/legacy  # compare against a released zip
VDJDB_REFERENCE_ZIP=reference.zip uv run pytest -q
```

Use the 2026-06-03 release as `reference.zip` for the release comparison and tests. The static
summary also requires R and pandoc; see [build requirements](docs/builds.md).

## Contributing

New records are submitted as chunks - see [the submission guide](docs/submission.md). A chunk pull
request is checked in under three minutes by `chunk-check`, against the specification documented
above.

## Citing

Please cite the most recent paper:

> Daniil V. Luppov, Anna E. Koneva, Dmitry V. Bagaev, Anastasiia V. Alexandrova, Elizaveta K.
> Vlasova, Dmitry M. Chudakov, Chihiro Motozono, Andrew K. Sewell & Mikhail Shugay. VDJdb in 2026:
> boosting T-cell receptor recognition evidence using paratope embeddings and AI-based structure
> prediction. *Nucleic Acids Research*, 2026.
> [doi:10.1093/nar/gkag904](https://doi.org/10.1093/nar/gkag904)

A release of the database itself is archived at
[doi:10.5281/zenodo.22104776](https://doi.org/10.5281/zenodo.22104776).
