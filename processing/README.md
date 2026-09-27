Auxiliary notebooks for handling input of different format, e.g. 10X files

## Raw inputs

`raw_data/` is untracked. Fetch what a notebook needs before running it:

| File | Origin |
|---|---|
| `GSE309696_filtered_contig_annotations.csv.gz` | GEO [GSE309696](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE309696), series supplementary |

```bash
mkdir -p raw_data && cd raw_data
curl -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE309nnn/GSE309696/suppl/GSE309696_filtered_contig_annotations.csv.gz
```

Experimental data as deposited, not derived. Nothing here is an input to the build -- the build
reads `chunks/` only.
