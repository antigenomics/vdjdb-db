# Standalone 10x processing

- [`dcode_chain_graph.py`](dcode_chain_graph.py) calls paired receptors and peptide labels from
  donor contigs and feature-barcode matrices.
- [`itrap_reprocess.py`](itrap_reprocess.py) converts deposited ITRAP count tables and invokes
  the same caller.

See [10x processing methods and commands](../docs/tenx-dextramer.md) for input formats,
formulas and output fields. The database build reads curated `chunks/`, not these raw inputs.

Keep downloaded inputs outside the tracked tree. The retired PMID 41315082 notebook used
`GSE309696_filtered_contig_annotations.csv.gz`, deposited as a series supplementary file at
[GEO GSE309696](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE309696). Retrieve it with:

```bash
mkdir -p processing/raw_data
curl --fail --location --output processing/raw_data/GSE309696_filtered_contig_annotations.csv.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE309nnn/GSE309696/suppl/GSE309696_filtered_contig_annotations.csv.gz
```

`processing/raw_data/` is ignored by Git.
