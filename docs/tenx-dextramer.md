# Process a 10x dextramer assay

A 10x Genomics dextramer experiment is one chunk. The caller is
`processing/dcode_chain_graph.py`. Call `process_donor(dcode, donor)` for one donor, or the
command below for the donors that belong to a chunk. Deposited ITRAP count tables use
`processing/itrap_reprocess.py`, which reshapes those tables and then calls the same model.
The negative-binomial mixture is not retuned there.

## Call the donor matrices

```bash
uv run python processing/dcode_chain_graph.py \
  --dcode /path/to/matrices \
  --donors 1 2 3 4 \
  --out /path/to/out
```

`--dcode` is the directory of per-donor matrices. `--donors` defaults to `1 2 3 4`. `--out` is
the directory that receives the five tables.

Each donor reads two files:

- `vdj_v1_hs_aggregated_donorN_all_contig_annotations.csv.gz`
- `vdj_v1_hs_aggregated_donorN_binarized_matrix.csv.gz`

`N` is the donor integer. The matrix file is numeric UMI counts. A column whose name ends in
`_binder` is a previous call and is not read. Surface-protein columns are not specificity evidence.

`process_donor` returns five frames: `chains`, `edges`, `clones`, `calls`, `checks`. The command
writes those frames as comma-separated tables named `chains.tsv`, `edges.tsv`, `clones.tsv`,
`calls.tsv` and `checks.tsv`. `calls.tsv` is the specificity table: one row per called label and
peptide, with the posterior of the specific hypothesis in `score`.

## What a call is

A node is one chain, named by locus, V, J and CDR3. An edge is two chains in one cell barcode
and stores both TCR contig UMI counts. A kept label is alpha-beta, alpha-alpha-beta, alpha only,
beta only, or two alphas. A second beta is not a label. An alpha-alpha-beta barcode is dropped
when its dextramer counts fit the sum of two parent clones. Contig UMI does not rescale the
dextramer UMI on a kept barcode.

Cells per label are the sample size. TCR contig UMI is the chain-graph weight. Dextramer UMI is
the multimer count. Specificity is a per-peptide negative-binomial mixture, and the call is the
maximum a posteriori hypothesis. A tie stays background. The module docstring of
`processing/dcode_chain_graph.py` is the specification of that model.

These matrices already use the AIRR junction, with both anchors included. A non-functional J
call and a non-canonical anchor are advisory findings. Leave the source sequence as reported.

## Corresponding chunks

| Chunk file | Rows | reference.id | Callable | Committed content |
|---|---:|---|---|---|
| `chunks/ISSUE_289.tsv` | 28125 | 10x application note | `process_donor` / the command above | chain-graph calls |
| `chunks/PMID_37133356.tsv` | 1440 | `PMID:37133356` | `itrap1_tables`, then `run_frames` | published ITRAP1 table |
| `chunks/DOI_10.64898_2025.12.19.695506.tsv` | 2468 | `doi:10.64898/2025.12.19.695506` | `viral_tables`, then `run_frames` | published ITRAP2 table |

The four-donor `reference.id` and `meta.study.id` are the application note
`https://www.10xgenomics.com/resources/application-notes/a-new-way-of-exploring-immunity-linking-highly-multiplexed-antigen-recognition-to-immune-repertoire-and-phenotype/#`.
The viral chunk's `meta.study.id` is `E-MTAB-13758`. Its `antigen.gene` and `antigen.species`
are blank on all 2468 rows, because the published table does not assign them.

`itrap1_tables` reads a deposited GEM table. `viral_tables` reads the deposited contig lists,
the 10x antibody matrix and the barcode-to-peptide key. Both return a contig frame and a
barcode-by-peptide count frame. `run_frames` runs one donor at a time through the chain graph
and the specificity model, and returns the same five keys. `write_result` writes those keys
under the same five filenames.
