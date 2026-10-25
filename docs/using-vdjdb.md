# Query receptors, inspect evidence and retrieve structures

VDJdb connects published receptor:pMHC observations to sequence motifs, experimental
methods, record confidence and structural models. Use the released data directly for
these tasks; rebuilding the database or recomputing motifs is unnecessary.

## Choose and identify your reference

Download a dated snapshot from the [release list](https://github.com/antigenomics/vdjdb-db/releases).
Where a release provides `manifest.json`, select its `primary` bundle for linked Parquet
tables, `legacy` for historical text readers, or `airr` for AIRR tools. Check the download
against `SHA256SUMS`. Keep the release tag, checksum and tool versions with your results.
Older snapshots can lack primary tables, score exports or TCREmp motifs; choose a snapshot
containing the required files rather than treating their absence as negative evidence.

The primary and legacy bundles include motif tables. Coordinates are downloaded separately
from [Hugging Face](https://huggingface.co/datasets/isalgo/vdjdb_structure_models).
The [output reference](outputs.md) lists the exact bundle members. A current source build,
the website and a dated release can describe different snapshots.

In the Python examples below, `release/` means the extracted directory containing
`records.parquet`, not the parent directory containing the ZIP. Install Polars in a
dedicated environment, or use this repository's existing environment:

```bash
uv venv .venv
uv pip install --python .venv/bin/python polars
```

Run the snippets with `.venv/bin/python`. They read downloaded tables and do not run a build.

## Annotate a receptor or repertoire with vdjmatch

[vdjmatch](https://docs.isalgo.dev/vdjmatch/) searches reference receptors and ranks
candidate peptide-MHC assignments. One input table can contain a single receptor or an
entire repertoire. Use AIRR `junction_aa`, `v_call` and `locus`, with `j_call` where
available. VDJdb's `cdr3` is also anchor-inclusive junction sequence; AIRR `cdr3_aa`
excludes the conserved anchors and is not interchangeable with it.

Install the published annotation interface separately from the database build environment:

```bash
uv tool install 'vdjmatch==0.6.2'
vdjmatch annotate --help
vdjmatch annotate sample.tsv --vdjdb reference.zip \
  --reference-population curated-benchmark --threads 4 -o annotation
```

Here `reference.zip` is your selected primary or legacy database bundle.
`curated-benchmark` is the recommended default population in version 0.6.2; it excludes
KLGGALQAK and PMID 40498839 source records. Choose `--reference-population all-valid`
when you need the loader's full eligible population or historical reproduction. Keep the
population in the analysis record. `--threads 4` is an example budget; set it to your allocated CPU count.

Omit `--vdjdb` for automatic reference acquisition, or pass the same local ZIP to every run
for a reproducible analysis. The older `vdjmatch match` interface has different flags and
outputs; do not run these `annotate` examples through it.

To use the latest completed build in a user-specified folder, point directly to its
`tables/` directory. This reads the current files there without a network request or
falling back to a dated release:

```bash
vdjmatch annotate sample.tsv --vdjdb /path/to/build/tables \
  --reference-population curated-benchmark --threads 4 -o local-annotation
```

The Python equivalent is `vdjmatch.db.vdjdb.load("/path/to/build/tables")`.
Select one completed build, not a parent containing multiple snapshots: duplicate table
basenames are rejected as ambiguous. Retain its `BUILD.json` and `SOURCES.json` where
available, plus the annotation manifest. A local folder is not refreshed from GitHub;
omitting `--vdjdb` selects automatic released-reference acquisition instead.

For paired receptors, supply linked AIRR TRA/TRB rows with a shared `cell_id` or `pair_id`,
or wide columns `junction_aa_alpha`, `junction_aa_beta`, `v_alpha`, `v_beta` and optional
`j_alpha`, `j_beta`. Preserve the reported pairing; two chains from different cells do not
form a paired query. Use a directory or a sample manifest for a cohort:

```bash
vdjmatch annotate --metadata samples.tsv --vdjdb reference.zip \
  --reference-population curated-benchmark --threads 4 -o cohort
```

`samples.tsv` contains `path` and optional `sample_id`; paths are relative to that TSV.
Batch the inputs in one invocation. The compact outputs are `annotation.queries.tsv`,
`annotation.selection.tsv` and `annotation.manifest.json` for the first command.
Inspect competing labels, ties and no-support dispositions as well as the selected target.
`--details` adds candidates and per-target evidence, with memory proportional to queries
times targets. Ordinary annotation retains exact hits; `--exclude-exact` is a benchmark option.

For a particular target, add `--epitope GILGFVFTL --mhc-a 'HLA-A*02:01'`.
This restricts the question; it does not test competing peptides outside that scope.
See vdjmatch's [workflow guide](https://docs.isalgo.dev/vdjmatch/how-to.html) for input
formats and [score explanation](https://docs.isalgo.dev/vdjmatch/explanation.html) for
the ranking semantics.

## Inspect the records supporting a target

Query the primary tables to recover the observations, methods and publications behind
a candidate assignment. This example selects human GILGFVFTL/HLA-A*02:01 reports and
prints their reported chains:

```python
from pathlib import Path
import polars as pl

root = Path("release")
records = pl.read_parquet(root / "records.parquet")
chains = pl.read_parquet(root / "chains.parquet")
target = records.filter(
    (pl.col("species") == "HomoSapiens")
    & (pl.col("antigen_epitope") == "GILGFVFTL")
    & (pl.col("mhc_a") == "HLA-A*02:01")
    & (pl.col("mhc_b") == "B2M")
    & (pl.col("mhc_class") == "MHCI")
)
observations = target.select(
    "record_id", "pmhc_id", "reference_id", "method_identification",
    "method_verification", "vdjdb_score",
).join(chains.select("record_id", "gene", "cdr3", "v_segm", "j_segm"),
       on="record_id", how="left", validate="1:m")
print(observations.sort("record_id", "gene").head(10))
```

To locate an exact query chain, filter `chains` by `gene`, `cdr3` and its available V/J
calls, then join `records` by `record_id`. For an exact paired lookup, require both
chains on the same record ID. A vdjmatch target selection is not itself a record ID.
Keep all supporting observations rather than choosing the first publication or assigning
one source record to a similarity-based hit. Submitted versus repaired junctions and
gene calls are available separately in [the chain table](outputs.md).

## Check reliability of a specificity record

Start with the publication and its identification, sequencing, pairing and verification
methods. Then inspect the legacy score, the two separate model forecasts and the underlying
evidence. For a v3.2 snapshot, continue the preceding example:

```python
scores = pl.read_parquet(root / "record_scores.parquet")
review = target.select(
    "record_id", "reference_id", "method_identification", "method_verification",
    "method_pairing", "vdjdb_score",
).join(scores.select(
    "record_id", "validation_score", "independent_validation_score",
    "paired_lower", "paired_upper", "positive", "score_version",
    "model_sha256", "independent_model_sha256",
), on="record_id", how="left", validate="1:1")
print(review.sort("record_id").head(10))

evidence = pl.read_parquet(root / "evidence.parquet")
print(evidence.join(target.select("record_id"), on="record_id", how="semi")
      .sort("record_id", "evidence_type", "evidence_id").head(10))
```

| Readout | Question it answers |
|---|---|
| `vdjdb_score`, integer 0-3 | Which legacy experimental-evidence criteria does the report meet? |
| `validation_score`, probability 0-1 | What is the v3.2a forecast of assay validation from integrated receptor:pMHC evidence? |
| `independent_validation_score`, probability 0-1 | What is the v3.2b forecast that either chain is reported by another publication at the same pMHC? |
| `paired_lower`, `paired_upper` | How does v3.2a vary over legal completions of missing predictors? |
| `positive` | Does the entire sensitivity range pass or fail the declared decision threshold, or is the decision missing? |
| `evidence_type`, `evidence_source`, `evidence_value`, `evidence_score` | Which supporting observations and supplied metrics contributed evidence? |

Use [the score specification](standards/record-scoring-v3.md) for model identities,
chain-only forecasts and evidence encoding. Do not average the two probabilities: they
predict different outcomes. Keep `record_id` when joining scores; joining on junction
alone mixes publications and can mix pMHC assignments. Count unique `reference_id` values
when describing publications, unique `record_id` values for reports, and chain rows for chains.

## Retrieve released sequence motifs

Use `cluster_members.txt` and `motif_pwms.txt` from the primary or legacy bundle for
TCRNET. Where supplied, the `_tcremp` counterparts describe TCREmp clusters. No background
download or motif recomputation is needed to inspect these released annotations.

This example selects human beta-chain motifs for GILGFVFTL and returns both membership
and position-specific amino-acid frequencies:

```python
members = pl.read_csv(root / "cluster_members.txt", separator="\t")
pwms = pl.read_csv(root / "motif_pwms.txt", separator="\t")
scope = ((pl.col("species") == "HomoSapiens")
         & (pl.col("gene") == "TRB")
         & (pl.col("antigen.epitope") == "GILGFVFTL"))
selected_members = members.filter(scope)
keys = ["species", "gene", "antigen.epitope", "cid"]
selected_pwms = pwms.join(selected_members.select(keys).unique(), on=keys, how="semi")
print(selected_members.select("cid", "cdr3aa", "v.segm", "j.segm").head(10))
print(selected_pwms.select("cid", "len", "pos", "aa", "count", "freq",
                           "freq.bg", "height.I.norm").head(10))
```

Keep producer, species, chain, epitope and `cid` together when identifying a motif.
For a fixed-width matrix or logo, select one `cid` **and one `len`** before pivoting
`pos` by `aa`. Positions are zero-based in anchor-inclusive junction space. `freq`
is the observed residue fraction; absent amino-acid rows have zero observed frequency.
`height.I.norm` is the background-adjusted logo height. The tables retain background
counts and frequencies, without distributing background receptors.

Motifs are computed within species/chain/peptide groups, pooling MHC restrictions.
MHC annotations on the exported rows do not turn them into allele-specific motif tests.
To find motifs containing a chain, compare `cdr3aa` and the relevant V/J fields within
that group. TCRNET membership is sequence-neighbourhood evidence; TCREmp membership uses
paratope embeddings. Read [motif statistics](denoising.md) for their definitions.

## Fetch coordinates for a record

Use the structure dataset's reviewed links, then fetch only the required archive.
Choose a dataset commit from its [revision history](https://huggingface.co/datasets/isalgo/vdjdb_structure_models/commits/main)
and set `HF_REV` to its full commit hash. The same revision must be used for metadata
and coordinates. Compare `metadata/current/verification.json` with your database snapshot;
stable record IDs do not guarantee unchanged sequences between snapshots.

```bash
# Set HF_REV to the chosen dataset commit before running this command.
: "${HF_REV:?Set HF_REV to a dataset commit hash}"
uvx --from huggingface-hub hf download isalgo/vdjdb_structure_models \
  --repo-type dataset --revision "$HF_REV" --local-dir structures \
  metadata/current/record_links.tsv.gz metadata/current/asset_sources.tsv.gz \
  metadata/current/asset_inventory.tsv.gz metadata/current/verification.json
```

Resolve a record to exact archive members. Replace the example ID with your record:

```python
from pathlib import Path
import polars as pl

metadata = Path("structures/metadata/current")
record_id = "VDJDB0000188038"
links = pl.read_csv(metadata / "record_links.tsv.gz", separator="\t")
sources = pl.read_csv(metadata / "asset_sources.tsv.gz", separator="\t")
assets = pl.read_csv(metadata / "asset_inventory.tsv.gz", separator="\t")
links = links.filter(
    (pl.col("record_id") == record_id)
    & pl.col("matches_repaired") & pl.col("exact_tcr_pmhc_match")
)
files = (assets.join(sources, on="asset_row", validate="1:m")
         .join(links.select("source_row").unique(), on="source_row", how="semi")
         .filter(pl.col("archive").is_in(["pdb_files.tgz", "pdb_files_native.tgz"]))
         .select("archive", "member", "size_bytes").unique().sort("archive", "member"))
print(files)
```

`matches_repaired` selects the built junction state; use `matches_original` for submitted
junctions. `exact_tcr_pmhc_match` additionally requires agreement on recorded species,
MHC chains and class. Inspect the full link before accepting an engineered construct.
The joins are many-to-many: several metadata rows can identify one coordinate file.
Do not count link rows as structures or independent experiments.

Copy `archive` and `member` from a selected output row into shell variables. Download
that archive, then stream only the chosen member to a new destination file:

```bash
: "${archive:?Set archive from the selected row}"
: "${member:?Set member from the selected row}"
uvx --from huggingface-hub hf download isalgo/vdjdb_structure_models \
  "data/$archive" --repo-type dataset --revision "$HF_REV" --local-dir structures
(set -C; tar -xOzf "structures/data/$archive" -- "$member" > selected-structure.pdb)
```

Check archive size and available disk before downloading. A compressed tar archive must
be downloaded and scanned even when selecting one member; the command avoids extracting
all its files. The no-clobber subshell prevents replacing an existing destination.
The same inventory maps contact tables and SVG complementarity maps: select their archive
types instead of the two PDB archive names above.

For deposited mmCIF sources and additional models, consult the
[native coverage metadata](https://huggingface.co/datasets/isalgo/vdjdb_structure_models/tree/main/metadata/native_review)
and [supplemental catalogue](https://huggingface.co/datasets/isalgo/vdjdb_structure_models/tree/main/metadata/supplemental).
In the supplemental catalogue, join `record_links.tsv.gz` to `structures.tsv.gz` on
`structure_id`, preserve `sequence_state`, and use its declared archive/member paths.
Do not apply the historical `source_row` join to that separate schema.

## Interpretation limits

- vdjmatch similarity and ranking scores are not assay-validation probabilities. A
  similar receptor supports a candidate assignment, not a new experimental observation.
- Record-score sensitivity ranges are not confidence intervals. Missing evidence or a
  missing decision is not experimental non-binding. Fitted-cohort performance does not
  establish universal production accuracy; retain the score version and model digests.
- Motif absence does not invalidate a record. Motif membership is supporting group-level
  sequence evidence and does not establish a query's binding outcome or MHC restriction.
- Coordinate availability, predicted-model confidence and experimental validation are
  different facts. Retain native versus predicted provenance. Exact paired junction/pMHC
  matching does not prove identical full variable domains or experimental constructs.
- Historical archives include withdrawn assignments and excluded modified ligands. Select
  through reviewed record links; neither an archive member nor a source positive/negative
  flag supplies an assay outcome. An empty strict lookup means no matching file was found
  in that catalogue and snapshot, not that the receptor cannot bind.
