# Downstream tasks across Antigenomics

Start with the requested biological output and reuse a released reference. Keep receptor,
record, cell, sample and publication identifiers through every join. A match can link one
query to several observations: never multiply clone counts by the number of matching
records. Decide whether the requested unit is cells, reads, unique clonotypes or reports.

## Tool routing

Public repository entry points checked on 2026-10-25. Read the selected version's README,
documentation and installed help before running a command; repository master can be ahead
of its published package. Install published versions in an isolated environment. Historical
study repositories are examples with pinned inputs, not interchangeable current APIs.

| Task | Repository | Input and useful output |
|---|---|---|
| Receptor/pMHC candidate annotation | [vdjmatch](https://github.com/antigenomics/vdjmatch) | One batched receptor/cohort input plus pinned VDJdb; query selections, competing labels and annotation manifest |
| Count-aware summaries, overlap, dynamics and neighbourhood enrichment | [vdjtools](https://github.com/antigenomics/vdjtools) | Canonical counted clonotype tables plus sample metadata; clone/sample summaries or neighbourhood evidence |
| Receptor embeddings and repertoire representations | [mirpy](https://github.com/antigenomics/mirpy) | Harmonized V/J/junction tables and fixed model; coordinates or sample features with model identity |
| Published prototype-embedding workflow | [tcremp](https://github.com/antigenomics/tcremp) | Study-compatible receptor inputs; embedding/clustering workflow described in the TCREMP paper |
| Read annotation or junction/domain checking | [arda](https://github.com/antigenomics/arda) | Receptor reads or supported sequence records; AIRR fields and germline/domain annotations |
| UMI preprocessing, when starting from raw reads | [migec](https://github.com/antigenomics/migec) | Protocol-specific raw reads and barcode layout; molecular consensus inputs for subsequent annotation |
| Batched exact/fuzzy sequence or proteome search | [seqtree](https://github.com/antigenomics/seqtree) | Immutable sequence index plus query batch; neighbour IDs, distances or peptide locations |
| Peptide presentation, origin and mimic context | [mhcmatch](https://github.com/antigenomics/mhcmatch) | Peptides plus MHC alleles and relevant reference; presentation/restriction, source or mimic results |
| Structural inspection and candidate ranking | [tcren](https://github.com/antigenomics/tcren) | Identity-checked TCR:pMHC coordinates and optional generator confidence; contact/geometry evidence or task-matched ranking |
| Interactive exploration | [vdjdb-web](https://github.com/antigenomics/vdjdb-web) | Web filters and record views; retain the queried snapshot and export IDs for reproducibility |

For existing receptor tables, skip raw-read processing. Low-level seqtree distance is
not a replacement for vdjmatch's pMHC ranking. Presentation prediction is a separate
question from recognition by a particular receptor.

## Annotate a repertoire and summarize support

1. Preserve species, chain, anchor-inclusive `junction_aa`, available V/J calls, clone
   counts and sample/cell IDs. Keep raw negatives and reported pairs. Do not pair bulk
   chains merely because they occur in the same sample.
2. Run one vdjmatch batch with the selected reference population and CPU allocation:

   ```bash
   vdjmatch annotate --metadata samples.tsv --vdjdb /path/to/build/tables \
     --reference-population curated-benchmark --threads 4 -o annotation
   ```

   `samples.tsv` supplies `path` and optional `sample_id`. A pinned release ZIP is also
   accepted. Keep `annotation.manifest.json` and all no-support/competing dispositions.
3. Resolve supporting VDJdb observations by chain/pMHC keys, then inspect methods,
   publications and the separate confidence forecasts as described in the usage guide.
   Retain exact-match versus neighbour status; a target selection alone has no unique
   source record ID.
4. If sample summaries were requested, use vdjtools to retain clonotype abundance and
   denominator definitions. Define how multi-target queries contribute before summing;
   retain an ambiguous category rather than double-counting their abundance. Report
   unique clonotype counts separately from their summed cell/read/frequency mass.

## Discover neighbourhoods, then ask what they resemble

Use released motif lookup first when the question is whether a receptor belongs to an
existing VDJdb motif. For a requested new repertoire-enrichment analysis, use vdjtools'
batched `overlap.tcrnet` workflow and its documented background interface. Select a
background compatible with species, chain and study design; keep its identity and the
reported enrichment quantities. Enrichment identifies an unusual neighbourhood and
does not assign a pMHC by itself.

Annotate the resulting members with vdjmatch and inspect supporting publications and
HLA context. Preserve highly expanded single clonotypes even if they have no enriched
neighbourhood. Do not pool sample-level counts as independent donor observations or
introduce new threshold scans and control analyses without a request.

For new embeddings, mirpy exposes `mir.embedding.tcremp.TCREmp`; its documented route
uses `TCREmp.from_defaults(species, chain)` and one `model.embed(frame)` batch. Use the
same fixed prototypes, gene handling and transformation for reference and query rows.
Fit a shared projection only when requested; separate PCAs cannot be overlaid as the
same coordinate system. Keep projection/model identifiers and queried record IDs.
Embedding proximity supports a neighbourhood interpretation, not a binding probability.

## Structures and candidate ranking

Retrieve coordinates through the pinned HF record links described in the main skill.
Verify junction state, both receptor chains, peptide, MHC, construct and coordinate
provenance before carrying record labels into structural analysis. Annotate domains
with arda when required by the consuming tool. Inspect deposited versus predicted
coordinates separately; do not substitute a model confidence value for an assay outcome.

Use tcren's installed CLI/documentation for contacts, geometry and a requested candidate
ranking. Current source documentation exposes `tcren screen`; check availability in the
installed release rather than assuming a repository-only command has been published.
Choose one ranking task: hold receptor fixed while varying peptide, or hold pMHC fixed
while varying receptor. A catalogue containing one cognate structure is suitable for
inspection, but does not provide a candidate-ranking cohort. New modelling is a separate
requested computation. Keep generator confidence and structural scores separate from
VDJdb assay/independent-report forecasts; none is an affinity measurement.

For peptide context, mhcmatch provides restriction, presentation, source and mimic
workflows. Preserve the observed MHC restriction alongside predicted alternatives.
Sequence resemblance or presentation of another peptide does not establish TCR
cross-reactivity; look for experiments or propose a specific follow-up.

## Evaluate or reproduce a published workflow

Only run benchmarks when requested. Declare the target population, chain/pairing scope,
reference snapshot, actual negative measurements, split unit and comparison metrics.
Separate within-known-pMHC annotation from transfer to unseen peptides, publications
or donors. Do not call untested receptor:peptide pairs experimental negatives.

For vdjmatch exact-hit exclusion, build the complete reference once and keep counts,
priors, motifs, graph structure and fitted parameters fixed. Exclude matching hits from
the returned hit list; do not puncture or refit the reference per query. Ordinary
annotation retains exact hits. Preserve source observations and report record, unique
receptor/pair and publication counts separately.

Historical examples include [annotation and motif discovery](https://github.com/antigenomics/tcr-annotation-methodology),
[precursor frequencies and rearrangement probabilities](https://github.com/antigenomics/public-epitope),
and [alpha/beta pairing preferences](https://github.com/antigenomics/trab-pairing-study).
Read their pinned data and notebook environments before reproducing them. Legacy Java/R
commands are not current Python API instructions, and reproducing an old analysis does
not silently replace its reference with today's VDJdb.

## Deliver the result

Return the requested table, motif/logo, coordinates or analysis with a small provenance
record: database and software versions, reference population, filters, matching/pairing
keys, model identifiers, raw-input origins and output paths. Report retained/dropped and
unmatched counts with denominators, and distinguish observation from inference. Keep
interpretation limits together. A newly discovered data defect belongs in a separately
authorized curation workflow; consumer analysis does not mutate VDJdb or publish data.
