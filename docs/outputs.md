# VDJdb build outputs

Every file the build produces, what it contains, whether it ships, and who consumes it.

`myst-parser` renders this file on the documentation site from the Markdown source, which is why it
stays Markdown rather than becoming `.rst`. Its section numbers are cited from `ROADMAP.md`, from the
package docstrings and from `docs/tuning/`, so treat them as stable identifiers and do not renumber a
section.

Status key: **shipped** in the release zip · **artifact** produced and uploaded by CI but not zipped
· **internal** produced during a build, not published.

---

## 1. Release bundles

Three zips per release. `manifest.json` names them by role, so consumers select by role rather than by
asset order (ROADMAP §3.1).

| Asset | Role | Contains |
|---|---|---|
| `vdjdb-<version>.zip` | `primary` | the new VDJdb format (§3) |
| `vdjdb-legacy-<version>.zip` | `legacy` | byte-layout-compatible with the historical release (§2) |
| `vdjdb-airr-<version>.zip` | `airr` | AIRR Rearrangement + Receptor/Reactivity (§4) |

Alongside them, as release assets rather than zip members: `manifest.json`, `SHA256SUMS`.

---

## 2. Legacy bundle - `vdjdb-legacy-<version>.zip`

Exactly ten members under a single `vdjdb-<version>/` directory. The member basenames are a contract:
`vdjmatch` looks inside the zip for `vdjdb.txt`, `vdjdb.slim.txt` and `vdjdb_full.txt` by basename,
and `vdjdb-web` resolves the rest as `<database.path>/<name>`.

| File | Rows (2026-06-03) | Cols | Consumer |
|---|---|---|---|
| `vdjdb.txt` | 284,546 | 22 | vdjdb-web; its schema is built from `vdjdb.meta.txt` |
| `vdjdb.meta.txt` | 21 + header | 8 | vdjdb-web; must match `vdjdb.txt` column-for-column, in order |
| `vdjdb.slim.txt` | 197,729 | 17 | vdjmatch, standalone R/Python users |
| `vdjdb.slim.meta.txt` | 16 + header | 2 | standalone users |
| `vdjdb_full.txt` | 192,753 | 35 | vdjmatch, standalone users |
| `cluster_members.txt` | 55,636 | 19 | vdjdb-web; parsed positionally, no header check |
| `motif_pwms.txt` | 40,061 | 27 | vdjdb-web; parsed positionally, no header check |
| `vdjdb_summary_embed.html` | - | - | vdjdb-web `/overview`, injected as an HTML fragment |
| `LICENSE` | - | - | - |
| `latest-version.txt` | 40 lines | 1 | legacy self-update clients; line 1 must point at this zip |

Optionally also `cluster_members_tcremp.txt` and `motif_pwms_tcremp.txt`, same positional schemas.

### Column orders (positional contracts)

`vdjdb.txt` (22): `complex.id gene cdr3 v.segm j.segm species mhc.a mhc.b mhc.class antigen.epitope
antigen.gene antigen.species reference.id vdjdb.score TCR_hash method meta cdr3fix web.method
web.method.seq web.cdr3fix.nc web.cdr3fix.unmp`

Production `vdjdb-web` additionally serves five appended `evidence.*` columns (§3.4).

`vdjdb.slim.txt` (17): `gene cdr3 species antigen.epitope antigen.gene antigen.species complex.id
v.segm j.segm mhc.a mhc.b mhc.class reference.id vdjdb.score TCR_hash j.start v.end`

`vdjdb_full.txt` (35): the 31 chunk columns + `cdr3fix.alpha cdr3fix.beta vdjdb.score TCR_hash`

`cluster_members.txt` (19): `species antigen.epitope antigen.gene antigen.species mhc.a mhc.b
mhc.class gene cdr3aa x y cid csz v.segm j.segm v.end j.start v.segm.repr j.segm.repr`

`motif_pwms.txt` (27): `species antigen.epitope gene aa pos len v.segm.repr j.segm.repr cid csz count
count.bg total.bg count.bg.i total.bg.i need.impute freq freq.bg I I.norm height.I height.I.norm
antigen.gene antigen.species mhc.a mhc.b mhc.class`

`vdjdb-web` parses `cluster_members.txt` and `motif_pwms.txt` positionally, against a fixed 19-entry
and 27-entry column-type array and with no header check. An inserted, removed or reordered column
mistypes or shifts the table without an error, so the column order of both files is a requirement,
not a convention.

### Preserved legacy quirks

- `vdjdb_full.txt`'s `cdr3fix.*` cells are Python `dict` repr, not JSON (single quotes). All 122,930
  non-empty cells in 2026-06-03 are this form.
- `vdjdb.txt`'s `method` / `meta` / `cdr3fix` use Python `json.dumps` defaults: `", "` separators and
  `ensure_ascii=True`, so the en dash in `M158–66` is escaped: it is written as `\u2013`.
  `struct.json_encode()` is not byte-compatible.
- No field is ever quoted; all three tables contain zero `"` characters.
- The `web.method` row of `vdjdb.meta.txt` has a space where a tab belongs, and all four `web.*` rows
  have a field shift. Reproduced in legacy, fixed in the new format.

---

## 3. New VDJdb format - `vdjdb-<version>.zip`

A normalised star schema rather than one denormalised table. Three fact tables plus the joined view.
Parquet, with a TSV projection of each for users without a parquet reader.

### 3.1 `records.parquet` - one row per submitted record

Primary key `record_id`, unique. One chunk row is one record: a chunk is one paper, a row is that
paper's report on one clone, and the row reports both chains. The table therefore has exactly as many
rows as the build reads, 192,753. `method.*` and `meta.*` sit here because they describe what the
publication reports about the record.

The receptor is not here: a chain is an observation, so it is a row of `chains`, while
`vdjdb_full.txt` folds both chains into paired columns and leaves half of them blank.

35 columns:

| Group | Columns |
|---|---|
| identity | `record_id`, `pmhc_id`, `epitope_id` |
| antigen | `species`, `mhc.a`, `mhc.b`, `mhc.class`, `antigen.epitope`, `antigen.gene`, `antigen.species` |
| provenance | `reference.id` |
| sample | `meta.study.id`, `.cell.subset`, `.subject.cohort`, `.subject.id`, `.replica.id`, `.clone.id`, `.tissue` - the id fields that are part of identity |
| annotation | `meta.epitope.id`, `.donor.MHC`, `.donor.MHC.method`, `.structure.id`, `.subset.frequency` |
| method | `method.identification`, `.frequency`, `.singlecell`, `.sequencing`, `.verification`, `.pairing` |
| score | `vdjdb.score` |
| curation | `chunk.file`, `chunk.row`, `chunk.id`, `submitter`, `comment` |

`submitter`, `comment`, `chunk.id`, `meta.subset.frequency` and `method.pairing` are kept here; the
legacy build discards all five.

`content_hash`, the record state and the release/commit provenance are in the registry (§6), not
here: they describe the record's history rather than the record.

### 3.2 `chains.parquet` - one row per TCR chain of a record

Primary key `(record_id, gene)`. This is the level `vdjdb.txt` is written at. Chains are a separate
table so that record fields are not duplicated per chain, as in `vdjdb.txt`, and not folded into
paired alpha/beta columns, as in `vdjdb_full.txt`.

28 columns: `record_id`, `gene` (`TRA`/`TRB`), `clonotype_id`, `clone_id`, `cdr3`, `v.segm`,
`d.segm`, `j.segm`,
`v.end`, `j.start`, `cdr3nt`, `cdr3nt.pgen`, `cdr3nt.margin`, `v.inferred`, `j.inferred`,
`d.inferred`, `d.start`, `d.end`, `d.posterior`, `d.entropy`, `cdr3.original`, `fix.needed`,
`fix.good`, `v.fix.type`, `j.fix.type`, `v.canonical`, `j.canonical`, `TCR_hash`.

`cdr3nt` is inferred, not observed (#461): it is the most plausible nucleotide junction behind the
amino-acid one, from the recombination model. 261,097 of 286,047 chains have one and each
back-translates to its junction, but two models agree on only 7.2 % of the sequences, so the column is
a representative history rather than evidence. `cdr3nt.pgen` is its generation probability and
`cdr3nt.margin` its margin over the runner-up; 9.4 % of those chains have a margin below 1.1, where
the choice was near-arbitrary. Filter on the margin rather than treating the sequence as observed
(ROADMAP §19).

`cdr3fix` is not a column here. Every member of the legacy JSON blob is its own variable:
`cdr3.original` is the sequence as submitted, the four `fix.*` / `*.fix.type` columns say what was
done to it, and `v.canonical` / `j.canonical` say whether the anchors are the expected ones. In the
legacy release the same field is a JSON number on one row and a string on the next.
`emit/legacy.py` reassembles the blob on the way out, and nothing else may.

`clonotype_id` is `CT` plus 16 hex digits of a sha256 over `(species, gene, cdr3, v.segm, j.segm)`.
Records reporting the same receptor chain share it, and it is the level motif evidence and the
independent-study support count attach at. `clone_id` is `CX` plus 16 hex digits over the record's two
sorted `clonotype_id`s, and is the empty string on the 99,459 records reporting a single chain.

Both are hashes of their own keys rather than counters, so adding or removing a chunk cannot renumber
anything, and sha256 rather than a library hash function so that a dependency upgrade cannot renumber
the database either. `records.pmhc_id` and `records.epitope_id` are the antigen-side equivalents.
[Identifiers](standards/identity.md) is the reference: the five levels, the two mechanisms, the
lifecycle fields and the seven invariants `vdjdb identity check` asserts.

`d.segm` is the curated D call, as the publication reported it. `d.inferred`, `d.start` and `d.end`
describe the D of the recombination scenario that produced `cdr3nt`, so the coordinates index that
sequence (0-based, half-open). They agree with the curated call at gene level on 78.4 % of the 40,892
beta chains that have one.

`v.inferred` and `j.inferred` hold a model-proposed call only where the curator named none (#462):
686 of the 711 chains with no V, 298 of the 596 with no J. They never sit beside a curated call. Read
their accuracy before using them: recovering a hidden V from the junction alone works on 23.8 % of
human TRB and 50.1 % of TRA, because the junction contains little V sequence; the J side is 95–98 %
(ROADMAP §21).

`d.posterior` is the probability of the gene `d.inferred` names, and `d.entropy` how decidable the D
was at all. A third of beta chains have a posterior below 0.6 and an entropy above 0.9, because TRBD1
and TRBD2 are short, heavily trimmed and similar, so the junction often cannot choose between them.
Filter on `d.posterior`; do not read `d.inferred` alone (ROADMAP §20).

### 3.3 `evidence.parquet` - long format, one row per piece of evidence

Primary key `(record_id, evidence_id)`. The table is long rather than wide because a record may have
any number of pieces of evidence of any number of kinds: a wide table would be mostly empty and would
gain a column per producer.

`record_id`, `gene` (empty when the evidence is record-level), `evidence_id`, `evidence_type`,
`evidence_source`, `evidence_value`, `evidence_score`, `first_seen_release`.

| `evidence_type` | `evidence_source` | `evidence_value` | `evidence_score` |
|---|---|---|---|
| `independent_study` | - | the other `reference.id`s, sorted | count of distinct references |
| `motif_tcrnet` | release tag | cluster id | cluster size |
| `motif_tcremp` | release tag | cluster id | cluster size |
| `structure_native` | PDB | PDB id | - |
| `structure_model` | model set id | structure hash | model confidence |

`independent_study` is the only producer today: 53,913 rows over 48,893 records, scores 2 to 41. It is
the same computation as the ROADMAP §11.1 tuning objective, in one implementation, so the shipped
column and the objective cannot disagree.

Structure evidence is keyed on the legacy `TCR_hash` today and moves to `record_id` when the structure
store is re-keyed.

**No held-out validation data is ever an evidence row.** See §7.

### 3.3a `epitopes.parquet` and `restriction.parquet` - the antigen catalogue

VDJdb's own list of epitopes and the MHCs that present them.

| Table | Key | Rows |
|---|---|---|
| `epitopes` | `(antigen.epitope, antigen.species)` | 2,132 |
| `restriction` | `(antigen.epitope, antigen.species, mhc.a, mhc.b)` | 2,373 |

`epitopes` has `antigen.gene`, `epitope.length`, `mhc.class`, and the support counts `records`,
`chains`, `clonotypes` and `references`. 379 epitopes are reported by two or more publications.

The key is the epitope and the species. A peptide is not unique to one organism: 13 epitopes are
reported under two species, and none of those is an error.
`patches/antigen_epitope_species_gene.dict` is keyed on the peptide alone and cannot express them, so
this is a table rather than a view over the patch.

`restriction` checks each allele against IPD-IMGT/HLA (<https://www.ebi.ac.uk/ipd/imgt/hla/>) by
prefix, since a VDJdb call is two-field and the authority stores four, and against
`proofreading/mhc_nonhuman.tsv` for the names that database does not cover (murine `H2-`, macaque
`Mamu-`, and `B2M`). So `mhc.a.status` and `mhc.b.status` read `known`, `declared` or `unknown`.

**A build carrying an `unknown` or blank call fails**, naming the value, the column, the cell count and
the chunks that report it, so only `known` and `declared` reach a release. The fix is an entry in
`patches/mhc.dict` when the call is wrong, or a row in `proofreading/mhc_nonhuman.tsv` when it is a
species IPD-IMGT/HLA does not cover. Phase 9e adds `mhcmatch` validation: whether the allele could present that peptide,
not only whether the allele name exists (ROADMAP §12, §27).

### 3.4 `vdjdb.parquet` - the joined view

`records ⋈ chains ⋈ evidence`, one row per chain, with each evidence type pivoted to a boolean. It is
the denormalised convenience table, derived on every build and never authored by hand.

Six boolean columns, all six always present: `evidence.motif.tcrnet`, `evidence.motif.tcremp`,
`evidence.structure.native`, `evidence.structure.model`, `evidence.validation.independent`,
`evidence.validation.same.study`. An evidence type with no producer yet is `false`, meaning no
evidence of that kind, rather than a missing column, so the shape of the view does not change as
phases 8–11 land.

The names match the five `evidence.*` columns production `vdjdb-web` already serves; nothing in this
repo produced them before this view.

### 3.5 Generated metadata

`vdjdb.schema.json` is the field registry in machine-readable form: per column, its `vdjdb.meta.txt`
attributes, its position in every table that includes it, and its physical dtype. Every other schema
artifact (`vdjdb.meta.txt`, the legacy column orders, the AIRR mapping in phase 7, the docs tables) is
a projection of the same registry, so none of them can drift.

The dtype is read off the written frame rather than declared, so the schema cannot claim a type the
shipped files do not have.

### 3.6 `clusters.parquet` and `motifs.parquet` - the motif tables

The new-format counterpart of `cluster_members.txt` and `motif_pwms.txt` (§2). They resolve the
cluster ids that `evidence.parquet` stores in `evidence_value` on its `motif_tcrnet` and
`motif_tcremp` rows; in the legacy bundle those ids resolve only into a positionally-parsed text file.
Both methods go into one pair of tables, separated by the `method` column rather than by a second pair
of files.

`clusters` has one row per `(method, cid, clonotype_id)`, i.e. a cluster's membership:

| Column | Note |
|---|---|
| `method` | `tcrnet` or `tcremp` |
| `cid` | `<species-initial>.<chain-initial>.<epitope>.<n>`, TCREMP appending `L<len>` |
| `clonotype_id` | joins to `chains.parquet`; the level motif evidence attaches at |
| `species`, `gene`, `antigen.epitope` | the scope the clustering ran in |
| `csz` | cluster size, in clonotypes |
| `v.segm.repr`, `j.segm.repr` | modal allele over the cluster |
| `x`, `y` | graph layout coordinates, for `vdjdb-web` |

Membership is keyed on `clonotype_id`, not on the CDR3/V/J triple. The legacy file repeats `cdr3aa
v.segm j.segm` plus seven annotation columns on every member row, which is why a 55,636-row file has
19 columns of mostly-duplicated epitope metadata. Here the annotation appears once, in `records` and
`epitopes`, and the membership row stores a key, so `cluster_members.txt` is one join away and the
normalisation rule of §3 holds.

`motifs` has one row per `(method, cid, pos, aa)`, the position weight matrix:

`method`, `cid`, `pos`, `aa`, `len`, `count`, `freq`, `count.bg`, `total.bg`, `count.bg.i`,
`total.bg.i`, `level.bg`, `freq.bg`, `I`, `I.norm`, `height.I`, `height.I.norm`.

A residue a cluster never shows has no row, rather than a zero-count one: the logo has no letter
there, and the legacy file omits it too.

`level.bg` names which background stratum supplied `count.bg`: the `(v.gene, j.gene, len)` cell when
it has support, otherwise the coarser `len`-only cell. The legacy schema records the imputation as a
bare `need.impute` boolean, which says that a fallback happened but not what it fell back to. The
legacy projection derives that boolean from `level.bg`.

Backgrounds never ship (CLAUDE.md hard rule 5). `count.bg` and `total.bg` are derived statistics
computed against a background streamed at build time; no background row reaches any output.

The per-epitope diagnostic for both tables is `reports/motifs_per_epitope.tsv` (§5), a report rather
than a shipped table, derived entirely from `clusters` and the corpus.

---

## 4. AIRR bundle - `vdjdb-airr-<version>.zip`

| File | Level | Rows (current corpus) |
|---|---|---|
| `vdjdb.rearrangement.tsv` | one row per chain; AIRR Rearrangement | 286,047 |
| `vdjdb.receptor.tsv` | one row per paired record; AIRR Receptor | 81,003 |
| `vdjdb.reactivity.tsv` | one row per record; AIRR Reactivity | 192,753 |
| `airr.yaml` | the AIRR schema version the files conform to (2.0) | - |

The files are linked by `cell_id`, which is the `record_id`: a VDJdb record is one publication's
report on one T-cell clone, and a clone is what AIRR's `Cell` names. Both are standard AIRR fields.

Reactivity has `ligand_type` (`MHC:peptide`), `antigen_type` (`peptide`), `antigen`,
`antigen_source_species`, `peptide_sequence_aa`, `mhc_class`, `mhc_allele_1`, `mhc_allele_2`,
`reactivity_method`, `reactivity_readout`, `reactivity_value`, `reactivity_unit` and
`reactivity_refs`.

`reactivity_readout` is `confidence` and `reactivity_value` is `vdjdb.score`, which is what the AIRR
spec asks of a non-physical assay: *"for inferred and annotated methods this should indicate a
confidence/quality level"*. `reactivity_method` is `MHC_peptide_multimer`, `native_protein` or
`annotated`; ROADMAP §18 gives the classification and its counts.

Nucleotide fields are present and empty until phase 8 (#461): `sequence`, `junction`, the two
alignments and the three cigars. The schema requires the column, not a value, and
`airr.validate_rearrangement` passes on the full table today.

`Receptor` requires `receptor_variable_domain_{1,2}_aa`, the complete mature variable domain and
non-nullable, so both are rebuilt by stitching germline V and J around the inferred nucleotide
junction and translating (ROADMAP §22). Domain 1 is the beta chain, domain 2 the alpha, as the
schema's controlled vocabularies require. `receptor_hash` is a sha256 over the two concatenated
domains and is not VDJdb's `TCR_hash`, which hashes CDR3s, segments, MHC and epitope and is the key
the structure store uses.

A receptor is a two-domain object, so only paired records appear in that file: 81,003 of the 93,294
paired records have both domains rebuilt. An unpaired record is not dropped; its chain is in the
Rearrangement file and it has a Reactivity row of its own, which is where AIRR puts a single
rearranged sequence. Overall the AIRR export keeps 1,141 records and 1,501 chains that the legacy
build discards (ROADMAP §18).

CI validates the Rearrangement file with the `airr` package's own schema validator. `airr` 2.0.0 has
no `Receptor` or `Reactivity` validator, so the Reactivity file is gated against the field list read
from the package's own `airr-schema.yaml` instead.

### Converting an older release

`vdjdb convert airr --legacy vdjdb.txt` produces the same two files from a legacy release zip, for
users who hold one. It runs the same emitter, since legacy `vdjdb.txt` already uses VDJdb's column
names, so there is no second mapping to drift. The output has no `d_call`, because the legacy file has
no D column, and the 1,501 chains the legacy build dropped are absent.

---

## 5. Produced but not shipped

Written to `out/reports/` on every build, uploaded as CI artifacts, excluded from every zip. Whether
they ship in future is open (ROADMAP §9).

| File | What it is |
|---|---|
| `vdjdb_full_filtered.txt` | records passing the gene, allele and canonical-CDR3 masks |
| `vdjdb_full_gene_broken.txt` | records whose V/J gene name fails the IMGT check |
| `vdjdb_full_allele_broken.txt` | records whose allele number exceeds the gene's allele count |
| `vdjdb_full_cdr3aa_broken.txt` | records whose CDR3 is not biologically valid |
| `vdjdb_full_scored.txt`, `vdjdb.slim.scored.txt`, `vdjdb.scored.txt` | the three tables with a `cluster.member` column |
| `qc.tsv` | every chunk QC finding: file, row, column, rule, value |
| `qc-summary.tsv` | one row per rule that fired: level, rule, findings, chunks, whether it is advisory. Written by `vdjdb qc --report`. Gated against `rules/qc_advisories.tsv` with a per-rule tolerance (`tests/unit/test_qc_advisories.py`), because an advisory is advisory since only a curator can resolve it, not because its size does not matter: 10,555 within-chunk duplicates and 209 segment calls with no CDR3 were printed to a log line and recorded nowhere |
| `chunk-lint.tsv` | text-level findings: encoding, BOM, CRLF, header shape |
| `records.diff.tsv` | added / amended / retired records against the previous release |
| `diff-report.md` | the comparison against a reference release (ROADMAP §5) |
| `motifs_per_epitope.tsv` | one row per (species, gene, method, epitope): clonotypes, clustered, retention, clusters, largest cluster, mean cluster size, singleton clusters, percolation, replicated, tp, precision, lift. 678 rows. Written by `vdjdb motifs` on every run; the pooled motif scorecard averages a strongly bimodal distribution and must not be reported without this table (`docs/clustering.md` §6) |
| `motif-timings.tsv` | one row per motif stage: the same columns as `build-timings.tsv`, written by `vdjdb motifs`. **This is where the pipeline's memory peak is** - 6,898 MiB against the assemble stage's 1,577, set by the two PWM-and-emit steps, so 43 % of a 16 GB runner. Five measurements span 6,786 to 7,531 MiB over two hosts, an 11 % spread on the runner alone, and the 10,240 MiB budget is set against the largest of them. Share gated against `rules/motif_timings.tsv`, peak against a 10,240 MiB budget. Unlike the assemble report, the shares here are **not** comparable across host classes: per-stage slowdown from a 16-core laptop to the 4-vCPU runner ranges from 1.46x to 24.0x, the outlier being the background fetch, so the baseline carries the core count it was recorded on and the comparison runs on that host class only |
| `lookalikes.tsv` | one row per spelling of a value that some other value differs from only in case or in a `-`, `_`, `.` or space: the column, the folded form both share, the record count, and `same.species`. **Advisory, gated by nothing.** Measured on the current corpus: 20 groups over `antigen.gene`, `method.identification` and `reference.id`, 10 of them within one species. A value one character from another may be a typo or may be two stains, two serotypes, or one gene under two species' symbol conventions - HGNC capitalises `MBP` for human and MGI title-cases `Mbp` for mouse - so `same.species` is the column to sort on and the curator is the only one who can decide. The fold deliberately keeps `+` and `-`: applied to `meta.cell.subset` a stripping fold reports 11 groups and every one is false, because `CD95-` and `CD95+` are two populations. Edit distance was measured and rejected - at ratio 0.85 it returns real gene families (`MAGE-A1`/`A2`/`A3`/`A4`, `PPM1`/`PPM1F`) and distinct organisms (`CMV`/`MCMV`/`LCMV`) |
| `build-timings.tsv` | one row per build stage: wall seconds, share of the recorded total, peak RSS in MiB, the record count and the core count. Written by `vdjdb build` on every run. The share is gated against `rules/build_timings.tsv` (`tests/release/test_build_timings.py`); the seconds are recorded and not gated, because they are a property of the host, so a uniform slowdown is visible in the artifact rather than caught by a bar. Peak RSS is gated absolutely at 4,096 MiB for this stage, because memory is a property of the data and the code rather than of the host; measured 1,577 MiB on a laptop and 1,064 MiB on the runner, which allocates less because polars chunks to fewer threads. The share baseline is recorded on the runner and carries its core count: measured, the nine shares here agree to within 1.3 points between 4 and 16 cores, which is why a share gate works for this report. ⚠ The assemble stage is **not** the pipeline's memory peak - see `motif-timings.tsv`. `annotate.junction.add_junction_nt` is 87.2 % of the wall time (`antigenomics/vdjtools#181`); it runs as one `vdjdb infer-nt` process per slice, see `docs/builds.md` |
| `motif-metrics.tsv` | one row per (species, gene, source, axis): twelve axes for four sources -- the last legacy release, the latest release, this build's two methods, and the partition that clusters nothing. 130 rows. Written by `vdjdb motif-metrics`, which also gates them: `current-*` against `latest` catches a code regression, every source against `rules/motif_metrics.tsv` catches a corpus one. The metrics used to live only inside test assertions, so a corpus change moved them inside the slack and nobody learned the new values |
| `motif-metrics.md` | the same table as markdown, with each axis against its baseline and against `latest`, written into the CI step summary so the values that did **not** trip a gate are still read |
| `motifs.debug/` | per-clonotype enrichment statistics, embeddings, cluster labels, eps sweeps, the pooled cross-epitope confusion matrix |
| `contact_sheet.png` | the eight dashboard panels tiled, for visual review |

The three `*_scored.txt` files are currently incorrect: `MotifsScoresAssembler.py` sets
`cluster.member` for every record of an `(epitope, species, gene)` rather than matching on CDR3, which
the rewrite fixes.

---

## 6. Record registry and the identity lifecycle

`identity-lifecycle.tsv` is one row per id ever published, at any level: `id`, `level`, `state`,
`first_release`, `last_release`, `replaced_by`. No key columns, because a key is recomputable from
`chunks/` and an id's history is not. 15.9 MB for 274,683 derived ids, so it ships as a release asset
beside the zips and is listed in `SHA256SUMS`. It is written by `vdjdb release` and never by a build:
a curation branch that adds a clonotype and removes it again has retired nothing. A build with no
previous copy produces exactly the same ids and reports only that the history is unknown.

`records.registry.tsv` maps `record_id` to its state, hashes, provenance and amendment history, so an
id survives a curator fixing a typo. It is not committed: at 72.7 MB for 192,753 records (19.8 MB
gzipped) it would add ~20 MB to the repo per curation pull request, against a `chunks/` corpus of
42 MB. It ships as a release asset, and the build fetches the previous release's copy to reconcile
against, so it is reviewed in the release diff rather than the pull-request diff (ROADMAP §17,
phase 14).

## 6a. Committed, not produced

| File | Role |
|---|---|
| `summary/reference_years.tsv` | publication years, so the dashboard render is offline and deterministic. A committed, reviewed input refreshed by its own pull request, never written by a build |
| `summary/annotations.tsv` | dashboard event callouts, with no hardcoded coordinates |
| `rules/expected_diffs.toml` | the differences `vdjdb diff` accepts, each declared with a measured row count |

The pinned motif parameters are not a data file. They are the `TUNED` constants in
`vdjdb.motifs.tcrnet` and `vdjdb.motifs.tcremp`, each stated with the scorecard it was chosen from
and the one cost it pays, in the module that uses it.

---

## 7. Never produced, never shipped

**TCRvdb / MATCHMAKERS** (Messemaker et al., doi:10.1101/2025.04.28.651095) is proprietary:
academic, non-commercial, **no redistribution, in whole or in part**.

It is used only as a held-out validation set, read from a path given by `VDJDB_TCRVDB` and never from
inside this repository. No file in any bundle, artifact or report may contain its rows, its per-record
labels, or any value derived from them at record granularity. Aggregate validation metrics (counts,
AUROC, recall at a threshold) may be reported; per-record verdicts may not.

Motif clustering is tuned on the independent-study support count: 4,129 human clonotype-epitope pairs
with ≥2 distinct `reference.id`, across 18 epitopes with ≥20 pairs. TCRvdb is read once, at the end,
to validate. Tuning on it would invalidate it as a validation set, and its 614 labels span only two
epitopes (YLQPRTFLL, GLCTLVAML), too narrow a basis for a global hyperparameter.

CI enforces this: a guard fails the build if any TCRvdb-derived file appears in the repository or in a
release bundle.
