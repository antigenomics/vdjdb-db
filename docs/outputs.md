# VDJdb build outputs — specification

Every file the build produces, what it contains, whether it ships, and who consumes it.

Becomes `docs/standards/database-outputs.rst` when the Sphinx site lands (ROADMAP phase 13); kept as
Markdown until then so it is useful now.

Status key: **shipped** in the release zip · **artifact** produced and uploaded by CI but not zipped
· **internal** produced during a build, not published.

---

## 1. Release bundles

Three zips per release. `manifest.json` names them by role so consumers select by role rather than by
asset order — see ROADMAP §3.1 for why that matters.

| Asset | Role | Contains |
|---|---|---|
| `vdjdb-<version>.zip` | `primary` | the new VDJdb format (§3) |
| `vdjdb-legacy-<version>.zip` | `legacy` | byte-layout-compatible with the historical release (§2) |
| `vdjdb-airr-<version>.zip` | `airr` | AIRR Rearrangement + Receptor/Reactivity (§4) |

Alongside them, as release assets rather than zip members: `manifest.json`, `SHA256SUMS`.

---

## 2. Legacy bundle — `vdjdb-legacy-<version>.zip`

Exactly ten members under a single `vdjdb-<version>/` directory. **The member basenames are a
contract**: `vdjmatch` looks inside the zip for `vdjdb.txt`, `vdjdb.slim.txt` and `vdjdb_full.txt` by
basename, and `vdjdb-web` resolves the rest as `<database.path>/<name>`.

| File | Rows (2026-06-03) | Cols | Consumer |
|---|---|---|---|
| `vdjdb.txt` | 284,546 | 22 | **vdjdb-web** — the schema is built from `vdjdb.meta.txt`; load-bearing |
| `vdjdb.meta.txt` | 21 + header | 8 | **vdjdb-web** — must match `vdjdb.txt` column-for-column, in order |
| `vdjdb.slim.txt` | 197,729 | 17 | vdjmatch, standalone R/Python users |
| `vdjdb.slim.meta.txt` | 16 + header | 2 | standalone users |
| `vdjdb_full.txt` | 192,753 | 35 | vdjmatch, standalone users |
| `cluster_members.txt` | 55,636 | 19 | **vdjdb-web** — parsed *positionally*, no header check |
| `motif_pwms.txt` | 40,061 | 27 | **vdjdb-web** — parsed *positionally*, no header check |
| `vdjdb_summary_embed.html` | — | — | **vdjdb-web** `/overview`, injected as an HTML fragment |
| `LICENSE` | — | — | — |
| `latest-version.txt` | 40 lines | 1 | legacy self-update clients; line 1 must point at **this** zip |

Optionally also `cluster_members_tcremp.txt` and `motif_pwms_tcremp.txt`, same positional schemas.

### Column orders (positional contracts)

`vdjdb.txt` (22): `complex.id gene cdr3 v.segm j.segm species mhc.a mhc.b mhc.class antigen.epitope
antigen.gene antigen.species reference.id vdjdb.score TCR_hash method meta cdr3fix web.method
web.method.seq web.cdr3fix.nc web.cdr3fix.unmp`

Production `vdjdb-web` additionally serves five appended `evidence.*` columns — see §3.4.

`vdjdb.slim.txt` (17): `gene cdr3 species antigen.epitope antigen.gene antigen.species complex.id
v.segm j.segm mhc.a mhc.b mhc.class reference.id vdjdb.score TCR_hash j.start v.end`

`vdjdb_full.txt` (35): the 31 chunk columns + `cdr3fix.alpha cdr3fix.beta vdjdb.score TCR_hash`

`cluster_members.txt` (19): `species antigen.epitope antigen.gene antigen.species mhc.a mhc.b
mhc.class gene cdr3aa x y cid csz v.segm j.segm v.end j.start v.segm.repr j.segm.repr`

`motif_pwms.txt` (27): `species antigen.epitope gene aa pos len v.segm.repr j.segm.repr cid csz count
count.bg total.bg count.bg.i total.bg.i need.impute freq freq.bg I I.norm height.I height.I.norm
antigen.gene antigen.species mhc.a mhc.b mhc.class`

### Legacy quirks preserved deliberately

- `vdjdb_full.txt`'s `cdr3fix.*` cells are **Python `dict` repr**, not JSON (single quotes). All
  122,930 non-empty cells in 2026-06-03 are this form.
- `vdjdb.txt`'s `method` / `meta` / `cdr3fix` use Python `json.dumps` defaults — `", "` separators
  and `ensure_ascii=True`, so `M158–66` is written `M158–66`. `struct.json_encode()` is *not*
  byte-compatible.
- No field is ever quoted; all three tables contain zero `"` characters.
- The `web.method` row of `vdjdb.meta.txt` has a space where a tab belongs, and all four `web.*` rows
  carry a field shift. Reproduced in legacy, fixed in the new format.

---

## 3. New VDJdb format — `vdjdb-<version>.zip`

A normalised star schema rather than one denormalised table. Three fact tables plus the joined view.
Parquet, with a TSV projection of each for users without a parquet reader.

### 3.1 `records.parquet` — one row per submitted record

Primary key `record_id`, unique. **One chunk row is one record**: a chunk is one paper, a row is its
report on one clone, and that row reports both chains. So this table has exactly as many rows as the
build reads — 192,753 — and `method.*` / `meta.*` sit here because the README defines them as what
the *publication* reports about the record.

The receptor is **not** here: a chain is an observation, so it is a row of `chains`. That is the
whole difference from `vdjdb_full.txt`, which folds both chains into paired columns and leaves half of
them blank.

33 columns:

| Group | Columns |
|---|---|
| identity | `record_id` |
| antigen | `species`, `mhc.a`, `mhc.b`, `mhc.class`, `antigen.epitope`, `antigen.gene`, `antigen.species` |
| provenance | `reference.id` |
| sample | `meta.study.id`, `.cell.subset`, `.subject.cohort`, `.subject.id`, `.replica.id`, `.clone.id`, `.tissue` — the id fields that are part of identity |
| annotation | `meta.epitope.id`, `.donor.MHC`, `.donor.MHC.method`, `.structure.id`, `.subset.frequency` |
| method | `method.identification`, `.frequency`, `.singlecell`, `.sequencing`, `.verification`, `.pairing` |
| score | `vdjdb.score` |
| curation | `chunk.file`, `chunk.row`, `chunk.id`, `submitter`, `comment` |

`submitter`, `comment`, `chunk.id`, `meta.subset.frequency` and `method.pairing` are **carried**, not
dropped — the legacy build discards all five. They are what makes a curation problem debuggable.

`content_hash`, the record state and the release/commit provenance live in the registry
(§6), not here: they describe the record's history rather than the record, and duplicating them into
every release would make the table a changelog.

### 3.2 `chains.parquet` — one row per TCR chain of a record

Primary key `(record_id, gene)`. This is the level `vdjdb.txt` is written at. Splitting it out is
what keeps the schema non-redundant: folding chains into records forces either duplicated record
fields (as `vdjdb.txt` does) or paired alpha/beta columns (as `vdjdb_full.txt` does).

27 columns: `record_id`, `gene` (`TRA`/`TRB`), `clonotype_id`, `cdr3`, `v.segm`, `d.segm`, `j.segm`,
`v.end`, `j.start`, `cdr3nt`, `cdr3nt.pgen`, `cdr3nt.margin`, `v.inferred`, `j.inferred`,
`d.inferred`, `d.start`, `d.end`, `d.posterior`, `d.entropy`, `cdr3.original`, `fix.needed`,
`fix.good`, `v.fix.type`, `j.fix.type`, `v.canonical`, `j.canonical`, `TCR_hash`.

**`cdr3nt` is inferred, not observed** (#461): the most plausible nucleotide junction behind the
amino-acid one, from the recombination model. 261,097 of 286,047 chains have one and every one of
them back-translates to its junction, but two models agree on only 7.2 % of the sequences, so it is a
representative history rather than evidence. `cdr3nt.pgen` is its generation probability and
`cdr3nt.margin` how far it beat the runner-up — 9.4 % are below 1.1, where the choice was
near-arbitrary. Filter on the margin rather than trusting the sequence (ROADMAP §19).

**`cdr3fix` is not a column here.** Every member of the legacy JSON blob is its own variable —
`cdr3.original` is the sequence as submitted, the four `fix.*` / `*.fix.type` columns say what was
done to it, and `v.canonical` / `j.canonical` say whether the anchors are the expected ones. A JSON
column cannot be filtered, grouped or joined without parsing, and in the release the same field is a
JSON number on one row and a string on the next. `emit/legacy.py` reassembles the blob on the way
out, which is the only place it belongs.

`clonotype_id` is a seeded hash of `(species, gene, cdr3, v.segm, j.segm)` — records reporting the
same receptor chain share it, and it is the level motif evidence and the independent-study support
count attach at. A hash rather than a counter, because a counter renumbers every clonotype the moment
a chunk is added.

`d.segm` is the **curated** D call, as the publication reported it. `d.inferred`, `d.start` and
`d.end` describe the D of the recombination scenario that produced `cdr3nt`, so the coordinates index
that sequence (0-based, half-open). They agree with the curated call at gene level on 78.4 % of the
40,892 beta chains that have one.

`v.inferred` and `j.inferred` carry a model-proposed call **only where the curator named none**
(#462) — 686 of the 711 chains with no V, 298 of the 596 with no J. They never sit beside a curated
call. Read their accuracy before using them: recovering a hidden V from the junction alone works on
23.8 % of human TRB and 50.1 % of TRA, because the junction carries little V; the J side is 95–98 %
(ROADMAP §21).

`d.posterior` is the probability of the gene `d.inferred` names, and `d.entropy` how decidable the D
was at all. **A third of beta chains have a posterior below 0.6 and an entropy above 0.9** — TRBD1
and TRBD2 are short, heavily trimmed and similar, so the junction often cannot choose between them.
Filter on `d.posterior`; do not read `d.inferred` alone (ROADMAP §20).

### 3.3 `evidence.parquet` — long format, one row per piece of evidence

Primary key `(record_id, evidence_id)`. Long rather than wide because a record may carry any number
of pieces of evidence of any number of kinds, and a wide table would be mostly empty — and would grow
a column per producer.

`record_id`, `gene` (empty when the evidence is record-level), `evidence_id`, `evidence_type`,
`evidence_source`, `evidence_value`, `evidence_score`, `first_seen_release`.

| `evidence_type` | `evidence_source` | `evidence_value` | `evidence_score` |
|---|---|---|---|
| `independent_study` | — | the *other* `reference.id`s, sorted | count of distinct references |
| `motif_tcrnet` | release tag | cluster id | cluster size |
| `motif_tcremp` | release tag | cluster id | cluster size |
| `structure_native` | PDB | PDB id | — |
| `structure_model` | model set id | structure hash | model confidence |

`independent_study` is the only producer today: **53,913 rows over 48,893 records**, scores 2 to 41.
It is the same computation as the ROADMAP §11.1 tuning objective, deliberately — one implementation,
so the shipped column and the objective cannot disagree.

Structure evidence is keyed on the legacy `TCR_hash` today and moves to `record_id` when the
structure store is re-keyed.

**No held-out validation data is ever an evidence row.** See §7.

### 3.3a `epitopes.parquet` and `restriction.parquet` — the antigen catalogue

VDJdb's own list of epitopes and the MHCs that present them.

| Table | Key | Rows |
|---|---|---|
| `epitopes` | `(antigen.epitope, antigen.species)` | 2,132 |
| `restriction` | `(antigen.epitope, antigen.species, mhc.a, mhc.b)` | 2,373 |

`epitopes` carries `antigen.gene`, `epitope.length`, `mhc.class`, and the support behind it —
`records`, `chains`, `clonotypes` and `references`. 379 epitopes are reported by two or more
publications.

**The key is the epitope and the species.** A peptide is not unique to one organism: 13 epitopes are
reported under two, and none is an error. `patches/antigen_epitope_species_gene.dict` is keyed on the
peptide alone and cannot express them, which is why this is a table rather than a view over the patch.

`restriction` checks each allele against IPD-IMGT/HLA (<https://www.ebi.ac.uk/ipd/imgt/hla/>) by
prefix, since a VDJdb call is two-field and the authority stores four: `mhc.a.status` and
`mhc.b.status` are `known`, `unknown`, or `unchecked` where no authority exists (murine and macaque
names, and `B2M`). Phase 9e adds `mhcmatch` validation — whether the allele *could present* that
peptide, not only whether its name is real (ROADMAP §12, §27).

### 3.4 `vdjdb.parquet` — the joined view

`records ⋈ chains ⋈ evidence`, one row per chain, with each evidence type pivoted to a boolean. The
convenient denormalised table, derived on every build and never authored — a consumer who edits it is
editing a cache.

Six boolean columns, **all six always present**: `evidence.motif.tcrnet`, `evidence.motif.tcremp`,
`evidence.structure.native`, `evidence.structure.model`, `evidence.validation.independent`,
`evidence.validation.same.study`. Everything without a producer yet is `false` — an honest "no
evidence of this kind", not a missing column, so the view's shape does not change as phases 8–11 land.

Those names are deliberate: production `vdjdb-web` already serves five `evidence.*` columns that
nothing in this repo produced. This is where they start being produced.

### 3.5 Generated metadata

`vdjdb.schema.json` — the field registry dumped machine-readably: per column, its `vdjdb.meta.txt`
attributes, its position in every table that carries it, and its physical dtype. Every other schema
artifact (`vdjdb.meta.txt`, the legacy column orders, the AIRR mapping in phase 7, the docs tables)
is a projection of the same registry, so none of them can drift.

The dtype is **read off the written frame**, not declared, so the schema cannot claim a type the
shipped files do not have.

### 3.6 `clusters.parquet` and `motifs.parquet` — the motif tables

The new-format counterpart of `cluster_members.txt` and `motif_pwms.txt` (§2). They exist because
`evidence.parquet` records a `motif_tcrnet` / `motif_tcremp` row whose `evidence_value` is a cluster
id: **that id has to resolve to something**, and in the legacy bundle it resolves only into a
positionally-parsed text file. One row per method per build; the `method` column separates them
rather than a second pair of files.

`clusters` — one row per `(method, cid, clonotype_id)`, i.e. a cluster's membership:

| Column | Note |
|---|---|
| `method` | `tcrnet` or `tcremp` |
| `cid` | `<species-initial>.<chain-initial>.<epitope>.<n>`, TCREMP appending `L<len>` |
| `clonotype_id` | joins to `chains.parquet`; the level motif evidence attaches at |
| `species`, `gene`, `antigen.epitope` | the scope the clustering ran in |
| `csz` | cluster size, in clonotypes |
| `v.segm.repr`, `j.segm.repr` | modal allele over the cluster |
| `x`, `y` | graph layout coordinates, for `vdjdb-web` |

**`clonotype_id`, not the CDR3/V/J triple.** The legacy file repeats `cdr3aa v.segm j.segm` plus
seven annotation columns on every member row, which is how a 55,636-row file carries 19 columns of
mostly-duplicated epitope metadata. Here the annotation lives once in `records`/`epitopes` and the
membership row carries a key — so `cluster_members.txt` is a join away, and the spec's own
normalisation rule holds (§3).

`motifs` — one row per `(method, cid, pos, aa)`, the position weight matrix:

`method`, `cid`, `pos`, `aa`, `len`, `count`, `freq`, `count.bg`, `total.bg`, `count.bg.i`,
`total.bg.i`, `level.bg`, `freq.bg`, `I`, `I.norm`, `height.I`, `height.I.norm`.

A residue a cluster never shows has **no row**, rather than a zero-count one — the logo has no letter
there, and the legacy file omits it too.

`level.bg` names which background stratum supplied `count.bg`: the `(v.gene, j.gene, len)` cell when
it has support, otherwise the coarser `len`-only cell. The legacy schema carries the imputation as a
bare `need.impute` boolean, which says that a fallback happened but not what it fell back **to**;
this names it, and the legacy projection derives the boolean from it.

**Backgrounds never ship** (CLAUDE.md hard rule 5). `count.bg` / `total.bg` are derived statistics
computed against a background streamed at build time; no background row reaches any output.

---

## 4. AIRR bundle — `vdjdb-airr-<version>.zip`

| File | Level | Rows (current corpus) |
|---|---|---|
| `vdjdb.rearrangement.tsv` | one row per chain — AIRR Rearrangement | 286,047 |
| `vdjdb.receptor.tsv` | one row per paired record — AIRR Receptor | 81,003 |
| `vdjdb.reactivity.tsv` | one row per record — AIRR Reactivity | 192,753 |
| `airr.yaml` | the AIRR schema version the files conform to (2.0) | — |

The two are linked by `cell_id`, which is the `record_id`: a VDJdb record is one publication's report
on one T-cell clone, and a clone is what AIRR's `Cell` names. Both are standard AIRR fields.

Reactivity carries `ligand_type` (`MHC:peptide`), `antigen_type` (`peptide`), `antigen`,
`antigen_source_species`, `peptide_sequence_aa`, `mhc_class`, `mhc_allele_1`, `mhc_allele_2`,
`reactivity_method`, `reactivity_readout`, `reactivity_value`, `reactivity_unit` and
`reactivity_refs`.

`reactivity_readout` is `confidence` and `reactivity_value` is `vdjdb.score`, which is what the spec
asks a non-physical assay for: *"for inferred and annotated methods this should indicate a
confidence/quality level"*. `reactivity_method` is `MHC_peptide_multimer`, `native_protein` or
`annotated` — see ROADMAP §18 for the classification and its counts.

**Nucleotide fields are present and empty** until phase 8 (#461): `sequence`, `junction`, the two
alignments and the three cigars. The schema requires the column, not a value, and
`airr.validate_rearrangement` passes on the full table today.

`Receptor` carries `receptor_variable_domain_{1,2}_aa` — the **complete mature variable domain**,
non-nullable — so both are rebuilt by stitching germline V and J around the inferred nucleotide
junction and translating (ROADMAP §22). Domain 1 is the beta chain, domain 2 the alpha, as the
schema's controlled vocabularies require. `receptor_hash` is a sha256 over the two concatenated
domains and is **not** VDJdb's `TCR_hash`, which hashes CDR3s, segments, MHC and epitope and is what
the structure store keys on.

**A receptor is a two-domain object**, so only paired records appear there: 81,003 of the 93,294
paired records have both domains rebuilt. An unpaired record is not dropped — its chain is in the
Rearrangement file and it has a Reactivity row of its own, which is where AIRR puts a single
rearranged sequence. Overall the AIRR export keeps **1,141 records and 1,501 chains that the legacy
build discards** (ROADMAP §18).

Validated in CI with the `airr` package's own schema validator — for Rearrangement. `airr` 2.0.0 has
no `Receptor` or `Reactivity` validator, so the Reactivity file is gated against the field list read
from the package's own `airr-schema.yaml` instead.

### Converting an older release

`vdjdb convert airr --legacy vdjdb.txt` produces the same two files from a legacy release zip, for
users who hold one. It runs the *same* emitter — legacy `vdjdb.txt` already speaks VDJdb's column
names — so there is no second mapping to drift. It carries no `d_call` (the file has no D column) and
the 1,501 chains the legacy build dropped are simply not there.

---

## 5. Produced but not shipped

Written to `out/reports/` on every build, uploaded as CI artifacts, **excluded from every zip**.
Their fate is an open question — ROADMAP §9.

| File | What it is |
|---|---|
| `vdjdb_full_filtered.txt` | records passing the gene, allele and canonical-CDR3 masks |
| `vdjdb_full_gene_broken.txt` | records whose V/J gene name fails the IMGT check |
| `vdjdb_full_allele_broken.txt` | records whose allele number exceeds the gene's allele count |
| `vdjdb_full_cdr3aa_broken.txt` | records whose CDR3 is not biologically valid |
| `vdjdb_full_scored.txt`, `vdjdb.slim.scored.txt`, `vdjdb.scored.txt` | the three tables with a `cluster.member` column |
| `qc.tsv` | every chunk QC finding: file, row, column, rule, value |
| `chunk-lint.tsv` | text-level findings: encoding, BOM, CRLF, header shape |
| `records.diff.tsv` | added / amended / retired records against the previous release |
| `diff-report.md` | the difference ledger against a reference release (ROADMAP §5) |
| `motifs.debug/` | per-clonotype enrichment statistics, embeddings, cluster labels, eps sweeps, the pooled cross-epitope confusion matrix |
| `contact_sheet.png` | the eight dashboard panels tiled, for visual review |

The `*_scored.txt` trio is currently **wrong** — `MotifsScoresAssembler.py` sets `cluster.member` for
every record of an `(epitope, species, gene)` rather than matching on CDR3. Fixed in the rewrite.

---

## 6. The record registry — a release asset

`records.registry.tsv` maps `record_id` to its state, hashes, provenance and amendment history, and
is what makes an id survive a curator fixing a typo. It is **not committed**: at 72.7 MB for 192,753
records (19.8 MB gzipped) it would add ~20 MB to the repo per curation PR, against a `chunks/` corpus
of 42 MB. It ships as a release asset and the build fetches the previous release's copy to reconcile
against, so it is reviewed in the release diff rather than the PR diff (ROADMAP §17, phase 14).

## 6a. Committed, not produced

| File | Role |
|---|---|
| `summary/reference_years.tsv` | publication years, so the dashboard render is offline and deterministic. A committed, reviewed input refreshed by its own PR -- never written by a build |
| `summary/annotations.tsv` | dashboard event callouts, with no hardcoded coordinates |
| `rules/expected_diffs.toml` | the declared differences the ledger accepts, each with a measured row count |

**The pinned motif parameters are not a data file.** They are `TUNED` in
`vdjdb.motifs.tcrnet` and `vdjdb.motifs.tcremp`, each carrying the scorecard it was chosen from and
the one cost it pays, in the module that uses them. A TOML file would separate a number from the
measurement that justifies it, and the measurement is the part that has to survive review.



---

## 7. Never produced, never shipped

**TCRvdb / MATCHMAKERS** (Messemaker et al., doi:10.1101/2025.04.28.651095) is proprietary:
academic, non-commercial, **no redistribution, in whole or in part**.

It is used **only** as a held-out validation set, read from a path given by `VDJDB_TCRVDB` and never
from inside this repository. No file in any bundle, artifact or report may contain its rows, its
per-record labels, or any value derived from them at record granularity. Aggregate validation metrics
(counts, AUROC, recall at a threshold) may be reported; per-record verdicts may not.

Motif clustering is **tuned** on the independent-study support count — 4,129 human clonotype-epitope
pairs with ≥2 distinct `reference.id`, across 18 epitopes with ≥20 pairs. TCRvdb is touched once, at
the end, to validate. Tuning on it would invalidate it, and its 614 labels span only two epitopes
(YLQPRTFLL, GLCTLVAML) — too narrow a basis for a global hyperparameter.

CI enforces this: a guard fails the build if any TCRvdb-derived file appears in the repository or in
a release bundle.
