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
build reads, and `method.*` / `meta.*` sit here because the README defines them as what the
*publication* reports about the record.

| Group | Columns |
|---|---|
| identity | `record_id`, `content_hash`, `record_state` |
| complex | the 15 identifying fields: `species`, `cdr3.alpha`, `v.alpha`, `j.alpha`, `cdr3.beta`, `v.beta`, `d.beta`, `j.beta`, `mhc.a`, `mhc.b`, `mhc.class`, `antigen.epitope`, `antigen.gene`, `antigen.species`, `reference.id` |
| method | `method.identification`, `.frequency`, `.singlecell`, `.sequencing`, `.verification`, `.pairing` |
| meta | `meta.study.id`, `.cell.subset`, `.subset.frequency`, `.subject.cohort`, `.subject.id`, `.replica.id`, `.clone.id`, `.epitope.id`, `.tissue`, `.donor.MHC`, `.donor.MHC.method`, `.structure.id` |
| curation | `submitter`, `comment`, `chunk.id` |
| score | `vdjdb.score` |
| provenance | `chunk.file`, `chunk.row`, `first_seen_release`, `first_seen_commit`, `last_modified_release`, `last_modified_commit`, `amendment_count` |

`submitter`, `comment`, `chunk.id`, `meta.subset.frequency` and `method.pairing` are **carried**, not
dropped — the legacy build discards all five. They are what makes a curation problem debuggable.

### 3.2 `chains.parquet` — one row per TCR chain of a record

Primary key `(record_id, gene)`. This is the level `vdjdb.txt` is written at. Splitting it out is
what keeps the schema non-redundant: folding chains into records forces either duplicated record
fields (as `vdjdb.txt` does) or paired alpha/beta columns (as `vdjdb_full.txt` does).

`record_id`, `gene` (`TRA`/`TRB`), `cdr3`, `v.segm`, `j.segm`, `d.segm`, `v.end`, `j.start`,
`d.start`, `d.end`, `cdr3nt`, `cdr3nt.pgen`, `cdr3nt.margin`, `cdr3fix` (JSON), `TCR_hash`,
`clonotype_id`.

`clonotype_id` collapses records that describe the same receptor chain, which is the level motif
evidence attaches at.

### 3.3 `evidence.parquet` — long format, one row per piece of evidence

Primary key `(record_id, evidence_id)`. Long rather than wide because a record may carry any number
of pieces of evidence of any number of kinds, and a wide table would be mostly null.

`record_id`, `gene` (null when the evidence is record-level), `evidence_id`, `evidence_type`,
`evidence_source`, `evidence_value`, `evidence_score`, `first_seen_release`.

| `evidence_type` | `evidence_source` | `evidence_value` | `evidence_score` |
|---|---|---|---|
| `motif_tcrnet` | release tag | cluster id | cluster size |
| `motif_tcremp` | release tag | cluster id | cluster size |
| `structure_native` | PDB | PDB id | — |
| `structure_model` | model set id | structure hash | model confidence |
| `independent_study` | — | the other `reference.id` | count of distinct references |

Structure evidence is keyed on the legacy `TCR_hash` today and moves to `record_id` when the
structure store is re-keyed.

**No held-out validation data is ever an evidence row.** See §6.

### 3.4 `vdjdb.parquet` — the joined view

`records ⋈ chains ⋈ evidence`, pivoted so each evidence type becomes a boolean or count column. This
is the convenient denormalised table, derived and never authored:

`evidence.motif.tcrnet`, `evidence.motif.tcremp`, `evidence.structure.native`,
`evidence.structure.model`, `evidence.validation.independent`, `evidence.validation.same.study`.

Those last names are deliberate: production `vdjdb-web` already serves five `evidence.*` columns that
nothing in this repo produces. This is where they start being produced.

### 3.5 Generated metadata

`vdjdb.schema.json` — the field registry dumped machine-readably: per column, its dtype, which tables
carry it, its position in each, its `vdjdb.meta.txt` attributes and its AIRR mapping. Every other
schema artifact (`vdjdb.meta.txt`, the legacy column orders, the AIRR mapping, the docs tables) is a
projection of it, so none of them can drift.

---

## 4. AIRR bundle — `vdjdb-airr-<version>.zip`

| File | Level |
|---|---|
| `vdjdb.rearrangement.tsv` | one row per chain — AIRR Rearrangement |
| `vdjdb.receptor.tsv` | one row per record — AIRR Receptor + Reactivity |
| `airr.yaml` | the AIRR schema version the files conform to |

Reactivity carries `antigen`, `antigen_type`, `peptide_sequence_aa`, `mhc_class`, `mhc_gene_1`,
`mhc_allele_1`, `mhc_gene_2`, `mhc_allele_2`, `reactivity_method`, `reactivity_readout`.

Unpaired records get a Receptor row with domain 2 empty — documented, not dropped.

Validated in CI with the `airr` package's own schema validator.

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

## 6. Committed, not produced

| File | Role |
|---|---|
| `records.registry.tsv` | the record-identity registry: `record_id` → state, hashes, provenance, amendment history. Committed and reviewed as part of a curation PR |
| `summary/reference_years.tsv` | publication-year cache, so the dashboard render is offline and deterministic |
| `summary/annotations.tsv` | dashboard event callouts, with no hardcoded coordinates |
| `rules/expected_diffs.toml` | the declared differences the ledger accepts, each with a measured row count |
| `config/motifs.toml` | every pinned motif parameter |

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
