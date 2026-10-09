# VDJdb build outputs

Every file the build produces, what it contains, whether it ships, and who consumes it.

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

Ten required members under a single `vdjdb-<version>/` directory. The member basenames are a contract:
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

**Column names here are `underscore_case`**, not the dotted names the legacy tables use: `mhc_a`,
`antigen_epitope`, `v_segm`, `cdr3nt_pgen`. A dot is a table qualifier in SQL and blocks attribute
access in most dataframe libraries, and these tables already mixed the two conventions -
`record_id` and `pmhc_id` beside `mhc.a`. The rule is a plain dot substitution with case left alone,
so `TCR_hash` and `meta.donor.MHC` ship as `TCR_hash` and `meta_donor_MHC`. The legacy tables do not
move: they are a positional contract `vdjdb-web` parses.

`vdjdb.schema.json` states both names per column - `name` is what the legacy tables ship and what
the registry indexes, `ships_as` is what these tables ship - and
[Columns](standards/columns.md) is the generated table for each, so neither can drift from the
build. The chunk format is unaffected: a submitted chunk still uses the dotted names.

### 3.1 `records.parquet` - one row per submitted record

Primary key `record_id`, unique. One chunk row is one record: a chunk is one paper, a row is that
paper's report on one clone, and the row reports both chains. `method.*` and `meta.*` sit here
because they describe what the publication reports about the record.

The table has one row per retained publication report. Exact within-publication duplicates
are merged; reports from different publications remain separate. When a publication appears
in two chunk files, compatible rows merge and fill each other's blanks. Different observations,
such as a solved structure and the sort that found the receptor, remain separate.
`out/reports/repeated-references.tsv` records these decisions. See the
[current master summary](dashboard.md) for database sizes.

The receptor is not here: a chain is an observation, so it is a row of `chains`, while
`vdjdb_full.txt` folds both chains into paired columns and leaves half of them blank.

The [generated column reference](standards/columns.md) lists names, types and descriptions.
Record fields cover identity, antigen, publication, sample, method, score and curation.

`submitter`, `comment`, `chunk_id`, `meta_subset_frequency` and `method_pairing` are kept here; the
legacy build discards all five.

**`method_frequency`, `method_frequency_count` and `method_frequency_total` are three independent
columns, and none is derived from another** (#696). A study that reports only a float has no count
behind it, so deriving the float would blank it exactly where it is the only measurement; deriving
the pair from a float is impossible. Measured: 43,231 records carry a count and total, 17,700 a
percentage, 2,601 a float, 129,109 nothing. The count and total are chunk columns, so a submitter
with a read count writes it as a number; where they do, the submitted value wins and nothing is
parsed. Where all three are present they must agree, which `vdjdb qc` reports and does not repair.

**`meta_subset_frequency` is not filled from `method_frequency`.** It is populated on 2,412 records
and left as submitted, because the records that do carry it are using it for a different quantity -
`method_frequency = 17/52` beside `meta_subset_frequency = 0.70%` is a clonotype's count within a
sorted subset beside that subset's share of the sample, and both are real. A consumer that wants
"the frequency of this clonotype in its subset" should coalesce the two:
`pl.coalesce("meta_subset_frequency", "method_frequency")`. Filling it in the build would have
written 61,253 cells across 117 chunks, every one a copy of the column beside it, and mixed the two
readings with nothing to tell them apart.

`content_hash`, the record state and the release/commit provenance are in the registry (§6), not
here: they describe the record's history rather than the record.

### 3.2 `chains.parquet` - one row per TCR chain of a record

Primary key `(record_id, gene)`. This is the level `vdjdb.txt` is written at. Chains are a separate
table so that record fields are not duplicated per chain, as in `vdjdb.txt`, and not folded into
paired alpha/beta columns, as in `vdjdb_full.txt`.

The [generated column reference](standards/columns.md) lists every chain field. Submitted
sequences and V/D/J calls remain available beside repaired junctions, canonical calls and
inferred segments, so consumers can distinguish the source observation from the build result.

`j_segm` is the J that ships, and it is not always the J the paper reported. Three values sit side by
side: `j_segm_submitted` is the call as submitted, `j_segm_arda` is the allele `arda.cdr3fix` aligned
the junction against, and `j_segm` is what ships. A J is **not used when it misses the junction's last
3 residues**, anchor excluded and one mismatch beside the anchor tolerated: it is replaced by the one
other gene of the chain's locus that matches 3 or more, or by the one gene of its own family that
matches 2 or more when nothing reaches 3. Ties are not broken, a call no other gene explains stands,
and the junction itself is never altered. `j_start` and `j_canonical` are then about the gene that
ships (#681; `vdjdb.curate.jcalls.recall`).

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

`v.inferred` and `j.inferred` hold a call proposed from the junction, and only where the publication
named no segment (#462, #658): 706 of the 745 chains with no V, and 3,272 chains with no submitted J.
They never sit beside a curated call.

The two do not have the same standing, and the split is deliberate. **`j.inferred` also ships**, in
`j.segm`, as the equivalent has in every release; `v.inferred` ships nowhere, and `v.segm` stays blank
on all 745 chains whose publication left it blank. The reason is measured: recovering a hidden J from
the junction alone works on 93.6–97.5 % of chains at gene level, because a J germline templates a
distinctive 3′ motif, while a V works on 23.8–50.1 %, because the junction contains little V sequence
and most of what it does contain templates the same `CAS`. So the J proposal is good enough to be the
record's J and the V proposal is not; it is offered beside the blank instead.

Two sources answer, in order: the recombination model (`vdjtools.model.infer_nt_batch`, with the blank
side marginalised), then arda's germline anchor table where the model declines - which is also the only
source for a species no model covers, and is what fills the 74 macaque `j.alpha` blanks that nothing
filled before. Until #658 this was a k-mer scan over `res/segments.txt`, a by-product of a 2023 IMGT
import whose V half never worked at all: 3 non-empty guesses in 4,000 sequences.

**`v.end.inferred` and `j.start.inferred` are the two sources of a V/J boundary, kept apart on
purpose (#631).** `v.end` and `j.start` are `arda.cdr3fix`'s answer, read off its protein alignment
against the germline of the segment the record names, and `-1` where it declined - 4,163 chains for
`v.end` and 1,164 for `j.start`. The `.inferred` pair is the boundary a second germline alignment
supports, from `vdjtools.model.germline_boundary`, which decides how far into the boundary codon the
germline reaches rather than rounding to a residue. It is filled **only where the first declined**:
2,442 and 488 chains. Everywhere else it is `-1`, so the two are distinguishable by column and
coalescing them cannot overwrite the shipped answer.

Two boundaries rather than one filled column, because they are answers to the same question with
different precision and the shipped column is what `vdjdb-web` reads. Against the external nucleotide
truth in `tests/release/test_cdr3fix_accuracy.py` the codon decision is exact on **92.9 %** of `v.end`
and **98.0 %** of `j.start`, where a protein alignment rounding to whole residues - VDJdb's k-mer
scanner and `arda.cdr3fix` alike - is 71.8 % and 97.1 %. It is not promoted into `v.end` anyway: that
is a change to what the database says about every record, which belongs to a curation decision and to
its own declared rules, not to a build improvement. What it does do is replace `-1`, which carries
nothing at all.

The `.inferred` pair used to come from the argmax recombination history behind `cdr3nt` instead. That
was 11 points worse on `v.end`, because maximising P(sequence) explains N-region nucleotides as
templated whenever it can; `antigenomics/vdjtools#182` was opened from this measurement and fixed it.
The legacy tables and the `cdr3fix` JSON `vdjdb-web` parses keep the alignment's `-1` untouched
throughout.

`d.posterior` is the probability of the gene `d.inferred` names, **from the same recombination
scenario weights that named it** - so the number beside the call is the probability of that call.
Naming the D and placing it are separate questions and one estimator answers each: the model names
the gene, the aligner places it. Measured on 4,000 real human TRB rearrangements whose D and
coordinates come from the nucleotide sequence, the gene is right on 74.35 % of all rows and 99.80 %
carry coordinates. TRBD1 and TRBD2 are short, heavily trimmed and similar, so the junction often
cannot choose between them: filter on `d.posterior`; do not read `d.inferred` alone (ROADMAP §20).

There is no `d.entropy`. It came from a separate estimator (`arda.dpost`) that named a different D
gene from the one it annotated on 21.5 % of chains, and that estimator is retired: the posterior now
comes from the scenario weights and there is no second distribution to take an entropy over.

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

| Table | Key |
|---|---|
| `epitopes` | `(antigen.epitope, antigen.species)` |
| `restriction` | `(antigen.epitope, antigen.species, mhc.a, mhc.b)` |

`epitopes` has `antigen.gene`, `epitope.length`, `mhc.class`, and the support counts `records`,
`chains`, `clonotypes` and `references`. Counts are recomputed on each build.

**`proteome_peptide` and `proteome_substitution`** link an epitope to the host-proteome peptide it is
one substitution from, where the epitope is not itself in the proteome (#632). In the
2026-09-29 measurement, 211 of the 2,131
epitopes carry the link, and **67 of those have the proteome form curated in VDJdb as a separate
row** - 35,346 records, 18.3% of the database, on rows nothing else says are two forms of one
peptide:

| epitope | records | `proteome_peptide` | records | `proteome_substitution` | `antigen_gene` |
|---|--:|---|--:|---|---|
| `SLLMWITQV` | 29,729 | `SLLMWITQC` | 13 | `9C>V` | NY-ESO-1 |
| `ELAGIGILTV` | 2,401 | `EAAGIGILTV` | 140 | `2A>L` | MLANA |
| `VEALYLVSG` | 2,495 | `VEALYLVCG` | 5,048 | `8C>S` | INS |
| `IMDQVPFSV` | 100 | `ITDQVPFSV` | 19 | `2T>M` | PMEL |

A query for NY-ESO-1 responses returns the 29,729 and never learns the 13 exist. That is what the
columns are for.

**Neither row is wrong, and this is not a defect report.** The epitope sequence is the ground truth -
it is the peptide the experiment used - and the difference from the proteome is almost always
deliberate: an anchor-optimised vaccine peptide, a designed altered-peptide ligand, a heteroclitic
variant, or a structure solved with a modified peptide. Of the 36,496 records on these 211 epitopes,
58 carry a `meta_structure_id` and 54 come from `PDB_Database.tsv`, so the crystallography case is
real and small. The columns are named for what was measured and not for why, because the sequence
tells those causes apart from none of the others.

Empty on the other 1,920 rows, and on every viral or bacterial epitope, because only the human and
mouse proteomes are read - a pathogen epitope's source is the pathogen's proteome, which is a
per-pathogen fetch and a different question. Both columns are read from
`proofreading/epitope_proteome.tsv`, a committed reviewed input refreshed by `vdjdb antigens` in its
own pull request, never resolved during a build: `mhcmatch` fetches the proteome from HuggingFace,
so doing it here would put a network call in the critical path (hard rule 9).

**The key is the epitope and the species, not the peptide.** A reader who assumes one row per peptide
- which the table's name invites - joins the records of 12 peptides twice.
`patches/antigen_epitope_species_gene.dict` is keyed on the peptide alone and cannot express them, so
this is a table rather than a view over the patch.

Those 12 are three different things, and `out/reports/epitope-sources.tsv` separates them (#633):

* **6 are conserved peptides** - the same sequence in two organisms' proteomes, so two publications are
  two independent reports. `VEALYLVCG` is in human `INS` and mouse `Ins2`, `RPIIRPATL` in influenza A
  and B `NP`, `LRVMMLAPF` in *E. coli* and *S.* Typhi `yeiH`. The report marks them `conserved`, which
  is the column to sort on, so a *new* collision is the only thing to look at.
* **6 are vocabulary gaps** - one organism written two ways (`AdV` beside `HAdV5`), or `Synthetic` in
  `antigen.species`, which is a provenance and not a species. These want an alias table and which way
  each folds is a curator's call (#632, #637).
* a mis-curation, now repaired: `RGPGRAFVTI` was `HomoSapiens` / `P18-I10` on one row against `HIV-1` /
  `GP160` on 85, where `P18-I10` is the laboratory name of the HIV-1 V3-loop peptide itself.

The same report carries the other half of the ambiguity, which nothing reported before. `epitopes` keeps
the **modal** `antigen.gene` for each `(peptide, species)`, and 30 of them have more than one label - so
the report counts them in a `genes` column beside the label the catalogue kept.

Its first run found `FVVKAYLPVNESFAFTADLRSNTGGQA` carrying **187** labels, `Eef2` to `Eef188`, one per
record and consecutive: a spreadsheet autofill had incremented one gene name down the column. That is
repaired (#397, #694) - `mhcmatch`'s mouse proteome assigns the peptide to `Eef2` and to nothing else -
so the largest remaining is 2.

`restriction` checks each allele against IPD-IMGT/HLA (<https://www.ebi.ac.uk/ipd/imgt/hla/>) by
prefix, since a VDJdb call is two-field and the authority stores four, and against
`proofreading/mhc_nonhuman.tsv` for the names that database does not cover (murine `H2-`, macaque
`Mamu-`, and `B2M`). So `mhc.a.status` and `mhc.b.status` read `known`, `unconfirmed`, `declared` or
`unknown`.

`unconfirmed` is a **refinement of `known`, not a rejection** (#634): the call resolves in IPD and every
allele under it is `Unconfirmed` - it rests on a single submission rather than an independent
observation. The existence check cannot make that distinction, because every one of these carries a name
IPD has. Measured 2026-09-29: 14 calls over 105 records, and 80 of those are one call,
`HLA-A*02:01:48` - a third-field allele resting on one cell from one submitting group, where
`HLA-A*02:01` has 169 Confirmed alleles under it. What to do about that is a curation question about
what the submitters meant, so it is reported and never fatal.

**An unknown call fails the build**, naming the value, column and affected chunks.
A reported class-II chain can have an unreported partner; the build preserves that partial
restriction. See [chunk format](standards/chunk-format.md) for permitted missing values. The fix is an entry in
`patches/mhc.dict` when the call is wrong, or a row in `proofreading/mhc_nonhuman.tsv` when it is a
species IPD-IMGT/HLA does not cover. Whether the allele could *present* that peptide, rather than only
whether its name exists, is phase 9e: `presentation.tsv` above for the offline half, and the six
columns below for the model's ranking.

`restriction` also carries six promiscuity columns (ROADMAP §10.6). An epitope is often presented by
several alleles and the curated one is not always the best binder; both facts belong in the database
and neither belongs in a key, because a model upgrade would otherwise renumber `pmhc_id` and break
every external reference while the build passed.

| Column | What it is |
|---|---|
| `alleles.reported` | distinct `mhc.a` values VDJdb records for this epitope. **Curation, not prediction**, so it is counted from this table and is present on a class II row where the other four are blank. Two or more means different publications restricted the same peptide differently - the question #372 asks |
| `mhc.a.top` | the panel allele that presents this epitope best |
| `mhc.a.rank` | where the curated allele sits in that ranking, 1 being the top |
| `mhc.a.percentile` | the curated allele's `%Rank_EL` against the human proteome background |
| `promiscuity` | panel alleles in the strong band for this epitope |
| `mhcmatch.version` | the model that produced the five above, per row. Empty on a row scored before the column existed |

All six are a **join against the committed `proofreading/epitope_promiscuity.tsv`** (13,510 rows over
1,729 epitopes and 107 alleles). Nothing is predicted during a build: `mhcmatch` fetches its reference
data from HuggingFace, and a build that downloads a model is neither offline nor deterministic (hard
rule 9), so the table is a reviewed input refreshed by `vdjdb promiscuity` through its own pull
request.

**A curated allele the prediction outranks is not a defect.** The curated allele is the one a
publication typed a donor for, which a proteome-background ranking has no access to, and most of these
epitopes are promiscuous. 1,793 of 1,989 class I pairs get a rank; the 196 that do not divide into
four causes, each a different statement: 86 pairs (16,316 records) name an allele outside the panel,
86 (10,395) an epitope the table has not been refreshed to cover since it was written, 24 (381) an
epitope outside the 8-11mer class I range, and the rest resolve only at a depth the panel does not
name. A deeper spelling such as `HLA-A*02:01:48` is scored at its two-field molecule, because the
panel is named at two fields and there is no deeper groove.

### 3.3b `epitope_assessment.parquet` - reported peptides and predicted binding cores

Also shipped as `epitope_assessment.tsv` in the primary bundle. This extends the presentation
annotation in section 3.3a and ROADMAP section 10.6 (#1315). It records predictions separately from
the `epitopes` provenance catalogue and the reported MHC pairs in `restriction`.

The key is `(antigen_epitope, antigen_species, antigen_gene, species, mhc_species, mhc_class,
mhc_a, mhc_b, prediction_allele)`. Every reported pair has a row with `reported = true`, including
unsupported peptides and molecules. Competing parent-gene labels are retained. `species` describes
the receptor; `mhc_species` selects the presentation model from the reported molecule; neither is
the peptide's `antigen_species`. Additional predicted pairings have `reported = false`, blank
`mhc_a`/`mhc_b` and zero `records`/`references`. Support counts belong to reported pairs only.

`prediction_allele` is the mhcmatch panel key. Class-II keys identify a molecule, including both
polymorphic chains for DP/DQ; an absent DP/DQ partner is not imputed. `allele_resolution` distinguishes
an exact resolution from prefix completion. The lowest presentation percentile is marked
`prediction_best`, with ties broken by allele name. All weak/strong predicted presenters and every
reported pairing are retained. The best panel allele is retained even when its band is non-binder.
That flag is a ranking within this panel, not evidence that the peptide is presented.

Class-II percentile ranks use a random-peptide background of the same length as the scored
peptide, so maximizing over additional registers does not inflate the presentation band.
Class I retains mhcmatch's marginal background and length preference. `presentation_p_present`
is the published scorer's separate isotonic probability; it is not a length-conditioned percentile.

`prediction_peptide` is the scored sequence. A reported class-I peptide longer than eleven residues
is assessed as binding-length windows, retaining the best window per allele and its 0-based
`prediction_offset` in the reported sequence. Class II is scored as the reported sequence, with its
allele-dependent register. `core` and `core_offset` come from the same model register, not an
allele-independent register guess. Class-I cores follow mhcmatch's footprint: eight residues for
an 8-mer, nine for a 9-11-mer, with central insertions omitted. Such a core need not be a contiguous
substring, and should not replace the full peptide in a structure model.

`tcr_facing` is the scored peptide with mhcmatch's class-default anchors masked by `X`.
`core_tcr_facing` applies the core residue mapping to that sequence, excluding class-II flanks.
These are predicted representations, not measured minimal recognition epitopes or measured contact
maps. A shared core under the same molecule supports a comparison between reported peptides; it
does not establish equivalent TCR recognition. Flanks and alternative registers can still matter.
No record, pMHC id, motif group or curation field is changed by this table.

`assessment_status` states coverage. `scored` rows contain the presentation percentile, calibrated
probability and mhcmatch class-specific band. Other rows name the absent reference, unsupported
peptide/species/class, absent panel allele, empty reference panel, or unscorable pair. Optional
measurements are text, with empty string as
the missing value; select scored rows and cast to numeric types for calculations. The generated
[column reference](standards/columns.md) declares every field.

Build-time predictions use the published mhcmatch version and revision/checksum-pinned reference
declared in `rules/epitope_assessment.toml`. Human and mouse class I/II are supported. Fetching the
reference is a separate input step; assembly is offline, recomputes predictions and calibration on
every run, and disables persisted calibration results. The fitted models bundled with mhcmatch are
immutable model inputs. Each row records the software version, reference SHA256, calibration seed,
background and footprint. CI runs:

```bash
uv run vdjdb epitope-reference --out out/inputs
uv run vdjdb build --out out/ --pmhc-reference out/inputs/pmhc/pmhc_full.tsv.gz --epitope-jobs 4
```

The same table can be recomputed independently with `vdjdb assess-epitopes`, over explicit chunk
files or `--tables`. See [check scopes](builds.md#choose-the-scope-of-a-check); selected-input
support counts are not full-corpus measurements.

Without `--pmhc-reference`, the table still catalogs all reported pairs and marks
`reference_not_supplied`; it makes no predictions. `--epitope-jobs` budgets processes over distinct
MHC species/class groups, each with one native thread. Results are sorted independently of worker
count. The older reviewed class-I promiscuity input and its `restriction` columns retain their
existing meaning; this table provides detailed freshly computed assessment alongside them.

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

### 3.6 Motif files

The primary bundle currently includes the same `cluster_members.txt` and `motif_pwms.txt`
files as the legacy bundle (§2), plus their optional TCREmp counterparts. It does not ship
`clusters.parquet` or `motifs.parquet`.

Run `vdjdb motifs --tables out/tables` to produce these files under `out/motifs`.
The column orders in §2 apply to both bundles. Background data are inputs; only derived
counts and statistics ship. The per-epitope diagnostic is
`reports/motifs_per_epitope.tsv` (§5).

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
| `motif-timings.tsv` | one row per motif stage: the same columns as `build-timings.tsv`, written by `vdjdb motifs`. `parent` is what makes it add up here: `motifs.tcrnet.background.*` runs inside `motifs.tcrnet.enrichment`, and while both rows carried their inclusive duration the report summed 268.7 s over a 198 s step and every `share` was understated by 36 %. **This is where the pipeline's memory peak is** - 6,898 MiB against the assemble stage's 1,577, set by the two PWM-and-emit steps, so 43 % of a 16 GB runner. Five measurements span 6,786 to 7,531 MiB over two hosts, an 11 % spread on the runner alone, and the 10,240 MiB budget is set against the largest of them. Share gated against `rules/motif_timings.tsv`, peak against a 10,240 MiB budget. Unlike the assemble report, the shares here are **not** comparable across host classes: per-stage slowdown from a 16-core laptop to the 4-vCPU runner ranges from 1.46x to 24.0x, the outlier being the background fetch, so the baseline carries the core count it was recorded on and the comparison runs on that host class only |
| `lookalikes.tsv` | one row per spelling of a value that some other value differs from only in case or in a `-`, `_`, `.` or space: the column, the folded form both share, the record count, and `same.species`. **Advisory, gated by nothing.** Measured on the current corpus: 19 groups over `antigen.gene` and `method.identification`, 9 of them within one species - the `reference.id` pair went when `PMID: 34433824` lost its space (#637). A value one character from another may be a typo or may be two stains, two serotypes, or one gene under two species' symbol conventions - HGNC capitalises `MBP` for human and MGI title-cases `Mbp` for mouse - so `same.species` is the column to sort on and the curator is the only one who can decide. The fold deliberately keeps `+` and `-`: applied to `meta.cell.subset` a stripping fold reports 11 groups and every one is false, because `CD95-` and `CD95+` are two populations. Edit distance was measured and rejected - at ratio 0.85 it returns real gene families (`MAGE-A1`/`A2`/`A3`/`A4`, `PPM1`/`PPM1F`) and distinct organisms (`CMV`/`MCMV`/`LCMV`) |
| `epitope-sources.tsv` | one row per `(antigen.epitope, antigen.species)` whose source is not single-valued: the modal `antigen.gene` the catalogue kept, `sources` (how many species carry this peptide), `genes` (how many labels this peptide-and-species carries), the record count, and `conserved`. **Advisory, gated by nothing.** Written by `vdjdb build` every run. Measured on the current corpus: 12 peptides under more than one species, 6 of them conserved between two proteomes and marked as such so a new collision is the only thing to read, and 30 rows carrying more than one gene label, the largest carrying 2. `build_epitopes` keeps the modal label and its comment claimed a function listed the rest; that function was never written, so until this report the discarded labels were visible nowhere (#633) |
| `anchors.tsv` | one row per chain whose `cdr3` contradicts the germline anchor of the V or J it names. VDJdb's `cdr3` is junction space - Cys104 through Phe/Trp118, both included - so a submission in IMGT CDR3 space, one carrying V or J framework past an anchor, or one with a mis-read anchor residue is none of those, and no rule in `vdjdb qc` tests for it: those check the residue alphabet and a minimum length. Written by `vdjdb build` on every run and never a gate. Columns: the chunk and row, the record, the shipped and submitted junctions, the two calls, arda's `v.canonical`/`j.canonical`, the germline anchor residue at each end, the named `defect`, a proposed `repair` sequence where the germline supports one, and a proposed `repair.call` where the junction matches a functional sibling allele instead. Measured on the 2026-09-29 corpus, after #646 repaired 4,838 junctions in `chunks/` and #647 corrected the mouse `TRAJ47` allele: **1,037 of 285,950 chains (0.36 %), 261 with a sequence repair and none with an allele repair** - 672 unexplained, 256 carrying framework past an anchor, 118 `unanchored`, 10 short a J anchor and 5 with a mis-read one. A chain can carry a defect at each end, so those count more than 1,037 between them. The anchor residue is read per segment and never assumed to be F or W: mouse `TRAJ47*01` is `HYANKMIC`, so "ends in Phe or Trp" calls 481 correct junctions broken. `unanchored` is the one case where that universal pair is still the test: with no J call, or one the reference does not have, there is no germline to read, and the check used to decline silently - those 118 chains were reported by the retired build's `vdjdb_full_cdr3aa_broken.txt` and by nothing here. No repair is proposed for them, because there is no germline to propose from |
| `harmonisation.tsv` | one row per value the build **rewrote**, across the four harmonisation passes `build_master` runs before identity is assigned (#700). The complement of `nomenclature.tsv` below, which names the calls it could *not* resolve: this one names the ones it could, and until #700 the four reports were computed on every build and discarded. Columns: `stage` (`segments`, `alleles`, `mhc`, `references`), the `issue` or patch table the rewrite cites, the `column`, the `species` where the pass scopes by organism, `from`, `to` and the record count. Written by `vdjdb build` on every run, **never a gate** - a rewrite is a declared correction, and the gate on those is `rules/expected_diffs.toml`'s generated `[[rename]]` block. Measured 2026-09-29: **270 rewrites over 7,450 records** - 235 segment respellings (4,146 records), 28 MHC corrections (1,494, of which the class II chain order is 149 records with a beta gene in `mhc.a`), 4 allele calls read off the junction (1,142, the TRAJ24 and mouse TRAJ47 signatures) and 3 `reference.id` values replaced by their PubMed id (668) |
| `j-calls.tsv` | one row per chain whose J call some **other** J gene explains better (#681). A junction's 3' end is templated by the J germline, so a record's J call can be checked against the sequence the same record reports: match the end - the anchor excluded, so a corrupt anchor cannot vote on its own diagnosis - against every J germline of that species and take the longest run. **This is not the question `anchors.tsv` asks.** That one reads the anchor residue of the segment a record *names*; this one asks whether a different gene explains the whole end better, and measured 2026-09-29 the two overlap on **19 of 718** chains. Columns: the chunk and row's record, the chain's gene and species, the junction, the call as written with the run its own germline supports, and the gene the sequence names with its run. Written by `vdjdb build` on every run, **never a gate** - #681 says re-calling a J from its junction is a curator's decision and what was missing is the list. The rule is stated rather than tuned: the best-matching gene must be unique (a tie is a fact about the germlines, not about the record), its run at least 5 residues (a chance run of 5 on a 20-letter alphabet is ~3e-7 per germline), and the called gene's own run strictly shorter. Measured: **715 chains over 96 chunks** against 285,545 checkable ones, so 99.75 % of J calls are not contradicted by their own junction. The largest single chunk is 106 chains of `PMID_32184241.tsv` |
| `nomenclature.tsv` | one row per segment call IMGT has at neither allele nor gene level, after harmonisation (#389). This is the report the retired build wrote as `vdjdb_full_gene_broken.txt` and `vdjdb_full_allele_broken.txt`, and nothing replaced it: `harmonise_segments` reports what it *rewrote*, `build_master` discards even that, and a call naming a gene no authority carries reached every shipped table with no report anywhere. Columns: the species, the chunk column, the call as written, the `part` that failed (a curator recording two candidates writes `TRBD1,TRBD2`, and each member is resolved separately), the chain count, and `family.members` with up to four `candidates`. Written by `vdjdb build` on every run, **never a gate**. Measured 2026-09-29: **88 rows over 69 distinct names and 3,457 chain-calls**, of which 2,823 are an under-specified *family* - `TRBV6` is nine IMGT genes and the record chose none of them, which only a curator can resolve - and the rest have no IMGT candidate at all, which is a spelling defect or a gene that species does not have. 850 are not human, and the retired check never saw one of those: its table was a 741-row human immunoglobulin list, so the driver ORed both masks with `species != 'HomoSapiens'`. It also compared `int(allele)` against a per-gene allele count rather than asking IMGT, which is a range check wearing the clothes of a membership check - its single finding on the whole corpus, `TRBV28*02`, is an allele IMGT lists. `tests/release/test_legacy_proofreading_parity.py` partitions every retired finding against this report and `anchors.tsv` with no remainder |
| `presentation.tsv` | one row per `(antigen.epitope, mhc.a, mhc.b, mhc.class)` pair that fails at least one of four checks on whether the recorded MHC could present the recorded epitope at all (ROADMAP phase 9e), and `presentation-summary.tsv` the same counted per finding. `proofreading/mhc_alleles.tsv.gz` answers whether IPD-IMGT/HLA lists a name; this asks whether the name reaches a **binding groove**, which is the 34 residues every presentation model reasons over. `mhcmatch` bundles those pseudosequences - 20,082 class I keys and 11,048 class II, loaded in 0.01 s with no network - so unlike `vdjdb promiscuity`, which fetches a model, this runs inside the build without breaking its offline determinism. The four: the call reaches a pseudosequence key; the key's own class agrees with `mhc.class`; a class I record's epitope fits a class I groove; and one molecule is not filed under two classes. **Advisory, never a gate** - a model is evidence about a pair, never authority over a publication. Measured 2026-09-29: **93 pairs over 1,213 records** of 2,343 and 192,641. 71 reach no groove, of which 64 are murine class II (`mhcmatch`'s class II pseudosequences are HLA, so this is a coverage statement rather than a finding against the record) and 2 are `H2-Qa-1b` and a four-field HLA spelling; 22 are a class I record whose epitope runs 12 to 20 residues; and 21 are `H2-IAb`, filed `MHCII` on 20 pairs and `MHCI` on the 77-record `QVYSLIRPNENPAH`, which is the one all three other checks agree on. An allele resolved by prefix is carried in `mhc.resolution` and is **not** a finding: 87 pairs over 17,836 records are an allele *group* the specification allows, `HLA-A*02` completed to its first member, so `mhcmatch` is guessing rather than the record being wrong |
| `functionality.tsv` | one row per chain-segment IMGT does not call functional, and `functionality-summary.tsv` the same counted per verdict (#634). `proofreading/imgt_alleles.tsv.gz` has carried IMGT's own F / ORF / P column since phase 9 and nothing read it: `vdjdb qc` asks whether a call *looks* like a TRBV name and `curate/nomenclature.py` asks whether IMGT *has* it, and neither asks whether IMGT thinks the gene is functional. Columns: the record, the chain, `V` or `J`, the call, IMGT's verdict as IMGT spells it, and `level` - `allele` where IMGT names that exact allele, `gene` where the verdict is inherited from the gene's alleles. Written by `vdjdb build` on every run, **never a gate**: a pseudogene V call is not automatically wrong, because a P gene can rearrange - `TRBV21-1` is 303 chains and turns up in real repertoires - and IMGT reclassifies genes between releases, so a gate would fail on a reference update rather than on a curation error. Measured 2026-09-29: **2,608 chain-segments, 1,172 V and 1,436 J**, largest `TRAJ58*01` ORF on 656 chains, `TRBJ1-6*01` ORF on 327 and `TRBV21-1*01` P on 303. The same verdict drives four advisory `non-functional *` rules in `vdjdb qc`. A gene's verdict is every verdict among its alleles and **one functional allele is enough**: `TRBJ2-7` reads `F/ORF` because `*02` is an ORF and it is one of the commonest J calls in VDJdb, so the stricter reading flagged 17,891 chunk rows with nothing actionable in the difference |
| `build-timings.tsv` | one row per build stage: the parent stage it sits inside, wall seconds **exclusive of any stage timed inside it**, share of the recorded total, parent peak RSS in MiB, optional `peak_tree_rss_mb` sampled across parent and descendants, the record count and the core count. Seconds sum to the wall clock and shares sum to 1, which they did not while a nested stage was counted both in its own row and in its parent's - the motif report summed 268.7 s over a 198 s step and understated every share by 36 %. Written by `vdjdb build` on every run. The share is gated against `rules/build_timings.tsv` (`tests/release/test_build_timings.py`); the seconds are recorded and not gated, because they are a property of the host, so a uniform slowdown is visible in the artifact rather than caught by a bar. Peak RSS is gated absolutely at 4,096 MiB for this stage, because memory is a property of the data and the code rather than of the host; measured 1,577 MiB on a laptop and 1,064 MiB on the runner, which allocates less because polars chunks to fewer threads. The share baseline is recorded on the runner and carries its core count: measured, the nine shares here agree to within 1.3 points between 4 and 16 cores, which is why a share gate works for this report. ⚠ The assemble stage is **not** the pipeline's memory peak - see `motif-timings.tsv`. `annotate.junction.add_junction_nt` was 87.2 % of the wall time (`antigenomics/vdjtools#181`); since `vdjtools` 4.5 it is one `infer_nt_batch` call per (species, locus) and **65.2 %** - 108.54 s of a 166.54 s stage on the runner, against 407.64 s of 467.72 s. See `docs/builds.md` |
| `motif-metrics.tsv` | one row per (species, gene, source, axis): eighteen axes for four sources -- the last legacy release, the latest release, this build's two methods, and the partition that clusters nothing. 180 rows. Thirteen axes score a clustering on this build's cohort (two of them, `percolation_median` and `percolation_excess`, are recorded and never gated, because the largest cluster's share depends on how many prominent motifs the epitope has); the five `partition_*` axes compare it against **the shipped file**, asking whether a record still has the cluster-mates the last release gave it. Those exist because nothing compared those tables: `vdjdb diff` keys `cluster_members.txt` on `cid`, a cid carries a position in a sorted list, so one renumbered cluster reads as the entire file replaced. Only `partition_neighbours_preserved` is gated, and the reason is measured: 19,971 of the released TRB clustering's 36,906 clonotypes sit in one cluster holding 94.7 % of the file's co-clustered pairs, so a pair-weighted score measures that one blob -- the do-nothing partition reads 0.9991 on it. Written by `vdjdb motif-metrics`, which also gates them: `current-*` against `latest` catches a code regression, every source against `rules/motif_metrics.tsv` catches a corpus one. The metrics used to live only inside test assertions, so a corpus change moved them inside the slack and nobody learned the new values |
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

`registry/records.tsv.gz` maps `record_id` to its state, hashes, provenance and amendment history, so an
id survives a curator fixing a typo. **It is committed** and is an input to the build, not an output of
it. Git LFS stores a deterministic gzip of the TSV. Decompression preserves the complete table, including every identifier and its history. Gzip metadata uses a zero timestamp and no filename, so identical inputs produce identical compressed bytes.
Run `git lfs pull` after cloning; CI checks out the LFS contents before building. It was going to ship as a release asset the build fetched, and #672 changed
that - without a committed copy the registry went stale between releases and landing one 40-record
chunk moved `record_id` on 168,723 of 192,753 records (#638). It is written **only** by
`vdjdb identity update`, which a chunk branch runs and commits alongside the chunk; `vdjdb build` reads
it.

Fifteen columns. Thirteen are the state, the two hashes, the chunk provenance, the release and commit
at first and last sighting, and the amendment count with the key hash it came from. The other two are
the lifecycle a consumer follows:

| Column | What it answers |
|---|---|
| `amended_from_key_hash` | backwards, and only for an amendment: which key this record used to have |
| `replaced_by` | forwards, and only for a retirement the amendment pass refused: which id took over. Two key fields moving is a new record by the rule the registry is built on, so the retirement is right and the pointer is what makes it diagnosable rather than a disappearance (#693, ROADMAP §10.4). Written when one id retires from a line of a chunk, exactly one is allocated against that same line, and the two natural keys name the same receptor; empty otherwise, including on every retirement that is a genuine deletion |

## 6a. Committed, not produced

| File | Role |
|---|---|
| `summary/reference_years.tsv` | publication years, so the dashboard render is offline and deterministic. A committed, reviewed input refreshed by its own pull request, never written by a build |
| `summary/annotations.tsv` | dashboard event callouts, with no hardcoded coordinates |
| `rules/expected_diffs.toml` | the differences `vdjdb diff` accepts, each declared with a measured row count |
| `proofreading/epitope_proteome.tsv` | where each **self** epitope sits in its own species' reference proteome, and how exactly (#632). Written by `vdjdb antigens`, never by a build: `mhcmatch` fetches the proteome from HuggingFace, so this is a committed, reviewed input refreshed by its own pull request. Three verdicts, spelled as what was measured. **`exact`** - the peptide is in the proteome at a named protein and position, and its `GN=` field gives the gene symbol, which is the authority `antigen.gene` has never had (315 human epitopes, 15 mouse). Where that symbol disagrees with the curated value, 99 epitopes over 2,195 records, that is worth a curator's time. **`one_substitution`** - one residue differs from a peptide that is in the proteome, named as `9C>V` (195 human epitopes, 16 mouse). **Not a defect and not a count to read as one.** A reference proteome is one genome and a patient cohort is not, so at least five unrelated things produce this shape and all are correct data: a tumour neoantigen, a germline or allelic variant, a cross-species homolog, a modified or hybrid peptide, and an anchor-optimised screening reagent. Read the epitope count and never the record total - one peptide, `SLLMWITQV`, is 81.5 % of that total, and the median epitope carries two records. Read `by_gene` too: `PMEL` carries 21 peptides, `INS` 13, `KRAS` 12 over 50 records, which is what a mutation panel or an antigen screen across a cohort looks like. The row's use is the **cross-reference** - `SLLMWITQV` and native `SLLMWITQC` are unrelated rows, so nothing can currently ask whether a response was found with the native peptide or a modified one - plus the `reference.id` list, because only the paper tells the five causes apart and `corpus/pubmed.tsv` has the title for 610 references. **`not_found`** - neither within one substitution; a splice junction, a fusion, a longer modification and a transcription error all look like this (236 human, 13 mouse). Nothing is ever rewritten. Not to be confused with `out/reports/epitope-sources.tsv` (#633), which asks whether the *corpus* agrees with itself about a peptide's source and needs no authority at all |

The pinned motif parameters are not a data file. They are the `TUNED` constants in
`vdjdb.motifs.tcrnet` and `vdjdb.motifs.tcremp`, each stated with the scorecard it was chosen from
and the one cost it pays, in the module that uses it.

---

## 7. Never produced, never shipped

**TCRvdb / MATCHMAKERS** (Messemaker et al., doi:10.1101/2025.04.28.651095) is used only as a held-out validation set, read from a path given by `VDJDB_TCRVDB` and never from
inside this repository. No file in any bundle, artifact or report may contain its rows, its per-record
labels, or any value derived from them at record granularity. Aggregate validation metrics (counts,
AUROC, recall at a threshold) may be reported; per-record verdicts may not.

Motif clustering is tuned on the independent-study support count: 4,129 human clonotype-epitope pairs
with ≥2 distinct `reference.id`, across 18 epitopes with ≥20 pairs. TCRvdb is read once, at the end,
to validate. Tuning on it would invalidate it as a validation set, and its 614 labels span only two
epitopes (YLQPRTFLL, GLCTLVAML), too narrow a basis for a global hyperparameter.

CI enforces this: a guard fails the build if any TCRvdb-derived file appears in the repository or in a
release bundle.
