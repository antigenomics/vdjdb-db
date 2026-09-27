# Identifiers

A VDJdb record is a receptor against a presented peptide. Both halves recur across records, and
consumers reference both halves, so each level carries an identifier of its own.

| Level | Identifies | Id | Where it appears |
|---|---|---|---|
| clonotype | one receptor chain | `CT` + 16 hex digits | `chains.clonotype_id` |
| clone | one alpha/beta pair | `CX` + 16 hex digits | `chains.clone_id` |
| pMHC | one peptide as presented | `PM` + 16 hex digits | `records.pmhc_id` |
| epitope | the peptide alone | `EP` + 16 hex digits | `records.epitope_id` |
| record | one curated line | `VDJDB` + 10 digits | `records.record_id` |

Counts on the 2026-09 build: 187,935 clonotypes over 286,047 chain rows, 82,266 clones over the
93,294 records reporting two chains, 2,364 pMHCs, 2,118 epitopes and 192,753 records.

## What each one keys on

| Level | Key |
|---|---|
| clonotype | `species`, `gene`, `cdr3`, `v.segm`, `j.segm` |
| clone | the record's two `clonotype_id`s, sorted |
| pMHC | `antigen.epitope`, `mhc.a`, `mhc.b` |
| epitope | `antigen.epitope` |
| record | the record natural key: the complex-information columns, the id fields, and `chunk.file` |

`mhc.class` is not in the pMHC key, because the alleles determine it: keying on it as well gives the
same 2,364 distinct values.

A record reporting one chain has no clone, and `clone_id` is the **empty string** rather than a null.
Empty is the only missing marker in VDJdb, so a null in a string column would be the ambiguity that
several shipped bugs came from. 99,459 of 192,753 records are in this state.

## Derived and allocated

Four of the five ids are **derived**: the id is a hash of its own key and of nothing else. One,
`record_id`, is **allocated** from a counter against a registry.

The difference answers the question a curator asks first, which is what happens when a chunk is
added or removed. For a derived id, nothing:

- the order `chunks/` is read in cannot reach it;
- adding a chunk names that chunk's new clonotypes and changes no existing id;
- removing a chunk retires only the ids nothing else supported;
- two hosts, two core counts and two release dates produce the same ids.

A counter has none of those properties. `record_id` is allocated anyway, because its purpose is to
survive a content change: a curator fixing a CDR3 typo has to keep the record's id, and no hash of
the content can do that. Every other level keys on content with no separate existence, so a change
there is not an amendment but a different clonotype.

One consequence: exactly one registry is consulted during a build, the record registry, and the four
derived levels need no history in order to be correct.

## The hash

sha256 over the key fields joined by `\x1f`, truncated to 16 hex digits, with the level's prefix in
front. The separator is the unit separator, so a field containing a tab cannot forge a key boundary.

sha256 and not a library hash function. polars' `Expr.hash` is xxhash and its output is not specified
across versions, so a dependency upgrade could renumber every clonotype with nothing failing
anywhere. The algorithm has to be one this repository owns, and
`tests/unit/test_identity_levels.py` freezes one id against a literal so it cannot move.

Truncating to 16 digits is a size decision rather than a security one. At 187,935 clonotypes the
chance of one collision is around 1 in 10⁹, and the invariants below fail a build rather than letting
one corrupt a join. Measured cost on the current build: 83 ms to hash the 187,935 distinct clonotype
keys and 7 ms to join them back onto 286,047 chain rows.

## Lifecycle

An id that vanishes is the failure a consumer cannot diagnose: a reference that used to resolve
returns nothing, and nothing distinguishes a correction from a deletion. So every level keeps a
lifecycle row, and the row holds only what a build cannot recompute.

| Field | Meaning |
|---|---|
| `id` | the identifier |
| `level` | `clonotype`, `clone`, `pmhc`, `epitope` or `record` |
| `state` | `active` or `retired` |
| `first_release` | the release tag the id first appeared in |
| `last_release` | the newest release carrying it, frozen at retirement |
| `replaced_by` | the id that took over, where the retirement was an amendment |

```bash
uv run vdjdb identity resolve CT3efa364e7f84aac9 --lifecycle identity-lifecycle.tsv
```

```
id              CT3efa364e7f84aac9
level           clonotype
state           active
first_release   v2026.09.1
last_release    v2026.09.1
replaced_by
```

**Retirement is per release, not per build.** A curation branch may add a clonotype and remove it
again before anything ships, and neither event is a lifecycle event: only the release job writes
these rows. That is also why a curation pull request never touches them.

**An id that comes back is the same id.** A derived id is the hash of its key, so the same key
returning has to produce the same id. It goes back to `active` keeping its original `first_release`,
and `vdjdb identity diff` reports it under `returned` rather than under `added`, because the two mean
different things to a reader: one has been published before and one has not.

The lifecycle table is a release asset rather than a committed file, at 15.9 MB for 274,683 derived
ids. The record registry is separate and larger, at 72.7 MB, because it carries each record's natural
key and content hash. The retired rows alone are committed, since they are the part a consumer needs
and the part that is otherwise lost.

A build with no lifecycle file still runs and produces exactly the same derived ids. Only the history
is unknown, which is the correct state for a fork and for a first build.

## Invariants

`vdjdb identity check` asserts six of the seven, exits 1 on any failure, and runs in `build.yml`
after the assembly step.

| # | Invariant | Needs |
|---|---|---|
| 1 | no id is shared by two distinct keys | the build |
| 2 | every stored id equals the hash of the key on its own row, recomputed | the build |
| 3 | a permuted chunk order, and a chunk added then removed, change no id | two builds |
| 4 | every `record_id` whose natural key is unchanged is unchanged | the previous registry |
| 5 | no id carries a prefix belonging to no level, and no published level empties silently | the previous lifecycle |
| 6 | every id referenced resolves: a clone covers exactly two chains of different `gene`, every chain names a record, every evidence row names a chain, `pmhc_id` and `epitope_id` are never empty | the build |
| 7 | `TCR_hash` is byte-identical wherever the clonotype key is unchanged | the previous tables |

Invariant 3 is a property of two builds rather than of one directory, so it is a unit test over a
fixture in `tests/unit/test_identity_levels.py` and costs milliseconds instead of two full builds.

Invariants 2 and 5 together are what rules out reuse. Reuse would be an id resolving to a different
key than the one it was retired under, and invariant 2 recomputes every id from the key beside it, so
an id pointing anywhere else fails. An id merely reappearing is not reuse.

Checks needing history are skipped rather than failed when the history is absent.

```bash
uv run vdjdb identity check --tables out/tables
uv run vdjdb identity check --tables out/tables --previous identity-lifecycle.tsv \
    --previous-tables previous/tables
uv run vdjdb identity diff identity-lifecycle.tsv --tables out/tables --release v2026.10.1
uv run vdjdb identity lifecycle --release v2026.10.1 --tables out/tables \
    --previous identity-lifecycle.tsv --out out/release/identity-lifecycle.tsv
```

## MHC promiscuity is an annotation, never part of a key

An epitope is often presented by several alleles, and the curated allele is not always the one that
binds it best. Both facts belong in the database; neither belongs in a key.

`pmhc_id` keys on the allele **as curated**. A prediction is a moving target: `mhcmatch` ships new
weights, and if a predicted allele were in the key then a model upgrade would renumber pMHC ids and
break every external reference while the build passed. Promiscuity is therefore measured into columns
on `restriction`, each row carrying the `mhcmatch` version that produced them, so re-running under a
new model rewrites those columns and no id.

## Legacy structure ids

`TCR_hash` is the identifier linking a record to a generated structure, and it ships exactly as
curated. It is stored as is: never recomputed, never renamed, never derived from the levels above.
Structure evidence keys on it. The identifiers on this page are additive, so a consumer holding a
`TCR_hash` is never asked to migrate, and invariant 7 fails the build if one moves on a clonotype
whose key did not.
