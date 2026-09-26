# Column reference

Every table below is built from `vdjdb.schema` **while this page renders**. There is no
generated file to refresh, so these tables cannot disagree with the build.

## `vdjdb.txt` -- the release table

```{vdjdb-schema} vdjdb
:columns: name, title, comment
```

## `vdjdb.txt` -- as production serves it

Five `evidence.*` columns and a `TCR_hash` row that the release table does not carry.

```{vdjdb-schema} vdjdb-web
:columns: name, title, comment
```

## `vdjdb.slim.txt`

```{vdjdb-schema} slim
:columns: name, title, comment
```

## `vdjdb_full.txt`

```{vdjdb-schema} full
:columns: name, title
```

## The definitive tables

`records` and `chains` are the database; every shipped file is a join and a pivot away
from them.

### `records`

```{vdjdb-schema} records
:columns: name, title, comment
```

### `chains`

```{vdjdb-schema} chains
:columns: name, title, comment
```

### `epitopes`

```{vdjdb-schema} epitopes
:columns: name, title
```

### `restriction`

```{vdjdb-schema} restriction
:columns: name, title
```

### `evidence`

```{vdjdb-schema} evidence
:columns: name, title
```

## The motif files

⚠ **Parsed positionally by `vdjdb-web`**, which hands Tablesaw a fixed column-type array
with no header check. Any inserted, removed or reordered column silently mistypes or
shifts the whole table: column order is a contract, not a convention.

### `cluster_members.txt`

```{vdjdb-schema} cluster_members
:columns: name, title
```

### `motif_pwms.txt`

```{vdjdb-schema} motif_pwms
:columns: name, title
```
