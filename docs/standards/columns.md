# Column reference

Every table below is built from `vdjdb.schema` while this page renders. There is no
generated file to refresh, so no copy of these tables can disagree with the build.

## `vdjdb.txt` - the release table

```{vdjdb-schema} vdjdb
:columns: name, title, comment
```

## `vdjdb.txt` - the vdjdb-web variant

The variant served by `vdjdb-web` adds five `evidence.*` columns and a `TCR_hash` row that
the release table does not have.

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

`records` and `chains` are the database. Every shipped file is a join and a pivot away
from them.

### `records`

```{vdjdb-schema} records
:columns: ships_as, title, comment
```

### `chains`

```{vdjdb-schema} chains
:columns: ships_as, title, comment
```

### `epitopes`

```{vdjdb-schema} epitopes
:columns: ships_as, title
```

### `restriction`

```{vdjdb-schema} restriction
:columns: ships_as, title
```

### `evidence`

```{vdjdb-schema} evidence
:columns: ships_as, title
```

### `epitope_assessment`

```{vdjdb-schema} epitope_assessment
:columns: ships_as, title
```

## The motif files

⚠ **These two files are parsed positionally by `vdjdb-web`**, which hands Tablesaw a fixed
column-type array with no header check. An inserted, removed or reordered column mistypes
or shifts every column after it, without an error. Column order is a contract.

### `cluster_members.txt`

```{vdjdb-schema} cluster_members
:columns: name, title
```

### `motif_pwms.txt`

```{vdjdb-schema} motif_pwms
:columns: name, title
```
