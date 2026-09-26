# AIRR mapping

`vdjdb convert airr` emits AIRR Rearrangement, Receptor and Reactivity from the definitive
tables. The mapping below is declared once, on the field registry, so both converters --
new-format and legacy -- fall out of one statement rather than two hand-written lists.

**An empty mapping is not "no counterpart exists"** but "none that is the *same*
quantity". VDJdb's `v.end` is an amino-acid offset in junction space and AIRR's
`v_sequence_end` a nucleotide offset in sequence space; declaring them equal would be a
lie. `vdjdb.convert.coords` is the one module allowed to convert between the four
coordinate spaces this codebase meets.

```{vdjdb-vocabulary} airr
```

## Species vocabulary

```{vdjdb-vocabulary} species
```

## Declared tables

```{vdjdb-vocabulary} tables
```
