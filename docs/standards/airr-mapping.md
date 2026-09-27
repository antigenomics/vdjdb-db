# AIRR mapping

`vdjdb convert airr` emits AIRR Rearrangement, Receptor and Reactivity from the definitive
tables. The mapping below is declared once, on the field registry, and both converters -
new-format and legacy - read it from there.

An empty mapping means that no AIRR field holds the same quantity, not that no related
field exists. VDJdb's `v.end` is an amino-acid offset in junction space, while AIRR's
`v_sequence_end` is a nucleotide offset in sequence space. `vdjdb.convert.coords` is the
only module that converts between the four coordinate spaces used in this codebase.

```{vdjdb-vocabulary} airr
```

## Species vocabulary

```{vdjdb-vocabulary} species
```

## Declared tables

```{vdjdb-vocabulary} tables
```
