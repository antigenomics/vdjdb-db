# AIRR mapping

`vdjdb convert airr` emits AIRR Rearrangement, Receptor and Reactivity from the definitive
tables. The mapping below is declared once, on the field registry, and both converters -
new-format and legacy - read it from there.

An empty mapping means that no AIRR field holds the same quantity, not that no related
field exists. VDJdb's `v.end` is an amino-acid offset in junction space, while AIRR's
`v_sequence_end` is a nucleotide offset in sequence space. `vdjdb.convert.coords` is the
only module that converts between the four coordinate spaces used in this codebase.

The explicit `v_end_nt`, `j_start_nt`, `d_start_nt` and `d_end_nt` fields are also
junction-relative, with 0-based half-open bounds. They cannot be copied directly
into AIRR sequence fields, which use 1-based closed bounds in the full sequence.
For a junction starting at 0-based sequence offset `o`, a junction interval `[s, e)`
becomes the AIRR interval `[o + s + 1, o + e]`. The origin must be known from the
sequence being exported; an inferred junction alone does not supply it.

The [CDR3 markup guide](cdr3-fixing.md) explains anchor inclusion, partial codons,
all eight nucleotide/amino-acid fields, missing values and legacy compatibility.

```{vdjdb-vocabulary} airr
```

## Species vocabulary

```{vdjdb-vocabulary} species
```

## Declared tables

```{vdjdb-vocabulary} tables
```
