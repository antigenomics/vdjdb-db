---
name: vdjdb-format
description: Normalise the controlled-vocabulary fields of a raw VDJdb chunk TSV to the form the database records - IMGT V/D/J gene and allele names including Adaptive ImmunoSEQ and Arden conversions, IPD-IMGT/HLA and murine H2 MHC allele names, species CamelCase, method vocabulary, reference id prefixes - and write a format log naming the authority behind every change. Use after vdjdb-extract and before vdjdb-proofread, or on any TSV whose gene or allele spellings need bringing to IMGT.
---

# Validate submission formatting

Read [AUTHORITIES.md](../AUTHORITIES.md). Use the current
[chunk specification](../../docs/standards/chunk-format.md) and a shipping header rather than
positional column numbers or copied column lists.

Run `vdjdb qc <file>` first. Resolve malformed headers, wrong separators and shifted fields from
the original source. An invalid layout cannot be repaired by renaming values in the shifted cells.

## Inspect assembly before editing

Preserve supplied sequence and segment cells when the build already resolves them. Inspect the
harmonisation, nomenclature, junction and submitted-call reports. Do not manually repeat allele
selection, old-name conversion or terminal-flank repair at import.

| Field | Review |
|---|---|
| Species | Declared vocabulary, including `BosTaurus`; no borrowed germline model |
| V/D/J | IMGT membership and locus; Arden/Adaptive conversions from reviewed tables |
| Multiple segment calls | Preserve the reported alternatives; check each call |
| Human MHC | IPD-IMGT/HLA membership, reported resolution and expression suffix |
| Nonhuman MHC | `proofreading/mhc_nonhuman.tsv` and declared `patches/mhc.dict` corrections |
| Method tokens | Declared vocabulary and actual experimental evidence |
| Reference | Verified publication/submission identifier and reviewed aliases |

Check suspicious allele-number runs against IMGT and the source. A numeric cutoff is not evidence
that an allele is wrong. Nonfunctional calls are advisory; do not substitute a functional allele
without sequence or source support.

Use `H2-` murine names where the authority maps them. A chain name and a haplotype molecule name
are different levels of reporting. Low-resolution human calls remain low resolution unless the
source provides more detail.

Class I requires `mhc.b=B2M`. Class II cannot use `B2M`; a reported class-II chain may have its
unreported partner blank. Resolve supported names and inspect partner inference where available;
never infer restriction from peptide length or a peptide-binding prediction.

## Preserve frequency and assay meaning

A percentage for a clonotype is not automatically a subset frequency. Keep `method.frequency`,
its count/total and `meta.subset.frequency` separate. Do not fill one from another.
A platform name is not sequencing evidence by itself. A sorting term does not become verification
unless the paper reports the reconstructed receptor being tested.

Record every proposed cell change with file/row, old/new value and authority. A missing authority
entry requires a reviewed addition, not a guessed nearest spelling. Apply data changes on the
chunk or data-issue branch, separately from build code. Re-run QC and
[proofread](../vdjdb-proofread/SKILL.md).
