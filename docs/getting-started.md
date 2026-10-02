# Explore your first records

This tutorial downloads a release and reads a few observations. You need an archive extractor
and Python 3; no database server is required. To search without downloading anything, use
[vdjdb.com](https://vdjdb.com).

## 1. Download and extract a release

Open the [release list](https://github.com/antigenomics/vdjdb-db/releases), choose a release,
and download its database zip. Extract it into a directory you can open in a terminal.
Keep the release tag with your analysis so another reader can identify the data you used.

## 2. Choose a table

For this example, use `vdjdb.txt`: it has one row per reported chain. Alpha and beta chains of a
paired receptor can therefore occupy two rows. `vdjdb_full.txt` puts paired chains in one row;
`vdjdb.slim.txt` is a reduced projection. See the [file reference](outputs.md) before treating
rows from different formats as interchangeable observations.

## 3. Read the first five rows

Run this command in the directory containing `vdjdb.txt`:

```bash
python3 - <<'PYTHON'
import csv
from itertools import islice

with open("vdjdb.txt", newline="") as handle:
    rows = csv.DictReader(handle, delimiter="\t")
    for row in islice(rows, 5):
        print(row["gene"], row["cdr3"], row["antigen.epitope"], row["reference.id"])
PYTHON
```

Each output line identifies a chain, its junction sequence, the tested epitope and the publication.
The exact values depend on the release. A repeated receptor sequence can describe distinct donors,
assays or publications; do not discard rows on sequence alone.

## Next steps

- Look up [columns](standards/columns.md) and [confidence scores](standards/confidence-score.md).
- [Submit records](submission.md) from a publication.
- [Build the database](builds.md) from the source chunks.

## Citing

Cite the most recent paper: Daniil V. Luppov, Anna E. Koneva, Dmitry V. Bagaev, Anastasiia V. Alexandrova, Elizaveta K. Vlasova, Dmitry M. Chudakov, Chihiro Motozono, Andrew K. Sewell & Mikhail Shugay. VDJdb in 2026: boosting T-cell receptor recognition evidence using paratope embeddings and AI-based structure prediction. *Nucleic Acids Research*, 2026. [doi:10.1093/nar/gkag904](https://doi.org/10.1093/nar/gkag904)
