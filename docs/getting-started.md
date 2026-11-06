# Explore your first records

This tutorial downloads a release and reads a few observations. You need an archive extractor
and Python 3; no database server is required. To search without downloading anything, use
[vdjdb.com](https://vdjdb.com).

## 1. Download and extract a release

Open the [release list](https://github.com/antigenomics/vdjdb-db/releases), choose a release,
and download the bundle whose manifest role is `legacy`. Extract it into a directory you can open in a terminal.
Keep the release tag with your analysis so another reader can identify the data you used.

## 2. Choose a table

For this example, select the **legacy** bundle and use `vdjdb.txt` or `vdjdb.txt.gz`: it has one row per reported chain. Alpha and beta chains of a
paired receptor can therefore occupy two rows. `vdjdb_full.txt` puts paired chains in one row;
`vdjdb.slim.txt` is a reduced projection. See the [file reference](outputs.md) before treating
rows from different formats as interchangeable observations.

## 3. Read the first five rows

Run this command in the extracted directory containing `vdjdb.txt` or `vdjdb.txt.gz`.
It streams either form without creating an uncompressed copy:

```bash
python3 - <<'PYTHON'
import csv
import gzip
from pathlib import Path
from itertools import islice

path = Path("vdjdb.txt.gz")
if not path.exists():
    path = Path("vdjdb.txt")
opener = gzip.open if path.suffix == ".gz" else open
with opener(path, "rt", encoding="utf-8", newline="") as handle:
    rows = csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE)
    for row in islice(rows, 5):
        print(row["gene"], row["cdr3"], row["antigen.epitope"], row["reference.id"])
PYTHON
```

Each output line identifies a chain, its junction sequence, the tested epitope and the publication.
The exact values depend on the release. A repeated receptor sequence can describe distinct donors,
assays or publications; do not discard rows on sequence alone.

## Next steps

- [Query receptors, retrieve motifs, inspect reliability and fetch coordinates](using-vdjdb.md).
- Look up [columns](standards/columns.md) and [confidence scores](standards/confidence-score.md).
- [Submit records](submission.md) from a publication.
- [Build the database](builds.md) from the source chunks.

## Citing

Cite the most recent paper: Daniil V. Luppov, Anna E. Koneva, Dmitry V. Bagaev, Anastasiia V. Alexandrova, Elizaveta K. Vlasova, Dmitry M. Chudakov, Chihiro Motozono, Andrew K. Sewell & Mikhail Shugay. VDJdb in 2026: boosting T-cell receptor recognition evidence using paratope embeddings and AI-based structure prediction. *Nucleic Acids Research*, 2026. [doi:10.1093/nar/gkag904](https://doi.org/10.1093/nar/gkag904)
