# Compute authoritative cdr3fix (vEnd/jStart/...) for PDB_Database.txt native chains using the
# vdjdb-db Cdr3Fixer, so the vdjdb-web reconcile step can populate generated native rows with a
# real cdr3fix JSON. An empty cdr3fix cell crashes the search-table stream (Json.parse("")).
# Run from py_src/ (the fixer resolves ../res and ../patches relatively).
#   python compute_pdb_cdr3fix.py ../chunks/PDB_Database.txt pdb_cdr3fix.tsv
# 2026-07-19
import csv
import json
import sys

from Cdr3Fixer import Cdr3Fixer

CHAINS = [("cdr3.alpha", "v.alpha", "j.alpha"), ("cdr3.beta", "v.beta", "j.beta")]


def main() -> int:
    src, out = sys.argv[1], sys.argv[2]
    fixer = Cdr3Fixer("../res/segments.txt", "../res/segments.aaparts.txt")

    seen = {}  # (species, cdr3, v, j) -> cdr3fix json
    with open(src, encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            species = (row.get("species") or "").strip()
            for cc, vc, jc in CHAINS:
                cdr3 = (row.get(cc) or "").strip()
                v = (row.get(vc) or "").strip()
                j = (row.get(jc) or "").strip()
                if not cdr3:
                    continue
                key = (species, cdr3, v, j)
                if key in seen:
                    continue
                res = fixer.fix_both(cdr3, v, j, species)
                seen[key] = json.dumps(res.results_to_dict() if res else {})

    with open(out, "w", encoding="utf-8", newline="") as w:
        wr = csv.writer(w, delimiter="\t")
        wr.writerow(["species", "cdr3", "v", "j", "cdr3fix"])
        for (sp, cdr3, v, j), cf in seen.items():
            wr.writerow([sp, cdr3, v, j, cf])
    # self-check: every emitted cdr3fix parses and carries vEnd/jStart
    for cf in seen.values():
        d = json.loads(cf)
        assert "vEnd" in d and "jStart" in d, d
    print(f"wrote {len(seen)} cdr3fix entries to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
