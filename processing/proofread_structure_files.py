"""Proofread the pinned structure dataset against completed current VDJdb tables.

Run via the repository environment. --fetch stages immutable hosted inputs; the
default verifies an existing local mirror. No source records or coordinates change.
"""
from pathlib import Path
import argparse
import json

from vdjdb.validate.structure_files import run

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--structures', type=Path, required=True)
    parser.add_argument('--tables', type=Path, default=Path('out/tables'))
    parser.add_argument('--out', type=Path, default=Path('out/reports/structures'))
    parser.add_argument('--manifest', type=Path, default=Path(__file__).resolve().parents[1] /
                        'proofreading/structure_inputs.tsv')
    parser.add_argument('--fetch', action='store_true')
    args = parser.parse_args()
    print(json.dumps(run(args.structures, args.tables, args.out, args.manifest, fetch=args.fetch), indent=2))
