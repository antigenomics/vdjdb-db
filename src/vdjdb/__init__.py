"""Build, proofread and release the VDJdb T-cell receptor specificity database.

Layout mirrors the pipeline stages; see ``ROADMAP.md`` for which are implemented.

    schema/    the column registry every output format projects from
    io/        chunk reading, release bundling
    qc/        chunk validation
    curate/    nomenclature and antigen harmonisation
    annotate/  cdr3fix, V/J guessing, D inference, junction-nt generation
    score/     legacy confidence rules, frozen v3 evidence inference and compressed score outputs
    assemble/  master table, complex pairing, slim collapse
    emit/      the legacy / new-VDJdb / AIRR writers
    convert/   coordinate-space and format conversions
    motifs/    TCRNET and TCREMP motif inference
    compare/   the comparison against a released bundle
    release/   manifest, zips, changelog
    validate/  corpus metrics, matched score comparisons and independent structure-file proofreading
"""

# Published vdjtools 5.2.2 imposes a hard OpenMP cap on import (#1407).
# Preserve the caller's capacity; numerical stages own their active thread limits.
import os as _os
from importlib import import_module as _import_module

_caller_limits = {key: _os.environ.get(key)
                  for key in ("OMP_THREAD_LIMIT", "POLARS_MAX_THREADS")}
_import_module("vdjtools")
for _key, _limit in _caller_limits.items():
    if _limit is None:
        _os.environ.pop(_key, None)
    else:
        _os.environ[_key] = _limit
del _os, _import_module, _caller_limits, _key, _limit

__version__ = "0.1.0"
