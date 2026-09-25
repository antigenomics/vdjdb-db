"""Build, proofread and release the VDJdb T-cell receptor specificity database.

Layout mirrors the pipeline stages; see ``ROADMAP.md`` for which are implemented.

    schema/    the column registry every output format projects from
    io/        chunk reading, release bundling
    qc/        chunk validation
    curate/    nomenclature and antigen harmonisation
    annotate/  cdr3fix, V/J guessing, D inference, junction-nt generation
    score/     the VDJdb confidence score
    assemble/  master table, complex pairing, slim collapse
    emit/      the legacy / new-VDJdb / AIRR writers
    convert/   coordinate-space and format conversions
    motifs/    TCRNET and TCREMP motif inference
    compare/   the release difference ledger
    release/   manifest, zips, changelog
"""

__version__ = "0.1.0"
