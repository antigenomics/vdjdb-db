"""Identifiers: the five levels, their lifecycle, and the invariants that hold them together.

:mod:`~vdjdb.identity.levels` derives the four content-keyed ids, :mod:`~vdjdb.identity.ids`
allocates ``record_id`` against a registry, :mod:`~vdjdb.identity.lifecycle` records what happened to
an id, and :mod:`~vdjdb.identity.checks` asserts the invariants. ``docs/standards/identity.md`` is the
reference page; ``ROADMAP.md`` section 10 is the design.
"""
from .checks import Finding, check
from .ids import (
    IdentityRegistry,
    RecordState,
    canonical_content_hash,
    natural_key,
    reconcile,
)
from .levels import (
    CLONE,
    CLONOTYPE,
    EPITOPE,
    LEVELS,
    PMHC,
    Level,
    clone_ids,
    derive,
    hash_id,
    level_of,
)
from .lifecycle import advance, compare, present, resolve

__all__ = [
    "CLONE",
    "CLONOTYPE",
    "EPITOPE",
    "LEVELS",
    "PMHC",
    "Finding",
    "IdentityRegistry",
    "Level",
    "RecordState",
    "advance",
    "canonical_content_hash",
    "check",
    "clone_ids",
    "compare",
    "derive",
    "hash_id",
    "level_of",
    "natural_key",
    "present",
    "reconcile",
    "resolve",
]
