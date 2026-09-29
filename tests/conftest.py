"""Puts ``tests/`` on ``sys.path`` so both tiers can import :mod:`legacy_qc`.

pytest inserts a conftest's own directory under the default ``prepend`` import mode, which is the
whole reason this file exists: ``tests/unit/`` and ``tests/release/`` both measure the new rules
against the same transcription of the retired build, and one copy of that transcription is the point.
"""
