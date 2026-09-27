# attic — retired, kept because something still depends on it

The Groovy implementation of the database build, retired in phase 4. It is **not** dead weight and
must not be deleted: `BuildDatabase.groovy` is the one correct specification for the two `*.meta.txt`
files and for deriving the `vdjdb.txt` header from the metadata, and `tests/unit/test_schema.py`
**parses it at test time** rather than copying it. A copy would drift exactly the way the nine
duplicated column lists did.

`AlignBestSegments.py` is the Groovy era's stage-II segment realignment. It sat in `src/` until
2026-09-27 and was never reachable: `BuildDatabase.groovy` was its only caller, the Python port
never invoked it, and it imports `Bio`, which this project does not depend on. It is here rather
than deleted because it is the only written form of what stage II did, and `arda` covers only part
of that -- arda repairs a CDR3 against a named germline but does not *name* a missing segment, which
is why `_legacy_fixer/` and `res/` are still live (issue #462, ROADMAP phase 8).

Nothing here runs. `uv run vdjdb build` replaced it; `docs/` specifies what it used to.
