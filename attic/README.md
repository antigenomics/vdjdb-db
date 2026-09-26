# attic — retired, kept because something still depends on it

The Groovy implementation of the database build, retired in phase 4. It is **not** dead weight and
must not be deleted: `BuildDatabase.groovy` is the one correct specification for the two `*.meta.txt`
files and for deriving the `vdjdb.txt` header from the metadata, and `tests/unit/test_schema.py`
**parses it at test time** rather than copying it. A copy would drift exactly the way the nine
duplicated column lists did.

Nothing here runs. `uv run vdjdb build` replaced it; `docs/` specifies what it used to.
