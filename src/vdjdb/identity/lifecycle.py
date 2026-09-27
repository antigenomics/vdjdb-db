"""What happened to the id you had.

An id that vanishes is the failure a consumer cannot diagnose: a reference that used to resolve
returns nothing, and nothing distinguishes a correction from a deletion. So every level keeps a
lifecycle row, and the row holds only what a build cannot recompute (``ROADMAP.md`` section 10.4).

Two properties shape this module.

**The lifecycle never decides what id anything gets.** The active set and its keys are recomputed
from ``chunks/`` on every build (hard rule 9), and :mod:`vdjdb.identity.levels` derives each id from
its own key. This file answers one question after the fact, and a build with no previous lifecycle
file still produces exactly the same ids.

**Retirement is per release, not per build.** A curation branch may add a clonotype and remove it
again before anything ships, and neither event is a lifecycle event: only the release job calls
:func:`advance`. That is also what keeps a curation pull request readable, since nothing here moves
when a chunk changes.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import polars as pl

from .levels import LEVELS, level_of

#: One row per id ever seen. No key columns: a key is recomputable and an id's history is not.
LIFECYCLE_COLUMNS: tuple[str, ...] = (
    "id", "level", "state", "first_release", "last_release", "replaced_by",
)

ACTIVE = "active"
RETIRED = "retired"

_SCHEMA = dict.fromkeys(LIFECYCLE_COLUMNS, pl.Utf8)


def empty() -> pl.DataFrame:
    """An empty lifecycle table with the right schema, which is what a first build reads."""
    return pl.DataFrame(schema=_SCHEMA)


def read(path: Path | None) -> pl.DataFrame:
    """A lifecycle table, or an empty one when the path is absent.

    Missing is not an error. A fork and a first build have no history, and they are expected to run:
    the ids come out the same and only the history is unknown.
    """
    if path is None or not Path(path).exists():
        return empty()
    # `empty_string_is_null=False`: an empty `replaced_by` means "no replacement", and a null would be
    # the pandas three-way ambiguity re-entering through the file (hard rule 6).
    return pl.read_csv(path, separator="\t", schema=_SCHEMA, empty_string_is_null=False)


def write(df: pl.DataFrame, path: Path) -> None:
    """Sorted by ``(level, id)``, so a diff between two releases is reviewable."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.select(list(LIFECYCLE_COLUMNS)).sort("level", "id").write_csv(
        path, separator="\t", quote_style="never")


def present(tables: dict[str, pl.DataFrame]) -> pl.DataFrame:
    """``(id, level)`` for every derived id in this build, one row per distinct id.

    Reads whichever of ``records`` and ``chains`` it is given, so it works on a fixture as well as on
    a full build. A level whose column is absent contributes nothing rather than raising, because
    ``vdjdb identity`` has to run against a partial build directory.
    """
    rows = []
    for level in LEVELS:
        for frame in tables.values():
            if level.column not in frame.columns:
                continue
            ids = (frame.select(pl.col(level.column).alias("id"))
                        .filter(pl.col("id").is_not_null() & (pl.col("id") != ""))
                        .unique())
            rows.append(ids.with_columns(pl.lit(level.name).alias("level")))
    if not rows:
        return pl.DataFrame(schema={"id": pl.Utf8, "level": pl.Utf8})
    return pl.concat(rows, how="vertical").unique().sort("level", "id")


@dataclass(frozen=True, slots=True)
class LifecycleReport:
    """Per level, how the id set moved between two releases."""

    release: str
    added: pl.DataFrame
    returned: pl.DataFrame
    retired: pl.DataFrame

    def counts(self) -> pl.DataFrame:
        """One row per level: added, returned, retired. Zero-filled, so every level appears."""
        parts = []
        for name, frame in (("added", self.added), ("returned", self.returned),
                            ("retired", self.retired)):
            parts.append(frame.group_by("level").len().rename({"len": name}))
        out = pl.DataFrame({"level": [lv.name for lv in LEVELS]})
        for frame in parts:
            out = out.join(frame, on="level", how="left")
        return out.with_columns(pl.col("added", "returned", "retired").fill_null(0))

    def as_markdown(self) -> str:
        lines = [f"# Identity lifecycle, release {self.release}", "",
                 "| Level | Added | Returned | Retired |", "|---|---|---|---|"]
        for row in self.counts().iter_rows(named=True):
            lines.append(f"| {row['level']} | {row['added']} | {row['returned']} "
                         f"| {row['retired']} |")
        return "\n".join(lines) + "\n"


def compare(previous: pl.DataFrame, now: pl.DataFrame, *, release: str) -> LifecycleReport:
    """Added, returned and retired, without writing anything.

    ``returned`` is separate from ``added`` because the two mean different things to a consumer: an
    id that comes back resolves to the key it always had, since a derived id *is* the hash of its
    key, while an added one has never been published.
    """
    prev_active = previous.filter(pl.col("state") == ACTIVE).select("id", "level")
    prev_retired = previous.filter(pl.col("state") == RETIRED).select("id", "level")
    known = previous.select("id", "level")
    return LifecycleReport(
        release=release,
        added=now.join(known, on=["id", "level"], how="anti").sort("level", "id"),
        returned=now.join(prev_retired, on=["id", "level"], how="semi").sort("level", "id"),
        retired=prev_active.join(now, on=["id", "level"], how="anti").sort("level", "id"),
    )


def advance(previous: pl.DataFrame, now: pl.DataFrame, *, release: str,
            replaced_by: dict[str, str] | None = None) -> pl.DataFrame:
    """The lifecycle table for ``release``.

    An id in this build is active with ``last_release`` moved to ``release``, keeping the
    ``first_release`` it already had. An id in ``previous`` and not in this build is retired with its
    ``last_release`` left where it was, because that is the release that last carried it.

    ``replaced_by`` names the id that took over, for the retirements where one did. It is supplied by
    the caller rather than guessed: for ``record_id`` the registry's amendment step knows the answer,
    and for a derived level a replacement is a curation judgement that no diff can infer.
    """
    repl = replaced_by or {}
    kept = (now.join(previous, on=["id", "level"], how="left")
               .with_columns(
                   pl.lit(ACTIVE).alias("state"),
                   pl.col("first_release").fill_null(release),
                   pl.lit(release).alias("last_release"),
                   pl.lit("").alias("replaced_by")))
    gone = (previous.filter(pl.col("state") == ACTIVE)
                    .join(now, on=["id", "level"], how="anti")
                    .with_columns(pl.lit(RETIRED).alias("state"),
                                  pl.col("replaced_by").fill_null("")))
    if repl:
        gone = gone.with_columns(
            pl.col("id").replace_strict(repl, default=None)
              .fill_null(pl.col("replaced_by")).alias("replaced_by"))
    # Rows already retired in `previous` and still absent keep their history untouched: a retirement
    # is a fact about a release that has shipped, so a later build must not rewrite it.
    frozen = (previous.filter(pl.col("state") == RETIRED)
                      .join(now, on=["id", "level"], how="anti"))
    out = pl.concat([kept.select(list(LIFECYCLE_COLUMNS)),
                     gone.select(list(LIFECYCLE_COLUMNS)),
                     frozen.select(list(LIFECYCLE_COLUMNS))], how="vertical")
    return out.sort("level", "id")


def resolve(lifecycle: pl.DataFrame, identifier: str) -> dict[str, str] | None:
    """One id's row, or ``None``. The level is read from the prefix, so a typo says so."""
    hit = lifecycle.filter(pl.col("id") == identifier)
    if hit.is_empty():
        return None
    return hit.row(0, named=True)


def unknown_prefixes(now: pl.DataFrame) -> list[str]:
    """Ids whose prefix belongs to no level. Empty is the passing state.

    Catches an id column that was populated by something other than
    :mod:`vdjdb.identity.levels`, which is the way a hand-written fixture drifts from the schema.
    """
    return [i for i in now["id"].to_list() if level_of(i) is None]

