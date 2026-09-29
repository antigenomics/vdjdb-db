"""The curation skills, checked against the build they instruct a curator to drive.

A skill is an instruction document an agent loads and follows. Nothing executes it, so nothing fails
when it goes stale, and these six went stale in a way that mattered: they were written against
``py_src/``, which phase 4 deleted, and two of them told a curator to normalise murine MHC names to
``H-2Db``. That is backwards. ``patches/mhc.dict`` declares the opposite, ``H2-`` is the MGI symbol
prefix, and the split those two spellings caused cost 768 records their motif badge on the deployed
site. Following the skill would have undone 3,001 rows of corpus curation.

So every factual claim a skill makes about this repository is asserted here: the paths it names, the
QC rule names it quotes, the commands it tells a curator to run. A skill may still be wrong about
immunology - no test catches that - but it can no longer be wrong about the code.

The structural checks come from the Agent Skills format's own guidance, which the two surveys of the
area agree on (arXiv:2607.25032 §4, arXiv:2608.29596 §3.2-3.4): metadata is always in context so the
description has to be specific enough to route on, the body loads whole on selection so it stays
small, and references are linked one level deep so a partial read cannot lose them.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

from vdjdb.qc.rules import RULES
from vdjdb.qc.runner import ADVISORY

ROOT = Path(__file__).resolve().parents[2]
SKILLS = ROOT / "skills"
SHARED = SKILLS / "AUTHORITIES.md"

#: Body length above which a skill is doing too much in one document. Not a spec limit - the format's
#: own guidance is "a few thousand tokens" - but a ceiling that keeps a skill reviewable and forces
#: long reference material into the files that own it. `vdjdb-proofread` was 822 lines and carried its
#: own copies of the IMGT and MHC tables, both of which had drifted from `proofreading/`.
MAX_BODY_LINES = 300

#: Top-level directories a backticked path may name. Anything else in backticks is prose, a command,
#: a column name or a value, and is not checked as a path.
REPO_DIRS = ("chunks", "patches", "proofreading", "registry", "res", "rules", "skills", "src",
             "tests", "docs", "summary", "pending", "withheld", "attic", ".github")

#: Retired paths. A skill naming one is instructing a curator to open a file that is not there.
RETIRED = ("py_src/", "src/*.groovy", "chunks2/", "vdjdb-motifs/", "summary/pubmed_years.tsv")

PATH_IN_BACKTICKS = re.compile(r"`([^`\s]+/[^`\s]*)`")
MARKDOWN_LINK = re.compile(r"\]\((?!https?:)([^)#]+)")
#: `uv run vdjdb <sub>` or `vdjdb <sub>`, capturing the subcommand.
VDJDB_COMMAND = re.compile(r"\bvdjdb\s+([a-z][a-z-]+)")
FRONTMATTER = re.compile(r"\A---\n(.*?)\n---\n", re.S)


def _skills() -> list[Path]:
    return sorted(SKILLS.glob("*/SKILL.md"))


def _frontmatter(path: Path) -> dict[str, str]:
    m = FRONTMATTER.match(path.read_text())
    assert m, f"{path.relative_to(ROOT)} has no YAML frontmatter"
    out: dict[str, str] = {}
    for line in m.group(1).splitlines():
        if ":" in line and not line.startswith((" ", "\t")):
            k, v = line.split(":", 1)
            out[k.strip()] = v.strip()
    return out


def _body(path: Path) -> str:
    return FRONTMATTER.sub("", path.read_text())


ALL = _skills()
DOCUMENTS = [*ALL, SHARED]


def test_there_are_skills_to_check() -> None:
    """A rename or a move that emptied `skills/` would otherwise make every test below vacuous."""
    assert len(ALL) >= 6, f"found {len(ALL)} skills under {SKILLS}"
    assert SHARED.is_file(), f"{SHARED.relative_to(ROOT)} is what every skill links to"


# ---------------------------------------------------------------------------------------------
# Frontmatter: the part that is always in context, and decides whether the skill is chosen at all
# ---------------------------------------------------------------------------------------------

@pytest.mark.parametrize("path", ALL, ids=lambda p: p.parent.name)
def test_the_name_in_the_frontmatter_is_the_directory_name(path: Path) -> None:
    """The directory name is the slash command. A mismatch means `/name` does not resolve, which is
    how the bodies came to say `/extract` while the command was `/vdjdb-extract`."""
    assert _frontmatter(path)["name"] == path.parent.name


@pytest.mark.parametrize("path", ALL, ids=lambda p: p.parent.name)
def test_the_description_says_what_and_when_in_enough_words_to_route_on(path: Path) -> None:
    """Selection is a match against this one line and nothing else, so a vague one is not chosen and
    an overlapping one is chosen wrongly. It also has to say *when*, not only what."""
    desc = _frontmatter(path)["description"]
    assert 120 <= len(desc) <= 900, f"description is {len(desc)} characters"
    assert re.search(r"\bUse (when|after|for|as)\b", desc), (
        "the description does not say when to use the skill")
    assert not re.match(r"\s*(I |This skill |You )", desc), "write the description in third person"


@pytest.mark.parametrize("path", ALL, ids=lambda p: p.parent.name)
def test_the_body_is_small_enough_to_load_whole(path: Path) -> None:
    n = len(_body(path).splitlines())
    assert n <= MAX_BODY_LINES, (
        f"{n} lines. Move reference material into the file that owns it - `proofreading/imgt.md`, "
        f"`proofreading/mhc.md`, `docs/standards/` - and link to it.")


@pytest.mark.parametrize("path", ALL, ids=lambda p: p.parent.name)
def test_every_skill_links_the_shared_authorities(path: Path) -> None:
    """The invariants and the authority table are in one file so six bodies cannot disagree about
    them, which is precisely what happened on murine MHC."""
    assert "AUTHORITIES.md" in _body(path)


# ---------------------------------------------------------------------------------------------
# Claims about this repository
# ---------------------------------------------------------------------------------------------

def _paths(doc: Path) -> set[str]:
    return {p for p in PATH_IN_BACKTICKS.findall(doc.read_text())
            if p.startswith(REPO_DIRS) and not p.endswith(("/", ":"))
            and "$" not in p and "<" not in p}      # a shell template is not a path claim


@pytest.mark.parametrize("doc", DOCUMENTS, ids=lambda p: p.parent.name)
def test_every_repository_path_a_skill_names_exists(doc: Path) -> None:
    missing = []
    for rel in sorted(_paths(doc)):
        if "*" in rel:
            if not list(ROOT.glob(rel)):
                missing.append(rel)
        elif not (ROOT / rel).exists():
            missing.append(rel)
    assert not missing, f"{doc.relative_to(ROOT)} names paths that do not exist: {missing}"


@pytest.mark.parametrize("doc", DOCUMENTS, ids=lambda p: p.parent.name)
def test_every_relative_link_resolves(doc: Path) -> None:
    """A reference is linked one level deep from the skill so a partial read cannot lose it. A link
    that 404s is the same defect with none of the excuse."""
    broken = [t for t in MARKDOWN_LINK.findall(doc.read_text())
              if not (doc.parent / t.strip()).resolve().exists()]
    assert not broken, f"{doc.relative_to(ROOT)} has broken links: {broken}"


@pytest.mark.parametrize("doc", DOCUMENTS, ids=lambda p: p.parent.name)
def test_no_skill_names_a_retired_path(doc: Path) -> None:
    text = doc.read_text()
    named = [r for r in RETIRED if r in text]
    assert not named, (
        f"{doc.relative_to(ROOT)} still refers to {named}. `py_src/` was the pandas pipeline phase 4 "
        "deleted; the build is `src/vdjdb/` and the CLI is `vdjdb`.")


def _quoted_rules(doc: Path) -> set[str]:
    """Backticked strings shaped like a QC finding name: `bad x.y`, `no.cdr3`, `non-functional v.beta`."""
    text = doc.read_text()
    return {m for m in re.findall(r"`((?:bad|no|non-functional) [a-z0-9.]+|no\.[a-z0-9.]+)`", text)}


@pytest.mark.parametrize("doc", DOCUMENTS, ids=lambda p: p.parent.name)
def test_every_qc_rule_name_a_skill_quotes_is_a_rule(doc: Path) -> None:
    """A skill telling a curator what `bad mhc.a` means, for a rule that has been renamed, sends them
    looking through a report for a string that is not in it."""
    unknown = sorted(_quoted_rules(doc) - set(RULES))
    assert not unknown, (
        f"{doc.relative_to(ROOT)} quotes findings that `vdjdb.qc.rules.RULES` does not define: "
        f"{unknown}")


def test_the_advisory_rules_a_skill_lists_are_the_advisory_ones() -> None:
    """`vdjdb-proofread` tells a curator which findings never fail a build and why. A rule that moved
    between the two lists makes that table wrong in the direction that costs a curator most: fixing
    something the build deliberately tolerates, in `chunks/`, which is the data."""
    body = _body(SKILLS / "vdjdb-proofread" / "SKILL.md")
    fatal, advisory = body.split("**Advisory rules", 1)
    for rule in sorted(set(RULES) & set(ADVISORY)):
        assert f"`{rule}`" in advisory, f"{rule} is advisory and the skill does not list it as one"
    for rule in sorted(set(RULES) - set(ADVISORY)):
        assert f"`{rule}`" in fatal, f"{rule} is fatal and the skill does not list it as one"


def _subcommands() -> set[str]:
    """Every group and command name the CLI exposes, one level deep."""
    from vdjdb.cli import app

    names: set[str] = set()
    for info in app.registered_commands:
        names.add(info.name or info.callback.__name__.replace("_", "-"))
    for group in app.registered_groups:
        names.add(group.name or "")
        inner = group.typer_instance
        names |= {i.name or i.callback.__name__.replace("_", "-")
                  for i in inner.registered_commands}
    return names - {""}


@pytest.mark.parametrize("doc", DOCUMENTS, ids=lambda p: p.parent.name)
def test_every_vdjdb_command_a_skill_shows_exists(doc: Path) -> None:
    """The skills now drive the CLI rather than carrying their own copy of what it does, which is only
    an improvement while the command names are right."""
    known = _subcommands() | {"run"}          # `uv run vdjdb ...`
    shown = set(VDJDB_COMMAND.findall(doc.read_text())) - {"qc"} | (
        {"qc"} if "vdjdb qc" in doc.read_text() else set())
    unknown = sorted(shown - known)
    assert not unknown, f"{doc.relative_to(ROOT)} shows `vdjdb {unknown}`, which the CLI has not"


# ---------------------------------------------------------------------------------------------
# The two claims that were actively wrong
# ---------------------------------------------------------------------------------------------

def test_no_skill_tells_a_curator_to_write_the_classical_murine_prefix() -> None:
    """`H2-` is the MGI gene symbol prefix and is what VDJdb records; `H-2` is the immunology spelling
    and is not a symbol. Two skills said to convert *toward* `H-2`, against `patches/mhc.dict` and
    against 3,001 rows of corpus curation. The corpus has 26 `H-2` cells left, all defects.

    `H-2` may appear as the source of a normalisation. It may not appear as the target, which is what
    this asserts in the two shapes a target takes: after an arrow, and in the second cell of a
    from/to table row.
    """
    arrow = re.compile(r"(?:\u2192|->|\bto\b)\s*`?H-2")
    for doc in DOCUMENTS:
        for n, line in enumerate(doc.read_text().splitlines(), start=1):
            if "H-2" not in line:
                continue
            where = f"{doc.relative_to(ROOT)}:{n}"
            assert not arrow.search(line), f"{where} normalises toward `H-2`:\n  {line.strip()}"
            cells = [c.strip(" `*") for c in line.split("|")]
            assert not (len(cells) >= 4 and "H2-" in cells[1] and "H-2" in cells[2]), (
                f"{where} has `H-2` in the target cell of a from/to row:\n  {line.strip()}")


def test_no_skill_states_the_fixed_junction_anchor_rule() -> None:
    """"Starts with C, ends with F or W" was the retired build's `is_qq_seq_biologically_valid`. Over
    the corpus it calls 481 correct chains broken - `TRAJ35*01` templates `IGFGNVLHC` - and misses 125
    that are not consistent. The anchor is read from the germline of the segment the record names.
    """
    pattern = re.compile(r"end(?:s|ing)?\s+(?:with|in)\s+`?(?:F|Phe)`?(?:/|\s+or\s+)`?(?:W|Trp)`?",
                         re.I)
    for doc in (*DOCUMENTS, ROOT / "proofreading" / "cdr3_repair.md"):
        for n, line in enumerate(doc.read_text().splitlines(), start=1):
            if not pattern.search(line):
                continue
            # Naming the retired rule in order to reject it is the point of several of these lines.
            assert re.search(r"\bno\b|not|never|calls|misses|assumed|retired|only\s+14|fixed", line, re.I), (
                f"{doc.relative_to(ROOT)}:{n} states the fixed anchor rule as if it held:\n"
                f"  {line.strip()}")
