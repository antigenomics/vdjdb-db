---
name: vdjdb-publish
description: Land a proofread VDJdb chunk on a chunk branch - one commit per chunk, on a branch named for that chunk, targeting dev - after finding or creating its PMID issue on antigenomics/vdjdb-db, refreshing registry/records.tsv with vdjdb identity update, and writing the message the chunk-change rule requires (files, row counts, reason, what the build shows, who decided). Use when a chunk is ready to commit, or when git shows new or modified files under chunks/.
---

# vdjdb-publish

Commit each new or changed chunk on its own branch, against its own issue, with a message a reviewer
can check. One chunk at a time, asking before every issue and every commit.

Read [`skills/AUTHORITIES.md`](../AUTHORITIES.md) first. Invariant 1 is what this skill enforces.

## Invocation

```
/vdjdb-publish
```

No arguments. Run from the repository root.

## The branch, before anything else

Gitflow here is `master` → `dev` → branch → `dev` → `master`, and `branch-policy.yml` fails any pull
request into `master` from anything but `dev` or `hotfix/*`.

A chunk branch **starts from `dev` and targets `dev`**. Never branch a chunk off `master`. Never put a
chunk edit on a branch whose subject is the build, the tests, the docs or CI, however mechanical the
edit looks: once every line of a file has changed in a commit nobody reviewed as a data change, a
curation edit and a whitespace edit are indistinguishable in `git log` forever.

```bash
git switch dev && git pull --ff-only
git switch -c chunk/PMID_<id>          # one chunk
git switch -c proofread/<issue>-<slug> # one data issue across several chunks
```

If the working tree already has chunk changes on `dev` or on an unrelated branch, say so and move them
before committing anything.

## Step 1 - collect the changed chunks

```bash
git restore --staged .
git diff --name-only HEAD -- chunks/
git ls-files --others --exclude-standard chunks/
```

Combine, deduplicate, sort. Empty: say "No new or changed chunks found in git" and stop.

## Step 2 - per chunk, in order

Work through the list one file at a time. Do not skip any. Ask before each issue and each commit.

### 2a. The PMID

`PMID_(\d+)\.txt` gives `$pubmedid`. A name that does not match the pattern
(`10xgenomics-2019-07-09.tsv`, `PDB_Database.tsv`) has no PMID: show the filename, say so, and ask
whether to skip it or commit it against a user-supplied issue and message.

### 2b. Check it is ready

```bash
uv run vdjdb qc chunks/PMID_$pubmedid.tsv
uv run vdjdb submission chunks/PMID_$pubmedid.tsv
```

`vdjdb qc` must exit 0. If it does not, stop and run [`/vdjdb-proofread`](../vdjdb-proofread/SKILL.md).
Keep the `submission` output - the commit message needs its numbers.

A chunk that cannot pass goes to `pending/` or `withheld/` instead of onto a branch, with the blocker
written on its issue and the issue left **open**. A branch is invisible: two submissions sat unlanded
for eight and ten years and were found only by checking every unmerged branch against the tracker.

### 2c. Find or create the issue

```bash
gh issue list --repo antigenomics/vdjdb-db --search "PMID:$pubmedid in:title" \
  --state all --json number,title,state,url,body --limit 5
git log --oneline --all -- "chunks/PMID_$pubmedid.tsv" | head -5
```

**If it exists**, show the number, title, state, URL and the first lines of the body. For a modified
tracked file also show the diff: lines added and removed, the row-count delta, and any change in the
column set. Where the new version has cleared metadata the old version had, offer to merge - keep the
old rows and append only rows the new version adds, matched on `cdr3.beta` and `antigen.epitope` - and
do the merge in Python if the user agrees.

Then ask: "Issue #N exists for PMID:$pubmedid. Commit `chunks/PMID_$pubmedid.tsv` with `Fixes #N`?
[y/n/skip]"

**If it does not exist**, fetch the citation:

```bash
curl -s "https://api.ncbi.nlm.nih.gov/lit/ctxp/v1/pubmed/?format=apa&id=$pubmedid"
```

On an error or empty body, fall back to `esummary.fcgi?db=pubmed&id=$pubmedid&retmode=json` and build
the citation from `authors`, `title`, `source` and `pubdate`. Propose title `PMID:$pubmedid` and body
`[<citation>](https://pubmed.ncbi.nlm.nih.gov/$pubmedid/)`, show both, and ask before calling
`gh issue create`.

### 2d. Refresh the registry

**Every chunk branch runs this, and commits the result with the chunk:**

```bash
uv run vdjdb identity update
```

`registry/records.tsv` is the record identity registry, one row per `record_id`, committed. Without
this step it goes stale and the next build retires every record of the chunk it has not seen. The
command is idempotent - on an unchanged corpus it rewrites the file byte for byte - so the diff in the
pull request is exactly the records this branch adds, amends or retires.

Before the registry was committed, landing one 40-record chunk moved `record_id` on 168,723 of 192,753
records.

### 2e. Commit

Stage the chunk and the registry, and nothing else:

```bash
git restore --staged .
git add chunks/PMID_$pubmedid.tsv registry/records.tsv
git status --short
```

The message needs five things, because the only instrument that would otherwise notice a chunk edit is
the comparison against the last release, which reports it as rows appearing and disappearing with no
reason attached:

1. **which files**, and per file how many rows this adds, removes or changes;
2. **why** - the paper, the tracker issue, or the `proofreading/` table the change comes from;
3. **what the build shows** - the row-count delta and the score distribution from step 2b, and the
   entry in `rules/expected_diffs.toml` the difference is declared under;
4. **who decided**, whenever the edit is a curation judgement rather than a mechanical repair. "A
   curator chose X over Y because Z" is the part no diff reconstructs later;
5. `Fixes #$issue_id` on the last line.

A mechanical repair across many files says so and states the invariant that makes it safe. `b0a479d`
is the precedent: 103 files, 137,538 insertions and 137,538 deletions, and the message says the content
is unchanged and names the four header defects it fixed alongside. Without that sentence the diff is
indistinguishable from rewriting the database.

Show the user the full message and the staged file list, then commit. Do not use `--no-verify` if a
hook fails: show the error and wait.

## Step 3 - the release comparison will go red, and that is expected

A new chunk's records are rows the reference release cannot contain, so the `added`, `removed` and
`row_delta` declarations in `rules/expected_diffs.toml` all move at once.

```bash
uv run vdjdb build --out out/
uv run vdjdb diff reference.zip out/legacy --report out/reports/release-diff.md
```

`vdjdb diff --report` prints declared against measured and names the file. Re-measure, update the
block, and **extend that block's `note` with the chunk and its record count** - the `PMID_18025130`
entry is the worked example. Nothing automates this on purpose: a command that re-froze the counts
would turn the gate into a rubber stamp.

## Step 4 - open the pull request

```bash
git push -u origin HEAD
gh pr create --repo antigenomics/vdjdb-db --base dev \
  --title "<chunk>: <what it adds>" --body "<the commit message, plus the submission report>"
```

`chunk-check` is the required check on a chunk pull request. It reports records added and removed, the
score histogram and QC findings by rule, and posts them as a sticky comment. Read it rather than
re-deriving it.

The chunk reaches `master` with the next `dev` → `master` merge, after the full build has run green on
`dev`. Not before, and not by a pull request of its own.

## Step 5 - report

Per chunk: committed (and to which issue), skipped, or quarantined (and where, with the blocker).
Then the branch name and the pull request URL.

## Errors

- `gh` not authenticated: stop, tell the user to run `gh auth login`.
- A `curl` fetch fails: show the error, ask for the citation.
- `git commit` fails on a hook: show the error, wait for guidance.
- `vdjdb identity update` changes more rows than the chunk has records: stop. Either the registry was
  already stale or the branch is not based on current `dev`. Do not commit it.
