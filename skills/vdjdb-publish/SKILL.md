---
name: vdjdb-publish
description: For each new or changed chunk in chunks/ (by git), find or create a GitHub issue for its PMID, then commit the chunk with "Closes #issue_id". Processes one chunk at a time, always asking user before creating issues or committing.
---

# /vdjdb-publish — Publish VDJdb Chunks to GitHub

## Purpose

Walk through every new or modified file in `chunks/` according to git, and for each one:
1. Determine the PubMed ID from the filename.
2. Find or create the matching GitHub issue (`PMID:$pubmedid`).
3. Commit only that chunk with `Closes #$issue_id`.

The skill processes chunks **one at a time**, always asking the user before creating an issue or committing.

## Invocation

```
/vdjdb-publish
```

No arguments. Run from the root of `vdjdb-db`.

---

## Step-by-step procedure

### 1. Collect changed/added chunks

```bash
git diff --name-only HEAD -- chunks/
git ls-files --others --exclude-standard chunks/
```

Combine both lists (modified tracked files + untracked new files). Deduplicate and sort.

**If the list is empty**, inform the user: "No new or changed chunks found in git." and stop.

**Important**: Before starting, unstage everything so you are working from a clean index:
```bash
git restore --staged .
```

---

### 2. For each chunk — one at a time

Work through the list sequentially. **Do not skip any file. Always ask the user before each commit.**

#### 2a. Extract PubMed ID

- If the filename matches `PMID_(\d+)\.txt`, extract the numeric ID as `$pubmedid`.
- If the filename does **not** match the `PMID_` pattern (e.g. `10xgenomics-2019-07-09.txt`, `PDB_Database.txt`), inform the user that this chunk has no PMID, show the filename, and ask how to proceed. Options: skip it, or commit it manually with a user-supplied message. Then move on.

#### 2b. Check for an existing GitHub issue

Search GitHub issues for the title `PMID:$pubmedid`:

```bash
gh issue list --repo antigenomics/vdjdb-db --search "PMID:$pubmedid in:title" --state all --json number,title,state,url,body --limit 5
```

**Also** check git log for any prior commits that reference this file (useful when a file was already tracked and modified):

```bash
git log --oneline --all -- "chunks/PMID_$pubmedid.txt" | head -5
```

#### 2c-A. Issue already exists

Display to the user:
- Issue number and title
- Issue state (open/closed)
- Issue URL
- First ~200 chars of the body

Ask: "Issue #$number already exists for PMID:$pubmedid. Do you want to commit `chunks/PMID_$pubmedid.txt` with `Closes #$number`? [y/n/skip]"

- **y**: proceed to step 3 (commit).
- **n / skip**: move to the next chunk without committing.

#### 2c-B. Issue does not exist

Fetch the citation from PubMed via the NCBI API:

```bash
curl -s "https://api.ncbi.nlm.nih.gov/lit/ctxp/v1/pubmed/?format=apa&id=$pubmedid"
```

Store the result as `$pubmedid_citation`. If the API returns an error or empty body, fall back to:

```bash
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=$pubmedid&retmode=json" | python3 -c "
import sys, json
d = json.load(sys.stdin)
r = d['result'][$pubmedid if isinstance(list(d['result'].keys())[0], str) else str($pubmedid)]
authors = ', '.join(a['name'] for a in r.get('authors', [])[:3])
title = r.get('title','')
source = r.get('source','')
year = r.get('pubdate','')[:4]
print(f'{authors} ({year}). {title} {source}.')
"
```

Construct the proposed issue:
- **Title**: `PMID:$pubmedid`
- **Body**: `[$pubmedid_citation](https://pubmed.ncbi.nlm.nih.gov/$pubmedid/)`

Show the user the proposed title and body, then ask:
"OK to create this issue on antigenomics/vdjdb-db? [y/n/skip]"

- **y**: create the issue and capture its number:
  ```bash
  gh issue create --repo antigenomics/vdjdb-db \
    --title "PMID:$pubmedid" \
    --body "[$pubmedid_citation](https://pubmed.ncbi.nlm.nih.gov/$pubmedid/)"
  ```
  Capture `$issue_id` from the output URL (the number at the end).
  Then proceed to step 3 (commit).
- **n**: move to the next chunk without committing.
- **skip**: same as n.

---

### 3. Stage and commit the single chunk

Make sure only this one chunk is staged:

```bash
git restore --staged .          # unstage everything first
git add "chunks/PMID_$pubmedid.txt"
```

Confirm to the user what will be committed (one line: filename + issue number + commit message), then:

```bash
git commit -m "Closes #$issue_id"
```

After a successful commit, move on to the next chunk in the list.

---

### 4. Finish

After all chunks have been processed, report a summary:
- How many chunks were committed (and to which issues).
- How many were skipped.

---

## Error handling

- If `gh` is not authenticated, stop immediately and tell the user to run `gh auth login`.
- If a `curl` fetch fails, show the error and ask the user to supply the citation manually before proceeding.
- If `git commit` fails (e.g. pre-commit hook), show the error and wait for user guidance — do **not** use `--no-verify`.
