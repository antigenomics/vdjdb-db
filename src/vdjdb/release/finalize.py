"""Commit only the release pointer on current master, preserving concurrent source changes."""
from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path

from .bundle import prepare_latest, version_of


def finalize_latest(repo: Path, tag: str, *, attempts: int = 3) -> None:
    """Fetch master and retry normal pushes after races; never move the pinned release tag.

    Each attempt has a fresh private worktree. Commit failures propagate immediately; only
    failed pushes are retried. A rerun is idempotent when the pointer is already current.
    """
    version_of(tag)
    if attempts < 1:
        raise ValueError("attempts must be positive")
    repo = repo.resolve()
    for _ in range(attempts):
        subprocess.run(["git", "fetch", "origin", "refs/heads/master"], cwd=repo, check=True)
        with tempfile.TemporaryDirectory(prefix="vdjdb-finalize-") as temporary:
            checkout = Path(temporary) / "checkout"
            subprocess.run(["git", "worktree", "add", "--detach", str(checkout), "FETCH_HEAD"],
                           cwd=repo, check=True)
            try:
                prepare_latest(tag, checkout / "latest-version.txt")
                subprocess.run(["git", "add", "latest-version.txt"], cwd=checkout, check=True)
                diff = subprocess.run(["git", "diff", "--cached", "--quiet"], cwd=checkout)
                if diff.returncode == 0:
                    return
                if diff.returncode != 1:
                    diff.check_returncode()
                subprocess.run(["git", "commit", "-m", f"latest-version.txt -> {tag}"],
                               cwd=checkout, check=True)
                push = subprocess.run(["git", "push", "origin", "HEAD:refs/heads/master"],
                                      cwd=checkout)
                if push.returncode == 0:
                    return
            finally:
                subprocess.run(["git", "worktree", "remove", "--force", str(checkout)],
                               cwd=repo, check=True)
    raise RuntimeError(f"could not finalize latest-version.txt after {attempts} normal pushes")
