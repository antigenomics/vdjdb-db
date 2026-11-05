"""Local Git remotes reproduce finalization after concurrent master updates."""
import subprocess
from pathlib import Path

import pytest

from vdjdb.release.bundle import legacy_url
from vdjdb.release.finalize import finalize_latest


def git(path, *args):
    return subprocess.run(["git", *args], cwd=path, check=True, text=True,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE).stdout.strip()


@pytest.fixture
def history(tmp_path):
    remote = tmp_path / "remote.git"
    remote.mkdir()
    git(remote, "init", "--bare", "--initial-branch=master")
    pinned = tmp_path / "pinned"
    git(tmp_path, "clone", str(remote), str(pinned))
    git(pinned, "config", "user.name", "Test")
    git(pinned, "config", "user.email", "test@example.invalid")
    (pinned / "latest-version.txt").write_text("old-url\n")
    git(pinned, "add", ".")
    git(pinned, "commit", "-m", "base")
    git(pinned, "push", "origin", "master")
    source = git(pinned, "rev-parse", "HEAD")
    git(pinned, "tag", "v2026.11.0")
    git(pinned, "checkout", "--detach", source)
    other = tmp_path / "other"
    git(tmp_path, "clone", str(remote), str(other))
    git(other, "config", "user.name", "Test")
    git(other, "config", "user.email", "test@example.invalid")
    (other / "concurrent.txt").write_text("preserve this source change\n")
    git(other, "add", ".")
    git(other, "commit", "-m", "concurrent source change")
    git(other, "push", "origin", "master")
    return pinned, remote, source


def test_finalize_preserves_advanced_master_and_pinned_source_and_is_repeatable(history):
    pinned, remote, source = history
    finalize_latest(pinned, "v2026.11.0")
    assert git(remote, "show", "master:concurrent.txt") == "preserve this source change"
    assert git(remote, "show", "master:latest-version.txt").splitlines()[0] == legacy_url("v2026.11.0")
    assert git(remote, "diff-tree", "--no-commit-id", "--name-only", "-r", "master") == "latest-version.txt"
    assert git(pinned, "rev-parse", "HEAD") == source
    assert git(pinned, "rev-parse", "v2026.11.0") == source
    first = git(remote, "rev-parse", "master")
    finalize_latest(pinned, "v2026.11.0")
    assert git(remote, "rev-parse", "master") == first


def test_commit_failures_are_not_masked(history):
    pinned, remote, _ = history
    hook = pinned / ".git/hooks/pre-commit"
    hook.write_text("#!/bin/sh\nexit 1\n")
    hook.chmod(0o755)
    before = git(remote, "rev-parse", "master")
    with pytest.raises(subprocess.CalledProcessError):
        finalize_latest(pinned, "v2026.11.0")
    assert git(remote, "rev-parse", "master") == before


def test_finalize_retries_a_normal_push_race(history, monkeypatch):
    pinned, remote, _ = history
    other = pinned.parent / "other"
    original_run = subprocess.run
    raced = False

    def run(command, *args, **kwargs):
        nonlocal raced
        if command[:3] == ["git", "push", "origin"] and Path(kwargs.get("cwd", ".")).name == "checkout" and not raced:
            raced = True
            (other / "second.txt").write_text("arrived during finalization")
            git(other, "add", "second.txt")
            git(other, "commit", "-m", "second concurrent source change")
            git(other, "push", "origin", "master")
        return original_run(command, *args, **kwargs)

    monkeypatch.setattr(subprocess, "run", run)
    finalize_latest(pinned, "v2026.11.0")
    assert raced
    assert git(remote, "show", "master:second.txt") == "arrived during finalization"
    assert git(remote, "show", "master:latest-version.txt").splitlines()[0] == legacy_url("v2026.11.0")
