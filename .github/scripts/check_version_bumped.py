"""Fail unless this branch declares a version above the base branch and above every tag.

A version is bumped early in the work, while the base branch is still where it was. By the
time the branch merges, the base branch has usually moved, and a rebase onto a base branch
that has itself been released leaves two trees declaring one version: git raises no conflict,
because both sides hold the same text, and the tag job on the base branch then declines to
create a tag that already exists. The released tag points at one tree and the merged code sits
in another, and nothing says so out loud.

This check states the invariant instead: the version this branch declares must be strictly
greater than the version the base branch declares, and strictly greater than every version
already tagged. Run on every pull request, so a rebase that overtakes a release is a red check
rather than a silent duplicate.

The version is read from `pyproject.toml`, following setuptools' dynamic `attr` directive to
whichever module actually holds it, so this script can be copied into another repository
without editing.

Reads git refs that are already in the repository; pass `--fetch` to update them first.

    python .github/scripts/check_version_bumped.py --base main --fetch
"""

from __future__ import annotations

import argparse
import ast
import re
import subprocess
import sys
from collections.abc import Callable
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

# A bare version, no v prefix, matching every tag this repository has and the pattern the tag
# job in .github/workflows/tests.yml enforces before it tags anything.
VERSION_PATTERN = re.compile(r"^(\d+)\.(\d+)\.(\d+)(?:[-.]([0-9A-Za-z.]+))?$")

# Returns the text of a repository-relative path, from the working tree or from a git ref.
ReadPath = Callable[[str], str]


class CheckFailed(Exception):
    """A version invariant does not hold, or the version could not be read."""


def git(*args: str) -> str:
    """Run git in the repository and return its stdout, raising CheckFailed on failure."""
    result = subprocess.run(
        ["git", *args],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise CheckFailed(f"git {' '.join(args)} failed: {result.stderr.strip()}")
    return result.stdout


def read_from_worktree(path: str) -> str:
    """Return the text of a repository-relative path in the working tree."""
    full_path = REPO_ROOT / path
    if not full_path.is_file():
        raise CheckFailed(f"{path} does not exist in the working tree.")
    return full_path.read_text()


def read_from_ref(ref: str) -> ReadPath:
    """Return a reader for repository-relative paths as they are at a git ref."""

    def read(path: str) -> str:
        result = subprocess.run(
            ["git", "show", f"{ref}:{path}"],
            cwd=REPO_ROOT,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise CheckFailed(f"{path} could not be read at {ref}: {result.stderr.strip()}")
        return result.stdout

    return read


def version_module_candidates(attr: str) -> list[tuple[str, str]]:
    """Return (path, attribute) candidates for a setuptools dynamic `attr` directive.

    `attr = "prodes.version.VERSION"` names the attribute VERSION in the module
    prodes.version, which may sit at prodes/version.py, at prodes/version/__init__.py,
    or under a src/ layout.
    """
    module_path, _, attribute = attr.rpartition(".")
    if not module_path or not attribute:
        raise CheckFailed(f'pyproject.toml declares version = {{attr = "{attr}"}}, which names no attribute.')
    parts = module_path.split(".")
    stems = ["/".join(parts), "src/" + "/".join(parts)]
    return [(f"{stem}{suffix}", attribute) for stem in stems for suffix in (".py", "/__init__.py")]


def assignment_value(text: str, attribute: str, path: str) -> str:
    """Return the string a module assigns to an attribute, read without importing it."""
    for node in ast.parse(text).body:
        if isinstance(node, ast.Assign) and any(getattr(target, "id", None) == attribute for target in node.targets):
            if isinstance(node.value, ast.Constant) and isinstance(node.value.value, str):
                return node.value.value
            raise CheckFailed(f"{path} assigns {attribute} something other than a string literal.")
    raise CheckFailed(f"{path} does not assign {attribute}.")


def declared_version(read_path: ReadPath) -> str:
    """Return the version a tree declares, wherever pyproject.toml says it is kept.

    pyproject.toml is parsed with regular expressions rather than tomllib, which is only in the
    standard library from Python 3.11, so this script runs under whatever Python is to hand.
    """
    pyproject = read_path("pyproject.toml")
    dynamic = re.search(r'^\s*version\s*=\s*\{\s*attr\s*=\s*"([^"]+)"', pyproject, re.MULTILINE)
    if dynamic:
        failures = []
        for path, attribute in version_module_candidates(dynamic.group(1)):
            try:
                return assignment_value(read_path(path), attribute, path)
            except CheckFailed as failure:
                failures.append(str(failure))
        raise CheckFailed("pyproject.toml declares version = {attr = " f'"{dynamic.group(1)}"}} and the version could not be read: ' + "; ".join(failures))
    static = re.search(r'^\s*version\s*=\s*"([^"]+)"', pyproject, re.MULTILINE)
    if static:
        return static.group(1)
    raise CheckFailed("pyproject.toml declares no version, statically or dynamically.")


def sort_key(version: str) -> tuple[int, int, int, int, str]:
    """Return an ordering key for a version, rejecting anything that is not a version number.

    A pre-release suffix sorts below the release it precedes, so 4.0.0-rc1 < 4.0.0, matching
    https://semver.org/. Suffixes are otherwise compared as text, which orders rc1 < rc2 and is
    all this check needs.
    """
    match = VERSION_PATTERN.fullmatch(version)
    if not match:
        raise CheckFailed(f"'{version}' is not a MAJOR.MINOR.PATCH version number. " "Tags in this repository are bare versions with no v prefix.")
    major, minor, patch, suffix = match.groups()
    return (int(major), int(minor), int(patch), 0 if suffix else 1, suffix or "")


def tagged_versions() -> list[str]:
    """Return every tag in the repository that is a version number."""
    return [tag for tag in git("tag", "--list").split() if VERSION_PATTERN.fullmatch(tag)]


def check(base_ref: str) -> list[str]:
    """Return the reasons the declared version is not acceptable, empty if it is."""
    version = declared_version(read_from_worktree)
    sort_key(version)  # Rejects a version that is not a version number at all.
    reasons = []

    base_version = declared_version(read_from_ref(f"origin/{base_ref}"))
    if sort_key(version) <= sort_key(base_version):
        reasons.append(
            f"This branch declares version {version} and origin/{base_ref} declares "
            f"{base_version}. Every branch has to raise the version, and a rebase onto a "
            f"base branch that has moved needs the version raised again: git does not report "
            f"a conflict when both sides declare the same version. Bump it to above "
            f"{base_version}."
        )

    tags = tagged_versions()
    if tags:
        highest_tag = max(tags, key=sort_key)
        if sort_key(version) <= sort_key(highest_tag):
            reasons.append(
                f"This branch declares version {version} and {highest_tag} is already tagged. "
                f"A released version cannot be declared a second time by different code: the "
                f"tag would keep pointing at the tree that was released and this one would "
                f"never get a tag. Bump it to above {highest_tag}."
            )
    return reasons


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base",
        default="main",
        help="the branch this work will merge into (default: main).",
    )
    parser.add_argument(
        "--fetch",
        action="store_true",
        help="fetch the base branch and the tags from origin first, rather than comparing " "against whatever the repository last saw.",
    )
    args = parser.parse_args()

    try:
        current_branch = git("rev-parse", "--abbrev-ref", "HEAD").strip()
        if current_branch == args.base:
            print(f"On {args.base} itself; there is no base branch to compare against.")
            return 0

        if args.fetch:
            git("fetch", "--quiet", "--tags", "origin", args.base)

        reasons = check(args.base)
    except CheckFailed as failure:
        print(str(failure), file=sys.stderr)
        return 1

    if reasons:
        for reason in reasons:
            print(reason, file=sys.stderr)
        return 1

    print(f"Version {declared_version(read_from_worktree)} is above origin/{args.base} and above every tag.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
