"""Compares the features of two sets of prodes output bundles, column by column.

Deliberately not a pytest test. Computing the full feature set for every
structure in tests/data takes minutes, which is too slow for CI, and the run it
compares against has to be made on a different commit. It is meant to be started
by hand by someone changing the way structures are read or described::

    python scripts/write_features.py before        # on the commit before
    python scripts/write_features.py after         # on the commit after
    python scripts/compare_features.py before after

The regression test in tests/test_sasa.py compares one structure against a
committed reference within a tolerance, which is what catches a change in the
science. This catches a change anywhere, on every shipped structure, at the last
written digit, which is what a change that is supposed to move nothing needs.

Exits non-zero if anything differs, so it can be used as a gate.
"""

import sys
from pathlib import Path

import pandas as pd

from prodes.output import read_features


def differing_columns(before, after):
    """Returns one description per column whose value changed, comparing exactly."""

    differences = []
    for column in before.columns:
        left, right = before[column].iloc[0], after[column].iloc[0]

        if isinstance(left, str) or isinstance(right, str):
            if left != right:
                differences.append(f"{column}: {left!r} -> {right!r}")
            continue

        # A feature that is absent in both runs is not a difference. Comparing
        # NaN with == would call every one of them one.
        if pd.isna(left) and pd.isna(right):
            continue

        if left != right:
            differences.append(f"{column}: {left!r} -> {right!r}")

    return differences


def compare(before_directory, after_directory):
    """Compares every bundle in one directory against its namesake in the other.

    Returns the descriptions of everything that differs, empty if nothing does.
    """

    failures = []
    for bundle in sorted(Path(before_directory).glob("*.zip")):
        counterpart = Path(after_directory) / bundle.name
        if not counterpart.exists():
            failures.append(f"{bundle.name}: no counterpart in {after_directory}")
            continue

        before, after = read_features(bundle), read_features(counterpart)
        if list(before.columns) != list(after.columns):
            failures.append(f"{bundle.name}: the columns themselves differ")
            continue

        differences = differing_columns(before, after)
        print(f"{bundle.stem:16s} {len(before.columns):4d} columns, {len(differences):3d} differing")
        for description in differences[:10]:
            print("   ", description)
        if differences:
            failures.append(f"{bundle.name}: {len(differences)} columns differ")

    return failures


def main():
    """Compares the two directories named on the command line."""

    if len(sys.argv) != 3:
        raise SystemExit("usage: python scripts/compare_features.py <before directory> <after directory>")

    failures = compare(sys.argv[1], sys.argv[2])

    print()
    if failures:
        print("DIFFERENT:")
        for failure in failures:
            print("  ", failure)
        raise SystemExit(1)

    print("identical in every column of every structure")


if __name__ == "__main__":
    main()
