"""Writes the full feature set of every shipped structure, for a before and after comparison.

Deliberately not a pytest test, for the same reason as scripts/compare_features.py:
it takes minutes, and the point of it is to be run twice on two different
commits::

    python scripts/write_features.py before        # on the commit before
    python scripts/write_features.py after         # on the commit after
    python scripts/compare_features.py before after

Run with the full 105 feature set and with screening off, so that the comparison
covers every column prodes can write and does not depend on the ionic strength
default.
"""

import sys
from pathlib import Path

import prodes

STRUCTURES = Path("tests/data")


def main():
    """Writes one bundle per shipped structure into the directory named on the command line."""

    if len(sys.argv) != 2:
        raise SystemExit("usage: python scripts/write_features.py <output directory>")

    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)

    for structure in sorted(STRUCTURES.glob("*.pdb.zip")):
        name = structure.name.split(".")[0]
        prodes.run_prodes(str(structure), str(out / f"{name}.zip"), full_features=True, ionic_strength_molar=0)
        print("wrote", name, flush=True)


if __name__ == "__main__":
    main()
