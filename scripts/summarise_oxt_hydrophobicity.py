"""Prints the tables of docs/oxt_hydrophobicity.md from the measurement in OUT_FILE.

Run after scripts/measure_oxt_hydrophobicity.py::

    python scripts/summarise_oxt_hydrophobicity.py

Prints Markdown, to be pasted into the doc.
"""

from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
OUT_FILE = REPO / "docs" / "oxt_hydrophobicity" / "oxt_features.csv"

SCALES = ["mj_scaled", "kd", "cw", "mj"]

# The features shown per structure: the surface mean and sum, the two counts,
# and the two default patch features that read hydrophobicity.
SHOWN = ["SurfMhpMean", "SurfMhpSum", "NSurfPosMhp", "NSurfNegMhp", "SurfMhpPatchArea", "SurfMhpLocalPosSum"]

METADATA = ["structure", "n_oxt", "scale", "variant", "ID"]

# 1GDW_h is 1GDW with hydrogens added. Prodes describes heavy atoms only, so its
# features equal those of 1GDW to the last digit, and counting it would count
# one structure twice.
DUPLICATES = ["1GDW_h"]


def variant(table, scale, name):
    """Returns the features of one variant on one scale, indexed by structure."""

    rows = table[(table.scale == scale) & (table.variant == name)]

    return rows.set_index("structure").drop(columns=["n_oxt", "scale", "variant", "ID"])


def moved_columns(table):
    """Prints, per scale, how many columns the OXT value moves and whether any is outside the Mhp block."""

    print("| scale | columns that move | of which not hydrophobicity | structures that move |")
    print("| --- | --- | --- | --- |")
    for scale in SCALES:
        old, new = variant(table, scale, "old"), variant(table, scale, "new")
        numeric = [column for column in old.columns if pd.api.types.is_numeric_dtype(old[column])]
        changed = (new[numeric] != old[numeric]) & ~(new[numeric].isna() & old[numeric].isna())
        columns = [column for column in numeric if changed[column].any()]
        structures = [structure for structure in changed.index if changed.loc[structure].any()]
        outside = [column for column in columns if "Mhp" not in column]
        print(f"| `{scale}` | {len(columns)} | {len(outside)} | {len(structures)} of {len(changed.index)} |")
    print()


def per_structure(table, scale):
    """Prints, for the structures with an OXT, the shown features without OXT, with the old value and with the new."""

    old, new, no_oxt = variant(table, scale, "old"), variant(table, scale, "new"), variant(table, scale, "no_oxt")

    print(f"`{scale}`, each cell no OXT / old / new:")
    print()
    print("| structure | " + " | ".join(f"`{feature}`" for feature in SHOWN) + " |")
    print("| --- |" + " --- |" * len(SHOWN))
    for structure in no_oxt.index:
        cells = [f"{no_oxt.at[structure, f]:.4g} / {old.at[structure, f]:.4g} / {new.at[structure, f]:.4g}" for f in SHOWN]
        print(f"| {structure} | " + " | ".join(cells) + " |")
    print()


def shift_from_no_oxt(table):
    """Prints, per scale and feature, the direction of the shift an OXT causes under each value."""

    print("| scale | feature | OXT raises it, old | OXT lowers it, old | OXT raises it, new | OXT lowers it, new | shift smaller with new |")
    print("| --- | --- | --- | --- | --- | --- | --- |")
    for scale in SCALES:
        old, new, no_oxt = variant(table, scale, "old"), variant(table, scale, "new"), variant(table, scale, "no_oxt")
        old, new = old.loc[no_oxt.index], new.loc[no_oxt.index]
        for feature in SHOWN:
            old_shift, new_shift = old[feature] - no_oxt[feature], new[feature] - no_oxt[feature]
            if old_shift.isna().all():
                continue
            smaller = (new_shift.abs() < old_shift.abs()).sum()
            print(
                f"| `{scale}` | `{feature}` | {(old_shift > 0).sum()} | {(old_shift < 0).sum()} | "
                f"{(new_shift > 0).sum()} | {(new_shift < 0).sum()} | {smaller} of {len(no_oxt.index)} |"
            )
    print()


def relative_change(table):
    """Prints, per scale and feature, the median and largest change from the old value to the new, in per cent."""

    print("| scale | feature | median change, % | largest change, % |")
    print("| --- | --- | --- | --- |")
    for scale in SCALES:
        old, new, no_oxt = variant(table, scale, "old"), variant(table, scale, "new"), variant(table, scale, "no_oxt")
        old, new = old.loc[no_oxt.index], new.loc[no_oxt.index]
        for feature in SHOWN:
            change = ((new[feature] - old[feature]) / old[feature].abs() * 100).dropna()
            if change.empty:
                continue
            largest = change.loc[change.abs().idxmax()]
            print(f"| `{scale}` | `{feature}` | {change.median():+.2f} | {largest:+.2f} |")
    print()


def main():
    """Prints every table of the doc."""

    table = pd.read_csv(OUT_FILE)
    table = table[~table.structure.isin(DUPLICATES)]

    moved_columns(table)
    relative_change(table)
    for scale in SCALES:
        per_structure(table, scale)
    shift_from_no_oxt(table)


if __name__ == "__main__":
    main()
