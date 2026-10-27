"""Measures how the hydrophobicity given to the C-terminal oxygen OXT moves the features.

Deliberately not a pytest test: it runs the full feature set of every shipped
structure on four scales, which takes tens of minutes. Run by hand::

    python scripts/measure_oxt_hydrophobicity.py

Each structure in tests/data is described three ways on each scale:

    old       as deposited, OXT given 1.0, which is what versions before 12.0 did
    new       as deposited, OXT given prodes.data.oxt_hydrophobicity
    no_oxt    with every OXT atom removed, where the OXT value cannot matter

The old value is reproduced by replacing prodes.data.oxt_hydrophobicity for the
run, which is the only place the value comes from, so the two runs differ in that
one number and in nothing else. The no_oxt runs measure what a structure source
that does not write OXT gets, against which both values can be compared.

Writes one row per structure, scale and variant to OUT_FILE, appended as soon as
it is measured. A run that is stopped can be started again with the same command:
every structure, scale and variant already in OUT_FILE is skipped. Delete OUT_FILE
to measure everything again. To leave the machine usable while it runs::

    PRODES_N_WORKERS=1 nice -n 19 python scripts/measure_oxt_hydrophobicity.py
"""

import tempfile
import zipfile
from pathlib import Path
from unittest import mock

import pandas as pd

from prodes.calculations.patches import PatchSettings
from prodes.data import hydrophobic_scale
from prodes.output import read_features
from prodes.run import calculate

REPO = Path(__file__).resolve().parent.parent
STRUCTURES = REPO / "tests" / "data"
OUT_FILE = REPO / "docs" / "oxt_hydrophobicity" / "oxt_features.csv"

SCALES = ["mj_scaled", "kd", "cw", "mj"]

# The value versions before 12.0 gave OXT on every scale.
OLD_OXT_VALUE = 1.0


def patch_settings(scale):
    """Returns the patch settings for this scale.

    mj has no negative values, so the default threshold of 0 is refused for it.
    The mean of the scale is used instead, which splits its residues roughly in
    half the way 0 does on the centred scales.
    """

    if min(hydrophobic_scale(scale).values()) >= 0:
        values = hydrophobic_scale(scale).values()
        return PatchSettings(hydrophobic_threshold=sum(values) / len(values))

    return PatchSettings()


def pdb_text(structure):
    """Returns the PDB text inside a shipped .pdb.zip."""

    with zipfile.ZipFile(structure) as archive:
        (member,) = archive.namelist()
        return archive.read(member).decode()


def is_oxt(line):
    """True for an ATOM record of an atom named OXT.

    HETATM records are not read by prodes, so an OXT there, such as the one of
    the acetate ligand in 4NZU, is neither counted nor removed.
    """

    return line.startswith("ATOM") and line[12:16].strip() == "OXT"


def features(pdb_file, scale, oxt_value=None):
    """Runs the full feature set, with OXT given oxt_value, or the shipped value when None."""

    with tempfile.TemporaryDirectory(prefix="prodes_oxt_") as directory:
        out_file = str(Path(directory) / "bundle.zip")
        settings = patch_settings(scale)

        if oxt_value is None:
            calculate(pdb_file, out_file, hydro_scale=scale, full_features=True, patches=settings)
        else:
            with mock.patch("prodes.data.oxt_hydrophobicity", lambda scale_dict: oxt_value):
                calculate(pdb_file, out_file, hydro_scale=scale, full_features=True, patches=settings)

        return read_features(out_file).iloc[0]


LEADING_COLUMNS = ["structure", "n_oxt", "scale", "variant"]


def measured():
    """Returns the (structure, scale, variant) triples already in OUT_FILE."""

    if not OUT_FILE.exists():
        return set()

    table = pd.read_csv(OUT_FILE, usecols=["structure", "scale", "variant"])

    return set(zip(table.structure, table.scale, table.variant, strict=True))


def append(row):
    """Appends one measured row to OUT_FILE, writing the header if the file is new."""

    OUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame([row])
    frame = frame[LEADING_COLUMNS + [column for column in frame.columns if column not in LEADING_COLUMNS]]

    if OUT_FILE.exists():
        header = pd.read_csv(OUT_FILE, nrows=0).columns.tolist()
        if header != frame.columns.tolist():
            raise ValueError(f"{OUT_FILE} has different columns from this run; delete it and measure everything again")
        frame.to_csv(OUT_FILE, mode="a", header=False, index=False)
    else:
        frame.to_csv(OUT_FILE, index=False)


def check_old_differs_from_new(name, scale, n_oxt):
    """Raises if a structure with an OXT came out the same under the old and the new value.

    That would mean the replacement of prodes.data.oxt_hydrophobicity did not
    reach the calculation, for instance because set_lipo stopped looking the
    function up at call time, and every old row would be a copy of the new one.
    """

    if not n_oxt:
        return

    table = pd.read_csv(OUT_FILE)
    rows = table[(table.structure == name) & (table.scale == scale)].set_index("variant")
    if rows.loc["old", "SurfMhpSum"] == rows.loc["new", "SurfMhpSum"]:
        raise RuntimeError(f"{name} on {scale}: the old and new OXT values gave the same SurfMhpSum, so the old value was not applied")


def main():
    """Measures the old, new and no_oxt features of every shipped structure on every scale not yet in OUT_FILE."""

    done = measured()
    with tempfile.TemporaryDirectory(prefix="prodes_oxt_structures_") as directory:
        for structure in sorted(STRUCTURES.glob("*.pdb.zip")):
            name = structure.name.split(".")[0]
            lines = pdb_text(structure).splitlines()
            n_oxt = sum(is_oxt(line) for line in lines)

            deposited = Path(directory) / f"{name}.pdb"
            deposited.write_text("\n".join(lines) + "\n")
            stripped = Path(directory) / f"{name}_no_oxt.pdb"
            stripped.write_text("\n".join(line for line in lines if not is_oxt(line)) + "\n")

            for scale in SCALES:
                variants = {"old": (deposited, OLD_OXT_VALUE), "new": (deposited, None)}
                if n_oxt:
                    variants["no_oxt"] = (stripped, None)

                for variant, (pdb_file, oxt_value) in variants.items():
                    if (name, scale, variant) in done:
                        print("already measured", name, scale, variant, flush=True)
                        continue

                    row = features(str(pdb_file), scale, oxt_value)
                    row["structure"], row["n_oxt"], row["scale"], row["variant"] = name, n_oxt, scale, variant
                    append(row)
                    print("measured", name, scale, variant, flush=True)

                check_old_differs_from_new(name, scale, n_oxt)

    print("complete:", OUT_FILE)


if __name__ == "__main__":
    main()
