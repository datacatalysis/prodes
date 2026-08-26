"""Compares the Prodes surface electrostatic potential against an APBS Poisson-Boltzmann solution.

Deliberately not a pytest test. It shells out to APBS and takes minutes per
structure, and is meant to be started by hand::

    source activate prodes
    python scripts/apbs_comparison.py

The question it answers is not whether the Prodes potential is physically
correct, which it is not and does not try to be. Prodes evaluates Coulomb's law
with a uniform relative permittivity of 4, with no solvent screening, no ionic
strength and no dielectric boundary. The question is whether that cheap
approximation ranks surface positions the same way a real Poisson-Boltzmann
solution does, and in particular whether it agrees about where the positive
patches are.

Three models are compared, all sampled at the exact same coordinates, namely the
Prodes surface points. Sampling at identical positions is the crux of the
design: it removes any difference in surface definition, leaving the difference
between the charge models and the solvent physics.

    P  Prodes         formal charges on 7 residue types plus termini, Coulomb, eps_r = 4
    C  matched Coulomb the PQR force field charges, same Coulomb kernel, eps_r = 4
    A  APBS           the same PQR charges, Poisson-Boltzmann with salt and a boundary

P against C isolates the charge model. C against A isolates the solvent physics.

Requires the two binaries from environment_apbs.yml. Nothing from that
environment is imported; both are invoked as subprocesses, located on PATH or
via PRODES_APBS_BIN and PRODES_PDB2PQR_BIN.

scipy is deliberately not used. It is not installed in the prodes env, and the
two things it would provide here, trilinear interpolation on a regular grid and a
tie aware rank correlation, are implemented below. They were cross checked
against scipy directly, agreeing to 4e-15 on interpolation and exactly on
Spearman over heavily tied data, and tests/test_apbs_comparison_math.py pins them
against analytic values so the check runs without the dependency.

Raw output goes to apbs_results/, which is git ignored, and is flushed after
every structure so a run killed part way through still leaves usable data. Each
structure's directory holds the run record: the resolved constants, the exact
APBS input used, the tool versions and the git commit.
"""

import json
import os
import platform
import shutil
import socket
import subprocess
import time
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd

from prodes.run import construct_surface_grid, prepare_structure

REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = REPO_ROOT / "apbs_results"

# The structures and their PROPKA output live in a sibling repository, which is
# not public and is not a dependency of Prodes itself. Override with
# PRODES_STRUCTURE_DIR and PRODES_PKA_DIR to point this at your own set.
DATA_ROOT = Path(os.environ.get("PRODES_DATA_ROOT", REPO_ROOT.parent / "biochai" / "data"))
STRUCTURE_DIR = Path(os.environ.get("PRODES_STRUCTURE_DIR", DATA_ROOT / "08_neijenhuis" / "structures"))
PKA_DIR = Path(os.environ.get("PRODES_PKA_DIR", DATA_ROOT / "09_neijenhuis_structural_features"))

# APBS and pdb2pqr come from a separate environment, built from
# environment_apbs.yml, because their dependency stacks have no business
# constraining the solve for the package Prodes ships. Found on PATH by default;
# set PRODES_APBS_BIN and PRODES_PDB2PQR_BIN, or activate that environment, if
# they are somewhere this cannot see.
APBS_BIN = Path(os.environ.get("PRODES_APBS_BIN") or shutil.which("apbs") or "apbs")
PDB2PQR_BIN = Path(os.environ.get("PRODES_PDB2PQR_BIN") or shutil.which("pdb2pqr") or "pdb2pqr")

# The 13 commercial proteins of Neijenhuis (2025). Their retention data is
# published, and their structures and PROPKA output already exist.
#
# convert_propka keys the pKa dictionary on residue number alone and ignores the
# chain column, so in a multi chain structure every chain silently receives chain
# A's pKa. The "protonation matched" premise is therefore false for the four
# multi chain entries regardless of what pdb2pqr does, and they are excluded from
# the headline and reported separately. That is a Prodes bug, tracked separately.
STRUCTURES = {
    "1AO6": ("serum albumin", 2),
    "1BSQ": ("beta-lactoglobulin", 1),
    "1CF3": ("glucose oxidase", 1),
    "1F6R": ("catalase-like", 6),
    "1F8N": ("lipoxygenase", 1),
    "1OVT": ("transferrin", 1),
    "1TRH": ("trypsin-like", 1),
    "3LA4": ("urease-like", 1),
    "4F5S": ("amylase", 2),
    "4PEP": ("pepsin", 1),
    "6FRV": ("hydrogenase", 1),
    "6PO0": ("catalase", 4),
    "AF-P01070": ("trypsin inhibitor A", 1),
}

PH = 7.0
R_PROBE = 1.4
FORCEFIELD = "AMBER"

# APBS conditions. pdie 4 matches the uniform relative permittivity hidden in the
# Prodes denominator. The generated pdb2pqr template carries no ion statement at
# all and uses the linearised equation, so the elec block is written from these
# constants rather than patched, or the headline would silently be measured at
# zero ionic strength.
#
# Scoped deliberately to one condition: pH 7, pdie 4, 150 mM. A pilot showed pdie
# 2 against 4 moves Spearman by 0.005, so a sensitivity sweep buys almost nothing
# for the runtime and disk it costs. Add 2.0 back here to run it.
PDIE_VALUES = (4.0,)
PDIE_MAIN = 4.0
SDIE = 78.54
ION_CONC_M = 0.150
ION_RADIUS = 2.0
TEMP_K = 298.15

# k_B T / e at 298.15 K. Used only to put both series in volts for readability;
# it cannot change a rank correlation.
KT_PER_E_VOLTS = 0.025693

# Prodes uses 1.6e-19 rather than the CODATA 1.602176634e-19, a 0.136 per cent
# low bias in every value it has ever reported. Reproduced here deliberately, so
# that model C differs from model P only in the charges and not in a constant.
ELEMENTARY_CHARGE = 1.6e-19
VACUUM_PERMITTIVITY = 8.854e-12
PRODES_EPS_R = 4.0

# Cap on the transient distance matrix of the Coulomb sum. The naive dense form
# allocates 1.37 GB for catalase, which contradicts the memory budget the rest of
# the package is built around.
CHUNK_BYTES = 256e6

# The one structure whose potential grid is kept, so the sampler can be
# revalidated later. Every other .dx is deleted as soon as it has been sampled.
REFERENCE_STRUCTURE = "1BSQ"


def strip_waters(pdb_file, out_file):
    """Writes a copy of a PDB with every water record removed, and returns the count dropped.

    pdb2pqr's own --drop-water cannot be relied on: on 6PO0 it logs "Dropping
    water from structure" and then writes 3,816 HOH atoms into the PQR anyway.
    Waters left in place sit exactly where the Prodes surface points are, punch
    low dielectric holes in the solvent and add hundreds of point dipoles that
    Prodes does not have at all, since PDBparser reads only ATOM records. Left
    in, they move Pearson from 0.89 to 0.16 and manufacture a -411 kT/e artefact.
    """

    kept, dropped = [], 0
    for line in pdb_file.read_text().splitlines(keepends=True):
        if line[17:20].strip() == "HOH" or " HOH " in line[:30]:
            dropped += 1
            continue
        kept.append(line)

    out_file.write_text("".join(kept))
    return dropped


def run_pdb2pqr(pdb_file, pqr_file, template_file, log_file):
    """Runs pdb2pqr with PROPKA titration states, returning its completed process.

    --keep-chain is off by default and without it the PQR carries no chain ID to
    join the per residue charge check on. --nodebump freezes heavy atom
    positions; without it pdb2pqr rotates side chains by up to 3.1 Angstrom,
    which would make model C a geometry change as well as a charge change.
    --apbs-input takes a path, it is not a bare flag.
    """

    command = [
        str(PDB2PQR_BIN),
        f"--ff={FORCEFIELD}",
        "--keep-chain",
        "--drop-water",
        "--nodebump",
        "--titration-state-method=propka",
        f"--with-ph={PH}",
        f"--apbs-input={template_file}",
        str(pdb_file),
        str(pqr_file),
    ]
    result = subprocess.run(command, capture_output=True, text=True, cwd=pqr_file.parent)
    log_file.write_text(result.stdout + result.stderr)
    return result


def grid_settings_from_template(template_file):
    """Reads the grid geometry pdb2pqr chose out of its generated APBS input.

    psize's grid sizing is the one part of the generated template worth keeping:
    it pads the fine grid by 10 Angstrom on every face and holds the spacing near
    0.5 Angstrom. Everything else in the template is replaced.
    """

    wanted = ("dime", "cglen", "fglen", "cgcent", "fgcent")
    settings = {}
    for line in template_file.read_text().splitlines():
        parts = line.split()
        if parts and parts[0] in wanted:
            settings[parts[0]] = " ".join(parts[1:])

    missing = [key for key in wanted if key not in settings]
    if missing:
        raise ValueError(f"{template_file} is missing {missing}; pdb2pqr template format changed")

    return settings


def write_apbs_input(in_file, pqr_name, grid, pdie, ion_conc, out_stem):
    """Writes an APBS input file from the module constants and pdb2pqr's grid geometry.

    Written rather than patched. The template pdb2pqr generates has no ion
    statement, uses lpbe and sets pdie 2, so running it unedited measures the
    headline comparison at zero ionic strength with the linearised equation, in
    exactly the direction that flatters Prodes.

    calcenergy is off: the total electrostatic energy the template prints is grid
    self energy dominated and meaningless without a matched reference, and a
    reader shown that number will misread it.
    """

    ions = ""
    if ion_conc > 0:
        ions = f"    ion charge  1 conc {ion_conc:.3f} radius {ION_RADIUS:.1f}\n    ion charge -1 conc {ion_conc:.3f} radius {ION_RADIUS:.1f}\n"

    in_file.write_text(f"""read
    mol pqr {pqr_name}
end
elec
    mg-auto
    dime {grid["dime"]}
    cglen {grid["cglen"]}
    fglen {grid["fglen"]}
    cgcent {grid["cgcent"]}
    fgcent {grid["fgcent"]}
    mol 1
    npbe
    bcfl sdh
{ions}    pdie {pdie:.4f}
    sdie {SDIE:.4f}
    srfm smol
    chgm spl2
    sdens 10.00
    srad 1.40
    swin 0.30
    temp {TEMP_K}
    calcenergy no
    calcforce no
    write pot dx {out_stem}
end
quit
""")


def run_apbs(in_file, log_file):
    """Runs APBS in the directory of its input file and returns the completed process.

    The working directory matters twice over: APBS unconditionally drops an
    io.mc capture file into it, and the paths inside the input are relative.
    """

    result = subprocess.run([str(APBS_BIN), in_file.name], capture_output=True, text=True, cwd=in_file.parent)
    log_file.write_text(result.stdout + result.stderr)
    return result


def read_dx(dx_file):
    """Reads an OpenDX potential grid, returning (origin, spacing, values).

    values is a 3D array indexed [x, y, z]. The data are written z fastest, so
    numpy's default C order reshape is already correct and no transpose is
    needed; the natural wrong guess, x fastest, produces a plausible looking
    scrambled grid rather than an error.

    The three spacings are not equal in general, so each axis takes its own
    delta. A reader that assumes a cubic voxel is wrong by up to 10 per cent per
    axis and looks entirely plausible.
    """

    origin, deltas, counts, n_values, header_lines = None, [], None, None, 0
    with dx_file.open() as handle:
        for line in handle:
            header_lines += 1
            parts = line.split()
            if line.startswith("object 1"):
                counts = tuple(int(value) for value in parts[-3:])
            elif line.startswith("origin"):
                origin = np.array([float(value) for value in parts[1:4]])
            elif line.startswith("delta"):
                deltas.append([float(value) for value in parts[1:4]])
            elif "data follows" in line:
                n_values = int(parts[parts.index("items") + 1])
                break

    if counts is None or origin is None or n_values is None or len(deltas) != 3:
        raise ValueError(f"{dx_file} is not a grid this reader understands")

    unit = dx_file.read_text(errors="ignore")[:400]
    if "kT/e" not in unit:
        raise ValueError(f"{dx_file} does not declare kT/e units; check the APBS settings")

    spacing = np.array([deltas[0][0], deltas[1][1], deltas[2][2]])
    if not np.all(spacing > 0):
        raise ValueError(f"{dx_file} has a non diagonal or degenerate delta matrix")

    # The C parser is what makes a 300 MB grid readable in seconds.
    frame = pd.read_csv(dx_file, skiprows=header_lines, nrows=-(-n_values // 3), sep=r"\s+", header=None, dtype=np.float64)
    values = frame.to_numpy().ravel()[:n_values]

    return origin, spacing, values.reshape(counts)


def trilinear_sample(origin, spacing, grid, points):
    """Interpolates a regular 3D grid at arbitrary points, returning NaN outside it.

    Equivalent to scipy's RegularGridInterpolator with method="linear", which is
    not available in the prodes environment. Checked against it in the tests.
    """

    fractional = (points - origin) / spacing
    lower = np.floor(fractional).astype(np.int64)
    weight = fractional - lower

    inside = np.all((lower >= 0) & (lower < np.array(grid.shape) - 1), axis=1)
    result = np.full(len(points), np.nan)
    if not np.any(inside):
        return result

    index, offset = lower[inside], weight[inside]
    total = np.zeros(len(index))
    for corner in range(8):
        dx, dy, dz = (corner >> 2) & 1, (corner >> 1) & 1, corner & 1
        share = (offset[:, 0] if dx else 1 - offset[:, 0]) * (offset[:, 1] if dy else 1 - offset[:, 1]) * (offset[:, 2] if dz else 1 - offset[:, 2])
        total += share * grid[index[:, 0] + dx, index[:, 1] + dy, index[:, 2] + dz]

    result[inside] = total
    return result


def read_pqr(pqr_file):
    """Reads a PQR into a dataframe of coordinates, charges, radii and residue identity."""

    rows = []
    for line in pqr_file.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        parts = line.split()
        # PQR is whitespace delimited from the coordinates onward, but the
        # leading columns keep PDB fixed widths, so slice those and split the rest.
        rows.append(
            {
                "name": line[12:16].strip(),
                "resname": line[17:20].strip(),
                "chain": line[21:22].strip(),
                "resnum": int(line[22:26]),
                "x": float(parts[-5]),
                "y": float(parts[-4]),
                "z": float(parts[-3]),
                "charge": float(parts[-2]),
                "radius": float(parts[-1]),
            }
        )

    return pd.DataFrame(rows)


def coulomb_potential(points, charge_coords, charges, eps_r=PRODES_EPS_R):
    """Volts at each point from point charges in a uniform medium, chunked to bound memory.

    This is the Prodes kernel: q / (4 pi eps0 eps_r d), summed over every charged
    atom with no distance cutoff, which is what Property_point.set_ep does with
    its default cutoff of 10000 Angstrom.
    """

    charge_coulombs = np.asarray(charges) * ELEMENTARY_CHARGE
    chunk = max(1, int(CHUNK_BYTES // max(len(charge_coords) * 3 * 8, 1)))

    potential = np.empty(len(points))
    for start in range(0, len(points), chunk):
        block = points[start : start + chunk]
        distances = np.linalg.norm(block[:, None, :] - charge_coords[None, :, :], axis=2) * 1e-10
        potential[start : start + chunk] = (charge_coulombs[None, :] / (4 * np.pi * VACUUM_PERMITTIVITY * eps_r * distances)).sum(axis=1)

    return potential


def average_ranks(values):
    """Ranks values, giving tied values their average rank.

    Tie handling is load bearing rather than a detail: set_ep rounds to two
    decimals, so a 21,129 point cloud carries only about 400 distinct values and
    a naive argsort ranking would bias the correlation badly.
    """

    values = np.asarray(values, dtype=np.float64)
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=np.float64)
    ranks[order] = np.arange(1, len(values) + 1, dtype=np.float64)

    sorted_values = values[order]
    start = 0
    for end in range(1, len(values) + 1):
        if end == len(values) or sorted_values[end] != sorted_values[start]:
            if end - start > 1:
                ranks[order[start:end]] = ranks[order[start:end]].mean()
            start = end

    return ranks


def pearson(first, second):
    """Pearson product moment correlation of two equal length arrays."""

    first, second = np.asarray(first, dtype=np.float64), np.asarray(second, dtype=np.float64)
    first, second = first - first.mean(), second - second.mean()
    denominator = np.sqrt((first**2).sum() * (second**2).sum())
    return float((first * second).sum() / denominator) if denominator else np.nan


def spearman(first, second):
    """Tie aware Spearman rank correlation."""

    return pearson(average_ranks(first), average_ranks(second))


def top_decile_auc(score, reference):
    """ROC AUC for retrieving the top decile of reference using score.

    Offset immune and range restriction immune, unlike Spearman, which is
    invariant to the additive offset that turns out to be the whole story when
    Prodes reports no positive surface point at all on a net negative protein.
    Computed from the rank sum identity, so ties are handled correctly.
    """

    positive = reference >= np.quantile(reference, 0.9)
    n_positive, n_negative = int(positive.sum()), int((~positive).sum())
    if n_positive == 0 or n_negative == 0:
        return np.nan

    ranks = average_ranks(score)
    return float((ranks[positive].sum() - n_positive * (n_positive + 1) / 2) / (n_positive * n_negative))


def positive_recall(prodes_ep, reference_ep):
    """Fraction of the reference's positive surface that Prodes also calls positive.

    This is the number the experiment exists to produce, and the one a rank
    correlation cannot see. Spearman is invariant to an additive offset, and an
    additive offset is exactly what a Coulomb sum with no ionic screening adds to
    a net negative protein: every surface point is pushed below zero together,
    the ordering survives, and the positive patches vanish.
    """

    positive = reference_ep > 0
    if not np.any(positive):
        return np.nan

    return float((prodes_ep[positive] > 0).mean())


def positive_precision(prodes_ep, reference_ep):
    """Fraction of the points Prodes calls positive that the reference agrees are positive."""

    claimed = prodes_ep > 0
    if not np.any(claimed):
        return np.nan

    return float((reference_ep[claimed] > 0).mean())


def sign_agreement(first, second):
    """Fraction of points on which two potentials agree in sign, ignoring exact zeros.

    Prodes rounds to two decimals, so any point below 0.005 V becomes exactly
    0.00 and belongs to neither sign. Those are excluded rather than assigned.
    """

    usable = (first != 0) & (second != 0)
    if not np.any(usable):
        return np.nan

    return float((np.sign(first[usable]) == np.sign(second[usable])).mean())


def per_residue_charges(pqr, structure):
    """Compares PQR and Prodes formal charge residue by residue, keyed on chain and number.

    A scalar total is invariant under any compensating pair of disagreements, and
    the likeliest mismatch here is exactly such a pair: an N terminal +1 against a
    C terminal -1 cancels exactly. pdb2pqr also cannot accept a precomputed pKa
    file, so it re-runs propka on its own structure and divergence from the
    Prodes pKa file is guaranteed by construction.
    """

    pqr_totals = pqr.groupby(["chain", "resnum"])["charge"].sum()

    prodes_totals = {}
    for residue in structure.residues:
        key = (residue.chain.name if residue.chain else "", int(residue.number))
        prodes_totals[key] = prodes_totals.get(key, 0.0) + sum(atom.charge(ph=PH) for atom in residue.atoms)

    rows = []
    for key in sorted(set(pqr_totals.index) | set(prodes_totals)):
        pqr_charge = float(pqr_totals.get(key, 0.0))
        prodes_charge = float(prodes_totals.get(key, 0.0))
        rows.append(
            {
                "chain": key[0],
                "resnum": key[1],
                "pqr_charge": round(pqr_charge, 4),
                "prodes_charge": round(prodes_charge, 4),
                "difference": round(pqr_charge - prodes_charge, 4),
            }
        )

    return pd.DataFrame(rows)


def tool_versions():
    """Returns the versions of the external tools, for the run record."""

    versions = {}
    for name, command in (("apbs", [str(APBS_BIN), "--version"]), ("pdb2pqr", [str(PDB2PQR_BIN), "--version"])):
        try:
            result = subprocess.run(command, capture_output=True, text=True)
            versions[name] = (result.stdout + result.stderr).strip().splitlines()[-1][:120]
        except OSError as error:
            versions[name] = f"unavailable: {error}"

    return versions


def write_run_record(out_dir, extra):
    """Writes the resolved constants, tool versions and provenance into the output directory.

    A shell one liner in a markdown file is not a record of a run.
    """

    commit = subprocess.run(["git", "-C", str(REPO_ROOT), "rev-parse", "HEAD"], capture_output=True, text=True)
    record = {
        "utc": datetime.now(UTC).isoformat(),
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "git_commit": commit.stdout.strip(),
        "ph": PH,
        "r_probe": R_PROBE,
        "forcefield": FORCEFIELD,
        "pdie_values": list(PDIE_VALUES),
        "sdie": SDIE,
        "ion_conc_M": ION_CONC_M,
        "temp_K": TEMP_K,
        "kt_per_e_volts": KT_PER_E_VOLTS,
        "tools": tool_versions(),
    }
    record.update(extra)
    (out_dir / "run_record.json").write_text(json.dumps(record, indent=2, default=str))


def prodes_surface(pdb_file, pka_file):
    """Returns (structure, coords, ep) for the Prodes surface grid of one structure.

    The points are the centres of the occupied cells of a 1 Angstrom grid laid
    over the Shrake-Rupley cloud, not the Shrake-Rupley points themselves, so
    they lie on an exact cubic lattice.
    """

    structure = prepare_structure(str(pdb_file), str(pka_file) if pka_file and pka_file.exists() else None)
    points = construct_surface_grid(structure, R_PROBE)
    coords = np.array([[point.x, point.y, point.z] for point in points])

    charged = [atom for atom in structure.atoms if atom.charge(ph=PH) != 0]
    if not charged:
        raise ValueError(f"{pdb_file} has no charged atoms at pH {PH}")

    charge_coords = np.array([[atom.x, atom.y, atom.z] for atom in charged])
    charges = np.array([atom.charge(ph=PH) for atom in charged])

    return structure, coords, np.round(coulomb_potential(coords, charge_coords, charges), 2)


def compare_structure(code, out_dir):
    """Runs the whole comparison for one structure and returns its summary rows."""

    out_dir.mkdir(parents=True, exist_ok=True)
    pdb_file = STRUCTURE_DIR / f"{code}.pdb"
    pka_file = PKA_DIR / f"{code}_pka.json"
    if not pdb_file.exists():
        raise FileNotFoundError(f"{pdb_file} not found; this script needs the biochai data checked out")

    stripped = out_dir / f"{code}_nw.pdb"
    dropped = strip_waters(pdb_file, stripped)

    pqr_file, template_file = out_dir / f"{code}.pqr", out_dir / f"{code}_template.in"
    result = run_pdb2pqr(stripped, pqr_file, template_file, out_dir / f"{code}_pdb2pqr.log")
    if result.returncode != 0 or not pqr_file.exists():
        raise RuntimeError(f"pdb2pqr failed on {code}, see {out_dir / f'{code}_pdb2pqr.log'}")

    pqr_text = pqr_file.read_text()
    if " HOH " in pqr_text:
        raise RuntimeError(f"{code}: waters survived into the PQR, which invalidates the comparison")

    structure, coords, ep_prodes = prodes_surface(pdb_file, pka_file)
    pqr = read_pqr(pqr_file)

    charge_diff = per_residue_charges(pqr, structure)
    charge_diff.to_csv(out_dir / f"{code}_residue_charges.csv", index=False)
    mismatched = int((charge_diff["difference"].abs() > 0.05).sum())

    charged_pqr = pqr[pqr["charge"] != 0]
    ep_coulomb = coulomb_potential(coords, charged_pqr[["x", "y", "z"]].to_numpy(), charged_pqr["charge"].to_numpy())

    grid = grid_settings_from_template(template_file)
    rows = []
    for pdie in PDIE_VALUES:
        stem = f"{code}_pdie{int(pdie)}"
        in_file = out_dir / f"{stem}.in"
        write_apbs_input(in_file, pqr_file.name, grid, pdie, ION_CONC_M, stem)

        started = time.time()
        apbs_result = run_apbs(in_file, out_dir / f"{stem}_apbs.log")
        dx_file = out_dir / f"{stem}.dx"
        log_text = (out_dir / f"{stem}_apbs.log").read_text()
        if apbs_result.returncode != 0 or not dx_file.exists() or dx_file.stat().st_mtime < in_file.stat().st_mtime:
            raise RuntimeError(f"APBS produced no usable grid for {stem}; exit {apbs_result.returncode}")
        if "Vio_ctor2" in log_text:
            raise RuntimeError(f"APBS reported Vio_ctor2 for {stem}, which it does while still exiting 0")

        origin, spacing, values = read_dx(dx_file)
        ep_apbs_kte = trilinear_sample(origin, spacing, values, coords)
        outside = int(np.isnan(ep_apbs_kte).sum())
        if outside:
            raise RuntimeError(f"{stem}: {outside} surface points fell outside the fine grid; grow dime, never fglen alone")

        ep_apbs_volts = ep_apbs_kte * KT_PER_E_VOLTS
        pd.DataFrame(
            {
                "x": coords[:, 0],
                "y": coords[:, 1],
                "z": coords[:, 2],
                "ep_prodes_volts": ep_prodes,
                "ep_coulomb_pqr_volts": ep_coulomb,
                "ep_apbs_kte": ep_apbs_kte,
                "ep_apbs_volts": ep_apbs_volts,
            }
        ).to_csv(out_dir / f"{stem}_points.csv.gz", index=False, compression="gzip")

        slope = float(np.polyfit(ep_prodes, ep_apbs_volts, 1)[0])
        rows.append(
            {
                "id": code,
                "protein": STRUCTURES[code][0],
                "chains": STRUCTURES[code][1],
                "pdie": pdie,
                "ion_conc_M": ION_CONC_M,
                "n_points": len(coords),
                "waters_dropped": dropped,
                "residues_mismatched": mismatched,
                "prodes_min_V": round(float(ep_prodes.min()), 3),
                "prodes_max_V": round(float(ep_prodes.max()), 3),
                "apbs_min_kte": round(float(ep_apbs_kte.min()), 2),
                "apbs_max_kte": round(float(ep_apbs_kte.max()), 2),
                "prodes_fraction_positive": round(float((ep_prodes > 0).mean()), 4),
                "apbs_fraction_positive": round(float((ep_apbs_volts > 0).mean()), 4),
                "sign_agreement": round(sign_agreement(ep_prodes, ep_apbs_volts), 4),
                "apbs_positive_recall": round(positive_recall(ep_prodes, ep_apbs_volts), 4),
                "prodes_positive_precision": round(positive_precision(ep_prodes, ep_apbs_volts), 4),
                "top_decile_auc": round(top_decile_auc(ep_prodes, ep_apbs_volts), 4),
                "spearman_P_A": round(spearman(ep_prodes, ep_apbs_volts), 4),
                "spearman_P_C": round(spearman(ep_prodes, ep_coulomb), 4),
                "spearman_C_A": round(spearman(ep_coulomb, ep_apbs_volts), 4),
                "pearson_P_A": round(pearson(ep_prodes, ep_apbs_volts), 4),
                "slope_volts_per_volt": round(slope, 5),
                "apbs_seconds": round(time.time() - started, 1),
            }
        )
        print(f"  pdie {pdie:.0f}: spearman {rows[-1]['spearman_P_A']:.3f}  sign {rows[-1]['sign_agreement']:.3f}  auc {rows[-1]['top_decile_auc']:.3f}")

        # The grid is 43 to 300 MB and its useful content is now sampled. Keep
        # only the reference structure's, so the sampler can be revalidated.
        if code != REFERENCE_STRUCTURE or pdie != PDIE_MAIN:
            dx_file.unlink()

    write_run_record(out_dir, {"structure": code, "waters_dropped": dropped, "residues_mismatched": mismatched})
    return rows


def main():
    """Runs the comparison over every structure and writes the summary table."""

    if not APBS_BIN.exists() or not PDB2PQR_BIN.exists():
        raise FileNotFoundError(
            f"apbs and pdb2pqr not found ({APBS_BIN}, {PDB2PQR_BIN}). Create the environment from "
            "environment_apbs.yml and activate it, or set PRODES_APBS_BIN and PRODES_PDB2PQR_BIN."
        )

    RESULTS_DIR.mkdir(exist_ok=True)
    summary_file = RESULTS_DIR / "apbs_vs_prodes_summary.csv"
    all_rows, failures = [], {}

    for code in STRUCTURES:
        print(f"{code} ({STRUCTURES[code][0]}, {STRUCTURES[code][1]} chain(s))")
        started = time.time()
        try:
            all_rows.extend(compare_structure(code, RESULTS_DIR / code))
        except Exception as error:  # noqa: BLE001 - one bad structure must not lose the rest
            failures[code] = str(error)
            print(f"  FAILED: {error}")
            continue

        print(f"  done in {time.time() - started:.0f}s")
        # Flushed per structure so a killed run still leaves usable data.
        pd.DataFrame(all_rows).to_csv(summary_file, index=False)

    pd.DataFrame(all_rows).to_csv(summary_file, index=False)
    if failures:
        (RESULTS_DIR / "failures.json").write_text(json.dumps(failures, indent=2))
        print(f"\n{len(failures)} structure(s) failed: {sorted(failures)}")

    print(f"\nwrote {summary_file}")
    if shutil.which("column"):
        print(pd.DataFrame(all_rows).to_string(index=False))


if __name__ == "__main__":
    main()
