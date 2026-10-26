"""Tests for the hydrophobicity projected onto a surface point, and the value given to OXT.

The expected values are read from hydrophobicity.json directly rather than through
prodes.data, so a test cannot pass by asking the code under test what it does.
"""

import json
import math
from pathlib import Path

import pytest

import prodes
from prodes.core.atom import Atom
from prodes.core.point import Property_point
from tests.pdb_records import atom_line, write_structure

SCALES = json.loads((Path(__file__).resolve().parent.parent / "src" / "prodes" / "data" / "hydrophobicity.json").read_text())

DISTANCE = 2.0


def oxt_atom(residue_name, chain="A", x=DISTANCE):
    """Returns an OXT atom of the given residue, DISTANCE angstrom from the origin by default."""

    return Atom("ATOM", "OXT", residue_name, chain, 1, x, 0, 0, element="O")


def lipo_at_origin(atoms, scale):
    """Returns the hydrophobicity projected onto a point at the origin."""

    point = Property_point(0, 0, 0)
    point.set_lipo(atoms, scale=scale)

    return point.lipo


@pytest.mark.parametrize("scale", sorted(SCALES))
def test_oxt_takes_the_asp_value_of_the_selected_scale(scale):
    """OXT on a C-terminal PHE contributes the ASP value of the scale, not PHE's and not a constant.

    PHE is the most hydrophobic residue on mj and mj_scaled, so a regression to
    "the residue's own value" fails here. The four ASP values differ from each
    other and from 1.0, so a regression to one constant fails on at least three
    of the scales.
    """

    expected = SCALES[scale]["ASP"] * math.exp(-DISTANCE)

    assert lipo_at_origin([oxt_atom("PHE")], scale) == pytest.approx(expected)
    assert SCALES[scale]["ASP"] != SCALES[scale]["PHE"]


@pytest.mark.parametrize("scale", sorted(SCALES))
def test_oxt_value_is_on_the_hydrophilic_side_of_the_scale(scale):
    """The ASP value given to OXT lies below the scale's midpoint, on every shipped scale."""

    values = SCALES[scale].values()
    midpoint = (min(values) + max(values)) / 2

    assert lipo_at_origin([oxt_atom("PHE")], scale) < midpoint * math.exp(-DISTANCE)


def test_oxt_value_does_not_depend_on_its_residue():
    """The carboxylate of a C-terminal LYS and of a C-terminal PHE contribute the same."""

    assert lipo_at_origin([oxt_atom("LYS")], "mj_scaled") == lipo_at_origin([oxt_atom("PHE")], "mj_scaled")


def test_other_atoms_of_the_terminal_residue_keep_the_residue_value():
    """Only OXT is special: the backbone O of the same PHE still takes PHE's value."""

    o_atom = Atom("ATOM", "O", "PHE", "A", 1, DISTANCE, 0, 0, element="O")

    assert lipo_at_origin([o_atom], "mj_scaled") == pytest.approx(SCALES["mj_scaled"]["PHE"] * math.exp(-DISTANCE))


def residue_with_oxt(serial, residue_name, chain, x):
    """Returns the backbone and OXT of one C-terminal residue, its OXT at (x, 0, 0)."""

    atoms = [
        ("N", x - 3.0, 1.0, 0.0, "N"),
        ("CA", x - 2.0, 1.0, 0.0, "C"),
        ("C", x - 1.0, 1.0, 0.0, "C"),
        ("O", x - 1.0, 2.2, 0.0, "O"),
        ("OXT", x, 0.0, 0.0, "O"),
    ]

    return [
        atom_line(serial + offset, name, residue_name, chain, 1, at_x, at_y, at_z, element) for offset, (name, at_x, at_y, at_z, element) in enumerate(atoms)
    ]


def test_every_chain_terminus_gets_the_oxt_value(tmp_path):
    """The parser keeps the OXT of every chain, and each one, projected on its own, contributes the ASP value.

    Chain A ends in PHE and chain B in LEU, so neither OXT passes by taking its
    residue's value.
    """

    path = write_structure(tmp_path / "two_chains.pdb", residue_with_oxt(1, "PHE", "A", 0.0) + residue_with_oxt(6, "LEU", "B", 100.0))
    structure = prodes.read(path)

    oxt_atoms = [atom for atom in structure.heavy_atoms if atom.name == "OXT"]
    assert sorted(atom.chain_name for atom in oxt_atoms) == ["A", "B"]

    for atom in oxt_atoms:
        point = Property_point(atom.x, atom.y - DISTANCE, atom.z)
        point.set_lipo([atom], scale="kd")

        assert point.lipo == pytest.approx(SCALES["kd"]["ASP"] * math.exp(-DISTANCE))
