"""Residues that are not one of the 20 standard amino acids, issue #37.

A structure predicted from a sequence holding an X has a residue named UNK,
which the residue table has no radius, pKa or hydrophobicity for. By default
it is left out with a warning and recorded, and a chain end it occupied
carries no terminal charge; with drop_unknown_residues=False the structure is
refused when it is read, naming every such residue.
"""

import json
import logging
import zipfile

import pytest

import prodes
from prodes.calculations.patches import PatchSettings
from prodes.io import parser
from prodes.io.pka_converter import ANY_CHAIN
from prodes.output import run_metadata
from prodes.run import calculate, without_removed_residues
from tests.pdb_records import atom_line, write_structure

# Backbone and CB, the atoms Boltz-2 builds for an UNK. GLY has no CB, but an
# extra atom changes nothing the parser decides, so every residue gets one.
BACKBONE = [("N", 0.0, 0.0, 0.0, "N"), ("CA", 1.2, 0.5, 0.6, "C"), ("C", 2.4, 0.0, 0.0, "C"), ("O", 2.4, -1.2, -0.6, "O"), ("CB", 1.2, 2.0, 1.4, "C")]


def residue_lines(serial, residue_name, chain, number, x, y, altloc="", occupancy=1.00):
    """Returns the five heavy atoms of one residue, placed at x along its chain."""

    return [
        atom_line(serial + offset, name, residue_name, chain, number, x + dx, y + dy, dz, element, altloc=altloc, occupancy=occupancy)
        for offset, (name, dx, dy, dz, element) in enumerate(BACKBONE)
    ]


def chain_lines(chain, names, y=0.0, serial=1):
    """Returns the records of one chain, residue n named names[n - 1], 3.8 A apart."""

    lines = []
    for number, residue_name in enumerate(names, start=1):
        lines += residue_lines(serial + len(lines), residue_name, chain, number, 3.8 * number, y)
    return lines


def pdb_text(*chains):
    """Returns PDB text for the given chain records, chains 20 A apart."""

    lines = []
    for offset, (chain, names) in enumerate(chains):
        lines += chain_lines(chain, names, y=20.0 * offset, serial=len(lines) + 1)
    return "\n".join(lines) + "\nEND\n"


# Unknown at the start, in the middle and at the end of chain A; chain B has none.
MIXED = pdb_text(("A", ["UNK", "VAL", "ALA", "UNK", "GLY", "UNK"]), ("B", ["SER", "LYS", "ASP"]))


def parse(text, drop):
    return parser.parse_pdb_text(text, "test", drop_unknown_residues=drop)


def residue(structure, chain, number):
    return next(r for r in structure.residues if r.chain.name == chain and r.number == number)


def terminal_groups(residue):
    return [key for group in (residue.pkas or []) for key in group if key in ("N+", "C-")]


def test_by_default_an_unknown_residue_is_left_out_with_a_warning(caplog):
    """An unattended run finishes: the structure is read, and the log names what was left out."""

    with caplog.at_level(logging.WARNING, logger="prodes.io.parser"):
        structure = parser.parse_pdb_text(MIXED, "test")

    assert structure.removed_residues == ["UNK A1", "UNK A4", "UNK A6"]
    assert "left out 3 residues that are not one of the 20 standard amino acids: UNK A1, UNK A4, UNK A6" in caplog.text


def test_without_the_option_an_unknown_residue_is_refused_when_the_structure_is_read():
    """Every unknown residue is named by chain and number, and the setting that removes them is named."""

    with pytest.raises(ValueError) as raised:
        parse(MIXED, drop=False)

    message = str(raised.value)
    assert "UNK A1, UNK A4, UNK A6" in message
    assert "drop_unknown_residues" in message
    assert "--no-drop-unknown-residues" in message


def test_one_unknown_residue_is_counted_in_the_singular():
    """One X in a sequence is the common case, so its message has to read right."""

    with pytest.raises(ValueError, match="holds 1 residue that is not one of the 20 standard amino acids, so Prodes has no data for it: UNK A1"):
        parse(pdb_text(("A", ["UNK", "VAL", "SER"])), drop=False)


def test_prodes_read_has_the_same_default_and_the_same_option(tmp_path):
    """The library entry point reads through the same parser, so it needs the same option."""

    path = tmp_path / "mixed.pdb"
    path.write_text(MIXED)

    assert [r.label for r in prodes.read(str(path)).residues] == ["A2", "A3", "A5", "B1", "B2", "B3"]
    with pytest.raises(ValueError, match="UNK A1"):
        prodes.read(str(path), drop_unknown_residues=False)


def test_the_option_removes_them_and_lists_them_in_file_order(caplog):
    """The removed residues are kept on the structure and named in the log, in the order the file holds them."""

    structure = parse(pdb_text(("B", ["UNK", "VAL"]), ("A", ["SER", "UNK", "LYS"])), drop=True)
    assert structure.removed_residues == ["UNK B1", "UNK A2"]

    with caplog.at_level(logging.WARNING, logger="prodes.io.parser"):
        structure = parse(MIXED, drop=True)

    assert [r.label for r in structure.residues] == ["A2", "A3", "A5", "B1", "B2", "B3"]
    assert all(atom.residue_name != "UNK" for atom in structure.atoms)
    assert structure.removed_residues == ["UNK A1", "UNK A4", "UNK A6"]
    assert "UNK A1, UNK A4, UNK A6" in caplog.text


def test_a_removed_first_residue_takes_the_n_terminus_with_it():
    """The new first residue is bonded to the removed one, so it has no free amine."""

    structure = parse(MIXED, drop=True)
    first = residue(structure, "A", 2)

    assert first.chain.n_terminus_removed is True
    assert first.terminus is None
    assert "N+" not in terminal_groups(first)
    assert next(atom for atom in first.atoms if atom.name == "N").charge(ph=7) == 0


def test_a_removed_last_residue_takes_the_c_terminus_with_it():
    """The same for the carboxylate at the other end."""

    structure = parse(MIXED, drop=True)
    last = residue(structure, "A", 5)

    assert last.chain.c_terminus_removed is True
    assert last.terminus is None
    assert "C-" not in terminal_groups(last)
    assert next(atom for atom in last.atoms if atom.name == "C").charge(ph=7) == 0


def test_a_chain_with_nothing_removed_keeps_both_termini():
    """Chain B is untouched, termini, groups and charges alike."""

    structure = parse(MIXED, drop=True)
    first, last = residue(structure, "B", 1), residue(structure, "B", 3)

    assert (first.terminus, last.terminus) == ("N", "C")
    assert terminal_groups(first) == ["N+"]
    assert terminal_groups(last) == ["C-"]
    assert next(atom for atom in first.atoms if atom.name == "N").charge(ph=7) == 1
    assert next(atom for atom in last.atoms if atom.name == "C").charge(ph=7) == -1


def test_a_removed_middle_residue_changes_no_terminus():
    """Removing A4 leaves A3 and A5 as ordinary residues."""

    structure = parse(pdb_text(("A", ["SER", "VAL", "ALA", "UNK", "GLY", "LYS"])), drop=True)

    assert structure.removed_residues == ["UNK A4"]
    assert [r.terminus for r in structure.residues] == ["N", None, None, None, "C"]
    assert terminal_groups(residue(structure, "A", 1)) == ["N+"]
    assert terminal_groups(residue(structure, "A", 6)) == ["C-"]


@pytest.mark.parametrize(
    ("names", "kept", "terminus", "group", "atom_name", "charge"),
    [
        (["UNK", "VAL"], 2, "C", "C-", "C", -1),
        (["VAL", "UNK"], 1, "N", "N+", "N", 1),
    ],
)
def test_a_chain_left_with_one_residue_keeps_the_end_that_was_not_removed(names, kept, terminus, group, atom_name, charge):
    """The terminus attribute, the pKa group and the atom charge agree on which end survived."""

    structure = parse(pdb_text(("A", names)), drop=True)
    remaining = residue(structure, "A", kept)

    assert remaining.terminus == terminus
    assert terminal_groups(remaining) == [group]
    assert next(atom for atom in remaining.atoms if atom.name == atom_name).charge(ph=7) == charge


def test_a_structure_of_nothing_but_unknown_residues_is_refused_either_way():
    """There is nothing left to describe, with or without the option."""

    text = pdb_text(("A", ["UNK", "UNK"]))

    with pytest.raises(ValueError, match="UNK A1, UNK A2"):
        parse(text, drop=False)
    with pytest.raises(ValueError, match="no residue that is one of the 20 standard amino acids"):
        parse(text, drop=True)


def test_a_structure_without_unknown_residues_is_the_same_with_and_without_the_option():
    """The option changes nothing for a structure it has nothing to do for."""

    text = pdb_text(("A", ["SER", "VAL", "LYS"]), ("B", ["ASP", "GLY"]))
    plain, dropped = parse(text, drop=False), parse(text, drop=True)

    def summary(structure):
        return [(r.label, r.name, r.terminus, r.pkas, len(r.atoms)) for r in structure.residues]

    assert summary(plain) == summary(dropped)
    for ph in (4, 7, 10):
        assert [atom.charge(ph=ph, formal=False) for atom in plain.atoms] == [atom.charge(ph=ph, formal=False) for atom in dropped.atoms]
    assert dropped.removed_residues == []
    assert not any(chain.n_terminus_removed or chain.c_terminus_removed for chain in dropped.chains)


@pytest.mark.parametrize(("winner", "occupancies"), [("ALA", (0.6, 0.4)), ("UNK", (0.4, 0.6))])
def test_the_elected_conformation_decides_whether_a_residue_is_unknown(winner, occupancies):
    """A position modelled as ALA and as UNK is judged by the amino acid the election keeps."""

    lines = chain_lines("A", ["SER"])
    lines += residue_lines(len(lines) + 1, "ALA", "A", 2, 7.6, 0.0, altloc="A", occupancy=occupancies[0])
    lines += residue_lines(len(lines) + 1, "UNK", "A", 2, 7.6, 0.0, altloc="B", occupancy=occupancies[1])
    lines += residue_lines(len(lines) + 1, "LYS", "A", 3, 11.4, 0.0)
    text = "\n".join(lines) + "\nEND\n"

    if winner == "ALA":
        assert residue(parse(text, drop=False), "A", 2).name == "ALA"
    else:
        with pytest.raises(ValueError, match="UNK A2"):
            parse(text, drop=False)
        dropped = parse(text, drop=True)
        assert dropped.removed_residues == ["UNK A2"]
        # Removed, so not also reported as renamed in the run record.
        assert dropped.alternate_conformers.summary()["residues_renamed"] == []


def test_a_hetatm_parse_is_not_checked():
    """Ligands and waters are never in the residue table, and are read separately."""

    text = MIXED.replace("END\n", "") + atom_line(999, "O", "HOH", "W", 1, 0.0, 50.0, 0.0, "O").replace("ATOM  ", "HETATM", 1) + "\nEND\n"

    structure = parser.parse_pdb_text(text, "test", identifier="HETATM")

    assert [r.name for r in structure.residues] == ["HOH"]


def test_a_predicted_terminal_pka_for_the_new_end_residue_is_ignored_and_reported(caplog):
    """A pKa file can give the residue after a removed one an N+ value; it has no such group."""

    structure = parse(MIXED, drop=True)

    with caplog.at_level(logging.WARNING, logger="prodes.core.structure"):
        structure.redo_pkas({"A": {2: [{"N+": 7.5}]}, ANY_CHAIN: {}})

    assert "N+" not in terminal_groups(residue(structure, "A", 2))
    assert "VAL A2 a N+ pKa, which that residue does not have" in caplog.text


def test_the_run_record_lists_the_removed_residues():
    """Empty when none were removed, so the key is always there to read."""

    import numpy as np

    assert run_metadata("x.pdb", {}, 0, np.array([]))["removed_residues"] == []
    assert run_metadata("x.pdb", {}, 0, np.array([]), removed_residues=["UNK A1"])["removed_residues"] == ["UNK A1"]


def test_a_residue_set_loses_the_removed_residues(caplog):
    """A residue with no atoms has no surface to restrict a patch to; the rest of the set is kept."""

    patches = PatchSettings(residue_sets={"Fv": ["A1", "A2", "A3"]})

    with caplog.at_level(logging.WARNING, logger="prodes.run"):
        kept = without_removed_residues(patches, ["UNK A1"])

    assert kept.residue_sets == {"Fv": ["A2", "A3"]}
    assert "residue set Fv loses the residues that were left out as unknown: A1" in caplog.text
    assert without_removed_residues(patches, []) is patches


def test_a_residue_set_of_nothing_but_removed_residues_is_refused():
    """An empty set would describe nothing, which is not what whoever named it asked for."""

    with pytest.raises(ValueError, match="residue set CDR names only residues that were left out as unknown"):
        without_removed_residues(PatchSettings(residue_sets={"CDR": ["A1"]}), ["UNK A1"])


def test_calculate_runs_by_default_and_records_what_it_left_out(tmp_path):
    """End to end: the bundle is written and its run record says what was left out."""

    pdb = write_structure(tmp_path / "mixed.pdb", MIXED.removesuffix("END\n").splitlines())
    out = tmp_path / "mixed.zip"

    with pytest.raises(ValueError, match="UNK A1"):
        calculate(pdb, str(out), drop_unknown_residues=False)

    calculate(pdb, str(out), patches=PatchSettings(residue_sets={"A": ["A1", "A2", "A3"]}))

    with zipfile.ZipFile(out) as bundle:
        record = json.loads(bundle.read(next(name for name in bundle.namelist() if name.endswith("prodes_run.json"))))
    assert record["removed_residues"] == ["UNK A1", "UNK A4", "UNK A6"]
    assert record["settings"]["drop_unknown_residues"] is True
    assert record["settings"]["patches"]["residue_sets"] == {"A": ["A2", "A3"]}
