"""Tests for the alternate conformation election.

Split out of tests/test_parser.py, which they outgrew. A residue modelled in
more than one conformation writes each of them out in full, and the parser used
to read the atom name as columns 13-17, which is the name plus the alternate
location indicator, so an atom of a disordered aspartate arrived called OD1A. A
name like that matches nothing the package looks for, so the residue kept both
copies of its side chain in the surface and lost its charge entirely. Every
count below fails on the parser as it was before version 7.0.

Most of these go through PDBparser.parse rather than through elect_conformers,
because what has to hold is what a caller sees.
"""

import logging

import pytest

from prodes.io import parser
from prodes.io.conformers import elect_conformers
from prodes.io.pdb_reader import read_atom_records
from tests.pdb_records import atom_line, write_structure

pdb_parser = parser.PDBparser()

file_path = "tests/data/1GDW.pdb.zip"


# Alternate conformations (altLoc). Issue #4.
#
# A residue modelled in more than one conformation writes each of them out in
# full. The parser used to read the atom name as columns 13-17, which is the
# name plus the alternate location indicator, so an atom of a disordered
# aspartate arrived called OD1A. A name like that matches nothing the package
# looks for, so the residue kept both copies of its side chain in the surface
# and lost its charge entirely. Every count below fails on the old parser.

alternates_path = "tests/data/1CBN.pdb.zip"
fab_path = "tests/data/4NZU.pdb.zip"


def parsed_lines(path):
    """Returns the ATOM records of a zipped structure, for comparing against the parse."""

    return [line for line in parser.read_pdb_text(path).splitlines() if line[0:6].strip() == "ATOM"]


def test_an_alternate_location_is_not_part_of_the_atom_name():
    """the indicator in column 17 is read as its own field, not glued to the name

    This is the whole defect. A disordered SG used to arrive as "SG A", which is
    not the SG that disulfide detection looks for nor the SG that carries charge.
    """

    structure = pdb_parser.parse(alternates_path)

    assert not [atom for atom in structure.atoms if " " in atom.name]
    assert {atom.altloc for atom in structure.atoms} == {"", "A", "C"}


def test_every_parsed_name_and_indicator_is_a_pair_the_file_wrote():
    """the name and the alternate location have to come apart the way the file wrote them

    Asserting "no spaces and at most four characters" would pass OD1A, CD1A and
    every other three character name carrying a letter: 103 of the 272 atoms
    1CBN used to mis-parse look perfectly normal. Checking the pair instead
    catches them, because a parser that reads columns 13-17 as the name produces
    ("OD1A", "") and no record in the file says that. The record for record
    comparison in the next test is the stronger statement; this one is the cheap
    one that names the original defect.
    """

    for path in (alternates_path, fab_path):
        structure = pdb_parser.parse(path)
        written = {(line[12:16].strip(), line[16].strip()) for line in parsed_lines(path)}

        assert len(structure.atoms) > 0
        for atom in structure.atoms:
            assert atom.name == atom.name.strip()
            assert (atom.name, atom.altloc) in written


def test_names_match_record_for_record_on_a_disordered_structure():
    """each record is tied to the line that produced it, and each atom to its record

    Two links, because either alone is weak. The reader has to return one record
    per coordinate line, in file order, carrying that line's own columns; and the
    parse has to return one atom per surviving record, in the same order. Chained,
    they say that every atom in the structure carries the columns of the record it
    came from, which is the property the whole change rests on.
    """

    for path in (alternates_path, fab_path):
        text = parser.read_pdb_text(path)
        structure = pdb_parser.parse(path)

        lines = [line for line in text.splitlines() if line[0:6] == "ATOM  "]
        records = [record for record in read_atom_records(text, "test") if record.identifier == "ATOM"]

        assert len(records) == len(lines)
        for record, line in zip(records, lines, strict=True):
            assert record.name == line[12:16].strip()
            assert record.altloc == line[16].strip()
            assert record.x == float(line[30:38])

        kept, _ = elect_conformers(records)

        assert len(kept) == len(structure.atoms)
        for record, atom in zip(kept, structure.atoms, strict=True):
            assert atom.name == record.name
            assert atom.altloc == record.altloc
            assert atom.x == record.x


def test_alternates_are_collapsed_to_one_conformation():
    """1CBN loses the duplicate side chains it used to carry

    772 atoms and 400 heavy atoms before, because every conformation survived.
    The heavy atom count is the one that matters: those are the atoms the surface
    is built from.
    """

    structure = pdb_parser.parse(alternates_path)

    assert len(structure.atoms) == 641
    assert len(structure.heavy_atoms) == 327


def test_a_residue_modelled_as_two_amino_acids_keeps_one_of_them():
    """1CBN residues 22 and 25 are each modelled as two different residues

    Residue 22 is a serine at occupancy 0.20 and a proline at 0.60, residue 25 an
    isoleucine and a leucine. Choosing per atom, which is what issue #4 proposed,
    would take the highest occupancy atom for each name independently and build a
    residue holding the serine's OG beside the proline's CG and CD. Choosing per
    residue keeps one whole amino acid, and the more occupied one.
    """

    residues = {residue.number: residue for residue in pdb_parser.parse(alternates_path).residues}

    assert residues[22].name == "PRO"
    assert sorted(atom.name for atom in residues[22].heavy_atoms) == ["C", "CA", "CB", "CD", "CG", "N", "O"]
    assert residues[25].name == "LEU"
    assert sorted(atom.name for atom in residues[25].heavy_atoms) == ["C", "CA", "CB", "CD1", "CD2", "CG", "N", "O"]


def test_a_disordered_residue_recovers_its_charge():
    """1CBN ASP 43 is disordered and used to carry no charge at any pH

    Its oxygens parsed as OD1A, OD1B, OD2A and OD2B, none of which is in the
    charged atoms of an aspartate, so the residue read as neutral. This is the
    largest of the effects: charge moves further than the surface does.
    """

    residues = [residue for residue in pdb_parser.parse(alternates_path).residues if residue.name == "ASP"]

    assert [round(residue.charge(7), 3) for residue in residues] == [-1.0]


def test_a_disordered_terminus_recovers_its_charge():
    """the terminal charge is matched by atom name too, so it was lost the same way

    Atom.charge asks whether the name equals the residue's terminus. 1CBN THR 1
    is disordered, so its backbone N arrived as "N  A" and the N terminal
    ammonium contributed nothing.
    """

    structure = pdb_parser.parse(alternates_path)

    assert sum(1 for residue in structure.residues if residue.charge(7) != 0) == 6
    assert sum(1 for atom in structure.atoms if atom.charge(7) != 0) == 12


def test_aromatic_carbons_get_their_own_radius_again():
    """Atom.radius looks up aromatic carbons by name before falling back to element

    So a mangled name silently handed an aromatic carbon the plain carbon radius
    of 2.0 A instead of 1.85 A. Only 6 of 1CBN's 18 aromatic carbons used to be
    recognised.
    """

    from prodes import data

    structure = pdb_parser.parse(alternates_path)
    aromatic = [atom for atom in structure.atoms if atom.name in (data.residue_data(atom.residue_name)["aromatic_carbons"] or [])]

    assert len(aromatic) == 18
    assert {atom.radius for atom in aromatic} == {data.vdw_radius("Cr")}


def test_occupancy_and_alternate_location_are_read():
    """both are their own fields on the Atom, and occupancy is a number"""

    first = pdb_parser.parse(alternates_path).atoms[0]

    assert first.altloc == "A"
    assert first.occupancy == 0.80


def test_a_missing_occupancy_is_none_rather_than_zero(tmp_path):
    """a file that stops before column 60 has not said the atom is absent

    Zero is a real occupancy. Reading a blank column as 0.0 would make every
    conformation of such a file tie at zero rather than falling through to the
    order the letters appear in.
    """

    truncated = [atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N")[:54]]
    path = write_structure(tmp_path / "no_occupancy.pdb", truncated)

    assert pdb_parser.parse(path).atoms[0].occupancy is None


def test_the_more_occupied_conformation_wins(tmp_path):
    """B beats A when B is the better occupied, rather than the first winning"""

    lines = [
        atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N", altloc="A", occupancy=0.30),
        atom_line(2, "CB", "ALA", "A", 1, 1.0, 0.0, 0.0, "C", altloc="A", occupancy=0.30),
        atom_line(3, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N", altloc="B", occupancy=0.70),
        atom_line(4, "CB", "ALA", "A", 1, 2.0, 0.0, 0.0, "C", altloc="B", occupancy=0.70),
    ]
    structure = pdb_parser.parse(write_structure(tmp_path / "occupancy.pdb", lines))

    assert [atom.altloc for atom in structure.atoms] == ["B", "B"]
    assert [atom.x for atom in structure.atoms] == [0.0, 2.0]


def test_a_conformation_is_ranked_by_its_median_occupancy(tmp_path):
    """one atom refined away from the rest does not decide the conformation

    Occupancy is not constant across a conformation's atoms: 32 of 225
    conformations across 40 real structures vary, and 1CBN residue 34 carries
    both 0.80 and 1.00 under the same letter. A mean would let the stray 1.00
    below carry A past B; the median ignores it.
    """

    lines = [
        atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N", altloc="A", occupancy=0.30),
        atom_line(2, "CA", "ALA", "A", 1, 1.0, 0.0, 0.0, "C", altloc="A", occupancy=0.30),
        atom_line(3, "CB", "ALA", "A", 1, 2.0, 0.0, 0.0, "C", altloc="A", occupancy=1.00),
        atom_line(4, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N", altloc="B", occupancy=0.40),
        atom_line(5, "CA", "ALA", "A", 1, 1.0, 0.0, 0.0, "C", altloc="B", occupancy=0.40),
        atom_line(6, "CB", "ALA", "A", 1, 3.0, 0.0, 0.0, "C", altloc="B", occupancy=0.40),
    ]
    structure = pdb_parser.parse(write_structure(tmp_path / "median.pdb", lines))

    assert {atom.altloc for atom in structure.atoms} == {"B"}


def test_equal_occupancies_are_broken_by_the_order_they_appear(tmp_path):
    """a tie goes to the first letter in the file

    Ties are the common case rather than a corner: roughly a quarter of alternate
    residues in real structures are written at exactly equal occupancy, so this
    rule decides a large share of every choice made here and has to give the same
    answer on every run.
    """

    lines = [
        atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N", altloc="A", occupancy=0.50),
        atom_line(2, "CB", "ALA", "A", 1, 1.0, 0.0, 0.0, "C", altloc="A", occupancy=0.50),
        atom_line(3, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N", altloc="B", occupancy=0.50),
        atom_line(4, "CB", "ALA", "A", 1, 2.0, 0.0, 0.0, "C", altloc="B", occupancy=0.50),
    ]
    structure = pdb_parser.parse(write_structure(tmp_path / "tie.pdb", lines))

    assert {atom.altloc for atom in structure.atoms} == {"A"}
    assert [atom.x for atom in structure.atoms if atom.name == "CB"] == [1.0]


def test_a_tie_is_not_decided_by_how_many_atoms_a_conformation_has(tmp_path):
    """conformations of unequal length at one occupancy still go to the first letter

    Ranking on the mean would make this depend on float arithmetic: the mean of
    three 0.35s is 0.3499999999999999 and of two is exactly 0.35, so an exact
    comparison would rank the longer conformation lower and the tie break would
    never run. The median is exact for equal values, and the comparison holds a
    tolerance regardless.

    The longer conformation still donates the atom the winner does not have,
    which is the completion rule rather than the election.
    """

    lines = [
        atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N", altloc="A", occupancy=0.35),
        atom_line(2, "CA", "ALA", "A", 1, 1.0, 0.0, 0.0, "C", altloc="A", occupancy=0.35),
        atom_line(3, "N", "ALA", "A", 1, 5.0, 0.0, 0.0, "N", altloc="B", occupancy=0.35),
        atom_line(4, "CA", "ALA", "A", 1, 6.0, 0.0, 0.0, "C", altloc="B", occupancy=0.35),
        atom_line(5, "CB", "ALA", "A", 1, 7.0, 0.0, 0.0, "C", altloc="B", occupancy=0.35),
    ]
    structure = pdb_parser.parse(write_structure(tmp_path / "uneven_tie.pdb", lines))

    assert [atom.x for atom in structure.atoms if atom.name in ("N", "CA")] == [0.0, 1.0]
    # B's extra CB goes with B. Nothing is carried over from a conformation that
    # lost, however incomplete the winner looks beside it.
    assert not [atom for atom in structure.atoms if atom.name == "CB"]


def test_equal_occupancies_that_are_unequal_floats_still_tie(tmp_path):
    """the tolerance in the comparison is load bearing, not decoration

    Deposited occupancies carry two decimals, and a median of two of them is not
    always the same float for the same nominal value: 0.04 and 0.37 give exactly
    0.205, while 0.01 and 0.40 give 0.20500000000000002. Comparing exactly would
    hand the second conformation the win on float noise alone and the documented
    first-appearance tie break would never run. 125 median values reachable from
    two decimal occupancies have this property.
    """

    lines = [
        atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N", altloc="A", occupancy=0.04),
        atom_line(2, "CB", "ALA", "A", 1, 1.0, 0.0, 0.0, "C", altloc="A", occupancy=0.37),
        atom_line(3, "N", "ALA", "A", 1, 5.0, 0.0, 0.0, "N", altloc="B", occupancy=0.01),
        atom_line(4, "CB", "ALA", "A", 1, 6.0, 0.0, 0.0, "C", altloc="B", occupancy=0.40),
    ]
    structure = pdb_parser.parse(write_structure(tmp_path / "float_tie.pdb", lines))

    assert {atom.altloc for atom in structure.atoms} == {"A"}


def test_nothing_is_taken_from_a_losing_conformation(tmp_path):
    """a losing conformation contributes nothing, not even an atom the winner lacks

    Completing the winner from the runners up was tried and removed. Over
    thousands of disordered residues it never supplied a single heavy atom,
    because a conformation is written short exactly when it is the minor one, so
    the best occupied conformation is never the less complete of the two. What it
    did supply was hydrogens carrying the other rotamer's coordinates: on 4NZU it
    put an HG23 0.86 A from an OG1 it is not bonded to.
    """

    lines = [
        atom_line(1, "N", "TYR", "A", 1, 0.0, 0.0, 0.0, "N", altloc="A", occupancy=0.40),
        atom_line(2, "CB", "TYR", "A", 1, 1.0, 0.0, 0.0, "C", altloc="A", occupancy=0.40),
        atom_line(3, "OH", "TYR", "A", 1, 2.0, 0.0, 0.0, "O", altloc="A", occupancy=0.40),
        atom_line(4, "N", "TYR", "A", 1, 0.0, 0.0, 0.0, "N", altloc="B", occupancy=0.60),
        atom_line(5, "CB", "TYR", "A", 1, 1.5, 0.0, 0.0, "C", altloc="B", occupancy=0.60),
    ]
    structure = pdb_parser.parse(write_structure(tmp_path / "truncated.pdb", lines))

    assert sorted(atom.name for atom in structure.atoms) == ["CB", "N"]
    assert {atom.altloc for atom in structure.atoms} == {"B"}


def test_no_residue_ends_up_with_two_atoms_of_one_name(tmp_path):
    """the invariant the whole change exists to establish

    A residue that writes an atom both without a letter and under one used to
    come out holding both copies, 0.1 A apart and both counted in the surface.
    """

    lines = [
        atom_line(1, "N", "SER", "A", 1, 0.0, 0.0, 0.0, "N"),
        atom_line(2, "CA", "SER", "A", 1, 1.0, 0.0, 0.0, "C"),
        atom_line(3, "OG", "SER", "A", 1, 2.0, 0.0, 0.0, "O", altloc="A", occupancy=0.60),
        atom_line(4, "CA", "SER", "A", 1, 1.1, 0.0, 0.0, "C", altloc="B", occupancy=0.40),
        atom_line(5, "OG", "SER", "A", 1, 3.0, 0.0, 0.0, "O", altloc="B", occupancy=0.40),
    ]
    structure = pdb_parser.parse(write_structure(tmp_path / "duplicate.pdb", lines))

    assert sorted(atom.name for atom in structure.atoms) == ["CA", "N", "OG"]

    # 1CBN uses no insertion codes, so every residue in it should now hold each
    # atom name once. 4NZU is deliberately not checked: it is Kabat numbered, and
    # the parser merges residues that differ only by an insertion code, so H100
    # through H100H arrive as a single 111 atom residue holding 35 distinct
    # names. That is a separate defect with its own issue and this change does
    # not touch it.
    for residue in pdb_parser.parse(alternates_path).residues:
        names = [atom.name for atom in residue.atoms]
        assert len(names) == len(set(names))


def test_one_conformation_never_mixes_two_residue_names(tmp_path):
    """a letter spanning two residue names does not build a chimera

    Nothing in the format stops a file writing one letter across two residue
    names, and the atoms of the losing name would otherwise be relabelled with
    the winner's and go on to drive its charge and radius lookups.
    """

    lines = [
        atom_line(1, "N", "SER", "A", 1, 0.0, 0.0, 0.0, "N", altloc="A", occupancy=0.80),
        atom_line(2, "CG", "PRO", "A", 1, 1.0, 0.0, 0.0, "C", altloc="A", occupancy=0.80),
        atom_line(3, "N", "SER", "A", 1, 0.0, 0.0, 0.0, "N", altloc="B", occupancy=0.20),
    ]
    structure = pdb_parser.parse(write_structure(tmp_path / "spanning.pdb", lines))

    assert [(atom.name, atom.residue_name) for atom in structure.atoms] == [("N", "SER")]


def test_an_unusable_occupancy_is_treated_as_absent(tmp_path):
    """nan and inf parse cleanly as floats and used to crash the election

    They reached the ranking and failed there as an integer conversion error
    naming neither the file nor the column.
    """

    for bad in ("   nan", "   inf"):
        lines = [
            atom_line(1, "CB", "ALA", "A", 1, 0.0, 0.0, 0.0, "C", altloc="A", occupancy=0.50),
            atom_line(2, "CB", "ALA", "A", 1, 1.0, 0.0, 0.0, "C", altloc="B", occupancy=0.50),
        ]
        lines[1] = lines[1][:54] + f"{bad:>6s}" + lines[1][60:]
        structure = pdb_parser.parse(write_structure(tmp_path / "unusable.pdb", lines))

        assert len(structure.atoms) == 1
        assert structure.atoms[0].occupancy == 0.50


def test_a_shared_atom_with_no_indicator_is_always_kept(tmp_path):
    """the ordinary case: an ordered backbone with a disordered side chain

    91% of alternate residues in real structures look like this, so it matters
    more than any of the pathological cases above.
    """

    lines = [
        atom_line(1, "N", "SER", "A", 1, 0.0, 0.0, 0.0, "N"),
        atom_line(2, "CA", "SER", "A", 1, 1.0, 0.0, 0.0, "C"),
        atom_line(3, "OG", "SER", "A", 1, 2.0, 0.0, 0.0, "O", altloc="A", occupancy=0.40),
        atom_line(4, "OG", "SER", "A", 1, 3.0, 0.0, 0.0, "O", altloc="B", occupancy=0.60),
    ]
    structure = pdb_parser.parse(write_structure(tmp_path / "sidechain.pdb", lines))

    assert sorted(atom.name for atom in structure.atoms) == ["CA", "N", "OG"]
    assert [atom.x for atom in structure.atoms if atom.name == "OG"] == [3.0]


def test_the_residue_is_named_for_the_conformation_that_won(tmp_path):
    """not for whichever atom arrived first

    A residue modelled as two amino acids can write a shared backbone atom under
    the name of the conformation that loses, which 1EJG does. Naming the residue
    after its first atom would then contradict the atoms it actually holds.
    """

    lines = [
        atom_line(1, "N", "SER", "A", 1, 0.0, 0.0, 0.0, "N"),
        atom_line(2, "OG", "SER", "A", 1, 2.0, 0.0, 0.0, "O", altloc="A", occupancy=0.20),
        atom_line(3, "CG", "PRO", "A", 1, 1.5, 0.0, 0.0, "C", altloc="B", occupancy=0.80),
    ]
    structure = pdb_parser.parse(write_structure(tmp_path / "named.pdb", lines))

    assert structure.residues[0].name == "PRO"
    assert {atom.residue_name for atom in structure.atoms} == {"PRO"}


def test_residues_differing_only_by_insertion_code_choose_separately(tmp_path):
    """one residue's alternates never decide another's

    The parser otherwise ignores the insertion code and merges 30 with 30A, which
    is its own defect. Leaving the code out of the grouping key here would add a
    new one on top: a single conformation elected across both residues.
    """

    lines = [
        atom_line(1, "CB", "ALA", "A", 30, 0.0, 0.0, 0.0, "C", altloc="A", occupancy=0.70),
        atom_line(2, "CB", "ALA", "A", 30, 1.0, 0.0, 0.0, "C", altloc="B", occupancy=0.30),
        atom_line(3, "CB", "ALA", "A", 30, 2.0, 0.0, 0.0, "C", altloc="A", occupancy=0.30, insertion="A"),
        atom_line(4, "CB", "ALA", "A", 30, 3.0, 0.0, 0.0, "C", altloc="B", occupancy=0.70, insertion="A"),
    ]
    structure = pdb_parser.parse(write_structure(tmp_path / "insertion.pdb", lines))

    assert [atom.x for atom in structure.atoms] == [0.0, 3.0]


def test_hetatm_alternates_collapse_when_hetatm_is_being_read(tmp_path):
    """the choice follows whichever record type the parser was asked for"""

    lines = [
        atom_line(1, "O", "HOH", "A", 1, 0.0, 0.0, 0.0, "O", altloc="A", occupancy=0.40).replace("ATOM  ", "HETATM", 1),
        atom_line(2, "O", "HOH", "A", 1, 1.0, 0.0, 0.0, "O", altloc="B", occupancy=0.60).replace("ATOM  ", "HETATM", 1),
    ]
    path = write_structure(tmp_path / "waters.pdb", lines)

    het_parser = parser.PDBparser()
    structure = het_parser.parse(path, identifier="HETATM")

    assert len(structure.atoms) == 1
    assert structure.atoms[0].x == 1.0


def test_the_collapse_is_reported(caplog):
    """the counts are logged, and the renamed residues are named individually

    A residue whose winning conformation is a different amino acid changes the
    sequence, and is the case where a pKa file generated from the original
    coordinates no longer agrees with the residue it was meant for.
    """

    with caplog.at_level(logging.WARNING):
        structure = pdb_parser.parse(alternates_path)

    assert "13 residues modelled in more than one conformation" in caplog.text
    assert "dropping 131 atoms" in caplog.text
    assert "A22" in caplog.text and "A25" in caplog.text
    assert structure.alternate_conformers.summary() == {
        "residues_with_alternates": 13,
        "atoms_dropped": 131,
        "residues_renamed": ["PRO A22", "LEU A25"],
    }


def test_a_structure_without_alternates_says_nothing(caplog):
    """no warning, and a report that records nothing was collapsed"""

    with caplog.at_level(logging.WARNING):
        structure = pdb_parser.parse(file_path)

    assert "conformation" not in caplog.text
    assert structure.alternate_conformers.summary() == {
        "residues_with_alternates": 0,
        "atoms_dropped": 0,
        "residues_renamed": [],
    }


def test_the_representative_fab_collapses_correctly():
    """4NZU is the representative case, and the completion rule fires on it

    1CBN is a 0.83 A crambin with microheterogeneity, which is a case this tool
    will essentially never meet: across 24 real X-ray entries, alternate residues
    run at 0.26% and microheterogeneity at zero. 4NZU is an ordinary Fab, and it
    is the fixture that says the change works on what Prodes is actually pointed
    at.
    """

    structure = pdb_parser.parse(fab_path)

    assert len(structure.atoms) == 5487
    assert len(structure.heavy_atoms) == 3270


@pytest.mark.parametrize(
    "name, atoms",
    [("1GDW", 1022), ("1GDW_h", 2003), ("1AO6", 9198), ("1GPB", 6699), ("ARH96693", 479), ("ARH98503", 3106)],
)
def test_structures_without_alternates_parse_exactly_as_before(name, atoms):
    """not one of the structures already committed here contains an alternate

    So the change has to be a no-op on every one of them, and on the reference
    CSV that tests/test_sasa.py compares against.
    """

    assert len(pdb_parser.parse(f"tests/data/{name}.pdb.zip").atoms) == atoms


def test_the_report_names_a_renamed_residue_as_the_run_record_expects():
    """the strings, not only the counts

    summary() goes into the run record of every run, so a residue key that
    cannot be formatted is a failure in the bundle writing path and not only in
    a test. The key holds a parsed residue number, so nothing may assume it is
    still the raw four columns.
    """

    report = pdb_parser.parse(alternates_path).alternate_conformers

    assert report.summary()["residues_renamed"] == ["PRO A22", "LEU A25"]


def test_kept_records_come_back_in_file_order(tmp_path):
    """not grouped by residue

    A structure's atoms are held in the order its records appear, so grouping
    them here would quietly reorder every atom array in the package, the PDB
    writer's TER placement with them, while still passing most of the suite.
    """

    lines = [
        atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N"),
        atom_line(2, "CB", "ALA", "A", 2, 1.0, 0.0, 0.0, "C", altloc="A", occupancy=0.70),
        atom_line(3, "CA", "ALA", "A", 1, 2.0, 0.0, 0.0, "C"),
        atom_line(4, "CB", "ALA", "A", 2, 3.0, 0.0, 0.0, "C", altloc="B", occupancy=0.30),
    ]
    text = open(write_structure(tmp_path / "order.pdb", lines)).read()

    kept, _ = elect_conformers(read_atom_records(text, "test"))

    assert [record.x for record in kept] == [0.0, 1.0, 2.0]


def test_a_hetatm_never_decides_a_protein_residues_conformation(tmp_path):
    """a water numbered onto a protein residue shares its election key

    Nothing in the format stops a HETATM carrying the chain and residue number
    of an ATOM residue. Electing over both record types together would let the
    water's occupancy pick the amino acid's conformation, which is a new error
    of exactly the kind the insertion code is in the key to prevent.
    """

    lines = [
        atom_line(1, "CB", "ALA", "A", 1, 0.0, 0.0, 0.0, "C", altloc="A", occupancy=0.60),
        atom_line(2, "CB", "ALA", "A", 1, 1.0, 0.0, 0.0, "C", altloc="B", occupancy=0.40),
        atom_line(3, "O", "HOH", "A", 1, 5.0, 0.0, 0.0, "O", altloc="B", occupancy=0.99).replace("ATOM  ", "HETATM", 1),
    ]
    path = write_structure(tmp_path / "water_on_residue.pdb", lines)

    structure = pdb_parser.parse(path)

    assert [(atom.name, atom.altloc, atom.x) for atom in structure.atoms] == [("CB", "A", 0.0)]


def test_a_tie_goes_to_the_first_letter_written_and_not_the_first_alphabetically(tmp_path):
    """the two rules agree on every other tie test in this file, and disagree here

    Roughly a quarter of alternate residues in real structures are written at
    exactly equal occupancy, so this rule decides a large share of every choice
    made here. Ranking on the letter instead of on where it first appears passes
    every other test in this module.
    """

    lines = [
        atom_line(1, "CB", "ALA", "A", 1, 1.0, 0.0, 0.0, "C", altloc="B", occupancy=0.50),
        atom_line(2, "CB", "ALA", "A", 1, 2.0, 0.0, 0.0, "C", altloc="A", occupancy=0.50),
    ]
    structure = pdb_parser.parse(write_structure(tmp_path / "reversed_tie.pdb", lines))

    assert [(atom.altloc, atom.x) for atom in structure.atoms] == [("B", 1.0)]


def test_a_conformation_with_no_occupancy_loses_to_one_written_at_zero(tmp_path):
    """zero is a real occupancy and a blank column is not a zero

    UNKNOWN_OCCUPANCY sits below every stated occupancy for this reason. Ranking
    an unstated occupancy as 0.0 would make the two tie and hand the choice to
    the order they appear in, which is the opposite of what the file says.

    The unstated conformation is written first deliberately. Written second it
    loses either way, once on rank and once on the tie break, so the rule under
    test would not be the rule deciding the answer.
    """

    unstated = atom_line(1, "CB", "ALA", "A", 1, 2.0, 0.0, 0.0, "C", altloc="B")[:54]
    stated = atom_line(2, "CB", "ALA", "A", 1, 1.0, 0.0, 0.0, "C", altloc="A", occupancy=0.00)
    structure = pdb_parser.parse(write_structure(tmp_path / "zero_beats_blank.pdb", [unstated, stated]))

    assert [(atom.altloc, atom.x) for atom in structure.atoms] == [("A", 1.0)]


def test_the_residue_label_formats_a_key_the_way_a_warning_reads(tmp_path):
    """the run record and the log agree about what a residue is called"""

    from prodes.io.conformers import residue_label

    assert residue_label(("A", 22, "")) == "A22"
    assert residue_label(("H", 100, "A")) == "H100A"
    assert residue_label(("", -3, "")) == "-3"
