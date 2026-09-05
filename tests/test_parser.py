import numpy as np
import pytest

from prodes import data
from prodes.io import parser

file_path = "tests/data/1GDW.pdb.zip"
pdb_parser = parser.PDBparser()
pdb_parser.identifier = "ATOM"


def test_main_fail():
    """tests if a incorect file will give an error"""

    with pytest.raises(ValueError):
        pdb_parser.parse(file_path[:-2])


def test_main_accept():
    """tests if a file is parsed correctly"""

    global structure
    structure = pdb_parser.parse(file_path)

    assert len(structure.atoms) == 1022
    assert len(structure.chains) == 1
    assert len(structure.residues) == 130


def test_read_text_builds_the_same_structure_as_reading_the_file():
    """the text of a file parses to what the file itself parses to

    parse_pdb_text is the whole of the parse below the archive handling, so a
    caller holding the text of a structure never has to write it to disk first.
    """

    from_text = parser.parse_pdb_text(parser.read_pdb_text(file_path), "test")
    from_file = pdb_parser.parse(file_path)

    assert from_text.name == "test"
    assert len(from_text.atoms) == len(from_file.atoms)
    assert [(atom.name, atom.x, atom.y, atom.z) for atom in from_text.atoms] == [(atom.name, atom.x, atom.y, atom.z) for atom in from_file.atoms]


def test_structure_name_comes_from_the_file_stem():
    """the structure name, which becomes the output ID, is the bare file name

    Guards the Windows case: a backslash path used to leave the whole path in
    the name, which matters because zipped structures are parsed from a
    temporary directory rather than from tests/data.
    """

    assert pdb_parser.parse(file_path).name == "1GDW"


def test_structure_name_keeps_every_extension_but_the_last(tmp_path):
    """a name with more than one dot keeps all of it, since only the suffix is dropped

    1abc.ent.pdb gives 1abc.ent, not 1abc. Pinned because the name becomes the
    output ID, and because splitting on the first dot, which is what the code
    used to do, gave the other answer.
    """

    extracted = parser.extract_pdb(file_path, tmp_path)
    multi_dot = tmp_path / "1abc.ent.pdb"
    multi_dot.write_text(open(extracted).read())

    assert pdb_parser.parse(str(multi_dot)).name == "1abc.ent"


def test_zipped_structure_is_named_after_the_archive(tmp_path):
    """an archive names the structure, not the member file inside it

    A user pointing at bar.pdb.zip expects to get a row called bar, whatever the
    file inside happens to be called.
    """

    import zipfile

    extracted = parser.extract_pdb(file_path, tmp_path)
    archive = tmp_path / "renamed.pdb.zip"
    with zipfile.ZipFile(archive, "w") as zipped:
        zipped.write(extracted, arcname="something_else.pdb")

    assert pdb_parser.parse(str(archive)).name == "renamed"


def test_zipped_and_plain_parse_identically(tmp_path):
    """a zipped structure parses to the same atoms as the extracted plain file"""

    extracted = parser.extract_pdb(file_path, tmp_path)

    from_archive = pdb_parser.parse(file_path)
    from_plain = pdb_parser.parse(extracted)

    assert len(from_archive.atoms) == len(from_plain.atoms)
    assert [(a.name, a.x, a.y, a.z) for a in from_archive.atoms] == [(a.name, a.x, a.y, a.z) for a in from_plain.atoms]


def test_write_pdb_round_trips(tmp_path):
    """a structure written back out parses to the same atoms it came from

    write_pdb had no test and, once the old scratch scripts were removed, no
    caller either. This keeps its column layout honest: the parser reads by
    fixed column position, so any drift in the writer shows up here.
    """

    structure = pdb_parser.parse(file_path)
    written = tmp_path / "round_trip.pdb"

    parser.write_pdb(structure, str(written))
    reparsed = pdb_parser.parse(str(written))

    assert len(reparsed.atoms) == len(structure.atoms)
    assert [atom.name for atom in reparsed.atoms] == [atom.name for atom in structure.atoms]
    assert [atom.chain_name for atom in reparsed.atoms] == [atom.chain_name for atom in structure.atoms]
    assert [atom.element for atom in reparsed.atoms] == [atom.element for atom in structure.atoms]


def test_write_pdb_rejects_a_non_pdb_name(tmp_path, capsys):
    """asking for a file without a pdb extension writes nothing and says so"""

    target = tmp_path / "structure.txt"

    parser.write_pdb(pdb_parser.parse(file_path), str(target))

    assert "can only make files" in capsys.readouterr().out
    assert not target.exists()


def test_write_pdb_rejects_residue_numbers_that_do_not_fit(tmp_path):
    """a 5 digit residue number is refused rather than silently truncated

    Columns 23-26 are all the format gives the residue number, so a wider one
    would overflow into the chain field and read back as a different residue.
    """

    structure = pdb_parser.parse(file_path)
    structure.atoms[0].residue_number = 12345
    target = tmp_path / "overflow.pdb"

    with pytest.raises(ValueError, match="more than the 4 columns"):
        parser.write_pdb(structure, str(target))

    assert not target.exists()


def test_builder_makes_a_placeholder_atom():
    """the Builder produces a dummy atom at the requested position"""

    atom = parser.Builder().build_dummy_atom(1.5, -2.5, 3.5, chain_name="B")

    assert (atom.x, atom.y, atom.z) == (1.5, -2.5, 3.5)
    assert atom.chain_name == "B"
    assert atom.element == "X"


def test_archive_with_no_pdb_is_rejected(tmp_path):
    """an archive that does not hold exactly one pdb file raises a clear error"""

    import zipfile

    empty = tmp_path / "empty.pdb.zip"
    with zipfile.ZipFile(empty, "w") as archive:
        archive.writestr("readme.txt", "no structure here")

    with pytest.raises(ValueError, match="expected exactly one .pdb file"):
        pdb_parser.parse(str(empty))


def test_parsed_atom():
    """tests parsed atom"""

    atom = structure.atoms[0]

    assert atom.name == "N"
    assert atom.identifier == "ATOM"
    assert atom.chain_name == "A"
    assert atom.residue_number == 1
    assert atom.residue_name == "LYS"
    assert atom.x == 1.134
    assert atom.y == 19.824
    assert atom.z == 22.575
    assert atom.structure == structure


def test_parsed_chain():
    """tests parsed chain"""

    chain = structure.chains[0]
    assert chain.name == "A"


def test_parsed_residue():
    """tests parsed residue"""

    residue = structure.residues[0]

    assert residue.name == "LYS"
    assert residue.number == 1
    assert len(residue.atoms) == 9


def test_a_file_with_no_records_of_the_requested_type_says_so():
    """rather than failing on an empty array

    1GDW holds no HETATM records at all, and asking for them used to end in an
    IndexError from the terminus assignment, naming neither the file nor what
    was missing.
    """

    with pytest.raises(ValueError, match="holds no HETATM records"):
        pdb_parser.parse(file_path, identifier="HETATM")


def test_a_coordinate_record_after_an_end_record_is_refused(tmp_path):
    """the reader stops at END, so the atoms after one would be lost silently

    Biopython treats the six column END of the specification, and any CONECT, as
    the end of the coordinates and hands back everything after it unparsed. A
    concatenated file, or one that writes CONECT for a ligand before the rest of
    its atoms, reaches this. prodes reads such a file today, so the choice is
    between losing atoms quietly and saying so.

    Written with the padded form deliberately: a bare three character END, which
    is what write_pdb and tests.pdb_records produce, is not recognised and does
    not truncate anything.
    """

    from tests.pdb_records import atom_line, write_structure

    for terminator in ("END" + " " * 77, "CONECT    1    2"):
        lines = [
            atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N"),
            terminator,
            atom_line(2, "CA", "ALA", "A", 2, 1.0, 0.0, 0.0, "C"),
        ]
        path = write_structure(tmp_path / "truncated.pdb", lines)

        with pytest.raises(ValueError, match="written after an END or CONECT record"):
            pdb_parser.parse(path)


def test_a_bare_end_record_does_not_truncate_anything(tmp_path):
    """which is what makes every structure written by write_pdb still readable"""

    from tests.pdb_records import atom_line, write_structure

    lines = [
        atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N"),
        "END",
        atom_line(2, "CA", "ALA", "A", 2, 1.0, 0.0, 0.0, "C"),
    ]
    path = write_structure(tmp_path / "bare_end.pdb", lines)

    assert len(pdb_parser.parse(path).atoms) == 2


def test_a_coordinate_written_to_more_decimals_than_the_format_holds_is_reported(tmp_path, caplog):
    """the eight column field is three decimals and a fourth is lost

    A fair trade for a file that does not conform, and a bad one to make in
    silence, because nothing else in the parse would ever mention it.
    """

    import logging

    from tests.pdb_records import pdb_line, write_structure

    line = pdb_line(
        [
            (1, "ATOM"),
            (7, "    1"),
            (13, " N  "),
            (18, "ALA"),
            (22, "A"),
            (23, "   1"),
            (31, " 1.12345"),
            (39, "   0.000"),
            (47, "   0.000"),
            (55, "  1.00"),
            (77, " N"),
        ]
    )
    path = write_structure(tmp_path / "over_precise.pdb", [line])

    with caplog.at_level(logging.WARNING):
        structure = pdb_parser.parse(path)

    assert structure.atoms[0].x == 1.123
    assert "a precision the PDB format does not hold" in caplog.text


def test_a_coordinate_that_loses_nothing_is_not_reported(tmp_path, caplog):
    """a trailing zero and an exponent are read exactly, so they must not warn

    Asking whether a column has a fourth decimal would flag both. Asking whether
    rounding changes the number is the question that was meant.
    """

    import logging

    from tests.pdb_records import pdb_line, write_structure

    line = pdb_line(
        [
            (1, "ATOM"),
            (7, "    1"),
            (13, " N  "),
            (18, "ALA"),
            (22, "A"),
            (23, "   1"),
            (31, "  1.1230"),
            (39, " 1.0e+02"),
            (47, "   0.000"),
            (55, "  1.00"),
            (61, "  0.00"),
            (77, " N"),
        ]
    )
    path = write_structure(tmp_path / "exact.pdb", [line])

    with caplog.at_level(logging.WARNING):
        structure = parser.PDBparser().parse(path)

    assert (structure.atoms[0].x, structure.atoms[0].y) == (1.123, 100.0)
    assert "precision" not in caplog.text


def test_only_the_first_model_of_an_ensemble_is_described(caplog):
    """an ensemble read whole is not a protein, so one model of it is described

    Every model used to be read into one structure, and because a residue number
    that reappears does not start a new residue, the later models piled onto the
    last residue made: 1PIT came out as 58 residues holding 17 780 atoms, 16 901
    of them in one residue, at a formal charge of -1096.
    """

    import logging

    from tests.pdb_records import atom_line

    lines = []
    for model in (1, 2, 3):
        lines += [f"MODEL     {model:>4d}", atom_line(1, "N", "ALA", "A", 1, float(model), 0.0, 0.0, "N"), "ENDMDL"]

    with caplog.at_level(logging.WARNING):
        structure = parser.parse_pdb_text("\n".join(lines) + "\nEND\n", "ensemble")

    assert [atom.x for atom in structure.atoms] == [1.0]
    assert "holds 3 models" in caplog.text
    assert "ignoring the 2 records of the rest" in caplog.text


def test_an_ensemble_gives_the_structure_its_first_model_alone_would(tmp_path):
    """atom for atom, coordinates included, which is the whole claim of the fix"""

    from tests.pdb_records import atom_line, write_structure

    def alanine(x):
        return [
            atom_line(1, "N", "ALA", "A", 1, x, 0.0, 0.0, "N"),
            atom_line(2, "CA", "ALA", "A", 1, x + 1.0, 0.0, 0.0, "C"),
            atom_line(3, "CB", "ALA", "A", 1, x + 2.0, 0.0, 0.0, "C"),
        ]

    ensemble = []
    for model, x in enumerate((0.0, 10.0, 20.0), start=1):
        ensemble += [f"MODEL     {model:>4d}", *alanine(x), "ENDMDL"]

    from_ensemble = pdb_parser.parse(write_structure(tmp_path / "ensemble.pdb", ensemble))
    from_one_model = pdb_parser.parse(write_structure(tmp_path / "single.pdb", alanine(0.0)))

    assert [(atom.name, atom.x, atom.y, atom.z) for atom in from_ensemble.atoms] == [(atom.name, atom.x, atom.y, atom.z) for atom in from_one_model.atoms]
    assert [(residue.name, residue.number, len(residue.atoms)) for residue in from_ensemble.residues] == [("ALA", 1, 3)]


def test_the_model_is_chosen_before_the_conformation_is(tmp_path):
    """or a discarded model's occupancies decide which conformation of model 1 is kept

    The election keys residues on chain, number and insertion code and carries
    no model, so every model's letters go into one ballot if they reach it
    together. Here model 1 writes A at 0.60 against B at 0.40 and model 2 the
    other way round by a wider margin: electing first keeps B, which is not
    even model 1's better conformation and sits at model 2's coordinates.
    """

    from tests.pdb_records import atom_line, write_structure

    def residue(x, a_occupancy, b_occupancy):
        return [
            atom_line(1, "N", "SER", "A", 1, x, 0.0, 0.0, "N"),
            atom_line(2, "OG", "SER", "A", 1, x + 1.0, 0.0, 0.0, "O", altloc="A", occupancy=a_occupancy),
            atom_line(3, "OG", "SER", "A", 1, x + 2.0, 0.0, 0.0, "O", altloc="B", occupancy=b_occupancy),
        ]

    lines = ["MODEL        1", *residue(0.0, 0.60, 0.40), "ENDMDL", "MODEL        2", *residue(10.0, 0.10, 0.90), "ENDMDL"]
    structure = pdb_parser.parse(write_structure(tmp_path / "ordering.pdb", lines))

    assert [(atom.name, atom.altloc, atom.x) for atom in structure.atoms] == [("N", "", 0.0), ("OG", "A", 1.0)]


def test_the_run_record_gets_the_model_count_from_the_parse():
    """the count has to survive the whole way, or every ensemble bundle claims to be one structure

    The bundle ships the input file unchanged, so a record saying one model when
    the file holds twenty is the one statement nothing else in the bundle can
    contradict.
    """

    from prodes.output import run_metadata

    structure = pdb_parser.parse("tests/data/1PIT.pdb.zip")
    record = run_metadata("1PIT.pdb.zip", {}, 0, np.array([]), models=structure.models)

    assert structure.models == 20
    assert record["models_in_file"] == 20


def test_a_single_model_structure_records_one_model():
    """the ordinary case, which is what makes the count above worth reading"""

    assert pdb_parser.parse(file_path).models == 1


def test_the_first_model_holding_records_of_the_requested_type_is_the_one_kept(tmp_path):
    """not model 1 outright, or a file whose ligand appears later holds nothing to describe"""

    from tests.pdb_records import atom_line, pdb_line, write_structure

    def water(serial, model):
        return pdb_line(
            [
                (1, "HETATM"),
                (7, f"{serial:5d}"),
                (13, " O  "),
                (18, "HOH"),
                (22, "A"),
                (23, f"{model:4d}"),
                (31, f"{float(model):8.3f}"),
                (39, "   0.000"),
                (47, "   0.000"),
                (55, "  1.00"),
                (61, "  0.00"),
                (77, " O"),
            ]
        )

    lines = [
        "MODEL        1",
        atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N"),
        "ENDMDL",
        "MODEL        2",
        atom_line(2, "N", "ALA", "A", 1, 1.0, 0.0, 0.0, "N"),
        water(3, 2),
        "ENDMDL",
        "MODEL        3",
        atom_line(4, "N", "ALA", "A", 1, 2.0, 0.0, 0.0, "N"),
        water(5, 3),
        "ENDMDL",
    ]
    path = write_structure(tmp_path / "late_ligand.pdb", lines)

    assert [atom.x for atom in pdb_parser.parse(path).atoms] == [0.0]
    assert [atom.x for atom in pdb_parser.parse(path, identifier="HETATM").atoms] == [2.0]


def test_a_single_model_structure_says_nothing_about_models(caplog):
    """the ordinary case has to stay quiet, or the warning is worthless"""

    import logging

    with caplog.at_level(logging.WARNING):
        pdb_parser.parse(file_path)

    assert "models" not in caplog.text


def test_an_ssbond_written_after_the_coordinates_is_still_read(caplog):
    """it is read from the text, so where it sits in the file does not matter

    Biopython does not parse SSBOND anywhere, and its header stops at the first
    coordinate record, so a bond written at the end would be lost by anything
    that took the reader's word for what the file contains. A lost SSBOND is a
    charge change and not a cosmetic one.
    """

    import logging

    from tests.pdb_records import cysteine_lines

    lines = cysteine_lines(1, "A", 1, 0.0) + cysteine_lines(7, "A", 2, 2.0) + ["SSBOND   1 CYS A    1    CYS A    2"]

    with caplog.at_level(logging.WARNING):
        structure = parser.parse_pdb_text("\n".join(lines) + "\nEND\n", "trailing_ssbond")

    assert len(structure.disulfides) == 1
    assert "SSBOND" not in caplog.text


# What a residue is, which from version 8.0 is its chain, its number and its
# insertion code. Until then the code was left out and H100 and H100A were read
# as one residue carrying both side chains.


def test_the_fab_keeps_its_kabat_insertions_apart():
    """4NZU is the real case, and it is an ordinary therapeutic antibody

    Three of its residue numbers carry insertion codes, H100 running to H100H,
    and reading each set as one residue lost 12 residues of the 434 the file
    holds and 1.3 kDa of its mass. Kabat and Chothia numbering write every CDR
    insertion this way, so this is the common case for the structures prodes is
    aimed at rather than a corner of the format.
    """

    structure = pdb_parser.parse("tests/data/4NZU.pdb.zip")
    insertions = [residue.label for residue in structure.residues if residue.insertion_code]

    assert len(structure.residues) == 434
    assert structure.mw == pytest.approx(46409.70)
    assert insertions == ["H52A", "H82A", "H82B", "H82C", "H100A", "H100B", "H100C", "H100D", "H100E", "H100F", "H100G", "H100H"]


def test_the_nmr_ensemble_is_described_by_its_first_model():
    """1PIT is bovine pancreatic trypsin inhibitor, 20 models of 58 residues

    Read whole it came out as 58 residues holding 17 780 atoms, 16 901 of them
    piled into one residue, with a formal charge of -1096 at pH 7 for a protein
    whose real charge there is about +6. One model gives 889 atoms and that +6.
    """

    structure = pdb_parser.parse("tests/data/1PIT.pdb.zip")

    assert len(structure.residues) == 58
    assert len(structure.atoms) == 889
    assert max(len(residue.atoms) for residue in structure.residues) == 24
    assert structure.charge(7) == pytest.approx(6.0)


def test_residues_differing_only_by_an_insertion_code_are_separate_residues(tmp_path):
    """H100 and H100A are two residues, and reading them as one is wrong in four ways at once

    The merged residue had one name for two residues, so its mass and its
    surface normalisation were wrong; one set of pKas, so every charged atom of
    the second residue titrated against the first residue's group; and one
    residue name, so a cysteine hidden inside a residue named something else was
    invisible to the disulfide detection.
    """

    from tests.pdb_records import atom_line, write_structure

    lines = [
        atom_line(1, "N", "ASP", "H", 100, 0.0, 0.0, 0.0, "N"),
        atom_line(2, "CA", "ASP", "H", 100, 1.0, 0.0, 0.0, "C"),
        atom_line(3, "N", "LYS", "H", 100, 2.0, 0.0, 0.0, "N", insertion="A"),
        atom_line(4, "NZ", "LYS", "H", 100, 3.0, 0.0, 0.0, "N", insertion="A"),
    ]
    structure = pdb_parser.parse(write_structure(tmp_path / "insertion.pdb", lines))

    assert [(residue.name, residue.number, residue.insertion_code, len(residue.atoms)) for residue in structure.residues] == [
        ("ASP", 100, "", 2),
        ("LYS", 100, "A", 2),
    ]
    assert [residue.label for residue in structure.residues] == ["H100", "H100A"]


def test_an_insertion_takes_its_charge_from_its_own_pka(tmp_path):
    """the merge did not only mislabel the residue, it titrated its atoms against the wrong group

    An atom keeps the residue name from its own record, so a lysine NZ merged
    into a residue named ASP was still recognised as chargeable and then given
    aspartate's pKa of 3.86: at pH 7 it lost its +1 entirely. The equivalent
    cysteine came out at a full -1.
    """

    from tests.pdb_records import atom_line, write_structure

    lines = [
        atom_line(1, "N", "ASP", "H", 100, 0.0, 0.0, 0.0, "N"),
        atom_line(2, "OD1", "ASP", "H", 100, 1.0, 0.0, 0.0, "O"),
        atom_line(3, "OD2", "ASP", "H", 100, 2.0, 0.0, 0.0, "O"),
        atom_line(4, "N", "LYS", "H", 100, 3.0, 0.0, 0.0, "N", insertion="A"),
        atom_line(5, "NZ", "LYS", "H", 100, 4.0, 0.0, 0.0, "N", insertion="A"),
    ]
    structure = pdb_parser.parse(write_structure(tmp_path / "insertion_charge.pdb", lines))

    lysine = structure.residues[1]

    assert lysine.side_chain_pka == data.residue_data("LYS")["pka"]
    assert lysine.charge(7) == pytest.approx(1.0)


def test_write_pdb_round_trips_an_insertion_code(tmp_path):
    """column 27 is part of which residue an atom belongs to

    A writer that dropped it would hand back a file whose two residues are one
    again, which is the defect this release fixes arriving by another route.
    """

    from tests.pdb_records import atom_line, write_structure

    lines = [
        atom_line(1, "N", "ASP", "H", 100, 0.0, 0.0, 0.0, "N"),
        atom_line(2, "N", "ALA", "H", 100, 1.0, 0.0, 0.0, "N", insertion="A"),
        atom_line(3, "N", "GLY", "H", 100, 2.0, 0.0, 0.0, "N", insertion="B"),
    ]
    structure = pdb_parser.parse(write_structure(tmp_path / "insertion_source.pdb", lines))

    written = tmp_path / "insertion_round_trip.pdb"
    parser.write_pdb(structure, str(written))
    reparsed = pdb_parser.parse(str(written))

    assert [(residue.name, residue.number, residue.insertion_code) for residue in reparsed.residues] == [
        (residue.name, residue.number, residue.insertion_code) for residue in structure.residues
    ]


def test_write_pdb_leaves_column_27_blank_when_there_is_no_insertion_code(tmp_path):
    """which is every atom of almost every structure, and the dummy surface points

    Checked on the columns rather than by reparsing, because a code written into
    the wrong column would shift the coordinates and be caught somewhere else,
    while a stray character in the right one would not be caught at all.
    """

    structure = pdb_parser.parse(file_path)
    written = tmp_path / "no_insertion.pdb"

    parser.write_pdb(structure, str(written))
    coordinates = [line for line in written.read_text().splitlines() if line.startswith(("ATOM", "HETATM"))]

    assert coordinates
    assert {line[26] for line in coordinates} == {" "}


# The grouping rules build_structure has to reproduce. Both are quirks rather
# than intentions and no shipped structure reaches either. Without these tests
# the next person to edit build_structure would be told by a passing suite that
# they are free to choose, so each was checked against the parser as it was
# before the reader changed and pins that answer.


def test_a_chain_that_reappears_is_the_chain_it_already_was(tmp_path):
    """not a second chain of the same name

    Which also decides where the C terminus goes, and the answer is not "the
    last residue of each chain". It is written on the chain being left at the
    moment a new chain name appears, so chain A's first residue takes a C the
    instant chain B starts, overwriting the N it had as the first of its chain;
    chain B never gets one, because nothing follows it; and the structure's last
    residue takes one at the end.
    """

    from tests.pdb_records import atom_line, write_structure

    lines = [
        atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N"),
        atom_line(2, "N", "GLY", "B", 1, 1.0, 0.0, 0.0, "N"),
        atom_line(3, "N", "SER", "A", 2, 2.0, 0.0, 0.0, "N"),
    ]
    structure = pdb_parser.parse(write_structure(tmp_path / "reappearing_chain.pdb", lines))

    assert [chain.name for chain in structure.chains] == ["A", "B"]
    assert [len(chain.residues) for chain in structure.chains] == [2, 1]
    assert [(residue.name, residue.terminus) for residue in structure.residues] == [("ALA", "C"), ("GLY", "N"), ("SER", "C")]


def test_a_residue_key_that_reappears_does_not_start_a_new_residue(tmp_path):
    """its atoms join whichever residue was made most recently

    A quirk and not a decision: the parser tests membership in the keys already
    seen in the chain, so a key that comes back after another has intervened
    finds itself already known. It was the mechanism behind the damage an NMR
    ensemble took, and selecting one model has put that out of reach rather than
    changing this, so a file that writes one residue's atoms in two places still
    reaches it.
    """

    from tests.pdb_records import atom_line, write_structure

    lines = [
        atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N"),
        atom_line(2, "N", "GLY", "A", 2, 1.0, 0.0, 0.0, "N"),
        atom_line(3, "CA", "ALA", "A", 1, 2.0, 0.0, 0.0, "C"),
    ]
    structure = pdb_parser.parse(write_structure(tmp_path / "reappearing_number.pdb", lines))

    assert [(residue.number, residue.name, len(residue.atoms)) for residue in structure.residues] == [(1, "ALA", 1), (2, "GLY", 2)]


def test_a_one_residue_chain_ends_up_a_c_terminus(tmp_path):
    """the C written at the chain transition overwrites the N written on the same residue

    Which is what the terminal charge is then read from, so it is not cosmetic.
    """

    from tests.pdb_records import atom_line, write_structure

    lines = [
        atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N"),
        atom_line(2, "N", "GLY", "B", 1, 1.0, 0.0, 0.0, "N"),
    ]
    structure = pdb_parser.parse(write_structure(tmp_path / "one_residue_chain.pdb", lines))

    assert [residue.terminus for residue in structure.residues] == ["C", "C"]


def test_a_coordinate_that_is_not_a_number_is_reported(tmp_path, caplog):
    """nan and inf parse cleanly and then poison every distance measured from them

    Read as they always were, because refusing them would be a change to what
    parses. Said out loud, because nothing downstream will mention it: a nan
    coordinate propagates through the surface and the grids as more nans.
    """

    import logging

    from tests.pdb_records import pdb_line, write_structure

    line = pdb_line(
        [
            (1, "ATOM"),
            (7, "    1"),
            (13, " N  "),
            (18, "ALA"),
            (22, "A"),
            (23, "   1"),
            (31, "     nan"),
            (39, "   0.000"),
            (47, "   0.000"),
            (55, "  1.00"),
            (61, "  0.00"),
            (77, " N"),
        ]
    )
    path = write_structure(tmp_path / "unusable.pdb", [line])

    with caplog.at_level(logging.WARNING):
        parser.PDBparser().parse(path)

    assert "not numbers" in caplog.text
    assert "precision" not in caplog.text
