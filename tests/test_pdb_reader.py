"""Tests for the Biopython-backed record reader.

The reader is the layer that replaced the hand-written column arithmetic. What
has to hold is that it returns one record per coordinate line, in file order,
carrying that line's own columns and nothing inferred, because everything above
it — the election, the residue grouping, every feature — is built on that.
"""

import logging

import pytest

from prodes.io.parser import read_pdb_text
from prodes.io.pdb_reader import raised_inside_prodes, read_atom_records, readable_occupancy
from tests.pdb_records import atom_line, write_structure

STRUCTURES = ["1GDW", "1GDW_h", "1CBN", "1AO6", "1GPB", "4NZU", "ARH96693", "ARH98503"]


def coordinate_lines(text):
    """Returns the ATOM and HETATM lines of a PDB file, matched on the six columns."""

    return [line for line in text.splitlines() if line[0:6] in ("ATOM  ", "HETATM")]


@pytest.mark.parametrize("name", STRUCTURES)
def test_one_record_per_coordinate_line_in_file_order(name):
    """no record is dropped, none is invented, and the order is the file's

    Order is not decoration: a structure's atoms are held in the order its
    records appear, and the PDB writer, the TER placement and every atom-level
    comparison depend on it.
    """

    text = read_pdb_text(f"tests/data/{name}.pdb.zip")
    lines = coordinate_lines(text)

    records = read_atom_records(text, name)

    assert len(records) == len(lines)
    for record, line in zip(records, lines, strict=True):
        assert record.identifier == line[0:6].strip()
        assert record.name == line[12:16].strip()
        assert record.altloc == line[16].strip()
        assert record.residue_name == line[17:20].strip()
        assert record.chain_name == line[21].strip()
        assert record.residue_number == int(line[22:26])
        assert record.insertion_code == line[26].strip()
        assert record.element == line[76:78].strip()
        assert record.segment_id == line[72:76].strip()
        assert record.occupancy == float(line[54:60])


@pytest.mark.parametrize("name", STRUCTURES)
def test_coordinates_are_the_columns_exactly(name):
    """not approximately, and not to a tolerance

    Biopython hands coordinates over as float32, which would perturb every one
    of them in the seventh significant digit. They are rounded back to the three
    decimals the format holds, which recovers the file's own number exactly for
    every value the eight column field can carry.
    """

    text = read_pdb_text(f"tests/data/{name}.pdb.zip")

    for record, line in zip(read_atom_records(text, name), coordinate_lines(text), strict=True):
        assert (record.x, record.y, record.z) == (float(line[30:38]), float(line[38:46]), float(line[46:54]))


def test_a_missing_occupancy_is_none_rather_than_zero(tmp_path):
    """a file that stops before column 60 has not said the atom is absent

    Zero is a real occupancy. Reading a blank column as 0.0 would make every
    conformation of such a file tie at zero rather than falling through to the
    order the letters appear in.
    """

    truncated = [atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N")[:54]]
    text = write_structure(tmp_path / "no_occupancy.pdb", truncated)

    assert read_atom_records(open(text).read(), "test")[0].occupancy is None


def test_a_missing_occupancy_survives_an_alternate_location(tmp_path):
    """the case that crashes Biopython's own entity builder

    DisorderedAtom.disordered_add compares the occupancy against the best seen
    so far, which raises TypeError when the column is blank. Nothing here builds
    a DisorderedAtom, so the blank column stays a blank column.
    """

    lines = [
        atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N", altloc="A")[:54],
        atom_line(2, "N", "ALA", "A", 1, 1.0, 0.0, 0.0, "N", altloc="B")[:54],
    ]
    text = open(write_structure(tmp_path / "no_occupancy_alt.pdb", lines)).read()

    records = read_atom_records(text, "test")

    assert [record.occupancy for record in records] == [None, None]
    assert [record.altloc for record in records] == ["A", "B"]


@pytest.mark.parametrize("unusable", ["   nan", "   inf", "  -inf"])
def test_an_unusable_occupancy_is_treated_as_absent(tmp_path, unusable):
    """nan and inf parse cleanly as floats and then poison every comparison

    They used to reach the ranking and fail there as an integer conversion error
    naming neither the file nor the column.
    """

    line = atom_line(1, "CB", "ALA", "A", 1, 0.0, 0.0, 0.0, "C", altloc="A")
    text = open(write_structure(tmp_path / "unusable.pdb", [line[:54] + unusable + line[60:]])).read()

    assert read_atom_records(text, "test")[0].occupancy is None


def test_a_blank_element_column_is_not_guessed_at(tmp_path):
    """Biopython's Atom class infers the element from the atom name; this does not

    Inferring it would be an improvement and is its own issue: with a blank
    element every atom raises from the van der Waals radius lookup, so a file
    with no element column produces no features at all. Improving it here would
    change numbers in a release that promises not to.
    """

    lines = [atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "")[:76], atom_line(2, "HB1", "ALA", "A", 1, 1.0, 0.0, 0.0, "")[:76]]
    text = open(write_structure(tmp_path / "no_element.pdb", lines)).read()

    assert [record.element for record in read_atom_records(text, "test")] == ["", ""]


def test_an_element_symbol_is_uppercased(tmp_path):
    """which the old reader did not do

    The format specifies a right justified uppercase symbol and no shipped
    structure writes anything else, so nothing real moves. It is pinned because
    it is a difference, and because prodes' own pseudo-element for aromatic
    carbons, Cr, would not survive the trip if a file ever wrote it.
    """

    lines = [atom_line(1, "CB", "ALA", "A", 1, 0.0, 0.0, 0.0, " h"), atom_line(2, "FE", "ALA", "A", 1, 1.0, 0.0, 0.0, "Fe")]
    text = open(write_structure(tmp_path / "element_case.pdb", lines)).read()

    assert [record.element for record in read_atom_records(text, "test")] == ["H", "FE"]


def test_an_atom_name_keeps_no_internal_spaces(tmp_path):
    """columns 13-16 stripped, which is what the package looks names up by

    Biopython leaves a name with internal spaces unstripped, and a name like
    "N B" matches nothing: not the charged atoms of its residue, not the
    aromatic carbons that carry their own radius, not the SG a disulfide is
    found from.
    """

    from tests.pdb_records import pdb_line

    line = pdb_line(
        [(1, "ATOM"), (7, "    1"), (13, " N B"), (18, "ALA"), (22, "A"), (23, "   1"), (31, "   0.000"), (39, "   0.000"), (47, "   0.000"), (77, " N")]
    )
    text = open(write_structure(tmp_path / "spaced_name.pdb", [line])).read()

    assert read_atom_records(text, "test")[0].name == "N B"


def test_the_hetero_flag_becomes_the_record_type(tmp_path):
    """a water and a ligand are HETATM records; a protein atom is an ATOM record"""

    lines = [
        atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N"),
        atom_line(2, "O", "HOH", "A", 2, 1.0, 0.0, 0.0, "O").replace("ATOM  ", "HETATM", 1),
        atom_line(3, "C1", "GOL", "A", 3, 2.0, 0.0, 0.0, "C").replace("ATOM  ", "HETATM", 1),
    ]
    text = open(write_structure(tmp_path / "hetero.pdb", lines)).read()

    assert [record.identifier for record in read_atom_records(text, "test")] == ["ATOM", "HETATM", "HETATM"]


def test_insertion_codes_are_read_as_their_own_field(tmp_path):
    """which the old reader never carried on an atom at all

    Not used to group residues yet: 4NZU's Kabat insertions still merge, which
    is its own issue. Carried so that fixing it is a change to one grouping key.
    """

    lines = [
        atom_line(1, "N", "ALA", "H", 100, 0.0, 0.0, 0.0, "N"),
        atom_line(2, "N", "GLY", "H", 100, 1.0, 0.0, 0.0, "N", insertion="A"),
    ]
    text = open(write_structure(tmp_path / "insertion.pdb", lines)).read()

    records = read_atom_records(text, "test")

    assert [(record.residue_number, record.insertion_code) for record in records] == [(100, ""), (100, "A")]


def test_every_model_is_read_and_tagged_with_its_own_number(tmp_path):
    """an NMR ensemble is still read whole, which is the behaviour that has to hold

    prodes has no notion of a model and merges them all, which is its own issue.
    The number is carried so that fixing it is a filter rather than a rewrite.
    It counts models rather than repeating the serial the MODEL record carries,
    because nothing stops a file writing one serial twice and a key that can
    collide is not a key.
    """

    lines = [
        "MODEL        1",
        atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N"),
        "ENDMDL",
        "MODEL        2",
        atom_line(2, "N", "ALA", "A", 1, 9.0, 0.0, 0.0, "N"),
        "ENDMDL",
    ]
    text = open(write_structure(tmp_path / "ensemble.pdb", lines)).read()

    records = read_atom_records(text, "test")

    assert [record.model for record in records] == [0, 1]
    assert [record.x for record in records] == [0.0, 9.0]


def test_a_structure_with_no_model_records_reads_as_one_model(tmp_path):
    """the ordinary case, where the parser's own counter supplies the number"""

    text = open(write_structure(tmp_path / "single.pdb", [atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N")])).read()

    assert read_atom_records(text, "test")[0].model == 0


def test_anisou_records_produce_no_atoms():
    """4NZU carries 6190 of them, and prodes uses none

    They are skipped rather than read. The builder has to say so explicitly:
    Biopython hands anisotropic factors to the atom it last built, and this
    builder never builds one, so an unguarded reader dies on the first record.
    """

    text = read_pdb_text("tests/data/4NZU.pdb.zip")

    assert any(line.startswith("ANISOU") for line in text.splitlines())
    assert len(read_atom_records(text, "4NZU")) == len(coordinate_lines(text))


def test_a_malformed_anisou_record_names_the_file(tmp_path):
    """a new failure on a record type prodes does not use

    Biopython reads those columns with a bare float() before the builder is
    reached, so the record cannot simply be ignored the way it used to be. What
    can be helped is the error, which otherwise names neither the file nor the
    column.
    """

    lines = [atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N"), "ANISOU    1  N   ALA A   1"]
    text = open(write_structure(tmp_path / "anisou.pdb", lines)).read()

    with pytest.raises(ValueError, match="could not read the coordinate records of anisou"):
        read_atom_records(text, "anisou")


@pytest.mark.parametrize(
    "description, line",
    [
        ("truncated before the insertion code", atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N")[:26]),
        (
            "blank residue number",
            atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N")[:22] + "    " + atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N")[26:],
        ),
        (
            "unreadable coordinates",
            atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N")[:30] + "  ??????" + atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N")[38:],
        ),
    ],
)
def test_an_unreadable_record_raises_a_value_error_naming_the_file(tmp_path, description, line):
    """three of Biopython's own column reads sit outside its exception handling

    Left alone they escape as a bare IndexError or ValueError from inside a
    third party parser, naming neither the file nor what went wrong.
    """

    text = open(write_structure(tmp_path / "broken.pdb", [line])).read()

    with pytest.raises(ValueError, match="could not read the coordinate records of broken"):
        read_atom_records(text, "broken")


def test_an_empty_file_names_itself(tmp_path):
    """rather than raising Biopython's bare "Empty file." """

    with pytest.raises(ValueError, match="could not read the coordinate records of empty"):
        read_atom_records("", "empty")


def test_the_readers_complaints_go_to_the_log(tmp_path, caplog):
    """not to stderr, and not two per atom

    Biopython warns once for an unreadable occupancy column and once for an
    unreadable temperature factor, so a file trimmed at column 54 would put two
    warnings per atom on the terminal. They are collected and counted instead.
    """

    lines = [atom_line(serial, "N", "ALA", "A", serial, 0.0, 0.0, 0.0, "N")[:54] for serial in range(1, 6)]
    text = open(write_structure(tmp_path / "trimmed.pdb", lines)).read()

    with caplog.at_level(logging.WARNING):
        records = read_atom_records(text, "trimmed")

    assert len(records) == 5
    assert "occupancy" in caplog.text.lower()


def test_a_bare_end_record_is_not_complained_about(tmp_path, caplog):
    """prodes writes one itself, in write_pdb, and 1GDW_h carries one

    Biopython recognises only the six column form of the specification, so a
    bare END is reported as an unrecognised record. Passing that on would warn
    about a marker prodes produced and correctly ignored.
    """

    with caplog.at_level(logging.WARNING):
        read_atom_records(read_pdb_text("tests/data/1GDW_h.pdb.zip"), "1GDW_h")

    assert "END" not in caplog.text


def test_a_header_in_another_language_does_not_stop_the_parse():
    """the header is dropped rather than parsed, and this is why

    Biopython's header parser looks for a date anywhere in a HEADER or REVDAT
    record and looks the month up in a list of English abbreviations, so a file
    stamped in German raises ValueError: 'Mai' is not in list and takes the
    whole parse with it. prodes reads nothing from the header.
    """

    text = "HEADER    HYDROLASE                               02-Mai-99   1ABC\n" + atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N") + "\nEND\n"

    assert len(read_atom_records(text, "german_header")) == 1


def test_a_carriage_return_ends_a_line():
    """a file written on a classic Mac separates its records with \\r alone

    Opening such a file translates the line endings, so prodes read it before
    this change. Handing its text to a reader that splits on \\n only would see
    one enormous line and no records at all.
    """

    text = "\r".join([atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N"), atom_line(2, "CA", "ALA", "A", 1, 1.0, 0.0, 0.0, "C")]) + "\r"

    assert [record.name for record in read_atom_records(text, "carriage_return")] == ["N", "CA"]


def test_windows_line_endings_do_not_reach_the_columns():
    """a trailing carriage return must not end up in the element or anywhere else"""

    text = "\r\n".join([atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N"), atom_line(2, "CA", "ALA", "A", 1, 1.0, 0.0, 0.0, "C")]) + "\r\n"

    assert [record.element for record in read_atom_records(text, "windows")] == ["N", "C"]


def test_a_failure_in_prodes_own_code_is_not_blamed_on_the_file():
    """the reader runs a builder of prodes' own, so its bugs surface as read failures

    Everything the builder does happens inside Biopython's call and is caught by
    the handler that catches an unreadable column. Judging on the innermost
    frame is what keeps a defect in prodes from being reported to the user as a
    defect in their file.
    """

    try:
        readable_occupancy("not a number")
    except TypeError as error:
        assert raised_inside_prodes(error)
    else:
        raise AssertionError("readable_occupancy accepted a string")

    try:
        float("not a number")
    except ValueError as error:
        assert not raised_inside_prodes(error)


def test_the_line_numbers_the_reader_quotes_are_the_files_own(tmp_path, caplog):
    """the header is blanked rather than dropped, so nothing shifts

    A complaint that names line 1 when the record is on line 4 is worse than no
    complaint, because it sends the reader to the wrong place in the file.
    """

    lines = ["HEADER    TEST", "REMARK   1 SOMETHING", "TITLE     X", atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N")[:54]]

    with caplog.at_level(logging.WARNING):
        read_atom_records("\n".join(lines) + "\nEND\n", "offset")

    assert "at line 4" in caplog.text


def test_the_complaints_counter_is_bounded(tmp_path, caplog):
    """every complaint names its own line, so counting distinct ones bounds nothing

    A large file trimmed at column 60 would otherwise build a counter with two
    entries per atom in order to log five of them.
    """

    lines = [atom_line(serial, "N", "ALA", "A", serial, 0.0, 0.0, 0.0, "N")[:54] for serial in range(1, 201)]

    with caplog.at_level(logging.WARNING):
        records = read_atom_records("\n".join(lines) + "\nEND\n", "many")

    assert len(records) == 200
    assert "further complaints" in caplog.text


def test_two_models_written_with_the_same_serial_are_still_two_models(tmp_path):
    """nothing in the format stops a file repeating a MODEL serial

    Numbering by the serial would report such a file as a single model, and the
    warning that exists to stop prodes describing a merged ensemble in silence
    would never fire on it.
    """

    lines = []
    for coordinate in (0.0, 9.0):
        lines += ["MODEL        1", atom_line(1, "N", "ALA", "A", 1, coordinate, 0.0, 0.0, "N"), "ENDMDL"]
    text = open(write_structure(tmp_path / "same_serial.pdb", lines)).read()

    assert [record.model for record in read_atom_records(text, "test")] == [0, 1]


def test_a_warning_that_is_not_the_readers_is_passed_on(tmp_path, recwarn):
    """a DeprecationWarning from inside the reader is somebody else's problem

    Collecting every category would report it to the user as something wrong
    with their file, and would swallow it on the way.
    """

    import warnings as warnings_module

    from prodes.io import pdb_reader

    original = pdb_reader.readable_occupancy

    def warn_and_read(occupancy):
        warnings_module.warn("not the reader's", DeprecationWarning, stacklevel=2)
        return original(occupancy)

    pdb_reader.readable_occupancy = warn_and_read
    try:
        read_atom_records(atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N") + "\nEND\n", "passed_on")
    finally:
        pdb_reader.readable_occupancy = original

    assert [str(warning.message) for warning in recwarn.list] == ["not the reader's"]


def test_the_number_of_further_complaints_is_the_file_s_and_not_the_counter_s(tmp_path, caplog):
    """the cap bounds what is logged, not what is counted

    Every complaint names its own line, so capping distinct messages caps every
    message. Reporting the tail of the cap rather than the tail of the file
    would understate the damage by whatever the file happens to be worth, which
    is a number this package would have invented.
    """

    lines = [atom_line(serial, "N", "ALA", "A", serial, 0.0, 0.0, 0.0, "N")[:54] for serial in range(1, 201)]

    with caplog.at_level(logging.WARNING):
        read_atom_records("\n".join(lines) + "\nEND\n", "many")

    # two complaints per atom, one for the occupancy column and one for the
    # temperature factor, less the five that are logged individually
    assert "and 395 further complaints" in caplog.text


def test_a_segment_id_is_read(tmp_path):
    """columns 73-76, which every shipped structure leaves blank

    Which is why this needs its own file: asserting the field against the column
    on the shipped structures compares an empty string with an empty string
    twenty-eight thousand times and would pass a builder that hard-coded it.
    """

    from tests.pdb_records import pdb_line

    line = pdb_line(
        [
            (1, "ATOM"),
            (7, "    1"),
            (13, " N  "),
            (18, "ALA"),
            (22, "A"),
            (23, "   1"),
            (31, "   0.000"),
            (39, "   0.000"),
            (47, "   0.000"),
            (55, "  1.00"),
            (61, "  0.00"),
            (73, "SEG1"),
            (77, " N"),
        ]
    )
    text = open(write_structure(tmp_path / "segid.pdb", [line])).read()

    assert read_atom_records(text, "test")[0].segment_id == "SEG1"


def test_a_failure_in_prodes_own_code_reaches_the_caller_unwrapped(tmp_path):
    """the guard is wired in, not merely present

    Everything the record builder does happens inside the reader's call, so a
    defect in prodes arrives at the same handler as an unreadable column. Broken
    here on purpose, by replacing the finiteness test readable_occupancy uses
    with something that raises from a frame inside prodes.
    """

    from prodes.io import pdb_reader

    original = pdb_reader.isfinite
    pdb_reader.isfinite = {}.__getitem__
    try:
        with pytest.raises(KeyError):
            read_atom_records(atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, "N") + "\nEND\n", "prodes_bug")
    finally:
        pdb_reader.isfinite = original
