import logging
import zipfile
from math import isfinite
from pathlib import Path

import numpy as np

from prodes.calculations.disulfides import IDENTITY_SYMMETRY_OPERATOR, assign_disulfides
from prodes.core.atom import Atom
from prodes.core.chain import Chain
from prodes.core.residue import Residue
from prodes.core.structure import Structure
from prodes.io.conformers import elect_conformers, report_alternate_conformers
from prodes.io.pdb_reader import COORDINATE_DECIMALS, pdb_lines, read_atom_records

logger = logging.getLogger(__name__)

# The record types Biopython reads coordinates from, spelled as it spells them:
# six columns, not a stripped name. Anything in this module that counts records
# has to agree with the reader about what one is.
COORDINATE_RECORDS = ("ATOM  ", "HETATM")


def pdb_name_in_archive(archive: zipfile.ZipFile, source):
    """Returns the name of the single PDB member of an open zip archive.

    Args:
        archive: an open ZipFile.
        source: the archive path, used only to make the error message useful.

    Raises:
        ValueError: if the archive does not hold exactly one PDB file.
    """

    names = [name for name in archive.namelist() if name.endswith(".pdb")]
    if len(names) != 1:
        raise ValueError(f"expected exactly one .pdb file inside {source}, found {len(names)}: {names}")

    return names[0]


def extract_pdb(archive, directory):
    """Extracts the PDB file held in a zip archive and returns the path to it.

    Args:
        archive: path to a zip holding exactly one PDB file.
        directory: directory to extract into, normally a temporary one.
    """

    with zipfile.ZipFile(archive) as zipped:
        return zipped.extract(pdb_name_in_archive(zipped, archive), directory)


def read_pdb_text(file):
    """Returns the text of a PDB file, reading straight through a zip archive if there is one.

    Lets a caller size or inspect a structure without extracting it to disk.
    """

    path = Path(file)
    if path.suffix != ".zip":
        return path.read_text()

    with zipfile.ZipFile(path) as zipped:
        return zipped.read(pdb_name_in_archive(zipped, path)).decode()


def read_ssbond_line(line):
    """Reads one SSBOND record and returns (chain1, number1, chain2, number2, sym1, sym2).

    The record names the two cysteines of a disulfide bond and the symmetry
    operator each is under. A blank or missing symmetry field is read as the
    identity operator, because hand-written and trimmed PDB files routinely stop
    before column 72 and a bond would otherwise be discarded for being
    incompletely written rather than for being crystallographic.

    Returns None for a line that cannot be read, so that one malformed record
    costs its own bond rather than the whole file.
    """

    try:
        chain1 = line[15].strip()
        number1 = int(line[17:21])
        chain2 = line[29].strip()
        number2 = int(line[31:35])
    except (IndexError, ValueError):
        logger.warning("could not read an SSBOND record, ignoring it: %s", line.rstrip())
        return None

    symmetry1 = line[59:65].strip() or IDENTITY_SYMMETRY_OPERATOR
    symmetry2 = line[66:72].strip() or IDENTITY_SYMMETRY_OPERATOR

    return chain1, number1, chain2, number2, symmetry1, symmetry2


def coordinate_lines(lines):
    """Returns the ATOM and HETATM records among the lines of a PDB file.

    Matched on the six columns rather than on the stripped record name, which is
    how the reader matches them, so that the count below cannot disagree with
    the reader about what a coordinate record is and then raise for the wrong
    reason.
    """

    return [line for line in lines if line[0:6] in COORDINATE_RECORDS]


def check_records_were_all_read(lines, records, source):
    """Raises if the reader returned fewer records than the file holds coordinate lines.

    Not a defensive assertion about a library bug. Biopython stops reading
    coordinates at the first END or CONECT record and returns everything after
    it as an unparsed trailer, so a file that writes coordinates after either
    one loses them silently. Concatenated structures, files that write CONECT
    for a ligand before the rest of its atoms, and anything that writes END per
    model all reach it. prodes reads such a file today, so the choice is between
    losing atoms quietly and saying so, and for a tool whose output is a set of
    numbers there is no choice.

    Note that only the six column form the specification defines is recognised.
    A bare three character END, which prodes' own write_pdb produces and which
    tests/pdb_records writes, does not truncate anything.

    Raises:
        ValueError: if any coordinate record was not read.
    """

    written = coordinate_lines(lines)
    if len(written) == len(records):
        return

    raise ValueError(
        f"{source} holds {len(written)} coordinate records but only {len(records)} could be read; "
        "the usual cause is an ATOM or HETATM record written after an END or CONECT record, "
        "which the reader treats as the end of the coordinates"
    )


def warn_about_unwritable_coordinates(lines, source):
    """Warns once if a coordinate column carries more precision than the format holds.

    The eight column coordinate field is three decimals, and prodes reads the
    columns through a reader that hands them over as float32 and rounds them
    back to those three decimals, which is exact for every value the field can
    hold. A file writing a fourth decimal loses it. That is a fair trade for a
    non-conforming file and a bad one to make in silence, so it is said out loud
    on the pass that counts the records.
    """

    unwritable = 0
    unusable = 0
    for line in coordinate_lines(lines):
        for column in (line[30:38], line[38:46], line[46:54]):
            try:
                written = float(column)
            except ValueError:
                continue

            # nan and inf parse cleanly and then poison every distance computed
            # from them, so they are their own complaint rather than a rounding
            # one: nan != nan would otherwise report them as imprecise.
            if not isfinite(written):
                unusable += 1
                continue

            # Asked as "does rounding change it" rather than as "does it have a
            # fourth decimal", so that a trailing zero and an exponent, which
            # lose nothing, do not warn.
            if round(written, COORDINATE_DECIMALS) != written:
                unwritable += 1

    if unwritable:
        logger.warning(
            "%s writes %d coordinate columns to a precision the PDB format does not hold; they are read to %d decimals",
            source,
            unwritable,
            COORDINATE_DECIMALS,
        )

    if unusable:
        logger.warning(
            "%s writes %d coordinate columns that are not numbers, and every distance measured from those atoms will be meaningless",
            source,
            unusable,
        )


def warn_about_models(records, source):
    """Warns if the file holds more than one model, because prodes merges them.

    prodes has no notion of a model. Every model's atoms are read into one
    structure, and because a residue number that reappears does not start a new
    residue, the later models' atoms pile onto whichever residue was made last:
    a 20 model NMR ensemble comes out as one impossible structure with tens of
    thousands of atoms in a single residue. That is its own defect with its own
    issue, and it is not fixed here, so the least this can do is say that the
    number it is about to report describes something that is not a protein.
    """

    models = {record.model for record in records}
    if len(models) > 1:
        logger.warning(
            "%s holds %d models and prodes has no notion of one; every model's atoms are read into a single structure, "
            "which is not a structure of anything, so take one model out of the file before describing it",
            source,
            len(models),
        )


def ssbond_records(lines):
    """Returns every SSBOND record of a PDB file, as read by read_ssbond_line.

    Read from the text rather than from the reader. Biopython does not parse
    SSBOND at all: its header stops at the first coordinate record and carries
    the title, the resolution and the journal, and nothing about disulfides. A
    lost SSBOND is a charge change rather than a cosmetic one, so this is read
    the way it always was, from anywhere in the file.
    """

    records = []
    for line in lines:
        if line[0:6].strip() == "SSBOND":
            record = read_ssbond_line(line)
            if record is not None:
                records.append(record)

    return records


def as_object_array(items):
    """Returns a list of prodes objects as the object array the entities hold them in."""

    array = np.empty(len(items), dtype=object)
    array[:] = items
    return array


def build_structure(name, records):
    """Builds a Structure, its chains, its residues and its atoms from elected records.

    The grouping rules are the ones prodes has always used, quirks included,
    because this change is about where the fields come from and not about what a
    residue is:

    - A chain is its name. A name seen before switches back to the chain that
      already has it rather than starting a second one.
    - A residue is a (chain, residue number) pair, and the insertion code is not
      in it, so a Kabat numbered H100 and H100A are one residue. That is issue
      #14 and it is deliberately not fixed here.
    - The test is membership in the numbers already seen in that chain, so a
      number that reappears after another has intervened does not start a new
      residue: its atoms join whichever residue was made most recently. That is
      the mechanism behind the damage an NMR ensemble takes, issue #13, and it
      is also deliberately not fixed here.
    - The N terminus is the first residue of a chain. The C terminus is written
      on the chain's last residue at the moment a new chain name appears, and on
      the structure's last residue at the end, which is not the same as "the
      last residue of every chain" for a file whose chains are interleaved.

    Collections are built as lists and converted to object arrays once. Growing
    them an element at a time with np.concatenate, which is what this replaced,
    cost more than the whole of the rest of the parse.
    """

    structure = Structure(name)
    chains = {}
    chain_order = []
    residue_numbers = {}
    chain_atoms = {}
    chain_residues = {}
    residue_atoms = {}
    atoms = []
    residues = []
    current_chain = None
    current_residue = None

    for record in records:
        atom = Atom(
            record.identifier,
            record.name,
            record.residue_name,
            record.chain_name,
            record.residue_number,
            record.x,
            record.y,
            record.z,
            segment_id=record.segment_id,
            element=record.element,
            altloc=record.altloc,
            occupancy=record.occupancy,
        )
        atom.structure = structure
        atoms.append(atom)

        if record.chain_name not in chains:
            if chain_order:
                chain_residues[current_chain.name][-1].terminus = "C"
            current_chain = Chain(record.chain_name, structure)
            chains[record.chain_name] = current_chain
            chain_order.append(record.chain_name)
            residue_numbers[record.chain_name] = set()
            chain_atoms[record.chain_name] = []
            chain_residues[record.chain_name] = []
        elif record.chain_name != current_chain.name:
            current_chain = chains[record.chain_name]

        atom.chain = current_chain
        chain_atoms[current_chain.name].append(atom)

        if record.residue_number not in residue_numbers[current_chain.name]:
            current_residue = Residue(record.residue_name, structure, record.residue_number, current_chain)
            if not chain_residues[current_chain.name]:
                current_residue.terminus = "N"
            residue_numbers[current_chain.name].add(record.residue_number)
            chain_residues[current_chain.name].append(current_residue)
            residue_atoms[current_residue] = []
            residues.append(current_residue)

        atom.residue = current_residue
        residue_atoms[current_residue].append(atom)

    structure.atoms = as_object_array(atoms)
    structure.residues = as_object_array(residues)
    structure.chains = as_object_array([chains[chain_name] for chain_name in chain_order])
    for chain_name in chain_order:
        chains[chain_name].atoms = as_object_array(chain_atoms[chain_name])
        chains[chain_name].residues = as_object_array(chain_residues[chain_name])
    for residue in residues:
        residue.atoms = as_object_array(residue_atoms[residue])

    return structure


def parse_pdb_text(text, name, identifier="ATOM"):
    """Reads the text of a PDB file into a Structure.

    The records are read in full before anything is built, because which
    conformation of a disordered residue to keep cannot be decided until all of
    them have been seen.

    SSBOND records are collected on the same text and handed to the disulfide
    detection at the end, so that any structure which has been parsed from a
    file already knows which of its cysteines are bonded. Anything that builds a
    Structure some other way has to call assign_disulfides itself. They are only
    read while reading protein atoms: under a different identifier the parser is
    not building residues those records could refer to.

    Args:
        text: the text of a PDB file.
        name: the name to give the structure, which becomes the output ID.
        identifier: the record type to build atoms from, ATOM or HETATM.

    Raises:
        ValueError: if the file holds no records of the requested type, or if
            any coordinate record could not be read.
    """

    lines = pdb_lines(text)
    records = read_atom_records(text, name)
    check_records_were_all_read(lines, records, name)
    warn_about_unwritable_coordinates(lines, name)

    # One record type at a time. Electing over both would let a water that
    # happens to share a chain and residue number with a protein residue decide
    # that residue's conformation.
    records = [record for record in records if record.identifier == identifier]
    if not records:
        raise ValueError(f"{name} holds no {identifier} records")

    warn_about_models(records, name)

    kept, report = elect_conformers(records)

    structure = build_structure(name, kept)
    structure.alternate_conformers = report
    report_alternate_conformers(report, name)

    structure.residues[-1].terminus = "C"
    assign_disulfides(structure, ssbond_records(lines) if identifier == "ATOM" else [])

    return structure


class PDBparser:

    def parse(self, file, identifier="ATOM"):
        """Parses pdb files, either plain or held in a zip archive

        A zipped structure is read straight out of the archive, so that callers
        can point at a .pdb.zip without unpacking it themselves.

        The structure is named after the file the caller named, which is also
        what ends up in the ID column of the output. For an archive that is the
        archive itself and not its member, so bar.pdb.zip holding foo.pdb gives
        a structure called bar.

        Raises:
            ValueError: if the file is not a .pdb or a .pdb.zip.
        """

        self.identifier = identifier
        path = Path(file)

        # Checked before anything is read, so that a name the parser does not
        # handle is refused for that reason rather than for not existing.
        if path.suffix == ".zip":
            name = Path(path.stem).stem
        elif path.suffix == ".pdb":
            name = path.stem
        else:
            raise ValueError("File extension not recognized by parser")

        return parse_pdb_text(read_pdb_text(path), name, identifier)


def read_pka(file):
    """Reads a pka json file and will return a dict containing residue numbers and pkas
    pka json files can be generated using the pka_converter in prodes.io"""

    import json

    with open(file) as f:
        pka_dict = {}
        pkas = json.loads(f.read())
        for residue, pka_list in pkas.items():
            residue = int(residue)
            pka_dict[residue] = []
            for pka in pka_list:
                for identifier, pka_value in pka.items():
                    pka_dict[residue].append({identifier: float(pka_value)})

    return pka_dict


def write_pdb(structure, filename, chain="all"):
    """
    takes a Structure object to generate a PDB formatted file using the Atoms information

    arguments:
    structure = Structure object
    filename = PATH + name of the to be generated file
    """

    from prodes import data

    if ".pdb" not in filename:
        print("can only make files with pdb extension")

    else:
        if chain == "all":
            atoms = structure.atoms

        else:
            atoms = np.empty([0])
            for struct_chain in structure.chains:
                if struct_chain.name in chain:
                    atoms = np.concatenate((atoms, struct_chain.atoms))

        # The residue number gets columns 23-26 and nothing more, so a 5 digit
        # number would overflow into the chain field and come back as a
        # different residue when the file is read again. Checked before the file
        # is opened, so that an unwritable structure leaves no half written file.
        too_wide = sorted({atom.residue_number for atom in atoms if len(str(atom.residue_number)) > 4})
        if too_wide:
            raise ValueError(f"residue numbers need more than the 4 columns the PDB format gives them: {too_wide[:5]}")

        with open(filename, "w") as f:
            viable_residues = data.all_residues().keys()
            atom_nmbr = 0
            col4 = ""
            col5 = ""
            col6 = ""
            for atom in atoms:

                if atom_nmbr > 0:
                    if col4.strip() in viable_residues and atom.residue_name in viable_residues and col5.strip() != atom.chain_name:

                        atom_nmbr += 1

                        col2 = str(atom_nmbr)
                        for _ in range(5 - len(col2)):
                            col2 = " " + col2

                        f.write(f"TER   {col2}     {col4}{col5}{col6}\n")
                    elif col4.strip() in viable_residues and atom.residue_name not in viable_residues:

                        atom_nmbr += 1

                        col2 = str(atom_nmbr)
                        for _ in range(5 - len(col2)):
                            col2 = " " + col2

                        f.write(f"TER   {col2}     {col4}{col5}{col6}\n")

                atom_nmbr += 1

                col1 = atom.identifier
                for _ in range(6 - len(col1)):
                    col1 = col1 + " "

                col2 = str(atom_nmbr)
                for _ in range(5 - len(col2)):
                    col2 = " " + col2

                col3 = atom.name
                for i in range(5 - len(col3)):
                    if i < 2:
                        col3 = " " + col3
                    else:
                        col3 = col3 + " "

                col4 = atom.residue_name
                for i in range(5 - len(col4)):
                    if i < 1:
                        col4 = " " + col4
                    else:
                        col4 = col4 + " "

                col5 = atom.chain_name

                # Columns 23-26 hold the residue number, which is what the parser
                # reads back from here; writing the residue name made the output
                # unreadable by PDBparser
                col6 = str(atom.residue_number)
                for _ in range(4 - len(col6)):
                    col6 = " " + col6

                col7 = str(atom.x)
                if len(col7) > 11:
                    col7 = col7[:7]
                for _ in range(11 - len(col7)):
                    col7 = " " + col7

                col8 = str(atom.y)
                if len(col8) > 8:
                    col8 = col8[:7]
                for _ in range(8 - len(col8)):
                    col8 = " " + col8

                col9 = str(atom.z)
                if len(col9) > 8:
                    col9 = col9[:7]
                for _ in range(8 - len(col9)):
                    col9 = " " + col9

                col10 = ""
                for _ in range(6 - len(col10)):
                    col10 = " " + col10

                col11 = ""
                for _ in range(6 - len(col11)):
                    col11 = " " + col11

                col12 = atom.segment_id
                for _ in range(4 - len(atom.segment_id)):
                    col12 = " " + col12

                col13 = atom.element
                for _ in range(2 - len(col13)):
                    col13 = " " + col13

                f.write(f"{col1}{col2}{col3}{col4}{col5}{col6} {col7}{col8}{col9}{col10}{col11}      {col12}{col13}\n")
            f.write("END")


class Builder:

    def build_dummy_atom(self, x, y, z, chain_name="A", ep="", size=""):
        identifier = "ATOM"
        name = "X"
        residue_name = "DUM"
        chain_name = chain_name
        residue_number = 0
        segid = ""
        element = "X"

        # Passed by keyword: Atom no longer takes occupancy and temperature_factor,
        # so the old positional call silently put them in segment_id and element
        atom = Atom(
            identifier,
            name,
            residue_name,
            chain_name,
            residue_number,
            x,
            y,
            z,
            segment_id=segid,
            element=element,
        )

        return atom
