import logging
import tempfile
import zipfile
from collections import defaultdict
from dataclasses import dataclass, field
from math import isfinite
from pathlib import Path
from statistics import median

import numpy as np

from prodes.calculations.disulfides import IDENTITY_SYMMETRY_OPERATOR, assign_disulfides
from prodes.core.atom import Atom
from prodes.core.chain import Chain
from prodes.core.residue import Residue
from prodes.core.structure import Structure

logger = logging.getLogger(__name__)

# Ranks a conformation whose file gives no occupancy at all. Below zero so that
# it loses to every stated occupancy, zero included, and the choice between two
# unstated ones falls to the order they appear in.
UNKNOWN_OCCUPANCY = -1.0

# Occupancies closer together than this are the same occupancy. Deposited values
# carry two decimals, so the gap is enormous compared with anything real and
# small compared with the float noise it exists to absorb.
OCCUPANCY_TOLERANCE = 1e-6


@dataclass
class AlternateConformerReport:
    """What collapsing the alternate conformations of one structure came to.

    Recorded in the run record as well as logged, because a warning is gone as
    soon as the run scrolls past and this is what is left to explain the numbers
    to somebody looking at them a month later.
    """

    residues: int = 0
    dropped: int = 0
    renamed: list = field(default_factory=list)

    def summary(self):
        """Returns the counts as plain types, for the run record."""

        return {
            "residues_with_alternates": self.residues,
            "atoms_dropped": self.dropped,
            "residues_renamed": [f"{name} {key[0]}{key[1].strip()}{key[2].strip()}" for key, name, _ in self.renamed],
        }


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


def residue_key(line):
    """Returns (chain, residue sequence number, insertion code) for one atom record.

    The insertion code is in the key even though the parser otherwise ignores
    column 27 and merges residues that differ only by it. That is not a fix for
    the merging, which is its own defect: it only stops one alternate location
    being elected across two residues that happen to share a number, which would
    be a new error rather than an old one.
    """

    return line[21], line[22:26], line[26]


def read_occupancy(line):
    """Returns the occupancy of one atom record, or None where the file does not give one.

    None rather than 0.0, because zero is a real occupancy and a file that stops
    before column 60 has not said the atom is absent. Keeping the two apart lets
    a structure with no occupancies fall through to the tie break instead of
    having every conformer tie at zero.
    """

    try:
        occupancy = float(line[54:60])
    except (IndexError, ValueError):
        return None

    # nan and inf parse cleanly and then poison every comparison downstream,
    # where the failure surfaces as an integer conversion error naming neither
    # the file nor the column. A value that is not a real occupancy is no more
    # usable than a blank one.
    if not isfinite(occupancy):
        return None

    return occupancy


def conformer_rank(occupancies):
    """Returns the occupancy a conformer is ranked by, given its atoms' occupancies.

    The median, not the mean. Occupancy is usually constant across a conformer's
    atoms, and the median recovers that constant while ignoring a single atom
    refined away from the rest; 1CBN residue 34 carries both 0.80 and 1.00 under
    the same letter, and 3NIR refines every alternate atom separately. The mean
    is dragged by those strays and the sum would prefer whichever conformer has
    more atoms rather than the one that is there more of the time.

    A conformer whose file gives no occupancy at all ranks below every conformer
    that has one, so that the choice falls to the order the letters appear in.
    """

    known = [occupancy for occupancy in occupancies if occupancy is not None]
    if not known:
        return UNKNOWN_OCCUPANCY

    return median(known)


def elect_conformers(lines, identifier):
    """Chooses one alternate location per residue and returns the atom records to keep.

    A residue modelled in more than one conformation writes each of them out in
    full, and every one of those atoms is otherwise read as a real atom, so the
    same side chain is counted twice in the surface and its charge is counted
    twice with it. Exactly one conformation is kept.

    The choice is made per residue rather than per atom. Choosing per atom looks
    equivalent and is not: a residue can be modelled as two *different* amino
    acids, and 1CBN residue 22 is a serine and a proline at one position.
    Choosing the highest occupancy atom for each name independently then
    assembles the serine's OG beside the proline's CG and CD, which is a side
    chain no protein has.

    Atoms with no alternate location belong to every conformation and are always
    kept, which is what makes the ordinary case, an ordered backbone with a
    disordered side chain, come out whole.

    Nothing is taken from a losing conformation to complete the winning one.
    That rule was tried and removed: measured over thousands of disordered
    residues it never once supplied a heavy atom, because a conformation is
    written short precisely when it is the minor one, so the best occupied
    conformation is never the less complete. What it did supply was hydrogens
    carrying the other rotamer's coordinates, several of them landing a fraction
    of an Angstrom from an atom they are not bonded to, and a duplicate of any
    atom the residue also wrote without a letter.

    Args:
        lines: every line of the file, in order.
        identifier: the record name being read, ATOM or HETATM.

    Returns:
        (keep, names, report). keep is the set of line indices to parse. names
        maps a residue key to the residue name of the winning conformation,
        which is not always the name its first atom carries. report is an
        AlternateConformerReport for the log and the run record.
    """

    conformers = defaultdict(lambda: defaultdict(list))
    for index, line in enumerate(lines):
        if line[0:6].strip() != identifier:
            continue
        conformers[residue_key(line)][line[16].strip()].append(index)

    keep = set()
    names = {}
    dropped = 0
    affected = 0
    renamed = []
    for key, by_letter in conformers.items():
        letters = [letter for letter in by_letter if letter]
        if not letters:
            keep.update(by_letter[""])
            continue

        affected += 1
        winner = elect_letter(lines, by_letter, letters)
        winning_name = lines[by_letter[winner][0]][17:20].strip()
        names[key] = winning_name

        # Filtered on the residue name as well as the letter. A letter is
        # normally one conformation of one residue, but nothing in the format
        # requires it, and a letter spanning two residue names would otherwise
        # produce the chimera that choosing per residue exists to prevent.
        kept = [index for index in by_letter[winner] if lines[index][17:20].strip() == winning_name]

        keep.update(by_letter[""])
        keep.update(kept)
        dropped += sum(len(by_letter[letter]) for letter in letters) - len(kept)

        # Reported against the name the residue used to end up with, which was
        # that of its first atom, so the warning marks residues whose reported
        # identity actually changes rather than every residue modelled as more
        # than one amino acid.
        first = min(index for indices in by_letter.values() for index in indices)
        if lines[first][17:20].strip() != winning_name:
            renamed.append((key, winning_name, [lines[first][17:20].strip()]))

    return keep, names, AlternateConformerReport(residues=affected, dropped=dropped, renamed=renamed)


def rank_letters(lines, by_letter, letters):
    """Returns one residue's alternate location letters, best first.

    Ranked on occupancy, and on the order the letters appear in the file where
    occupancy cannot separate them. Ties are not a corner worth glossing over:
    across a survey of real structures roughly a quarter of alternate residues
    are written at exactly equal occupancy, so the tie break decides a large
    share of every choice made here and has to give the same answer on each run.

    Occupancies are compared at a tolerance rather than exactly. A median of
    equal values is not the same float for every atom count, so an exact
    comparison would rank a three atom conformation below a two atom one at the
    same nominal occupancy and the tie break would never run at all.
    """

    ranks = {letter: conformer_rank([read_occupancy(lines[index]) for index in by_letter[letter]]) for letter in letters}

    return sorted(letters, key=lambda letter: (-round(ranks[letter] / OCCUPANCY_TOLERANCE), by_letter[letter][0]))


def elect_letter(lines, by_letter, letters):
    """Returns the winning alternate location letter for one residue."""

    return rank_letters(lines, by_letter, letters)[0]


def report_alternate_conformers(report, name):
    """Warns about the alternate conformations collapsed out of one structure.

    One line for the structure rather than one per residue, which on a real
    structure would be hundreds. A structure with no alternates says nothing.

    A residue whose winning conformation carries a different amino acid is named
    individually, because that changes the sequence and not only the geometry,
    and because it is the case where a pKa file generated from the original
    coordinates no longer agrees with the residue it was meant for.
    """

    if report.residues == 0:
        return

    logger.warning(
        "%s has %d residues modelled in more than one conformation; keeping the best occupied of each and dropping %d atoms",
        name,
        report.residues,
        report.dropped,
    )

    for key, winning_name, alternatives in report.renamed:
        logger.warning(
            "%s %s%s%s is modelled as %s as well; keeping %s, so the sequence and any pKa read for it differ from the file",
            name,
            key[0],
            key[1].strip(),
            key[2].strip(),
            " and ".join(alternatives),
            winning_name,
        )


class PDBparser:

    def parse(self, file, identifier="ATOM"):
        """Parses pdb files, either plain or held in a zip archive

        A zipped structure is extracted to a temporary file and parsed from
        there, so that callers can point at a .pdb.zip without unpacking it
        themselves.

        The structure is named after the file the caller named, which is also
        what ends up in the ID column of the output. For an archive that is the
        archive itself and not its member, so bar.pdb.zip holding foo.pdb gives
        a structure called bar.
        """

        self.identifier = identifier
        path = Path(file)

        if path.suffix == ".zip":
            with tempfile.TemporaryDirectory() as directory:
                return self._parse_pdb(Path(extract_pdb(path, directory)), Path(path.stem).stem)

        return self._parse_pdb(path, path.stem)

    def _parse_pdb(self, path: Path, name: str):
        """Parses one plain PDB file under a name chosen by the caller.

        Args:
            path: path to a .pdb file.
            name: name to give the structure, which becomes the output ID.

        Raises:
            ValueError: if the file is not a .pdb.
        """

        if path.suffix != ".pdb":
            raise ValueError("File extension not recognized by parser")

        with open(path) as pdb:
            return self._read_pdb(pdb, name)

    def _read_pdb(self, pdb, name):
        """Function takes a PDB formatted file and returns a Structure object which contains all atom information of the file

        SSBOND records are collected on the same pass and handed to the
        disulfide detection at the end, so that any structure which has been
        parsed from a file already knows which of its cysteines are bonded.
        Anything that builds a Structure some other way has to call
        assign_disulfides itself.
        """

        current_structure = self._create_structure(name)

        # Read in full before building anything: which conformation of a
        # disordered residue to keep cannot be decided until all of them have
        # been seen, and _add_atom commits an atom the moment it is read. The
        # text is small next to the Atom objects the parse is about to make.
        lines = list(pdb)
        keep, winning_names, report = elect_conformers(lines, self.identifier)
        current_structure.alternate_conformers = report
        report_alternate_conformers(report, name)

        ssbond_records = []
        for index, line in enumerate(lines):
            identifier = line[0:6].strip()
            if identifier == self.identifier:
                if index not in keep:
                    continue
                information = self._read_line(line)
                # The residue is named for the conformation that won, not for
                # whichever atom happened to arrive first. A residue modelled as
                # two different amino acids can write a shared backbone atom
                # under the name of the conformation that lost.
                information["residue_name"] = winning_names.get(residue_key(line), information["residue_name"])
                current_atom = Atom(identifier, **information)
                self._add_atom(current_atom, current_structure)

            # Only while reading protein atoms. Under a different identifier the
            # parser is not building residues these records could refer to.
            elif identifier == "SSBOND" and self.identifier == "ATOM":
                record = read_ssbond_line(line)
                if record is not None:
                    ssbond_records.append(record)

        current_structure.residues[-1].terminus = "C"
        assign_disulfides(current_structure, ssbond_records)

        return current_structure

    def _create_structure(self, name):
        """initializes the parser"""

        self.current_chain = ""
        self.current_residue = ""
        self.residue_numbers = {}
        self.chain_names = []
        return Structure(name)

    def _read_line(self, line):
        """Reads one atom record into the fields an Atom is built from.

        Columns 13-16 are the atom name and column 17 is the alternate location.
        Reading the two together, which this did, gave a disordered cysteine the
        name "SG A", and a name like that matches nothing the rest of the package
        looks for: not the charged atoms of its residue, not the aromatic carbons
        that carry their own radius, not the SG a disulfide is found from. The
        residue kept its atoms and quietly lost its charge.

        Returned as a dict and passed to Atom by keyword. Positional was how
        occupancy and temperature factor once ended up in segment_id and element.
        """

        return {
            "name": line[12:16].strip(),
            "altloc": line[16].strip(),
            "residue_name": line[17:20].strip(),
            "chain_name": line[21].strip(),
            "residue_number": int(line[22:26].strip()),
            "x": float(line[30:38].strip()),
            "y": float(line[38:46].strip()),
            "z": float(line[46:54].strip()),
            "occupancy": read_occupancy(line),
            "segment_id": line[72:76].strip(),
            "element": line[76:78].strip(),
        }

    def _add_atom(self, atom, structure):
        """Adds a atom object to a structure object"""

        addarray = np.array([atom])
        structure.atoms = np.concatenate((structure.atoms, addarray))
        atom.structure = structure

        # Find the chain of the atom
        if atom.chain_name not in self.chain_names:
            if len(self.chain_names) != 0:
                self.current_chain.residues[-1].terminus = "C"
            self.current_chain = self._make_new_chain(structure, atom)

        elif atom.chain_name != self.current_chain.name:
            self.current_chain = structure.chains[self.chain_names.index(atom.chain_name)]

        atom.chain = self.current_chain
        self._add_atom_to_chain(atom, self.current_chain)

        # Find the residue of the atom
        if atom.residue_number not in self.residue_numbers[self.current_chain.name]:
            self.current_residue = self._make_new_residue(structure, atom)

        atom.residue = self.current_residue
        self._add_atom_to_residue(atom, self.current_residue)

    def _make_new_chain(self, structure, atom):
        """Makes a new chain object"""

        self.current_chain = Chain(atom.chain_name, atom.structure)
        self._add_chain_to_structure(self.current_chain, structure)
        self.chain_names.append(atom.chain_name)
        self.residue_numbers[self.current_chain.name] = []
        return self.current_chain

    def _add_atom_to_chain(self, atom, chain):
        """Adds an atom object to a chain object"""

        addarray = np.array([atom])
        chain.atoms = np.concatenate([chain.atoms, addarray])

    def _add_chain_to_structure(self, chain, structure):
        """Adds a chain object to a structure object"""

        addarray = np.array([chain])
        structure.chains = np.concatenate((structure.chains, addarray))

    def _make_new_residue(self, structure, atom):
        """Makes a new residue object"""

        self.current_residue = Residue(atom.residue_name, atom.structure, atom.residue_number, atom.chain)
        if len(self.current_chain.residues) == 0:
            self.current_residue.terminus = "N"
        self._add_residue_to_structure(self.current_residue, structure)
        self.residue_numbers[atom.chain_name].append(atom.residue_number)
        self._add_residue_to_chain(self.current_residue, self.current_chain)
        return self.current_residue

    def _add_residue_to_structure(self, residue, structure):
        """ "adds residue to a structure"""

        addarray = np.array([residue])
        structure.residues = np.concatenate((structure.residues, addarray))

    def _add_atom_to_residue(self, atom, residue):
        """adds atom to a residue"""

        addarray = np.array([atom])
        residue.atoms = np.concatenate([residue.atoms, addarray])

    def _add_residue_to_chain(self, residue, chain):
        "adds residue to a chain"

        addarray = np.array([residue])
        chain.residues = np.concatenate([chain.residues, addarray])


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
