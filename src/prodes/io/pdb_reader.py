"""Reads the coordinate records of a PDB file into flat AtomRecords.

Biopython does the column arithmetic. It is used as a record reader rather than
as an object model: ``PDBParser`` is driven with the ``StructureBuilder`` below,
which collects one record per coordinate line instead of assembling Biopython's
``Structure > Model > Chain > Residue > Atom`` tree.

That is a deliberate choice and not the obvious one. Biopython's entities
represent disorder rather than resolving it, and its defaults for resolving it
are weaker than the rules prodes settled on in version 7.0, so building on the
tree would silently undo them. Worse, the tree cannot express what prodes reads:

- A residue modelled as two amino acids that share an atom written with a blank
  alternate location makes ``StructureBuilder.init_residue`` raise, and in
  permissive mode the whole residue is then dropped with a warning and no error.
- An atom with no alternate location arriving after a disordered residue has
  formed is discarded as "defined twice", which is exactly the atom the 7.0
  rules always keep.
- A disordered atom whose occupancy column is blank crashes the comparison in
  ``DisorderedAtom.disordered_add``, and a blank occupancy is a case prodes
  handles on purpose.
- Coordinates are stored as float32, which perturbs every one of them.

Collecting records avoids all four, because no ``DisorderedAtom`` or
``DisorderedResidue`` is ever built. It also hands prodes the fields its own
identity rules need, one per record: the insertion code, the hetero flag and the
model number, none of which the old hand-written reader carried, and the first
and last of which are what version 8.0 uses to tell one residue from another and
one model from the next.

``MMCIFParser`` drives the same builder interface, so an mmCIF reader is a
second entry point onto the same records rather than a second parser.
"""

import logging
import warnings
from collections import Counter
from dataclasses import dataclass
from io import StringIO
from math import isfinite
from pathlib import Path

from Bio.Data.IUPACData import atom_weights
from Bio.PDB.PDBExceptions import PDBConstructionException, PDBConstructionWarning
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.StructureBuilder import StructureBuilder

logger = logging.getLogger(__name__)

# How many distinct complaints from the reader are worth logging before the rest
# are only counted. A file trimmed at column 60 produces two per atom, each
# naming its own line, so "distinct" is no protection on its own.
REPORTED_WARNINGS = 5

# The complaints that are never worth passing on, because prodes reads or writes
# the records they are about.
#
# Biopython recognises the end of a file only as the six column "END   " of the
# specification, so a file ending in a bare "END" is an unrecognised record to
# it. prodes writes that form itself in write_pdb, and 1GDW_h and every file
# round tripped through the writer carry it.
#
# SSBOND is unrecognised to Biopython everywhere, because it does not parse it
# at all. prodes reads it from the text, so a warning that it was ignored is not
# only noise but false.
BENIGN_WARNINGS = ("Ignoring unrecognized record 'END'", "Ignoring unrecognized record 'SSBOND'")

# How many distinct complaints to hold on to at all. Beyond this the counts of
# the ones already seen still rise, so the total stays right, but a file that
# complains differently about every one of its atoms cannot grow the counter
# without bound.
DISTINCT_WARNINGS = 50

# The records that begin the coordinate section, matched as Biopython matches
# them. Everything before the first of them is header, which prodes does not
# read and which is not handed to the reader at all.
COORDINATE_SECTION_RECORDS = ("ATOM  ", "HETATM", "MODEL ")

# Coordinates are three decimals in the format and Biopython hands them over as
# float32, whose spacing at the top of the eight column field is 0.000977. Half
# of that is 0.000488, under the 0.0005 that would round to a different three
# decimal value, so rounding recovers the number the file wrote, exactly, for
# every value the columns can hold. Checked over the whole grid.
COORDINATE_DECIMALS = 3


@dataclass(frozen=True)
class AtomRecord:
    """One coordinate record, with every field read separately.

    A record, not an atom: every conformation of a disordered residue is here,
    and choosing between them is the election's job rather than the reader's.

    The identity fields are what the old reader did not have. ``residue_number``
    and ``insertion_code`` are separate, so a Kabat numbered residue can be told
    from its neighbour, and ``model`` is the model the record was written under,
    so an NMR ensemble can be told from a single structure. From version 8.0
    both decide what is built: the first two thirds of a residue's key in
    ``build_structure``, and which records reach it at all.
    """

    model: int
    identifier: str
    name: str
    altloc: str
    residue_name: str
    chain_name: str
    residue_number: int
    insertion_code: str
    x: float
    y: float
    z: float
    occupancy: float | None
    segment_id: str
    element: str
    # True where the column was blank and the element above is a guess from the
    # atom name rather than something the file said. Carried on the record and
    # not just logged, because build_structure keeps only the elected records
    # and the run record wants the count for the structure it actually built,
    # not for every record the file held.
    element_inferred: bool = False


def infer_element(name, fullname):
    """Guesses the element from the atom name, for a record whose element column is blank.

    Follows the rule ``Bio.PDB.Atom._assign_element`` applies rather than calling it, since
    calling it would mean constructing a ``Bio.PDB.Atom`` per record, which this builder
    deliberately avoids. ``name`` and ``fullname`` are the same two fields ``PDBParser``
    already passes to ``init_atom``: ``fullname`` is columns 13-16 verbatim, and ``name`` is
    the same field with internal spaces stripped, but only where there are none to lose.

    The first character of ``fullname`` is what tells an inorganic element from an organic
    one: a metal is written from column 13, for example ``"FE  "``, while a carbon,
    nitrogen, oxygen, sulphur, phosphorus or hydrogen leaves it blank, for example
    ``" CA "`` for an alpha carbon rather than a calcium ion. A name starting with a digit,
    the convention for a stereo-specific hydrogen such as ``1HB1``, is read from its second
    character instead.

    Returns ``"X"`` where nothing sensible comes out, matching ``Bio.PDB.Atom``'s own
    fallback for an atom neither rule can place. That includes a name that is a single
    digit and nothing else: real files never write one, but a blank element column tends
    to arrive on files with other blank or malformed columns too, and guessing is not
    worth an ``IndexError`` naming neither the file nor the atom.
    """

    if fullname[0].isalpha() and not fullname[2:].isdigit():
        putative_element = name.strip()
    elif name[0].isdigit():
        putative_element = name[1] if len(name) > 1 else ""
    else:
        putative_element = name[0]

    return putative_element.upper() if putative_element.capitalize() in atom_weights else "X"


def readable_occupancy(occupancy):
    """Returns an occupancy prodes can rank by, or None where the file gives none.

    None rather than 0.0, because zero is a real occupancy and a file that stops
    before column 60 has not said the atom is absent. Keeping the two apart lets
    a structure with no occupancies fall through to the tie break instead of
    having every conformer tie at zero.

    nan and inf parse cleanly and then poison every comparison downstream, where
    the failure surfaces as an integer conversion error naming neither the file
    nor the column. A value that is not a real occupancy is no more usable than
    a blank one.
    """

    if occupancy is None or not isfinite(occupancy):
        return None

    return float(occupancy)


class PdbRecordBuilder(StructureBuilder):
    """Collects one AtomRecord per coordinate record instead of building entities.

    ``PDBParser`` calls this the way it calls Biopython's own builder, so it
    receives every field already extracted from the columns: the residue name,
    hetero flag, sequence number and insertion code from ``init_residue``, the
    chain from ``init_chain``, the model from ``init_model``, and the rest from
    ``init_atom``. Nothing is assembled and nothing is discarded.

    The ``set_anisou``, ``set_siguij`` and ``set_sigatm`` overrides are not
    decoration. Biopython calls them on ``self.atom``, which this builder never
    sets, so without them the first ANISOU record of a real structure raises
    ``AttributeError: 'NoneType' object has no attribute 'set_anisou'``. 4NZU has
    6190 of them.

    ``set_line_counter`` is deliberately not overridden: the inherited one
    assigns an attribute and touches nothing else.
    """

    def __init__(self):
        super().__init__()
        self.records = []
        self.model_number = 0
        self.chain_name = ""
        self.residue = ("", " ", 0, " ")
        self.segment_id = ""

    def init_structure(self, structure_id):
        """Ignores the structure name, which the caller already has."""

    def init_model(self, model_id, serial_num=None):
        """Starts a new model.

        Numbered by the reader's own count of models rather than by the serial
        the MODEL record carries. The serial is what a person would quote, but
        nothing in the format stops a file writing the same one twice, and a
        model number that can collide is not a key: counting distinct serials
        would report a two model file as one. The count cannot collide.
        """

        self.model_number = model_id

    def init_chain(self, chain_id):
        """Starts a new chain."""

        self.chain_name = chain_id.strip()

    def init_seg(self, segid):
        """Starts a new segment."""

        self.segment_id = segid.strip()

    def init_residue(self, resname, field, resseq, icode):
        """Starts a new residue.

        Held rather than built. The parser calls this again whenever the residue
        name changes, so every record is filed under the name its own columns
        carry, which is what the election's residue name filter depends on.
        """

        self.residue = (resname, field, resseq, icode)

    def init_atom(self, name, coord, b_factor, occupancy, altloc, fullname, serial_number=None, element=None, **keywords):
        """Appends one record.

        The name comes from ``fullname``, which is columns 13-16 verbatim, and
        not from ``name``, which is the same string put through a splitting rule
        that belongs to Biopython and could change. Columns 13-16 stripped is
        what prodes has always read, and a name carrying a space matches nothing
        the rest of the package looks for: not the charged atoms of its residue,
        not the aromatic carbons that carry their own radius, not the SG a
        disulfide is found from.

        A blank element column is filled in by ``infer_element`` rather than
        left blank. Left blank, it reaches ``data.vdw_radius("")`` from the
        surface area calculation and raises, and reaches ``Residue.heavy_atoms``
        and ``Residue.protons`` before that, where it counts every hydrogen as
        heavy and no atom as a proton. A column that is not blank is kept
        exactly as written, guesswork being a fallback and not a correction.
        """

        resname, field, resseq, icode = self.residue
        element = (element or "").strip()
        element_inferred = not element
        if element_inferred:
            element = infer_element(name, fullname)
        self.records.append(
            AtomRecord(
                model=self.model_number,
                identifier="ATOM" if field == " " else "HETATM",
                name=fullname.strip(),
                altloc=altloc.strip(),
                residue_name=resname.strip(),
                chain_name=self.chain_name,
                residue_number=resseq,
                insertion_code=icode.strip(),
                x=round(float(coord[0]), COORDINATE_DECIMALS),
                y=round(float(coord[1]), COORDINATE_DECIMALS),
                z=round(float(coord[2]), COORDINATE_DECIMALS),
                occupancy=readable_occupancy(occupancy),
                segment_id=self.segment_id,
                element=element,
                element_inferred=element_inferred,
            )
        )

    def get_structure(self):
        """Returns the records, which is what this builder builds."""

        return self.records

    def set_header(self, header):
        """Ignores the header. prodes reads SSBOND from the text and nothing else."""

    def set_symmetry(self, spacegroup, cell):
        """Ignores the crystal symmetry, which prodes does not use."""

    def set_anisou(self, anisou_array):
        """Ignores anisotropic temperature factors, which prodes does not use."""

    def set_siguij(self, siguij_array):
        """Ignores anisotropic standard deviations, which prodes does not use."""

    def set_sigatm(self, sigatm_array):
        """Ignores coordinate standard deviations, which prodes does not use."""


def report_reader_warnings(complaints, total, source):
    """Logs what the reader complained about, as counted messages.

    Biopython writes its complaints to the warnings machinery, which puts them
    on stderr and out of reach of anything that reads prodes' logs. They are
    collected instead and reported here, so that a file with an unreadable
    column says so once through the same channel as everything else rather than
    twice per atom on the terminal.
    """

    reported = 0
    for message, count in list(complaints.items())[:REPORTED_WARNINGS]:
        logger.warning("%s: %s%s", source, message, f" ({count} times)" if count > 1 else "")
        reported += count

    if total > reported:
        logger.warning("%s: and %d further complaints from the reader", source, total - reported)


def raised_inside_prodes(error):
    """True if the exception came from prodes' own code rather than from the reader.

    The reader is handed a builder of prodes' own, so anything that builder gets
    wrong is raised from inside Biopython's call and caught by the same handler
    that catches an unreadable column. Without this, a defect in prodes is
    reported to the user as a defect in their file.

    Judged on the innermost frame, which is where the exception actually came
    from. Every frame above it is in prodes by construction, because prodes is
    what called the reader.
    """

    traceback = error.__traceback__
    while traceback is not None and traceback.tb_next is not None:
        traceback = traceback.tb_next

    if traceback is None:
        return False

    return traceback.tb_frame.f_code.co_filename.startswith(str(Path(__file__).parent.parent))


def pdb_lines(text):
    """Returns the lines of a PDB file, split the way opening the file would split them.

    A file written on a classic Mac, or by a tool that ends its lines with a
    carriage return, separates records with \\r alone. Reading such a file with
    open() splits it into records, because Python translates line endings in
    text mode, and splitting the same text with a naive split on \\n does not.
    Every reader of the text has to agree about how many lines it has, or the
    count that checks nothing was lost compares one file against another.
    """

    return StringIO(text, newline=None).readlines()


def coordinate_section(text):
    """Returns the lines of a PDB file with the header blanked out.

    The header is blanked rather than parsed. prodes reads nothing from it, and
    Biopython's header parser is a liability on files prodes reads today: it
    looks for a date anywhere in a HEADER or REVDAT record and looks the month
    up in a list of English abbreviations, so a file stamped in any other
    language raises ``ValueError: 'Mai' is not in list`` and takes the whole
    parse with it. Blanking also silences the reader's complaints about header
    records it does not recognise, SSBOND among them, which prodes does read,
    from the text and not from here.

    Blanked and not dropped, so that the line numbers the reader quotes in its
    complaints are still the file's own.
    """

    lines = pdb_lines(text)
    for index, line in enumerate(lines):
        if line[0:6] in COORDINATE_SECTION_RECORDS:
            return [""] * index + lines[index:]

    return [""] * len(lines)


def read_atom_records(text, source):
    """Returns every ATOM and HETATM record of a PDB file, in file order.

    Both record types come back together and carry their own ``identifier``;
    choosing between them is the caller's business.

    Args:
        text: the text of a PDB file.
        source: the name of the file it came from, for error messages.

    Not thread safe: the reader's complaints are captured by replacing the
    warnings module's global handler for the length of the read. prodes parses
    in the parent process before its workers start, so nothing in the package
    reaches this, and run.py already documents the package as process-parallel
    rather than thread-parallel.

    Raises:
        ValueError: if the file cannot be read. Three of Biopython's own column
            reads sit outside its exception handling and would otherwise escape
            as a bare IndexError or ValueError naming neither the file nor the
            column: the insertion code, the residue sequence number, and the
            ANISOU records, which prodes does not use and which cannot break a
            parse today.
    """

    builder = PdbRecordBuilder()
    complaints = Counter()
    total = 0

    def collect(message, category, filename, lineno, file=None, line=None):
        # Only the reader's own complaints. A DeprecationWarning raised
        # somewhere inside it is somebody else's problem and is passed on
        # rather than counted as something wrong with the file.
        nonlocal total
        if not issubclass(category, PDBConstructionWarning):
            passed_on(message, category, filename, lineno, file, line)
            return

        text_of_message = str(message)
        if any(text_of_message.startswith(benign) for benign in BENIGN_WARNINGS):
            return

        total += 1
        # Every complaint names its own line, so the number of distinct messages
        # is no bound at all: a 9000 atom file trimmed at column 60 would build a
        # counter of 18000 strings in order to print five of them. The total is
        # counted separately, so capping the messages does not cost the count.
        if len(complaints) < DISTINCT_WARNINGS or text_of_message in complaints:
            complaints[text_of_message] += 1

    with warnings.catch_warnings():
        passed_on = warnings.showwarning
        warnings.simplefilter("always")
        warnings.showwarning = collect
        try:
            PDBParser(structure_builder=builder, PERMISSIVE=True, QUIET=False).get_structure(
                source, StringIO("\n".join(line.rstrip("\n") for line in coordinate_section(text)))
            )
        except (PDBConstructionException, ValueError, IndexError, KeyError) as error:
            if raised_inside_prodes(error):
                raise
            raise ValueError(f"could not read the coordinate records of {source}: {type(error).__name__}: {error}") from error

    report_reader_warnings(complaints, total, source)

    inferred = sum(1 for record in builder.records if record.element_inferred)
    if inferred:
        logger.warning("%s: inferred the element for %d atom(s) with a blank element column", source, inferred)

    return builder.records
