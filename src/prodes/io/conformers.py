"""Chooses one alternate conformation per residue.

A residue modelled in more than one conformation writes each of them out in
full, and every one of those atoms is otherwise read as a real atom, so the same
side chain is counted twice in the surface and its charge is counted twice with
it. This module decides which conformation to keep. The rules were settled in
version 7.0 and are documented in docs/alternate_conformations.md; they are
prodes' own and are deliberately not the ones a PDB reading library would apply.
"""

import dataclasses
import logging
from collections import defaultdict
from dataclasses import dataclass, field
from statistics import median

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
            "residues_renamed": [f"{name} {residue_label(key)}" for key, name, _ in self.renamed],
        }


def residue_label(key):
    """Returns a residue key as it is written in a warning, for example A22 or H100A.

    The key holds a parsed residue number rather than the raw columns, so the
    number is formatted here rather than stripped. Doing it in one place keeps
    the run record and the log agreeing about what a residue is called.
    """

    chain, number, insertion_code = key

    return f"{chain}{number}{insertion_code}"


def residue_key(record):
    """Returns (chain, residue number, insertion code) for one atom record.

    The insertion code is in the key even though the parser otherwise ignores it
    and merges residues that differ only by it. That is not a fix for the
    merging, which is its own defect: it only stops one alternate location being
    elected across two residues that happen to share a number, which would be a
    new error rather than an old one.
    """

    return record.chain_name, record.residue_number, record.insertion_code


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


def rank_letters(records, by_letter, letters):
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

    ranks = {letter: conformer_rank([records[index].occupancy for index in by_letter[letter]]) for letter in letters}

    return sorted(letters, key=lambda letter: (-round(ranks[letter] / OCCUPANCY_TOLERANCE), by_letter[letter][0]))


def elect_letter(records, by_letter, letters):
    """Returns the winning alternate location letter for one residue."""

    return rank_letters(records, by_letter, letters)[0]


def elect_conformers(records):
    """Chooses one alternate location per residue and returns the records to keep.

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
        records: the AtomRecords of one record type, in file order. One type
            only: electing over ATOM and HETATM together would let a water that
            happens to share a chain and residue number with a protein residue
            decide that residue's conformation.

    Returns:
        (kept, report). kept is the records to build atoms from, in file order,
        each carrying the residue name of the conformation that won, which is
        not always the name its first atom carries. report is an
        AlternateConformerReport for the log and the run record.
    """

    conformers = defaultdict(lambda: defaultdict(list))
    for index, record in enumerate(records):
        conformers[residue_key(record)][record.altloc].append(index)

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
        winner = elect_letter(records, by_letter, letters)
        winning_name = records[by_letter[winner][0]].residue_name
        names[key] = winning_name

        # Filtered on the residue name as well as the letter. A letter is
        # normally one conformation of one residue, but nothing in the format
        # requires it, and a letter spanning two residue names would otherwise
        # produce the chimera that choosing per residue exists to prevent.
        kept = [index for index in by_letter[winner] if records[index].residue_name == winning_name]

        keep.update(by_letter[""])
        keep.update(kept)
        dropped += sum(len(by_letter[letter]) for letter in letters) - len(kept)

        # Reported against the name the residue used to end up with, which was
        # that of its first atom, so the warning marks residues whose reported
        # identity actually changes rather than every residue modelled as more
        # than one amino acid.
        first = min(index for indices in by_letter.values() for index in indices)
        if records[first].residue_name != winning_name:
            renamed.append((key, winning_name, [records[first].residue_name]))

    # Returned in file order. The order of a structure's atoms is the order of
    # the file's records, and grouping them by residue here would quietly change
    # it for every consumer, the PDB writer included.
    kept_records = [rename_record(records[index], names) for index in sorted(keep)]

    return kept_records, AlternateConformerReport(residues=affected, dropped=dropped, renamed=renamed)


def rename_record(record, names):
    """Returns the record under the residue name of the conformation that won.

    A residue modelled as two different amino acids can write a shared backbone
    atom under the name of the conformation that lost, so the name an atom
    carries is not always the name of the residue it ends up in.
    """

    winning_name = names.get(residue_key(record))
    if winning_name is None or winning_name == record.residue_name:
        return record

    return dataclasses.replace(record, residue_name=winning_name)


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
            "%s %s is modelled as %s as well; keeping %s, so the sequence and any pKa read for it differ from the file",
            name,
            residue_label(key),
            " and ".join(alternatives),
            winning_name,
        )
