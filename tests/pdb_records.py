"""Builders for PDB records, so a test can write exactly the structure it needs.

Shared rather than local to one test module: the parser reads by fixed column
position, so every test that wants to exercise a column has to write that column
correctly, and one builder that gets it right is better than several that each
might not.
"""


def pdb_line(fields):
    """Returns an 80 column PDB record built from (start column, text) pairs.

    Columns are given as the format's own 1-based numbers, so they can be read
    straight off the specification rather than counted out in an f-string.
    """

    line = [" "] * 80
    for column, text in fields:
        line[column - 1 : column - 1 + len(text)] = text

    return "".join(line).rstrip()


def atom_line(serial, name, residue_name, chain, number, x, y, z, element, altloc="", occupancy=1.00, insertion=""):
    """Returns one ATOM record.

    The atom name occupies columns 13-16 and the alternate location column 17,
    which is the distinction this builder exists to keep straight. A name of four
    characters starts in column 13; a shorter one is written from column 14, the
    convention that leaves room for a two letter element symbol.
    """

    placed = f"{name:<4s}" if len(name) == 4 else f" {name:<3s}"

    return pdb_line(
        [
            (1, "ATOM"),
            (7, f"{serial:5d}"),
            (13, placed),
            (17, altloc),
            (18, f"{residue_name:>3s}"),
            (22, chain),
            (23, f"{number:4d}"),
            (27, insertion),
            (31, f"{x:8.3f}"),
            (39, f"{y:8.3f}"),
            (47, f"{z:8.3f}"),
            (55, f"{occupancy:6.2f}"),
            (61, "  0.00"),
            (77, f"{element:>2s}"),
        ]
    )


def cysteine_lines(serial, chain, number, x, altloc="", occupancy=1.00, insertion=""):
    """Returns the six heavy atoms of one cysteine, placed with its SG at x.

    The other atoms only have to be somewhere sensible: nothing under test reads
    them, but a residue with no backbone is not a residue the parser would ever
    produce.
    """

    atoms = [
        ("N", x - 3.0, 0.0, 0.0, "N"),
        ("CA", x - 2.0, 0.0, 0.0, "C"),
        ("C", x - 2.0, 1.5, 0.0, "C"),
        ("O", x - 2.0, 2.5, 0.0, "O"),
        ("CB", x - 1.0, 0.0, 0.0, "C"),
        ("SG", x, 0.0, 0.0, "S"),
    ]

    return [
        atom_line(serial + offset, name, "CYS", chain, number, at_x, at_y, at_z, element, altloc=altloc, occupancy=occupancy, insertion=insertion)
        for offset, (name, at_x, at_y, at_z, element) in enumerate(atoms)
    ]


def write_structure(path, lines):
    """Writes the given records out as a PDB file and returns the path as a string."""

    path.write_text("\n".join(lines) + "\nEND\n")

    return str(path)
