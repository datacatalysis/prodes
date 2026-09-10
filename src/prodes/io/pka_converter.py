import argparse
import json
import logging

from prodes import data

logger = logging.getLogger(__name__)

# PROPKA writes this in place of a pKa for a group it has decided cannot
# titrate, which in practice means a cysteine it found bonded into a disulfide.
# It is not a prediction, and it is passed through unchanged: it puts the charge
# at zero at every reachable pH, which is the same answer prodes reaches from
# the structure itself.
PROPKA_NOT_TITRATABLE = 99.99

# A pKa file keys its predictions by chain, {chain: {residue number: [...]}},
# except when the predictor's own output does not distinguish chains at all.
# ANY_CHAIN marks such an entry: applied to every chain that has a residue of
# that number and type, exactly as the whole file was applied before version
# 9.0, when there was no chain in the format to key on. See issue #12.
ANY_CHAIN = "*"


def is_legacy_pka_dict(pka_dict):
    """True if pka_dict is the pre-9.0 flat {residue number: [...]} shape.

    Used by both prodes.io.parser.read_pka, for a file loaded from disk, and
    Structure.redo_pkas, for a dict built by hand or passed straight from a
    converter's return value rather than through read_pka. A current-format
    dict is {chain: {residue number: [...]}}, so its values are themselves
    dicts; a pre-9.0 dict's values are the [...] lists directly.
    """

    return any(isinstance(value, list) for value in pka_dict.values())


def parse_arguments():
    """parses the arguments"""

    parser = argparse.ArgumentParser(description="reads pka output files and returns a dictionary or json formatted file")
    parser.add_argument("pka_file", help="Path to the pka file to parse")
    parser.add_argument("source", help="software used to generate the pka file")
    parser.add_argument("-o", "--output", help="Path to the output file ")
    arg = parser.parse_args()

    return arg.pka_file, arg.source, arg.output


def write_json(dictionary, outputfile):
    """writes dictionary as a json"""

    written_json = json.dumps(dictionary, indent=4)

    with open(outputfile, "w") as f:
        f.write(written_json)


def convert_hpp(pka_file):
    """Converts H++ format to a dictionary, keyed {ANY_CHAIN: {residue number: [...]}}.

    H++'s own output has no chain column: a residue is named only by its type
    and number, e.g. ``ASP-45``. So unlike convert_propka and convert_pypka,
    this converter has no chain to read, and every entry it produces is
    wrapped under ANY_CHAIN, applied by Structure.redo_pkas to every chain
    that has a residue of that number and type. See issue #12.
    """

    pkas = {}
    with open(pka_file) as f:
        for i, line in enumerate(f):

            if i == 0:
                pass

            else:
                if line[0:4] == "Site":
                    break

                else:
                    splitted_line = line.split()
                    res_numb = int(splitted_line[0].split("-")[1])
                    identifier = splitted_line[0].split("-")[0]
                    res_pka = splitted_line[2]
                    if ">" in res_pka:
                        res_pka = 14.0

                    elif "<" in res_pka:
                        res_pka = 0

                    res_pka = float(res_pka)

                    if identifier[0:2] == "NT":
                        identifier = "N+"

                    elif identifier[0:2] == "CT":
                        identifier = "C-"

                    else:
                        identifier = identifier[0:3].upper()

                    pkas.setdefault(res_numb, []).append({identifier: res_pka})

    return {ANY_CHAIN: pkas}


def convert_propka(pka_file):
    """Converts PROPKA format to a dictionary, keyed {chain: {residue number: [...]}}.

    PROPKA's summary line carries the chain in a column of its own, between
    the residue number and the pKa value, e.g. ``ASP  18 A     3.92``. Before
    version 9.0 that column was skipped over: every entry was keyed by
    residue number alone, so a value predicted for one chain was offered to
    every chain in the structure. See issue #12.
    """

    known_residues = data.all_residues()
    summary = False
    line_numb = 1
    pkas = {}
    with open(pka_file) as f:
        for line in f:

            if "----" in line:
                summary = False

            if summary:
                if line_numb == 1:
                    pass
                else:

                    stripped = line.strip()
                    identifier = stripped[0:4].strip()
                    if identifier in known_residues or identifier in ["N+", "C-"]:
                        res_numb = int(stripped[4:7].strip())
                        chain_id = stripped[7:12].strip()
                        res_pka = float(stripped[12:18].strip())

                        pkas.setdefault(chain_id, {}).setdefault(res_numb, []).append({identifier: res_pka})
                line_numb += 1

            if "SUMMARY OF THIS PREDICTION" in line:
                summary = True

    return pkas


def convert_pypka(pka_file):
    """Converts pypka output to a dictionary, keyed {chain: {residue number: [...]}}.

    pypka writes a "Chain: <id>" header before each chain's block of residue
    rows. Before version 9.0 that header was read only as a one-shot signal
    that data had started, never for the chain id it names, and never
    re-checked once seen: a second chain's own header line fell through into
    the residue-row parser instead and raised an IndexError. Every header is
    now read for its chain id, and the current chain is updated each time one
    is seen. See issue #12.

    A line that looks like a header or a residue row but is not one — no ':'
    on a "Chain" line, fewer than four fields on a data line — is logged and
    skipped rather than raising. A blank line is skipped silently.
    """

    from prodes.data import residue_data

    current_chain = None
    pkas = {}
    with open(pka_file) as f:
        for line in f:

            if line[:3] == "API":
                break

            if line[:5] == "Chain":
                if ":" not in line:
                    logger.warning("%r looks like a chain header but has no ':'; skipping it", line.rstrip("\n"))
                    continue
                current_chain = line.split(":", 1)[1].strip()
                continue

            if current_chain is None:
                continue

            split_line = line.split()
            if not split_line:
                continue
            if len(split_line) < 4:
                logger.warning("%r does not look like a pypka residue row (fewer than 4 fields); skipping it", line.rstrip("\n"))
                continue

            res_numb = int(split_line[1])
            identifier = split_line[2]
            res_pka = split_line[3]
            if identifier != "SER" and identifier != "THR":
                if res_pka == "Not":
                    potential_charge = residue_data(identifier)["potential_charge"]
                    if potential_charge > 0:
                        res_pka = 14
                    else:
                        res_pka = 0

                res_pka = float(res_pka)

                if identifier == "NTR":
                    identifier = "N+"

                elif identifier == "CTR":
                    identifier = "C-"

                pkas.setdefault(current_chain, {}).setdefault(res_numb, []).append({identifier: res_pka})

    return pkas


def main():
    pka_file, source, output_file = parse_arguments()

    convert = f"convert_{source}"
    pkas = globals()[convert](pka_file)
    print(output_file)
    if output_file:
        write_json(pkas, output_file)
    else:
        print(pkas)


if __name__ == "__main__":
    main()
