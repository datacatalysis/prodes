import importlib.resources
import json


def hydrophobic_scale(scale):
    """Extracts the hydrophobicity scale"""

    with importlib.resources.open_text("prodes.data", "hydrophobicity.json") as f:
        data = json.load(f)

    try:
        return data[scale]

    except KeyError as err:
        raise KeyError(f"Scale {scale} is not found in the list of hydrophobicity scales") from err


def oxt_hydrophobicity(scale_dict):
    """Returns the hydrophobicity given to the C-terminal oxygen OXT on this scale.

    OXT is the second oxygen of the C-terminal carboxylate, the same chemical
    group as the side chain of aspartate, so it takes the scale's ASP value.
    The atom-level schemes type it that way: the Ghose-Crippen typing of the
    surface_analyses package gives OXT the acid oxygen type O12 of ASP OD2, and
    the Eisenberg-McLachlan and EEF1 solvation models give it the value of ASP
    OD1 and OD2. Every other heavy atom takes the value of its own residue,
    which for OXT would put a hydrophobic value on the carboxylate of a
    C-terminal PHE. See docs/oxt_hydrophobicity.md.

    Args:
        scale_dict: residue name to hydrophobicity, as returned by hydrophobic_scale.
    """

    return scale_dict["ASP"]


def residue_data(residue):
    """loads the data relating to a specific amino acid residue"""
    aas = all_residues()

    try:
        return aas[residue]
    except KeyError as err:
        raise KeyError(f"{residue} is not found in the list with residues " "Please make sure you provide the three letter code all caps") from err


def all_residues():
    """loads all residue data"""

    return data["residues"]


def vdw_radius(element):
    """returns the van der waals radius of a particular element"""

    try:
        return data["vdw_radius"][element]

    except KeyError as err:
        raise KeyError(f"{element} is not found in the list of elements with") from err


def heavy_atom_mass(atom):
    """returns the mass of a heavy atom"""
    return data["heavy_atom_masses"][atom]


def c_term_pka():
    """returns the standard C terminus pka"""

    return data["C_term_pka"]


def n_term_pka():
    """returns the standard N terminus pka"""

    return data["N_term_pka"]


def _open_general():
    """opens the general data"""
    with importlib.resources.open_text("prodes.data", "general_data.json") as f:
        data = json.load(f)

    return data


data = _open_general()
