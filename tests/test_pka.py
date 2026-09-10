"""Tests for loading externally predicted pKa values into a structure.

prodes carries one default pKa per titratable residue type, and a
structure-based predictor such as PROPKA, H++ or pypka can override them per
residue. Nothing exercised that path: convert_propka and its two siblings had
no test, Structure.redo_pkas had no test, and the only tests that mention a pKa
file pass a made up filename to a mocked calculate(), which pins the argument
plumbing and nothing else.

It is not a corner of the package that can be left to look after itself. On
1GDW a real PROPKA prediction moves 19 of the 54 reduced features, formal
charge included, so a refactor that quietly stopped applying the values would
change more than a third of the output while every other test still passed.
Downstream code already depends on it: biochai runs PROPKA and imports
convert_propka to read the result, so this is a public API with an outside
caller.

The chain has four links and one test each below::

    PROPKA output -> convert_propka -> JSON -> read_pka -> redo_pkas -> features

tests/data/1GDW.pka is a genuine propka3.4.0 run of 1GDW, committed by the
original author, so the whole chain is covered from the repository alone. None
of this needs PROPKA, or any other predictor, to be installed.
"""

import json
import logging

import pytest

from prodes.core.structure import Structure
from prodes.io.parser import PDBparser, read_pka
from prodes.io.pka_converter import ANY_CHAIN, PROPKA_NOT_TITRATABLE, convert_hpp, convert_propka, convert_pypka, write_json
from tests.pipeline_output import pipeline_output

PDB_PATH = "tests/data/1GDW.pdb.zip"
PROPKA_PATH = "tests/data/1GDW.pka"
FAB_PATH = "tests/data/4NZU.pdb.zip"
FAB_PKA_PATH = "tests/data/4NZU.pka"
HPP_PATH = "tests/data/hpp_1GDW.pkout"
HPP_JSON_PATH = "tests/data/1GDW_hpp_pka.json"


def test_convert_propka_reads_a_real_propka_file():
    """The summary table of a propka3 run becomes a chain to residue number to pKa mapping.

    Values are read out of the file by column position, so the assertions name
    residues at both ends of the chain: a shifted slice would still parse and
    still return plausible looking numbers.
    """

    pkas = convert_propka(PROPKA_PATH)

    assert list(pkas) == ["A"], "1GDW.pka is a single-chain file, update this test if the fixture changes"
    assert len(pkas["A"]) == 45
    assert pkas["A"][18] == [{"ASP": 3.92}]
    assert pkas["A"][35] == [{"GLU": 7.28}]

    # The termini are the two entries propka names differently from the residue
    # types, and the N terminus shares its residue number with a titratable side
    # chain, so both have to land under the same key.
    assert pkas["A"][1] == [{"LYS": 11.34}, {"N+": 8.35}]
    assert pkas["A"][130] == [{"C-": 3.22}]


def test_converted_pkas_survive_the_json_round_trip(tmp_path):
    """What convert_propka writes is what read_pka gives back.

    JSON has no integer keys, so the residue numbers go out as strings and have
    to come back as ints. read_pka is the only thing that converts them, and a
    structure looks up its residues by integer number, so a mapping keyed by
    strings would silently apply to nothing at all.
    """

    pkas = convert_propka(PROPKA_PATH)
    json_path = tmp_path / "1GDW_pka.json"
    write_json(pkas, str(json_path))

    assert read_pka(str(json_path)) == pkas


def test_redo_pkas_overrides_only_the_residues_it_is_given():
    """Predicted values replace the defaults, and untouched residues keep theirs."""

    structure = PDBparser().parse(PDB_PATH)
    residues = {residue.number: residue for residue in structure.residues}

    assert residues[18].pkas == [{"ASP": 3.86}], "default ASP pKa changed, update this test"
    assert residues[1].pkas == [{"LYS": 10.5}, {"N+": 9.69}]

    structure.redo_pkas(convert_propka(PROPKA_PATH))

    assert residues[18].pkas == [{"ASP": 3.92}]
    assert residues[1].pkas == [{"LYS": 11.34}, {"N+": 8.35}]

    # 45 of the 130 residues are in the prediction. The rest are not titratable
    # and must come through untouched rather than being blanked.
    assert residues[2].pkas is None


def test_redo_pkas_is_a_public_method_of_structure():
    """calculate() reaches the override through this name, so it cannot be renamed quietly."""

    assert callable(Structure.redo_pkas)


def test_a_pka_file_changes_the_calculated_features(tmp_path):
    """End to end: predicted pKas reach the output and move the charge features.

    The test is that the values are actually used, not what they come out as.
    Asserting the numbers themselves would pin PROPKA's predictions rather than
    prodes' handling of them, and would have to be rewritten every time the
    charge calculation is legitimately improved.
    """

    json_path = tmp_path / "1GDW_pka.json"
    write_json(convert_propka(PROPKA_PATH), str(json_path))

    default = pipeline_output(PDB_PATH, full_features=False)
    predicted = pipeline_output(PDB_PATH, full_features=False, pkas_file=str(json_path))

    assert list(default.columns) == list(predicted.columns)
    assert default["ID"].iloc[0] == predicted["ID"].iloc[0]

    numeric = default.select_dtypes("number").columns
    changed = [column for column in numeric if default[column].iloc[0] != predicted[column].iloc[0]]

    # Formal charge is the most direct consequence: a predicted pKa either side
    # of pH 7 flips whether a residue counts as charged at all.
    assert "Formal charge" in changed
    assert "Isoelectric point" in changed
    assert len(changed) > 10, f"only {len(changed)} features responded to the pKa file: {changed}"


# The pKa file and the structure can disagree about a residue, and until issue #6
# the file always won by default: redo_pkas replaced a residue's whole list, so a
# group the file did not mention was deleted rather than left at its default.


def test_a_group_the_file_omits_keeps_its_default():
    """A file naming a side chain but not the terminus must not delete the terminus.

    PROPKA always emits both, so this never bit on a PROPKA file, but the
    deletion was silent and the terminal atom then took its charge from the
    side-chain value.
    """

    structure = PDBparser().parse(PDB_PATH)
    first = structure.residues[0]

    structure.redo_pkas({first.chain.name: {1: [{"LYS": 11.34}]}})

    assert first.pkas == [{"LYS": 11.34}, {"N+": 9.69}]
    assert first.group_pka("N+") == 9.69


def test_omitted_groups_are_reported(caplog):
    """Silence used to mean "use the textbook value", which is rarely what was meant."""

    structure = PDBparser().parse(PDB_PATH)
    first = structure.residues[0]

    with caplog.at_level(logging.WARNING):
        structure.redo_pkas({first.chain.name: {1: [{"LYS": 11.34}]}})

    assert "are not in the pKa file" in caplog.text


def test_a_complete_pka_file_is_reported_silently(caplog):
    """PROPKA covers all 46 titratable groups of 1GDW, so nothing should be flagged."""

    structure = PDBparser().parse(PDB_PATH)

    with caplog.at_level(logging.WARNING):
        structure.redo_pkas(convert_propka(PROPKA_PATH))

    assert caplog.text == ""


def test_a_predicted_pka_cannot_titrate_a_bonded_cysteine():
    """The structure outranks the file: the group the file predicts does not exist."""

    structure = PDBparser().parse(PDB_PATH)
    cysteine = next(residue for residue in structure.residues if residue.number == 6)
    assert cysteine.disulfide_partner is not None, "1GDW cysteine 6 is bonded to 128"

    structure.redo_pkas({cysteine.chain.name: {6: [{"CYS": 9.0}]}})

    assert cysteine.pkas is None
    assert round(cysteine.charge(11), 3) == 0


def test_propka_reports_bonded_cysteines_as_not_titratable():
    """PROPKA writes 99.99 for a bridged cysteine rather than leaving it out.

    The issue that prompted this expected an omission, and the safety net was
    designed around that. It is worth pinning what the file really contains, since
    the two mechanisms have to agree rather than fight.
    """

    pkas = convert_propka(PROPKA_PATH)
    bonded = [6, 30, 65, 77, 81, 95, 116, 128]

    assert [pkas["A"][number] for number in bonded] == [[{"CYS": PROPKA_NOT_TITRATABLE}]] * len(bonded)


def test_a_pka_for_a_group_the_residue_does_not_have_is_dropped(caplog):
    """Such an entry used to reach charged_atoms and raise a TypeError there."""

    structure = PDBparser().parse(PDB_PATH)
    second = structure.residues[1]

    with caplog.at_level(logging.WARNING):
        structure.redo_pkas({second.chain.name: {second.number: [{"SER": 10.0}]}})

    assert "which that residue does not have" in caplog.text
    assert second.pkas is None
    assert second.charge(7) == 0


# Issue #12: redo_pkas and the converters used to key on residue number alone,
# so a value predicted for one chain was offered to every chain. Where the
# residue types differed the value had nowhere to go and was dropped; where
# they matched, as below, it was silently applied to the wrong chain.


def test_a_chain_specific_pka_is_not_applied_to_another_chains_same_numbered_residue():
    """The core issue #12 bug: chain H's prediction must not leak onto chain L.

    4NZU is a Fab, chains H and L, both numbered from 1. Residue 91 is a
    tyrosine on both chains, exactly the case that used to be silently
    misapplied: the group-name safety net in redo_pkas has nothing to catch
    two residues of the same type sharing a number on different chains.
    """

    structure = PDBparser().parse(FAB_PATH)
    residues = {(residue.chain.name, residue.number): residue for residue in structure.residues}

    heavy_91 = residues[("H", 91)]
    light_91 = residues[("L", 91)]
    assert heavy_91.name == "TYR"
    assert light_91.name == "TYR", "4NZU L91 changed type, update this test"

    default_light_pka = light_91.group_pka("TYR")

    structure.redo_pkas({"H": {91: [{"TYR": 2.5}]}})

    assert heavy_91.group_pka("TYR") == 2.5
    assert light_91.group_pka("TYR") == default_light_pka


def test_a_wildcard_pka_dict_applies_to_every_chain():
    """ANY_CHAIN is the explicit, documented form of the pre-9.0 behavior.

    A pre-9.0 file, or a convert_hpp conversion, carries no chain information
    at all, so redo_pkas still has to apply it to every chain that has a
    matching residue.
    """

    structure = PDBparser().parse(FAB_PATH)
    residues = {(residue.chain.name, residue.number): residue for residue in structure.residues}
    heavy_91 = residues[("H", 91)]
    light_91 = residues[("L", 91)]

    structure.redo_pkas({ANY_CHAIN: {91: [{"TYR": 2.5}]}})

    assert heavy_91.group_pka("TYR") == 2.5
    assert light_91.group_pka("TYR") == 2.5


def test_a_chain_specific_entry_takes_priority_over_a_wildcard_entry():
    """A residue's own chain is checked before falling back to ANY_CHAIN."""

    structure = PDBparser().parse(FAB_PATH)
    residues = {(residue.chain.name, residue.number): residue for residue in structure.residues}
    heavy_91 = residues[("H", 91)]
    light_91 = residues[("L", 91)]

    structure.redo_pkas({"H": {91: [{"TYR": 2.5}]}, ANY_CHAIN: {91: [{"TYR": 8.0}]}})

    assert heavy_91.group_pka("TYR") == 2.5
    assert light_91.group_pka("TYR") == 8.0


def test_convert_propka_and_redo_pkas_end_to_end_on_a_real_multi_chain_file():
    """The full chain, on a real multi-chain PROPKA run, not a hand-built dict.

    tests/data/4NZU.pka is a genuine propka3.5.1 run of the 4NZU Fab (chains H
    and L). PROPKA predicts a different pKa for chain H's tyrosine 91 and
    chain L's tyrosine 91, so this is also a real-world check that the two
    values land on the right residue rather than only checking that a
    hand-built dict is looked up correctly.
    """

    pkas = convert_propka(FAB_PKA_PATH)

    assert set(pkas) == {"H", "L"}
    assert pkas["H"][91] == [{"TYR": 11.83}]
    assert pkas["L"][91] == [{"TYR": 14.12}], "PROPKA prediction changed, update this test"

    structure = PDBparser().parse(FAB_PATH)
    residues = {(residue.chain.name, residue.number): residue for residue in structure.residues}

    structure.redo_pkas(pkas)

    assert residues[("H", 91)].group_pka("TYR") == 11.83
    assert residues[("L", 91)].group_pka("TYR") == 14.12


def test_redo_pkas_accepts_a_legacy_flat_dict_directly():
    """A pre-9.0 dict passed straight to redo_pkas, not through read_pka, still works.

    Every example before version 9.0 built a flat {number: [...]} dict by
    hand and passed it straight to redo_pkas, exactly as the tests above this
    section used to. Only read_pka gained the pre-9.0 fallback in this
    change; redo_pkas needs the same one, or a caller that never went through
    read_pka would have their dict silently applied to nothing at all.
    """

    structure = PDBparser().parse(FAB_PATH)
    residues = {(residue.chain.name, residue.number): residue for residue in structure.residues}
    heavy_91 = residues[("H", 91)]
    light_91 = residues[("L", 91)]

    structure.redo_pkas({91: [{"TYR": 2.5}]})

    assert heavy_91.group_pka("TYR") == 2.5
    assert light_91.group_pka("TYR") == 2.5


def test_redo_pkas_warns_when_given_a_legacy_flat_dict(caplog):
    """The pre-9.0 fallback in redo_pkas is explicit, like read_pka's own."""

    structure = PDBparser().parse(FAB_PATH)

    with caplog.at_level(logging.WARNING):
        structure.redo_pkas({91: [{"TYR": 2.5}]})

    assert "pre-9.0 flat pKa dict" in caplog.text


def test_convert_hpp_reads_a_real_hpp_file():
    """H++'s own output has no chain column, so its predictions are wrapped under ANY_CHAIN.

    tests/data/1GDW_hpp_pka.json is a genuine H++ conversion committed
    alongside tests/data/hpp_1GDW.pkout but never wired into a test; this pins
    that a fresh conversion still matches it exactly.
    """

    with open(HPP_JSON_PATH) as f:
        expected = {int(number): pka_list for number, pka_list in json.load(f).items()}

    pkas = convert_hpp(HPP_PATH)

    assert list(pkas) == [ANY_CHAIN]
    assert len(pkas[ANY_CHAIN]) == 213
    assert pkas[ANY_CHAIN] == expected


def test_convert_pypka_reads_chains_from_the_header_lines(tmp_path):
    """pypka writes a "Chain: <id>" header before each chain's residues.

    There is no real pypka output committed in the repository, so this
    fixture is hand-written to match the format the parser assumes: a header
    line naming the chain, then whitespace-separated rows of index, residue
    number, identifier and pKa, repeated per chain, ending at a line starting
    "API". This also covers the two bugs the chain-nesting fix required
    fixing: before this change a second chain's header was never recognized
    as a header at all (it fell into the residue-row parser and raised
    IndexError), only the first chain's block was ever read.
    """

    pypka_output = tmp_path / "structure.out"
    pypka_output.write_text(
        "some preamble\n"
        "Chain: A\n"
        "0     1   NTR    7.50\n"
        "1     1   ASP    3.80\n"
        "2     2   SER    Not\n"
        "Chain: B\n"
        "0     5   GLU    4.20\n"
        "1     8   CTR    3.10\n"
        "API\n"
        "unrelated trailer\n"
    )

    pkas = convert_pypka(str(pypka_output))

    assert set(pkas) == {"A", "B"}
    assert pkas["A"][1] == [{"N+": 7.5}, {"ASP": 3.8}]
    assert 2 not in pkas["A"], "SER is skipped, same as THR"
    assert pkas["B"][5] == [{"GLU": 4.2}]
    assert pkas["B"][8] == [{"C-": 3.1}]


def test_convert_pypka_skips_a_malformed_row_and_logs_it(tmp_path, caplog):
    """A truncated data row is skipped and reported, not silently dropped.

    A blank line is expected between blocks and is skipped without comment;
    a non-blank line with too few fields to be a residue row is not expected,
    and is logged so it does not vanish without a trace.
    """

    pypka_output = tmp_path / "structure.out"
    pypka_output.write_text("Chain: A\n\n0     1   truncated\n1     5   ASP    3.80\nAPI\n")

    with caplog.at_level(logging.WARNING):
        pkas = convert_pypka(str(pypka_output))

    assert pkas["A"][5] == [{"ASP": 3.8}]
    assert 1 not in pkas["A"]
    assert "does not look like a pypka residue row" in caplog.text


def test_convert_pypka_skips_a_chain_header_with_no_colon(tmp_path, caplog):
    """A line starting "Chain" with no ':' is not a real header; it is reported and skipped."""

    pypka_output = tmp_path / "structure.out"
    pypka_output.write_text("ChainXfoo bar baz qux\nAPI\n")

    with caplog.at_level(logging.WARNING):
        pkas = convert_pypka(str(pypka_output))

    assert pkas == {}
    assert "looks like a chain header but has no" in caplog.text


def test_read_pka_rejects_a_file_mixing_the_two_formats(tmp_path):
    """A file with both a flat legacy entry and a chain-keyed one is corrupt, not ambiguous."""

    pka_file = tmp_path / "structure.pka"
    pka_file.write_text(json.dumps({"A": {"5": [{"N+": 7.541}]}, "10": [{"ARG": 14.0}]}))

    with pytest.raises(ValueError, match="mixes the pre-9.0"):
        read_pka(str(pka_file))
