"""Tests for the surface patch step, prodes.calculations.patches.

The synthetic tests build points directly on the 1 A lattice, where the right
answer can be counted by hand. The tests on a real structure check what the
synthetic ones cannot: that the step is wired into calculate without changing
any default feature, that the point areas add up to the surface area, and that
a residue set selects the points it names.
"""

import tempfile
from functools import cache
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from prodes import run as run_module
from prodes.calculations import grid_wizard
from prodes.calculations.patches import (
    DEFAULT_HYDROPHOBIC_THRESHOLD,
    PatchSettings,
    calculate_patch_features,
    kept_patches,
    label_patches,
    lattice_indices,
    local_average,
    neighbourhood_offsets,
    patch_statistics,
    residue_set_masks,
)
from prodes.feature_dictionary import PATCH_FEATURE_CODES, REDUCED_FEATURE_CODES, FeatureDictionary
from prodes.output import read_features, read_surface_points, write_bundle
from prodes.run import calculate, construct_surface_grid, prepare_structure

STRUCTURE = "tests/data/ARH96693.pdb.zip"

# The first twenty residues of ARH96693, a single chain named A.
N_TERMINUS = [f"A{number}" for number in range(1, 21)]


def block(origin, size):
    """Returns the lattice points of a cube of side size starting at origin."""

    steps = range(size)
    return [(origin[0] + x, origin[1] + y, origin[2] + z) for x in steps for y in steps for z in steps]


# --- the lattice --------------------------------------------------------------


def test_lattice_indices_recover_whole_number_offsets():
    """Cell centres carry float rounding; the indices must not."""

    coords = np.array([[10.1, 5.3, -2.7], [11.1, 5.3, -2.7], [10.1 + 1e-9, 7.3, -1.7]])

    assert lattice_indices(coords).tolist() == [[0, 0, 0], [1, 0, 0], [0, 2, 1]]


def test_points_off_the_lattice_are_refused():
    with pytest.raises(ValueError, match="1 A lattice"):
        lattice_indices(np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]))


def test_points_sharing_a_position_are_refused():
    with pytest.raises(ValueError, match="share a lattice position"):
        lattice_indices(np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]))


def test_the_neighbourhood_boundary_is_exact():
    """(3, 4, 0) lies exactly 5 A away and is inside; (3, 4, 1) is not."""

    offsets = {tuple(offset) for offset in neighbourhood_offsets(5.0)}

    assert (3, 4, 0) in offsets
    assert (0, 0, 5) in offsets
    assert (3, 4, 1) not in offsets
    assert (0, 0, 0) in offsets


# --- the neighbourhood average ------------------------------------------------


def test_local_average_on_a_line():
    """One hydrophobic point at the end of a line of eleven, radius 2."""

    coords = np.array([[x, 0.0, 0.0] for x in range(11)])
    values = np.zeros(11)
    values[0] = 1.0

    local = local_average(lattice_indices(coords), values, np.ones(11), 2.0)

    # Point 0 sees points 0, 1, 2; point 2 sees 0 to 4; point 3 no longer sees 0.
    assert local[0] == pytest.approx(1 / 3)
    assert local[2] == pytest.approx(1 / 5)
    assert local[3] == 0


def test_local_average_is_weighted_by_area():
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])

    local = local_average(lattice_indices(coords), np.array([1.0, 0.0]), np.array([3.0, 1.0]), 1.0)

    assert local == pytest.approx([0.75, 0.75])


# --- patches ------------------------------------------------------------------


def synthetic_surface():
    """Three separated cubes: 27 points, 8 points and 1 point, every point 1 A^2."""

    points = block((0, 0, 0), 3) + block((10, 0, 0), 2) + [(20, 0, 0)]
    coords = np.array(points, dtype=float)
    return lattice_indices(coords), np.ones(len(points))


def test_separated_regions_are_separate_patches():
    indices, _ = synthetic_surface()

    labels = label_patches(indices, np.ones(len(indices), dtype=bool))

    assert sorted(np.bincount(labels).tolist()) == [1, 8, 27]


def test_diagonal_neighbours_are_linked():
    """The surface is built on the 26-neighbourhood, so a corner contact joins two points."""

    indices = lattice_indices(np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [3.0, 3.0, 3.0]]))

    assert label_patches(indices, np.ones(3, dtype=bool)).tolist() == [0, 0, 1]


def test_points_outside_the_mask_are_in_no_patch():
    indices, _ = synthetic_surface()
    mask = np.zeros(len(indices), dtype=bool)
    mask[:27] = True

    labels = label_patches(indices, mask)

    assert set(labels[27:]) == {-1}
    assert set(labels[:27]) == {0}


def test_patches_under_the_minimum_area_are_discarded():
    indices, areas = synthetic_surface()
    labels = kept_patches(label_patches(indices, np.ones(len(indices), dtype=bool)), areas, 8.0)

    statistics = patch_statistics("SurfMhp", labels, areas, np.ones(len(indices), dtype=bool))

    assert statistics == {"SurfMhpPatchArea": 35.0, "SurfMhpPatchMaxArea": 27.0, "NSurfMhpPatch": 2}


def test_restriction_counts_the_in_set_share_of_whole_patches():
    """A set holding 9 points of the 27-point patch and nothing of the 1-point one."""

    indices, areas = synthetic_surface()
    labels = kept_patches(label_patches(indices, np.ones(len(indices), dtype=bool)), areas, 8.0)
    in_set = np.zeros(len(indices), dtype=bool)
    in_set[:9] = True
    in_set[-1] = True

    statistics = patch_statistics("SurfMhp", labels, areas, in_set)

    # The minimum area was judged on the whole patch, so the 9-point share
    # counts although it is itself above 8, and the 1-point patch, discarded
    # on the whole surface, does not come back.
    assert statistics == {"SurfMhpPatchArea": 9.0, "SurfMhpPatchMaxArea": 9.0, "NSurfMhpPatch": 1}


def test_no_patch_gives_zeros():
    indices, areas = synthetic_surface()
    labels = np.full(len(indices), -1)

    assert patch_statistics("SurfEpPos", labels, areas, np.ones(len(indices), dtype=bool)) == {
        "SurfEpPosPatchArea": 0.0,
        "SurfEpPosPatchMaxArea": 0.0,
        "NSurfEpPosPatch": 0,
    }


# --- settings -----------------------------------------------------------------


def test_the_default_threshold_is_refused_on_a_scale_with_no_negative_values():
    with pytest.raises(ValueError, match="no negative values"):
        PatchSettings().resolved_hydrophobic_threshold("mj")


def test_an_explicit_threshold_is_accepted_on_any_scale():
    assert PatchSettings(hydrophobic_threshold=6.0).resolved_hydrophobic_threshold("mj") == 6.0
    assert PatchSettings().resolved_hydrophobic_threshold("mj_scaled") == DEFAULT_HYDROPHOBIC_THRESHOLD


@pytest.mark.parametrize(
    ("arguments", "message"),
    [
        ({"radius": -1.0}, "radius"),
        ({"radius": 0.5}, "lattice spacing"),
        ({"radius": 50.0}, "between"),
        ({"radius": "5"}, "radius"),
        ({"hydrophobic_threshold": "0.1"}, "hydrophobic_threshold"),
        ({"residue_sets": {"CDR": "H95"}}, "list of residue labels"),
        ({"residue_sets": None}, "dict"),
        ({"min_area": float("nan")}, "min_area"),
        ({"ep_threshold": -0.1}, "ep_threshold"),
        ({"residue_sets": {"CDR-H3": ["H95"]}}, "letters and digits"),
        ({"residue_sets": {"CDR": []}}, "empty"),
    ],
)
def test_bad_settings_are_refused(arguments, message):
    with pytest.raises(ValueError, match=message):
        PatchSettings(**arguments)


def test_a_residue_the_structure_does_not_have_is_refused():
    with pytest.raises(ValueError, match="H999"):
        residue_set_masks({"CDR": ["H95", "H999"]}, np.array(["H95"]), ["H95", "H96"])


def test_a_residue_set_with_no_surface_points_is_refused():
    with pytest.raises(ValueError, match="owns a surface point"):
        residue_set_masks({"Core": ["H96"]}, np.array(["H95"]), ["H95", "H96"])


# --- the residue of a surface cell --------------------------------------------


def dot(residue_label, radius):
    """A Shrake-Rupley dot on an atom of the given radius, in the given residue."""

    atom = SimpleNamespace(
        residue=SimpleNamespace(label=residue_label),
        max_area=lambda probe_r: 4 * np.pi * (radius + probe_r) ** 2,
        max_n_points=lambda probe_r: int(4 * np.pi * (radius + probe_r) ** 2 * 2),
    )
    return SimpleNamespace(atom=atom)


def cell_of(*dots):
    return SimpleNamespace(filtered_content=lambda *_: list(dots))


def test_a_cell_belongs_to_the_residue_with_most_area():
    """Two dots of A5 outweigh one of A4."""

    area, residue = grid_wizard.surface_cell_area_and_residue(cell_of(dot("A4", 1.7), dot("A5", 1.7), dot("A5", 1.7)), 1.4, {})

    assert residue == "A5"
    assert area == pytest.approx(3 * 4 * np.pi * 3.1**2 / int(4 * np.pi * 3.1**2 * 2))


def test_a_tie_goes_to_the_first_label_whatever_the_order():
    first = grid_wizard.surface_cell_area_and_residue(cell_of(dot("A5", 1.7), dot("A4", 1.7)), 1.4, {})
    second = grid_wizard.surface_cell_area_and_residue(cell_of(dot("A4", 1.7), dot("A5", 1.7)), 1.4, {})

    assert first[1] == second[1] == "A4"


# --- on a real structure ------------------------------------------------------


@cache
def _run(full_features, with_set):
    with tempfile.TemporaryDirectory(prefix="prodes_patch_tests_") as directory:
        out_file = Path(directory) / "bundle.zip"
        patches = PatchSettings(residue_sets={"Nterm": N_TERMINUS}) if with_set else None
        calculate(STRUCTURE, str(out_file), full_features=full_features, patches=patches)
        return read_features(out_file), read_surface_points(out_file)


def run(full_features=False, with_set=False):
    """Returns copies of the feature and point tables of one run per setting.

    With with_set, the run is given one residue set, Nterm.
    """

    features, points = _run(full_features, with_set)
    return features.copy(), points.copy()


def test_point_areas_add_up_to_the_surface_area():
    """ARH96693 has no internal cavity, so the flood fill reaches every dot."""

    structure = prepare_structure(STRUCTURE, None)
    points = construct_surface_grid(structure, 1.4)

    assert sum(point.area for point in points) == pytest.approx(structure.surface_area(1.4), rel=1e-9)
    assert all(point.residue for point in points)


def test_a_residue_set_only_appends_columns():
    """Patches are found on the whole surface, so a residue set cannot change a whole-surface value."""

    plain, _ = run()
    restricted, _ = run(with_set=True)

    assert list(plain.columns) == ["ID", *REDUCED_FEATURE_CODES]
    assert list(restricted.columns[: len(plain.columns)]) == list(plain.columns)
    assert restricted[plain.columns].equals(plain)


def test_the_columns_are_the_ones_the_dictionary_documents():
    restricted, _ = run(with_set=True)

    assert list(restricted.columns[1:]) == FeatureDictionary().get_feature_codes(residue_sets=["Nterm"])


def test_the_patch_step_reads_the_formal_potential_under_full_features():
    """The partial-charge phase overwrites the potential on the points; the patch step must run before it."""

    reduced, _ = run(with_set=True)
    full, _ = run(full_features=True, with_set=True)
    codes = [code for code in FeatureDictionary().get_feature_codes(residue_sets=["Nterm"]) if code in PATCH_FEATURE_CODES or code.endswith("_Nterm")]

    assert full[codes].equals(reduced[codes])


def test_the_residue_set_selects_the_points_it_names():
    features, points = run(with_set=True)
    in_set = points["residue"].isin(N_TERMINUS)

    assert features["SurfArea_Nterm"].iloc[0] == pytest.approx(points.loc[in_set, "area_A2"].sum(), abs=1e-3)
    assert set(points.loc[in_set, "residue"]) <= set(N_TERMINUS)
    assert 0 < in_set.sum() < len(points)


def test_the_point_table_agrees_with_the_hydrophobic_patch_area():
    features, points = run(with_set=True)
    in_patch = points["mhp_patch_id"] >= 0

    assert features["SurfMhpPatchArea"].iloc[0] == pytest.approx(points.loc[in_patch, "area_A2"].sum(), abs=1e-3)
    assert features["NSurfMhpPatch"].iloc[0] == points.loc[in_patch, "mhp_patch_id"].nunique()


def test_restricted_values_do_not_exceed_whole_surface_ones():
    features = run(with_set=True)[0].iloc[0]

    for code in ("SurfMhpPatchArea", "SurfEpPosPatchArea", "SurfEpNegPatchArea", "SurfMhpLocalPosSum"):
        assert features[f"{code}_Nterm"] <= features[code]
    assert features["SurfMhpLocalMax_Nterm"] <= features["SurfMhpLocalMax"]


def test_calculate_patch_features_refuses_an_unknown_residue_before_any_work():
    with pytest.raises(ValueError, match="H999"):
        calculate_patch_features(
            np.zeros((1, 3)),
            np.ones(1),
            np.array(["A1"]),
            np.zeros(1),
            np.zeros(1),
            PatchSettings(residue_sets={"CDR": ["H999"]}),
            "mj_scaled",
            ["A1"],
        )


def test_settings_with_residue_sets_can_be_hashed():
    """Frozen settings are hashable even though the residue sets are a dict."""

    assert hash(PatchSettings(residue_sets={"CDR": ["H95"]})) == hash(PatchSettings(residue_sets={"CDR": ["H96"]}))


def test_a_point_with_no_area_is_refused():
    with pytest.raises(ValueError, match="positive area"):
        calculate_patch_features(np.zeros((1, 3)), np.zeros(1), np.array(["A1"]), np.zeros(1), np.zeros(1), PatchSettings(), "mj_scaled", ["A1"])


def test_calculate_refuses_bad_patch_settings_before_building_the_surface(monkeypatch, tmp_path):
    """An unknown residue or a refused threshold fails in the first second, not after the surface phase."""

    def surface_was_built(*args, **kwargs):
        raise AssertionError("the surface was built before the patch settings were checked")

    monkeypatch.setattr(run_module, "construct_surface_grid", surface_was_built)

    with pytest.raises(ValueError, match="A999"):
        calculate(STRUCTURE, str(tmp_path / "out.zip"), patches=PatchSettings(residue_sets={"X": ["A999"]}))
    with pytest.raises(ValueError, match="no negative values"):
        calculate(STRUCTURE, str(tmp_path / "out.zip"), hydro_scale="mj", patches=PatchSettings())


def test_residue_labels_read_back_as_text(tmp_path):
    """A chainless label such as 100 must not come back as the integer 100, nor a blank one as NaN."""

    import pandas as pd

    bundle = write_bundle(
        tmp_path / "out.zip",
        "x",
        pd.DataFrame([{"ID": "x", "Area": 1.0}]),
        np.zeros((2, 3)),
        np.zeros(2),
        np.zeros(2),
        STRUCTURE,
        {"prodes_version": "test"},
        "mj_scaled",
        point_columns={"area_A2": np.ones(2), "residue": np.array(["100", ""])},
    )

    assert read_surface_points(bundle)["residue"].tolist() == ["100", ""]


# --- the feature dictionary ---------------------------------------------------


def test_the_dictionary_resolves_a_residue_set_suffix_on_patch_codes_only():
    fd = FeatureDictionary()

    assert fd.get_entry("SurfMhpPatchArea_CDR") == fd.get_entry("SurfMhpPatchArea")
    assert "SurfArea_CDR" in fd
    assert "SurfArea" not in fd
    assert "Area_CDR" not in fd
    assert "SurfMhpPatchArea_CDR_H3" not in fd
    assert "SurfMhpPatchArea_CDR²" not in fd
    with pytest.raises(KeyError, match="residue set suffix"):
        fd.get_entry("Area_CDR")


def test_the_patch_codes_are_the_last_of_the_default_set():
    fd = FeatureDictionary()

    assert fd.get_feature_codes()[-11:] == fd.get_patch_feature_codes() == list(PATCH_FEATURE_CODES)
    assert set(PATCH_FEATURE_CODES) <= set(fd.get_feature_codes(full=True))
    assert fd.get_feature_codes(residue_sets=["CDR"])[65:67] == ["SurfArea_CDR", "SurfMhpPatchArea_CDR"]
