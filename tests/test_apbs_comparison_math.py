"""Pins the numerical helpers in scripts/apbs_comparison.py against analytic values.

The comparison script deliberately implements trilinear interpolation and a tie
aware rank correlation itself rather than depending on scipy, which is not in the
prodes environment. Those two functions decide every number the comparison
reports, so they are pinned here against values that can be worked out by hand,
which means the check runs anywhere without the dependency.

They were also cross checked directly against scipy once, agreeing to 4e-15 on
interpolation and exactly on Spearman over heavily tied data. These tests are
what keeps them right.
"""

import sys
import zipfile
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from apbs_comparison import (  # noqa: E402
    average_ranks,
    coulomb_potential,
    pearson,
    read_dx,
    sign_agreement,
    spearman,
    top_decile_auc,
    trilinear_sample,
)


def test_trilinear_is_exact_on_a_linear_field():
    """Trilinear interpolation reproduces a linear function exactly, by construction."""

    origin, spacing, shape = np.array([-3.0, 1.5, 0.25]), np.array([0.4689, 0.5044, 0.4553]), (9, 7, 8)
    coefficients = np.array([1.7, -0.9, 2.3])

    axes = np.meshgrid(*[origin[axis] + spacing[axis] * np.arange(shape[axis]) for axis in range(3)], indexing="ij")
    grid = coefficients[0] * axes[0] + coefficients[1] * axes[1] + coefficients[2] * axes[2] + 4.2

    rng = np.random.default_rng(0)
    low = origin + 1e-6
    high = origin + spacing * (np.array(shape) - 1) - 1e-6
    points = rng.uniform(low, high, size=(500, 3))

    expected = points @ coefficients + 4.2
    np.testing.assert_allclose(trilinear_sample(origin, spacing, grid, points), expected, rtol=1e-12)


def test_trilinear_uses_a_separate_spacing_per_axis():
    """A cubic voxel assumption would be wrong by up to 10 per cent per axis."""

    origin, spacing = np.zeros(3), np.array([1.0, 2.0, 4.0])
    grid = np.zeros((2, 2, 2))
    grid[1, 0, 0], grid[0, 1, 0], grid[0, 0, 1] = 1.0, 1.0, 1.0

    # Half a cell along each axis is a different distance on each axis.
    assert trilinear_sample(origin, spacing, grid, np.array([[0.5, 0.0, 0.0]]))[0] == pytest.approx(0.5)
    assert trilinear_sample(origin, spacing, grid, np.array([[0.0, 1.0, 0.0]]))[0] == pytest.approx(0.5)
    assert trilinear_sample(origin, spacing, grid, np.array([[0.0, 0.0, 2.0]]))[0] == pytest.approx(0.5)


def test_trilinear_returns_nan_outside_the_grid():
    """Points outside must be detectable, not silently extrapolated.

    The comparison raises on any NaN here, because a surface point outside the
    fine grid would otherwise be given a plausible looking wrong potential.
    """

    origin, spacing, grid = np.zeros(3), np.ones(3), np.zeros((4, 4, 4))
    points = np.array([[1.5, 1.5, 1.5], [-0.1, 1.0, 1.0], [1.0, 99.0, 1.0]])

    sampled = trilinear_sample(origin, spacing, grid, points)
    assert not np.isnan(sampled[0])
    assert np.isnan(sampled[1]) and np.isnan(sampled[2])


def test_average_ranks_averages_ties():
    """Tied values share their average rank, which is what makes Spearman correct."""

    # 10 is smallest; the three 20s occupy ranks 2, 3, 4 and so all get 3.
    np.testing.assert_allclose(average_ranks([20, 10, 20, 30, 20]), [3.0, 1.0, 3.0, 5.0, 3.0])
    np.testing.assert_allclose(average_ranks([5, 5]), [1.5, 1.5])


def test_spearman_is_one_for_any_monotone_transform():
    """Spearman is invariant to every monotone transform, which is why it cannot see an offset."""

    values = np.array([1.0, 2.0, 3.0, 4.0, 9.0])
    assert spearman(values, np.exp(values)) == pytest.approx(1.0)
    assert spearman(values, -values) == pytest.approx(-1.0)
    # The property that matters for this experiment: a large additive shift is invisible.
    assert spearman(values, values - 1000) == pytest.approx(1.0)


def test_spearman_matches_a_hand_worked_tied_case():
    """Worked by hand: ranks [1, 2.5, 2.5, 4] against [1, 2, 3, 4]."""

    assert spearman([10, 20, 20, 30], [1, 2, 3, 4]) == pytest.approx(0.9486832980505138)


def test_pearson_matches_a_hand_worked_case():
    """Perfect linear relations are exactly 1 and -1 regardless of scale and offset."""

    assert pearson([1, 2, 3], [2, 4, 6]) == pytest.approx(1.0)
    assert pearson([1, 2, 3], [10, 8, 6]) == pytest.approx(-1.0)
    assert pearson([1, 2, 3, 4], [1, 3, 2, 4]) == pytest.approx(0.8)


def test_top_decile_auc_is_one_when_the_ranking_is_perfect():
    """A score that orders the reference perfectly retrieves its top decile perfectly."""

    reference = np.arange(1000.0)
    assert top_decile_auc(reference, reference) == pytest.approx(1.0)
    assert top_decile_auc(-reference, reference) == pytest.approx(0.0)
    # A constant score cannot separate anything, so it sits at chance.
    assert top_decile_auc(np.zeros(1000), reference) == pytest.approx(0.5)


def test_sign_agreement_ignores_exact_zeros():
    """Prodes rounds to two decimals, so sub-millivolt points become exactly 0.00.

    Those belong to neither sign and must not be counted as agreement or disagreement.
    """

    first = np.array([1.0, -1.0, 0.0, 2.0])
    second = np.array([1.0, 1.0, 5.0, 3.0])
    # Only three usable points, two of which agree.
    assert sign_agreement(first, second) == pytest.approx(2 / 3)


def test_coulomb_potential_is_unchanged_by_chunking(monkeypatch):
    """Chunking bounds memory and must not change a single value.

    The sum is over charges within one point, so chunking over points is exactly
    associative and the result must be bit identical, not merely close.
    """

    import apbs_comparison

    rng = np.random.default_rng(1)
    points = rng.uniform(-20, 20, size=(500, 3))
    charge_coords = rng.uniform(-10, 10, size=(60, 3))
    charges = rng.choice([-1.0, -0.5, 0.5, 1.0], size=60)

    monkeypatch.setattr(apbs_comparison, "CHUNK_BYTES", 1e9)
    whole = apbs_comparison.coulomb_potential(points, charge_coords, charges)
    monkeypatch.setattr(apbs_comparison, "CHUNK_BYTES", 2000.0)
    chunked = apbs_comparison.coulomb_potential(points, charge_coords, charges)

    assert np.array_equal(whole, chunked)


def test_coulomb_potential_signs_are_the_right_way_round():
    """A positive charge raises the potential nearby, a negative one lowers it.

    Cheap, and it catches the one error that would silently invert every figure.
    """

    point = np.array([[5.0, 0.0, 0.0]])
    origin = np.array([[0.0, 0.0, 0.0]])

    assert coulomb_potential(point, origin, np.array([1.0]))[0] > 0
    assert coulomb_potential(point, origin, np.array([-1.0]))[0] < 0


def test_coulomb_potential_reproduces_the_prodes_kernel():
    """Pins the comparison to Prodes itself rather than to a reimplementation of it.

    Property_point.set_ep is the reference: whatever it computes is what the
    Prodes features are built from, so model P must match it exactly.
    """

    from prodes.core.point import Property_point
    from prodes.io.parser import PDBparser

    archive = REPO_ROOT / "tests" / "data" / "ARH96693.pdb.zip"
    with zipfile.ZipFile(archive) as zipped:
        name = next(member for member in zipped.namelist() if member.endswith(".pdb"))
        extracted = zipped.extract(name, REPO_ROOT / "tests" / "data")

    structure = PDBparser().parse(extracted)
    Path(extracted).unlink()

    charged = [atom for atom in structure.atoms if atom.charge(ph=7) != 0]
    charge_coords = np.array([[atom.x, atom.y, atom.z] for atom in charged])
    charges = np.array([atom.charge(ph=7) for atom in charged])

    rng = np.random.default_rng(2)
    probes = rng.uniform(-15, 15, size=(25, 3)) + np.array([structure.x, structure.y, structure.z])

    expected = []
    for x, y, z in probes:
        point = Property_point(x, y, z)
        point.set_ep(structure.atoms, ph=7)
        expected.append(point.ep)

    mine = np.round(coulomb_potential(probes, charge_coords, charges), 2)
    np.testing.assert_allclose(mine, expected, atol=0.01)


def test_read_dx_parses_a_grid_written_in_the_apbs_layout(tmp_path):
    """Values are written three per line, z fastest, so a C order reshape is correct.

    The natural wrong guess, x fastest, produces a plausible looking scrambled
    grid rather than an error, so this pins the axis order explicitly.
    """

    counts = (2, 3, 4)
    values = np.arange(np.prod(counts), dtype=float)
    lines = [" ".join(f"{value:.6e}" for value in values[start : start + 3]) for start in range(0, len(values), 3)]

    dx_file = tmp_path / "test.dx"
    dx_file.write_text(
        "# Data from a test\n#\n# POTENTIAL (kT/e)\n#\n"
        f"object 1 class gridpositions counts {counts[0]} {counts[1]} {counts[2]}\n"
        "origin -1.000000e+00 2.000000e+00 3.000000e+00\n"
        "delta 5.000000e-01 0.000000e+00 0.000000e+00\n"
        "delta 0.000000e+00 2.500000e-01 0.000000e+00\n"
        "delta 0.000000e+00 0.000000e+00 1.250000e-01\n"
        f"object 2 class gridconnections counts {counts[0]} {counts[1]} {counts[2]}\n"
        f"object 3 class array type double rank 0 items {len(values)} data follows\n" + "\n".join(lines) + '\nattribute "dep" string "positions"\n'
    )

    origin, spacing, grid = read_dx(dx_file)

    np.testing.assert_allclose(origin, [-1.0, 2.0, 3.0])
    np.testing.assert_allclose(spacing, [0.5, 0.25, 0.125])
    assert grid.shape == counts
    np.testing.assert_allclose(grid, values.reshape(counts))
    # z fastest: the second value written is the neighbour along z, not along x.
    assert grid[0, 0, 1] == 1.0


def test_read_dx_rejects_a_grid_that_is_not_in_kt_per_e(tmp_path):
    """The unit is parsed and asserted rather than assumed."""

    dx_file = tmp_path / "wrong.dx"
    dx_file.write_text(
        "# Data\n#\n# POTENTIAL (kJ/mol)\n#\n"
        "object 1 class gridpositions counts 2 2 2\n"
        "origin 0 0 0\ndelta 1 0 0\ndelta 0 1 0\ndelta 0 0 1\n"
        "object 2 class gridconnections counts 2 2 2\n"
        "object 3 class array type double rank 0 items 8 data follows\n"
        "0 0 0\n0 0 0\n0 0\n"
    )

    with pytest.raises(ValueError, match="kT/e"):
        read_dx(dx_file)
