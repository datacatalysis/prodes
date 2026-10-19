"""Finds contiguous patches of hydrophobic, positive and negative surface, and reports them over named residue sets.

The per-point values the surface grid phase computes (projected hydrophobicity
and the formal-charge electrostatic potential) are summarised elsewhere over the
whole surface by sign alone. This module adds the spatial step: it averages each
value over a neighbourhood of the point, marks the points whose average passes a
threshold, links marked points that touch on the 1 A surface lattice into
patches, and reports the area of those patches. Every descriptor can also be
reported restricted to a named set of residues, such as the CDRs of an antibody.

The method and the reasons for each setting are in docs/surface_patches.md.

Patches are found on the whole surface and only then restricted. A patch that
crosses the border of a residue set is one patch, and the minimum area applies
to the whole of it; the restricted columns count only its share inside the set.
"""

import re
from collections import deque
from dataclasses import dataclass, field
from itertools import product

import numpy as np

#: Neighbourhood radius in angstrom over which each point's value is averaged.
DEFAULT_PATCH_RADIUS = 5.0

#: Patches smaller than this, in square angstrom of solvent accessible surface,
#: are discarded: about half the solvent accessible area of one fully exposed
#: phenylalanine. Checked against antibody HIC data and kept; see docs/surface_patches.md.
DEFAULT_PATCH_MIN_AREA = 100.0

#: The largest radius accepted. The neighbourhood average reads the whole point
#: set once per lattice offset within the radius, and the offsets grow with its
#: cube: 515 at 5 A, 14,147 at 15 A.
MAX_PATCH_RADIUS = 15.0

#: A point is hydrophobic when its neighbourhood-averaged hydrophobicity is
#: above this. 0 is the sign of the projected hydrophobicity, which only means
#: something on a scale that has negative values.
DEFAULT_HYDROPHOBIC_THRESHOLD = 0.0

#: A point is positive when its neighbourhood-averaged potential is above this
#: many volts, and negative when below minus this. See docs/surface_patches.md
#: for how the value was chosen.
DEFAULT_EP_THRESHOLD = 0.5

RESIDUE_SET_NAME = re.compile(r"^[A-Za-z0-9]+$")


@dataclass(frozen=True)
class PatchSettings:
    """The settings of the patch step, which runs on every structure. calculate uses PatchSettings() when given none.

    Attributes:
        radius: neighbourhood radius in angstrom for the average.
        min_area: smallest patch kept, in square angstrom.
        hydrophobic_threshold: neighbourhood-averaged hydrophobicity above which
            a point is hydrophobic. None means DEFAULT_HYDROPHOBIC_THRESHOLD,
            which is refused for a scale with no negative values.
        ep_threshold: neighbourhood-averaged potential in volts above which a
            point is positive and below the negative of which it is negative.
        residue_sets: name to residue labels (such as H100A), each reported as
            its own set of restricted columns.
    """

    radius: float = DEFAULT_PATCH_RADIUS
    min_area: float = DEFAULT_PATCH_MIN_AREA
    hydrophobic_threshold: float | None = None
    ep_threshold: float = DEFAULT_EP_THRESHOLD
    # Left out of the hash, since a dict cannot be hashed; equal settings still
    # hash equal, which is all a hash has to guarantee.
    residue_sets: dict = field(default_factory=dict, hash=False)

    def __post_init__(self):
        for name, value in (("radius", self.radius), ("min_area", self.min_area), ("ep_threshold", self.ep_threshold)):
            if isinstance(value, bool) or not isinstance(value, int | float) or not np.isfinite(value) or value < 0:
                raise ValueError(f"patch {name} must be a finite number and not negative, got {value!r}")
        if not 1 <= self.radius <= MAX_PATCH_RADIUS:
            raise ValueError(f"patch radius must be between the 1 A lattice spacing and {MAX_PATCH_RADIUS} A, got {self.radius}")
        threshold = self.hydrophobic_threshold
        if threshold is not None and (isinstance(threshold, bool) or not isinstance(threshold, int | float) or not np.isfinite(threshold)):
            raise ValueError(f"patch hydrophobic_threshold must be a finite number or None, got {threshold!r}")
        if not isinstance(self.residue_sets, dict):
            raise ValueError(f"residue_sets must be a dict of name to residue labels, got {type(self.residue_sets).__name__}")
        for name, labels in self.residue_sets.items():
            if not isinstance(name, str) or not RESIDUE_SET_NAME.match(name):
                raise ValueError(f"residue set name {name!r} must be letters and digits only, since it becomes part of a column name")
            if isinstance(labels, str) or not all(isinstance(label, str) for label in labels):
                raise ValueError(f"residue set {name} must be a list of residue labels such as ['H95', 'H100A'], got {labels!r}")
            if len(labels) == 0:
                raise ValueError(f"residue set {name} is empty")

    def resolved_hydrophobic_threshold(self, hydro_scale):
        """Returns the hydrophobic threshold to use with this scale.

        Raises:
            ValueError: when no threshold was given and the scale has no
                negative values, so that 0 would mark every point hydrophobic.
        """

        if self.hydrophobic_threshold is not None:
            return self.hydrophobic_threshold

        from prodes.data import hydrophobic_scale

        if min(hydrophobic_scale(hydro_scale).values()) >= 0:
            raise ValueError(
                f"the hydrophobicity scale {hydro_scale} has no negative values, so the default patch threshold of "
                f"{DEFAULT_HYDROPHOBIC_THRESHOLD} would mark every point hydrophobic; pass one explicitly, with "
                "--patch-hydrophobic-threshold or PatchSettings(hydrophobic_threshold=...)"
            )

        return DEFAULT_HYDROPHOBIC_THRESHOLD

    def record(self, hydro_scale):
        """Returns the settings as written into prodes_run.json, with the threshold resolved."""

        return {
            "radius_A": self.radius,
            "min_area_A2": self.min_area,
            "hydrophobic_threshold": self.resolved_hydrophobic_threshold(hydro_scale),
            "ep_threshold_volts": self.ep_threshold,
            "residue_sets": {name: list(labels) for name, labels in self.residue_sets.items()},
        }


def lattice_indices(coords):
    """Returns the integer lattice position of each point of the 1 A surface grid.

    The property points are cell centres of a grid of 1 A cells, so their
    offsets from the lowest corner are whole numbers up to rounding. Working on
    the integers makes every neighbour and radius test exact.

    Raises:
        ValueError: if the points are not on a 1 A lattice, or two share a position.
    """

    offsets = coords - coords.min(axis=0)
    indices = np.rint(offsets).astype(np.int64)
    if not np.allclose(offsets, indices, atol=1e-6):
        raise ValueError("the surface points are not on a 1 A lattice")
    if len(np.unique(indices, axis=0)) != len(indices):
        raise ValueError("two surface points share a lattice position")

    return indices


def neighbourhood_offsets(radius):
    """Returns the integer offsets whose length is at most radius, the origin included."""

    reach = int(np.floor(radius))
    steps = np.arange(-reach, reach + 1)
    offsets = np.array(list(product(steps, steps, steps)))

    return offsets[(offsets**2).sum(axis=1) <= radius**2]


def local_average(indices, values, areas, radius):
    """Returns each point's area-weighted mean value over the points within radius of it.

    The point itself is included. Built by placing area and area x value on a
    dense array and reading it back at every offset within the radius, which is
    a few hundred gathers of the whole point set rather than a pairwise distance
    matrix.
    """

    offsets = neighbourhood_offsets(radius)
    pad = int(np.floor(radius))
    shape = tuple(indices.max(axis=0) + 1 + 2 * pad)

    area_grid = np.zeros(shape)
    weighted_grid = np.zeros(shape)
    x, y, z = (indices + pad).T
    area_grid[x, y, z] = areas
    weighted_grid[x, y, z] = areas * values

    area_sum = np.zeros(len(indices))
    weighted_sum = np.zeros(len(indices))
    for dx, dy, dz in offsets:
        area_sum += area_grid[x + dx, y + dy, z + dz]
        weighted_sum += weighted_grid[x + dx, y + dy, z + dz]

    return weighted_sum / area_sum


def label_patches(indices, mask):
    """Returns the patch index of each point, or -1 for a point outside the mask.

    Two marked points belong to the same patch when they touch on the lattice,
    which is any of the 26 neighbours, the same connectivity the flood fill that
    builds the surface uses. Patches are numbered in the order their first point
    appears.
    """

    position = {tuple(index): row for row, index in enumerate(indices) if mask[row]}
    steps = [step for step in product((-1, 0, 1), repeat=3) if step != (0, 0, 0)]

    labels = np.full(len(indices), -1, dtype=np.int64)
    patch = 0
    for row in range(len(indices)):
        if not mask[row] or labels[row] >= 0:
            continue

        labels[row] = patch
        queue = deque([row])
        while queue:
            ix, iy, iz = indices[queue.popleft()]
            for dx, dy, dz in steps:
                neighbour = position.get((ix + dx, iy + dy, iz + dz))
                if neighbour is not None and labels[neighbour] < 0:
                    labels[neighbour] = patch
                    queue.append(neighbour)
        patch += 1

    return labels


def kept_patches(labels, areas, min_area):
    """Returns labels with every patch smaller than min_area set to -1."""

    if not np.any(labels >= 0):
        return labels

    patch_areas = np.bincount(labels[labels >= 0], weights=areas[labels >= 0])
    small = patch_areas < min_area
    kept = labels.copy()
    kept[(labels >= 0) & small[np.maximum(labels, 0)]] = -1

    return kept


def patch_statistics(name, labels, areas, in_set):
    """Returns the summed area, largest area and count of the kept patches inside a set.

    Args:
        name: the population prefix, such as SurfMhp.
        labels: kept patch index per point, -1 for none.
        areas: area per point.
        in_set: boolean per point, True for the points counted.
    """

    counted = (labels >= 0) & in_set
    if not np.any(counted):
        return {f"{name}PatchArea": 0.0, f"{name}PatchMaxArea": 0.0, f"N{name}Patch": 0}

    per_patch = np.bincount(labels[counted], weights=areas[counted])
    present = per_patch > 0

    return {
        f"{name}PatchArea": round(float(per_patch.sum()), 3),
        f"{name}PatchMaxArea": round(float(per_patch.max()), 3),
        f"N{name}Patch": int(np.count_nonzero(present)),
    }


def hydrophobic_strength(local_mhp, areas, in_set):
    """Returns the maximum and the area-weighted sum of the positive part of the local hydrophobicity."""

    return {
        "SurfMhpLocalMax": round(float(local_mhp[in_set].max()), 3),
        "SurfMhpLocalPosSum": round(float(np.sum(np.maximum(local_mhp[in_set], 0) * areas[in_set])), 3),
    }


def check_residue_sets(residue_sets, structure_labels):
    """Raises ValueError if a residue set names a residue the structure does not have.

    Cheap, so calculate runs it as soon as the structure is read rather than
    after the surface has been built.
    """

    known = set(structure_labels)
    for name, labels in residue_sets.items():
        missing = sorted(set(labels) - known)
        if missing:
            raise ValueError(f"residue set {name} names residues the structure does not have: {', '.join(missing)}")


def residue_set_masks(residue_sets, point_residues, structure_labels):
    """Returns a boolean point mask per residue set.

    A surface point belongs to the residue whose dots carry the most area in
    its 1 A cell, so a residue with only a sliver of exposed surface can own
    no point at all.

    Raises:
        ValueError: if a set names a residue the structure does not have, or
            none of its residues owns a point of the surface.
    """

    check_residue_sets(residue_sets, structure_labels)

    masks = {}
    for name, labels in residue_sets.items():
        mask = np.isin(point_residues, list(labels))
        if not np.any(mask):
            raise ValueError(
                f"none of the residues of set {name} owns a surface point: each is buried, exposed too little to own one "
                "of the 1 A surface cells, or on a part of the surface that does not touch the rest; see docs/surface_patches.md"
            )
        masks[name] = mask

    return masks


def calculate_patch_features(coords, areas, residues, lipo, ep, settings, hydro_scale, structure_labels):
    """Returns the patch features and the per-point columns for the bundle.

    Args:
        coords: (N, 3) property point coordinates on the 1 A lattice.
        areas: (N,) area per point in square angstrom.
        residues: (N,) residue label per point.
        lipo: (N,) projected hydrophobicity per point.
        ep: (N,) formal-charge electrostatic potential per point, in volts.
        settings: a PatchSettings.
        hydro_scale: the hydrophobicity scale lipo was projected with.
        structure_labels: every residue label in the structure, to validate the sets.

    Returns:
        features: whole-surface columns, then per residue set SurfArea_<NAME>
            and the same columns with _<NAME> appended.
        point_columns: local_mhp, local_ep and mhp_patch_id per point.
    """

    hydrophobic_threshold = settings.resolved_hydrophobic_threshold(hydro_scale)
    masks = residue_set_masks(settings.residue_sets, residues, structure_labels)
    if np.any(areas <= 0):
        raise ValueError("every surface point needs a positive area, or its neighbourhood average is undefined")

    indices = lattice_indices(coords)
    local_mhp = local_average(indices, lipo, areas, settings.radius)
    local_ep = local_average(indices, ep, areas, settings.radius)

    labels = {
        "SurfMhp": kept_patches(label_patches(indices, local_mhp > hydrophobic_threshold), areas, settings.min_area),
        "SurfEpPos": kept_patches(label_patches(indices, local_ep > settings.ep_threshold), areas, settings.min_area),
        "SurfEpNeg": kept_patches(label_patches(indices, local_ep < -settings.ep_threshold), areas, settings.min_area),
    }

    # Column order is the order of src/prodes/data/features_reduced.yaml, which
    # documents each one; tests/test_patches.py checks the two agree.
    def block(in_set):
        return {
            **patch_statistics("SurfMhp", labels["SurfMhp"], areas, in_set),
            **hydrophobic_strength(local_mhp, areas, in_set),
            **patch_statistics("SurfEpPos", labels["SurfEpPos"], areas, in_set),
            **patch_statistics("SurfEpNeg", labels["SurfEpNeg"], areas, in_set),
        }

    features = block(np.ones(len(coords), dtype=bool))
    for name, mask in masks.items():
        features[f"SurfArea_{name}"] = round(float(areas[mask].sum()), 3)
        features.update({f"{code}_{name}": value for code, value in block(mask).items()})

    point_columns = {"local_mhp": np.round(local_mhp, 4), "local_ep_volts": np.round(local_ep, 4), "mhp_patch_id": labels["SurfMhp"]}

    return features, point_columns
