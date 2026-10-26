# Surface patches

Most Prodes features summarise the surface by the sign of each point's value: every hydrophobic point on the protein is pooled, whether the points form one region or twenty scattered ones. The patch step adds the spatial grouping. It reports how much of the surface is covered by contiguous hydrophobic, positive and negative regions, how large the largest one is, and the same quantities restricted to a named set of residues, such as the CDRs of an antibody.

The patch step runs on every structure. Its 11 columns are the last of the 65 default features, and the residue-set columns, when a run is given residue sets, come after them. It adds about 1 s to a run of 57,000 surface points.

Added to a random-forest model of antibody hydrophobic interaction chromatography (HIC) retention built on sequence features, the patch columns raised the mean per-fold Spearman R on the five published folds of the GDPa1 panel (236 antibodies) from 0.533 to 0.573 with ABodyBuilder3 Fv models and from 0.529 to 0.589 with Boltz-2 models of the whole antibody, improving every fold, and raised it on the held-out GDPa3 panel (79 antibodies) from 0.527 to 0.620. Added to models of anion and cation exchange retention, they changed the Spearman R by less than 0.01.

## Running it

Command line:

```bash
prodes Fv.pdb Fv.zip --residue-set CDR=cdr_residues.txt --residue-set Fv=fv_residues.txt
```

Python, where `patches` is optional and only needed to name residue sets or change a setting:

```python
from prodes.run import calculate
from prodes.calculations.patches import PatchSettings

calculate("Fv.pdb", "Fv.zip", patches=PatchSettings(residue_sets={"CDR": ["H27", "H28", "H100A", "L91"]}))
```

A residue set file holds one residue label per line. A label is the chain, the residue number and the insertion code, written as `H100A` or `L91`; anything after a `#` is ignored. A label the structure does not contain stops the run with an error naming it, before the surface is built. Which residues form the CDRs depends on the numbering scheme of the structure and is decided by the caller, not by Prodes.

## Method

| Step | What happens |
| --- | --- |
| Per-point value | Each point of the 1 A surface grid carries the projected hydrophobicity and the formal-charge electrostatic potential that the default features are calculated from. |
| Per-point area and residue | Each point is one occupied 1 A cell of the Shrake-Rupley solvent accessible surface. Its area is the summed area of the surface dots in the cell, weighted as `Atom.surface_area` weights them, so the point areas add up to `Area` apart from internal cavities (see Limitations). Its residue is the residue whose dots carry the most area in the cell. |
| Local value | Each point's value is replaced by the area-weighted mean over every surface point within the radius, the point itself included. |
| Threshold | A point is hydrophobic when its local hydrophobicity is above the hydrophobic threshold, positive when its local potential is above the electrostatic threshold, and negative when it is below the negative of the electrostatic threshold. |
| Patch | Marked points that touch on the surface lattice (any of the 26 neighbours, the connectivity the surface itself is built with) form one patch. Patches smaller than the minimum area are discarded. |
| Restriction | Patches are found on the whole surface. A restricted column counts only the area of patch points owned by a residue in the set. |

The local average is needed because the per-point projected hydrophobicity is short-ranged. Each atom contributes with a weight of `exp(-d)`, a decay length of 1 A: 0.05 at 3 A and 4.5e-5 at 10 A. A surface point sits about 3 A from the nearest atom centres, so its value mostly reflects the one or two atoms closest to it. Averaging over 5 A gives the value the length scale of a surface region. Waibl et al. (Front. Mol. Biosci. 2022) compared averaging radii for hydrophobicity descriptors of antibodies and found 5 A better than 10 A in almost every case.

Finding patches on the whole surface and then restricting them means a patch that crosses the border of a residue set is one patch, and the minimum area is judged on the whole of it. Restricting first would cut patches at the border and make the minimum area depend on where the border falls.

## Settings

| Setting | CLI | Default | Notes |
| --- | --- | --- | --- |
| `radius` | `--patch-radius` | 5.0 A | Radius of the local average, at most 15 A. The cost of the average grows with the cube of the radius. |
| `min_area` | `--patch-min-area` | 100.0 A^2 | Smallest patch kept; see below. |
| `hydrophobic_threshold` | `--patch-hydrophobic-threshold` | 0.0 | The sign of the local hydrophobicity. Refused for a scale with no negative values, where 0 would mark every point hydrophobic; pass a threshold explicitly there. |
| `ep_threshold` | `--patch-ep-threshold` | 0.5 V | Applied symmetrically to positive and negative patches. |
| `residue_sets` | `--residue-set NAME=FILE` | none | Repeatable. `NAME` is letters and digits only. |

With `--hydro mj`, whose values are all positive, the default hydrophobic threshold of 0 would mark every point hydrophobic, so a run with that scale is refused unless `--patch-hydrophobic-threshold` is given.

The hydrophobicity scale is the run's `--hydro` setting. The thresholds are in the units of the local value, so a threshold chosen for one scale does not carry over to another.

The defaults were checked against antibody HIC retention time on the GDPa1 panel (236 antibodies). Settings tuned on that panel raised the correlation there but not on the held-out GDPa3 panel (79 antibodies), so the defaults were kept. Any re-tuning needs an independent panel to check it on.

The tuning searched radius 3 to 12 A, minimum area 0 to 500 A^2, hydrophobic threshold -0.08 to 0.06 and electrostatic threshold 0.1 to 1.0 V, with patches on ABodyBuilder3 Fv models and residue sets for the CDRs and the Fv. The best setting on GDPa1 (radius 3 A, minimum area 0, hydrophobic threshold -0.06, electrostatic threshold 0.1 V) raised the Spearman R of `SurfMhpPatchArea_CDR` against retention time from 0.34 to 0.46 there, but lowered it from 0.355 to 0.328 on GDPa3. A random-forest HIC model on sequence features plus the Prodes columns, fitted on GDPa1, scored a Spearman R of 0.589 on GDPa3 with the tuned setting against 0.615 with the defaults, lower on all five random-forest seeds. Choosing the setting inside each training fold of GDPa1 instead gave 0.579 against 0.581 on its own five published folds.

How each default was arrived at:

- **Minimum area.** Areas here are solvent accessible surface areas with a 1.4 A probe, which are larger than molecular surface areas of the same region, and by an amount that depends on curvature, so a minimum area quoted by a tool that measures on a molecular surface does not transfer to this one. 100 A^2 is about half the solvent accessible area of one fully exposed phenylalanine residue, 210 to 240 A^2 (Tien et al., PLoS ONE 2013): a patch has to be larger than the exposed half of one large hydrophobic side chain.
- **Electrostatic threshold.** At 0 V, the surface with the same sign as the net charge becomes one or two patches covering most of the protein, and the patch area then says little beyond `NSurfPosEpFormal` or the total area. 0.5 V was chosen from the local potential at the default ionic strength of 0.15 M on the eight structures in `tests/data` (1AO6, 1CBN, 1GDW, 1GPB, 1PIT, 4NZU, ARH96693 and ARH98503): over those structures the median of the 75th percentile of the absolute local potential was 0.39 V and the median of the 90th percentile 0.61 V. At 0.5 V the marked points covered between 13 and 26 per cent of the surface area on six of the eight structures, 2 per cent on 1CBN and 32 per cent on 1PIT. The potential is a screened Coulomb sum at a relative permittivity of 4, so thresholds from Poisson-Boltzmann tools, usually given in kT/e, do not carry over.
- **Hydrophobic threshold.** 0 on `mj_scaled` separates hydrophobic from hydrophilic residues. It is lower than the 0.03 the hydrophobicity viewer uses to colour points, because the viewer shows the per-point value and the patches use the local value, whose range is narrower. The local value on `mj_scaled` is negative over most of a typical surface, so at 0 only the few most hydrophobic regions pass, and the hydrophobic patch area changes steeply with the threshold.

All settings, with the hydrophobic threshold resolved and every residue set's labels, are recorded in `prodes_run.json` under `settings.patches`.

## Columns

The whole-surface columns are, in output order:

| Code | Meaning |
| --- | --- |
| `SurfMhpPatchArea` | summed area of the hydrophobic patches, A^2 |
| `SurfMhpPatchMaxArea` | area of the largest hydrophobic patch, A^2 |
| `NSurfMhpPatch` | number of hydrophobic patches |
| `SurfMhpLocalMax` | highest local hydrophobicity |
| `SurfMhpLocalPosSum` | sum over points of area times the positive part of the local hydrophobicity, A^2 |
| `SurfEpPosPatchArea`, `SurfEpPosPatchMaxArea`, `NSurfEpPosPatch` | the same three for positive patches |
| `SurfEpNegPatchArea`, `SurfEpNegPatchMaxArea`, `NSurfEpNegPatch` | the same three for negative patches |

For each residue set `NAME` the run then writes `SurfArea_NAME`, the surface area owned by the set, followed by every code above with `_NAME` appended. Within a set, `NSurf*Patch_NAME` counts the patches with at least one point in the set, and `*PatchMaxArea_NAME` is the largest share of one patch that lies in the set, not the largest patch lying entirely inside it. A framework patch that only touches the CDRs contributes its small CDR share.

`SurfMhpLocalMax` and `SurfMhpLocalPosSum` need no threshold and no minimum area, so they are not affected by either choice. They are similar in spirit to spatial aggregation propensity (Chennamsetty et al., PNAS 2009), which also averages a hydrophobicity scale over a neighbourhood and sums the positive part, although that method works per atom and weights each residue by its exposed side-chain fraction. `SurfMhpLocalPosSum` is in square angstrom times the units of the hydrophobicity scale, so its value depends on the scale.

Every code is documented in `src/prodes/data/features_reduced.yaml`, where the patch features are marked `per_residue_set`. `FeatureDictionary` looks up a code with or without a residue set suffix, and `get_feature_codes(residue_sets=[...])` lists the columns of a run with residue sets.

The point table, `*_surface_points.csv`, carries `area_A2`, `residue`, `local_mhp`, `local_ep_volts` and `mhp_patch_id`, the hydrophobic patch the point belongs to (-1 for none), for each point, so a patch can be traced back to its residues.

## Limitations

- **Truncated structures.** On an Fv or Fab cut out of a larger molecule, the surface that faces the removed domains is exposed only because of the cut, and can form or enlarge whole-surface patches. Columns restricted to a residue set away from the cut, such as the CDRs, are not affected.
- **Clefts.** The neighbourhood is a Euclidean sphere, so at 5 A it can reach across a narrow cleft to the facing surface.
- **Internal cavities and separate surfaces.** The surface grid is built by a flood fill from one point, so a surface component that does not touch the rest, such as an internal cavity or a second molecule that does not touch the first, is not part of it. This applies to every surface feature, not only the patches, and is why the point areas can add up to less than `Area`.
- **Residues with little exposed surface.** A surface point belongs to the residue whose dots carry the most area in its 1 A cell, so a residue with only a sliver of exposed surface can own no point. A residue set none of whose residues owns a point stops the run with an error rather than reporting zeros.
- **Memory.** The local average places the points on a dense array the size of the structure's bounding box, two arrays of 8 bytes per cubic angstrom. This is a few MB for an Fv and a few tens of MB for a whole antibody, and it is not counted against `--mem-limit`.
- **Electrostatic patches and chromatography.** In hydrophobic interaction chromatography the salt concentration is high and the Debye length is a few angstrom, against the 0.15 M Prodes uses by default. The electrostatic patches describe charge distribution at the run's ionic strength and are not intended as hydrophobic interaction descriptors.
- **One structure.** As with every Prodes feature, the values come from one static structure. Patch areas are known to vary between structure models of the same molecule.

## References

- Chennamsetty et al., *PNAS* 2009, DOI: 10.1073/pnas.0904191106. Spatial aggregation propensity.
- Waibl et al., *Front. Mol. Biosci.* 2022. Hydrophobicity scales and averaging radii for antibody descriptors.
- Hoerschinger et al., *J. Chem. Inf. Model.* 2023, DOI: 10.1021/acs.jcim.3c01490. PEP-Patch, connected-component surface patches.
- Li et al., *mAbs* 2025, 17:2474521. PROPERMAB, patch descriptors on the Fv with CDR-localised versions.
- Tien et al., *PLoS ONE* 2013, DOI: 10.1371/journal.pone.0080635. Maximum solvent accessible surface areas of the amino acids.
