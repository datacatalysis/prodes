# The hydrophobicity of the terminal oxygen

What hydrophobicity Prodes gives the C-terminal oxygen `OXT`, why, and how much it moves the features. Behaviour from version 12.0.

## What `OXT` is and what it gets

The C terminus of a protein is a carboxylate, `-C(=O)-O`. The PDB format names the first oxygen `O`, which every residue has, and the second one `OXT`, which only the last residue of a chain has.

Prodes assigns hydrophobicity per residue: `Property_point.set_lipo` gives every heavy atom the value of its own residue on the selected scale. `OXT` is the one exception. It takes the **ASP value of the selected scale**, from `prodes.data.oxt_hydrophobicity`, because the C-terminal carboxylate is the same chemical group as the aspartate side chain.

| scale (`--hydro`) | `OXT` value | scale range | rank of ASP, most hydrophilic first |
| --- | --- | --- | --- |
| `mj_scaled` (default) | -0.80 | -1 (LYS) to 1 (PHE) | 2 of 20 |
| `mj` | 3.56 | 2.95 (LYS) to 9.03 (PHE) | 2 of 20 |
| `kd` | -3.5 | -4.5 (ARG) to 4.5 (ILE) | 3 of 20, tied with ASN, GLN, GLU |
| `cw` | -0.51 | -2.28 (HIS) to 1.81 (ILE) | 7 of 20 |

`cw` is the Cowan-Whittaker scale measured at pH 3.4, where the aspartate side chain, with a pKa near 3.9, is mostly protonated and HIS, LYS and ARG are the most hydrophilic residues. ASP therefore sits in the middle of that scale rather than at its hydrophilic end, and so does `OXT`. The C-terminal carboxyl has a lower pKa, near 3.3, so at pH 3.4 it is more ionised than an aspartate side chain, and on `cw` the ASP value probably understates how hydrophilic the C terminus is. GLU would serve as well as ASP: it is the same group, -0.77 against -0.80 on `mj_scaled` and equal on `kd`.

Versions before 12.0 gave `OXT` 1.0 on every scale. That value did not follow `--hydro`: it is the most hydrophobic value on `mj_scaled`, the value of PHE, near the middle of `kd`, and below every residue on `mj`.

`OXT` carries no charge in Prodes, so no electrostatic feature depends on this value.

## Why the ASP value

The atom-level hydrophobicity and solvation schemes all give `OXT` the value of an aspartate side-chain oxygen:

| scheme | `OXT` | ASP side-chain oxygen | source |
| --- | --- | --- | --- |
| Ghose-Crippen atom types, as used for molecular hydrophobicity potentials | `O12`, acid oxygen | `OD2` is `O12` | [`surface_analyses` `crippen.py`](https://github.com/liedllab/surface_analyses/blob/162a7c822acdedc9f281a7cc68efe88d8e33928f/surface_analyses/crippen.py), [`propermab` `surface.py`](https://github.com/regeneron-mpds/propermab/blob/main/propermab/features/surface.py) |
| Eisenberg-McLachlan atomic solvation parameters | `O-`, -24 cal/mol/Å² | `O-` | Eisenberg and McLachlan 1986, doi 10.1038/319199a0; [`surface_analyses` `eisenberg.py`](https://github.com/liedllab/surface_analyses/blob/162a7c822acdedc9f281a7cc68efe88d8e33928f/surface_analyses/eisenberg.py) |
| EEF1 solvation free energy | -10.0 | -10.0 | [`hipatch.py`](https://github.com/Dignon-Lab/Hydration_Free_Energy_Paper_Figures/blob/main/hipatch/hipatch.py) |

The tools that work from a residue scale, as Prodes does, either leave the backbone out of the hydrophobicity sum altogether, `OXT` with it (Chennamsetty's spatial aggregation propensity and its Rosetta implementation), or give `OXT` the value of its residue like any other atom (`surface_analyses` in residue-scale mode). Leaving `OXT` out alone would treat one backbone atom differently from `N`, `CA`, `C` and `O`. Giving it its residue's value would put a hydrophobic value on the carboxylate of a C-terminal PHE.

The atom-level schemes treat both oxygens of the C-terminal carboxylate alike. Prodes changes only `OXT`: the backbone `O` of the last residue keeps its residue's value, as the backbone `O` of every other residue does. Changing that would move every residue's backbone towards atom-level hydrophobicity, which is a larger question than this one atom.

## Measurement

**Method.** `scripts/measure_oxt_hydrophobicity.py` runs the full feature set, at the default ionic strength and without a pKa file, on every structure in `tests/data`, on each of the four scales, three ways:

| variant | structure | `OXT` value |
| --- | --- | --- |
| old | as deposited | 1.0, as before version 12.0 |
| new | as deposited | the ASP value |
| no `OXT` | every `OXT` in an `ATOM` record removed | none to give |

The old value is reproduced by replacing `prodes.data.oxt_hydrophobicity` for the run, so old and new differ in that one number only. On `mj`, which has no negative values, the patch threshold is the mean of the scale. The results are in [`oxt_hydrophobicity/oxt_features.csv`](oxt_hydrophobicity/oxt_features.csv), and `scripts/summarise_oxt_hydrophobicity.py` prints every table below from it.

Six of the eight structures contain one `OXT` each: 1CBN, 1GDW, 1GPB, 4NZU (on chain L), ARH96693 and ARH98503. 1AO6 and 1PIT contain none and are the controls. 1GDW_h, 1GDW with hydrogens added, gives the same values as 1GDW to the last digit and is not counted separately.

**Only hydrophobicity columns move, and only on structures with an `OXT`.** No column outside the `Mhp` block moves on any scale, and the two controls do not move at all:

| scale | columns that move | of which not hydrophobicity | structures that move |
| --- | --- | --- | --- |
| `mj_scaled` | 21 | 0 | 6 of 8 |
| `kd` | 23 | 0 | 6 of 8 |
| `cw` | 22 | 0 | 6 of 8 |
| `mj` | 11 | 0 | 6 of 8 |

**The change is about one per cent.** From the old value to the new, over the six structures with an `OXT`:

| scale | feature | median change, % | largest change, % |
| --- | --- | --- | --- |
| `mj_scaled` | `SurfMhpMean` | -1.09 | -6.06 |
| `mj_scaled` | `SurfMhpSum` | -0.68 | -3.58 |
| `mj_scaled` | `NSurfPosMhp` | -0.21 | -3.13 |
| `mj_scaled` | `NSurfNegMhp` | +0.08 | +0.85 |
| `mj_scaled` | `SurfMhpPatchArea` | -0.86 | -1.85 |
| `mj_scaled` | `SurfMhpLocalPosSum` | -0.82 | -36.63 |
| `kd` | `SurfMhpMean` | -0.42 | -13.04 |
| `kd` | `SurfMhpSum` | -0.39 | -12.21 |
| `kd` | `NSurfPosMhp` | -0.11 | -0.95 |
| `kd` | `NSurfNegMhp` | +0.03 | +0.28 |
| `kd` | `SurfMhpPatchArea` | -0.29 | -16.39 |
| `kd` | `SurfMhpLocalPosSum` | -0.28 | -12.83 |
| `cw` | `SurfMhpMean` | -0.63 | -8.00 |
| `cw` | `SurfMhpSum` | -0.59 | -6.04 |
| `cw` | `NSurfPosMhp` | -0.24 | -1.03 |
| `cw` | `NSurfNegMhp` | +0.19 | +0.70 |
| `cw` | `SurfMhpPatchArea` | -0.32 | -1.42 |
| `cw` | `SurfMhpLocalPosSum` | -0.42 | -4.22 |
| `mj` | `SurfMhpMean` | +0.10 | +0.31 |
| `mj` | `SurfMhpSum` | +0.08 | +0.27 |
| `mj` | `NSurfPosMhp` | +0.00 | +0.00 |
| `mj` | `SurfMhpLocalPosSum` | +0.08 | +0.26 |

The largest relative changes are on small values: `SurfMhpLocalPosSum` of 1GDW on `mj_scaled` falls from 2.886 to 1.829. On the centred scales every feature shown moves towards hydrophilic, and on `mj`, where 1.0 was below every residue, towards hydrophobic.

On the default `mj_scaled`, per structure, each cell giving the value without `OXT` / with the old value / with the new value:

| structure | `SurfMhpMean` | `SurfMhpSum` | `NSurfPosMhp` | `NSurfNegMhp` | `SurfMhpPatchArea` | `SurfMhpLocalPosSum` |
| --- | --- | --- | --- | --- | --- | --- |
| 1CBN | -0.034 / -0.033 / -0.035 | -108.4 / -106.5 / -110.1 | 953 / 954 / 951 | 2232 / 2233 / 2236 | 807.4 / 811.9 / 804.9 | 30.35 / 30.47 / 30.27 |
| 1GDW | -0.087 / -0.086 / -0.087 | -610.2 / -606.3 / -612.6 | 757 / 766 / 742 | 6283 / 6273 / 6297 | 0 / 0 / 0 | 2.046 / 2.886 / 1.829 |
| 1GPB | -0.086 / -0.086 / -0.086 | -2827 / -2824 / -2831 | 5194 / 5195 / 5190 | 2.775e+04 / 2.775e+04 / 2.775e+04 | 2060 / 2072 / 2054 | 116.4 / 119.4 / 113.3 |
| 4NZU | -0.098 / -0.098 / -0.099 | -1947 / -1944 / -1949 | 1769 / 1769 / 1769 | 1.8e+04 / 1.8e+04 / 1.8e+04 | 599.4 / 599.4 / 599.4 | 19.2 / 19.2 / 19.2 |
| ARH96693 | -0.052 / -0.051 / -0.053 | -251.2 / -245.7 / -254.5 | 1465 / 1490 / 1462 | 3326 / 3304 / 3332 | 1068 / 1082 / 1062 | 54.61 / 54.97 / 54.43 |
| ARH98503 | -0.077 / -0.077 / -0.077 | -1596 / -1593 / -1598 | 3562 / 3565 / 3562 | 1.718e+04 / 1.717e+04 / 1.717e+04 | 1669 / 1678 / 1668 | 78.13 / 78.27 / 78.04 |

**An `OXT` still moves the features, now the other way.** Compared with the same structure without its `OXT`, the old value raised the hydrophobicity features and the new value lowers them. How many of the six structures each feature rises or falls on, and on how many the shift from the structure without `OXT` is strictly smaller under the new value:

| scale | feature | OXT raises it, old | OXT lowers it, old | OXT raises it, new | OXT lowers it, new | shift smaller with new |
| --- | --- | --- | --- | --- | --- | --- |
| `mj_scaled` | `SurfMhpMean` | 3 | 0 | 0 | 3 | 1 of 6 |
| `mj_scaled` | `SurfMhpSum` | 6 | 0 | 0 | 6 | 5 of 6 |
| `mj_scaled` | `NSurfPosMhp` | 5 | 0 | 0 | 4 | 2 of 6 |
| `mj_scaled` | `NSurfNegMhp` | 1 | 4 | 4 | 2 | 2 of 6 |
| `mj_scaled` | `SurfMhpPatchArea` | 4 | 0 | 0 | 4 | 4 of 6 |
| `mj_scaled` | `SurfMhpLocalPosSum` | 5 | 0 | 0 | 5 | 4 of 6 |
| `kd` | `SurfMhpMean` | 2 | 0 | 0 | 5 | 0 of 6 |
| `kd` | `SurfMhpSum` | 6 | 0 | 0 | 6 | 0 of 6 |
| `kd` | `NSurfPosMhp` | 4 | 0 | 0 | 3 | 1 of 6 |
| `kd` | `NSurfNegMhp` | 2 | 4 | 4 | 2 | 1 of 6 |
| `kd` | `SurfMhpPatchArea` | 4 | 0 | 0 | 4 | 0 of 6 |
| `kd` | `SurfMhpLocalPosSum` | 4 | 0 | 0 | 4 | 0 of 6 |
| `cw` | `SurfMhpMean` | 3 | 0 | 0 | 0 | 3 of 6 |
| `cw` | `SurfMhpSum` | 6 | 0 | 0 | 6 | 5 of 6 |
| `cw` | `NSurfPosMhp` | 6 | 0 | 0 | 4 | 5 of 6 |
| `cw` | `NSurfNegMhp` | 0 | 5 | 4 | 2 | 4 of 6 |
| `cw` | `SurfMhpPatchArea` | 5 | 0 | 1 | 4 | 4 of 6 |
| `cw` | `SurfMhpLocalPosSum` | 6 | 0 | 0 | 6 | 5 of 6 |
| `mj` | `SurfMhpMean` | 1 | 0 | 4 | 0 | 0 of 6 |
| `mj` | `SurfMhpSum` | 3 | 3 | 5 | 1 | 1 of 6 |
| `mj` | `NSurfPosMhp` | 3 | 3 | 3 | 3 | 0 of 6 |
| `mj` | `NSurfNegMhp` | 0 | 0 | 0 | 0 | 0 of 6 |
| `mj` | `SurfMhpPatchArea` | 0 | 0 | 0 | 0 | 0 of 6 |
| `mj` | `SurfMhpLocalPosSum` | 3 | 3 | 5 | 1 | 1 of 6 |

On `mj_scaled` and `cw` the shift is smaller under the new value for most structures, for `SurfMhpSum` on 5 of 6. On `kd` it is larger, for `SurfMhpSum` on every structure: 1.0 lies near the middle of `kd`, so the old `OXT` contributed little, while the ASP value of -3.5 is near its hydrophilic end.

So the value of `OXT` does not by itself make structures with and without one comparable. A structure source that does not write `OXT`, such as Boltz2, still gives slightly different hydrophobicity features from one that does. Preparing every structure, as in [preparing a structure](structure_preparation.md), gives every chain its `OXT` and removes that difference.
