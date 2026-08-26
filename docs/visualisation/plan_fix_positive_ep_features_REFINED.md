# Fixing the surface electrostatic potential features: refined plan

**Status:** plan only. **Branch:** `visualisation`. **Date:** 2026-08-21.

**Supersedes** `plan_fix_positive_ep_features.md`, which should not be implemented from: its recommendation was backwards. **Provenance:** original plan reviewed independently by a statistics reviewer (`_REVIEW1.md`) and a computational chemistry reviewer (`_REVIEW2.md`), both of whom recomputed rather than read. Every contested claim was then re-checked directly. Verified on 2026-08-21 is tagged **[verified]**.

**What changed.** The original recommended patching the threshold (option A) and evaluating Debye screening (option C) later. That is exactly backwards: **screening wins outright, with zero fitted constants, and beats even an oracle threshold.** Option A is dropped. The original's claimed repair also turned out to be largely a protein-size correlation. The memory model was wrong, the two reviewers disagreed about it by a factor of two, and the disagreement is resolved below.

---

## 1. The decision

Add a Debye screening factor to the Prodes Coulomb kernel. One line in `Property_point.set_ep`:

```python
ep = sum( q * exp(-d / DEBYE_LENGTH) / (4 * pi * eps0 * 4 * d) )
```

Measured against APBS at identical surface coordinates on all 13 proteins **[verified, reproduced independently of the review]**:

| scheme | fitted constants | median Jaccard overlap with the APBS positive surface |
| --- | --- | --- |
| shipped, `ep > 0` | 0 | **0.058** |
| option A, recalibrated threshold | 2 | 0.586 |
| oracle, thresholded at the true APBS zero | per protein, from APBS | 0.601 |
| option B, `ep > mean(ep)` | 0 | 0.618 |
| **option C, screened, `> 0` unchanged** | **0** | **0.658** |

Screening with no fitted parameter beats the oracle threshold. Cross-protein agreement with APBS:

| feature | shipped | option A | **option C** |
| --- | --- | --- | --- |
| `NSurfPosEpFormal` | 0.641 | 0.973 | **0.995** |
| `SurfEpPosMeanFormal` | 0.234 | 0.335 | **0.813** |
| `SurfEpPosStdFormal` | 0.123 | 0.126 | **0.577** |
| `SurfEpStdFormal` (global) | 0.060 | unchanged | **0.516** |

Option A repairs one feature and leaves the rest; option C improves the negative-side and global features too, because it fixes the field rather than the cutoff.

**Insensitive to the parameter.** Median Jaccard is 0.651 at lambda = 5 A, 0.658 at 7.85 A, 0.640 at 12 A **[verified]**. The "it adds a parameter that must be chosen" objection in the original plan does not survive contact with the data.

**Control:** an unscreened reimplementation reproduces the stored Prodes potentials to **0.0000 V** on all 13 **[verified]**, which is what licenses the comparison.

## 2. Why option A was rejected

**Its headline repair was a size correlation.** `n_surface_points` alone, with no electrostatics whatsoever, scores **0.929** against the APBS positive count **[verified]**. Option A's 0.973 is a gain of +0.05 over that, bootstrap 95 % CI [-0.10, +0.24], P(no gain) = 0.22. On the size-free positive **fraction** the recalibrated rule scores 0.769 against the naive 0.780, that is slightly **worse**, while the oracle reaches 0.885 **[verified]**. The information exists; a two-parameter linear rule does not extract it.

Leave-one-out gives 0.956 against 0.973 in-sample, so overfitting was not the problem. The rule was measuring the wrong thing.

Worse, both fixes push the correlation of `NSurfPosEpFormal` with `NSurfPoints` from R2 0.242 to 0.743 under A, turning it into the size proxy that the repo's own comment in `calculate_surface_grid_features` deleted `NSurfNegEpFormal` for being.

## 3. Honesty about what the screening factor is

`lambda = 3.04 / sqrt(I)` Angstrom is the Debye length **in water**, at eps_r = 78.54. The Prodes kernel uses eps_r = 4, where the self-consistent Debye length would be **1.77 A**. Applying a water Debye length inside a kernel with a protein dielectric is therefore **an empirical correction that mimics screening, not a first-principles Poisson-Boltzmann calculation**, and the documentation must say so. It works, and the result is flat across lambda, but it should not be described as making Prodes physically correct.

The strictly correct linearised Debye-Huckel form is `exp(-kappa(r-a)) / ((1 + kappa*a) * r)`. The extra terms form a distance-independent prefactor `exp(kappa*a)/(1 + kappa*a)` = 1.0282 at a = 2 A, and the measured effect on every metric here is exactly 0.000. Include it for correctness; expect no accuracy change.

## 4. Corrections to the physical argument

The original justified the offset as the monopole term `Q/R`, using `sqrt(n_surface_points)` as a proxy for R. That proxy tracks the radius of gyration only to about 31 %, and the physically correct harmonic-mean radius actually fits **worse** (R2 0.944) than the proxy (0.991), which is a curve-fitting signature at n = 13.

The clean argument, which should replace it: for a harmonic function, the average over a sphere equals the value at its centre, so the mean of `1/d` over a spherical surface is exactly `1/R`. The monopole reading is right; the `sqrt(n)` proxy was doing more work than it deserved.

Also correct two overstatements. The Prodes-to-APBS map is affine with a **slope that varies twofold across proteins** (0.035 to 0.072), not a pure offset. And "an affine map reproduces 73 % of the APBS variance" is the **best** of the 13; the median is 54 %.

## 5. Scope: 10 features, not 3

Sign-thresholded EP features in the default 54:

- surface: `NSurfPosEpFormal`, `SurfEpPosMeanFormal`, `SurfEpPosStdFormal`, `SurfEpNegMeanFormal`, `SurfEpNegStdFormal`
- shell: `NShellPosEpFormal`, `ShellEpPosMeanFormal`, `ShellEpPosStdFormal`, and the negative pair

The shell features use the same unscreened kernel on far-field points, **where the monopole dominates more**, and no APBS ground truth was sampled out there at all. Sampling the shell points is a required addition to the measurement.

**The negative side is not inherently safe.** That asymmetry is a property of the 13-protein test set, all of which are net negative. `tests/data/1GDW` (net charge +7, pI 8.9) shows the mirror failure already, `NSurfPosEpFormal` 7020 of 7039 and `NSurfNegEpFormal` 18. Any claim that "anion exchange is fine" holds only for net-negative proteins and must be stated that way.

## 6. What the published models actually use

From Neijenhuis (2024), scored against APBS on the 13 **[verified]**:

| feature | Spearman vs APBS | role |
| --- | --- | --- |
| `SurfEpMinFormal` | **-0.236** | cation exchange model, top feature, coefficient 31.24 |
| `SurfEpNegMedianFormal` | 0.758 | anion exchange model, #1, coefficient -31 |
| `NSurfPosEpFormal` | 0.641 | anion exchange model, #2, coefficient 18.17 |
| *pure size baseline* | *0.929* | *no electrostatics* |

Both published models lean on features that disagree with their APBS counterparts across proteins. Caveat: n = 13, and these are not the same proteins the models were fitted on.

Separately, the cation exchange model reports a cross-validated R2 of 0.95 from **10 features on 72 datapoints drawn from 16 proteins at 5 pH values**. The paper invokes a "five datapoints per feature" rule to justify 3 features for the anion model and then uses 10 for the cation model. Structure-derived features barely change with pH, so the pH replicates are closer to pseudo-replication than to independent samples. The paper itself notes that total charge takes a negative coefficient in a cation exchange model, which is backwards, and attributes it to collinearity.

## 7. Data availability, which constrains everything

**There is no cation exchange retention data in this repository** **[verified]**. `08_neijenhuis/retention_times.csv` is 52 rows, all HiTrap Q Sepharose XL, an anion exchanger, 13 proteins at 4 pH values. The Disela set is anion-mode IEX with 94.6 % of proteins net negative at measurement pH. The SP Sepharose HP data behind the cation model would have to be obtained from the source the paper cites.

So the end-to-end test the original plan called "the only test that measures what anyone actually cares about" **cannot be run here**. It should be listed as blocked pending that data, not scheduled.

## 8. Validation

1. **Build the group split.** `11_protein_clusters` cannot support it: it holds 2 clusters covering 4 of 834 proteins and is a duplicate detector, not homology clustering **[verified by both reviewers]**. Use MMseqs2 or CD-HIT at 30 % identity.
2. **Report the size-only baseline next to every number.** `n_surface_points` scores 0.929 on the count feature. Any proposed fix that does not clearly beat that has not been shown to work.
3. **Score size-free variants** (positive *fraction*, not count) as the primary endpoint, since the count is size-confounded by construction.
4. **Sample the shell points too**, or the shell features are being changed without evidence.
5. **Include both charge signs.** Only 67 of 834 AlphaFold structures are net positive at pH 7, about 8 % **[verified]**, so a held-out basic subset is roughly 13 proteins, the same n this plan calls insufficient. **Run a subset at pH 5, where about 55 % are net positive.** That is the cheapest way to get coverage of the regime the fix exists to serve, and the original plan missed it.
6. **Include 1GDW** as the mirror-failure regression case alongside the four zero-positive proteins.
7. Screening changes every feature value, so **all published figures and models must be regenerated**, and this is a major version bump.

## 9. Compute cost, corrected

**The memory disagreement is resolved.** Reviewer 1 read APBS's own log counter and got 467.5 MB per million grid points; reviewer 2 measured resident memory externally and got 224.4. I measured it directly: on a 8.0 Mpt grid, APBS self-reports 3,755 MB high water while actual peak RSS is **1,845 MB**, that is **231 MB per Mgridpoint** **[verified]**. APBS's internal counter roughly doubles reality. Reviewer 2 is right.

Consequences: the largest structure (`ARH96700`, 577x385x481 = 106.9 Mpts) needs about **25 GB**, not the 81 GB implied by the log counter and not the 28 GB the original plan asserted. **Memory is not the binding constraint and the memory-aware scheduler in the original deliverables is unnecessary.** Ordinary N-way parallelism is fine on 113 GB, with the single largest structure given room.

Runtime stands: `seconds = 5.35 * Mgridpoints + 0.9`, R2 0.999, about **9.8 h single core, 1.2 h on 8 processes** for 834 structures.

**The expensive tail is AlphaFold disorder, not protein size.** Among the 25 largest grids the median fraction of residues below pLDDT 70 is 0.16 against 0.017 overall; one 156-residue model needs a 21.6 M grid. Dropping structures with more than 10 % low-confidence residues removes 106 of 834 (12.7 %), saves 21.6 % of runtime and cuts peak memory from 23.6 to 5.4 GiB. Do that, and record it as an exclusion rather than a silent filter.

Whether AlphaFold models and crystal structures belong in one calibration is a real question; stratify and check rather than pooling blind.

## 10. Implementation

No two-pass restructuring is needed. `calculate_surface_grid_features` already materialises the full `eps` array in the parent after `run_tasks` returns, so neither `process_surface_grid_cell` nor `SURFACE_GRID_STATE` changes.

The screening factor goes in `Property_point.set_ep` behind a module constant, defaulting to off until validated, with ionic strength as the exposed parameter.

**Trap to avoid:** with option A's negative threshold, `eps > t` and `eps < 0` overlap, breaking the documented `NSurfNeg = NSurfPoints - NSurfPos` partition invariant. Option C keeps the cutoff at 0 and preserves it. One more reason to prefer C.

**On `standard_features` returning 0 for an empty array:** leave it. All 210 tests pass with NaN, but exactly one call site in the suite ever sees an empty array and it would then pass for the wrong reason. Scope is 39 of 105 features, and NaN is wrong for the `Sum` family. More importantly, in the production data the censored rows have median net charge -26.7 against -7.5 overall, so 0 places them at the correct extreme; NaN would be worse.

## 11. Files touched

| File | Change |
| --- | --- |
| `src/prodes/core/point.py` | screening factor in `set_ep`, behind a constant, default off |
| `src/prodes/run.py` | plumb ionic strength through `calculate` |
| `scripts/apbs_comparison.py` | sample shell points; structure-list input; exclusion of low-pLDDT models |
| `scripts/validate_screening.py` | group split, size baseline, both charge signs, pH 5 subset |
| `docs/visualisation/screening_validation.md` | results |
| `tests/test_screening.py` | 1GDW mirror case, the four zero-positive proteins, lambda insensitivity |
