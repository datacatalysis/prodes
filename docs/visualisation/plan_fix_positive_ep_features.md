# Plan: fix the positive surface electrostatic potential features

**Status:** plan only, for review. **Branch:** `visualisation`. **Date:** 2026-08-21.

Follows from `apbs_vs_prodes_results.md`. Every number below tagged **[measured]** was produced on the analysis machine on 2026-08-20 or 2026-08-21 and is reproducible from `scripts/apbs_comparison.py`.

---

## 1. The problem, in one paragraph

Prodes computes the surface electrostatic potential as an unscreened Coulomb sum with a uniform dielectric of 4. Compared against an APBS Poisson-Boltzmann solution sampled at identical coordinates, the two agree **up to an affine transform**: on beta-lactoglobulin an affine map of the Prodes potential reproduces 73 % of the APBS variance **[measured]**, and after rescaling the colour ramp the two point clouds are visually near identical, patch for patch. So for any monotone feature (mean, min, max, standard deviation) the discrepancy is largely harmless.

It is not harmless for features defined by a **threshold**, because a threshold is the one operation an affine shift does not survive. Prodes selects its positive population with `ep > 0`, but the Prodes value at which APBS actually crosses zero ranges from **-0.40 V to -3.68 V** across the 13 test proteins **[measured]**. The consequence, measured as set overlap between the points Prodes calls positive and the points APBS calls positive:

| | median Jaccard overlap |
| --- | --- |
| Prodes, naive `ep > 0` | **0.058** |
| thresholded at the APBS zero crossing | 0.601 |

0.058 means the two sets are very nearly **disjoint**. Three features in the default 54-feature set are built on that population, `NSurfPosEpFormal`, `SurfEpPosFormalMean` and `SurfEpPosFormalStd`, and `standard_features` returns 0 for an empty array, so all three collapse to exactly 0 for four of thirteen proteins where APBS counts 489 to 7,754 positive points **[measured]**.

The damage is asymmetric, and the asymmetry is the point:

| feature | Spearman against the APBS equivalent, 13 proteins | relevant to |
| --- | --- | --- |
| `NSurfNegEp` | 0.956 | anion exchange |
| `SurfEpNegMean` | 0.571 | anion exchange |
| `NSurfPosEp` | 0.641 | **cation exchange** |
| `SurfEpPosMean` | 0.234 | **cation exchange** |
| `SurfEpPosStd` | 0.123 | **cation exchange** |

Anion exchange is driven by negative surface, which Prodes describes well. Cation exchange is driven by positive surface, which Prodes describes through a mis-placed cutoff.

## 2. How the threshold is picked, since this was asked directly

**There is one rule, fitted once, and it produces a different threshold for every protein.**

The rule currently reads

```
threshold_volts = 0.743 * (this protein's mean surface EP) - 0.121
```

The two constants `0.743` and `-0.121` are global: fitted once across the whole calibration set and then frozen. What varies per protein is the input, that protein's own mean surface potential, and therefore the output. Catalase gets -0.40 V, pepsin gets -3.68 V.

**Why it has to be per protein.** The offset between Prodes and APBS is the unscreened monopole term. The potential at the surface of a sphere of radius R carrying net charge Q is proportional to `Q/R`, and because Prodes applies no ionic screening, that term survives in full and shifts the entire surface of that protein by a single amount. A different protein has a different `Q` and a different `R`, so it gets a different shift. A single global threshold cannot work, and neither can `0`.

This is confirmed rather than assumed. The Prodes mean surface potential is **99.1 %** explained by `netQ / sqrt(n_surface_points)`, which is exactly `Q/R` **[measured]**; the same regression against APBS gives 68.3 %. And the correct threshold is predictable from quantities Prodes already has **[measured]**:

| predictor | R2 with the correct threshold |
| --- | --- |
| net charge Q | 0.923 |
| Q / sqrt(n_points) | 0.940 |
| the protein's own mean surface EP | **0.954** |

Using the last of those, with no APBS at run time:

| `NSurfPosEp`, Spearman against APBS | |
| --- | --- |
| naive `ep > 0` | 0.641 |
| at the true APBS zero crossing | 0.995 |
| **at the predicted threshold, no APBS needed** | **0.973** |

**The honest caveat, and the reason this document exists:** those constants were fitted and evaluated on the same 13 proteins, with no held-out data. That is not a calibration, it is a demonstration that the breakage is systematic and correctable.

## 3. Are 13 proteins enough? No.

Three independent reasons.

**Statistical.** A two parameter fit on 13 points is not a calibration. At n = 13 a Spearman of 0.6 has a 95 % interval of roughly 0.05 to 0.87, which is the width of the entire conclusion. Every cross-protein number in section 1 sits inside that uncertainty. Only the within-protein results and the 0.991 monopole regression are tight enough to stand on 13.

**Coverage.** All 13 are net negative, from -5 to -35 e. The proteins the rule most needs to work on are the ones with positive or near-neutral net charge, since those are the ones that bind in cation exchange, and **the calibration set contains none of them**. Extrapolating a linear rule past the edge of its training range is exactly the failure mode to expect.

**Composition.** The 13 are large, mostly multi-domain, four are multi-chain and two carry covalently modified residues that both tools discard. They are not representative of the host cell proteins the models are ultimately applied to.

The fix must therefore be fitted on a set that is larger, spans both signs of net charge, and is drawn from the same population as the application.

## 4. What it would take to run the Disela dataset

`../biochai/data/04a_alphafolddb_structures` holds **834 AlphaFold structures**, each with a matching `_pka.json` in `04b_compute_surface_pka` **[measured]**. With the 13 Neijenhuis proteins that gives **847**.

Sizes are modest compared to the Neijenhuis set: median 2,562 atoms, p95 6,082, max 11,962 **[measured]**.

**Runtime, calibrated on six real runs spanning the size range** **[measured]**:

| structure | atoms | surface points | grid | APBS | total |
| --- | --- | --- | --- | --- | --- |
| ARH98169 | 811 | 6,825 | 129x129x129 | 9 s | 11 s |
| ARH98911 | 2,562 | 22,389 | 257x161x225 | 41 s | 48 s |
| ARH97435 | 6,135 | 33,238 | 225x193x225 | 44 s | 56 s |
| ARH96666 | 10,348 | 128,695 | 385x321x321 | 176 s | 213 s |

Runtime tracks grid volume almost exactly, `seconds = 5.35 * (grid points in millions) + 0.9`, R2 = 0.999 on the six **[measured]**. Extrapolated across the real size distribution of all 834:

| | |
| --- | --- |
| median per structure | 28 s |
| p95 | 111 s |
| slowest | 928 s |
| **total, one core** | **9.8 h** |
| on 4 processes | 2.4 h |
| on 8 processes | 1.2 h |

APBS is single threaded, so parallelism is just N independent processes. The machine has 16 logical CPUs, 113 GB free RAM and 352 GB free disk **[measured]**.

**Two constraints that need designing around.**

*Memory, not CPU, is the binding limit.* The largest structure needs a 385x321x321 grid, roughly 28 GB of resident memory by the rough estimate of about 20 grid-sized arrays. Eight of those at once would exhaust the machine. Schedule by predicted grid volume, not by count: run the small ones many-at-once and the largest few one-at-a-time. Predicted grid volume is available before the run from the bounding box, so this can be a static schedule.

*Disk, if the grids are kept.* At a mean of about 100 MB per `.dx`, 847 grids is roughly 85 GB. The existing pipeline already samples and deletes immediately, keeping one reference grid, so peak disk stays at one grid. Do not change that.

**The 5 % largest account for 23 % of the total time** and only 12 structures have a bounding box over 150 Angstrom **[measured]**. Dropping those saves 11 % of runtime; not worth the loss of coverage, but useful if a first pass is wanted quickly.

**Estimate: one overnight run, or about 2 hours with a memory-aware scheduler.** This is comfortably affordable and there is no reason to settle for 13.

## 5. Options for the fix

Three, in increasing order of ambition. They are not exclusive; the validation in section 6 should score all three.

### Option A: recalibrate the threshold

Keep the physics, move the cutoff. Add a per-protein threshold computed from the rule in section 2, and select the positive population with `ep > threshold` instead of `ep > 0`.

*Pros.* Smallest possible change. Fixes the count feature nearly completely (0.641 to 0.973 on the pilot). Requires no external tool at run time. Leaves every other feature bit-identical, so nothing already published moves except the three positive features.

*Cons.* Only partly repairs `SurfEpPosMean` (0.234 to 0.335): recentring fixes *which* points are selected but the mean of a censored tail is still sensitive to the scale. Introduces a fitted empirical constant into a physics calculation, which needs documenting honestly. Changes three released feature values, so any model trained on them must be refitted.

### Option B: subtract the monopole term

Report `ep - mean(ep)`, or better, subtract the fitted `Q/R` term analytically, so the potential is centred the way APBS is. Then `> 0` becomes correct again by construction.

*Pros.* Fixes the offset for **every** feature at once rather than patching three. Makes `> 0` meaningful, so no empirical threshold constant is needed in the feature code. Cheap.

*Cons.* Destroys the absolute scale, which is currently the thing carrying the net-charge signal, and that signal is genuinely predictive. Risk of removing real information: a protein that is uniformly negative *is* different from one that is not. Mitigate by keeping the mean as its own explicit feature so the information is retained rather than discarded, which is arguably better practice anyway.

### Option C: add Debye screening to the kernel

Fix the physics at source. Multiply each term by `exp(-d / lambda)` with `lambda = 3.04 / sqrt(I)` Angstrom, 7.85 Angstrom at 150 mM. This is a one line change to `Property_point.set_ep`.

*Pros.* Addresses the actual cause rather than its symptom, and should improve every feature simultaneously, not just the thresholded ones. Screening is the single largest physical omission. Introduces ionic strength as a real, meaningful parameter, which matters because chromatography is run at defined salt. Still no external dependency and still O(charges) per point.

*Cons.* Changes **every** feature value Prodes has ever produced, so every published number and figure would need regenerating. Adds a parameter that must be chosen, and the papers' data span several buffer conditions. Does not reproduce a full dielectric boundary, so it will not match APBS exactly. Needs its own validation before anyone would trust it.

**Recommendation: implement A first**, because it is contained and the count feature is the one the repo's own redundancy analysis identified as carrying signal. **Evaluate C in the same validation run**, because if screening alone recovers most of the agreement it is the better long term answer and A becomes unnecessary.

## 6. Validation design

This is the part the pilot got wrong and must not repeat.

1. **Split before fitting.** Group split by sequence similarity, not random: `../biochai/data/11_protein_clusters` already exists and should be used, or the calibration will leak between homologues. Fit on the training clusters, report on held-out clusters only.
2. **Report the naive baseline alongside every number.** The claim is "0.641 becomes 0.973", so both must appear, on held-out data.
3. **Stratify by net charge.** The rule's whole justification is a `Q/R` dependence, so its residual must be examined against `Q`. If it degrades for `Q > 0` proteins, that is a blocking result and it will not be visible in an aggregate.
4. **Include the four zero-positive proteins as a named test case.** They are the failure that motivated this, and they must be shown to be fixed.
5. **Compare the three options on the same split**, scored by agreement with APBS on the same held-out clusters.
6. **The end-to-end test, if retention data allows it.** Recompute the surface features from the APBS potentials and refit the retention models. If APBS-derived features beat Prodes-derived features on **cation exchange** and tie on **anion exchange**, that confirms the mechanism exactly where this analysis predicts. This is the only test that measures what anyone actually cares about, and it should be attempted rather than deferred.

## 7. Deliverables

| File | Purpose |
| --- | --- |
| `scripts/apbs_comparison.py` | extend to take a structure list, add memory-aware scheduling |
| `scripts/calibrate_ep_threshold.py` | fit and validate on held-out clusters, write the constants |
| `src/prodes/calculations/...` | the chosen fix, behind a flag until validated |
| `tests/test_ep_threshold.py` | pin the rule and the zero-positive regression cases |
| `docs/visualisation/ep_threshold_calibration.md` | the results, with the held-out numbers |

## 8. Open questions for review

- Is option C the right primary target, making A a stopgap?
- Should the threshold rule use mean EP, `Q/sqrt(n)`, or `Q` directly? Mean EP fits best on 13 but is itself derived from the potential being corrected, which is a circularity worth a second opinion.
- Is `standard_features` returning 0 for an empty array the right behaviour at all, or should it be NaN so that a censored feature is visible to a model rather than silently equal to a real zero?
- Does changing three released feature values require a major version bump?
