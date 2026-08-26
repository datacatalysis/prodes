# Independent review: `plan_fix_positive_ep_features.md`

**Reviewer brief:** computational chemistry and software implementation. **Date:** 2026-08-21. **Machine:** the analysis machine. **Envs:** `prodes` (Prodes + numpy/pandas), `prodes-apbs` (APBS 3.4.1, pdb2pqr 3.6.1).

Everything below tagged **[measured]** was produced by running code during this review, not by reading the plan. Everything tagged **[verified]** is a number the plan already claims that I reproduced independently. Working scripts live in the session scratchpad (`screening.py`, `analyse.py`, `optb.py`, `gridpred.py`, `netq.py`, `nanpatch.py`) and are quoted inline where the exact form matters.

---

## Verdict

The diagnosis is right and reproduces exactly. Every `[measured]` number in section 2 came back to three decimal places when I recomputed it from the structures. The plan is honest about its own weaknesses and section 3 is correct that 13 proteins is not a calibration.

The **recommendation is backwards**, and I can show it rather than argue it. Option C, the one the plan calls least tested and defers, is the best of the three on every metric I could construct, needs no fitted constant, and beats the *oracle* threshold that uses the real APBS zero crossing. Option A, the one the plan wants to implement first, is beaten by the trivial one-line Option B (`ep > mean(ep)`, zero fitted parameters). Implementing A first would spend a release, a model refit and a version bump on the weakest of the three.

Separately, the headline evidence for A (`NSurfPosEp` Spearman 0.641 to 0.973) is **confounded with protein size**: on these 13 proteins, the surface point count alone scores 0.929 against APBS with no electrostatics involved at all. Both proposed fixes make `NSurfPosEpFormal` substantially *more* of a size proxy than it currently is, which by the repo's own redundancy criterion is the argument for deleting the feature, not repairing it.

The compute plan is affordable but its memory model is wrong in both directions at once, its named "largest structure" is not the largest, and the group split it depends on does not exist in the data directory it points at.

Recommended sequencing: **C first, behind a flag, validated on the 834; B as the cheap interim; A only if C fails.** Do not start the 847-structure run until findings F3 and F7 are settled.

---

## F1. The recommendation is backwards. Screening alone beats every thresholding option, including the oracle. [measured]

**What the plan says.** Section 5: "Recommendation: implement A first, because it is contained... Evaluate C in the same validation run, because if screening alone recovers most of the agreement it is the better long term answer and A becomes unnecessary." Option C is listed last, its cons are the longest, and it is described as needing "its own validation before anyone would trust it."

**What is actually true.** Screening alone does recover most of the agreement, and more than that: it beats the best possible thresholding, so A is already unnecessary on the evidence that exists today. I recomputed the Prodes potential from the structures with an added `exp(-d/lambda)` factor and scored the result against the same APBS values in `apbs_results/<id>/<id>_pdie4_points.csv.gz`.

Median Jaccard overlap of the positive-surface point set against APBS, 13 proteins:

| variant | fitted constants | median Jaccard |
| --- | --- | --- |
| Prodes today, `ep > 0` | 0 | 0.058 |
| **A**, threshold `0.743 * meanEP - 0.121` | 2 | 0.586 |
| **A oracle**, threshold at the true APBS zero crossing (needs APBS) | 13 | 0.601 |
| **B**, `ep > mean(ep)` | 0 | 0.618 |
| **C**, screened at lambda = 7.85 A, `ep > 0` unchanged | 0 | **0.658** |

Option A with two fitted constants scores *below* the parameter-free Option B. Option C with no fitted constants scores above the oracle threshold that is allowed to look at the APBS answer for each protein individually. That is the single most important result in this review.

Cross-protein Spearman against the APBS equivalent, the same table the plan uses to make its case, n = 13:

| feature | Prodes today | Option A (predicted threshold) | **Option C (lambda 7.85)** |
| --- | --- | --- | --- |
| `NSurfPosEp` | 0.641 | 0.973 | **0.995** |
| `SurfEpPosMean` | 0.234 | 0.335 | **0.813** |
| `SurfEpPosStd` | 0.123 | not reported | **0.577** |

Option C also repairs everything A leaves alone, which is the claim the plan makes for C but does not test:

| feature | Prodes today | Option C |
| --- | --- | --- |
| `NSurfNegEp` | 0.956 | 0.984 |
| `SurfEpNegMean` | 0.571 | 0.720 |
| `SurfEpNegStd` | 0.192 | 0.335 |
| `SurfEpMean` | 0.582 | 0.819 |
| `SurfEpStd` | 0.060 | 0.516 |
| `SurfEpMax` | 0.346 | 0.467 |
| `SurfEpMin` | -0.236 | -0.071 |

No feature gets worse. Median point-level Spearman rises from 0.779 to 0.816 (0.834 at lambda 15), and the median point-level affine R2 from 0.537 to 0.576. The four proteins that report zero positive surface (1CF3, 1F6R, 4PEP, 6FRV) all report a nonzero, roughly correct positive fraction under screening: 1CF3 0.264 against APBS 0.261, 1F6R 0.313 against 0.264, 4PEP 0.016 against 0.035, 6FRV 0.094 against 0.129. The median absolute error in the positive surface fraction drops from 0.271 to 0.031.

The lambda sweep shows the result is not sensitive to the choice, which removes the plan's "adds a parameter that must be chosen" objection to C:

| lambda (A) | median Jaccard | median Spearman | `NSurfPosEp` Spearman | `SurfEpPosMean` Spearman |
| --- | --- | --- | --- | --- |
| 5 | 0.651 | 0.795 | 0.989 | 0.808 |
| **7.85 (150 mM)** | **0.658** | 0.816 | 0.995 | 0.813 |
| 10 | 0.654 | 0.827 | 0.995 | 0.846 |
| 12 | 0.640 | 0.832 | 0.995 | 0.890 |
| 15 | 0.599 | 0.834 | 0.989 | 0.901 |
| 20 | 0.540 | 0.831 | 0.967 | 0.907 |

Anything from 5 to 12 A is within noise of the optimum on the thresholded metrics. The physically motivated value at the condition APBS was run at is inside that window.

**Evidence.** Prodes formal charges were regenerated per structure with `prodes.run.prepare_structure` at pH 7 using the same `_pka.json` files `scripts/apbs_comparison.py` uses, and the potential recomputed at the exact coordinates stored in the per-point CSVs. The unscreened reimplementation reproduces the stored `ep_prodes_volts` column exactly (`max |round(recomputed,2) - stored| = 0.0` for all 13 proteins), which is what licenses the comparison. Kernel:

```python
d = np.linalg.norm(points[:, None, :] - charge_coords[None, :, :], axis=2)
f = np.exp(-kappa * d)                       # kappa = 1/lambda, Angstrom^-1
v = (q * 1.6e-19 * f / (4*np.pi*8.854e-12*4.0*d*1e-10)).sum(axis=1)
```

**Concrete fix.** Reverse section 5. Make C the primary target, behind a flag, and validate it on the 834. Keep B as the interim for anyone who cannot regenerate published numbers. Drop A unless C fails on held-out data, and if A survives at all, note in section 5 that it is beaten by B.

---

## F2. "0.641 to 0.973" is largely a protein-size correlation, and both fixes turn the feature into a size proxy. [measured]

**What the plan says.** Sections 1 and 2 rest on `NSurfPosEp` Spearman 0.641 rising to 0.973. Option A's pros say it "fixes the count feature nearly completely" and that this is "the one the repo's own redundancy analysis identified as carrying signal."

**What is actually true.** On these 13 proteins, a feature with no electrostatics in it at all scores 0.929:

```
Spearman(n_surface_points, APBS positive count) = 0.929
Spearman(-netQ,            APBS positive count) = -0.132
R2      (n_surface_points, APBS positive count) = 0.938
```

The APBS positive count on this set is itself 94 per cent explained by how big the protein is. The 13 span 8,466 to 58,743 surface points, so any count-type feature is dominated by size, and 0.641 to 0.973 has to be read against a 0.929 floor, not against zero. At n = 13 the difference between 0.929 and 0.973 is not resolvable, which is the plan's own section 3 argument turned on its own headline.

Worse, both fixes push `NSurfPosEpFormal` toward being a pure size proxy:

| `NSurfPosEpFormal` | R2 with `NSurfPoints` |
| --- | --- |
| today, `ep > 0` | 0.242 |
| Option A, predicted threshold | 0.743 |
| Option C, screened | 0.902 |
| APBS reference | 0.938 |

The code comment in `calculate_surface_grid_features` keeps `NSurfPosEpFormal` and drops `NSurfNegEpFormal` on exactly this criterion: "NSurfNegEpFormal is also the size proxy of the pair (R2 0.87 with Area, vs 0.03 for NSurfPosEpFormal), so NSurfPosEpFormal is the one that carries signal." Repairing the feature moves it from R2 0.24 to 0.74 or 0.90 against `NSurfPoints`, which is the condition under which the same analysis deleted its sibling.

This does not mean the fix is wrong. It means the *justification* is wrong: the repaired feature is a better estimate of positively charged surface area, and positively charged surface area really is close to total surface area times a modest factor. The honest framing is that the fix makes the feature correct and simultaneously makes it redundant with `Area`, and the interesting descriptor is the *fraction* (`NSurfPosEpFormal / NSurfPoints`), not the count.

**Evidence.** Computed from `apbs_results/<id>/<id>_pdie4_points.csv.gz` and the recomputed screened potentials; ranks with the same tie-aware `average_ranks` implementation `scripts/apbs_comparison.py` uses, so the numbers are comparable to the published ones.

**Concrete fix.** Add the size baseline to every table in sections 1, 2 and 6: "Spearman(n_surface_points, APBS positive count) = 0.929". Add `NSurfPosEpFormal / NSurfPoints` as the feature actually being proposed, and re-run the repaired features through the existing feature reduction rather than assuming their old redundancy verdict still holds. Section 6 needs a seventh item: "re-run the redundancy analysis on the repaired features."

---

## F3. The group split the validation design depends on does not exist. [measured]

**What the plan says.** Section 6, item 1: "Group split by sequence similarity, not random: `../biochai/data/11_protein_clusters` already exists and should be used, or the calibration will leak between homologues."

**What is actually true.** `11_protein_clusters/similarity-clusters.json` defines **2 clusters covering 4 proteins**, out of 834. They are `lacI`/`lacI` at 100 per cent identity and `tufA`/`tufB` at 99.75 per cent. The directory is a near-duplicate detector built from ESM-C embedding cosine similarity, not a homology clustering, as `duplicate-candidates.json` and a 3-line `duplicate_pair_identities.csv` confirm. Using it as a group split would put 830 of 834 proteins in singleton groups, which is a random split with two pairs held together.

```
$ python3 -c "import json; d=json.load(open('.../11_protein_clusters/similarity-clusters.json')); \
  print(len(d), sum(c['member_count'] for c in d))"
2 4
```

The raw material does exist: `cosine-similarity-matrix.json` is a 27 MB full pairwise matrix and could be thresholded into real clusters, or mmseqs2 could be run on the sequences in `03_uniprot_sequence_lookup`.

**Concrete fix.** Change section 6 item 1 to: "build a homology clustering first, by thresholding `11_protein_clusters/cosine-similarity-matrix.json` or by running mmseqs2 at 30 per cent identity on `03_uniprot_sequence_lookup`; `similarity-clusters.json` is a duplicate detector and covers 4 of 834." Add it to the deliverables table. This is a prerequisite for the compute run being worth anything.

---

## F4. The screening kernel as specified is dielectrically incoherent, and the plan should say so. [measured]

**What the plan says.** Option C: "Multiply each term by `exp(-d / lambda)` with `lambda = 3.04 / sqrt(I)` Angstrom, 7.85 Angstrom at 150 mM. This is a one line change to `Property_point.set_ep`."

**What is actually true.** `lambda = 3.04 / sqrt(I)` is the Debye length of *water*, derived from `kappa^2 = 2 N_A e^2 I / (eps_r eps_0 k_B T)` with `eps_r = 78.54`. `Property_point.set_ep` uses `eps_r = 4`. Computing kappa at the dielectric the kernel actually uses gives:

```
Debye length at eps_r = 78.54 : 7.86 A     (the plan's 7.85, [verified])
Debye length at eps_r = 4.0   : 1.77 A
```

A self-consistent Debye-Huckel kernel at `eps_r = 4` would have a 1.77 A screening length, which annihilates the potential a couple of Angstrom from each charge and would not reproduce anything. So the proposal is not "fix the physics at source": it is an empirical solvent-screening kernel with a water Debye length bolted onto a protein-dielectric Coulomb sum. It works very well (F1), but the plan must not present it as first-principles Poisson-Boltzmann, because the next reader will try to justify `eps_r = 4` and `lambda = 7.85` in the same sentence and will not be able to.

There is a clean way out that costs nothing. The 20x magnitude error the results document reports is exactly `78.54 / 4 = 19.6`, and a uniform rescale changes no rank, no sign, no threshold and no Spearman. Setting `eps_r = 78.54` alongside the screening makes the kernel a coherent screened-Coulomb model of a charge in bulk water, puts the output in physically interpretable volts, and costs exactly one extra constant change with zero effect on any metric in this review.

**Concrete fix.** Rewrite Option C as: "screened Coulomb in a bulk-water dielectric: `eps_r = 78.54`, `exp(-kappa d)` with `kappa = 1/7.85 A^-1` at 150 mM. This is not a dielectric-boundary calculation and does not claim to be; it is the cheapest kernel that gets the sign of the surface potential right." Add a sentence noting the current `eps_r = 4` is what produces the documented 20x magnitude error and that fixing it is free.

---

## F5. Answering the ion-exclusion question directly: the `(1 + kappa a)` term is real, and it provably cannot matter here. [measured]

**What the plan says.** Nothing. It specifies `exp(-d / lambda)` without qualification.

**What is actually true.** The correct linearised Poisson-Boltzmann potential of a point charge `q` at the centre of an ion-free sphere of radius `a` is

```
phi(r) = q * exp(-kappa (r - a)) / (4 pi eps_0 eps_r (1 + kappa a) r)
```

so relative to `exp(-kappa r)` the correction is a factor `exp(kappa a) / (1 + kappa a)`. That factor **does not depend on r**. It is a single positive constant multiplying the entire sum, identical for every point and every protein.

Measured, with `a = 2.0 A` (the `ION_RADIUS` used in the APBS runs) and lambda = 7.85 A:

```
max |delta| in Jaccard         : 0.000e+00
max |delta| in Spearman        : 0.000e+00
max |delta| in fraction positive: 0.000e+00
max |delta| in NSurfPos        : 0.000e+00
mean(v_with_a) / mean(v_without_a) = 1.0282 for all 13 proteins
analytic exp(kappa a)/(1 + kappa a) = 1.0282
```

So the correct expression differs from the plan's by a uniform 2.8 per cent scale at `a = 2 A` (10.3 per cent at `a = 4 A`, which is closer to physical for ion-to-atom-centre closest approach). It changes absolute magnitudes and nothing else: not a rank, not a sign, not a count, not a mean-over-a-thresholded-set boundary.

**Concrete fix.** State the full expression in Option C for correctness, note in one sentence that the prefactor is distance-independent and therefore cannot affect any thresholded or rank-based feature, and pick a value of `a` explicitly rather than leaving it at zero by omission. This is a documentation fix, not a physics risk.

---

## F6. The memory model is wrong by 3.2x in the safe direction and 1.5x in the unsafe direction at the same time, and the named largest structure is not the largest. [measured]

**What the plan says.** Section 4: "The largest structure needs a 385x321x321 grid, roughly 28 GB of resident memory by the rough estimate of about 20 grid-sized arrays. Eight of those at once would exhaust the machine."

**What is actually true.** Three separate errors.

*First, 28 GB is 3.2x too high for that grid.* Measured peak RSS with `/usr/bin/time -v` on four real APBS `mg-auto` runs at this repo's exact `elec` settings, spanning 2.1 to 39.7 million grid points:

| dime | grid points | peak RSS | bytes / grid point |
| --- | --- | --- | --- |
| 129 x 129 x 129 | 2,146,689 | 601 MiB | 293.7 |
| 193 x 193 x 161 | 5,997,569 | 1,492 MiB | 260.8 |
| 225 x 225 x 257 | 13,010,625 | 3,172 MiB | 255.6 |
| **385 x 321 x 321** | **39,670,785** | **8,829 MiB (8.6 GiB)** | **239.0** |

Least squares over the four: `peak_RSS_MiB = 224.4 * (grid points in millions) + 164.6`, asymptotically 235 bytes per grid point. The plan's own named grid measures **8.6 GiB, not 28 GB**.

*Second, the stated rule of thumb does not produce the stated number, and is wrong in the crash direction.* "About 20 grid-sized arrays" is 20 x 8 bytes = 160 bytes per grid point, which for 385x321x321 gives 5.9 GiB, not 28 GB. The 28 GB figure is only reachable from a grid of about 173 million points, which is presumably the "slowest 928 s" structure rather than the one named. Whichever it is, 160 bytes per point is **32 per cent below measurement**, so anyone who applies the plan's rule of thumb to a new structure rather than copying its number will under-provision. `pdb2pqr`'s own `psize.py` uses `GMEMFAC = 200` bytes per grid point, which is much closer to the measured 235.

*Third, 385x321x321 is not the largest grid in the set.* It is the largest of the six structures they ran. I replicated `psize`'s sizing rules exactly (`CFAC = 1.7`, `FADD = 20.0`, `SPACE = 0.50`, and the `32k+1` rounding in `set_fine_grid_points`) and reproduced all four rows of the plan's own table to the grid point, including `ARH96666` at 385x321x321. Applied to all 834:

```
Mpts: median 4.17   p95 12.84   max 106.85
largest: ARH96700  577 x 385 x 481 = 106.9 M points  -> predicted peak 23.6 GiB
next:    ARH96666  385 x 321 x 321 =  39.7 M points  ->  measured peak  8.6 GiB
structures needing > 8 GiB: 4      > 4 GiB: 25      > 2 GiB: 148
8 processes at the median grid: 8.2 GiB total
```

The real conclusion is the opposite of the plan's: memory is **not** the binding limit. Eight parallel processes at the median grid need 8.2 GiB of the 113 GiB free. Only four structures in the whole set exceed 8 GiB, and the single worst needs 24 GiB. A static "run the four biggest alone, everything else 8-wide" rule is enough; the elaborate memory-aware scheduler the plan wants to build is not needed. See F7, which removes even that.

**Evidence.**

```bash
# in a scratch dir, with the pqr and .in copied out of apbs_results so nothing there is touched
sed 's/^    dime .*/    dime 385 321 321/' 6PO0_pdie4.in > 6PO0_big.in
/usr/bin/time -v apbs 6PO0_big.in
# Elapsed (wall clock) time: 4:15.19
# Maximum resident set size (kbytes): 9260996
```

**Concrete fix.** Replace the memory paragraph with the measured model `peak_RSS_MiB = 224 * Mpts + 165` (four runs, 2.1 to 39.7 M points, `/usr/bin/time -v`), name `ARH96700` (577x385x481, 23.6 GiB predicted) as the largest, and delete "memory, not CPU, is the binding limit" and the memory-aware scheduler from the deliverables. Say instead: four structures run alone, the rest 8-wide.

---

## F7. The expensive tail of the compute run is AlphaFold disorder, not protein size. Trimming it removes 22 per cent of the runtime and the entire memory problem. [measured]

**What the plan says.** Section 4 treats grid volume as a property of protein size and proposes scheduling around it. Section 3's "Composition" argument worries about the *Neijenhuis* set being unrepresentative and says nothing about the AlphaFold models.

**What is actually true.** Grid volume in this set is driven by low-confidence extended tails, not by mass. Per-residue pLDDT from the B-factor column of the 834 models, joined to the replicated grid sizing:

```
median mean-pLDDT across the 834          : 95.0
median fraction of residues with pLDDT<70 : 0.017
same, among the 25 largest grids          : 0.160   (9x)
Spearman(fraction pLDDT<70, grid points)  : 0.354
```

Individual cases:

| id | atoms | residues | mean pLDDT | frac pLDDT<70 | bbox | grid | Mpts |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ARH96700 | 11,962 | 1,486 | 81.5 | 0.11 | 259 A | 577x385x481 | 106.9 |
| **ARH96666** (the plan's "largest") | 10,348 | 1,329 | **66.7** | **0.48** | 168 A | 385x321x321 | 39.7 |
| ARH96687 | 4,305 | 557 | 69.9 | 0.48 | 163 A | 353x257x385 | 34.9 |
| ARH99399 | **1,211** | **156** | 90.5 | 0.06 | **182 A** | 417x161x321 | 21.6 |

`ARH96666`, the structure the plan's whole memory paragraph is built around, has **48 per cent of its residues below pLDDT 70**. `ARH99399` is a 156-residue protein that needs a 21.6 million point grid because AlphaFold laid it out as a string. Solving Poisson-Boltzmann to 0.5 A on the electrostatics of a coordinate set that AlphaFold itself says it does not believe is compute spent on an artefact, and the surface points Prodes generates on those tails are equally meaningless.

Dropping structures by low-confidence fraction:

| cut | structures dropped | runtime saved | largest remaining grid | peak memory |
| --- | --- | --- | --- | --- |
| frac pLDDT<70 > 0.5 | 2 (0.2 %) | 0.4 % | 106.9 Mpts | 23.6 GiB |
| frac pLDDT<70 > 0.3 | 13 (1.6 %) | 4.3 % | 106.9 Mpts | 23.6 GiB |
| frac pLDDT<70 > 0.2 | 29 (3.5 %) | 7.5 % | 106.9 Mpts | 23.6 GiB |
| **frac pLDDT<70 > 0.1** | **106 (12.7 %)** | **21.6 %** | **23.8 Mpts** | **5.4 GiB** |

The last row eliminates the memory problem entirely: nothing needs more than 5.4 GiB, so a flat 8-wide schedule fits in 44 GiB and the scheduler in the deliverables can be deleted.

The counterpoint, in the plan's favour and worth putting in the document: the disorder concern is *mild overall*. Median mean pLDDT is 95.0 and only 3.5 per cent of structures have more than a fifth of their residues below 70. It is concentrated almost entirely in the expensive tail, which is what makes trimming so cheap.

**Concrete fix.** Add a pre-filter step to section 4 and to `scripts/apbs_comparison.py`: read pLDDT from the B-factor column, and either drop structures with more than 10 per cent of residues below pLDDT 70 or trim those residues before computing the surface. Report the count dropped. Then state the schedule as flat 8-wide and remove the memory-aware scheduler from the deliverables.

---

## F8. The coverage argument, the strongest reason to run 834 structures, is only 8 per cent satisfied. [measured]

**What the plan says.** Section 3, "Coverage": "All 13 are net negative, from -5 to -35 e. The proteins the rule most needs to work on are the ones with positive or near-neutral net charge... the calibration set contains none of them." Section 4 then proposes the 834 as the remedy.

**What is actually true.** The 834 are also overwhelmingly acidic. Net charge at pH 7 computed by Henderson-Hasselbalch from the PROPKA values in the `_pka.json` files each model already carries, which is the same charge model Prodes uses:

```
n = 834, missing pKa files = 0            [verified: the plan's "each with a matching _pka.json"]
netQ:  min -49.4   p5 -21.8   median -6.5   p95 +1.8   max +34.2
net positive (netQ > 0): 67 of 834  =  8.0 %
netQ > -5              : 327 of 834 = 39.2 %
```

Going from 0 to 67 net-positive proteins is a real improvement and the run is worth doing for it. But after a homology group split (F3) the held-out net-positive population is on the order of 13 proteins, which is exactly the n the plan spends section 3 explaining is not a calibration. A random or size-stratified split will produce a validation that says nothing about the regime the fix exists for.

**Concrete fix.** Change section 6 item 3 from a post-hoc check to a design constraint: stratify the split on the sign of net charge so that a fixed number of net-positive clusters is held out, and report the held-out metrics separately for `Q > 0`, `-5 < Q <= 0` and `Q <= -5`. State the n in each stratum up front in section 3, because 67 net-positive proteins is the real limit on what this run can conclude. If the goal is a rule that works on cation exchange, consider deliberately sourcing basic proteins outside the Disela set rather than relying on the 8 per cent.

---

## F9. It is 10 features, not 3, and the "negative side is fine" asymmetry is a property of the test set, not of Prodes. [measured]

**What the plan says.** Section 1: "Three features in the default 54-feature set are built on that population, `NSurfPosEpFormal`, `SurfEpPosFormalMean` and `SurfEpPosFormalStd`." And the asymmetry table: "Anion exchange is driven by negative surface, which Prodes describes well."

**What is actually true, part one: the count.** Ten of the default 54 are sign-thresholded electrostatic populations, from `src/prodes/data/features_reduced.yaml`:

```
NSurfPosEpFormal  SurfEpPosMeanFormal  SurfEpPosStdFormal
                  SurfEpNegMeanFormal  SurfEpNegStdFormal
NShellPosEpFormal ShellEpPosMeanFormal ShellEpPosStdFormal
                  ShellEpNegMeanFormal ShellEpNegStdFormal
```

`calculate_shell_features` builds `pos_potentials` and `neg_potentials` with the same `> 0` / `< 0` split, over potentials from the same unscreened Coulomb kernel via `charged_atom_arrays` and `process_shell_plane`, so the shell trio carries the identical monopole offset. A further 27 in the full 105 are affected, including the whole `SurfEpPos*Average` and `SurfEpNeg*Average` suites from `calculate_average_chargesurface_grid_features`. Option A therefore needs its threshold applied in **three** functions with **three** separately fitted constant pairs, since the shell potential and the average-charge potential are different kernels with different scales. The plan's "smallest possible change" and "leaves every other feature bit-identical" are both wrong as written. (Option C, applied in `Property_point.set_ep`, fixes the surface and average-charge families in one place and needs a matching factor in `process_shell_plane`.)

**What is actually true, part two, and this is the part that will mislead a biotechnologist.** The clean asymmetry in section 1 exists because all 13 test proteins are acidic. For a basic protein the failure mirrors exactly, and the repo already contains the demonstration. `tests/data/1GDW.pdb.zip`, formal charge +7, pI 8.9:

```
NSurfPoints        7039
NSurfPosEpFormal   7020      (99.7 % of the surface)
NSurfNegEpFormal     18
SurfEpMinFormal    -0.27     SurfEpMaxFormal 2.86
```

and with the PROPKA pKa file from `tests/data/hpp_1GDW.pkout` applied, the negative population becomes **empty** (see F13, where this is the single empty-array call site in the whole test suite). So on a net-positive protein, Prodes calls essentially the entire surface positive and the anion-exchange-relevant features are the ones that collapse. A reader who takes away "the negative features are safe, only cation exchange is affected" will apply that to exactly the basic proteins where it is false.

**Concrete fix.** Correct the count in section 1 to ten of 54 and name the shell features. Add a paragraph: "the asymmetry is a property of this test set, not of the method. All 13 proteins are net negative. On a net-positive protein the failure mirrors: `tests/data/1GDW` has formal charge +7 and Prodes calls 7,020 of 7,039 surface points positive, so the negative-surface features collapse instead." Add 1GDW to section 6 item 4 as the mirror test case.

---

## F10. The section 2 physics is arguing for a true conclusion with a proxy that happens to fit better than the correct quantity. [measured]

**What the plan says.** "The offset between Prodes and APBS is the unscreened monopole term. The potential at the surface of a sphere of radius R carrying net charge Q is proportional to `Q/R`... Since the number of surface points scales with area, `netQ / sqrt(n_points)` is that term." Followed by "This is confirmed rather than assumed. The Prodes mean surface potential is 99.1 % explained by `netQ / sqrt(n_surface_points)`, which is exactly `Q/R`."

**First, the numbers all reproduce.** [verified] Recomputed independently from the structures and the per-point files:

| plan claim | my value |
| --- | --- |
| mean surface EP vs `Q/sqrt(n)`, Prodes R2 0.991 | 0.9907 |
| same, APBS R2 0.683 | 0.6827 |
| threshold range -0.40 to -3.68 V | -0.396 to -3.677 V |
| threshold R2 vs `Q` 0.923 | 0.923 |
| threshold R2 vs `Q/sqrt(n)` 0.940 | 0.940 |
| threshold R2 vs mean EP 0.954 | 0.954 |
| residual of the fitted rule 0.21 V | 0.214 V |
| median Jaccard 0.058 naive, 0.601 at the APBS zero | 0.058, 0.601 |
| `NSurfNegEp` Spearman 0.956 | 0.956 |

The conclusion of section 2 is right. The argument for it has three problems.

**Problem one: `sqrt(n_surface_points)` is a mediocre proxy for radius, and the plan asserts it without checking.** Measured against four different radius definitions over the 13:

| radius definition | Pearson with `sqrt(n)` | R2 through the origin | spread of the ratio `sqrt(n)/R` |
| --- | --- | --- | --- |
| radius of gyration | 0.976 | 0.943 | 5.14 to 6.75 (**31 %**) |
| mean surface point distance from centroid | 0.962 | 0.843 | 50 % |
| harmonic mean surface radius | 0.882 | 0.717 | 63 % |
| max atom distance from centroid | 0.595 | 0.144 | 181 % |
| longest bounding box axis | 0.872 | 0.735 | 77 % |

The best case is a 31 per cent spread in the constant of proportionality. That is not "exactly `Q/R`". It is good enough to carry a qualitative argument and not good enough to be described as an identity.

**Problem two: the physically correct radius fits worse than the proxy, which is a curve-fitting signature.** The quantity that actually enters a mean-of-`1/d` is the harmonic mean surface radius. Regressing the Prodes mean surface EP on `netQ / R` for each definition:

```
netQ / sqrt(n_points)         R2 = 0.9907
netQ / radius_of_gyration     R2 = 0.9873
netQ / mean_surface_radius    R2 = 0.9731
netQ / harmonic_mean_radius   R2 = 0.9442
netQ  alone                   R2 = 0.7950
```

When the ad hoc proxy beats the physically motivated quantity by 0.05 R2 on n = 13, the 0.991 is telling you about the fit, not about the physics. All the divisors do the same job, which is to remove a size trend; `sqrt(n)` wins because it is the smoothest of them.

**Problem three, and this is the fix.** There is a clean, exact statement available and the plan does not use it. By the mean value theorem for harmonic functions, the average of `1/|r - r_i|` over a *sphere* of radius R enclosing all the charges is exactly `1/R`, for every interior charge, whatever the distribution. So the mean surface potential of a sphere is *exactly* `Q / (4 pi eps_0 eps_r R)`, independent of where the charge sits inside. That is the real reason the mean surface EP tracks `Q/R` so tightly, and it is a much stronger argument than a regression on 13 points. It also tells you precisely where the 31 per cent spread comes from: proteins are not spheres.

**Problem four: "shifts the entire surface of that protein by a single amount" is not what the data show.** The Prodes-to-APBS map is affine with a *slope* that varies by a factor of two across the set. From `apbs_results/apbs_vs_prodes_summary.csv`, `slope_volts_per_volt` runs 0.0348 (1TRH) to 0.0719 (4F5S). The threshold is `-intercept / slope`, so it mixes an additive monopole with a per-protein scale change. Calling it "an offset" is a simplification that a non-specialist reader will take literally, and it is the reason Option A needs a fitted slope of 0.743 rather than 1.0.

**Concrete fix.** Replace the "since the number of surface points scales with area, `sqrt(n)` is R" sentence with the mean value theorem statement, then add: "proteins are not spheres, and `sqrt(n_surface_points)` tracks the radius of gyration only to about 31 per cent, so this is an argument for the form of the dependence, not a calibration." Replace "shifts the entire surface by a single amount" with "shifts and rescales", and cite the measured slope range 0.035 to 0.072.

---

## F11. The 0.21 V residual is presented as small. In feature units it is 32 to 228 per cent. [measured]

**What the plan says.** Section 2 and the results document: "A single linear rule... has a residual of 0.21 V over a threshold range of 3.3 V." The framing invites the reader to divide and conclude 6 per cent.

**What is actually true.** The relevant question is what a 0.21 V threshold error does to the feature, and the answer is that the surface EP distribution is dense right where the threshold sits. Measured swing in `NSurfPosEpFormal` between thresholds of `t - 0.21` and `t + 0.21` V, as a percentage of the count at `t`:

| protein | N at t | N at t-0.21 | N at t+0.21 | swing |
| --- | --- | --- | --- | --- |
| 1TRH | 6,658 | 7,711 | 5,586 | 32 % |
| 1BSQ | 2,795 | 3,366 | 2,391 | 35 % |
| 1AO6 | 28,439 | 34,968 | 22,583 | 44 % |
| 6PO0 | 25,312 | 32,931 | 17,754 | 60 % |
| 1CF3 | 4,983 | 6,500 | 3,062 | 69 % |
| 4PEP | 136 | 324 | 71 | **186 %** |
| 6FRV | 1,380 | 3,526 | 383 | **228 %** |

The worst cases are exactly the proteins the fix exists for: the ones with very little positive surface, where the threshold sits far out in the tail and the count is most sensitive. Median swing across the 13 is 60 per cent.

**Concrete fix.** Report the residual in feature units alongside volts: "a one sigma residual of 0.21 V changes `NSurfPosEpFormal` by a median of 60 per cent, and by 186 and 228 per cent on 4PEP and 6FRV, the two proteins with the least positive surface." This also strengthens the case for C, whose error is not concentrated at the decision boundary in the same way.

---

## F12. Answering the implementation question: no two-pass is needed, the parallel model is not in the way, and the change breaks a documented invariant. [measured, by reading the executed code path]

**What was asked.** Whether a per-protein threshold that depends on the mean over all points forces a two-pass structure, and whether `prodes/parallel.py` and `process_surface_grid_cell` make it awkward.

**What is actually true.** It does not, and they do not. `calculate_surface_grid_features` already materialises the complete potential array **in the parent process**, after `run_tasks` returns and before the positive population is built:

```python
cell_values = run_tasks(process_surface_grid_cell, len(cells), "surface grid features")
...
eps = np.array([point.ep for point in surface_points])     # every point, in the parent
...
positive_eps = np.array([ep for ep in eps if ep > 0])
```

`process_surface_grid_cell` returns `(ep, lipo)` tuples that the parent writes back with `Property_point.set_values`, so the workers never need to know the threshold. `SURFACE_GRID_STATE` needs no new key. The change is three lines in one function, single pass, zero effect on the fork-based worker model:

```python
eps = np.array([point.ep for point in surface_points])

# Per-protein positive-surface cutoff. The unscreened Coulomb sum carries a
# monopole term that shifts the whole surface, so `> 0` selects the wrong
# population; the constants are fitted, see docs/visualisation/....
threshold = EP_POSITIVE_THRESHOLD_SLOPE * eps.mean() + EP_POSITIVE_THRESHOLD_INTERCEPT
features["SurfEpPosThresholdFormal"] = round(threshold, 3)

positive_eps = eps[eps > threshold]
negative_eps = eps[eps <= threshold]        # NOT eps < 0, see below
```

**The trap the plan does not mention.** Keeping `negative_eps = eps[eps < 0]` while moving the positive cutoff to a negative threshold makes the two populations **overlap** on `threshold < ep < 0`. Every point in that band is counted in both. That silently breaks the invariant the code's own comment relies on:

```
# NSurfNegEpFormal = NSurfPoints - NSurfPosEpFormal, near-exact: the two counts
# partition the same point set, apart from points whose potential is exactly zero.
```

and it invalidates the redundancy argument that comment is used to justify. Either split at the same threshold on both sides (which changes `SurfEpNegMeanFormal` and `SurfEpNegStdFormal`, so the plan's "leaves every other feature bit-identical" becomes false) or document the overlap explicitly.

**Other implementation notes.**

- Option A must be repeated in `calculate_average_chargesurface_grid_features` and `calculate_shell_features` with their own fitted constants, or the `*Average` and `Shell*` positive features stay broken (F9).
- The threshold is a function of `eps.mean()`, which is `SurfEpMeanFormal` before rounding. Emit the threshold as its own feature so the fitted constant is auditable in the output, not just in the source.
- Option C is genuinely a one-line change in `Property_point.set_ep`, and lambda travels to the workers by the route `ph` already travels. One thing to check: `set_ep` ends with `round(float(ep), 2)`. Screening shrinks the spread, so the quantisation gets coarser. Measured, it stays acceptable: 258 to 422 distinct values per protein after rounding (against 321 to 586 unscreened), and 0.1 to 4.5 per cent of points land on exactly 0.00 and would be excluded from both populations by the strict `> 0` / `< 0` tests. Move to `round(..., 3)` when implementing C.
- Option C also permits a real distance cutoff for the first time. `set_ep` currently defaults to `cutoff=10000`, i.e. every charge for every point. With screening, 5 lambda (about 40 A) drops terms below `exp(-5)` and would cut the inner loop substantially on large structures. Worth listing in Option C's pros.

---

## F13. Answering the NaN question: no tests break, exactly one test touches the branch, and it would then pass for the wrong reason. Scope is 39 features, not 3. [measured]

**What the plan asks.** Open question 3: "Is `standard_features` returning 0 for an empty array the right behaviour at all, or should it be NaN?"

**What is actually true.**

*Nothing breaks.* I patched `prodes.run.standard_features` to return NaN for empty input via a pytest plugin loaded before collection, and ran the full suite in the `prodes` env:

```bash
PYTHONPATH=<scratch> python -m pytest -q --no-header -p nanpatch
# 210 passed in 32.09s     (baseline, unpatched: 210 passed in 33.62s)
```

*But the suite has essentially no coverage of the branch.* Instrumenting the patch to count calls:

```
NANPATCH: standard_features calls=120  empty-array calls=1  names=['SurfEpNeg']
```

One empty-array call in the entire suite. The stack traces to `tests/test_pka.py::test_a_pka_file_changes_the_calculated_features`, on `tests/data/1GDW.pdb.zip`, the net-positive protein from F9: with the PROPKA pKa file applied it has no negative surface points at all.

*And under NaN that test would pass for the wrong reason.* It counts changed features with

```python
changed = [column for column in numeric if default[column].iloc[0] != predicted[column].iloc[0]]
...
assert len(changed) > 10
```

`nan != nan` is True, so a feature that is censored in both runs would be counted as having changed. The `len(changed) > 10` assertion becomes silently easier to satisfy. It happens not to matter today because only one of the two runs is censored, but the test would stop measuring what it claims to.

*The scope is much wider than the plan implies.* Features produced by `standard_features` over a population that can be empty: **12 of the default 54** (`SurfEpPosMean/StdFormal`, `SurfEpNegMean/StdFormal`, `SurfPosMhpMean/Std`, `SurfNegMhpMean/Std`, `ShellEpPosMean/StdFormal`, `ShellEpNegMean/StdFormal`) and **27 of the 51 full-only** ones, so **39 of 105**. Changing the sentinel is a change to a third of the feature set, not to three features.

*And a blanket NaN would be wrong for one family.* The sum of an empty set is 0. `SurfEpPosSumFormal = 0` for a protein with no positive surface is the correct answer, not a censored one, and turning it into NaN would destroy real information. Same for the `N*` counts, which are computed by `len()` outside `standard_features` and correctly stay 0.

**Concrete fix.** Answer open question 3 as: NaN for `Mean`, `Trimean`, `Median` and `Std`; keep 0 for `Sum`; the `N*` counts are already right. Note that it touches 39 of 105 features. Add a real regression test on a net-positive structure, using the 1GDW fixture already in `tests/data`, that pins the empty-population behaviour on both signs, since the suite currently exercises it exactly once and incidentally. Fix `test_a_pka_file_changes_the_calculated_features` to compare with a NaN-aware equality before the sentinel changes under it.

---

## F14. The runtime model is sound for the 834 and does not transfer to the 13, and section 4 is internally inconsistent about which structure is the biggest. [measured]

**What the plan says.** "`seconds = 5.35 * (grid points in millions) + 0.9`, R2 = 0.999 on the six." Median 28 s, p95 111 s, slowest 928 s, total 9.8 h.

**What is actually true.**

*The model is defensible for AlphaFold monomers.* The plan's own four tabulated runs give 4.19, 4.40, 4.50 and 4.44 APBS seconds per million grid points across an 18x size range, which is flat to within 7 per cent. And I can confirm the grid sizes are right: replicating `psize`'s `CFAC/FADD/SPACE` rules and the `32k+1` rounding in `set_fine_grid_points` reproduces **all four rows exactly**, including 129x129x129, 257x161x225, 225x193x225 and 385x321x321. That is a strong independent check on the whole extrapolation.

*It does not transfer to the Neijenhuis 13.* Measured on this machine at the repo's own `elec` settings, wall clock: 1BSQ 10.05 s at 2.15 Mpts, 1TRH 27.71 s at 6.00 Mpts, 6PO0 94.67 s at 13.01 Mpts. A linear fit through the two smaller ones predicts 59.8 s for 6PO0; it took 94.7 s, 58 per cent over. `npbe` is a nonlinear solve and its iteration count grows with the magnitude of the potential, so multi-chain, highly charged structures are superlinear in a way monomers are not. If the 13 are re-run in the same batch, budget for that; it does not affect the 834.

*Section 4 contradicts itself about the biggest structure.* "Slowest 928 s" implies `(928 - 0.9) / 5.35 = 173` million grid points. "The largest structure needs a 385x321x321 grid" is 39.7 million. Both cannot be true. My reproduction of `psize` over all 834 gives a maximum of 106.9 million (`ARH96700`, 577x385x481) and a total of 7.3 h under the plan's own formula, against the plan's 9.8 h. The 28 GB memory figure is consistent with the 173 Mpts number and inconsistent with the 385x321x321 one, which is how the error in F6 arose.

*Two smaller checks that came back clean, worth stating so nobody re-does them.* Every structure has at least 5 A of fine-grid padding on every axis, so the `cfac` cap in `set_fine_grid_dims` never bites and the `"surface points fell outside the fine grid"` guard in `compare_structure` should not fire on any of the 834 [measured, 0 of 834 at risk]. And the plan's size statistics reproduce exactly: median 2,562 atoms, p95 6,081, max 11,962; 12 structures with a bounding box over 150 A; 834 structures with 0 missing `_pka.json` [verified].

**Concrete fix.** State which structure is the slowest and which is the largest, with its grid, and reconcile 928 s against 385x321x321. Add one sentence: "the linear model is calibrated on AlphaFold monomers and does not transfer to the multi-chain Neijenhuis structures, where 6PO0 at 13.0 M points took 94.7 s against a 59.8 s linear prediction." Show all six calibration runs or say six and show six; the table has four.

---

## F15. Section 1's "harmless up to an affine transform" is a within-protein result being read as a between-protein one, and it quotes the best of 13. [measured]

**What the plan says.** Section 1: "on beta-lactoglobulin an affine map of the Prodes potential reproduces 73 % of the APBS variance, and after rescaling the colour ramp the two point clouds are visually near identical, patch for patch. So for any monotone feature (mean, min, max, standard deviation) the discrepancy is largely harmless."

**What is actually true.** Two problems.

*1BSQ is the best of the 13, not a typical one.* From `apbs_results/apbs_vs_prodes_summary.csv`, `pearson_P_A` runs 0.468 to 0.857, and 0.857 (R2 = 0.734, the quoted 73 per cent) is the **maximum**. The median is 0.733, R2 = 0.537. Section 1 opens by quoting its single best case as the reassurance.

*The "monotone features are harmless" inference does not survive the move to between-protein.* A monotone feature is invariant to an affine shift *within one protein*. The models operate across proteins, where each protein has its own shift and its own slope (F10). Measured cross-protein Spearman against the APBS equivalent, n = 13:

| feature | Prodes today |
| --- | --- |
| `SurfEpMean` | 0.582 |
| `SurfEpStd` | 0.060 |
| `SurfEpMax` | 0.346 |
| `SurfEpMin` | -0.236 |

`SurfEpStd` at 0.060 and `SurfEpMin` at -0.236 are not "largely harmless"; they are noise and anti-correlation respectively, in features the plan explicitly names as safe. All four improve under Option C (0.819, 0.516, 0.467, -0.071), which is a further argument for F1.

**Concrete fix.** Change "on beta-lactoglobulin" to "on beta-lactoglobulin, the best of the 13; the median is 54 per cent." Replace "for any monotone feature the discrepancy is largely harmless" with the measured cross-protein table above, and say plainly that the affine invariance holds within a protein and does not hold across proteins, which is the level a retention model works at. This is the single change most likely to prevent a biotechnologist reader from drawing the wrong conclusion from this document.

---

## Answers to the plan's own open questions

**Is option C the right primary target, making A a stopgap?** Yes, and stronger than that. C beats A, beats B, and beats the oracle threshold on held-in data, with zero fitted constants and no sensitivity to lambda across 5 to 12 A. See F1. Fix the dielectric incoherence when you write it up (F4) and state the correct DH form (F5).

**Should the threshold rule use mean EP, `Q/sqrt(n)`, or `Q` directly?** Moot if C is adopted. If A is kept: `Q/sqrt(n)` (R2 0.940) rather than mean EP (0.954). The 0.014 difference is far inside the n = 13 uncertainty the plan itself invokes in section 3, and `Q/sqrt(n)` is not derived from the quantity being corrected, so it removes the circularity. The circularity is real: mean EP is 99.1 per cent predicted by `Q/sqrt(n)` (F10), so the two predictors are the same predictor and the fit is choosing between them on noise.

**Is 0 the right sentinel for an empty array?** No for `Mean`/`Trimean`/`Median`/`Std`, yes for `Sum`. It touches 39 of 105 features, breaks nothing in the current suite, and the one test that touches the branch would then pass for the wrong reason. See F13.

**Does changing three released feature values require a major version bump?** It is not three. Under A it is at least ten in the default set and 30-plus in the full set (F9); under C it is all 105, since every EP-derived feature changes. Either way this is a breaking change to the output contract of a package whose output feeds trained models. Major bump, and ship the flag default-off for one minor release so both feature sets can be produced from one install for the model refit.

---

## What I would do next, in order

1. Settle F3 (build a real homology clustering) and F7 (pLDDT pre-filter) **before** launching the 847-structure run. Both change what the run produces; neither is expensive.
2. Run the 834 with the flat 8-wide schedule and the four large structures serialised (F6). Budget one overnight, not a scheduler.
3. Score C, B and A on held-out clusters, stratified by the sign of net charge (F8), with the size baseline reported alongside every number (F2).
4. Only then decide the sentinel question (F13) and the version policy, since C changes every feature anyway.
