# Review 1: `plan_fix_positive_ep_features.md`

**Reviewer lens:** statistics, calibration, experimental validity. **Date:** 2026-08-21. **Branch:** `visualisation`.

Everything below was re-computed from `apbs_results/*/*_points.csv.gz`, `apbs_results/*/*_residue_charges.csv`, `docs/visualisation/apbs_vs_prodes_summary.csv`, `../biochai/data/04b_compute_surface_pka/`, `../biochai/data/04a_alphafolddb_structures/`, `../biochai/data/11_protein_clusters/` and `../biochai/data/13_combined_dataset/combined_dataset.csv`, using `spearman`, `pearson` and `average_ranks` from `scripts/apbs_comparison.py`. No file other than this one was touched.

---

## Verdict

**Every headline number in the plan reproduces exactly. [VERIFIED]** The arithmetic is clean and the document is unusually honest about its own limits. The problem is not the numbers, it is what they are being taken to mean.

Three findings are, in my judgement, blocking:

1. The `0.641 -> 0.973` improvement is **not evidence that the threshold rule works**. Simply counting surface points, with no electrostatics at all, scores 0.929 on the same target. The fitted rule beats that trivial baseline by +0.056 with a bootstrap 95 % CI of [-0.106, +0.249] and P(no improvement) = 0.22. On the size-free version of the same quantity (positive *fraction* rather than positive *count*) the recalibrated rule scores **0.769 against the naive rule's 0.780**, i.e. it is very slightly worse. The whole apparent repair is the count feature recovering surface area.
2. The retention data this is all meant to serve is **anion exchange at pH 7 to 10**. There is no cation exchange arm anywhere in `13_combined_dataset`, so section 6 item 6, the "only test that measures what anyone actually cares about", cannot be run as written. The plan's central motivating asymmetry points at a chromatographic mode the project has no data for.
3. `../biochai/data/11_protein_clusters` does **not** contain a clustering that can support the group split in section 6 item 1. It contains an ESM cosine-similarity matrix and a duplicate-detection result with **2 clusters covering 4 proteins**. The split the plan says "already exists" has to be built from scratch.

Two further findings materially change the plan: the fix touches **six** features in the 54-set, not three, and the second constant in the rule is an artefact of one high-leverage protein.

The overnight APBS run is still worth doing. Almost nothing else in the plan should be executed in its current form.

---

## F1 (blocking). The headline repair is a surface-area effect, not an electrostatics repair

**Plan says.** Section 2: `NSurfPosEp` Spearman against APBS goes from 0.641 naive to 0.995 at the true threshold and 0.973 at the predicted threshold. Section 5 Option A: "Fixes the count feature nearly completely (0.641 to 0.973 on the pilot)."

**Actually true.** All three values reproduce to the digit [VERIFIED: 0.6407, 0.9945, 0.9725]. But the comparison omits the only baseline that matters. `NSurfPosEp` is a **count**, and counts are dominated by how many surface points the protein has:

| predictor of the APBS positive-point count | Spearman, n=13 |
| --- | --- |
| `n_surface_points` alone, no electrostatics whatsoever | **0.9286** |
| naive `ep > 0` | 0.6407 |
| threshold rule, in sample | 0.9725 |
| threshold rule, leave-one-out | 0.9560 |
| oracle threshold (APBS zero crossing) | 0.9945 |

Paired bootstrap over the 13 proteins, 20 000 resamples:

- rule minus naive: **+0.350**, 95 % CI [+0.093, +0.817], P(diff <= 0) = 0.000. Real.
- rule minus `n_surface_points`: **+0.056**, 95 % CI [-0.106, +0.249], P(diff <= 0) = **0.222**. Not distinguishable from zero.
- LOO rule minus `n_surface_points`: **+0.033**, 95 % CI [-0.147, +0.230], P(diff <= 0) = **0.343**.

Removing size settles it. Comparing the *fraction* of surface each method calls positive:

| | Spearman of predicted positive fraction vs APBS positive fraction |
| --- | --- |
| naive `ep > 0` | **0.780** |
| threshold rule, in sample | 0.769 |
| threshold rule, LOO | 0.758 |
| oracle threshold | 0.885 |

The recalibration buys **nothing** on the size-free quantity. The oracle threshold does (0.885), which proves the information is there and the two-parameter rule fails to extract it. The naive rule's 0.641 is low because four proteins tie at exactly zero; among the nine proteins with a non-zero naive count the naive Spearman is 0.350 while `n_surface_points` alone scores 0.933 on the same nine.

**Evidence.**
```
conda activate prodes
# sys.path.insert(0,'scripts'); from apbs_comparison import spearman
# counts from apbs_results/*/*_points.csv.gz at thresholds 0, t_true, 0.743*meanEP-0.121
```

**Fix.** Report every count claim against `n_surface_points` as the null model, and report the **fraction** as the primary endpoint. Both belong in section 6 item 2 alongside the naive baseline. If the honest statement is "the rule recovers the count ranking no better than counting points does", say that; it changes the recommendation in section 5.

---

## F2 (blocking). There is no cation exchange data, so the motivating asymmetry cannot be tested

**Plan says.** Section 1: the damage is "specific to cation exchange". Section 6 item 6: recompute features from APBS and refit; "if APBS-derived features beat Prodes-derived features on cation exchange and tie on anion exchange, that confirms the mechanism exactly where this analysis predicts."

**Actually true.** `13_combined_dataset/combined_dataset.csv` has 886 rows, target `DRT IEX [-]`, and a `pH` column taking values **7, 8, 9, 10 only**. Disela 2024 contributes 836 rows all at pH 7; Neijenhuis 2025 contributes 13/12/13/12 rows at pH 7/8/9/10. `seq:Charge at measurement pH` is negative for **94.6 %** of rows (median -8.5); only 48 rows are net positive and all 48 are at pH 7. Disela 2024 is an E. coli HCP gradient run (`Disela (2024) data explanation.md`), and Neijenhuis at pH 7 to 10 with proteins at median charge -8 to -58 is anion exchange. Nothing in the dataset is labelled or plausibly cation exchange. [VERIFIED]

So the confirmatory experiment the plan calls the only one that matters has no cation exchange arm. It reduces to "does APBS tie with Prodes on anion exchange", which is not a discriminating test.

**Fix.** Either drop the cation exchange framing and justify the fix on anion exchange terms (see F3, where it survives, on different grounds), or state explicitly that confirming the mechanism requires new CEX retention data that does not exist and is out of scope. Do not write a validation plan whose decisive test cannot be executed.

---

## F3 (major). The "collapse to zero" is monotone-consistent for the actual target, and naively fixing it makes the feature a worse univariate predictor

**Plan says.** Section 1: the three features "collapse to exactly 0 for four of thirteen proteins", framed as pure damage. Section 8 open question 3 asks whether 0 should be NaN so a censored feature is visible.

**Actually true.** In the 880 model-ready production rows, `str:NSurfPosEpFormal` is 0 for **97 rows (11.0 %)** and `str:NShellPosEpFormal` for **114 rows (13.0 %)**. Those censored rows are not random: their median `seq:Charge at measurement pH` is **-26.7** against -7.5 for the rest. The feature's Spearman with `DRT IEX [-]` is **-0.444**, i.e. more positive surface means weaker anion exchange retention, exactly as physics predicts. Because the censored rows are the *most* negative proteins, collapsing them to 0 places them at the correct extreme of a monotone relationship. Zero is, by luck, the right answer for rank-based and tree-based models.

Replacing 0 with NaN, as open question 3 proposes, would therefore **destroy** information for a monotone model, not reveal it. And a proxy for the recalibrated feature (fraction above threshold predicted from `str:SurfEpMeanFormal` using a rule fitted on the 13, R2 0.881, applied to the 880) drops the marginal Spearman with `DRT IEX` from **-0.444 to -0.248**.

The counter-argument, which the plan does not make and should: partialling out net charge,

| feature | partial Spearman with `DRT IEX` given `seq:Charge at measurement pH` |
| --- | --- |
| `str:NSurfPosEpFormal` as shipped | **-0.104** |
| proxy recalibrated feature | **-0.243** |
| total surface points | -0.200 |
| `str:SurfEpMeanFormal` | -0.210 |

So the shipped feature's -0.444 is almost entirely net-charge collinearity, and the recalibrated feature does carry more charge-independent signal. That is the real case for the fix. Note however that "total surface points" alone gets -0.200 of that -0.243, which is F1 again: most of the conditional signal is area, not electrostatics.

**Evidence.** The proxy is a linear extrapolation from 13 proteins and is indicative only, but the marginal correlations and the 11 % / -26.7 censoring numbers are direct measurements. [VERIFIED for the measured parts]

**Fix.** Answer open question 3 with "0 is correct, NaN would be worse, and here is why". Add a section stating that the case for the fix is the **partial** correlation given net charge, not the marginal one, and make that the section 6 endpoint. Any validation that scores the fix by marginal correlation with retention will report a regression.

---

## F4 (major). `11_protein_clusters` cannot support the group split the plan specifies

**Plan says.** Section 6 item 1: "Group split by sequence similarity, not random: `../biochai/data/11_protein_clusters` already exists and should be used, or the calibration will leak between homologues."

**Actually true.** That directory holds `cosine-similarity-matrix.json` (836 x 836 ESM-C embedding cosines), `similarity-clusters.json` containing **2 clusters, each of member_count 2**, `duplicate-candidates.json` with **2 pairs**, and `duplicate_pair_identities.csv` with **2 data rows**. It is a near-duplicate detector, not a clustering of the 834. [VERIFIED]

Worse, the cosine matrix cannot be thresholded into clusters. Off-diagonal ESM cosines over the 836: min -0.208, p5 0.716, **median 0.891**, p95 0.949, max 1.000. **43.3 % of all 349 030 pairs exceed cosine 0.9.** A 0.9 cut would merge nearly half the dataset into one component. ESM cosine at these values is not a sequence-identity proxy.

**Fix.** Build a real clustering with MMseqs2 or CD-HIT at a stated identity threshold (30 % or 40 %) from `03_uniprot_sequence_lookup`, commit the cluster assignment, and cite it. Budget the time; it is minutes of compute but it is a step the plan currently assumes away. Until it exists, section 6 item 1 is not executable.

---

## F5 (major). Six features are affected, not three, and one whole block is outside the calibration's domain

**Plan says.** Section 1 and section 5: "Three features in the default 54-feature set are built on that population, `NSurfPosEpFormal`, `SurfEpPosFormalMean` and `SurfEpPosFormalStd`."

**Actually true.** `src/prodes/data/features_reduced.yaml` (feature_count 54) contains **eighteen** EP features, of which **ten are threshold-defined**:

- positive-threshold: `NSurfPosEpFormal`, `SurfEpPosMeanFormal`, `SurfEpPosStdFormal`, **`NShellPosEpFormal`**, **`ShellEpPosMeanFormal`**, **`ShellEpPosStdFormal`**
- negative-threshold: `SurfEpNegMeanFormal`, `SurfEpNegStdFormal`, `ShellEpNegMeanFormal`, `ShellEpNegStdFormal`

The Shell block is selected by the same `potential > 0` test in `run.py` (the `pos_potentials` line in the shell feature function), on the Sunflower far-field reference points, not the surface grid. Those points sit further from the protein, so the unscreened monopole dominates them *more*, and the calibration constants fitted on surface points do not transfer to them. `str:NShellPosEpFormal` is already 0 for 13.0 % of production rows, a higher censoring rate than the surface feature. The plan never mentions the shell block, and the APBS comparison never sampled the shell points, so there is no ground truth for them at all. [VERIFIED]

The plan's names are also wrong: the code emits `SurfEpPosMeanFormal` and `SurfEpPosStdFormal`, not `SurfEpPosFormalMean` / `SurfEpPosFormalStd`.

**Fix.** Correct the scope statement to six positive-threshold features. Add sampling of the Shell reference points to `apbs_comparison.py` in the same run, or state explicitly that the Shell block is out of scope and remains broken. Fix the feature names before anyone greps for them.

---

## F6 (major). The intercept is a single-protein artefact, and that protein has a known charge-model defect

**Plan says.** Section 2: "The two constants `0.743` and `-0.121` are global: fitted once across the whole calibration set and then frozen."

**Actually true.** The fit reproduces exactly [VERIFIED: slope 0.7434, intercept -0.1207, R2 0.9540, residual SD 0.206 V population / 0.224 V with ddof=2]. But 4PEP sits at meanEP = -5.10 while the next most extreme protein is at -3.26, giving it a **hat value of 0.545** against a mean leverage of 2/13 = 0.154. It is 3.5x the average leverage and dominates the fit.

| fit | slope | intercept | R2 |
| --- | --- | --- | --- |
| all 13 | 0.743 | **-0.121** | 0.954 |
| drop 4PEP | 0.820 | **-0.008** | 0.941 |
| drop 1F6R | 0.739 | -0.124 | 0.950 |

Dropping one protein moves the intercept from -0.121 to essentially zero. The second constant is not identified by this dataset. Compounding it, `apbs_vs_prodes_results.md` itself flags 4PEP as carrying a **phosphoserine that both tools discard**, so the single point setting the intercept is the one with a known missing charge. And 4PEP is simultaneously one of the four zero-positive proteins the fix exists to repair.

Note that a slope of ~0.82 with a zero intercept is close to a scaled mean-centring, which makes Option B and Option A far less distinct than the plan presents them (see F9).

**Fix.** Report leverage and a leave-one-out constant table in `docs/visualisation/ep_threshold_calibration.md`. Do not freeze a two-parameter rule until the intercept survives dropping any single protein. Consider constraining the intercept to 0 on physical grounds and fitting one parameter.

---

## F7 (moderate). Overfitting is small; the plan over-weights it and under-weights everything else. Corrected LOO numbers

**Plan says.** Section 3: "A two parameter fit on 13 points is not a calibration." Section 2: "those constants were fitted and evaluated on the same 13 proteins, with no held-out data."

**Actually true.** The optimism from the two-parameter fit is **small**:

| | in sample | leave-one-out |
| --- | --- | --- |
| Spearman, `NSurfPosEp` count | 0.9725 | **0.9560** |
| Spearman, positive fraction | 0.7692 | 0.7582 |
| R2 for the threshold itself | 0.9540 | 0.9282 |
| RMSE of the threshold | 0.208 V | **0.257 V** |

**LOO Spearman 0.956 against the in-sample 0.973 is the number the plan asked for.** The two-parameter fit is not the problem. F1 is the problem: 0.956 is still statistically indistinguishable from the 0.929 you get by counting points.

Per-protein LOO behaviour matters more than the aggregate, and it is bad exactly where it needs to be good. Ratio of LOO-predicted positive count to the APBS count, for the four proteins that currently report zero:

| | 1CF3 | 1F6R | 4PEP | 6FRV |
| --- | --- | --- | --- | --- |
| LOO count / APBS count | 1.01 | 0.81 | **2.37** | **0.42** |

Across all 13 the LOO ratio spans 0.42 to 2.37 (geometric mean 0.96, 10th to 90th percentile 0.60 to 1.22). So out of sample the "fixed" count is wrong by up to a factor 2.4 on the very cases in section 6 item 4.

**Confidence intervals at n=13** (Fisher z with the Bonett-Wright standard error 1.06/sqrt(n-3)):

| rho | 95 % CI |
| --- | --- |
| 0.641 (naive count) | [0.103, 0.889] |
| 0.929 (`n_surface_points` baseline) | [0.758, 0.980] |
| 0.956 (LOO rule) | [0.846, 0.988] |
| 0.973 (in-sample rule) | [0.903, 0.993] |
| 0.995 (oracle) | [0.982, 0.999] |
| 0.234 (`SurfEpPosMean` naive) | [-0.396, 0.714] |
| 0.123 (`SurfEpPosStd` naive) | [-0.488, 0.653] |

The plan's own statement, "at n = 13 a Spearman of 0.6 has a 95 % interval of roughly 0.05 to 0.87", computes to **[0.036, 0.874]**. [VERIFIED, correctly argued] Note that the CIs for 0.929 and 0.973 overlap almost entirely, which is the same conclusion as the paired bootstrap in F1, and that the `SurfEpPosMean` and `SurfEpPosStd` numbers in the section 1 table are not distinguishable from zero, so the "asymmetry" table should carry intervals.

**Fix.** Replace "the constants were fitted and evaluated on the same 13" with the measured LOO numbers and the per-protein ratio table. Add CIs to the section 1 asymmetry table. Then reframe section 3: 13 is not too few for the *fit*, it is too few for the *comparison against the baseline*, which is a different and more damning problem.

---

## F8 (moderate). The circularity is real but harmless, and irrelevant, because netQ works at least as well

**Plan says.** Section 8: "Mean EP fits best on 13 but is itself derived from the potential being corrected, which is a circularity worth a second opinion."

**Actually true.** It is not vicious circularity. The rule uses a *global* statistic of the field to correct a *local* threshold on the same field, which is ordinary within-sample normalisation, the same logic as subtracting a batch mean. It would be circular only if mean EP depended on the threshold, which it does not. So the plan can stop worrying about it on statistical grounds.

It is nonetheless the wrong choice, because the EP-free alternatives are not worse:

| predictor | R2 with the true threshold | in-sample Spearman, count | **LOO Spearman, count** | LOO RMSE of threshold |
| --- | --- | --- | --- | --- |
| Prodes mean surface EP | **0.9540** | 0.9725 | 0.9560 | **0.257 V** |
| netQ / sqrt(n_points) | 0.9404 | 0.9451 | 0.9286 | 0.296 V |
| netQ | 0.9232 | **0.9835** | **0.9725** | 0.331 V |

netQ alone gives the **best** LOO count Spearman of the three (0.9725 against mean EP's 0.9560), while mean EP gives the best threshold RMSE. The tradeoff is therefore about 0.07 V of threshold accuracy against a slightly better count ranking, and both sit inside the noise established in F1. There is no measurable cost to using netQ, and netQ is independent of the field being corrected, is available before any surface calculation, is one integer per protein, and extrapolates on a physically interpretable axis.

**Fix.** Answer open question 2 with "use netQ". It removes the circularity question, removes a dependence on the surface point count, and costs nothing that this dataset can resolve.

---

## F9 (moderate). Option B is statistically sound and does not throw the signal away, but it fails on the extreme protein

**Plan says.** Section 5 Option B: subtracting the mean "destroys the absolute scale, which is currently the thing carrying the net-charge signal", with a mitigation of keeping the mean as its own feature.

**Actually true, tested directly.** Setting `threshold = mean(ep)` per protein and comparing the selected point set to the APBS positive set:

| scheme | median Jaccard | min Jaccard | Spearman, count | Spearman, fraction | `SurfEpPosMean` Spearman | `SurfEpPosStd` Spearman |
| --- | --- | --- | --- | --- | --- | --- |
| naive `ep > 0` | 0.058 | 0.000 | 0.641 | 0.780 | 0.234 | 0.123 |
| Option A, `0.743*meanEP - 0.121` | 0.586 | 0.247 | 0.973 | 0.769 | 0.341 | **0.126** |
| Option B, mean-centre | **0.618** | **0.058** | 0.896 | **-0.154** | 0.341 | 0.291 |
| oracle, true APBS crossing | 0.601 | 0.221 | 0.995 | 0.885 | 0.335 | 0.407 |

Mean-centring gives the **best median Jaccard of any scheme, better than the oracle threshold**, so it emphatically does not destroy agreement on the thresholded point set. The plan's stated concern is wrong on the point-set question. With the mean retained as its own feature, no information is lost at all, exactly as the plan's own mitigation says.

Two real problems the plan does not identify:

- Mean-centring forces the positive fraction to be near-constant across proteins by construction, so `Spearman(fraction) = -0.154`, destroying the between-protein positive-surface ranking that the oracle threshold recovers at 0.885. If a model uses `NSurfPosEp` as anything other than an area proxy, B is worse than A.
- B's **minimum** Jaccard is 0.058, on 4PEP, where APBS calls only 3.5 % of the surface positive and mean-centring selects roughly 40 %. B fails hardest on the most extreme protein, the same one that dominates A's fit (F6).

Also note the constrained slope-1 model `threshold = meanEP + 0.387` reaches R2 0.840 against the two-parameter 0.954, and dropping 4PEP gives A a slope of 0.820 with a zero intercept. A and B are much closer to the same operation than section 5 implies. They should not be presented as independent options.

**Fix.** Rewrite the Option B "cons" to the two real failure modes above, delete the incorrect "destroys the signal" claim, and state that A and B differ only in one fitted slope so the comparison is a one-parameter question, not a three-way choice.

---

## F10 (moderate). Section 6 is not sufficient. What is missing

Item by item.

- **Item 1, group split.** Cannot be executed (F4). It also does not say what the split ratio is, how many folds, or what the seed policy is. A single held-out split at 20 % gives roughly 170 proteins, on which a Spearman of 0.9 has a 95 % CI of about [0.87, 0.92] by Fisher z, which is workable, but a single split is a lottery. **Use repeated grouped k-fold (5-fold, 5 repeats) and report the spread across folds, not one number.**
- **Item 3, stratify by net charge.** Stratifying is not enough, because subgroup Spearman at the sizes available has no power. Measured on the 834 pKa jsons at pH 7 the net-positive proteins are **52 of 834, 6.2 %** (F11). A 20 % held-out split contains about **10** of them; a Spearman on 10 has a 95 % CI roughly [-0.4, 0.9]. The plan's blocking criterion ("if it degrades for Q > 0 proteins, that is a blocking result") is untestable as specified. **Replace the subgroup comparison with a continuous test: regress the threshold residual on netQ over the full held-out range and test the slope and a quadratic term. That has power; a 10-protein subgroup does not.**
- **Item 5, three options on one split.** Yes, this invites selection bias. Picking the maximum of three correlated estimates on a single split inflates the reported score by roughly 0.5 to 0.9 standard errors of the fold-to-fold spread. **Nest it: choose among A, B and C on inner folds, report only the winner's outer-fold score, and report all three outer scores anyway.**
- **Item 6, end-to-end.** Not executable as written (F2).
- **Missing entirely.**
  - The `n_surface_points` null model (F1). Without it nothing in section 6 is interpretable.
  - The positive **fraction** as an endpoint alongside the count.
  - The Shell block (F5), which has no ground truth at all.
  - A pre-registered decision rule. Section 6 lists measurements but never says which number, at which value, causes the fix to ship or to be abandoned. Write that down before the run, not after.
  - Sensitivity to ionic strength and pH. One condition was run (pH 7, 150 mM); `13_combined_dataset` spans pH 7 to 10 and the papers span buffer conditions. A rule fitted at one condition and applied across four is an untested extrapolation.

---

## F11 (moderate). The coverage argument is right in direction, wrong in magnitude, and the fix is cheaper than the plan thinks

**Plan says.** Section 3: "All 13 are net negative, from -5 to -35 e... the calibration set contains none of them [near-neutral or positive]." Section 4: 834 structures with matching pKa jsons.

**Actually true.** Net formal charge at pH 7 from `apbs_results/*/*_residue_charges.csv`: the 13 span **-5 (1OVT) to -35 (4PEP)**, all negative. [VERIFIED]

834 pdb files in `04a_alphafolddb_structures` and 834 `_pka.json` in `04b_compute_surface_pka`. [VERIFIED] Net formal charge at pH 7 recomputed from the pKa jsons with the Prodes convention (positive types ARG/LYS/HIS/N+ charged when pKa > pH, negative types ASP/GLU/TYR/CYS/C- charged when pKa < pH; validated against the 13, exact on the 7 single-chain structures, off by 1 to 6 on the four multi-chain ones because of the known `convert_propka` chain bug):

| | 834 AlphaFold structures, pH 7 |
| --- | --- |
| min / median / max | -52 / **-7** / +33 |
| net positive, Q > 0 | **52 (6.2 %)** |
| Q >= 0 | 70 (8.4 %) |
| Q > -5, i.e. above the calibration range | **251 (30.1 %)** |
| Q < -35, below the calibration range | 13 |
| inside the calibrated [-35, -5] | 570 (68.3 %) |

So the plan's coverage argument is correct: **30 % of the application set sits outside the calibrated netQ range** and the whole positive tail is unrepresented. But 6.2 % net positive is a thin tail, not the fix the plan implies (F10 item 3).

**The cheap fix the plan misses.** Net charge is a function of pH, and the same jsons give:

| pH | net positive | median netQ |
| --- | --- | --- |
| 4 | 828 (99.3 %) | +19 |
| 5 | 456 (54.7 %) | +1 |
| 6 | 197 (23.6 %) | -3 |
| **7** | **52 (6.2 %)** | **-7** |
| 8 | 31 (3.7 %) | -9 |
| 9 | 29 (3.5 %) | -10 |

**Running a subset at pH 5 as well as pH 7 buys a balanced netQ distribution for a fraction of the cost of the full 834 at one pH.** It is the only way to get genuine Q > 0 coverage, it tests the rule's extrapolation on the axis that matters, and it directly tests the physical claim that the offset is a monopole term. If any single change is made to section 4, it should be this.

---

## F12 (moderate). The memory and largest-structure numbers in section 4 are wrong, in the direction that kills a 10-hour run at hour 8

**Plan says.** Section 4: "The largest structure needs a 385x321x321 grid, roughly 28 GB of resident memory by the rough estimate of about 20 grid-sized arrays."

**Actually true.** Three separate errors.

1. The plan's own arithmetic does not give 28 GB. 385 x 321 x 321 = 39.7 M grid points; 20 arrays x 39.7 M x 8 bytes = **6.3 GB**, not 28.
2. The real number is measurable, and neither figure is right. Fitting APBS high-water memory from the 13 `*_pdie4_apbs.log` files against grid volume gives **high-water MB = 467.5 x (grid points in millions) + 111.6, R2 = 0.9995** [VERIFIED], i.e. about **468 bytes per grid point**. For 39.7 M points that is **18.7 GB**.
3. 385x321x321 is not the largest grid. That is the grid for **ARH96666** (10 348 atoms), which the plan's own table lists. The largest structure in the set is **ARH96700** (`ARH96700-AF-P22523-F1-model_v6.pdb`, 11 962 atoms) with a **263 Angstrom** bounding box, an extended non-globular fold. Applying the APBS `psize` rule (cfac 1.7, fadd 20, 0.5 Angstrom spacing, dime = 32k+1) gives roughly **577 x 417 x 513 = 123 M grid points**, about **58 GB** at the measured 468 bytes per point. My psize estimator reproduces the plan's four tabulated grids to within one dime step, so this is not an artefact of the estimator.

Against 113 GB free RAM, one structure needing ~50 GB is survivable alone and fatal if it overlaps with anything substantial. Predicted memory above 16 GB applies to 4 structures and above 8 GB to 34.

**Runtime is sound.** Independently reconstructing the grid for all 834 from their bounding boxes reproduces the plan's atom statistics exactly (median 2 562, p95 6 082, max 11 962 [VERIFIED]) and gives a single-core total of **9.5 h** using the plan's `5.35 * Mgridpts + 0.9`, or 9.8 h using a fit to the 13 real APBS runs (`7.13 * Mgridpts - 11.34`, R2 0.907). Bounding box over 150 Angstrom: I count 13, the plan says 12. [VERIFIED, with the caveat that the per-structure maximum is 660 to 870 s, not the plan's 928 s, and that the two calibrations differ by 30 % in slope]

**Fix.** Replace "about 20 grid-sized arrays, 28 GB" with the measured **468 MB per million grid points, R2 0.9995**, cite the log files, and schedule on `sum(predicted MB) < 80 GB` with the top 5 structures forced to run alone. Name ARH96700 explicitly as the one that must be solo, and dry-run its `psize` before committing to the overnight schedule.

---

## F13 (minor). Numbers that are stated with the wrong provenance or are slightly off

- **`SurfEpPosMean` 0.234 to 0.335.** Section 5 attributes 0.335 to Option A. Measured, **0.335 is the oracle threshold** and Option A's predicted threshold gives **0.341**. Almost identical, but the plan mixes an oracle number into an Option A "pros/cons" list. [VERIFIED, mislabelled]
- **`SurfEpPosStd` is not repaired by Option A at all.** Naive 0.123, Option A **0.126**, oracle 0.407. The plan lists three features as being fixed; one of them moves by 0.003. Say so.
- **Total surface points.** `apbs_vs_prodes_results.md` says 390 427; the `n_points` column of `apbs_vs_prodes_summary.csv` sums to **391 027**, and the point files contain 391 027 rows. [VERIFIED] Correct the results doc.
- **Feature names.** `SurfEpPosFormalMean` / `SurfEpPosFormalStd` do not exist. The emitted names are `SurfEpPosMeanFormal` / `SurfEpPosStdFormal` (`run.py`, the `standard_features(positive_eps, "SurfEpPos", ...)` call).
- **"0.991 monopole regression is tight enough to stand on 13."** Correct [VERIFIED: Prodes 0.9907, APBS 0.6827; using the PQR charges instead of the Prodes formal charges gives 0.9918 and 0.7356]. This is the one claim in the document that is genuinely robust at n=13, and it should be leaned on harder.

---

## What I would change before running anything

1. Add `n_surface_points` as the null model to every count claim, and add the positive **fraction** as a primary endpoint. If the rule does not beat those, the plan's recommendation changes. **Do this before the 10-hour run**, because it may change what the run is for.
2. Drop the cation exchange framing or acknowledge that no data exists to test it. Rewrite section 6 item 6 around the partial correlation given net charge (F3), which is testable today on 880 rows.
3. Build a real sequence clustering. Section 6 item 1 is currently fiction.
4. Use netQ, not mean EP. Answers open question 2, removes the circularity question, costs nothing measurable.
5. Run a pH 5 subset alongside pH 7. It is the only affordable route to Q > 0 coverage and it tests the monopole explanation directly.
6. Extend `apbs_comparison.py` to sample the Shell reference points, or declare the Shell block out of scope.
7. Fix the memory model to 468 MB per million grid points and force ARH96700 to run alone.
8. Write the decision rule down before the run.
