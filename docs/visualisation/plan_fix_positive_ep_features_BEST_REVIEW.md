# Which review of the EP fix plan did the better work

**Date:** 2026-08-21. Inputs: `plan_fix_positive_ep_features.md` plus `_REVIEW1.md` (statistics and calibration) and `_REVIEW2.md` (computational chemistry and implementation). Output: `_REFINED.md`.

## Ranking

**1. Review 2. 2. Review 1.** Both are excellent and the gap is narrow.

**What this ranking is not.** It is not a filter on findings. Review 1 produced the single most damaging result in either review, and the refined plan is built on it as much as on review 2's. Every blocking finding from both was adopted.

## Why review 2 wins

It **ran the experiment that decided the plan**. The plan recommended patching the threshold first and evaluating Debye screening later. Review 2 implemented the screening, measured it, and showed the recommendation was backwards: screening with **zero fitted constants** reaches a median Jaccard of 0.658 against the APBS positive surface, beating the two-constant threshold patch (0.586) and even an oracle threshold taken from APBS itself (0.601). It also showed the result is flat across lambda from 5 to 12 Angstrom, which kills the plan's own stated objection that screening "adds a parameter that must be chosen".

Crucially it validated its own tooling first: an unscreened reimplementation reproduces the stored Prodes potentials to 0.0000 V, which is what licenses the comparison. I reproduced the whole result independently and it holds exactly.

Three further findings that only a reviewer working at the code and physics level would reach:

- **The negative side is not inherently safe.** The apparent asymmetry (anion exchange fine, cation exchange broken) is a property of a test set in which all 13 proteins are net negative. `tests/data/1GDW`, already in the repository, has net charge +7 and shows the exact mirror failure. This reframes the whole problem and neither the plan nor review 1 saw it.
- **The scope is 10 of the default 54 features, not 3**, because the shell features apply the same `> 0` test on far-field points where the monopole dominates more, and no APBS ground truth was sampled out there at all.
- **`lambda = 3.04/sqrt(I)` is water's Debye length**, but the kernel runs at `eps_r = 4`, where the self-consistent value is 1.77 A. The fix is an empirical correction that mimics screening, not first-principles Poisson-Boltzmann, and saying otherwise would have been a false claim in a released document.

It also answered both direct questions correctly and negatively, which is often more useful than a positive finding: no two-pass restructuring is needed, and the `standard_features` NaN proposal should be dropped.

## What review 1 uniquely contributed, and it is the best single finding

**The plan's headline repair was a protein-size correlation.** `n_surface_points` alone, with no electrostatics at all, scores 0.929 against the APBS positive count. The plan's celebrated "0.641 to 0.973" is a gain of +0.05 over that, with a bootstrap CI spanning zero. On the size-free positive *fraction* the recalibrated rule is slightly **worse** than the naive one, 0.769 against 0.780, while an oracle threshold reaches 0.885.

That is the finding that killed option A, and review 2 reached the same conclusion by a different route but did not frame it as sharply. It is also a general methodological lesson: the plan reported an absolute count against another absolute count, and both are dominated by size.

Review 1 also uniquely produced: the leave-one-out number that showed overfitting was *not* the problem (0.956 against 0.973 in-sample), so the rule was measuring the wrong thing rather than measuring it badly; the demonstration that circularity is a non-issue because netQ works as well as mean EP; the finding that the intercept is a one-protein artefact driven by 4PEP at hat value 0.545; and that mean-centring (option B) is the best zero-constant scheme apart from screening.

Both reviewers independently found that `11_protein_clusters` cannot support the group split, and both quantified the net-charge coverage problem.

## Re-verification table

| Claim | Raised by | Re-check result |
| --- | --- | --- |
| Screening beats the oracle threshold, median Jaccard 0.658 vs 0.601 | R2 | **CONFIRMED.** Reproduced independently: 0.058 naive, 0.618 mean-centred, 0.651/0.658/0.640 at lambda 5/7.85/12 |
| `NSurfPosEp` 0.641 to 0.995 under screening; `SurfEpPosMean` 0.234 to 0.813 | R2 | **CONFIRMED** exactly |
| Unscreened reimplementation reproduces stored values to 0.0 | R2 | **CONFIRMED**, max difference 0.0000 V |
| `n_surface_points` alone scores 0.929 on the count feature | R1 | **CONFIRMED.** Rule gain +0.05, CI [-0.10, +0.24], P(no gain) 0.22 |
| Recalibrated rule is worse than naive on positive *fraction* | R1 | **CONFIRMED**, 0.769 vs 0.780, oracle 0.885 |
| **Memory: 467.5 MB/Mgridpoint** | R1 | **REFUTED.** From APBS's own log counter, which overstates |
| **Memory: 224.4 MB/Mgridpoint** | R2 | **CONFIRMED.** Direct measurement: 8.0 Mpt grid, APBS self-reports 3,755 MB high water, actual peak RSS 1,845 MB = **231 MB/Mgridpoint**. Largest structure needs about 25 GB, not 81 GB. Memory is not the binding constraint and the plan's scheduler is unnecessary |
| `11_protein_clusters` holds 2 clusters covering 4 proteins | R1, R2 | **CONFIRMED** |
| No cation exchange retention data in the repository | R1 | **CONFIRMED.** All 52 Neijenhuis rows are HiTrap Q Sepharose, an anion exchanger |
| Coverage: 6.2 % (R1) vs 8.0 % (R2) net positive at pH 7 | both | **BOTH DEFENSIBLE**, different charge models (sequence vs PROPKA). Either way the held-out basic subset is about 13 proteins |
| 1GDW shows the mirror failure, `NSurfPos` 7020/7039 | R2 | **CONFIRMED** by inspection |

The memory row is the only outright contradiction between the two reviews, and it mattered: the original plan would have built a scheduler against a phantom constraint, using a number that was itself wrong in the crash direction.

## Scorecard

| Criterion | Review 1 | Review 2 |
| --- | --- | --- |
| Verified truth | Reproduced every number, then reinterpreted it | Reproduced every number, then ran the decisive new experiment |
| Blocking value | Killed option A; found the missing CEX data | Reversed the recommendation; found the mirror failure |
| Unique coverage | Size baseline, LOO, leverage, circularity, pH 5 idea | Screening result, shell features, lambda physics, memory truth |
| Actionability | Concrete alternative endpoints | Concrete one-line code change, measured |

## Meta lesson

The plan asked for a critique of a proposal and got something better from both: one reviewer tested the option the plan had **deferred**, and the other tested the baseline the plan had **never stated**.

Review 2's win came from refusing the plan's sequencing. The plan said "implement A first, evaluate C later"; the reviewer evaluated C anyway, and C turned out to dominate. A reviewer who had respected the plan's own priority ordering would have produced a careful critique of the wrong option.

Review 1's win came from asking what a null model would score. The plan compared a feature against another feature and never asked what pure protein size would achieve. That single baseline turned a celebrated result into a non-result.

For the next round, both belong in the instructions: **always test the option the plan defers, and always compute the dumbest possible baseline before believing an improvement.**
