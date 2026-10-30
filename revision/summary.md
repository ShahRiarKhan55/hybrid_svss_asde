# Summary of the revision experiments

All numbers: TEST FPGAs 27–50 (24 devices), primary protocol (held-out 30% test set of about 952 sites per path, additive bounded RBF GPR inside selection and for prediction, selection on FPGA 1 only, error = |mu − v| / d_spec averaged test sites → path → FPGA, mean ± SD over 5 seeds). Design choices were made on TUNE FPGAs 2–26 only; TEST was evaluated once. Details and secondary tables (FPGA 1 included, the paper's old test set): `stage4_report.md`. Stage reports: `stage0_audit.md`, `stage1_report.md`, `stage2_report.md`, `stage3_report.md`, `stage4_report.md`, `stage4_equivalence.md`, `stage5_report.md`.

## 1. Final numbers (TEST, primary protocol)
| method | 95 sites (3%) | 158 (5%) | 222 (7%) | 317 (10%) | first budget < 0.09 | first budget < 0.085 |
|---|---|---|---|---|---|---|
| Hybrid (frozen) | 0.0867 ± 0.0006 | 0.0798 ± 0.0006 | 0.0774 ± 0.0003 | 0.0758 ± 0.0002 | 95 | 158 |
| AL | 0.0897 ± 0.0009 | 0.0881 ± 0.0005 | 0.0867 ± 0.0006 | 0.0855 ± 0.0008 | 95 | -- |
| AL-max | 0.0840 ± 0.0004 | 0.0791 ± 0.0005 | 0.0773 ± 0.0003 | 0.0755 ± 0.0003 | 95 | 95 |
| SVTS | 0.0836 ± 0.0011 | 0.0799 ± 0.0005 | 0.0782 ± 0.0004 | 0.0760 ± 0.0005 | 95 | 95 |
| SDE | 0.0820 ± 0.0005 | 0.0782 ± 0.0003 | 0.0766 ± 0.0004 | 0.0753 ± 0.0002 | 95 | 95 |
| Random | 0.0823 ± 0.0007 | 0.0788 ± 0.0005 | 0.0769 ± 0.0005 | 0.0752 ± 0.0004 | 95 | 95 |
| Maximin | 0.0823 ± 0.0006 | 0.0787 ± 0.0002 | 0.0769 ± 0.0004 | 0.0754 ± 0.0004 | 95 | 95 |
| anchors 1% (exploratory) | 0.0818 ± 0.0010 | 0.0790 ± 0.0005 | 0.0772 ± 0.0005 | 0.0756 ± 0.0003 | 95 | 95 |
| maximin anchors (exploratory) | 0.0823 ± 0.0006 | 0.0789 ± 0.0004 | 0.0772 ± 0.0004 | 0.0756 ± 0.0005 | 95 | 95 |
| constant predictor (reference) | – | – | – | 0.0959 ± 0.0002 | – | – |
| GPR on 50% of the sites, 1587 sites (reference, seed 1 only) | – | – | – | 0.0715 | – | – |

Hybrid = SVSS (95 anchors, alpha = beta = 0.5) + ASDE, frozen unchanged after Stage 3 (no candidate met the pre-written rule). SVTS weights (0.3 similarity, 0.7 variance) were chosen on TUNE. The two exploratory rows were defined before any result and are not part of the headline comparison. Budgets include seed points.

## 2. Ablation (TEST, 317 sites)
| stage 1 (95 sites) \ stage 2 (222 sites) | random fill | fixed SDE (2,2) | adaptive ASDE |
|---|---|---|---|
| random | 0.0752 ± 0.0004 | 0.0753 ± 0.0003 | 0.0755 ± 0.0003 |
| SVSS | 0.0756 ± 0.0003 | 0.0760 ± 0.0003 | 0.0758 ± 0.0002 (hybrid) |

| other cells @ 317 sites | error |
|---|---|
| SVSS-only | 0.0785 ± 0.0007 |
| AL-max | 0.0755 ± 0.0003 |
| ASDE-only(rand seed) | 0.0756 ± 0.0003 |
| ASDE-only(maximin seed) | 0.0756 ± 0.0003 |
| SDE | 0.0753 ± 0.0002 |

Paired differences (path-level, 32 paths, 95% bootstrap CI; Holm over all 46 comparisons of Stage 4):

| contrast @ 317 sites (negative = first lower error) | mean difference | 95% CI | first / second wins (of 32 paths) | Holm p (46 comparisons) |
|---|---|---|---|---|
| Hybrid − rand+SDE | +0.00048 | [+0.00017, +0.00080] | 7 / 25 | 0.069 |
| Hybrid − rand+ASDE | +0.00036 | [+0.00004, +0.00065] | 9 / 23 | 0.45 |
| Hybrid − SVSS+rand | +0.00023 | [-0.00011, +0.00059] | 15 / 17 | 1 |
| Hybrid − SVSS+SDE | -0.00016 | [-0.00052, +0.00021] | 22 / 10 | 1 |
| Hybrid − SVSS-only | -0.00272 | [-0.00345, -0.00193] | 28 / 4 | 4e-05 |
| Hybrid − ASDE-only(rand seed) | +0.00019 | [-0.00018, +0.00052] | 10 / 22 | 1 |
| Hybrid − ASDE-only(maximin seed) | +0.00023 | [-0.00011, +0.00057] | 10 / 22 | 1 |
| AL-max − SVSS-only | -0.00307 | [-0.00376, -0.00240] | 31 / 1 | 4.3e-08 |
| ASDE-only(maximin seed) − ASDE-only(rand seed) | -0.00004 | [-0.00022, +0.00015] | 17 / 15 | 1 |
| Stage-1 main effect: SVSS − random (mean of 3 stage-2 choices) | +0.00045 | [+0.00028, +0.00062] | 4 / 28 | 0.0032 |
| Stage-2 main effect: ASDE − fixed SDE | -0.00002 | [-0.00025, +0.00022] | 18 / 14 | 1 |
| Stage-2 main effect: ASDE − random fill | +0.00023 | [-0.00007, +0.00055] | 13 / 19 | 1 |
| Stage-2 main effect: fixed SDE − random fill | +0.00025 | [-0.00011, +0.00060] | 16 / 16 | 1 |
| Interaction: (SVSS − random) under ASDE minus under fixed SDE | -0.00029 | [-0.00076, +0.00019] | 19 / 13 | 1 |

## 3. Cost (final designs)
| method | GPR fits per path | selection time per path, s (mean over 32 paths of the median of 3 repeats) | min – max over paths | relative to hybrid |
|---|---|---|---|---|
| AL | 285 | 45.033 | 34.275 – 62.137 | 234.43x |
| AL-max | 285 | 44.903 | 36.303 – 53.717 | 233.75x |
| SVTS | 285 | 46.322 | 38.996 – 58.544 | 241.14x |
| Hybrid | 2 | 0.192 | 0.157 – 0.308 | 1.00x |
| SDE | 0 | 0.042 | 0.039 – 0.077 | 0.22x |
| Maximin | 0 | 0.005 | 0.004 – 0.007 | 0.02x |
| Random | 0 | 0.000 (below 1 ms) | 0.000 – 0.000 | 0.00x |

Final GPR fit on 317 sites (same for every method): 0.323 s (min–max 0.221–0.608); prediction on 2856 sites: 0.055 s. One thread, sequential, 3 repeats, median. Selection times varied by about 30% between repeats (mean (max − min)/median over path and method pairs), the machine ran on the Balanced power plan at 1.8–2.3 GHz with other services active, so treat the absolute seconds as indicative; the fit counts and the ratio (about 234x for AL, 241x for SVTS) are the robust numbers. Details: `stage5_report.md`.

## 4. Equivalence (TOST, hybrid vs the others, TEST path-level differences)
**The margin of ±0.001 in error (about 1.3% of the hybrid's 10%-budget error) is a choice, not derived from data.** Paired-t TOST on the 32 path-level differences, Holm over the 15 tests; a Wilcoxon-based version gives the same 4 equivalences (`stage4_equivalence.md`).

| budget | hybrid vs | mean difference | 90% CI | Holm p (TOST) | equivalent within ±0.001 | smallest margin that would give equivalence |
|---|---|---|---|---|---|---|
| 158 | SVTS | -0.00005 | [-0.00073, +0.00064] | 0.102 | no | 0.00073 |
| 158 | AL-max | +0.00069 | [-0.00005, +0.00142] | 1 | no | 0.00142 |
| 158 | SDE | +0.00162 | [+0.00102, +0.00222] | 1 | no | 0.00222 |
| 158 | Random | +0.00101 | [+0.00053, +0.00149] | 1 | no | 0.00149 |
| 158 | Maximin | +0.00107 | [+0.00063, +0.00151] | 1 | no | 0.00151 |
| 222 | SVTS | -0.00080 | [-0.00146, -0.00013] | 1 | no | 0.00146 |
| 222 | AL-max | +0.00002 | [-0.00054, +0.00057] | 0.0315 | yes | 0.00057 |
| 222 | SDE | +0.00072 | [+0.00035, +0.00110] | 0.656 | no | 0.00110 |
| 222 | Random | +0.00048 | [+0.00022, +0.00074] | 0.0153 | yes | 0.00074 |
| 222 | Maximin | +0.00042 | [+0.00001, +0.00083] | 0.102 | no | 0.00083 |
| 317 | SVTS | -0.00020 | [-0.00073, +0.00034] | 0.0858 | no | 0.00073 |
| 317 | AL-max | +0.00035 | [-0.00009, +0.00080] | 0.102 | no | 0.00080 |
| 317 | SDE | +0.00050 | [+0.00024, +0.00076] | 0.0211 | yes | 0.00076 |
| 317 | Random | +0.00059 | [+0.00029, +0.00088] | 0.102 | no | 0.00088 |
| 317 | Maximin | +0.00042 | [+0.00011, +0.00073] | 0.0211 | yes | 0.00073 |

Equivalent within ±0.001 after Holm: hybrid vs AL-max at 222 sites, vs Random at 222, vs SDE at 317, vs Maximin at 317 (4 of 15). Not shown for SVTS at any budget (at 317 the 90% CI is [−0.00073, +0.00034], Holm p = 0.086), for the others at 158, or for Random and AL-max at 317 (Holm p = 0.10).

## 5. Claims
**Supported by the data**
1. The hybrid has lower error than AL at every budget (−0.0030 to −0.0097; wins on 28–31 of 32 paths; Holm p ≤ 1.3e-6).
2. Selection cost: the hybrid needs 2 GPR fits and about 0.19 s per path, against 285 fits and about 45–46 s for AL, AL-max and SVTS (roughly 234–241 times faster on this machine). SDE, Maximin and Random need no GPR fit and are cheaper still.
3. From 158 sites on, no difference between the hybrid and SVTS or AL-max was detected (|difference| ≤ 0.0008, Holm p = 1). The error of AL-max is statistically equivalent to the hybrid's within ±0.001 at 222 sites only.
4. At 222 sites the hybrid is equivalent to Random within ±0.001 (difference +0.00048), and at 317 sites it is equivalent to SDE and to Maximin (+0.00050 and +0.00042). At the other budgets or comparators equivalence was not established.
5. SVSS's similarity term adds nothing: SVSS-only is worse than the hybrid (+0.0027, Holm p = 4e-5) and plain highest-variance AL-max beats SVSS-only by 0.0031 on 31 of 32 paths.
6. SVSS anchors are slightly worse than random anchors (+0.00045, Holm p = 0.0032); the choice of stage-2 rule (adaptive ASDE, fixed SDE, random fill) makes no detectable difference at 317 sites; ASDE-only matches the hybrid.

**Not supported by the data**
1. Lowest error. At 317 sites: Random 0.0752, SDE 0.0753, Maximin 0.0754, AL-max 0.0755, hybrid 0.0758, SVTS 0.0760, AL 0.0855. The hybrid is the best of the seven on 0–1 of 32 paths at every budget.
2. Earlier convergence. At 95 sites the hybrid is worse than Random, SDE, Maximin, AL-max and SVTS (+0.0026 to +0.0047, Holm-significant); it first falls below 0.085 at 158 sites, the others (except AL) at 95.
3. Better than SVTS. Worse at 95 sites (+0.0031) and no difference later; equivalence with SVTS within ±0.001 is not shown either.
4. Equal to or better than SDE, Random or Maximin at small budgets: the hybrid is worse than SDE at 95, 158 and 222 sites and worse than Random and Maximin at 95 and 158 sites (Holm-significant). At 317 sites it is nominally 0.0004–0.0006 above them; the two-sided tests are not significant after Holm (p = 0.067–0.19) while equivalence is shown for SDE and Maximin only.
5. Tighter spread. SD across the 24 FPGAs at 317 sites 0.0041 for the hybrid, Random and SDE; SD across paths 0.0241 for the hybrid and Random.
6. A benefit of the adaptive ASDE windows, of the SVSS anchors, or of the cosine-without-centring similarity (Stage 3: no candidate met the rule).
7. That the exploratory variants (anchors 1%, maximin anchors) beat Random; they remove the hybrid's early deficit but are at the Random level.
8. Any inference beyond these 32 paths and this device family: FPGAs reuse the same site sets, so agreement across FPGAs (the hybrid is worse than Random on 24 of 24 FPGAs at 317 sites) is not independent evidence.

**Limits of this work.** One device family, one predictor, 5 seeds. The equivalence margin is a choice. Timings come from a non-idle laptop in the Balanced power plan with about 30% repeat-to-repeat variation. The 50%-GPR reference uses one seed. Stage 3 searched only the trimmed list of candidates agreed beforehand.

## 6. Predictor comparison on the frozen sites (Stage 6; all 24 TEST FPGAs)
Sites: the first 317 sites of the saved FPGA-1 orders, 5 seeds, for the frozen hybrid and for Random; same sites for every predictor; held-out 30% test set; error = |mu − v| / d_spec averaged test sites → path → FPGA → mean over the 24 FPGAs; mean ± SD over the 5 seeds. Paired differences: path-level (32 paths, mean over seeds and FPGAs), 95% bootstrap CI, Wilcoxon, Holm over 16 comparisons (10 versus the additive GPR, 6 site contrasts). Full report: `stage6_report.md`; wider-grid report: `stage6_wide_report.md`.

### 6.1 Hyperparameter tuning and cross-validation actually used
| predictor | tuned hyperparameters (grid) | cross-validation | fixed settings |
|---|---|---|---|
| SVR (RBF kernel) | C ∈ {0.1, 1, 10, 100}; gamma ∈ {0.01, 0.1, 1, 10}; epsilon ∈ {0.01, 0.1, 0.3} (48 settings) | 5-fold `KFold(shuffle=True, random_state=0)` on the 317 training sites, scoring = negative mean absolute error of the predictions in the original value units | inputs standardised (`StandardScaler`, fitted inside each fold), target standardised (`TransformedTargetRegressor`); refit on all 317 sites with the best setting |
| Random forest | min_samples_leaf ∈ {1, 3, 5, 10}; max_features ∈ {1.0, 0.5} (8 settings) | same | 200 trees, `random_state=0`; same input/target standardisation |
| MLP | hidden layers (16,16), (32,32), (64,64); alpha ∈ {0.001, 0.1, 10} (9 settings) | same | tanh, Adam (scikit-learn default solver), `max_iter=1000`, no early stopping, `random_state=0`; same standardisation |
| SVR, wider grid (FPGAs 27 and 29) | C ∈ {0.1, 1, 10, 100, 1000}; gamma ∈ {0.03, 0.1, 0.3, 1, 3, 10}; epsilon ∈ {0.001, 0.01, 0.1} (90 settings) | same | as above |
| MLP, wider grid (FPGAs 27 and 29) | widths (16,16), (32,32), (64,64), (128,128); alpha ∈ {1e-5, 1e-4, 1e-3, 1e-2, 1e-1}; activation ∈ {tanh, relu} (40 settings) | same | Adam, `max_iter=1000`, `random_state=0` |
| Additive and one-lengthscale GPR | none (hyperparameters by marginal likelihood) | – | same bounded GPR in all stages: lengthscale 0.5–300, kernel variance 1e-4–1e3, noise 1e-4–10, at most 200 iterations, no restarts |
| Deep GP | none (variational bound) | – | layers [1,3,2], 50 inducing points, 200 iterations, isotropic RBF kernels |

Tuning was done separately for every FPGA, path, seed and site type, on the 317 training sites only (no test site is used). The grids are my choice.

### 6.2 Error
| predictor | hybrid sites | Random sites |
|---|---|---|
| Additive RBF GPR | 0.0758 ± 0.0002 | 0.0752 ± 0.0004 |
| One-lengthscale RBF GPR | 0.0797 ± 0.0002 | 0.0798 ± 0.0002 |
| SVR (RBF) | 0.0757 ± 0.0003 | 0.0758 ± 0.0003 |
| Random forest | 0.0746 ± 0.0005 | 0.0743 ± 0.0002 |
| MLP (tanh, 2 hidden layers) | 0.0789 ± 0.0002 | 0.0787 ± 0.0001 |
| Deep GP [1,3,2] | 0.0798 ± 0.0002 | 0.0804 ± 0.0002 |

### 6.3 Paired differences
| predictor − additive GPR | hybrid sites: mean difference [95% CI], paths won (predictor / GPR), Holm p | Random sites |
|---|---|---|
| One-lengthscale RBF GPR | +0.00393 [+0.00303, +0.00490], 0/32 paths, Holm p = 7.5e-09 | +0.00460 [+0.00368, +0.00559], 0/32 paths, Holm p = 7.5e-09 |
| SVR (RBF) | -0.00011 [-0.00151, +0.00093], 11/21 paths, Holm p = 0.44 | +0.00059 [-0.00086, +0.00169], 5/27 paths, Holm p = 0.0032 |
| Random forest | -0.00124 [-0.00171, -0.00073], 27/5 paths, Holm p = 0.0012 | -0.00093 [-0.00130, -0.00055], 23/9 paths, Holm p = 0.00074 |
| MLP (tanh, 2 hidden layers) | +0.00306 [+0.00233, +0.00385], 1/31 paths, Holm p = 1.1e-08 | +0.00346 [+0.00274, +0.00424], 1/31 paths, Holm p = 1.1e-08 |
| Deep GP [1,3,2] | +0.00398 [+0.00320, +0.00480], 0/32 paths, Holm p = 7.5e-09 | +0.00521 [+0.00424, +0.00624], 0/32 paths, Holm p = 7.5e-09 |

Effect of the site choice for each predictor:

| predictor | hybrid sites − Random sites: mean difference [95% CI], paths where hybrid sites win/lose, Holm p |
|---|---|
| Additive RBF GPR | +0.00059 [+0.00027, +0.00094], 11/21 paths, Holm p = 0.014 |
| One-lengthscale RBF GPR | -0.00009 [-0.00044, +0.00026], 18/14 paths, Holm p = 0.95 |
| SVR (RBF) | -0.00011 [-0.00031, +0.00010], 19/13 paths, Holm p = 0.76 |
| Random forest | +0.00028 [-0.00014, +0.00069], 12/20 paths, Holm p = 0.58 |
| MLP (tanh, 2 hidden layers) | +0.00018 [-0.00013, +0.00052], 15/17 paths, Holm p = 0.95 |
| Deep GP [1,3,2] | -0.00065 [-0.00097, -0.00033], 22/10 paths, Holm p = 0.008 |

### 6.4 Wider-grid check (FPGAs 27 and 29 only, 32 paths, 5 seeds)
| predictor (FPGAs 27 and 29) | hybrid sites | Random sites |
|---|---|---|
| Additive RBF GPR | 0.0761 ± 0.0002 | 0.0756 ± 0.0004 |
| SVR, Stage 6 grid | 0.0762 ± 0.0003 | 0.0764 ± 0.0003 |
| SVR, wider grid | 0.0761 ± 0.0004 | 0.0762 ± 0.0003 |
| MLP, Stage 6 grid (tanh) | 0.0797 ± 0.0002 | 0.0796 ± 0.0002 |
| MLP, wider grid (tanh or relu) | 0.0788 ± 0.0003 | 0.0787 ± 0.0002 |

| wider-grid comparison | site | mean difference [95% CI] | paths won (first/second) | Holm p (10 tests) |
|---|---|---|---|---|
| SVR (wider grid) − additive GPR | hybrid | -0.00008 [-0.00160, +0.00114] | 13/19 | 0.28 |
| MLP (wider grid) − additive GPR | hybrid | +0.00271 [+0.00187, +0.00358] | 2/30 | 2.3e-08 |
| SVR wider grid − SVR Stage 6 grid | hybrid | -0.00015 [-0.00028, -0.00002] | 22/10 | 0.056 |
| MLP wider grid − MLP Stage 6 grid | hybrid | -0.00084 [-0.00114, -0.00054] | 27/5 | 8e-05 |
| SVR (wider grid) − MLP (wider grid) | hybrid | -0.00278 [-0.00411, -0.00175] | 30/2 | 2.3e-08 |
| SVR (wider grid) − additive GPR | random | +0.00064 [-0.00093, +0.00189] | 7/25 | 0.0062 |
| MLP (wider grid) − additive GPR | random | +0.00316 [+0.00231, +0.00408] | 2/30 | 2.3e-08 |
| SVR wider grid − SVR Stage 6 grid | random | -0.00017 [-0.00029, -0.00005] | 23/9 | 0.047 |
| MLP wider grid − MLP Stage 6 grid | random | -0.00083 [-0.00109, -0.00055] | 29/3 | 2.8e-05 |
| SVR (wider grid) − MLP (wider grid) | random | -0.00252 [-0.00373, -0.00159] | 31/1 | 2.3e-08 |

With the wider grids the order additive GPR < SVR < MLP is unchanged on Random sites; on hybrid sites SVR and the additive GPR are tied (difference −0.00008, Holm p = 0.28); the MLP remains last (+0.0027 on hybrid sites, +0.0032 on Random sites above the additive GPR). The wider SVR grid is not a strict superset of the Stage 6 grid (gamma 0.01 and epsilon 0.3 were dropped; the best epsilon is at the upper grid value 0.1 in 56% of fits). The MLP optimum is at a wider-grid edge in about half of the fits (width 128: 34%, alpha 0.1: 28%; relu chosen in 85%).

### 6.5 Claims about the predictors
**Supported:** (1) The random forest has the lowest error on both site types (−0.00124 on hybrid sites, −0.00093 on Random sites against the additive GPR; Holm p = 0.0012 and 0.00074; about 1.2–1.6% of the error). (2) SVR and the additive GPR are not separable on the hybrid sites; on Random sites the additive GPR is better by 0.00059 (Holm p = 0.0032); this does not change with the wider grids (checked on 2 FPGAs). (3) The additive structure matters: the one-lengthscale GPR is 0.0039–0.0046 worse than the additive GPR; the MLP and the deep GP are 0.003–0.005 worse, on 31–32 of 32 paths (Holm p ≤ 1.1e-8). (4) The GPR predictive intervals are close to nominal (90.8–91.1% and 94.4–94.6%). (5) Only the additive GPR (+0.00059) and the deep GP (−0.00065) show a significant difference between hybrid and Random sites; no predictor shows a hybrid advantage except the deep GP, which is the worst predictor.
**Not supported:** (1) That the additive GPR is the most accurate predictor (the random forest is better). (2) That the ranking of predictors is independent of the site choice in every detail (SVR and the additive GPR swap places). (3) That the MLP or the deep GP can match the additive GPR, even with wider MLP grids. (4) That the wider-grid result holds on all 24 FPGAs (it was run on 2). (5) That the hybrid sites are better than Random sites for the best predictor (random forest: +0.00028, Holm p = 0.58).

## 7. Exact numbers to quote in the paper (Stage 6; TEST FPGAs, 24 devices, 317 sites, held-out test set)
| predictor | error, hybrid sites | error, Random sites | difference to additive GPR, hybrid sites [95% CI] | difference to additive GPR, Random sites [95% CI] | fit incl. tuning, s | predict, s |
|---|---|---|---|---|---|---|
| Additive RBF GPR | 0.0758 ± 0.0002 | 0.0752 ± 0.0004 | – | – | 0.25 | 0.051 |
| One-lengthscale RBF GPR | 0.0797 ± 0.0002 | 0.0798 ± 0.0002 | +0.00393 [+0.00303, +0.00490] (Holm p = 7.5e-09) | +0.00460 [+0.00368, +0.00559] (Holm p = 7.5e-09) | 0.22 | 0.035 |
| SVR (RBF) | 0.0757 ± 0.0003 | 0.0758 ± 0.0003 | -0.00011 [-0.00151, +0.00093] (Holm p = 0.44) | +0.00059 [-0.00086, +0.00169] (Holm p = 0.0032) | 1.85 | 0.027 |
| Random forest | 0.0746 ± 0.0005 | 0.0743 ± 0.0002 | -0.00124 [-0.00171, -0.00073] (Holm p = 0.0012) | -0.00093 [-0.00130, -0.00055] (Holm p = 0.00074) | 5.10 | 0.019 |
| MLP (tanh, 2 hidden layers) | 0.0789 ± 0.0002 | 0.0787 ± 0.0001 | +0.00306 [+0.00233, +0.00385] (Holm p = 1.1e-08) | +0.00346 [+0.00274, +0.00424] (Holm p = 1.1e-08) | 7.28 | 0.003 |
| Deep GP [1,3,2] | 0.0798 ± 0.0002 | 0.0804 ± 0.0002 | +0.00398 [+0.00320, +0.00480] (Holm p = 7.5e-09) | +0.00521 [+0.00424, +0.00624] (Holm p = 7.5e-09) | 5.34 | 0.000 |

| quantity | value |
|---|---|
| 90% / 95% interval coverage, additive GPR (hybrid / Random sites) | 90.9% / 94.5% and 90.9% / 94.5% |
| 90% / 95% interval coverage, one-lengthscale GPR (hybrid / Random sites) | 90.8% / 94.4% and 91.1% / 94.6% |
| SVR, wider grid (FPGAs 27, 29): hybrid / Random sites | 0.0761 / 0.0762 (additive GPR 0.0761 / 0.0756) |
| MLP, wider grid (FPGAs 27, 29): hybrid / Random sites | 0.0788 / 0.0787 (additive GPR 0.0761 / 0.0756) |

Notes for quoting: errors are mean ± SD over 5 seeds of the FPGA-averaged error; differences are path-level (32 paths) with 95% bootstrap CIs and Holm-adjusted Wilcoxon p-values (16 comparisons); times are for one thread, sequential, median of 3 repeats on FPGA 27 (8 paths) and include cross-validation where applicable (the deep GP fits and predicts in one call, its time is listed under fit); the wider-grid rows use FPGAs 27 and 29 only.
