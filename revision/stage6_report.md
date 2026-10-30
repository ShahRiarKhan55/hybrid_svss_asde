# Stage 6 – Predictor comparison on the frozen sites (TEST FPGAs: 24 of the 24)

## Protocol
- Sites: the saved FPGA-1 selection orders, first 317 sites, 5 seeds: **hybrid** (the frozen current hybrid) and **Random**. Same sites for every predictor. Mean training sites per FPGA {'hybrid': 317.0, 'random': 317.0}.
- Test: the held-out 30% sites (about 952 per path), 24 of the 24 TEST FPGAs (27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50); error = |mu − v| / d_spec averaged test sites → path → FPGA → mean over the 24 FPGAs; mean ± SD over the 5 seeds.
- Predictors: additive RBF GPR (the bounded GPR used in all stages); one-lengthscale RBF GPR (same normaliser, bounds and 200 iterations); SVR (RBF; C ∈ {0.1, 1, 10, 100}, gamma ∈ {0.01, 0.1, 1, 10}, epsilon ∈ {0.01, 0.1, 0.3}; scaled inputs, standardised target); random forest (200 trees; min_samples_leaf ∈ {1, 3, 5, 10}, max_features ∈ {1.0, 0.5}); MLP (tanh, two hidden layers of width w, w ∈ {16, 32, 64}, alpha ∈ {0.001, 0.1, 10}, Adam, up to 1000 iterations); deep GP (layers [1,3,2], 50 inducing points, 200 iterations, isotropic RBF kernels, as in `predictor_comparison`). **Assumptions:** the grids above are my choice; SVR, RF and MLP are tuned by 5-fold CV (shuffled, fixed seed) with mean absolute error as the score, on the training sites only, separately for every FPGA, path, seed and site type, and the final model is refit on all training sites.
- Paired differences: path-level values (mean over seeds and the 24 FPGAs), 32 paths, 95% bootstrap CI (10,000 resamples, seed 0), Wilcoxon signed-rank (two-sided); Holm correction over all 16 comparisons below (10 predictor-vs-additive-GPR comparisons plus 6 hybrid-vs-Random site contrasts).

## Error on the hybrid sites
| predictor | error (mean ± SD over 5 seeds) | rank |
|---|---|---|
| Additive RBF GPR | 0.0758 ± 0.0002 | 3 |
| One-lengthscale RBF GPR | 0.0797 ± 0.0002 | 5 |
| SVR (RBF) | 0.0757 ± 0.0003 | 2 |
| Random forest | 0.0746 ± 0.0005 | 1 |
| MLP (tanh, 2 hidden layers) | 0.0789 ± 0.0002 | 4 |
| Deep GP [1,3,2] | 0.0798 ± 0.0002 | 6 |

## Error on Random sites
| predictor | error (mean ± SD over 5 seeds) | rank |
|---|---|---|
| Additive RBF GPR | 0.0752 ± 0.0004 | 2 |
| One-lengthscale RBF GPR | 0.0798 ± 0.0002 | 5 |
| SVR (RBF) | 0.0758 ± 0.0003 | 3 |
| Random forest | 0.0743 ± 0.0002 | 1 |
| MLP (tanh, 2 hidden layers) | 0.0787 ± 0.0001 | 4 |
| Deep GP [1,3,2] | 0.0804 ± 0.0002 | 6 |

Ranking by site type: the ranking differs between hybrid and Random sites (see the rank columns).

## Paired differences to the additive GPR
| comparison (negative = first lower error) | mean diff | 95% CI | first / second wins (of 32) | Wilcoxon p | Holm p (16 tests) | significant |
|---|---|---|---|---|---|---|
| One-lengthscale RBF GPR − additive GPR, hybrid sites | +0.00393 | [+0.00303, +0.00490] | 0 / 32 | 4.7e-10 | 7.5e-09 | yes |
| SVR (RBF) − additive GPR, hybrid sites | -0.00011 | [-0.00151, +0.00093] | 11 / 21 | 0.087 | 0.44 | no |
| Random forest − additive GPR, hybrid sites | -0.00124 | [-0.00171, -0.00073] | 27 / 5 | 0.00013 | 0.0012 | yes |
| MLP (tanh, 2 hidden layers) − additive GPR, hybrid sites | +0.00306 | [+0.00233, +0.00385] | 1 / 31 | 9.3e-10 | 1.1e-08 | yes |
| Deep GP [1,3,2] − additive GPR, hybrid sites | +0.00398 | [+0.00320, +0.00480] | 0 / 32 | 4.7e-10 | 7.5e-09 | yes |
| One-lengthscale RBF GPR − additive GPR, random sites | +0.00460 | [+0.00368, +0.00559] | 0 / 32 | 4.7e-10 | 7.5e-09 | yes |
| SVR (RBF) − additive GPR, random sites | +0.00059 | [-0.00086, +0.00169] | 5 / 27 | 0.0004 | 0.0032 | yes |
| Random forest − additive GPR, random sites | -0.00093 | [-0.00130, -0.00055] | 23 / 9 | 7.4e-05 | 0.00074 | yes |
| MLP (tanh, 2 hidden layers) − additive GPR, random sites | +0.00346 | [+0.00274, +0.00424] | 1 / 31 | 9.3e-10 | 1.1e-08 | yes |
| Deep GP [1,3,2] − additive GPR, random sites | +0.00521 | [+0.00424, +0.00624] | 0 / 32 | 4.7e-10 | 7.5e-09 | yes |

## Does the site choice matter for each predictor? (hybrid sites − Random sites)
| comparison (negative = first lower error) | mean diff | 95% CI | first / second wins (of 32) | Wilcoxon p | Holm p (16 tests) | significant |
|---|---|---|---|---|---|---|
| Additive RBF GPR: hybrid sites − Random sites | +0.00059 | [+0.00027, +0.00094] | 11 / 21 | 0.0024 | 0.014 | yes |
| One-lengthscale RBF GPR: hybrid sites − Random sites | -0.00009 | [-0.00044, +0.00026] | 18 / 14 | 0.6 | 0.95 | no |
| SVR (RBF): hybrid sites − Random sites | -0.00011 | [-0.00031, +0.00010] | 19 / 13 | 0.25 | 0.76 | no |
| Random forest: hybrid sites − Random sites | +0.00028 | [-0.00014, +0.00069] | 12 / 20 | 0.14 | 0.58 | no |
| MLP (tanh, 2 hidden layers): hybrid sites − Random sites | +0.00018 | [-0.00013, +0.00052] | 15 / 17 | 0.48 | 0.95 | no |
| Deep GP [1,3,2]: hybrid sites − Random sites | -0.00065 | [-0.00097, -0.00033] | 22 / 10 | 0.0011 | 0.008 | yes |

## Interval coverage of the GPRs (held-out test sites, predictive variance including noise)
| GPR | sites | 90% interval coverage | 95% interval coverage |
|---|---|---|---|
| Additive RBF GPR | hybrid | 90.9% | 94.5% |
| One-lengthscale RBF GPR | hybrid | 90.8% | 94.4% |
| Additive RBF GPR | random | 90.9% | 94.5% |
| One-lengthscale RBF GPR | random | 91.1% | 94.6% |

## Fit and predict time per path (means over all fits of the main run on hybrid sites; 8 jobs ran in parallel on a laptop, so these are indicative; see the sequential timing below if present)
| predictor | fit incl. tuning, s | predict, s |
|---|---|---|
| Additive RBF GPR | 0.46 | 0.080 |
| One-lengthscale RBF GPR | 0.49 | 0.058 |
| SVR (RBF) | 3.19 | 0.044 |
| Random forest | 8.65 | 0.031 |
| MLP (tanh, 2 hidden layers) | 15.38 | 0.006 |
| Deep GP [1,3,2] | 8.37 | 0.000 |

## Tuned hyperparameters at a grid edge
| predictor | hyperparameter | share of fits with the best value at a grid edge |
|---|---|---|
| SVR (RBF) | C (grid 0.1 … 100) | 30% |
| SVR (RBF) | gamma (grid 0.01 … 10) | 22% |
| SVR (RBF) | epsilon (grid 0.01 … 0.3) | 64% |
| Random forest | min_samples_leaf (grid 1 … 10) | 31% |
| Random forest | max_features (grid 1.0 … 0.5) | 100% |
| MLP (tanh, 2 hidden layers) | hidden_layer_sizes (grid [16, 16] … [64, 64]) | 72% |
| MLP (tanh, 2 hidden layers) | alpha (grid 0.001 … 10.0) | 61% |

## Secondary: paper's old test set (all non-training sites)
| predictor | hybrid sites | Random sites |
|---|---|---|
| Additive RBF GPR | 0.0764 ± 0.0002 | 0.0757 ± 0.0003 |
| One-lengthscale RBF GPR | 0.0803 ± 0.0002 | 0.0802 ± 0.0002 |
| SVR (RBF) | 0.0764 ± 0.0002 | 0.0763 ± 0.0003 |
| Random forest | 0.0751 ± 0.0003 | 0.0746 ± 0.0001 |
| MLP (tanh, 2 hidden layers) | 0.0794 ± 0.0002 | 0.0791 ± 0.0001 |
| Deep GP [1,3,2] | 0.0804 ± 0.0001 | 0.0810 ± 0.0002 |

## Sequential timing (one thread, nothing else started by me, FPGA 27, hybrid seed-1 sites, 8 paths, median of 3 repeats; measured in the first Stage 6 session)
| predictor | fit incl. tuning, s | predict on ~2,850 sites, s |
|---|---|---|
| Additive RBF GPR | 0.25 | 0.051 |
| One-lengthscale RBF GPR | 0.22 | 0.035 |
| SVR (RBF) | 1.85 | 0.027 |
| Random forest | 5.10 | 0.019 |
| MLP (tanh, 2 hidden layers) | 7.28 | 0.003 |
| Deep GP [1,3,2] | 5.34 | 0.000 |

The deep GP function fits and predicts in one call; its time is reported under fit.

## Reading the tables (plain words, all 24 TEST FPGAs)
- On the hybrid sites the random forest has the lowest error (0.0746), then SVR (0.0757) and the additive GPR (0.0758); the MLP (0.0789), the one-lengthscale GPR (0.0797) and the deep GP (0.0798) are clearly worse. On Random sites the order is random forest 0.0743, additive GPR 0.0752, SVR 0.0758, MLP 0.0787, one-lengthscale GPR 0.0798, deep GP 0.0804.
- Random forest − additive GPR: −0.00124 on hybrid sites (CI −0.00171 to −0.00073, better on 27 of 32 paths, Holm p = 0.0012) and −0.00093 on Random sites (23 of 32, Holm p = 0.00074); about 1.2–1.6% of the error.
- SVR and the additive GPR are not separable on the hybrid sites (−0.00011, CI −0.00151 to +0.00093, Holm p = 0.44). On Random sites the additive GPR is better than SVR by 0.00059 (Holm p = 0.0032), although the bootstrap CI of the mean includes 0 (the Wilcoxon test is driven by the median shift).
- The one-lengthscale GPR, the MLP and the deep GP are worse than the additive GPR by 0.003–0.005 on both site types, on 31–32 of 32 paths (all Holm p ≤ 1.1e-8).
- Does the ranking depend on the site choice? The random forest is first on both. The only change is that SVR and the additive GPR swap places between second and third, which is within noise on the hybrid sites. Hybrid sites − Random sites: additive GPR +0.00059 (hybrid sites slightly worse, Holm p = 0.014) and deep GP −0.00065 (hybrid sites slightly better, Holm p = 0.008); no significant site effect for the other four predictors.
- GPR interval coverage is close to nominal for both GPRs and both site types: 90.8–91.1% for the 90% interval and 94.4–94.6% for the 95% interval.
- The additive GPR error on the hybrid sites (0.0758) equals the Stage 4 value, as it should.
- The wider-grid check (below, FPGAs 27 and 29): with wider SVR and MLP grids SVR and the additive GPR are tied on hybrid sites (0.0761 vs 0.0761; the order of the two flips in the fifth decimal) and the additive GPR stays ahead of SVR on Random sites (0.0756 vs 0.0762). The MLP improves by 0.0008 but remains 0.0027–0.0032 above the additive GPR (Holm p ≤ 2.3e-8).

## Things to know
- The grids of Stage 6 are my choice and coarse (MLP best width at a grid edge in 72% of fits, alpha at its smallest value in 61%; SVR epsilon at an edge in 64%). The wider-grid check addresses this only for FPGAs 27 and 29.
- The wider SVR grid is not a strict superset of the Stage 6 grid: C gains 1000; gamma gains 0.03, 0.3 and 3 but loses 0.01; epsilon becomes {0.001, 0.01, 0.1} and loses 0.3. The best epsilon is the upper edge (0.1) in 56% of wider-grid fits and at an edge in 78%, so a larger epsilon might still help SVR. The MLP grid gains width 128, alpha 1e-5/1e-4/1e-2 and relu, and loses alpha = 10 (chosen in 1% of Stage 6 fits).
- The MLP optimum is at the edge of the wider grid in about half of the fits (width 128 in 34%, alpha 0.1 in 28%), and relu is chosen in 85% of fits.
- The deep GP is the isotropic-kernel version from the earlier stages (not ARD).
- Times in the main-run table are with 8 jobs in parallel (indicative); the sequential table is the one to cite. Holm correction covers 16 comparisons in the main table and 10 in the wider-grid check.
- The main run was done in two sessions (8 FPGAs, then the remaining 16) with identical code, seeds and sites.

## Wider-grid check (FPGAs 27 and 29 only)
### Protocol
Same sites (hybrid and Random, 5 seeds, 317 sites), held-out test set, error definition, 5-fold CV on the training sites (mean absolute error), standardised inputs/target, FPGAs 27 and 29, all 32 paths. Wider grids: SVR with C ∈ {0.1, 1, 10, 100, 1000}, gamma ∈ {0.03, 0.1, 0.3, 1, 3, 10}, epsilon ∈ {0.001, 0.01, 0.1}; MLP (Adam, up to 1000 iterations) with hidden widths (16,16), (32,32), (64,64), (128,128), alpha ∈ {1e-5, 1e-4, 1e-3, 1e-2, 1e-1} and activation ∈ {tanh, relu} chosen by CV. The "Stage 6 grid" rows and the additive GPR are the Stage 6 results for the same two FPGAs. Errors are the mean over the 2 FPGAs; SD over the 5 seeds; paired differences are path-level (32 paths, mean over seeds and the 2 FPGAs) with a 95% bootstrap CI (10,000 resamples) and Wilcoxon, Holm over the 10 comparisons below. Only 2 devices were used, so the numbers differ from the 24-FPGA table.

### Error
| predictor | hybrid sites | rank among the 3 families* | Random sites | rank* |
|---|---|---|---|---|
| Additive RBF GPR | 0.0761 ± 0.0002 | 1 | 0.0756 ± 0.0004 | 1 |
| SVR, Stage 6 grid | 0.0762 ± 0.0003 | 2 | 0.0764 ± 0.0003 | 2 |
| SVR, wider grid | 0.0761 ± 0.0004 | 1 | 0.0762 ± 0.0003 | 2 |
| MLP, Stage 6 grid (tanh) | 0.0797 ± 0.0002 | 3 | 0.0796 ± 0.0002 | 3 |
| MLP, wider grid (tanh or relu) | 0.0788 ± 0.0003 | 3 | 0.0787 ± 0.0002 | 3 |

\* rank among the three families (additive GPR, SVR, MLP), using the Stage 6 grids for the SVR and MLP rows marked "Stage 6 grid" and the wider grids for the rows marked "wider grid".

Order of the three families (lowest error first): hybrid sites — Stage 6 grids: Additive GPR < SVR < MLP; wider grids: SVR < Additive GPR < MLP. Random sites — Stage 6 grids: Additive GPR < SVR < MLP; wider grids: Additive GPR < SVR < MLP.

### Paired differences
| site | comparison (negative = first lower error) | mean diff | 95% CI | first / second wins (of 32) | Wilcoxon p | Holm p (10 tests) | significant |
|---|---|---|---|---|---|---|---|
| hybrid | SVR (wider grid) − additive GPR | -0.00008 | [-0.00160, +0.00114] | 13 / 19 | 0.28 | 0.28 | no |
| hybrid | MLP (wider grid) − additive GPR | +0.00271 | [+0.00187, +0.00358] | 2 / 30 | 3.3e-09 | 2.3e-08 | yes |
| hybrid | SVR wider grid − SVR Stage 6 grid | -0.00015 | [-0.00028, -0.00002] | 22 / 10 | 0.028 | 0.056 | no |
| hybrid | MLP wider grid − MLP Stage 6 grid | -0.00084 | [-0.00114, -0.00054] | 27 / 5 | 1.6e-05 | 8e-05 | yes |
| hybrid | SVR (wider grid) − MLP (wider grid) | -0.00278 | [-0.00411, -0.00175] | 30 / 2 | 2.3e-09 | 2.3e-08 | yes |
| random | SVR (wider grid) − additive GPR | +0.00064 | [-0.00093, +0.00189] | 7 / 25 | 0.0015 | 0.0062 | yes |
| random | MLP (wider grid) − additive GPR | +0.00316 | [+0.00231, +0.00408] | 2 / 30 | 2.3e-09 | 2.3e-08 | yes |
| random | SVR wider grid − SVR Stage 6 grid | -0.00017 | [-0.00029, -0.00005] | 23 / 9 | 0.016 | 0.047 | yes |
| random | MLP wider grid − MLP Stage 6 grid | -0.00083 | [-0.00109, -0.00055] | 29 / 3 | 4.6e-06 | 2.8e-05 | yes |
| random | SVR (wider grid) − MLP (wider grid) | -0.00252 | [-0.00373, -0.00159] | 31 / 1 | 2.3e-09 | 2.3e-08 | yes |

### Wider-grid hyperparameters
| predictor | hyperparameter | share at a grid edge | most frequent best value |
|---|---|---|---|
| SVR | C (grid 0.1 … 1000) | 22% | 1 (46%) |
| SVR | gamma (grid 0.03 … 10) | 20% | 1 (26%) |
| SVR | epsilon (grid 0.001 … 0.1) | 78% | 0.1 (56%) |
| MLP | hidden_layer_sizes (grid [16, 16] … [128, 128]) | 48% | [128, 128] (34%) |
| MLP | alpha (grid 1e-05 … 0.1) | 48% | 0.1 (28%) |
| MLP | activation (grid tanh … relu) | 100% | "relu" (85%) |

MLP activation chosen: relu 85%, tanh 15%.

Fit time including the wider grid search (mean over fits, 8 jobs in parallel, indicative): SVR 20.0 s, MLP 61.4 s.

