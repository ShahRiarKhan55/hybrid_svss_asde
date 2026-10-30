# Stage 6 – wider-grid check (TEST FPGAs 27 and 29 only)

## Protocol
Same sites (hybrid and Random, 5 seeds, 317 sites), held-out test set, error definition, 5-fold CV on the training sites (mean absolute error), standardised inputs/target, FPGAs 27 and 29, all 32 paths. Wider grids: SVR with C ∈ {0.1, 1, 10, 100, 1000}, gamma ∈ {0.03, 0.1, 0.3, 1, 3, 10}, epsilon ∈ {0.001, 0.01, 0.1}; MLP (Adam, up to 1000 iterations) with hidden widths (16,16), (32,32), (64,64), (128,128), alpha ∈ {1e-5, 1e-4, 1e-3, 1e-2, 1e-1} and activation ∈ {tanh, relu} chosen by CV. The "Stage 6 grid" rows and the additive GPR are the Stage 6 results for the same two FPGAs. Errors are the mean over the 2 FPGAs; SD over the 5 seeds; paired differences are path-level (32 paths, mean over seeds and the 2 FPGAs) with a 95% bootstrap CI (10,000 resamples) and Wilcoxon, Holm over the 10 comparisons below. Only 2 devices were used, so the numbers differ from the 24-FPGA table.

## Error
| predictor | hybrid sites | rank among the 3 families* | Random sites | rank* |
|---|---|---|---|---|
| Additive RBF GPR | 0.0761 ± 0.0002 | 1 | 0.0756 ± 0.0004 | 1 |
| SVR, Stage 6 grid | 0.0762 ± 0.0003 | 2 | 0.0764 ± 0.0003 | 2 |
| SVR, wider grid | 0.0761 ± 0.0004 | 1 | 0.0762 ± 0.0003 | 2 |
| MLP, Stage 6 grid (tanh) | 0.0797 ± 0.0002 | 3 | 0.0796 ± 0.0002 | 3 |
| MLP, wider grid (tanh or relu) | 0.0788 ± 0.0003 | 3 | 0.0787 ± 0.0002 | 3 |

\* rank among the three families (additive GPR, SVR, MLP), using the Stage 6 grids for the SVR and MLP rows marked "Stage 6 grid" and the wider grids for the rows marked "wider grid".

Order of the three families (lowest error first): hybrid sites — Stage 6 grids: Additive GPR < SVR < MLP; wider grids: SVR < Additive GPR < MLP. Random sites — Stage 6 grids: Additive GPR < SVR < MLP; wider grids: Additive GPR < SVR < MLP.

## Paired differences
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

## Wider-grid hyperparameters
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
