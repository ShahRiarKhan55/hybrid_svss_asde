# Stage 0 – Audit

Folder: `revision/` (sibling of `predictor_comparison/`; assumption: you meant the repo root). Script: `stage0_audit.py`. Evidence: `stage0/` (`audit_output.txt`, `code_lines.csv`, `budget_check.csv`, `asde_window_modes.csv`, `leakage_test.csv`). No existing file was changed. Line numbers refer to cell 0 of the named notebook.

## Summary

| # | Item | Result |
|---|---|---|
| 1 | AL = lowest predictive variance; AL-max; random 1% seed; refit rule | PASS on the rule, FAIL on the seed and the refit rule as you described them |
| 2 | SVTS alpha=0.3, beta=0.7; combination and normalisation | FAIL on the weights (code is 0.7 / 0.3); the similarity term is degenerate |
| 3 | SVSS alpha=beta=0.5, features [mu, sigma^2], refit per 1% batch, 3% | PASS; real GPR fits = 2 |
| 4 | ASDE K=2, tau=2.0, per-anchor windows; flooring; random fill | PASS against the paper text; the flooring is not a difference; 27.5% random fill confirmed |
| 5 | SDE fixed (2,2), 5x5 window; GPR fits | PASS; GPR fits in selection = 0 |
| 6 | Same number of unique, active, in-range sites | PASS unique/active/in range, FAIL on equal count (315–317) |
| 7 | Leakage | PASS (perturbation test) |
| 8 | GPR identical across methods | PASS identical; FAIL as "tuned": no bounds, no restarts |
| 9 | Timing | FAIL for final timings: one run, not idle, no median of 3 |

## 1. AL
| Check | Evidence | Result |
|---|---|---|
| Picks the lowest predictive variance | `raw_al` L49 `remaining[np.argmin(var.ravel())]` | PASS |
| Variance includes the noise term | `m.predict` returns variance with noise: difference to the noise-free variance equals the fitted noise variance (0.9601 in both) | as designed |
| Random 1% seed | `raw_al` L41 `train_idx = [rng.choice(idx_all)]`: one random site, not 1% | FAIL vs your description |
| Refit rule | `raw_al` L47 refit inside the loop after every added site: 314–316 fits per path, mean 315.9 (measured) | per site, not per 1% batch |
| Budget | L38 `int(n*0.10)` = 315–317 | floor, matches the paper's T |
| AL-max | not in the old notebooks; exists only in my `practical_cost` run | to be implemented as a separate baseline |

Another AL exists inside `hybrid_al_sde.ipynb` (L42–49): 1% random start, batches of 1%, `np.argsort(var)[:pick]` (also lowest variance). `point_find.ipynb` uses the highest variance. So the code base has three AL variants. The paper numbers for the AL baseline come from `raw_al`.

## 2. SVTS
| Check | Evidence | Result |
|---|---|---|
| Weights | `svts` L18 `ALPHA, BETA = 0.7, 0.3`; L70 `score = ALPHA*sim_n + BETA*var_n`. So 0.7 is on similarity and 0.3 on variance | FAIL vs your alpha=0.3, beta=0.7 (or your alpha/beta are the reverse roles) |
| Combination | features [mu, sigma^2] z-scored jointly over candidates and training sites; similarity = mean-centred cosine to the training centroid; `sim` and variance min-max normalised over the candidates; pick = argmax (L62–71) | as coded |
| Similarity term | Numeric check, path 1, 32 random sites: `feat_r_c` has rank 1 (every row is a multiple of (1,-1)); \|sim\| = 1 for 100% of candidates; after min-max, `sim_n` is exactly 0 for 50.0% and 1 for 50.0% of candidates | WARN: the "similarity" is a binary sign split |

Consequence: with weights 0.7/0.3 (and with SVSS 0.5/0.5) every candidate in the sign-agreeing half scores higher than (SVSS: at least as high as) every candidate in the other half, so the rule is "highest variance within the favoured half". The same code is used in SVSS.

## 3. SVSS
| Check | Evidence | Result |
|---|---|---|
| alpha = beta = 0.5 | `hybrid_svss_asde` L16–17, L87 | PASS |
| Features [mu, sigma^2], joint z-score, mean-centred cosine, min-max | L60–87, matches `Algorithms.txt` | PASS |
| Batch | L56 `step = round(n*0.01)` = 32 (paper text: floor, 31 for n = 3173); L57 `target = round(n*0.03)` = 95; L89 top-`pick` by score | PASS, minor floor/round difference |
| GPR fits | start = 32 random sites (no fit); fit at 32 adds 32; fit at 64 adds 31 → 95. Measured: 2.0 fits per path on all 32 paths | 2 fits, not 3 |

## 4. ASDE
| Check | Evidence | Result |
|---|---|---|
| K = 2, tau = 2.0 | L18–19 | PASS |
| Per-anchor window | L109 `ceil(mean dx) + bump`, bump = 1 if std of the 2 neighbours' values > tau; matches `Algorithms.txt` | PASS |
| Window of new sites | L110 `default_ab = int(np.mean(...))` floors. Paper text: `ᾱ ← mean({α_j})`, real-valued, and the acceptance test is `\|x_i − x_j\| ≤ α_j`. With integer coordinates a real window 2.33 and a floored window 2 give the same test | The flooring is NOT a difference from the paper: floor and real-valued mean accept exactly the same sites (paths 1–4, `asde_window_modes.csv`) |
| Measured windows | anchors' mean window (2.33, 2.42); new-site window (2.00, 1.97) = SDE's fixed (2,2). Rounding differs on 3 of 4 checked paths, ceiling on all 4 (path 1: 134 accepted additions with floor, 76 with ceil (3,3)) | the new-site window is effectively SDE's |
| Random fill | additions accepted by the rule: 134.8 of 222 (60.7%); random fill 87.2 = 27.5% of all hybrid sites (anchors 30.0%, accepted 42.5%). The paper text also fills with a shuffle of the remaining sites | confirmed |

## 5. SDE
| Check | Evidence | Result |
|---|---|---|
| Window | `sde` L11 `ALPHA, BETA = 2, 2` hard-coded; L40 reject when \|dx\| ≤ 2 and \|dy\| ≤ 2 (5×5); no threshold is derived from calibration (only the site positions come from FPGA 1) | PASS |
| Start and fill | L33 shuffle of all sites with `Random(42)`, first site starts; short targets are filled with rejected sites. 84.5% of sites pass the test, 15.5% are fill | as coded |
| GPR fits in selection | none in `sde_select`; GPR appears only in `run_path` (L73–74) for the final prediction. Measured: 0 fits | 0 |
| Budget | L12 `round(0.10*n)` | round |

## 6. Budget
| Method | sites per path | unique | active | in range |
|---|---|---|---|---|
| AL | 315–317 | yes | yes | yes |
| SVTS | 315–317 | yes | yes | yes |
| SDE | 316–317 | yes | yes | yes |
| hybrid | 316–317 | yes | yes | yes |

- Active sites per path: 3,155 on path 1, 3,173 on the other 31 paths.
- AL and SVTS use `int(0.10 n)`, SDE `round(0.10 n)`, hybrid `round(0.03 n) + round(0.07 n)`, so counts differ by up to 2 sites.
- On the other FPGAs the 3-sigma filter removes on average 0.04–0.07 selected sites.
- Assumption for the revision: fixed counts for every path, 3% / 5% / 7% / 10% = 95 / 158 / 222 / 317 sites (floor of the fraction × 3,173, as in the paper text), including the seed points.

## 7. Leakage
| Check | Result |
|---|---|
| Selection reads only FPGA 1 | `hybrid_svss_asde` L142 `load_dat(CALIB_FPGA, path)`; the other notebooks use `CALIB_FPGA = 1` the same way. PASS |
| No candidate value used | Perturbation test, path 1: add N(0, 5 sd) noise to the values of every non-selected FPGA-1 site and re-run selection. Order identical for hybrid, SVTS, AL and SDE (`leakage_test.csv`). PASS |
| Test sites in fitting | old runs: test = sites outside the training set, never fitted. PASS |
| 3-sigma filter and d_spec | computed from the full path, including test sites (offline, expected) |
| Aggregation order | test sites → path mean → FPGA mean → mean over FPGAs (`final_gpr`/`practical_cost` analysis: `groupby("FPGA")[col].mean().mean()` over per-path means). FPGA 1 excluded with `FPGA >= 2`. PASS |

New protocol note: the old runs let selection pick from all sites. The new fixed 30% held-out set shrinks the selectable pool to about 2,221 sites; 10% of 3,173 is then 14.3% of the pool.

## 8. GPR
| Check | Evidence | Result |
|---|---|---|
| Kernel, normalizer | `RBF(1, active_dims=[0]) + RBF(1, active_dims=[1])`, `normalizer=True` everywhere (`practical_cost`, `final_gpr`, `strategies`, `budget`, `why_sde`) | PASS |
| Init | lengthscales 1.0, variances 1.0, noise 1.0 (raw coordinates) | identical |
| Bounds, restarts | only GPy's positivity constraints; no bounds; single start; no restarts | FAIL vs "tuned" |
| Optimiser setting | selection loops: default optimiser (up to 1000 iterations). Predictor: `max_iters=200`. Same training set (path 1, hybrid): identical fit in both (converged in 19 iterations, LL −364.593, lengthscales 0.857 / 4.308) | different settings, same result here |
| Carry-over | a new `GPRegression` per fit, defaults every time; no warm start | PASS |
| Extreme fits | `final_gpr`: 4.2% of additive fits end with a lengthscale below 0.01 or above 100 training-std units | motivates bounds |

## 9. Timing
| Check | Evidence | Result |
|---|---|---|
| How `practical_cost` timed | `time.perf_counter()` and `time.process_time()` around the selection function (fits and ASDE included, data loading excluded); one process, sequential, one run per method and path | single run, no repeats |
| Threads | env vars set before importing NumPy; checked now: OpenBLAS 1, OpenMP 1 | PASS |
| Idle machine | not idle: Balanced power plan, other processes running (AWCC, Killer, the Claude app, explorer); CPU clock read 2.3 GHz (base), the AL-max run slowed from 40 s to over 100 s per path | FAIL for final timing |
| Repeat spread | hybrid selection path 1, three repeats now: 2.79, 2.94, 2.88 s (5% spread) | |
| For the final timings | OMP/MKL threads = 1, sequential, 3 repeats, median, High-performance power plan and other apps closed. I can set the first three; the power plan and closing apps are yours | |

## Differences the revision protocol will introduce (all numbers will change)
- Fixed 30% held-out test set that no method can select; old runs used all non-selected sites.
- TUNE (FPGAs 2–26) / TEST (27–50) split; primary results exclude FPGA 1.
- Fixed budgets 95 / 158 / 222 / 317 for every method and path.
- 5 seeds (random seed set, random fill, random baseline); the old runs had one fixed seed.
- One GPR function with bounds, used both inside selection and for prediction.

## Decisions I need before Stage 1 (defaults if you do not answer)
1. SVTS weights: the code is 0.7 on similarity and 0.3 on variance. Default: run the code's setting (0.7, 0.3) as "SVTS" and your setting (0.3, 0.7) as "SVTS-var", and report both.
2. AL seed and refit: default is a random 1% seed shared with the hybrid (same seed set per seed index), then one site per refit (lowest variance for AL, highest for AL-max). Alternative: the old single seed site.
3. GPR bounds and restarts: default is bounds only (lengthscale 0.5–300 cells, variances at least 1e-4), no restarts, because restarts would multiply every selection cost (AL: 316 fits per path) and a 5-restart test earlier changed nothing.
4. The ASDE window: floor, float and rounding give the same result only for floor and float. Stage 3c will test round and ceil as real changes, not as a bug fix.
5. Similarity term: the mean-centred cosine on two features is binary. Default: keep it unchanged in the baseline hybrid; Stage 3 may add a variant that uses a meaningful similarity (for example the cosine without centring) if you approve.

Stopped here, as requested.
