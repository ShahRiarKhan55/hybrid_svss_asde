# Stage 5 – Cost (final designs)

## Protocol
- One Python process, BLAS/OpenMP threads = 1 (set before importing NumPy; checked in Stage 0), sequential, nothing else started by me during the run. Selection on FPGA 1, all 32 paths, seed 1, held-out pool as in Stages 1–4; data loading and the held-out set are outside the timer. Each (path, method) was timed 3 times (repeats interleaved: repeat 1 over all paths, then repeat 2, then repeat 3); the table uses the median of the 3 repeats per path, then the mean over the 32 paths.
- Same additive bounded RBF GPR (`lib.fit`, up to 200 optimiser iterations) inside every selection loop and for the final fit. Final fit and prediction: hybrid training set (317 sites) on FPGAs 2, 14 and 26 (TUNE FPGAs, timing only) × 32 paths, 3 repeats, median; prediction on all non-training sites of the FPGA.
- Machine: i7-12700H (14 cores / 20 threads), 31.7 GB RAM, Windows 11, Python 3.14.7, NumPy 2.3.5, GPy 1.14.2.
- Power plan during the run: Balanced. I did not change it (the machine was not switched to High performance). Windows reported processor frequency 1776–2259 MHz (mean 1958), 41–91% of maximum frequency, whole-machine CPU utilisation 4–67% (mean 13%) at the start of each path block. Background services (vendor utilities, the Claude app) were running, so the machine was not idle.

## Selection cost
| method | GPR fits per path | selection time per path, s (median of 3 repeats; mean over 32 paths) | median over paths | min – max over paths | CPU time, s | relative to hybrid | mean repeat-to-repeat spread |
|---|---|---|---|---|---|---|---|
| AL | 285.0 | 45.033 | 42.362 | 34.275 – 62.137 | 44.758 | 234.4x | 36% |
| AL-max | 285.0 | 44.903 | 45.004 | 36.303 – 53.717 | 44.583 | 233.8x | 33% |
| SVTS (0.3 sim, 0.7 var) | 285.0 | 46.322 | 45.558 | 38.996 – 58.544 | 46.021 | 241.1x | 29% |
| SDE | 0.0 | 0.042 | 0.041 | 0.039 – 0.077 | 0.043 | 0.2x | 24% |
| Random | 0.0 | 0.000 | 0.000 | 0.000 – 0.000 | 0.000 | 0.0x | 32% |
| Maximin | 0.0 | 0.005 | 0.005 | 0.004 – 0.007 | 0.002 | 0.0x | 24% |
| Hybrid (SVSS 95 + ASDE 222) | 2.0 | 0.192 | 0.189 | 0.157 – 0.308 | 0.191 | 1.0x | 31% |

## Final fit and prediction (317 training sites)
| quantity | mean over 96 (FPGA, path) | min – max |
|---|---|---|
| final GPR fit, s | 0.323 | 0.221 – 0.608 |
| prediction on 2856 sites, s | 0.055 | 0.053 – 0.059 |

The final fit and prediction are the same for every method (same GPR, same number of training sites).

## Notes
- Times are wall-clock seconds on a laptop; the repeat-to-repeat spread column shows how reproducible they were. Absolute values depend on the machine; the ratios and the fit counts are the portable numbers.
- The hybrid's 2 fits are on 32 and 64 sites; the AL, AL-max and SVTS loops refit after every added site (285 fits, on 32 to 316 sites).
- SDE, Random and Maximin use no GPR; their times are for the Python selection loops only.
