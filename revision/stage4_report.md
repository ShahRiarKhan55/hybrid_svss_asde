# Stage 4 – Final TEST run (FPGAs 27–50, run once)

## Protocol
- Frozen design: the current hybrid (SVSS 95 anchors, alpha = beta = 0.5, then ASDE with random fill), unchanged; chosen by the Stage 3 rule (no candidate qualified). SVTS weights (0.3 similarity, 0.7 variance) were chosen on TUNE in Stage 1. The two exploratory rows (anchors 1%, maximin anchors) were defined in `stage3_rule_section.md` before any result.
- Selection orders are the saved FPGA-1 orders of Stages 1–3 (5 seeds, same sites on all FPGAs); selection never used TEST data. Evaluation: additive bounded RBF GPR (same function as in selection), held-out 30% test set (about 952 test sites per path; never selectable, never fitted), error = |mu − v| / d_spec, test sites → path → FPGA → mean over the 24 TEST FPGAs. Mean ± SD over 5 seeds of that FPGA-averaged error. FPGA 1 is excluded from the primary tables.
- Paired differences: path-level values (mean over the 5 seeds and 24 TEST FPGAs), 32 paths, 95% percentile bootstrap CI (10,000 resamples, seed 0), Wilcoxon signed-rank (two-sided, zero differences dropped) on the 32 path-level differences. Holm correction over all 46 comparisons reported in this stage (24 headline, 8 exploratory, 14 ablation) as one family; the headline-only Holm p-values are in `stage4/tests.csv`.
- Budgets 95 / 158 / 222 / 317 sites including seed points. Noise-floor references: constant predictor = mean of a random 317-site training set; "GPR 50%" = the same GPR trained on a random 1587 sites (50% of the 3,173 sites, taken from the selectable 70% pool), both on the same held-out test sets. The GPR 50% reference uses one seed only (each fit takes about 20 s).
- Training sites per TEST FPGA: the mean equals the nominal budget (95 / 158 / 222 / 317; 1,587 for the 50% reference) to one decimal, so the 3-sigma filter removed almost none of the selected sites.

## Headline table (TEST 27–50, held-out test set, mean error ± SD over 5 seeds)
| method | 95 | 158 | 222 | 317 |
|---|---|---|---|---|
| Hybrid (frozen) | 0.0867 ± 0.0006 | 0.0798 ± 0.0006 | 0.0774 ± 0.0003 | 0.0758 ± 0.0002 |
| AL | 0.0897 ± 0.0009 | 0.0881 ± 0.0005 | 0.0867 ± 0.0006 | 0.0855 ± 0.0008 |
| AL-max | 0.0840 ± 0.0004 | 0.0791 ± 0.0005 | 0.0773 ± 0.0003 | 0.0755 ± 0.0003 |
| SVTS | 0.0836 ± 0.0011 | 0.0799 ± 0.0005 | 0.0782 ± 0.0004 | 0.0760 ± 0.0005 |
| SDE | 0.0820 ± 0.0005 | 0.0782 ± 0.0003 | 0.0766 ± 0.0004 | 0.0753 ± 0.0002 |
| Random | 0.0823 ± 0.0007 | 0.0788 ± 0.0005 | 0.0769 ± 0.0005 | 0.0752 ± 0.0004 |
| Maximin | 0.0823 ± 0.0006 | 0.0787 ± 0.0002 | 0.0769 ± 0.0004 | 0.0754 ± 0.0004 |

Exploratory rows (not part of the headline comparison):

| method | 95 | 158 | 222 | 317 |
|---|---|---|---|---|
| anchors 1% (exploratory) | 0.0818 ± 0.0010 | 0.0790 ± 0.0005 | 0.0772 ± 0.0005 | 0.0756 ± 0.0003 |
| maximin anchors (exploratory) | 0.0823 ± 0.0006 | 0.0789 ± 0.0004 | 0.0772 ± 0.0004 | 0.0756 ± 0.0005 |

References: constant predictor 0.0959 ± 0.0002; GPR trained on 50% of the sites 0.0715 (seed 1 only, one fit of 1587 sites per path); headroom between them 0.0245.

## First budget below error 0.09 and 0.085
| method | first budget below 0.09 | first budget below 0.085 |
|---|---|---|
| Hybrid (frozen) | 95 | 158 |
| AL | 95 | -- |
| AL-max | 95 | 95 |
| SVTS | 95 | 95 |
| SDE | 95 | 95 |
| Random | 95 | 95 |
| Maximin | 95 | 95 |
| anchors 1% (exploratory) | 95 | 95 |
| maximin anchors (exploratory) | 95 | 95 |

## Paired differences: hybrid against each baseline (headline)
| comparison (negative = first lower error) | mean diff | 95% CI | first wins / second wins (of 32) | Wilcoxon p | Holm p (all 46 tests) | significant after Holm |
|---|---|---|---|---|---|---|
| Hybrid − AL @ 95 | -0.00301 | [-0.00387, -0.00223] | 28 / 4 | 3.3e-08 | 1.3e-06 | yes |
| Hybrid − AL-max @ 95 | +0.00263 | [+0.00156, +0.00373] | 6 / 26 | 3.5e-05 | 0.0012 | yes |
| Hybrid − SVTS @ 95 | +0.00308 | [+0.00187, +0.00429] | 7 / 25 | 2.8e-05 | 0.00099 | yes |
| Hybrid − SDE @ 95 | +0.00469 | [+0.00335, +0.00603] | 4 / 28 | 1.4e-07 | 5.7e-06 | yes |
| Hybrid − Random @ 95 | +0.00431 | [+0.00306, +0.00559] | 2 / 30 | 3e-07 | 1.2e-05 | yes |
| Hybrid − Maximin @ 95 | +0.00433 | [+0.00316, +0.00546] | 4 / 28 | 3.5e-07 | 1.3e-05 | yes |
| Hybrid − AL @ 158 | -0.00825 | [-0.00943, -0.00700] | 31 / 1 | 1.4e-09 | 6e-08 | yes |
| Hybrid − AL-max @ 158 | +0.00069 | [-0.00014, +0.00151] | 15 / 17 | 0.18 | 1 | no |
| Hybrid − SVTS @ 158 | -0.00005 | [-0.00081, +0.00075] | 19 / 13 | 0.66 | 1 | no |
| Hybrid − SDE @ 158 | +0.00162 | [+0.00096, +0.00233] | 7 / 25 | 4.9e-05 | 0.0016 | yes |
| Hybrid − Random @ 158 | +0.00101 | [+0.00046, +0.00156] | 7 / 25 | 0.0011 | 0.032 | yes |
| Hybrid − Maximin @ 158 | +0.00107 | [+0.00055, +0.00156] | 9 / 23 | 0.0004 | 0.012 | yes |
| Hybrid − AL @ 222 | -0.00935 | [-0.01067, -0.00800] | 31 / 1 | 9.3e-10 | 4.3e-08 | yes |
| Hybrid − AL-max @ 222 | +0.00002 | [-0.00061, +0.00067] | 14 / 18 | 0.96 | 1 | no |
| Hybrid − SVTS @ 222 | -0.00080 | [-0.00157, -0.00004] | 20 / 12 | 0.054 | 1 | no |
| Hybrid − SDE @ 222 | +0.00072 | [+0.00029, +0.00115] | 7 / 25 | 0.0014 | 0.042 | yes |
| Hybrid − Random @ 222 | +0.00048 | [+0.00019, +0.00078] | 9 / 23 | 0.005 | 0.13 | no |
| Hybrid − Maximin @ 222 | +0.00042 | [-0.00007, +0.00086] | 11 / 21 | 0.031 | 0.62 | no |
| Hybrid − AL @ 317 | -0.00968 | [-0.01114, -0.00816] | 31 / 1 | 9.3e-10 | 4.3e-08 | yes |
| Hybrid − AL-max @ 317 | +0.00035 | [-0.00014, +0.00087] | 14 / 18 | 0.22 | 1 | no |
| Hybrid − SVTS @ 317 | -0.00020 | [-0.00079, +0.00043] | 19 / 13 | 0.31 | 1 | no |
| Hybrid − SDE @ 317 | +0.00050 | [+0.00021, +0.00080] | 8 / 24 | 0.0029 | 0.077 | no |
| Hybrid − Random @ 317 | +0.00059 | [+0.00026, +0.00094] | 11 / 21 | 0.0024 | 0.067 | no |
| Hybrid − Maximin @ 317 | +0.00042 | [+0.00006, +0.00075] | 10 / 22 | 0.0083 | 0.19 | no |

## Paired differences: hybrid against the exploratory designs
| comparison (negative = first lower error) | mean diff | 95% CI | first wins / second wins (of 32) | Wilcoxon p | Holm p (all 46 tests) | significant after Holm |
|---|---|---|---|---|---|---|
| Hybrid − anchors 1% (exploratory) @ 95 | +0.00483 | [+0.00381, +0.00587] | 2 / 30 | 2.3e-09 | 9.8e-08 | yes |
| Hybrid − maximin anchors (exploratory) @ 95 | +0.00433 | [+0.00318, +0.00546] | 4 / 28 | 3.5e-07 | 1.3e-05 | yes |
| Hybrid − anchors 1% (exploratory) @ 158 | +0.00080 | [+0.00016, +0.00144] | 10 / 22 | 0.018 | 0.41 | no |
| Hybrid − maximin anchors (exploratory) @ 158 | +0.00094 | [+0.00036, +0.00152] | 9 / 23 | 0.0078 | 0.19 | no |
| Hybrid − anchors 1% (exploratory) @ 222 | +0.00017 | [-0.00037, +0.00066] | 14 / 18 | 0.15 | 1 | no |
| Hybrid − maximin anchors (exploratory) @ 222 | +0.00020 | [-0.00026, +0.00062] | 12 / 20 | 0.12 | 1 | no |
| Hybrid − anchors 1% (exploratory) @ 317 | +0.00019 | [-0.00017, +0.00052] | 10 / 22 | 0.084 | 1 | no |
| Hybrid − maximin anchors (exploratory) @ 317 | +0.00025 | [-0.00009, +0.00055] | 12 / 20 | 0.071 | 1 | no |

## Paths where each method is the best of the seven headline methods (of 32 paths)
| method | 95 | 158 | 222 | 317 |
|---|---|---|---|---|
| Hybrid (frozen) | 0 | 1 | 0 | 1 |
| AL | 0 | 0 | 0 | 0 |
| AL-max | 5 | 7 | 5 | 6 |
| SVTS | 2 | 1 | 5 | 4 |
| SDE | 9 | 13 | 7 | 4 |
| Random | 5 | 7 | 7 | 11 |
| Maximin | 11 | 3 | 8 | 6 |

## Noise floor
Share of the headroom (constant − GPR 50%) closed by each method:

| method | share of headroom closed @ 95 | share of headroom closed @ 158 | share of headroom closed @ 222 | share of headroom closed @ 317 |
|---|---|---|---|---|
| Hybrid (frozen) | 37.9% | 65.9% | 75.9% | 82.2% |
| AL | 25.6% | 32.2% | 37.7% | 42.6% |
| AL-max | 48.7% | 68.7% | 76.0% | 83.7% |
| SVTS | 50.5% | 65.7% | 72.7% | 81.4% |
| SDE | 57.1% | 72.5% | 78.9% | 84.3% |
| Random | 55.5% | 70.0% | 77.9% | 84.6% |
| Maximin | 55.6% | 70.3% | 77.6% | 83.9% |
| anchors 1% (exploratory) | 57.7% | 69.2% | 76.6% | 83.0% |
| maximin anchors (exploratory) | 55.6% | 69.7% | 76.7% | 83.2% |

The hybrid − random differences relative to the headroom: +17.6% @ 95, +4.1% @ 158, +2.0% @ 222, +2.4% @ 317.

## Ablation on TEST (Stage 2 repeated, 317 sites)
| stage 1 (95 sites) \ stage 2 (222 sites) | random fill | fixed SDE (2,2) | adaptive ASDE |
|---|---|---|---|
| random | 0.0752 ± 0.0004 | 0.0753 ± 0.0003 | 0.0755 ± 0.0003 |
| SVSS | 0.0756 ± 0.0003 | 0.0760 ± 0.0003 | 0.0758 ± 0.0002 (hybrid) |

| cell @ 317 sites | error (mean ± SD over seeds) |
|---|---|
| SVSS-only | 0.0785 ± 0.0007 |
| AL-max (1% seed) | 0.0755 ± 0.0003 |
| ASDE-only(rand seed) | 0.0756 ± 0.0003 |
| ASDE-only(maximin seed) | 0.0756 ± 0.0003 |
| SDE ref | 0.0753 ± 0.0002 |

| comparison (negative = first lower error) | mean diff | 95% CI | first wins / second wins (of 32) | Wilcoxon p | Holm p (all 46 tests) | significant after Holm |
|---|---|---|---|---|---|---|
| Hybrid − rand+SDE @ 317 | +0.00048 | [+0.00017, +0.00080] | 7 / 25 | 0.0026 | 0.069 | no |
| Hybrid − rand+ASDE @ 317 | +0.00036 | [+0.00004, +0.00065] | 9 / 23 | 0.022 | 0.45 | no |
| Hybrid − SVSS+rand @ 317 | +0.00023 | [-0.00011, +0.00059] | 15 / 17 | 0.43 | 1 | no |
| Hybrid − SVSS+SDE @ 317 | -0.00016 | [-0.00052, +0.00021] | 22 / 10 | 0.17 | 1 | no |
| Hybrid − SVSS-only @ 317 | -0.00272 | [-0.00345, -0.00193] | 28 / 4 | 1.1e-06 | 4e-05 | yes |
| Hybrid − ASDE-only(rand seed) @ 317 | +0.00019 | [-0.00018, +0.00052] | 10 / 22 | 0.084 | 1 | no |
| Hybrid − ASDE-only(maximin seed) @ 317 | +0.00023 | [-0.00011, +0.00057] | 10 / 22 | 0.12 | 1 | no |
| AL-max − SVSS-only @ 317 | -0.00307 | [-0.00376, -0.00240] | 31 / 1 | 9.3e-10 | 4.3e-08 | yes |
| ASDE-only(maximin seed) − ASDE-only(rand seed) @ 317 | -0.00004 | [-0.00022, +0.00015] | 17 / 15 | 0.56 | 1 | no |
| Stage-1 main effect: SVSS − random (mean of 3 stage-2 choices) @ 317 | +0.00045 | [+0.00028, +0.00062] | 4 / 28 | 0.0001 | 0.0032 | yes |
| Stage-2 main effect: ASDE − fixed SDE @ 317 | -0.00002 | [-0.00025, +0.00022] | 18 / 14 | 0.82 | 1 | no |
| Stage-2 main effect: ASDE − random fill @ 317 | +0.00023 | [-0.00007, +0.00055] | 13 / 19 | 0.38 | 1 | no |
| Stage-2 main effect: fixed SDE − random fill @ 317 | +0.00025 | [-0.00011, +0.00060] | 16 / 16 | 0.37 | 1 | no |
| Interaction: (SVSS − random) under ASDE minus under fixed SDE @ 317 | -0.00029 | [-0.00076, +0.00019] | 19 / 13 | 0.28 | 1 | no |

## Secondary: TEST + FPGA 1 (headline and exploratory rows)
| method | 95 | 158 | 222 | 317 |
|---|---|---|---|---|
| Hybrid (frozen) | 0.0871 ± 0.0006 | 0.0803 ± 0.0006 | 0.0778 ± 0.0003 | 0.0763 ± 0.0002 |
| AL | 0.0901 ± 0.0009 | 0.0885 ± 0.0005 | 0.0872 ± 0.0006 | 0.0860 ± 0.0008 |
| AL-max | 0.0845 ± 0.0004 | 0.0796 ± 0.0005 | 0.0778 ± 0.0003 | 0.0759 ± 0.0003 |
| SVTS | 0.0840 ± 0.0011 | 0.0803 ± 0.0005 | 0.0786 ± 0.0004 | 0.0765 ± 0.0005 |
| SDE | 0.0824 ± 0.0005 | 0.0786 ± 0.0003 | 0.0771 ± 0.0004 | 0.0758 ± 0.0002 |
| Random | 0.0828 ± 0.0007 | 0.0792 ± 0.0005 | 0.0773 ± 0.0005 | 0.0757 ± 0.0003 |
| Maximin | 0.0828 ± 0.0006 | 0.0792 ± 0.0002 | 0.0774 ± 0.0004 | 0.0758 ± 0.0003 |
| anchors 1% (exploratory) | 0.0823 ± 0.0010 | 0.0794 ± 0.0005 | 0.0776 ± 0.0005 | 0.0761 ± 0.0003 |
| maximin anchors (exploratory) | 0.0828 ± 0.0006 | 0.0793 ± 0.0004 | 0.0776 ± 0.0004 | 0.0760 ± 0.0005 |

## Secondary: paper's old test set (all non-training sites), TEST 27–50
| method | 95 | 158 | 222 | 317 |
|---|---|---|---|---|
| Hybrid (frozen) | 0.0874 ± 0.0005 | 0.0805 ± 0.0005 | 0.0781 ± 0.0002 | 0.0764 ± 0.0002 |
| AL | 0.0904 ± 0.0009 | 0.0892 ± 0.0005 | 0.0881 ± 0.0007 | 0.0872 ± 0.0009 |
| AL-max | 0.0846 ± 0.0005 | 0.0797 ± 0.0005 | 0.0780 ± 0.0003 | 0.0763 ± 0.0002 |
| SVTS | 0.0842 ± 0.0010 | 0.0807 ± 0.0005 | 0.0791 ± 0.0004 | 0.0771 ± 0.0005 |
| SDE | 0.0824 ± 0.0005 | 0.0786 ± 0.0002 | 0.0771 ± 0.0004 | 0.0757 ± 0.0003 |
| Random | 0.0828 ± 0.0007 | 0.0792 ± 0.0005 | 0.0773 ± 0.0005 | 0.0757 ± 0.0003 |
| Maximin | 0.0827 ± 0.0005 | 0.0792 ± 0.0002 | 0.0774 ± 0.0004 | 0.0759 ± 0.0003 |
| anchors 1% (exploratory) | 0.0823 ± 0.0009 | 0.0795 ± 0.0004 | 0.0777 ± 0.0004 | 0.0761 ± 0.0003 |
| maximin anchors (exploratory) | 0.0827 ± 0.0005 | 0.0793 ± 0.0003 | 0.0776 ± 0.0004 | 0.0761 ± 0.0004 |

## Hyperparameters at a bound (share of predictor fits, within 5% of the bound)
| method | column lengthscale lower / upper | row lengthscale lower / upper | column variance lower / upper | row variance lower / upper | noise lower / upper |
|---|---|---|---|---|---|
| AL | 21.5% / 0.0% | 5.8% / 0.1% | 1.9% / 0.0% | 4.8% / 0.0% | 0.0% / 0.0% |
| AL-max | 3.5% / 0.0% | 7.4% / 0.0% | 0.6% / 0.0% | 1.7% / 0.0% | 0.5% / 0.0% |
| ASDE-only(maximin seed) | 0.7% / 0.0% | 3.0% / 0.0% | 0.0% / 0.0% | 1.4% / 0.0% | 0.0% / 0.0% |
| ASDE-only(rand seed) | 0.9% / 0.0% | 2.9% / 0.0% | 0.0% / 0.0% | 1.1% / 0.0% | 0.0% / 0.0% |
| Hybrid (frozen) | 6.3% / 0.0% | 6.5% / 0.0% | 0.9% / 0.0% | 3.0% / 0.0% | 0.0% / 0.0% |
| Maximin | 4.2% / 0.0% | 8.2% / 0.0% | 0.2% / 0.0% | 4.1% / 0.0% | 0.0% / 0.0% |
| Random | 5.0% / 0.0% | 6.7% / 0.0% | 0.3% / 0.0% | 4.0% / 0.0% | 0.0% / 0.0% |
| SDE | 4.6% / 0.0% | 8.3% / 0.0% | 0.5% / 0.0% | 4.1% / 0.0% | 0.0% / 0.0% |
| SVSS+SDE | 1.1% / 0.0% | 3.3% / 0.0% | 0.0% / 0.0% | 1.1% / 0.0% | 0.0% / 0.0% |
| SVSS+rand | 1.1% / 0.0% | 2.8% / 0.0% | 0.0% / 0.0% | 1.0% / 0.0% | 0.0% / 0.0% |
| SVSS-only | 3.3% / 0.0% | 3.0% / 0.0% | 0.0% / 0.0% | 2.0% / 0.0% | 0.0% / 0.0% |
| SVTS | 2.8% / 0.0% | 8.3% / 0.0% | 0.3% / 0.0% | 1.8% / 0.0% | 0.2% / 0.0% |
| anchors 1% (exploratory) | 4.0% / 0.0% | 6.7% / 0.0% | 0.4% / 0.0% | 3.1% / 0.0% | 0.0% / 0.0% |
| gpr50 | 0.1% / 0.0% | 0.0% / 0.0% | 0.0% / 0.0% | 0.0% / 0.0% | 0.0% / 0.0% |
| maximin anchors (exploratory) | 3.3% / 0.0% | 8.3% / 0.0% | 0.1% / 0.0% | 3.2% / 0.0% | 0.0% / 0.0% |
| rand+ASDE | 1.1% / 0.0% | 1.9% / 0.0% | 0.0% / 0.0% | 1.2% / 0.0% | 0.0% / 0.0% |
| rand+SDE | 0.6% / 0.0% | 4.0% / 0.0% | 0.0% / 0.0% | 1.2% / 0.0% | 0.0% / 0.0% |

## Spread (descriptive, added after the main run; not part of the pre-specified tests)
Per-FPGA error (mean over paths and 5 seeds), then the spread across the 24 TEST FPGAs:

| method | SD across 24 FPGAs @ 95 | SD across 24 FPGAs @ 158 | SD across 24 FPGAs @ 222 | SD across 24 FPGAs @ 317 | min - max across FPGAs @ 317 |
|---|---|---|---|---|---|
| Hybrid (frozen) | 0.0043 | 0.0042 | 0.0041 | 0.0041 | 0.0638 - 0.0818 |
| AL | 0.0040 | 0.0041 | 0.0042 | 0.0041 | 0.0729 - 0.0911 |
| AL-max | 0.0047 | 0.0045 | 0.0043 | 0.0041 | 0.0633 - 0.0809 |
| SVTS | 0.0045 | 0.0045 | 0.0044 | 0.0042 | 0.0635 - 0.0816 |
| SDE | 0.0043 | 0.0042 | 0.0041 | 0.0041 | 0.0633 - 0.0811 |
| Random | 0.0042 | 0.0042 | 0.0042 | 0.0041 | 0.0632 - 0.0813 |
| Maximin | 0.0044 | 0.0043 | 0.0042 | 0.0041 | 0.0634 - 0.0814 |

Spread across paths (path-level error, mean over seeds and FPGAs):

| method | SD across 32 paths @ 95 | SD across 32 paths @ 158 | SD across 32 paths @ 222 | SD across 32 paths @ 317 |
|---|---|---|---|---|
| Hybrid (frozen) | 0.0221 | 0.0235 | 0.0238 | 0.0241 |
| AL | 0.0215 | 0.0216 | 0.0221 | 0.0222 |
| AL-max | 0.0227 | 0.0225 | 0.0226 | 0.0231 |
| SVTS | 0.0234 | 0.0229 | 0.0227 | 0.0231 |
| SDE | 0.0230 | 0.0234 | 0.0236 | 0.0238 |
| Random | 0.0231 | 0.0238 | 0.0240 | 0.0241 |
| Maximin | 0.0238 | 0.0239 | 0.0241 | 0.0239 |

Hybrid - random per FPGA (mean, min, max over the 24 FPGAs): 95: +0.00431 (+0.00332 to +0.00547); 158: +0.00101 (+0.00061 to +0.00147); 222: +0.00048 (+0.00018 to +0.00077); 317: +0.00059 (+0.00041 to +0.00073); the hybrid is worse than random on 24 of 24 FPGAs at 317 sites.

## Claims the data support and do not support (TEST FPGAs 27–50, held-out test set, Holm over 46 comparisons)

**Supported**
1. *The hybrid has lower error than AL (lowest variance) at every budget.* Differences −0.0030 (95 sites) to −0.0097 (317); the hybrid wins on 28–31 of 32 paths; Holm p ≤ 1.3e-6 at all four budgets.
2. *From 158 sites on, the hybrid is not distinguishable from SVTS (weights chosen on TUNE) or from AL-max.* Hybrid − SVTS: −0.00005 (158), −0.00080 (222), −0.00020 (317); hybrid − AL-max: +0.00069, +0.00002, +0.00035. All Holm p = 1. This is "no difference detected", not proof of equivalence.
3. *The hybrid needs far fewer GPR fits in selection.* 2 fits per path against 285 for AL, AL-max and SVTS (fit counts from Stage 1; the timing is Stage 5). SDE, Random and Maximin need no GPR fit.
4. *At 222 and 317 sites the hybrid is within 0.0004–0.0007 of Random, SDE and Maximin* (2–3% of the 0.0245 headroom between the constant predictor and a GPR on 50% of the sites). After Holm, none of these differences is significant at 222 and 317 sites except hybrid − SDE at 222 sites (+0.00072, Holm p = 0.042). At 317 sites the unadjusted CIs of hybrid − Random and hybrid − SDE exclude 0; the Holm p-values are 0.067 and 0.077.
5. *The similarity term of SVSS adds nothing.* SVSS-only is worse than the hybrid (+0.0027, Holm p = 4e-5) and AL-max (plain highest variance, same seed) beats SVSS-only by 0.0031 on 31 of 32 paths (Holm p = 4e-8).
6. *The SVSS anchors are slightly worse than random anchors* (Stage-1 main effect +0.00045, Holm p = 0.0032).

**Not supported**
1. *Lowest error.* At 317 sites the order is Random 0.0752, SDE 0.0753, Maximin 0.0754, AL-max 0.0755, hybrid 0.0758, SVTS 0.0760, AL 0.0855. The hybrid is the best of the seven headline methods on 0–1 of 32 paths at every budget.
2. *Earlier convergence.* At 95 sites the hybrid is worse than Random, SDE, Maximin, AL-max and SVTS (+0.0026 to +0.0047, all Holm-significant). It first falls below 0.085 at 158 sites; every other method except AL does so at 95. All methods are below 0.09 at 95 sites (AL: 0.0897).
3. *Better than SVTS.* Only AL-level differences are significant. Hybrid − SVTS is worse at 95 sites (+0.0031) and not different later.
4. *Better than, or equal to, SDE, Random and Maximin at small budgets.* The hybrid is worse than SDE at 95, 158 and 222 sites, and worse than Random and Maximin at 95 and 158 sites (all Holm-significant).
5. *Tighter spread.* The SD across the 24 FPGAs is 0.0041 for the hybrid and for Random and SDE at 317 sites; the SD across paths is 0.0241 for the hybrid and Random.
6. *A benefit of adaptive ASDE windows.* ASDE − fixed SDE −0.00002, ASDE − random fill +0.00023 (neither is significant); ASDE-only from a random or a maximin seed matches the hybrid.
7. *A benefit of the exploratory variants.* Anchors 1% and maximin anchors remove the early deficit of the hybrid (the hybrid is worse than both at 95 sites, Holm-significant) but are not better than Random (0.0756 against 0.0752 at 317).
8. *Any claim that an observed small difference is real beyond the 32 paths.* FPGAs are replicates of the same site sets, so agreement across FPGAs (the hybrid is worse than Random on 24 of 24 FPGAs at 317 sites) is not independent evidence; the path-level test is the relevant one.

## Things to know
- TEST FPGAs were evaluated once (768 path-jobs, one run, no failures). A dry run of this script on TUNE FPGA 2 had accidentally written its output files into the Stage 4 folder; they were moved out and deleted before the TEST run, and no TEST data was involved in the dry run.
- The mean number of training sites per FPGA equals the nominal budget to one decimal (95, 158, 222, 317; 1,587 for the 50% reference).
- GPR 50% reference: one seed only (about 20 s per fit); constant predictor: mean of a random 317-site training set, 5 seeds.
- Holm family: 46 unique comparisons (24 headline, 8 exploratory, 14 ablation; the hybrid − Random / SDE / AL-max comparisons at 317 sites appear once). Headline-only Holm p-values are in `stage4/tests.csv` and are smaller; the conclusions above use the stricter all-comparison correction. Several unadjusted CIs exclude 0 where the Holm p is above 0.05 (for example hybrid − Random at 317 sites).
- The hyperparameter bound shares on TEST are in the last table; the column lengthscale is at its lower bound in 21.5% of AL predictor fits and 3–6% of the others.
- Stage 5 (cost) has not been started.
