# Stage 1 – Fair baseline table (TUNE FPGAs)

## Protocol as run (assumptions in bold)
- Selection once per path on FPGA 1, same coordinates on all FPGAs, selection uses only FPGA 1 values of already selected sites. Predictor and the GPR inside selection: one function (`lib.fit`): additive RBF (column + row), raw coordinates, normalizer, initial lengthscales/variances/noise 1.0, bounds lengthscale 0.5–300 cells, kernel variance 1e-4–1e3, noise 1e-4–10, no restarts, at most 200 optimiser iterations, a fresh model for every fit.
- Held-out test set: a fixed random 30% of the active sites of each path (taken on FPGA 1, 946 sites on path 1, 951.8 test sites per path on average over FPGAs), never selectable and never fitted. The selectable pool is the other 70%. Secondary: the paper's old test set (all non-training sites).
- Budgets 95 / 158 / 222 / 317 sites including seed points, fixed for every method and path. 5 seeds (1–5): the random 1% seed set (32 sites) shared by AL, AL-max, SVTS, SVSS, maximin; the SDE shuffle; the ASDE scan order and fill; the random baseline.
- AL, AL-max, SVTS: 32-site random seed, then one site per refit (285 fits). SVSS: batches of 32 (32 / 32 / 31), 2 fits. Maximin starts from the same 32-site seed set (**assumption**: Euclidean distance on raw coordinates).
- Error = mean |mu − v| / d_spec: test sites → path → FPGA → mean over FPGAs. Mean ± SD over the 5 seeds of that FPGA-averaged error. Primary evaluation excludes FPGA 1; tables with FPGA 1 included are below.
- **Deviation from the Stage 1 text: TEST FPGAs (27–50) are not run here.** The protocol says TEST is touched once, at the end (Stage 4). Stage 1 and 2 use TUNE (FPGAs 2–26, 25 devices) only.
- **SVTS weights were chosen on TUNE only.** Mean TUNE error over the four budgets: (0.7 similarity, 0.3 variance) 0.0798; (0.3 similarity, 0.7 variance) 0.0786. The second is used as "SVTS" below; both are shown.
- Checks done before running: the copied ASDE gives the same sites as `hybrid_svss_asde.ipynb`'s `adaptive_sde` (same seed, same anchors); `sde_extend` gives the same sites as `sde.ipynb`'s `sde_select`; no held-out site is in any training set (asserted in every evaluation).

## Main table – TUNE FPGAs 2–26, held-out 30% test set (mean error ± SD over 5 seeds)
| method | 95 (3%) | 158 (5%) | 222 (7%) | 317 (10%) |
|---|---|---|---|---|
| AL (lowest variance) | 0.0888 ± 0.0008 | 0.0871 ± 0.0006 | 0.0858 ± 0.0005 | 0.0846 ± 0.0006 |
| AL-max (highest variance) | 0.0833 ± 0.0004 | 0.0784 ± 0.0004 | 0.0765 ± 0.0003 | 0.0747 ± 0.0003 |
| SVTS (0.3 sim, 0.7 var; chosen on TUNE) | 0.0828 ± 0.0011 | 0.0790 ± 0.0004 | 0.0774 ± 0.0005 | 0.0752 ± 0.0005 |
| SVTS (0.7 sim, 0.3 var) | 0.0840 ± 0.0009 | 0.0802 ± 0.0007 | 0.0783 ± 0.0006 | 0.0766 ± 0.0006 |
| SDE (fixed 5x5) | 0.0812 ± 0.0004 | 0.0774 ± 0.0003 | 0.0759 ± 0.0003 | 0.0746 ± 0.0002 |
| Random | 0.0815 ± 0.0007 | 0.0780 ± 0.0005 | 0.0761 ± 0.0005 | 0.0744 ± 0.0003 |
| Maximin (farthest point) | 0.0815 ± 0.0005 | 0.0779 ± 0.0002 | 0.0762 ± 0.0004 | 0.0746 ± 0.0003 |
| Hybrid SVSS+ASDE (current) | 0.0859 ± 0.0004 | 0.0790 ± 0.0007 | 0.0766 ± 0.0003 | 0.0751 ± 0.0003 |

## Same, FPGA 1 included (26 devices)
| method | 95 (3%) | 158 (5%) | 222 (7%) | 317 (10%) |
|---|---|---|---|---|
| AL (lowest variance) | 0.0893 ± 0.0008 | 0.0876 ± 0.0006 | 0.0863 ± 0.0005 | 0.0851 ± 0.0006 |
| AL-max (highest variance) | 0.0838 ± 0.0004 | 0.0788 ± 0.0004 | 0.0769 ± 0.0003 | 0.0751 ± 0.0003 |
| SVTS (0.3 sim, 0.7 var; chosen on TUNE) | 0.0832 ± 0.0011 | 0.0794 ± 0.0004 | 0.0778 ± 0.0005 | 0.0756 ± 0.0005 |
| SVTS (0.7 sim, 0.3 var) | 0.0845 ± 0.0009 | 0.0807 ± 0.0007 | 0.0788 ± 0.0006 | 0.0771 ± 0.0006 |
| SDE (fixed 5x5) | 0.0816 ± 0.0004 | 0.0779 ± 0.0003 | 0.0763 ± 0.0003 | 0.0750 ± 0.0002 |
| Random | 0.0820 ± 0.0007 | 0.0784 ± 0.0005 | 0.0765 ± 0.0005 | 0.0749 ± 0.0003 |
| Maximin (farthest point) | 0.0820 ± 0.0005 | 0.0784 ± 0.0002 | 0.0766 ± 0.0004 | 0.0750 ± 0.0003 |
| Hybrid SVSS+ASDE (current) | 0.0864 ± 0.0004 | 0.0795 ± 0.0007 | 0.0770 ± 0.0003 | 0.0755 ± 0.0003 |

## Secondary: paper's old test set (all non-training sites), TUNE 2–26
| method | 95 (3%) | 158 (5%) | 222 (7%) | 317 (10%) |
|---|---|---|---|---|
| AL (lowest variance) | 0.0896 ± 0.0008 | 0.0882 ± 0.0005 | 0.0872 ± 0.0006 | 0.0863 ± 0.0007 |
| AL-max (highest variance) | 0.0838 ± 0.0005 | 0.0790 ± 0.0004 | 0.0771 ± 0.0003 | 0.0755 ± 0.0002 |
| SVTS (0.3 sim, 0.7 var; chosen on TUNE) | 0.0834 ± 0.0011 | 0.0797 ± 0.0004 | 0.0783 ± 0.0004 | 0.0762 ± 0.0005 |
| SVTS (0.7 sim, 0.3 var) | 0.0846 ± 0.0009 | 0.0810 ± 0.0007 | 0.0792 ± 0.0005 | 0.0777 ± 0.0005 |
| SDE (fixed 5x5) | 0.0816 ± 0.0005 | 0.0779 ± 0.0003 | 0.0763 ± 0.0003 | 0.0749 ± 0.0003 |
| Random | 0.0820 ± 0.0007 | 0.0785 ± 0.0005 | 0.0765 ± 0.0005 | 0.0749 ± 0.0003 |
| Maximin (farthest point) | 0.0820 ± 0.0004 | 0.0784 ± 0.0002 | 0.0766 ± 0.0004 | 0.0751 ± 0.0003 |
| Hybrid SVSS+ASDE (current) | 0.0867 ± 0.0004 | 0.0797 ± 0.0005 | 0.0773 ± 0.0002 | 0.0757 ± 0.0002 |

## Selection cost (indicative only)
Measured while 8 selection tasks ran in parallel on a laptop that throttles; Stage 5 repeats it single-threaded. Hybrid: 168.1 of the 222 ASDE additions (75.7%) passed the exclusion rule, the rest is random fill.

| method | GPR fits per path | seconds per path (8 workers in parallel) |
|---|---|---|
| AL | 285.0 | 82.25 |
| AL-max | 285.0 | 83.70 |
| SVTS (chosen) | 285.0 | 83.62 |
| SVTS (0.7/0.3) | 285.0 | 81.96 |
| Hybrid | 2.0 | 0.32 |
| SDE | 0.0 | 0.07 |
| Maximin | 0.0 | 0.01 |
| Random | 0.0 | 0.00 |

## Hyperparameters at a bound (share of GPR fits; lower = at the lower bound, upper = at the upper bound; within 5% of the bound)
Selection-loop fits (all fits of the method's selection run; hybrid has only 2 fits per path, on 32 and 64 sites):

| method | column lengthscale lower / upper | row lengthscale lower / upper | column variance lower / upper | row variance lower / upper | noise lower / upper |
|---|---|---|---|---|---|
| al | 23.6% / 0.1% | 6.1% / 0.0% | 2.7% / 0.0% | 4.2% / 0.0% | 0.1% / 0.0% |
| almax | 5.5% / 0.2% | 7.7% / 0.0% | 2.3% / 0.0% | 2.2% / 0.0% | 2.2% / 0.0% |
| svts_b | 4.3% / 0.3% | 7.6% / 0.0% | 1.9% / 0.0% | 1.4% / 0.0% | 1.0% / 0.0% |
| svts_a | 4.6% / 0.3% | 6.1% / 0.0% | 1.8% / 0.0% | 2.0% / 0.0% | 1.0% / 0.0% |
| hybrid | 24.1% / 0.6% | 18.1% / 0.0% | 13.4% / 0.0% | 12.2% / 0.0% | 7.2% / 0.0% |

Predictor fits in the evaluation (all budgets, all FPGAs incl. FPGA 1):

| method | column lengthscale lower / upper | row lengthscale lower / upper | column variance lower / upper | row variance lower / upper | noise lower / upper |
|---|---|---|---|---|---|
| al | 21.6% / 0.0% | 5.5% / 0.1% | 1.9% / 0.0% | 4.7% / 0.0% | 0.0% / 0.0% |
| almax | 3.6% / 0.0% | 7.1% / 0.0% | 0.7% / 0.0% | 1.8% / 0.0% | 0.6% / 0.0% |
| svts_b | 3.0% / 0.1% | 7.5% / 0.0% | 0.4% / 0.0% | 1.5% / 0.0% | 0.1% / 0.0% |
| svts_a | 3.1% / 0.1% | 5.6% / 0.0% | 0.4% / 0.0% | 2.3% / 0.0% | 0.2% / 0.0% |
| sde | 4.7% / 0.0% | 7.6% / 0.0% | 0.5% / 0.0% | 3.3% / 0.0% | 0.0% / 0.0% |
| random | 5.4% / 0.0% | 6.0% / 0.0% | 0.4% / 0.0% | 3.5% / 0.0% | 0.0% / 0.0% |
| maximin | 4.3% / 0.0% | 7.3% / 0.0% | 0.2% / 0.0% | 3.4% / 0.0% | 0.0% / 0.0% |
| hybrid | 6.5% / 0.1% | 5.5% / 0.0% | 1.1% / 0.0% | 2.1% / 0.0% | 0.0% / 0.0% |

## Reading the table (plain words)
- AL (lowest variance) is the worst method at every budget (0.0846 at 317 sites; random 0.0744).
- At 95 sites (the SVSS anchors alone) the hybrid is worse than random, SDE, maximin, AL-max and SVTS (0.0859 against 0.0815 for random).
- At 317 sites the order is random 0.0744, SDE 0.0746, maximin 0.0746, AL-max 0.0747, hybrid 0.0751, SVTS 0.0752, AL 0.0846. The spread between random and the hybrid is 0.0007, about two seed SDs (0.0003); the spread between the best six methods is 0.0008.
- The hybrid is not the lowest-error method at any budget in this table. Random, SDE and maximin are at or below the hybrid at all four budgets. Significance tests are deferred to Stage 4.
- The seed SD is 0.0002–0.0011; it is largest for SVTS at 95 sites and for AL.
- The column lengthscale reaches its lower bound (0.5 cells) in 22% of AL predictor fits and in 3–7% of the predictor fits of the other methods; in the selection loops it is 24% for AL and for the hybrid's two small fits. Upper bounds are almost never reached.
- With FPGA 1 included, every mean moves by at most 0.0005; the order of the methods is unchanged at all four budgets.

## Differences from the old results
Different held-out test set (about 950 test sites per path; the old protocol tested on all non-selected sites, about 2,850), selection restricted to the 70% pool, fixed budgets for all methods, 1% random seed for AL / AL-max / SVTS (old: one site), bounded GPR, 5 seeds instead of one, TUNE devices instead of all 50.
