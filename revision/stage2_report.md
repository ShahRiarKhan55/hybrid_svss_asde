# Stage 2 – Factorial ablation for R1.3 (TUNE FPGAs, 10% budget = 317 sites)

## Protocol
Same as Stage 1 (held-out 30% test set, additive bounded RBF GPR, 5 seeds, fixed budgets, selection on FPGA 1, TUNE FPGAs 2–26, FPGA 1 only as a secondary table). TEST FPGAs are not used. Split 95 sites (stage 1) + 222 sites (stage 2).
- Stage-1 choice: **random** = the first 95 sites of the random order; **SVSS** = the current SVSS anchors (alpha = beta = 0.5, 32 random seed sites, 2 fits).
- Stage-2 choice: **random fill**; **fixed SDE** = the SDE rule (window 2,2) continued from the anchors, with random fill if short; **adaptive ASDE** = the current ASDE (K = 2, tau = 2.0, per-anchor windows) from the anchors, with random fill if short.
- Extra rows: SVSS-only (SVSS continued to 317 sites in batches of 32, 9 fits); AL-max (highest predictive variance, same 32-site seed, one site per refit, from Stage 1); ASDE-only from the random 32-site seed; ASDE-only from a maximin (farthest-point) 32-site seed started at the first random seed site (**assumption**: Euclidean distance, raw coordinates).
- Paired differences: path-level values (mean over 5 seeds and the 25 TUNE FPGAs), 32 paths, 95% bootstrap CI over paths (10,000 resamples, seed 0). Not corrected for multiple comparisons; no p-values are used.

## Error at 317 sites, factorial cells (mean ± SD over 5 seeds; TUNE 2–26; held-out 30% test set)
| stage 1 (95 sites) \ stage 2 (222 sites) | random fill | fixed SDE (2,2) | adaptive ASDE |
|---|---|---|---|
| random | 0.0744 ± 0.0003 | 0.0745 ± 0.0003 | 0.0747 ± 0.0003 |
| SVSS | 0.0748 ± 0.0003 | 0.0752 ± 0.0003 | 0.0751 ± 0.0003 (hybrid) |

## Other rows
| cell | TUNE 2–26, held-out (primary) | TUNE + FPGA 1, held-out | TUNE 2–26, old test set |
|---|---|---|---|
| SVSS-only | 0.0777 ± 0.0006 | 0.0783 ± 0.0006 | 0.0791 ± 0.0006 |
| AL-max (1% seed) | 0.0747 ± 0.0003 | 0.0751 ± 0.0003 | 0.0755 ± 0.0002 |
| ASDE-only(rand seed) | 0.0749 ± 0.0003 | 0.0753 ± 0.0003 | 0.0753 ± 0.0002 |
| ASDE-only(maximin seed) | 0.0748 ± 0.0002 | 0.0753 ± 0.0002 | 0.0753 ± 0.0003 |
| SDE (Stage 1 reference) | 0.0746 ± 0.0002 | 0.0750 ± 0.0002 | 0.0749 ± 0.0003 |
| rand+rand | 0.0744 ± 0.0003 | 0.0749 ± 0.0003 | 0.0749 ± 0.0003 |
| SVSS+ASDE (hybrid) | 0.0751 ± 0.0003 | 0.0755 ± 0.0003 | 0.0757 ± 0.0002 |

## Paired differences (TUNE, primary test set)
| contrast (first − second; negative = first better) | mean difference | 95% CI | paths first better (of 32) | CI excludes 0 |
|---|---|---|---|---|
| hybrid - rand+rand | +0.00063 | [+0.00032, +0.00097] | 9 | yes |
| hybrid - rand+SDE | +0.00054 | [+0.00024, +0.00083] | 5 | yes |
| hybrid - rand+ASDE | +0.00036 | [+0.00004, +0.00067] | 9 | yes |
| hybrid - SVSS+rand | +0.00027 | [-0.00006, +0.00061] | 12 | no |
| hybrid - SVSS+SDE | -0.00013 | [-0.00047, +0.00025] | 20 | no |
| hybrid - SVSS-only | -0.00264 | [-0.00334, -0.00188] | 28 | yes |
| hybrid - AL-max (1% seed) | +0.00037 | [-0.00014, +0.00088] | 14 | no |
| hybrid - ASDE-only(rand seed) | +0.00022 | [-0.00013, +0.00054] | 10 | no |
| hybrid - ASDE-only(maximin seed) | +0.00023 | [-0.00012, +0.00058] | 11 | no |
| hybrid - SDE (Stage 1 reference) | +0.00051 | [+0.00022, +0.00081] | 8 | yes |
| AL-max - SVSS-only | -0.00301 | [-0.00368, -0.00236] | 31 | yes |
| ASDE-only(maximin seed) - ASDE-only(rand seed) | -0.00002 | [-0.00020, +0.00017] | 16 | no |
| main effect stage 1: SVSS - random (mean of 3 stage-2 choices) | +0.00046 | [+0.00028, +0.00064] | 4 | yes |
| main effect stage 2: ASDE - SDE-fixed (mean of 2 stage-1 choices) | +0.00002 | [-0.00020, +0.00025] | 17 | no |
| main effect stage 2: ASDE - random fill | +0.00027 | [-0.00003, +0.00061] | 13 | no |
| main effect stage 2: SDE-fixed - random fill | +0.00025 | [-0.00009, +0.00060] | 16 | no |
| interaction: (SVSS - random) under ASDE minus under SDE-fixed | -0.00030 | [-0.00078, +0.00018] | 18 | no |

## Selection side facts
| cell | GPR fits per path | additions passing the ASDE rule |
|---|---|---|
| SVSS+ASDE (hybrid, Stage 1) | 2.0 | 168.1 of 222 (75.7%) |
| rand+ASDE | 0.0 | 93.5 of 222 (42.1%) |
| ASDE-only(rand seed) | 0.0 | 44.8 of 285 (15.7%) |
| ASDE-only(maximin seed) | 0.0 | 9.9 of 285 (3.5%) |
| SVSS-only | 9.0 | – |
| SVSS+SDE | 0.0 | – |
| SVSS+rand | 0.0 | – |
| rand+SDE | 0.0 | – |

Hyperparameters at a bound, predictor fits at 317 sites (share of fits, within 5% of the bound):

| cell | column lengthscale lower / upper | row lengthscale lower / upper | column variance lower / upper | row variance lower / upper | noise lower / upper |
|---|---|---|---|---|---|
| AL-max (1% seed) | 1.5% / 0.0% | 0.9% / 0.0% | 0.0% / 0.0% | 0.3% / 0.0% | 0.0% / 0.0% |
| ASDE-only(maximin seed) | 0.8% / 0.0% | 2.7% / 0.0% | 0.0% / 0.0% | 1.1% / 0.0% | 0.0% / 0.0% |
| ASDE-only(rand seed) | 1.0% / 0.0% | 2.3% / 0.0% | 0.0% / 0.0% | 1.0% / 0.0% | 0.0% / 0.0% |
| SDE (Stage 1 reference) | 1.5% / 0.0% | 3.1% / 0.0% | 0.0% / 0.0% | 0.8% / 0.0% | 0.0% / 0.0% |
| SVSS+ASDE (hybrid) | 1.2% / 0.0% | 1.3% / 0.0% | 0.0% / 0.0% | 0.4% / 0.0% | 0.0% / 0.0% |
| SVSS+SDE | 1.2% / 0.0% | 2.6% / 0.0% | 0.1% / 0.0% | 0.8% / 0.0% | 0.0% / 0.0% |
| SVSS+rand | 1.3% / 0.0% | 2.0% / 0.0% | 0.0% / 0.0% | 0.9% / 0.0% | 0.0% / 0.0% |
| SVSS-only | 2.8% / 0.0% | 2.9% / 0.0% | 0.0% / 0.0% | 2.0% / 0.0% | 0.0% / 0.0% |
| rand+ASDE | 1.4% / 0.0% | 1.7% / 0.0% | 0.0% / 0.0% | 1.1% / 0.0% | 0.0% / 0.0% |
| rand+SDE | 0.9% / 0.0% | 3.6% / 0.0% | 0.0% / 0.0% | 0.8% / 0.0% | 0.0% / 0.0% |
| rand+rand | 1.1% / 0.0% | 2.4% / 0.0% | 0.0% / 0.0% | 1.1% / 0.0% | 0.0% / 0.0% |

Selection-loop fits of SVSS-only (9 per path): column lengthscale at the lower bound in 10.1% of fits, row lengthscale 6.9%, variances 3.8% / 3.8%, noise 1.6%; upper bounds below 0.3%.

## Reading the tables (plain words)
- At 317 sites the plain random baseline has the lowest error of the factorial cells (0.0744). The hybrid is 0.0751: worse than random by 0.00063 (CI +0.00032 to +0.00097); it is better than random on 9 of 32 paths. It is also worse than SDE (+0.00051, CI excludes 0).
- Stage-1 choice: SVSS anchors are worse than random anchors by 0.00046 (CI +0.00028 to +0.00064), averaged over the three stage-2 choices; SVSS is better on 4 of 32 paths.
- Stage-2 choice: adaptive ASDE, fixed SDE and random fill are not separable at this budget. ASDE − fixed SDE +0.00002 (CI −0.00020 to +0.00025); ASDE − random fill +0.00027 (CI −0.00003 to +0.00061); fixed SDE − random fill +0.00025 (CI −0.00009 to +0.00060). The SVSS × stage-2 interaction is also not separable from zero.
- Does the similarity term add anything? No. AL-max (plain highest variance) beats SVSS-only by 0.00301 (CI −0.00368 to −0.00236) and is better on 31 of 32 paths. The hybrid is better than SVSS-only by 0.00264 (CI -0.00334 to -0.00188); 222 of the hybrid's 317 sites come from ASDE or random fill, not from SVSS.
- ASDE-only is not handicapped by its random seed: ASDE-only from the random seed (0.0749) and from the maximin seed (0.0748) are within 0.0003 of the hybrid (0.0751), CIs include 0. The hybrid's SVSS anchors therefore give no measurable gain over a random or space-filling seed.
- SVSS-only (0.0777) is the worst cell; AL-max (0.0747) is among the best.
- With FPGA 1 included and on the old test set, the values shift by about 0.0005 and the order of the cells is almost the same (see the second and third columns of the second table).

## Things to know
- The ASDE window is set by the anchor spacing. With the maximin seed the anchors are far apart, the windows are large, and only 3.5% of the 285 additions pass the rule; the rest is random fill, so that cell is mostly a random fill after a space-filling seed. With the random seed 15.7% pass; with the SVSS anchors 75.7%.
- Hyperparameters ending at a bound are rare at 317 sites (column lengthscale at the lower bound 0.8–2.8%); the shares are larger for smaller budgets (Stage 1).
- Stage 2 uses one budget (317 sites); the 95-site row of Stage 1 is the SVSS anchor set alone.
- Seeds: the SD over 5 seeds is 0.0002–0.0006; the CI over paths reflects path-to-path variation and uses the seed-averaged values.
