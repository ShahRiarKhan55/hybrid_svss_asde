
## Claims the data support and do not support (TEST FPGAs 27–50, held-out test set, Holm over 46 comparisons)

**Supported**
1. *The hybrid has lower error than AL (lowest variance) at every budget.* Differences −0.0030 (95 sites) to −0.0097 (317); the hybrid wins on 28–31 of 32 paths; Holm p ≤ 1.3e-6 at all four budgets.
2. *From 158 sites on, the hybrid is not distinguishable from SVTS (weights chosen on TUNE) or from AL-max.* Hybrid − SVTS: −0.00005 (158), −0.00080 (222), −0.00020 (317); hybrid − AL-max: +0.00069, +0.00002, +0.00035. All Holm p = 1. This is "no difference detected", not proof of equivalence.
3. *The hybrid needs far fewer GPR fits in selection.* 2 fits per path against 285 for AL, AL-max and SVTS (fit counts from Stage 1; the timing is Stage 5). SDE, Random and Maximin need no GPR fit.
4. *At 222 and 317 sites the hybrid is within 0.0004–0.0006 of Random, SDE and Maximin* (2–2.4% of the 0.0245 headroom between the constant predictor and a GPR on 50% of the sites), and none of these differences is significant after Holm at 222 or 317 sites. (At 317 the unadjusted CIs of hybrid − random and hybrid − SDE exclude 0; the Holm p-values are 0.067 and 0.077.)
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
