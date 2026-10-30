# Equivalence test (TOST) – TEST FPGAs 27–50

**The margin is a choice.** The equivalence margin is ±0.001 in error. I chose it as about 1.3% of the 10%-budget error of the hybrid (0.0758); it is not derived from the data or from a physical requirement. A different margin changes the answer; the last column shows the smallest margin at which each comparison would be called equivalent.

Method. For each comparison the 32 path-level differences (hybrid − other; mean over the 5 seeds and the 24 TEST FPGAs) are tested against H0: |mean difference| ≥ 0.001 with two one-sided tests at level 0.05 (paired t on the 32 paths; the TOST p-value is the larger of the two one-sided p-values). The 90% CI of the mean difference lying inside (−0.001, +0.001) is the same decision before correction. Holm correction over the 15 tests (3 budgets × 5 comparators). The "Wilcoxon version" shifts the differences by ±0.001 and uses one-sided Wilcoxon signed-rank tests (distribution-free sensitivity check, Holm over the same 15). Equivalence is declared only if the Holm-adjusted p < 0.05.

| budget | hybrid vs | mean difference | 90% CI | TOST p (t) | Holm p | equivalent within ±0.001 (t) | Holm p (Wilcoxon version) | equivalent (Wilcoxon version) | smallest margin that would give equivalence |
|---|---|---|---|---|---|---|---|---|---|
| 158 | SVTS | -0.00005 | [-0.00073, +0.00064] | 0.0126 | 0.102 | no | 0.125 | no | 0.00073 |
| 158 | AL-max | +0.00069 | [-0.00005, +0.00142] | 0.236 | 1 | no | 1 | no | 0.00142 |
| 158 | SDE | +0.00162 | [+0.00102, +0.00222] | 0.956 | 1 | no | 1 | no | 0.00222 |
| 158 | Random | +0.00101 | [+0.00053, +0.00149] | 0.513 | 1 | no | 1 | no | 0.00149 |
| 158 | Maximin | +0.00107 | [+0.00063, +0.00151] | 0.607 | 1 | no | 1 | no | 0.00151 |
| 222 | SVTS | -0.00080 | [-0.00146, -0.00013] | 0.303 | 1 | no | 1 | no | 0.00146 |
| 222 | AL-max | +0.00002 | [-0.00054, +0.00057] | 0.00262 | 0.0315 | yes | 0.0253 | yes | 0.00057 |
| 222 | SDE | +0.00072 | [+0.00035, +0.00110] | 0.109 | 0.656 | no | 0.672 | no | 0.00110 |
| 222 | Random | +0.00048 | [+0.00022, +0.00074] | 0.00102 | 0.0153 | yes | 0.0155 | yes | 0.00074 |
| 222 | Maximin | +0.00042 | [+0.00001, +0.00083] | 0.0114 | 0.102 | no | 0.0666 | no | 0.00083 |
| 317 | SVTS | -0.00020 | [-0.00073, +0.00034] | 0.0078 | 0.0858 | no | 0.096 | no | 0.00073 |
| 317 | AL-max | +0.00035 | [-0.00009, +0.00080] | 0.0102 | 0.102 | no | 0.0578 | no | 0.00080 |
| 317 | SDE | +0.00050 | [+0.00024, +0.00076] | 0.00151 | 0.0211 | yes | 0.0253 | yes | 0.00076 |
| 317 | Random | +0.00059 | [+0.00029, +0.00088] | 0.0124 | 0.102 | no | 0.0578 | no | 0.00088 |
| 317 | Maximin | +0.00042 | [+0.00011, +0.00073] | 0.00159 | 0.0211 | yes | 0.0253 | yes | 0.00073 |

Reading: 4 of 15 comparisons are equivalent within ±0.001 after Holm (t version); 4 of 15 with the Wilcoxon version. Equivalence at 95 sites was not tested (the hybrid is worse there).
