"""Equivalence test (TOST) on the Stage 4 TEST path-level differences: hybrid vs SVTS, AL-max, SDE, Random, Maximin at 158, 222, 317 sites. Margin +-0.001 (a choice)."""
from lib import *
from scipy import stats
MARGIN, ALPHA = 0.001, 0.05
R = pd.read_csv(f"{HERE}/stage4/eval_all.csv"); P = R.groupby(["method", "budget", "PATH"]).err.mean()
OTHERS = ["SVTS", "AL-max", "SDE", "Random", "Maximin"]; rows = []


def holm(p):
    p = np.asarray(p); o = np.argsort(p); m = len(p); adj = np.empty(m); run = 0
    for rank, i in enumerate(o): run = max(run, min(1, (m - rank) * p[i])); adj[i] = run
    return adj


for b in (158, 222, 317):
    for m in OTHERS:
        d = (P["Hybrid (frozen)", b] - P[m, b]).values; n = len(d); mean, se = d.mean(), d.std(ddof=1) / np.sqrt(n)
        p_lo, p_up = 1 - stats.t.cdf((mean + MARGIN) / se, n - 1), stats.t.cdf((mean - MARGIN) / se, n - 1); tcrit = stats.t.ppf(1 - ALPHA, n - 1)
        pw = max(stats.wilcoxon(d + MARGIN, alternative="greater").pvalue, stats.wilcoxon(d - MARGIN, alternative="less").pvalue)  # distribution-free sensitivity version
        rows.append(dict(budget=b, other=m, mean_diff=mean, ci90_low=mean - tcrit * se, ci90_high=mean + tcrit * se, p_tost_t=max(p_lo, p_up), p_tost_wilcoxon=pw, smallest_margin_for_equivalence=abs(mean) + tcrit * se))
D = pd.DataFrame(rows); D["p_holm_t"] = holm(D.p_tost_t); D["p_holm_wilcoxon"] = holm(D.p_tost_wilcoxon)
D["equivalent_t"] = D.p_holm_t < ALPHA; D["equivalent_wilcoxon"] = D.p_holm_wilcoxon < ALPHA; D.to_csv(f"{HERE}/stage4/equivalence.csv", index=False)
ref = R[(R.method == "Hybrid (frozen)") & (R.budget == 317)].groupby("seed").apply(lambda g: g.groupby("FPGA").err.mean().mean(), include_groups=False).mean()
t = ["| budget | hybrid vs | mean difference | 90% CI | TOST p (t) | Holm p | equivalent within ±0.001 (t) | Holm p (Wilcoxon version) | equivalent (Wilcoxon version) | smallest margin that would give equivalence |", "|---|---|---|---|---|---|---|---|---|---|"]
for _, r in D.iterrows(): t.append(f"| {r.budget} | {r.other} | {r.mean_diff:+.5f} | [{r.ci90_low:+.5f}, {r.ci90_high:+.5f}] | {r.p_tost_t:.3g} | {r.p_holm_t:.3g} | {'yes' if r.equivalent_t else 'no'} | {r.p_holm_wilcoxon:.3g} | {'yes' if r.equivalent_wilcoxon else 'no'} | {r.smallest_margin_for_equivalence:.5f} |")
md = f"""# Equivalence test (TOST) – TEST FPGAs 27–50

**The margin is a choice.** The equivalence margin is ±0.001 in error. I chose it as about 1.3% of the 10%-budget error of the hybrid ({ref:.4f}); it is not derived from the data or from a physical requirement. A different margin changes the answer; the last column shows the smallest margin at which each comparison would be called equivalent.

Method. For each comparison the 32 path-level differences (hybrid − other; mean over the 5 seeds and the 24 TEST FPGAs) are tested against H0: |mean difference| ≥ 0.001 with two one-sided tests at level 0.05 (paired t on the 32 paths; the TOST p-value is the larger of the two one-sided p-values). The 90% CI of the mean difference lying inside (−0.001, +0.001) is the same decision before correction. Holm correction over the 15 tests (3 budgets × 5 comparators). The "Wilcoxon version" shifts the differences by ±0.001 and uses one-sided Wilcoxon signed-rank tests (distribution-free sensitivity check, Holm over the same 15). Equivalence is declared only if the Holm-adjusted p < 0.05.

{chr(10).join(t)}

Reading: {int(D.equivalent_t.sum())} of 15 comparisons are equivalent within ±0.001 after Holm (t version); {int(D.equivalent_wilcoxon.sum())} of 15 with the Wilcoxon version. Equivalence at 95 sites was not tested (the hybrid is worse there).
"""
open(f"{HERE}/stage4_equivalence.md", "w", encoding="utf8").write(md); print(D.round(5).to_string(index=False))
