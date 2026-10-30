from lib import *
g = pd.read_csv(f"{HERE}/stage4/results_primary.csv").set_index(["method", "budget"]); D = pd.read_csv(f"{HERE}/stage4/tests.csv"); E = pd.read_csv(f"{HERE}/stage4/equivalence.csv")
S = pd.read_csv(f"{HERE}/stage5/selection_timing.csv"); F = pd.read_csv(f"{HERE}/stage5/final_timing.csv")
HEAD = ["Hybrid (frozen)", "AL", "AL-max", "SVTS", "SDE", "Random", "Maximin"]; EX = ["anchors 1% (exploratory)", "maximin anchors (exploratory)"]
cell = lambda m, b: f"{g.loc[(m, b), 'mean']:.4f} ± {g.loc[(m, b), 'sd']:.4f}"
th = lambda m, t: next((str(b) for b in BUDGETS if g.loc[(m, b), "mean"] < t), "--")
t1 = ["| method | 95 sites (3%) | 158 (5%) | 222 (7%) | 317 (10%) | first budget < 0.09 | first budget < 0.085 |", "|---|---|---|---|---|---|---|"]
for m in HEAD + EX: t1.append(f"| {m} | " + " | ".join(cell(m, b) for b in BUDGETS) + f" | {th(m, .09)} | {th(m, .085)} |")
t1.append(f"| constant predictor (reference) | – | – | – | {cell('constant', 317)} | – | – |"); t1.append(f"| GPR on 50% of the sites, 1587 sites (reference, seed 1 only) | – | – | – | {g.loc[('gpr50', 1587), 'mean']:.4f} | – | – |")
A = lambda m: cell(m, 317)
fac = ["| stage 1 (95 sites) \\ stage 2 (222 sites) | random fill | fixed SDE (2,2) | adaptive ASDE |", "|---|---|---|---|", f"| random | {A('Random')} | {A('rand+SDE')} | {A('rand+ASDE')} |", f"| SVSS | {A('SVSS+rand')} | {A('SVSS+SDE')} | {A('Hybrid (frozen)')} (hybrid) |"]
oth = ["| other cells @ 317 sites | error |", "|---|---|"] + [f"| {m} | {A(m)} |" for m in ["SVSS-only", "AL-max", "ASDE-only(rand seed)", "ASDE-only(maximin seed)", "SDE"]]
ct = ["| contrast @ 317 sites (negative = first lower error) | mean difference | 95% CI | first / second wins (of 32 paths) | Holm p (46 comparisons) |", "|---|---|---|---|---|"]
for _, r in D[D.family == "A"].iterrows(): ct.append(f"| {r.comparison.replace(' @ 317', '')} | {r.mean_diff:+.5f} | [{r.ci_low:+.5f}, {r.ci_high:+.5f}] | {r.first_wins} / {r.second_wins} | {r.p_holm_all:.2g} |")
M = S.groupby(["method", "PATH"]).agg(wall=("wall_s", "median"), fits=("gpr_fits", "first")).reset_index(); nm = {"al": "AL", "almax": "AL-max", "svts": "SVTS", "sde": "SDE", "random": "Random", "maximin": "Maximin", "hybrid": "Hybrid"}
hb = M[M.method == "hybrid"].wall.mean(); c3 = ["| method | GPR fits per path | selection time per path, s (mean over 32 paths of the median of 3 repeats) | min – max over paths | relative to hybrid |", "|---|---|---|---|---|"]
for m in ["al", "almax", "svts", "hybrid", "sde", "maximin", "random"]:
    x = M[M.method == m]; c3.append(f"| {nm[m]} | {x.fits.mean():.0f} | {x.wall.mean():.3f}" + (" (below 1 ms)" if x.wall.mean() < 0.001 else "") + f" | {x.wall.min():.3f} – {x.wall.max():.3f} | {x.wall.mean() / hb:.2f}x |")
FM = F.groupby(["FPGA", "PATH"]).agg(fit=("fit_s", "median"), pred=("predict_s", "median"), npr=("n_predict", "first"))
c3 += ["", f"Final GPR fit on 317 sites (same for every method): {FM.fit.mean():.3f} s (min–max {FM.fit.min():.3f}–{FM.fit.max():.3f}); prediction on {FM.npr.mean():.0f} sites: {FM.pred.mean():.3f} s. One thread, sequential, 3 repeats, median. Selection times varied by about 30% between repeats (mean (max − min)/median over path and method pairs), the machine ran on the Balanced power plan at 1.8–2.3 GHz with other services active, so treat the absolute seconds as indicative; the fit counts and the ratio (about {M[M.method == 'al'].wall.mean() / hb:.0f}x for AL, {M[M.method == 'svts'].wall.mean() / hb:.0f}x for SVTS) are the robust numbers. Details: `stage5_report.md`."]
eq = ["| budget | hybrid vs | mean difference | 90% CI | Holm p (TOST) | equivalent within ±0.001 | smallest margin that would give equivalence |", "|---|---|---|---|---|---|---|"]
for _, r in E.iterrows(): eq.append(f"| {r.budget} | {r.other} | {r.mean_diff:+.5f} | [{r.ci90_low:+.5f}, {r.ci90_high:+.5f}] | {r.p_holm_t:.3g} | {'yes' if r.equivalent_t else 'no'} | {r.smallest_margin_for_equivalence:.5f} |")
md = f"""# Summary of the revision experiments

All numbers: TEST FPGAs 27–50 (24 devices), primary protocol (held-out 30% test set of about 952 sites per path, additive bounded RBF GPR inside selection and for prediction, selection on FPGA 1 only, error = |mu − v| / d_spec averaged test sites → path → FPGA, mean ± SD over 5 seeds). Design choices were made on TUNE FPGAs 2–26 only; TEST was evaluated once. Details and secondary tables (FPGA 1 included, the paper's old test set): `stage4_report.md`. Stage reports: `stage0_audit.md`, `stage1_report.md`, `stage2_report.md`, `stage3_report.md`, `stage4_report.md`, `stage4_equivalence.md`, `stage5_report.md`.

## 1. Final numbers (TEST, primary protocol)
{chr(10).join(t1)}

Hybrid = SVSS (95 anchors, alpha = beta = 0.5) + ASDE, frozen unchanged after Stage 3 (no candidate met the pre-written rule). SVTS weights (0.3 similarity, 0.7 variance) were chosen on TUNE. The two exploratory rows were defined before any result and are not part of the headline comparison. Budgets include seed points.

## 2. Ablation (TEST, 317 sites)
{chr(10).join(fac)}

{chr(10).join(oth)}

Paired differences (path-level, 32 paths, 95% bootstrap CI; Holm over all 46 comparisons of Stage 4):

{chr(10).join(ct)}

## 3. Cost (final designs)
{chr(10).join(c3)}

## 4. Equivalence (TOST, hybrid vs the others, TEST path-level differences)
**The margin of ±0.001 in error (about 1.3% of the hybrid's 10%-budget error) is a choice, not derived from data.** Paired-t TOST on the 32 path-level differences, Holm over the 15 tests; a Wilcoxon-based version gives the same 4 equivalences (`stage4_equivalence.md`).

{chr(10).join(eq)}

Equivalent within ±0.001 after Holm: hybrid vs AL-max at 222 sites, vs Random at 222, vs SDE at 317, vs Maximin at 317 (4 of 15). Not shown for SVTS at any budget (at 317 the 90% CI is [−0.00073, +0.00034], Holm p = 0.086), for the others at 158, or for Random and AL-max at 317 (Holm p = 0.10).

## 5. Claims
**Supported by the data**
1. The hybrid has lower error than AL at every budget (−0.0030 to −0.0097; wins on 28–31 of 32 paths; Holm p ≤ 1.3e-6).
2. Selection cost: the hybrid needs 2 GPR fits and about {M[M.method == 'hybrid'].wall.mean():.2f} s per path, against 285 fits and about {M[M.method == 'al'].wall.mean():.0f}–{M[M.method == 'svts'].wall.mean():.0f} s for AL, AL-max and SVTS (roughly {M[M.method == 'al'].wall.mean() / hb:.0f}–{M[M.method == 'svts'].wall.mean() / hb:.0f} times faster on this machine). SDE, Maximin and Random need no GPR fit and are cheaper still.
3. From 158 sites on, no difference between the hybrid and SVTS or AL-max was detected (|difference| ≤ 0.0008, Holm p = 1). The error of AL-max is statistically equivalent to the hybrid's within ±0.001 at 222 sites only.
4. At 222 sites the hybrid is equivalent to Random within ±0.001 (difference +0.00048), and at 317 sites it is equivalent to SDE and to Maximin (+0.00050 and +0.00042). At the other budgets or comparators equivalence was not established.
5. SVSS's similarity term adds nothing: SVSS-only is worse than the hybrid (+0.0027, Holm p = 4e-5) and plain highest-variance AL-max beats SVSS-only by 0.0031 on 31 of 32 paths.
6. SVSS anchors are slightly worse than random anchors (+0.00045, Holm p = 0.0032); the choice of stage-2 rule (adaptive ASDE, fixed SDE, random fill) makes no detectable difference at 317 sites; ASDE-only matches the hybrid.

**Not supported by the data**
1. Lowest error. At 317 sites: Random 0.0752, SDE 0.0753, Maximin 0.0754, AL-max 0.0755, hybrid 0.0758, SVTS 0.0760, AL 0.0855. The hybrid is the best of the seven on 0–1 of 32 paths at every budget.
2. Earlier convergence. At 95 sites the hybrid is worse than Random, SDE, Maximin, AL-max and SVTS (+0.0026 to +0.0047, Holm-significant); it first falls below 0.085 at 158 sites, the others (except AL) at 95.
3. Better than SVTS. Worse at 95 sites (+0.0031) and no difference later; equivalence with SVTS within ±0.001 is not shown either.
4. Equal to or better than SDE, Random or Maximin at small budgets: the hybrid is worse than SDE at 95, 158 and 222 sites and worse than Random and Maximin at 95 and 158 sites (Holm-significant). At 317 sites it is nominally 0.0004–0.0006 above them; the two-sided tests are not significant after Holm (p = 0.067–0.19) while equivalence is shown for SDE and Maximin only.
5. Tighter spread. SD across the 24 FPGAs at 317 sites 0.0041 for the hybrid, Random and SDE; SD across paths 0.0241 for the hybrid and Random.
6. A benefit of the adaptive ASDE windows, of the SVSS anchors, or of the cosine-without-centring similarity (Stage 3: no candidate met the rule).
7. That the exploratory variants (anchors 1%, maximin anchors) beat Random; they remove the hybrid's early deficit but are at the Random level.
8. Any inference beyond these 32 paths and this device family: FPGAs reuse the same site sets, so agreement across FPGAs (the hybrid is worse than Random on 24 of 24 FPGAs at 317 sites) is not independent evidence.

**Limits of this work.** One device family, one predictor, 5 seeds. The equivalence margin is a choice. Timings come from a non-idle laptop in the Balanced power plan with about 30% repeat-to-repeat variation. The 50%-GPR reference uses one seed. Stage 3 searched only the trimmed list of candidates agreed beforehand.
"""
open(f"{HERE}/summary.md", "w", encoding="utf8").write(md); print(len(md.splitlines()), "lines")
