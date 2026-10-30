from lib import *
OUT = f"{HERE}/stage2"
G = pd.read_csv(f"{OUT}/results.csv"); Cn = pd.read_csv(f"{OUT}/contrasts.csv"); be = pd.read_csv(f"{OUT}/bounds_eval.csv", index_col=0)
C2 = pd.concat(pd.read_csv(f"{OUT}/cost_s{s}_p{p}.csv") for s in SEEDS for p in range(1, 33)).fillna(0); C1 = pd.read_csv(f"{HERE}/stage1/cost_all.csv")
hyb_acc = C1[C1.method == "hybrid"].asde_accepted.mean()
g = lambda fp, tn: G[(G.fpgas == fp) & (G.test == tn)].set_index("method")
f = lambda m, fp="TUNE 2-26", tn="held-out 30%": f"{g(fp, tn).loc[m, 'mean']:.4f} ± {g(fp, tn).loc[m, 'sd']:.4f}"
fac = ["| stage 1 (95 sites) \\ stage 2 (222 sites) | random fill | fixed SDE (2,2) | adaptive ASDE |", "|---|---|---|---|",
       f"| random | {f('rand+rand')} | {f('rand+SDE')} | {f('rand+ASDE')} |", f"| SVSS | {f('SVSS+rand')} | {f('SVSS+SDE')} | {f('SVSS+ASDE (hybrid)')} (hybrid) |"]
oth = ["| cell | TUNE 2–26, held-out (primary) | TUNE + FPGA 1, held-out | TUNE 2–26, old test set |", "|---|---|---|---|"]
for m in ["SVSS-only", "AL-max (1% seed)", "ASDE-only(rand seed)", "ASDE-only(maximin seed)", "SDE (Stage 1 reference)", "rand+rand", "SVSS+ASDE (hybrid)"]:
    oth.append(f"| {m} | {f(m)} | {f(m, 'TUNE + FPGA 1')} | {f(m, 'TUNE 2-26', 'old (all non-selected)')} |")
ct = ["| contrast (first − second; negative = first better) | mean difference | 95% CI | paths first better (of 32) | CI excludes 0 |", "|---|---|---|---|---|"]
for _, r in Cn.iterrows(): ct.append(f"| {r.contrast} | {r.mean_diff:+.5f} | [{r.ci95_low:+.5f}, {r.ci95_high:+.5f}] | {r.paths_first_better} | {'yes' if r.ci_excludes_0 else 'no'} |")
cc = C2.groupby("cell").agg(fits=("gpr_fits", "mean"), acc=("asde_accepted", "mean")); n_add = {"rand+ASDE": 222, "ASDE-only(rand seed)": 285, "ASDE-only(maximin seed)": 285}
cost = ["| cell | GPR fits per path | additions passing the ASDE rule |", "|---|---|---|", f"| SVSS+ASDE (hybrid, Stage 1) | 2.0 | {hyb_acc:.1f} of 222 ({100 * hyb_acc / 222:.1f}%) |"]
for c in ["rand+ASDE", "ASDE-only(rand seed)", "ASDE-only(maximin seed)", "SVSS-only", "SVSS+SDE", "SVSS+rand", "rand+SDE"]:
    cost.append(f"| {c} | {cc.loc[c, 'fits']:.1f} | " + (f"{cc.loc[c, 'acc']:.1f} of {n_add[c]} ({100 * cc.loc[c, 'acc'] / n_add[c]:.1f}%)" if c in n_add else "–") + " |")
P5 = ["ls0", "ls1", "var0", "var1", "noise"]; nm = {"ls0": "column lengthscale", "ls1": "row lengthscale", "var0": "column variance", "var1": "row variance", "noise": "noise"}
bt = ["| cell | " + " | ".join(f"{nm[q]} lower / upper" for q in P5) + " |", "|---|" + "---|" * 5]
for m in be.index: bt.append(f"| {m} | " + " | ".join(f"{100 * be.loc[m, q + '_lo']:.1f}% / {100 * be.loc[m, q + '_hi']:.1f}%" for q in P5) + " |")
sv = C2[C2.cell == "SVSS-only"]; prm = [f"{q}_{e}" for q in P5 for e in ("lo", "hi")]; sb = sv[prm].sum() / sv.gpr_fits.sum()
md = f"""# Stage 2 – Factorial ablation for R1.3 (TUNE FPGAs, 10% budget = 317 sites)

## Protocol
Same as Stage 1 (held-out 30% test set, additive bounded RBF GPR, 5 seeds, fixed budgets, selection on FPGA 1, TUNE FPGAs 2–26, FPGA 1 only as a secondary table). TEST FPGAs are not used. Split 95 sites (stage 1) + 222 sites (stage 2).
- Stage-1 choice: **random** = the first 95 sites of the random order; **SVSS** = the current SVSS anchors (alpha = beta = 0.5, 32 random seed sites, 2 fits).
- Stage-2 choice: **random fill**; **fixed SDE** = the SDE rule (window 2,2) continued from the anchors, with random fill if short; **adaptive ASDE** = the current ASDE (K = 2, tau = 2.0, per-anchor windows) from the anchors, with random fill if short.
- Extra rows: SVSS-only (SVSS continued to 317 sites in batches of 32, 9 fits); AL-max (highest predictive variance, same 32-site seed, one site per refit, from Stage 1); ASDE-only from the random 32-site seed; ASDE-only from a maximin (farthest-point) 32-site seed started at the first random seed site (**assumption**: Euclidean distance, raw coordinates).
- Paired differences: path-level values (mean over 5 seeds and the 25 TUNE FPGAs), 32 paths, 95% bootstrap CI over paths (10,000 resamples, seed 0). Not corrected for multiple comparisons; no p-values are used.

## Error at 317 sites, factorial cells (mean ± SD over 5 seeds; TUNE 2–26; held-out 30% test set)
{chr(10).join(fac)}

## Other rows
{chr(10).join(oth)}

## Paired differences (TUNE, primary test set)
{chr(10).join(ct)}

## Selection side facts
{chr(10).join(cost)}

Hyperparameters at a bound, predictor fits at 317 sites (share of fits, within 5% of the bound):

{chr(10).join(bt)}

Selection-loop fits of SVSS-only (9 per path): column lengthscale at the lower bound in {100 * sb['ls0_lo']:.1f}% of fits, row lengthscale {100 * sb['ls1_lo']:.1f}%, variances {100 * sb['var0_lo']:.1f}% / {100 * sb['var1_lo']:.1f}%, noise {100 * sb['noise_lo']:.1f}%; upper bounds below 0.3%.

## Reading the tables (plain words)
- At 317 sites the plain random baseline has the lowest error of the factorial cells (0.0744). The hybrid is 0.0751: worse than random by 0.00063 (CI +0.00032 to +0.00097); it is better than random on 9 of 32 paths. It is also worse than SDE (+0.00051, CI excludes 0).
- Stage-1 choice: SVSS anchors are worse than random anchors by 0.00046 (CI +0.00028 to +0.00064), averaged over the three stage-2 choices; SVSS is better on 4 of 32 paths.
- Stage-2 choice: adaptive ASDE, fixed SDE and random fill are not separable at this budget. ASDE − fixed SDE +0.00002 (CI −0.00020 to +0.00025); ASDE − random fill +0.00027 (CI −0.00003 to +0.00061); fixed SDE − random fill +0.00025 (CI −0.00009 to +0.00060). The SVSS × stage-2 interaction is also not separable from zero.
- Does the similarity term add anything? No. AL-max (plain highest variance) beats SVSS-only by 0.00301 (CI −0.00368 to −0.00236) and is better on 31 of 32 paths. The hybrid is better than SVSS-only by 0.00264 (CI -0.00334 to -0.00188); 222 of the hybrid's 317 sites come from ASDE or random fill, not from SVSS.
- ASDE-only is not handicapped by its random seed: ASDE-only from the random seed (0.0749) and from the maximin seed (0.0748) are within 0.0003 of the hybrid (0.0751), CIs include 0. The hybrid's SVSS anchors therefore give no measurable gain over a random or space-filling seed.
- SVSS-only (0.0777) is the worst cell; AL-max (0.0747) is among the best.
- With FPGA 1 included and on the old test set, the values shift by about 0.0005 and the order of the cells is almost the same (see the second and third columns of the second table).

## Things to know
- The ASDE window is set by the anchor spacing. With the maximin seed the anchors are far apart, the windows are large, and only {100 * cc.loc['ASDE-only(maximin seed)', 'acc'] / 285:.1f}% of the 285 additions pass the rule; the rest is random fill, so that cell is mostly a random fill after a space-filling seed. With the random seed {100 * cc.loc['ASDE-only(rand seed)', 'acc'] / 285:.1f}% pass; with the SVSS anchors {100 * hyb_acc / 222:.1f}%.
- Hyperparameters ending at a bound are rare at 317 sites (column lengthscale at the lower bound 0.8–2.8%); the shares are larger for smaller budgets (Stage 1).
- Stage 2 uses one budget (317 sites); the 95-site row of Stage 1 is the SVSS anchor set alone.
- Seeds: the SD over 5 seeds is 0.0002–0.0006; the CI over paths reflects path-to-path variation and uses the seed-averaged values.
"""
open(f"{HERE}/stage2_report.md", "w", encoding="utf8").write(md); print("ok")
