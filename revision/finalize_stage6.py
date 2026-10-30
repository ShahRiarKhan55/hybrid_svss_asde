import pandas as pd
HERE = r"C:\Users\sysdp\OneDrive - 島根大学\Desktop\Research\MY PAPERS\hybrid_svss_asde\revision"
r = open(f"{HERE}/stage6_report.md", encoding="utf8").read().split("\n## Sequential timing")[0].split("\n## Reading the tables")[0].rstrip()
T = pd.read_csv(f"{HERE}/stage6/sequential_timing.csv"); M = T.groupby(["predictor", "PATH"]).median().groupby("predictor").mean()
nm = {"add": "Additive RBF GPR", "iso": "One-lengthscale RBF GPR", "svr": "SVR (RBF)", "rf": "Random forest", "mlp": "MLP (tanh, 2 hidden layers)", "dgp": "Deep GP [1,3,2]"}
t = ["| predictor | fit incl. tuning, s | predict on ~2,850 sites, s |", "|---|---|---|"] + [f"| {nm[k]} | {M.loc[k, 'fit_s']:.2f} | {M.loc[k, 'predict_s']:.3f} |" for k in nm]
w = open(f"{HERE}/stage6_wide_report.md", encoding="utf8").read().split("\n", 2)[2].replace("\n## ", "\n### ").replace("## Protocol", "### Protocol").lstrip()
r += f"""

## Sequential timing (one thread, nothing else started by me, FPGA 27, hybrid seed-1 sites, 8 paths, median of 3 repeats; measured in the first Stage 6 session)
{chr(10).join(t)}

The deep GP function fits and predicts in one call; its time is reported under fit.

## Reading the tables (plain words, all 24 TEST FPGAs)
- On the hybrid sites the random forest has the lowest error (0.0746), then SVR (0.0757) and the additive GPR (0.0758); the MLP (0.0789), the one-lengthscale GPR (0.0797) and the deep GP (0.0798) are clearly worse. On Random sites the order is random forest 0.0743, additive GPR 0.0752, SVR 0.0758, MLP 0.0787, one-lengthscale GPR 0.0798, deep GP 0.0804.
- Random forest − additive GPR: −0.00124 on hybrid sites (CI −0.00171 to −0.00073, better on 27 of 32 paths, Holm p = 0.0012) and −0.00093 on Random sites (23 of 32, Holm p = 0.00074); about 1.2–1.6% of the error.
- SVR and the additive GPR are not separable on the hybrid sites (−0.00011, CI −0.00151 to +0.00093, Holm p = 0.44). On Random sites the additive GPR is better than SVR by 0.00059 (Holm p = 0.0032), although the bootstrap CI of the mean includes 0 (the Wilcoxon test is driven by the median shift).
- The one-lengthscale GPR, the MLP and the deep GP are worse than the additive GPR by 0.003–0.005 on both site types, on 31–32 of 32 paths (all Holm p ≤ 1.1e-8).
- Does the ranking depend on the site choice? The random forest is first on both. The only change is that SVR and the additive GPR swap places between second and third, which is within noise on the hybrid sites. Hybrid sites − Random sites: additive GPR +0.00059 (hybrid sites slightly worse, Holm p = 0.014) and deep GP −0.00065 (hybrid sites slightly better, Holm p = 0.008); no significant site effect for the other four predictors.
- GPR interval coverage is close to nominal for both GPRs and both site types: 90.8–91.1% for the 90% interval and 94.4–94.6% for the 95% interval.
- The additive GPR error on the hybrid sites (0.0758) equals the Stage 4 value, as it should.
- The wider-grid check (below, FPGAs 27 and 29): with wider SVR and MLP grids SVR and the additive GPR are tied on hybrid sites (0.0761 vs 0.0761; the order of the two flips in the fifth decimal) and the additive GPR stays ahead of SVR on Random sites (0.0756 vs 0.0762). The MLP improves by 0.0008 but remains 0.0027–0.0032 above the additive GPR (Holm p ≤ 2.3e-8).

## Things to know
- The grids of Stage 6 are my choice and coarse (MLP best width at a grid edge in 72% of fits, alpha at its smallest value in 61%; SVR epsilon at an edge in 64%). The wider-grid check addresses this only for FPGAs 27 and 29.
- The wider SVR grid is not a strict superset of the Stage 6 grid: C gains 1000; gamma gains 0.03, 0.3 and 3 but loses 0.01; epsilon becomes {{0.001, 0.01, 0.1}} and loses 0.3. The best epsilon is the upper edge (0.1) in 56% of wider-grid fits and at an edge in 78%, so a larger epsilon might still help SVR. The MLP grid gains width 128, alpha 1e-5/1e-4/1e-2 and relu, and loses alpha = 10 (chosen in 1% of Stage 6 fits).
- The MLP optimum is at the edge of the wider grid in about half of the fits (width 128 in 34%, alpha 0.1 in 28%), and relu is chosen in 85% of fits.
- The deep GP is the isotropic-kernel version from the earlier stages (not ARD).
- Times in the main-run table are with 8 jobs in parallel (indicative); the sequential table is the one to cite. Holm correction covers 16 comparisons in the main table and 10 in the wider-grid check.
- The main run was done in two sessions (8 FPGAs, then the remaining 16) with identical code, seeds and sites.

## Wider-grid check (FPGAs 27 and 29 only)
{w}
"""
open(f"{HERE}/stage6_report.md", "w", encoding="utf8").write(r); print(len(r.splitlines()), "lines")
