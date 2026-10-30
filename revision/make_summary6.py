import pandas as pd, numpy as np, re
H = r"C:\Users\sysdp\OneDrive - 島根大学\Desktop\Research\MY PAPERS\hybrid_svss_asde\revision"
R = pd.read_csv(f"{H}/stage6/res_all.csv"); D = pd.read_csv(f"{H}/stage6/tests.csv"); Q = pd.read_csv(f"{H}/stage6/sequential_timing.csv")
W = pd.read_csv(f"{H}/stage6_wide/res_all.csv"); WD = pd.read_csv(f"{H}/stage6_wide/tests.csv")
NM = {"add": "Additive RBF GPR", "iso": "One-lengthscale RBF GPR", "svr": "SVR (RBF)", "rf": "Random forest", "mlp": "MLP (tanh, 2 hidden layers)", "dgp": "Deep GP [1,3,2]"}
seedm = lambda d, by: d.groupby(by + ["seed", "FPGA"]).err.mean().groupby(by + ["seed"]).mean().groupby(by).agg(["mean", "std"])
G = seedm(R, ["site", "predictor"]); GW = seedm(W, ["site", "predictor"])
ce = lambda s, k: f"{G.loc[(s, k), 'mean']:.4f} ± {G.loc[(s, k), 'std']:.4f}"
T = Q.groupby(["predictor", "PATH"]).median().groupby("predictor").mean()
cov = R[R.predictor.isin(["add", "iso"])].groupby(["site", "predictor"])[["cov90", "cov95"]].mean()
nfp = R.FPGA.nunique()


def diff(site, k):
    r = D[D.comparison == f"{NM[k]} − additive GPR, {site} sites"].iloc[0]; return r
def site_eff(k): return D[D.comparison == f"{NM[k]}: hybrid sites − Random sites"].iloc[0]
fmt = lambda r: f"{r.mean_diff:+.5f} [{r.ci_low:+.5f}, {r.ci_high:+.5f}], {r.first_wins}/{r.second_wins} paths, Holm p = {r.p_holm:.2g}"

cv = ["| predictor | tuned hyperparameters (grid) | cross-validation | fixed settings |", "|---|---|---|---|",
      "| SVR (RBF kernel) | C ∈ {0.1, 1, 10, 100}; gamma ∈ {0.01, 0.1, 1, 10}; epsilon ∈ {0.01, 0.1, 0.3} (48 settings) | 5-fold `KFold(shuffle=True, random_state=0)` on the 317 training sites, scoring = negative mean absolute error of the predictions in the original value units | inputs standardised (`StandardScaler`, fitted inside each fold), target standardised (`TransformedTargetRegressor`); refit on all 317 sites with the best setting |",
      "| Random forest | min_samples_leaf ∈ {1, 3, 5, 10}; max_features ∈ {1.0, 0.5} (8 settings) | same | 200 trees, `random_state=0`; same input/target standardisation |",
      "| MLP | hidden layers (16,16), (32,32), (64,64); alpha ∈ {0.001, 0.1, 10} (9 settings) | same | tanh, Adam (scikit-learn default solver), `max_iter=1000`, no early stopping, `random_state=0`; same standardisation |",
      "| SVR, wider grid (FPGAs 27 and 29) | C ∈ {0.1, 1, 10, 100, 1000}; gamma ∈ {0.03, 0.1, 0.3, 1, 3, 10}; epsilon ∈ {0.001, 0.01, 0.1} (90 settings) | same | as above |",
      "| MLP, wider grid (FPGAs 27 and 29) | widths (16,16), (32,32), (64,64), (128,128); alpha ∈ {1e-5, 1e-4, 1e-3, 1e-2, 1e-1}; activation ∈ {tanh, relu} (40 settings) | same | Adam, `max_iter=1000`, `random_state=0` |",
      "| Additive and one-lengthscale GPR | none (hyperparameters by marginal likelihood) | – | same bounded GPR in all stages: lengthscale 0.5–300, kernel variance 1e-4–1e3, noise 1e-4–10, at most 200 iterations, no restarts |",
      "| Deep GP | none (variational bound) | – | layers [1,3,2], 50 inducing points, 200 iterations, isotropic RBF kernels |"]
rt = ["| predictor | hybrid sites | Random sites |", "|---|---|---|"] + [f"| {NM[k]} | {ce('hybrid', k)} | {ce('random', k)} |" for k in NM]
dt = ["| predictor − additive GPR | hybrid sites: mean difference [95% CI], paths won (predictor / GPR), Holm p | Random sites |", "|---|---|---|"] + [f"| {NM[k]} | {fmt(diff('hybrid', k))} | {fmt(diff('random', k))} |" for k in NM if k != "add"]
st = ["| predictor | hybrid sites − Random sites: mean difference [95% CI], paths where hybrid sites win/lose, Holm p |", "|---|---|"] + [f"| {NM[k]} | {fmt(site_eff(k))} |" for k in NM]
wk = {"add": "Additive RBF GPR", "svr_narrow": "SVR, Stage 6 grid", "svr_wide": "SVR, wider grid", "mlp_narrow": "MLP, Stage 6 grid (tanh)", "mlp_wide": "MLP, wider grid (tanh or relu)"}
wt = ["| predictor (FPGAs 27 and 29) | hybrid sites | Random sites |", "|---|---|---|"] + [f"| {n} | {GW.loc[('hybrid', k), 'mean']:.4f} ± {GW.loc[('hybrid', k), 'std']:.4f} | {GW.loc[('random', k), 'mean']:.4f} ± {GW.loc[('random', k), 'std']:.4f} |" for k, n in wk.items()]
wdt = ["| wider-grid comparison | site | mean difference [95% CI] | paths won (first/second) | Holm p (10 tests) |", "|---|---|---|---|---|"] + [f"| {r.comparison} | {r.site} | {r.mean_diff:+.5f} [{r.ci_low:+.5f}, {r.ci_high:+.5f}] | {r.first_wins}/{r.second_wins} | {r.p_holm:.2g} |" for _, r in WD.iterrows()]
q = ["| predictor | error, hybrid sites | error, Random sites | difference to additive GPR, hybrid sites [95% CI] | difference to additive GPR, Random sites [95% CI] | fit incl. tuning, s | predict, s |", "|---|---|---|---|---|---|---|"]
for k in NM:
    dh = "–" if k == "add" else (lambda r: f"{r.mean_diff:+.5f} [{r.ci_low:+.5f}, {r.ci_high:+.5f}] (Holm p = {r.p_holm:.2g})")(diff("hybrid", k)); dr = "–" if k == "add" else (lambda r: f"{r.mean_diff:+.5f} [{r.ci_low:+.5f}, {r.ci_high:+.5f}] (Holm p = {r.p_holm:.2g})")(diff("random", k))
    q.append(f"| {NM[k]} | {ce('hybrid', k)} | {ce('random', k)} | {dh} | {dr} | {T.loc[k, 'fit_s']:.2f} | {T.loc[k, 'predict_s']:.3f} |")
q2 = ["| quantity | value |", "|---|---|", f"| 90% / 95% interval coverage, additive GPR (hybrid / Random sites) | {100 * cov.loc[('hybrid', 'add'), 'cov90']:.1f}% / {100 * cov.loc[('hybrid', 'add'), 'cov95']:.1f}% and {100 * cov.loc[('random', 'add'), 'cov90']:.1f}% / {100 * cov.loc[('random', 'add'), 'cov95']:.1f}% |",
      f"| 90% / 95% interval coverage, one-lengthscale GPR (hybrid / Random sites) | {100 * cov.loc[('hybrid', 'iso'), 'cov90']:.1f}% / {100 * cov.loc[('hybrid', 'iso'), 'cov95']:.1f}% and {100 * cov.loc[('random', 'iso'), 'cov90']:.1f}% / {100 * cov.loc[('random', 'iso'), 'cov95']:.1f}% |",
      f"| SVR, wider grid (FPGAs 27, 29): hybrid / Random sites | {GW.loc[('hybrid', 'svr_wide'), 'mean']:.4f} / {GW.loc[('random', 'svr_wide'), 'mean']:.4f} (additive GPR {GW.loc[('hybrid', 'add'), 'mean']:.4f} / {GW.loc[('random', 'add'), 'mean']:.4f}) |",
      f"| MLP, wider grid (FPGAs 27, 29): hybrid / Random sites | {GW.loc[('hybrid', 'mlp_wide'), 'mean']:.4f} / {GW.loc[('random', 'mlp_wide'), 'mean']:.4f} (additive GPR {GW.loc[('hybrid', 'add'), 'mean']:.4f} / {GW.loc[('random', 'add'), 'mean']:.4f}) |"]
s = open(f"{H}/summary.md", encoding="utf8").read(); s = s.split("\n## 6. Predictor comparison")[0].rstrip()
s += f"""

## 6. Predictor comparison on the frozen sites (Stage 6; all {nfp} TEST FPGAs)
Sites: the first 317 sites of the saved FPGA-1 orders, 5 seeds, for the frozen hybrid and for Random; same sites for every predictor; held-out 30% test set; error = |mu − v| / d_spec averaged test sites → path → FPGA → mean over the {nfp} FPGAs; mean ± SD over the 5 seeds. Paired differences: path-level (32 paths, mean over seeds and FPGAs), 95% bootstrap CI, Wilcoxon, Holm over 16 comparisons (10 versus the additive GPR, 6 site contrasts). Full report: `stage6_report.md`; wider-grid report: `stage6_wide_report.md`.

### 6.1 Hyperparameter tuning and cross-validation actually used
{chr(10).join(cv)}

Tuning was done separately for every FPGA, path, seed and site type, on the 317 training sites only (no test site is used). The grids are my choice.

### 6.2 Error
{chr(10).join(rt)}

### 6.3 Paired differences
{chr(10).join(dt)}

Effect of the site choice for each predictor:

{chr(10).join(st)}

### 6.4 Wider-grid check (FPGAs 27 and 29 only, 32 paths, 5 seeds)
{chr(10).join(wt)}

{chr(10).join(wdt)}

With the wider grids the order additive GPR < SVR < MLP is unchanged on Random sites; on hybrid sites SVR and the additive GPR are tied (difference −0.00008, Holm p = 0.28); the MLP remains last (+0.0027 on hybrid sites, +0.0032 on Random sites above the additive GPR). The wider SVR grid is not a strict superset of the Stage 6 grid (gamma 0.01 and epsilon 0.3 were dropped; the best epsilon is at the upper grid value 0.1 in 56% of fits). The MLP optimum is at a wider-grid edge in about half of the fits (width 128: 34%, alpha 0.1: 28%; relu chosen in 85%).

### 6.5 Claims about the predictors
**Supported:** (1) The random forest has the lowest error on both site types (−0.00124 on hybrid sites, −0.00093 on Random sites against the additive GPR; Holm p = 0.0012 and 0.00074; about 1.2–1.6% of the error). (2) SVR and the additive GPR are not separable on the hybrid sites; on Random sites the additive GPR is better by 0.00059 (Holm p = 0.0032); this does not change with the wider grids (checked on 2 FPGAs). (3) The additive structure matters: the one-lengthscale GPR is 0.0039–0.0046 worse than the additive GPR; the MLP and the deep GP are 0.003–0.005 worse, on 31–32 of 32 paths (Holm p ≤ 1.1e-8). (4) The GPR predictive intervals are close to nominal (90.8–91.1% and 94.4–94.6%). (5) Only the additive GPR (+0.00059) and the deep GP (−0.00065) show a significant difference between hybrid and Random sites; no predictor shows a hybrid advantage except the deep GP, which is the worst predictor.
**Not supported:** (1) That the additive GPR is the most accurate predictor (the random forest is better). (2) That the ranking of predictors is independent of the site choice in every detail (SVR and the additive GPR swap places). (3) That the MLP or the deep GP can match the additive GPR, even with wider MLP grids. (4) That the wider-grid result holds on all 24 FPGAs (it was run on 2). (5) That the hybrid sites are better than Random sites for the best predictor (random forest: +0.00028, Holm p = 0.58).

## 7. Exact numbers to quote in the paper (Stage 6; TEST FPGAs, {nfp} devices, 317 sites, held-out test set)
{chr(10).join(q)}

{chr(10).join(q2)}

Notes for quoting: errors are mean ± SD over 5 seeds of the FPGA-averaged error; differences are path-level (32 paths) with 95% bootstrap CIs and Holm-adjusted Wilcoxon p-values (16 comparisons); times are for one thread, sequential, median of 3 repeats on FPGA 27 (8 paths) and include cross-validation where applicable (the deep GP fits and predicts in one call, its time is listed under fit); the wider-grid rows use FPGAs 27 and 29 only.
"""
open(f"{H}/summary.md", "w", encoding="utf8").write(s); print(len(s.splitlines()), "lines")
