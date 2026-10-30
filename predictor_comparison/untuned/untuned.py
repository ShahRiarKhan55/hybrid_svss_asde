"""Hybrid sites only. Setup A: untuned SVR/RF/MLP + paper GPR. Setup B: GPR kernel variants (train-only selection).
Usage: python untuned.py   (parallel over FPGAs, resumable, then writes the CSVs and prints the table)"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import time
from multiprocessing import Pool
import GPy, numpy as np, pandas as pd
from scipy.stats import wilcoxon
from sklearn.compose import TransformedTargetRegressor as TTR
from sklearn.ensemble import RandomForestRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
import run, stages  # stages: same hybrid sites; run: load_dat, gpr (paper function), same scaling recipe

R = f"{HERE}/results"
SK = {"svr": SVR(), "rf": RandomForestRegressor(random_state=0), "mlp": MLPRegressor(max_iter=2000, random_state=0)}
KERN = {  # GPy kernels, same fit recipe as the paper's gpr(): normalizer=True, optimize(max_iters=200)
    "rbf_iso": lambda: GPy.kern.RBF(2), "rbf_ard": lambda: GPy.kern.RBF(2, ARD=True),
    "m32_iso": lambda: GPy.kern.Matern32(2), "m32_ard": lambda: GPy.kern.Matern32(2, ARD=True),
    "rbf_iso_r5": lambda: GPy.kern.RBF(2),
}


def sk_fit(est, tr, te):  # same scaling as before: X and y standardised on the training sites
    m = TTR(make_pipeline(StandardScaler(), est), transformer=StandardScaler())
    return m.fit(tr[["X", "Y"]].values, tr.Value.values).predict(te[["X", "Y"]].values)


def gp_fit(name, tr, te):
    m = GPy.models.GPRegression(tr[["X", "Y"]].values, tr[["Value"]].values, KERN[name](), normalizer=True)
    if name.endswith("_r5"):  # run 0 starts from the default init, runs 1-5 from random inits; best kept
        m.optimize_restarts(5, robust=True, verbose=False, max_iters=200)
    else:
        m.optimize(max_iters=200, messages=False)
    return m.predict(te[["X", "Y"]].values)[0].ravel(), m.log_likelihood()  # LL on training sites only


def job(fpga):
    f = f"{R}/raw_fpga{fpga}.csv"
    if os.path.exists(f):
        return
    rows = []
    for path in range(1, 33):
        df = run.load_dat(fpga, path); dspec = df.Value.max() - df.Value.min()
        mk = stages.sites(fpga, path, df)["hybrid"]; tr, te = df[mk], df[~mk]
        d = lambda p: np.mean(np.abs(p - te.Value.values)) / dspec
        t = time.time(); rows.append((fpga, path, "gpr", d(run.gpr(tr, te)[0]), np.nan, time.time() - t))
        for n, est in SK.items():
            t = time.time(); rows.append((fpga, path, n, d(sk_fit(est, tr, te)), np.nan, time.time() - t))
        for n in KERN:
            t = time.time(); p, ll = gp_fit(n, tr, te); rows.append((fpga, path, n, d(p), ll, time.time() - t))
    pd.DataFrame(rows, columns=["FPGA", "PATH", "model", "Delta", "TrainLL", "sec"]).to_csv(f + ".tmp", index=False)
    os.replace(f + ".tmp", f); print(time.strftime("%H:%M:%S"), f"FPGA-{fpga} done", flush=True)


def fmt(g):  # per-FPGA mean, then mean over FPGAs (as in the existing delta CSVs)
    return g.groupby("FPGA").Delta.mean().mean()


if __name__ == "__main__":
    os.makedirs(R, exist_ok=True)
    with Pool(16) as p:
        p.map(job, [f for f in range(1, 51) if not os.path.exists(f"{R}/raw_fpga{f}.csv")], chunksize=1)
    raw = pd.concat(pd.read_csv(f"{R}/raw_fpga{i}.csv") for i in range(1, 51))
    w = raw.pivot(index=["FPGA", "PATH"], columns="model", values="Delta")
    ll = raw.pivot(index=["FPGA", "PATH"], columns="model", values="TrainLL")[list(KERN)]
    gps = ll.columns
    sel = w[gps].lookup(w.index, ll.idxmax(axis=1)) if hasattr(w, "lookup") else np.array([w.loc[i, c] for i, c in ll.idxmax(axis=1).items()])
    w["gp_best_LL_per_path"] = sel  # per-path pick by training log-likelihood (no test values)
    for m in w:  # existing CSV formats
        g = w[m].rename("Delta").reset_index(); g.to_csv(f"{R}/{m}_hybrid_delta_per_path.csv", index=False)
        pf = g.groupby("FPGA").Delta.mean().reset_index(); pf.to_csv(f"{R}/{m}_hybrid_delta_per_fpga.csv", index=False)
        pd.DataFrame({"Overall Delta": [pf.Delta.mean()]}).to_csv(f"{R}/{m}_hybrid_delta_overall.csv", index=False)
    old = pd.read_csv(f"{os.path.dirname(HERE)}/results/final/summary.csv").query("set=='all' and cond=='hybrid'").set_index("model")
    best = ll.mean().idxmax()
    out = [(m, "untuned" if m in ("gpr", *SK) else "GPR variant", w[m].reset_index().rename(columns={m: "Delta"}).pipe(fmt),
            w[m].std(), raw[raw.model == m].sec.mean(), wilcoxon(w[m], w["gpr"]).pvalue if m != "gpr" else np.nan,
            ll[m].mean() if m in ll else np.nan) for m in w if m != "gp_best_LL_per_path"]
    out.append(("gp_best_LL_per_path", "GPR selected by train LL", fmt(w["gp_best_LL_per_path"].rename("Delta").reset_index()),
                w["gp_best_LL_per_path"].std(), np.nan, wilcoxon(w["gp_best_LL_per_path"], w["gpr"]).pvalue, np.nan))
    out += [(f"{m} (tuned, earlier)", "tuned", old.loc[m, "mean_delta"], old.loc[m, "std_per_path"], old.loc[m, "sec_per_path"], np.nan, np.nan)
            for m in ("svr", "rf", "mlp", "deepgp")]
    t = pd.DataFrame(out, columns=["model", "group", "mean_delta", "std_per_path", "sec_per_path", "p_vs_gpr", "mean_train_LL"])
    t.to_csv(f"{R}/summary_table.csv", index=False)
    print(t.round(4).to_string(index=False))
    print("best mean train LL variant:", best, "| picked per path:", ll.idxmax(axis=1).value_counts().to_dict())
    print("rbf_iso vs paper gpr max |diff| per path:", (w["rbf_iso"] - w["gpr"]).abs().max())
