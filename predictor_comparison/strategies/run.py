"""Site-selection strategies under one predictor: additive RBF GPR (reference: paper's one-lengthscale GPR), raw coordinates.
Sites come from the existing runs (FPGA1T train files, not reselected): hybrid SVSS+ASDE, SDE, SVTS, AL (raw_al); random = 3 seeds, same size as hybrid.
Usage: python run.py  (8 workers, resumable per path, then the analysis; writes strategies/results/)"""
import os
for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[v] = "1"
import sys
HERE = os.path.dirname(os.path.abspath(__file__)); PC = os.path.dirname(HERE); ROOT = os.path.dirname(PC)
sys.path.insert(0, PC)  # `import run` is the earlier predictor_comparison/run.py (same file name as this script)
import time
from multiprocessing import Pool
import GPy, numpy as np, pandas as pd
from scipy.stats import wilcoxon
import run, stages

R = f"{HERE}/results"
FILES = {"hybrid": ("hybrid_svss_asde", "hybrid_svts_sde"), "sde": ("sde", "sde"), "svts": ("svts", "svts"), "al": ("raw_al", "raw_al")}
SETS = list(FILES) + ["rand1", "rand2", "rand3"]
KERN = {"add": lambda: GPy.kern.RBF(1, active_dims=[0]) + GPy.kern.RBF(1, active_dims=[1]), "iso": lambda: GPy.kern.RBF(2)}
STRATS = ["hybrid", "sde", "svts", "al", "random"]


def positions(strategy, path):  # saved training sites of the existing run (calibrated on FPGA 1, reused for every FPGA)
    d, s = FILES[strategy]
    f = f"{ROOT}/{d}/FPGA1T/fpga1_path{path}_train_{s}.csv"
    if not os.path.exists(f):
        return None
    t = pd.read_csv(f); return set(zip(t.X, t.Y))


def gp(kind, tr, te, ds):  # paper recipe: RBF, normalizer=True, optimize(max_iters=200), raw coordinates
    m = GPy.models.GPRegression(tr[["X", "Y"]].values, tr[["Value"]].values, KERN[kind](), normalizer=True)
    m.optimize(max_iters=200, messages=False)
    e = np.abs(m.predict(te[["X", "Y"]].values)[0].ravel() - te.Value.values) / ds
    return e, m.optimization_runs[-1].status == "Converged"


def job(a):
    f, p = a
    if os.path.exists(f"{R}/path_f{f}_p{p}.csv"):
        return
    df = run.load_dat(f, p); ds = df.Value.max() - df.Value.min()  # 3-sigma filter + d_spec as in all earlier runs
    masks = {}
    for s in FILES:
        pos = positions(s, p)
        if pos is not None: masks[s] = np.array([(x, y) in pos for x, y in zip(df.X, df.Y)])
    for i in (1, 2, 3):  # random sites: same count as hybrid, seeds 1-3
        m = np.zeros(len(df), bool); m[np.random.default_rng([i, f, p]).choice(len(df), masks["hybrid"].sum(), replace=False)] = True; masks[f"rand{i}"] = m
    row, errs = dict(FPGA=f, PATH=p), {}
    for s, mk in masks.items():
        tr, te = df[mk], df[~mk]; row[f"{s}_ntrain"] = len(tr)
        if s != "hybrid" and not s.startswith("rand"): row[f"{s}_overlap_hybrid"] = (mk & masks["hybrid"]).sum() / masks["hybrid"].sum()
        for k in KERN:
            np.random.seed(f * 100 + p)
            try:
                e, ok = gp(k, tr, te, ds)
                errs[f"{s}_{k}"] = e.astype(np.float32); row[f"{s}_{k}_delta"] = e.mean(); row[f"{s}_{k}_rmse"] = np.sqrt(np.mean(e ** 2)); row[f"{s}_{k}_conv"] = ok
            except Exception as ex:
                print("FAILED", a, s, k, repr(ex)[:100], flush=True)
    np.savez_compressed(f"{R}/err_f{f}_p{p}.npz", **errs)
    pd.DataFrame([row]).to_csv(f"{R}/path_f{f}_p{p}.csv.tmp", index=False); os.replace(f"{R}/path_f{f}_p{p}.csv.tmp", f"{R}/path_f{f}_p{p}.csv")
    print(time.strftime("%H:%M:%S"), a, "done", flush=True)


def analyse():
    P = pd.concat(pd.read_csv(f"{R}/path_f{f}_p{p}.csv") for f in range(1, 51) for p in range(1, 33)); P.to_csv(f"{R}/per_path_all.csv", index=False)
    for k in KERN:  # random = mean over the 3 seeds
        for m in ("delta", "rmse"): P[f"random_{k}_{m}"] = P[[f"rand{i}_{k}_{m}" for i in (1, 2, 3)]].mean(axis=1)
    names = lambda s: [s] if s != "random" else ["rand1", "rand2", "rand3"]
    cols, subsets = {}, {"all50": P.FPGA >= 1, "no_fpga1": P.FPGA >= 2}
    keys = {(f, p): i for i, (f, p) in enumerate(zip(P.FPGA, P.PATH))}
    pooled = {(s, k, t): [] for s in STRATS for k in KERN for t in subsets}
    for (f, p) in keys:
        z = np.load(f"{R}/err_f{f}_p{p}.npz")
        for s in STRATS:
            for k in KERN:
                for t, msk in subsets.items():
                    if f >= (1 if t == "all50" else 2): pooled[(s, k, t)] += [z[f"{n}_{k}"] for n in names(s) if f"{n}_{k}" in z]
    rows = []
    for t, msk in subsets.items():
        Q = P[msk]
        for s in STRATS:
            for k in KERN:
                e = np.concatenate(pooled[(s, k, t)])
                rows.append(dict(fpgas=t, strategy=s, predictor=k, delta_bar=Q.groupby("FPGA")[f"{s}_{k}_delta"].mean().mean(), rmse=Q.groupby("FPGA")[f"{s}_{k}_rmse"].mean().mean(),
                                 median_abs=np.median(e), p99_abs=np.quantile(e, .99), mean_ntrain=Q[f"{'rand1' if s == 'random' else s}_ntrain"].mean()))
    S = pd.DataFrame(rows); S.to_csv(f"{R}/summary_table.csv", index=False)
    rng, rows = np.random.default_rng(0), []
    for t, msk in subsets.items():
        Q = P[msk]
        for k in KERN:
            V = {s: Q.groupby("PATH")[f"{s}_{k}_delta"].mean() for s in STRATS}  # 32 path-number values (mean over FPGAs)
            for o in STRATS[1:]:
                d = (V["hybrid"] - V[o]).values; bs = [d[rng.integers(0, 32, 32)].mean() for _ in range(10000)]
                rows.append(dict(fpgas=t, predictor=k, comparison=f"hybrid - {o}", mean_diff=d.mean(), ci95_low=np.quantile(bs, .025), ci95_high=np.quantile(bs, .975),
                                 paths_hybrid_better=int((d < 0).sum()), of=32, wilcoxon_p=wilcoxon(d).pvalue))
    X = pd.DataFrame(rows)
    X["bonf_4"] = X.wilcoxon_p < 0.05 / 4  # family = the 4 comparisons of one (FPGA set, predictor)
    X["bonf_16"] = X.wilcoxon_p < 0.05 / len(X)  # family = all comparisons in this table
    X.to_csv(f"{R}/path_level_tests.csv", index=False)
    pd.set_option("display.width", 250)
    print(S.round(4).to_string(index=False)); print(); print(X.round(5).to_string(index=False)); print()
    cnt = P[[f"{s}_ntrain" for s in SETS]]
    print("train counts per set (min/mean/max):", {s: (int(cnt[f"{s}_ntrain"].min()), round(cnt[f"{s}_ntrain"].mean(), 1), int(cnt[f"{s}_ntrain"].max())) for s in SETS})
    print("paths where the count differs from hybrid:", {s: int((cnt[f"{s}_ntrain"] != cnt["hybrid_ntrain"]).sum()) for s in SETS if s != "hybrid"})
    print("mean overlap with hybrid sites:", {s: round(P[f"{s}_overlap_hybrid"].mean(), 3) for s in FILES if s != "hybrid"})
    print("non-converged fits:", {c: int((~P[c].astype(bool)).sum()) for c in P if c.endswith("_conv") and (~P[c].astype(bool)).sum()}, "| NaN deltas:", int(P.filter(like="_delta").isna().sum().sum()))


if __name__ == "__main__":
    os.makedirs(R, exist_ok=True)
    missing = [p for p in range(1, 33) for s in FILES if positions(s, p) is None]
    print("missing site files:", len(missing), flush=True)
    with Pool(8) as pool:
        pool.map(job, [(f, p) for f in range(1, 51) for p in range(1, 33) if not os.path.exists(f"{R}/path_f{f}_p{p}.csv")], chunksize=1)
    analyse()
