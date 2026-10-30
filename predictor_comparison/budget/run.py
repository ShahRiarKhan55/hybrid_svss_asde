"""Budget sweep: hybrid SVSS+ASDE vs SDE vs random at total budgets 3/5/7/10 %, additive RBF GPR, all 50 FPGAs x 32 paths.
Selection is calibrated on FPGA 1 per path (as in the existing runs) and the selection ORDER is saved: results/sel_path{p}.csv.
Stage 1 reuses the existing functions: hybrid_svss_asde.ipynb cell 0 (svts_select, adaptive_sde) and sde.ipynb cell 0 (sde_select).
Test set per (FPGA, path, budget) = sites outside ALL training sets of that budget (hybrid, SDE, 3 random seeds), identical for every strategy.
Usage: python run.py  (8 workers, resumable; selection first, then evaluation, then the analysis)"""
import os
for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[v] = "1"
import sys
HERE = os.path.dirname(os.path.abspath(__file__)); PC = os.path.dirname(HERE); ROOT = os.path.dirname(PC)
sys.path.insert(0, PC)  # `import run` is the earlier predictor_comparison/run.py (same file name as this script)
import json, random, time, types
from multiprocessing import Pool
import GPy, numpy as np, pandas as pd
from scipy.stats import wilcoxon
import run

R = f"{HERE}/results"
BUDGETS = [0.03, 0.05, 0.07, 0.10]
hy = run.hy  # hybrid_svss_asde.ipynb cell 0: load_dat, svts_select (SVSS, 3 %), adaptive_sde (ASDE)
sde = types.ModuleType("sde"); exec("".join(json.load(open(f"{ROOT}/sde.ipynb", encoding="utf8"))["cells"][0]["source"]), sde.__dict__)  # sde_select, SEED


def select(p):  # FPGA-1 selection for path p, orders saved; budgets are prefixes of the hybrid / SDE orders
    f = f"{R}/sel_path{p}.csv"
    if os.path.exists(f):
        return
    df = hy.load_dat(1, p); n = len(df)
    svss = hy.svts_select(df)  # SVSS anchors (3 %), in selection order
    hy.SDE_RATIO = 0.07; asde = hy.adaptive_sde(df, svss)  # ASDE additions, in acceptance order (a smaller target is a prefix)
    k = {b: round(n * (b - 0.03)) for b in BUDGETS}; nb = {b: len(svss) + k[b] for b in BUDGETS}
    full = sde.sde_select(df.X.values, df.Y.values, nb[0.10], random.Random(sde.SEED))
    for b in BUDGETS:  # SDE run separately per budget: must equal the prefix of the 10 % run
        assert sde.sde_select(df.X.values, df.Y.values, nb[b], random.Random(sde.SEED)) == full[:nb[b]]
    saved = pd.read_csv(f"{ROOT}/hybrid_svss_asde/FPGA1T/fpga1_path{p}_train_hybrid_svts_sde.csv")  # existing 10 % hybrid sites
    common = len(set(zip(saved.X, saved.Y)) & set(map(tuple, df.loc[svss + asde[:k[0.10]], ["X", "Y"]].values)))
    print(f"path {p}: regenerated 10% hybrid sites share {common} of {len(saved)} with the saved ones", flush=True)  # GPy numerics: not bit-reproducible on some paths
    o = pd.concat([df.loc[svss, ["X", "Y"]].assign(strategy="hybrid", phase="svss"), df.loc[asde, ["X", "Y"]].assign(strategy="hybrid", phase="asde"),
                   df.loc[full, ["X", "Y"]].assign(strategy="sde", phase="sde")])
    o["rank"] = o.groupby("strategy").cumcount() + 1
    o.assign(**{f"n_{int(b * 100)}": nb[b] for b in BUDGETS}, n_active=n).to_csv(f + ".tmp", index=False); os.replace(f + ".tmp", f)
    print(time.strftime("%H:%M:%S"), "selected path", p, flush=True)


def gp_err(tr, te, ds):  # additive RBF GPR, paper recipe, raw coordinates
    m = GPy.models.GPRegression(tr[["X", "Y"]].values, tr[["Value"]].values, GPy.kern.RBF(1, active_dims=[0]) + GPy.kern.RBF(1, active_dims=[1]), normalizer=True)
    m.optimize(max_iters=200, messages=False)
    e = np.abs(m.predict(te[["X", "Y"]].values)[0].ravel() - te.Value.values) / ds
    return e.mean(), np.sqrt(np.mean(e ** 2)), m.optimization_runs[-1].status == "Converged"


def evaluate(a):
    f, p = a
    out = f"{R}/path_f{f}_p{p}.csv"
    if os.path.exists(out):
        return
    S = pd.read_csv(f"{R}/sel_path{p}.csv"); df = run.load_dat(f, p); ds = df.Value.max() - df.Value.min()
    coords = list(zip(df.X, df.Y)); perm = {i: np.random.default_rng([i, f, p]).permutation(len(df)) for i in (1, 2, 3)}
    row = dict(FPGA=f, PATH=p)
    for b in BUDGETS:
        t = f"{int(b * 100)}"; nb = S[f"n_{t}"].iloc[0]
        hs = S[S.strategy == "hybrid"].iloc[:nb]; ss = S[S.strategy == "sde"].iloc[:nb]
        hset, sset = set(zip(hs.X, hs.Y)), set(zip(ss.X, ss.Y))
        mask = {"hybrid": np.array([c in hset for c in coords]), "sde": np.array([c in sset for c in coords])}
        for i in (1, 2, 3):
            m = np.zeros(len(df), bool); m[perm[i][:mask["hybrid"].sum()]] = True; mask[f"rand{i}"] = m  # same size as hybrid on this FPGA
        test = ~np.any(list(mask.values()), axis=0); te = df[test]; row[f"{t}_n_test"] = len(te)
        for s, mk in mask.items():
            np.random.seed(f * 100 + p)
            row[f"{t}_{s}_n"] = int(mk.sum()); row[f"{t}_{s}_delta"], row[f"{t}_{s}_rmse"], row[f"{t}_{s}_conv"] = gp_err(df[mk], te, ds)
    pd.DataFrame([row]).to_csv(out + ".tmp", index=False); os.replace(out + ".tmp", out)
    print(time.strftime("%H:%M:%S"), a, "done", flush=True)


def analyse():
    P = pd.concat(pd.read_csv(f"{R}/path_f{f}_p{p}.csv") for f in range(1, 51) for p in range(1, 33)); P.to_csv(f"{R}/per_path_all.csv", index=False)
    for b in BUDGETS:
        t = f"{int(b * 100)}"
        for m in ("delta", "rmse"): P[f"{t}_random_{m}"] = P[[f"{t}_rand{i}_{m}" for i in (1, 2, 3)]].mean(axis=1)
    subsets = {"all50": P.FPGA >= 1, "no_fpga1": P.FPGA >= 2}; rng, S, X = np.random.default_rng(0), [], []
    for tag, msk in subsets.items():
        Q = P[msk]
        for b in BUDGETS:
            t = f"{int(b * 100)}"
            for s in ("hybrid", "sde", "random"):
                S.append(dict(fpgas=tag, budget=f"{t}%", strategy=s, delta_bar=Q.groupby("FPGA")[f"{t}_{s}_delta"].mean().mean(), rmse=Q.groupby("FPGA")[f"{t}_{s}_rmse"].mean().mean(),
                              mean_n_train=Q[f"{t}_{'rand1' if s == 'random' else s}_n"].mean(), mean_n_test=Q[f"{t}_n_test"].mean()))
            V = {s: Q.groupby("PATH")[f"{t}_{s}_delta"].mean() for s in ("hybrid", "sde", "random")}
            for o in ("sde", "random"):
                d = (V["hybrid"] - V[o]).values; bs = [d[rng.integers(0, 32, 32)].mean() for _ in range(10000)]
                X.append(dict(fpgas=tag, budget=f"{t}%", comparison=f"hybrid - {o}", mean_diff=d.mean(), ci95_low=np.quantile(bs, .025), ci95_high=np.quantile(bs, .975),
                              paths_hybrid_better=int((d < 0).sum()), of=32, wilcoxon_p=wilcoxon(d).pvalue))
    S, X = pd.DataFrame(S), pd.DataFrame(X)
    X["bonf_budgets"] = X.wilcoxon_p < 0.05 / len(BUDGETS)  # 4 budgets per comparison and FPGA set
    X["bonf_all"] = X.wilcoxon_p < 0.05 / (len(X) // 2)  # 8 tests per FPGA set (4 budgets x 2 comparisons)
    S.to_csv(f"{R}/summary_table.csv", index=False); X.to_csv(f"{R}/path_level_tests.csv", index=False)
    pd.set_option("display.width", 250); print(S.round(4).to_string(index=False)); print(); print(X.round(5).to_string(index=False)); print()
    cv = [c for c in P if c.endswith("_conv")]
    print("non-converged fits:", int((~P[cv].astype(bool)).sum().sum()), "of", P[cv].size, "| NaN deltas:", int(P.filter(like="_delta").isna().sum().sum()))
    print("train-size differences vs hybrid (paths): sde", {b: int((P[f"{int(b * 100)}_sde_n"] != P[f"{int(b * 100)}_hybrid_n"]).sum()) for b in BUDGETS})


if __name__ == "__main__":
    os.makedirs(R, exist_ok=True)
    with Pool(8) as pool:
        pool.map(select, [p for p in range(1, 33) if not os.path.exists(f"{R}/sel_path{p}.csv")], chunksize=1)
        pool.map(evaluate, [(f, p) for f in range(1, 51) for p in range(1, 33) if not os.path.exists(f"{R}/path_f{f}_p{p}.csv")], chunksize=1)
    analyse()
