"""Final check on all 50 FPGAs x 32 paths, hybrid sites: does an RBF GPR with a better kernel structure match RF?
Models: iso (paper), ard, add (RBF on X + RBF on Y), sel (best training log-lik of the three), rf (tuned, stage-1 grid), const (train mean).
Usage: python run.py   (8 workers, resumable per path, then the analysis)"""
import os
for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[v] = "1"
import sys
HERE = os.path.dirname(os.path.abspath(__file__)); PC = os.path.dirname(HERE)
sys.path.insert(0, PC)  # `import run` below is the earlier predictor_comparison/run.py (same name as this script)
import time
from multiprocessing import Pool
import GPy, numpy as np, pandas as pd
from scipy.stats import kurtosis, wilcoxon
import run, stages

R = f"{HERE}/results"
KERN = {"iso": lambda: GPy.kern.RBF(2), "ard": lambda: GPy.kern.RBF(2, ARD=True),
        "add": lambda: GPy.kern.RBF(1, active_dims=[0]) + GPy.kern.RBF(1, active_dims=[1])}
MODELS = ["iso", "ard", "add", "sel", "rf", "const"]
GP4 = ["iso", "ard", "add", "sel"]


def gp(name, tr, te, ds):  # paper recipe: normalizer=True, optimize(max_iters=200), raw coordinates
    t = time.time()
    m = GPy.models.GPRegression(tr[["X", "Y"]].values, tr[["Value"]].values, KERN[name](), normalizer=True)
    m.optimize(max_iters=200, messages=False)
    mu, var = (a.ravel() for a in m.predict(te[["X", "Y"]].values))  # predictive mean / variance of y (noise included)
    r, sd = (te.Value.values - mu) / ds, np.sqrt(var) / ds  # delta scale
    parts = m.kern.parts if name == "add" else [m.kern]
    ls = [float(l) for p in parts for l in np.atleast_1d(p.lengthscale.values)]; ls = ls * 2 if len(ls) == 1 else ls
    sx = tr[["X", "Y"]].values.std(0)
    return dict(err=np.abs(r), rmse=np.sqrt(np.mean(r ** 2)), cov90=np.mean(np.abs(r) <= 1.645 * sd), cov95=np.mean(np.abs(r) <= 1.96 * sd),
                nlpd=np.mean(0.5 * np.log(2 * np.pi * sd ** 2) + r ** 2 / (2 * sd ** 2)), LL=float(m.log_likelihood()), ls_x=ls[0] / sx[0], ls_y=ls[-1] / sx[1],
                conv=m.optimization_runs[-1].status == "Converged", sec=time.time() - t)


def job(a):
    f, p = a
    if os.path.exists(f"{R}/path_f{f}_p{p}.csv"):
        return
    run.NJ = 1
    df = run.load_dat(f, p); ds = df.Value.max() - df.Value.min()  # 3-sigma filter + d_spec as in all earlier runs
    mk = stages.sites(f, p, df)["hybrid"]; tr, te = df[mk], df[~mk]
    res = {}
    for n in KERN:
        np.random.seed(f * 100 + p)
        try:
            res[n] = gp(n, tr, te, ds)
        except Exception as e:
            print("FAILED", a, n, repr(e)[:120], flush=True)
    pick = max(res, key=lambda k: res[k]["LL"]) if res else None
    res["sel"] = dict(res[pick], sec=sum(r["sec"] for r in res.values())) if pick else None  # sel runtime = fitting all three
    t = time.time(); pr = run.MODELS["rf"](tr, te); e = np.abs(pr - te.Value.values) / ds
    res["rf"] = dict(err=e, rmse=np.sqrt(np.mean(e ** 2)), sec=time.time() - t)
    e = np.abs(te.Value.values - tr.Value.mean()) / ds; res["const"] = dict(err=e, rmse=np.sqrt(np.mean(e ** 2)), sec=0.0)
    row = dict(FPGA=f, PATH=p, n_test=len(te), picked=pick)
    for n, r in res.items():
        if r:
            row[f"{n}_delta"] = float(r["err"].mean())
            row.update({f"{n}_{k}": v for k, v in r.items() if k != "err"})
    np.savez_compressed(f"{R}/err_f{f}_p{p}.npz", **{n: r["err"].astype(np.float32) for n, r in res.items() if r})
    pd.DataFrame([row]).to_csv(f"{R}/path_f{f}_p{p}.csv.tmp", index=False); os.replace(f"{R}/path_f{f}_p{p}.csv.tmp", f"{R}/path_f{f}_p{p}.csv")
    print(time.strftime("%H:%M:%S"), a, "done", flush=True)


def analyse():
    P = pd.concat(pd.read_csv(f"{R}/path_f{f}_p{p}.csv") for f in range(1, 51) for p in range(1, 33)); P.to_csv(f"{R}/per_path_all.csv", index=False)
    E = {m: [] for m in MODELS}
    for f in range(1, 51):
        for p in range(1, 33):
            z = np.load(f"{R}/err_f{f}_p{p}.npz")
            for m in MODELS:
                if m in z: E[m].append(z[m])
    E = {m: np.concatenate(v) for m, v in E.items()}
    fm = lambda col: P.groupby("FPGA")[col].mean().mean()  # path -> FPGA mean -> mean (paper averaging)
    w = P.n_test
    rows = []
    for m in MODELS:
        rows.append(dict(model=m, delta_bar=fm(f"{m}_delta"), rmse=fm(f"{m}_rmse"), median_abs=np.median(E[m]), p99_abs=np.quantile(E[m], .99),
                         cov90=np.average(P[f"{m}_cov90"], weights=w) if m in GP4 else np.nan, cov95=np.average(P[f"{m}_cov95"], weights=w) if m in GP4 else np.nan,
                         nlpd=np.average(P[f"{m}_nlpd"], weights=w) if m in GP4 else np.nan, sec_per_path=P[f"{m}_sec"].mean(), failed_paths=int(P[f"{m}_delta"].isna().sum())))
    T = pd.DataFrame(rows).set_index("model"); T.to_csv(f"{R}/summary_table.csv"); print(T.round(4).to_string())
    share = P.picked.value_counts(normalize=True).reindex(["iso", "ard", "add"]).fillna(0); share.to_csv(f"{R}/selected_share.csv")
    print("\nselected picks (share of 1600 paths):", share.round(4).to_dict(), "| share picking ard by path number (min/max):", P.groupby("PATH").picked.apply(lambda s: (s == "ard").mean()).agg(["min", "max"]).round(2).to_dict())
    V = {m: P.groupby("PATH")[f"{m}_delta"].mean() for m in MODELS}  # 32 path-number values, each averaged over the 50 FPGAs
    rng, rows = np.random.default_rng(0), []
    for o in ["iso", "ard", "add", "rf"]:
        d = (V["sel"] - V[o]).values; bs = [d[rng.integers(0, 32, 32)].mean() for _ in range(10000)]
        rows.append(dict(comparison=f"selected - {o}", mean_diff=d.mean(), ci95_low=np.quantile(bs, .025), ci95_high=np.quantile(bs, .975),
                         paths_selected_better=int((d < 0).sum()), of=32, wilcoxon_p=wilcoxon(d).pvalue))
    X = pd.DataFrame(rows); X.to_csv(f"{R}/path_level_tests.csv", index=False); print("\n", X.round(5).to_string(index=False))
    pd.DataFrame(V).to_csv(f"{R}/delta_by_path_number.csv")
    ref = pd.read_csv(f"{PC}/../hybrid_svss_asde/hybrid_svts_sde_delta_per_path.csv").set_index(["FPGA", "PATH"]).Delta
    print("\niso vs the paper's per-path delta file: max |diff| =", (P.set_index(["FPGA", "PATH"]).iso_delta - ref).abs().max())
    for n in ["iso", "ard", "add"]:
        ext = ((P[[f"{n}_ls_x", f"{n}_ls_y"]] < .01) | (P[[f"{n}_ls_x", f"{n}_ls_y"]] > 100)).any(axis=1)
        print(f"{n}: not converged {int((~P[f'{n}_conv']).sum())}, extreme lengthscale (<0.01 or >100 train-std) {ext.sum()} ({ext.mean():.1%}), failed {int(P[f'{n}_delta'].isna().sum())}")
    print("kurtosis of pooled abs errors:", {m: round(float(kurtosis(E[m])), 2) for m in MODELS})


if __name__ == "__main__":
    os.makedirs(R, exist_ok=True)
    with Pool(8) as pool:
        pool.map(job, [(f, p) for f in range(1, 51) for p in range(1, 33) if not os.path.exists(f"{R}/path_f{f}_p{p}.csv")], chunksize=1)
    analyse()
