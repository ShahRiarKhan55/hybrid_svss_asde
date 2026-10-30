"""Debug the PyDeepGP model. Hybrid sites, FPGAs 2,10,20,30,40,50 x paths 1,8,16,24,32 (30 paths).
Usage: python debug.py describe   -> prints/saves how the model is built (init vs learned values, one path)
       python debug.py            -> parallel, resumable run + table"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); PC = os.path.dirname(HERE)
sys.path.insert(0, PC); sys.path.insert(0, f"{PC}/why_dgpr")
import time
from multiprocessing import Pool
import GPy, deepgp, numpy as np, pandas as pd
import run, stages, why

R = f"{HERE}/results"
JOBS = [(f, p) for f in why.FPGAS for p in (1, 8, 16, 24, 32)]
MODELS = ["exact_ard", "sparse_raw", "sparse_std", "dgp1_base", "dgp1_fixed", "dgp1_fixed_lsinit", "dgp2_base", "dgp2_fixed", "dgp2_fixed_lsinit"]


def build(kind, Xs, ys, ind=50):
    K = [GPy.kern.RBF(2, ARD=True)] if kind == 1 else [GPy.kern.RBF(3, ARD=True), GPy.kern.RBF(2, ARD=True)]
    return deepgp.DeepGP([1, 2] if kind == 1 else [1, 3, 2], ys[:, None], X=Xs, kernels=K, num_inducing=ind, back_constraint=False)


def fix_inputs(m, kind):
    """make the layer that takes X use fixed inputs with variance 1e-6; returns True if it was a latent variable (a real change)"""
    lay = m.layers[0] if kind == 1 else m.layers[1]
    if not lay.uncertain_inputs:
        return False  # already plain fixed numbers (two-layer model): nothing to do
    lay.X.variance[:] = 1e-6
    lay.X.mean.fix(warning=False); lay.X.variance.fix(warning=False)
    return True


def snapshot(m, kind, sx):  # parameters worth printing
    top = m.layers[0] if kind == 1 else m.layers[1]
    d = {"top_inputs": type(top.X).__name__, "top_lengthscale_raw": np.round(top.kern.lengthscale.values * sx, 2).tolist(),
         "top_kern_var": round(float(top.kern.variance[0]), 3), "top_noise": round(float(top.likelihood.variance[0]), 4),
         "top_Kuu_sigma_mean": round(float(top.Kuu_sigma.values.mean()), 5)}
    if kind == 2:
        o = m.layers[0]
        d.update(obs_noise=round(float(o.likelihood.variance[0]), 4), obs_lengthscale=np.round(o.kern.lengthscale.values, 2).tolist(),
                 obs_kern_var=round(float(o.kern.variance[0]), 3), hidden_var_mean=round(float(o.X.variance.values.mean()), 5),
                 hidden_mean_std=np.round(o.X.mean.values.std(0), 3).tolist())
    else:
        d.update(input_var_mean=round(float(top.X.variance.values.mean()), 5), input_mean_std=np.round(top.X.mean.values.std(0), 3).tolist())
    return d


def fit(name, tr, te, seed):
    np.random.seed(seed)
    X = tr[["X", "Y"]].values.astype(float); y = tr.Value.values; Xt = te[["X", "Y"]].values.astype(float)
    mx, sx, my, sy = X.mean(0), X.std(0), y.mean(), y.std()
    if name in ("exact_ard", "sparse_raw"):
        p, obj, ls = why.fit("exact" if name == "exact_ard" else "sparse", 50, 200, tr, te, seed)
        return p, float(obj), ls, {}
    if name == "sparse_std":  # sparse GPR on the same standardised inputs the deep GP sees
        m = GPy.models.SparseGPRegression((X - mx) / sx, y[:, None], GPy.kern.RBF(2, ARD=True), num_inducing=50, normalizer=True)
        m.optimize(max_iters=200, messages=False)
        return m.predict((Xt - mx) / sx)[0].ravel(), float(m.log_likelihood()), m.kern.lengthscale.values * sx, {"noise": float(m.Gaussian_noise.variance[0])}
    kind = 1 if name.startswith("dgp1") else 2
    m = build(kind, (X - mx) / sx, (y - my) / sy)
    info = {}
    if "fixed" in name:
        info["was_latent"] = fix_inputs(m, kind)
    if "lsinit" in name:  # first layer (the kernel acting on X, Y) starts at the exact ARD GPR lengthscales
        ls = why.fit("exact", 0, 0, tr, te, seed)[2]
        (m.layers[0] if kind == 1 else m.layers[1]).kern.lengthscale = ls / sx
    info["init"] = snapshot(m, kind, sx)
    m.optimize(max_iters=200, messages=False)
    info["learned"] = snapshot(m, kind, sx)
    Z = (Xt - mx) / sx
    p = np.concatenate([m.predict(Z[i:i + 250])[0].ravel() for i in range(0, len(Z), 250)]) * sy + my
    top = m.layers[0] if kind == 1 else m.layers[1]
    return p, float(m.log_likelihood()), top.kern.lengthscale.values * sx, info


def job(a):
    f, p = a
    out = f"{R}/f{f}_p{p}.csv"
    if os.path.exists(out):
        return
    df = run.load_dat(f, p); dspec = df.Value.max() - df.Value.min()
    mk = stages.sites(f, p, df)["hybrid"]; tr, te = df[mk], df[~mk]
    rows, preds = [], {}
    for n in MODELS:
        t = time.time(); pr, obj, ls, info = fit(n, tr, te, f * 100 + p); preds[n] = pr
        rows.append((f, p, n, np.mean(np.abs(pr - te.Value.values)) / dspec, obj, ls[0], ls[1], time.time() - t,
                     info.get("noise", info.get("learned", {}).get("top_noise", np.nan)), np.corrcoef(pr, preds["exact_ard"])[0, 1]))
    pd.DataFrame(rows, columns=["FPGA", "PATH", "model", "Delta", "TrainObj", "ls_x", "ls_y", "sec", "noise", "corr_exact"]).to_csv(out + ".tmp", index=False)
    os.replace(out + ".tmp", out); print(time.strftime("%H:%M:%S"), a, "done", flush=True)


if __name__ == "__main__":
    os.makedirs(R, exist_ok=True)
    if sys.argv[1:] == ["describe"]:
        f, p = 2, 1
        df = run.load_dat(f, p); mk = stages.sites(f, p, df)["hybrid"]; tr, te = df[mk], df[~mk]
        for n in ("dgp1_base", "dgp2_base"):
            _, obj, _, info = fit(n, tr, te, 201)
            print(n, "ELBO", round(obj, 1)); [print("  ", k, v) for k, v in info.items()]
        m = build(2, np.random.randn(316, 2), np.random.randn(316))
        print([(type(l).__name__, "uncertain_inputs=%s" % l.uncertain_inputs, type(l.X).__name__, "kern dim %d" % l.kern.input_dim, "Z %s" % (l.Z.shape,)) for l in m.layers])
        sys.exit()
    with Pool(16) as pool:
        pool.map(job, [j for j in JOBS if not os.path.exists(f"{R}/f{j[0]}_p{j[1]}.csv")], chunksize=1)
    r = pd.concat(pd.read_csv(f"{R}/f{f}_p{p}.csv") for f, p in JOBS); r["ratio_y_over_x"] = r.ls_y / r.ls_x; r.to_csv(f"{R}/raw_all.csv", index=False)
    t = r.groupby("model", sort=False).agg(mean_delta=("Delta", "mean"), std=("Delta", "std"), obj=("TrainObj", "mean"), ls_x=("ls_x", "median"),
                                           ls_y=("ls_y", "median"), ratio=("ratio_y_over_x", "median"), noise=("noise", "median"), corr_exact=("corr_exact", "mean"), sec=("sec", "mean"))
    t.to_csv(f"{R}/table.csv"); print(t.round(4).to_string())
