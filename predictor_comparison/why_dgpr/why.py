"""Why is the deep GP not better than ARD GPR? Hybrid sites, FPGAs 2,10,20,30,40,50 x 32 paths, RBF ARD everywhere.
Usage: python why.py   (parallel over (FPGA, path), resumable, then writes the CSVs and the table)"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); PC = os.path.dirname(HERE)
sys.path.insert(0, PC)
import time
from multiprocessing import Pool
import GPy, deepgp, numpy as np, pandas as pd
import run, stages

R = f"{HERE}/results"
FPGAS = [2, 10, 20, 30, 40, 50]
# name -> (kind, inducing (None = all training sites), iterations)
FITS = {"sparse_gpr": ("sparse", 50, 200), "dgp1_ard_50_200": ("dgp1", 50, 200), "dgp2_ard_50_200": ("dgp2", 50, 200),
        "dgp2_ard_50_600": ("dgp2", 50, 600), "dgp2_ard_50_1500": ("dgp2", 50, 1500), "dgp2_ard_all_200": ("dgp2", None, 200)}


def fit(kind, ind, iters, tr, te, seed):
    """returns prediction, training objective (log-lik / ELBO), (ls_x, ls_y) of the kernel acting on the inputs, in raw coordinate units"""
    np.random.seed(seed)  # same init for every setting of a path (the earlier runs were unseeded)
    X = tr[["X", "Y"]].values.astype(float); y = tr.Value.values
    if kind == "exact":
        m = GPy.models.GPRegression(X, y[:, None], GPy.kern.RBF(2, ARD=True), normalizer=True)
        m.optimize(max_iters=200, messages=False)
        return m.predict(te[["X", "Y"]].values)[0].ravel(), m.log_likelihood(), m.kern.lengthscale.values
    if kind == "sparse":
        m = GPy.models.SparseGPRegression(X, y[:, None], GPy.kern.RBF(2, ARD=True), num_inducing=ind, normalizer=True)
        m.optimize(max_iters=iters, messages=False)
        return m.predict(te[["X", "Y"]].values)[0].ravel(), m.log_likelihood(), m.kern.lengthscale.values
    mx, sx, my, sy = X.mean(0), X.std(0), y.mean(), y.std()  # deep GP: standardised inputs and target (train only)
    kern = [GPy.kern.RBF(2, ARD=True)] if kind == "dgp1" else [GPy.kern.RBF(3, ARD=True), GPy.kern.RBF(2, ARD=True)]
    m = deepgp.DeepGP([1, 2] if kind == "dgp1" else [1, 3, 2], ((y - my) / sy)[:, None], X=(X - mx) / sx, kernels=kern,
                      num_inducing=ind or len(tr), back_constraint=False)
    m.optimize(max_iters=iters, messages=False)
    lay = [l for l in m.layers if l.kern.input_dim == 2][0]  # the layer that sees X and Y
    Xt = (te[["X", "Y"]].values - mx) / sx  # predict in chunks: 316 inducing points need N*M*M floats (2+ GB per array) otherwise
    pred = np.concatenate([m.predict(Xt[i:i + 250])[0].ravel() for i in range(0, len(Xt), 250)])
    return pred * sy + my, m.log_likelihood(), lay.kern.lengthscale.values * sx


def job(a):
    try:
        return job_(a)
    except Exception as e:  # one failed path must not kill the pool; rerun the script to retry it
        print("FAILED", a, repr(e)[:200], flush=True)


def job_(a):
    fpga, path = a
    f = f"{R}/f{fpga}_p{path}.csv"
    if os.path.exists(f):
        return
    df = run.load_dat(fpga, path); dspec = df.Value.max() - df.Value.min()
    mk = stages.sites(fpga, path, df)["hybrid"]; tr, te = df[mk], df[~mk]
    rows = []
    for name, (kind, ind, iters) in {"gpr_ard": ("exact", 0, 0), **FITS}.items():
        t = time.time(); p, obj, ls = fit(kind, ind, iters, tr, te, fpga * 100 + path)
        rows.append((fpga, path, name, np.mean(np.abs(p - te.Value.values)) / dspec, obj, ls[0], ls[1], time.time() - t))
    pd.DataFrame(rows, columns=["FPGA", "PATH", "model", "Delta", "TrainObj", "ls_x", "ls_y", "sec"]).to_csv(f + ".tmp", index=False)
    os.replace(f + ".tmp", f); print(time.strftime("%H:%M:%S"), f"FPGA-{fpga} path-{path} done", flush=True)


def table():
    raw = pd.concat(pd.read_csv(f"{R}/f{f}_p{p}.csv") for f in FPGAS for p in range(1, 33)); raw["TrainObj"] = raw.TrainObj.astype(str).str.strip("[] ").astype(float)  # sparse GPy returns a 1x1 array
    raw["ratio_y_over_x"] = raw.ls_y / raw.ls_x
    ref = lambda f, m, c=None: (lambda d: d[d.FPGA.isin(FPGAS) & (d.model == m)].assign(model=c or m))(pd.read_csv(f))
    refs = [ref(f"{PC}/results/final/raw_all.csv", "deepgp", "earlier_dgp2_iso_50_200").query("cond=='hybrid'"),
            ref(f"{PC}/dgpr_rbf/results/raw_all.csv", "dgp_ard", "earlier_dgp2_ard_50_200_unseeded")]
    allr = pd.concat([raw] + refs)
    num = ["TrainObj", "ls_x", "ls_y", "ratio_y_over_x"]  # absent in the reference runs -> NaN
    allr[num] = allr[num].apply(pd.to_numeric, errors="coerce")
    raw.to_csv(f"{R}/raw_all.csv", index=False)
    for m, g in allr.groupby("model"):  # existing delta-CSV formats
        g[["FPGA", "PATH", "Delta"]].to_csv(f"{R}/{m}_hybrid_delta_per_path.csv", index=False)
        pf = g.groupby("FPGA").Delta.mean().reset_index(); pf.to_csv(f"{R}/{m}_hybrid_delta_per_fpga.csv", index=False)
        pd.DataFrame({"Overall Delta": [pf.Delta.mean()]}).to_csv(f"{R}/{m}_hybrid_delta_overall.csv", index=False)
    t = allr.groupby("model").agg(n=("Delta", "size"), mean_delta=("Delta", "mean"), std=("Delta", "std"), sec_per_path=("sec", "mean"),
                                  train_obj=("TrainObj", "mean"), ls_x=("ls_x", "mean"), ls_y=("ls_y", "mean"), ratio_y_over_x=("ratio_y_over_x", "median"))
    t = t.reindex(["gpr_ard", "sparse_gpr", "dgp1_ard_50_200", "earlier_dgp2_iso_50_200", "earlier_dgp2_ard_50_200_unseeded",
                   "dgp2_ard_50_200", "dgp2_ard_50_600", "dgp2_ard_50_1500", "dgp2_ard_all_200"])
    t.to_csv(f"{R}/table.csv"); print(t.round(4).to_string())
    w = raw.pivot(index=["FPGA", "PATH"], columns="model", values="Delta")
    print("dgp2 200 better than gpr_ard on paths:", (w.dgp2_ard_50_200 < w.gpr_ard).mean(), "| dgp1 vs sparse corr:", w.dgp1_ard_50_200.corr(w.sparse_gpr))


if __name__ == "__main__":
    os.makedirs(R, exist_ok=True)
    with Pool(16) as p:
        p.map(job, [(f, q) for f in FPGAS for q in range(1, 33) if not os.path.exists(f"{R}/f{f}_p{q}.csv")], chunksize=1)
    table()
