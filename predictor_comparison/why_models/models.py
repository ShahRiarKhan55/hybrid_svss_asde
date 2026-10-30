"""Fit six models on hybrid sites for FPGAs 2,10,20,30,40,50 x 32 paths; save predictions at ALL sites (train rows flagged).
Usage: python models.py  (8 workers, resumable: one pred_*/meta_* file pair per path)"""
import os
for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[v] = "1"
import sys
HERE = os.path.dirname(os.path.abspath(__file__)); PC = os.path.dirname(HERE)
sys.path.insert(0, PC); sys.path.insert(0, f"{PC}/untuned")
import time
from multiprocessing import Pool
import GPy, numpy as np, pandas as pd
import run, stages, untuned

R = f"{HERE}/results"
FPGAS = [2, 10, 20, 30, 40, 50]
KERN = {"gpr_iso": lambda: GPy.kern.RBF(2),  # the paper's GPR
        "gpr_ard": lambda: GPy.kern.RBF(2, ARD=True),
        "gpr_add": lambda: GPy.kern.RBF(1, active_dims=[0]) + GPy.kern.RBF(1, active_dims=[1])}


def gp(name, tr, ev):  # paper recipe: normalizer=True, optimize(max_iters=200); raw coordinates
    m = GPy.models.GPRegression(tr[["X", "Y"]].values, tr[["Value"]].values, KERN[name](), normalizer=True)
    m.optimize(max_iters=200, messages=False)
    parts = m.kern.parts if name == "gpr_add" else [m.kern]
    ls = [float(l) for p in parts for l in np.atleast_1d(p.lengthscale.values)]
    var = float(sum(p.variance[0] for p in parts)); noise = float(m.Gaussian_noise.variance[0])
    return m.predict(ev[["X", "Y"]].values)[0].ravel(), dict(model=name, ls1=ls[0], ls2=ls[-1], kern_var=var, noise=noise, noise_frac=noise / (noise + var), LL=float(m.log_likelihood()))


def job(a):
    f, p = a
    if os.path.exists(f"{R}/pred_f{f}_p{p}.csv"):
        return
    run.NJ = 1
    df = run.load_dat(f, p); mk = stages.sites(f, p, df)["hybrid"]; tr, te = df[mk], df[~mk]
    ev = pd.concat([tr, te])  # train rows first, then test rows
    out = ev.assign(train=[1] * len(tr) + [0] * len(te)); meta = []
    for n in KERN:
        np.random.seed(f * 100 + p); out["pred_" + n], info = gp(n, tr, ev); meta.append(info)
    out["pred_svr"] = untuned.sk_fit(untuned.SK["svr"], tr, ev)  # untuned defaults, same scaling as stage 1
    out["pred_rf"] = run.MODELS["rf"](tr, ev)   # stage-1 grid search (CV on the training sites only)
    out["pred_mlp"] = run.MODELS["mlp"](tr, ev)
    pd.DataFrame(meta).assign(FPGA=f, PATH=p, dspec=df.Value.max() - df.Value.min(), vstd=df.Value.std(), n=len(df)).to_csv(f"{R}/meta_f{f}_p{p}.csv", index=False)
    out.to_csv(f"{R}/pred_f{f}_p{p}.csv.tmp", index=False); os.replace(f"{R}/pred_f{f}_p{p}.csv.tmp", f"{R}/pred_f{f}_p{p}.csv")
    print(time.strftime("%H:%M:%S"), a, "done", flush=True)


if __name__ == "__main__":
    os.makedirs(R, exist_ok=True)
    with Pool(8) as pool:
        pool.map(job, [(f, p) for f in FPGAS for p in range(1, 33) if not os.path.exists(f"{R}/pred_f{f}_p{p}.csv")], chunksize=1)
    print("all done", flush=True)
