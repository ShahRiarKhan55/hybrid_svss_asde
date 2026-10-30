"""Parallel, resumable driver. Usage: python stages.py {1|2|sanity}
stage 1: gpr/svr/rf/mlp x {hybrid,random}; stage 2: deepgp x hybrid; one raw CSV per FPGA."""
import os
for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[v] = "1"
import sys, time
from multiprocessing import Pool
import numpy as np, pandas as pd
import run

WORKERS = 16
STAGES = {"1": (["gpr", "svr", "rf", "mlp"], ["hybrid", "random"]), "2": (["deepgp"], ["hybrid"])}


def sites(fpga, path, df):  # same hybrid / random sites for every model (as in the pilot)
    pos = run.hybrid_sites(path)
    h = np.array([(x, y) in pos for x, y in zip(df.X, df.Y)])
    r = np.zeros(len(df), bool)
    r[np.random.default_rng(fpga * 100 + path).choice(len(df), h.sum(), replace=False)] = True
    return {"hybrid": h, "random": r}


def job(a):
    out, fpga, models, conds, paths, fits = a
    f = f"{out}/raw_fpga{fpga}.csv"
    if os.path.exists(f):
        return fpga
    run.NJ = 1
    rows = []
    for path in paths:
        df = run.load_dat(fpga, path)
        dspec = df.Value.max() - df.Value.min()
        masks = sites(fpga, path, df)
        for cond in conds:
            tr, te = df[masks[cond]], df[~masks[cond]]
            for m in models:
                t = time.time()
                pred = fits.get(m, run.MODELS[m])(tr, te)
                rows.append((cond, m, fpga, path, np.mean(np.abs(pred - te.Value.values)) / dspec, time.time() - t))
    pd.DataFrame(rows, columns=["cond", "model", "FPGA", "PATH", "Delta", "sec"]).to_csv(f + ".tmp", index=False)
    os.replace(f + ".tmp", f)  # atomic: a finished FPGA is never half-written
    print(f"{time.strftime('%H:%M:%S')} FPGA-{fpga} done", flush=True)
    return fpga


if __name__ == "__main__":
    which = sys.argv[1]
    if which == "sanity":  # deepgp settings on 5 paths of FPGA 2
        from functools import partial
        paths, out = [1, 8, 16, 24, 32], "results/sanity"
        os.makedirs(out, exist_ok=True)
        res = {}
        for name, kw in [("old_50ind_200it", dict(ind=50, iters=200)), ("new_100ind_500it", dict(ind=100, iters=500))]:
            dirs = [f"{out}/{name}_p{pa}" for pa in paths]
            for d in dirs:
                os.makedirs(d, exist_ok=True)
            with Pool(5) as p:  # one path per worker
                p.map(job, [(d, 2, ["deepgp"], ["hybrid"], [pa], {"deepgp": partial(run.deepgp_fit, **kw)}) for d, pa in zip(dirs, paths)])
            res[name] = pd.concat(pd.read_csv(f"{d}/raw_fpga2.csv") for d in dirs)
        t = pd.DataFrame({k: v.set_index("PATH")[["Delta", "sec"]].stack() for k, v in res.items()}).unstack()
        print(t.to_string()); print(pd.concat(res).groupby(level=0)[["Delta", "sec"]].mean())
    else:
        models, conds = STAGES[which]
        out = f"results/stage{which}"
        os.makedirs(out, exist_ok=True)
        todo = [f for f in range(1, 51) if not os.path.exists(f"{out}/raw_fpga{f}.csv")]
        print(f"stage {which}: {len(todo)} FPGAs to do", flush=True)
        with Pool(WORKERS) as p:
            list(p.imap_unordered(job, [(out, f, models, conds, range(1, 33), {}) for f in todo]))
        print("stage done", flush=True)
