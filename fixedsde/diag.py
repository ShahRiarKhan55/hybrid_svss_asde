"""Diagnosis of the non-reproducing paths. python diag.py A|B|C  (A instrumented, B plain, C plain+sequential on the 7 paths)"""
import sys, json, os, numpy as np, pandas as pd, GPy
from multiprocessing import Pool
import run as r
hy = r.hy; H = r.H; SEEDS = hy["SEED"]; BAD = [2, 15, 16, 19, 23, 30, 31]


def fit(tr, te, log=None):
    m = GPy.models.GPRegression(tr[["X", "Y"]].values, tr[["Value"]].values, GPy.kern.RBF(input_dim=2), normalizer=True)
    o = m.optimize(max_iters=hy["MAX_ITERS"], messages=False)
    if log is not None: log.append(dict(ls=float(m.kern.lengthscale[0]), var=float(m.kern.variance[0]), noise=float(m.Gaussian_noise.variance[0]), status=str(o.status), n_iter=int(getattr(o, "funct_eval", -1))))
    if te is None: return m
    mu, v = m.predict(te[["X", "Y"]].values); return mu.flatten(), v.flatten()


def svss(df, info):
    """copy of hy.svts_select with logging of fits, boundary ties; returns same list."""
    rng = np.random.default_rng(SEEDS); n = len(df); step = max(1, round(n * 0.01)); target = max(1, round(n * hy["AL_RATIO"]))
    X = df[["X", "Y"]].values.astype(float); Y = df["Value"].values.astype(float); train = list(rng.choice(n, size=step, replace=False))
    while len(train) < target:
        remaining = [i for i in range(n) if i not in train]
        m = GPy.models.GPRegression(X[train], Y[train].reshape(-1, 1), GPy.kern.RBF(input_dim=2), normalizer=True); o = m.optimize(messages=False)
        info["svss_fits"].append(dict(ls=float(m.kern.lengthscale[0]), var=float(m.kern.variance[0]), noise=float(m.Gaussian_noise.variance[0]), status=str(o.status)))
        mean_r, var_r = m.predict(X[remaining]); mean_s, var_s = m.predict(X[train])
        am = np.concatenate([mean_r.ravel(), mean_s.ravel()]); av = np.concatenate([var_r.ravel(), var_s.ravel()])
        mz = (am - am.mean()) / (am.std() + 1e-9); vz = (av - av.mean()) / (av.std() + 1e-9); nr = len(mean_r)
        fr = np.column_stack([mz[:nr], vz[:nr]]); fs = np.column_stack([mz[nr:], vz[nr:]]); c = fs.mean(axis=0); cc = c - c.mean(); frc = fr - fr.mean(axis=1, keepdims=True)
        sim = (frc @ cc) / (np.linalg.norm(frc, axis=1) * np.linalg.norm(cc) + 1e-9); sn = (sim - sim.min()) / (sim.max() - sim.min() + 1e-9)
        vn = (var_r.ravel() - var_r.min()) / (var_r.max() - var_r.min() + 1e-9); sc = hy["SVTS_ALPHA"] * sn + hy["SVTS_BETA"] * vn
        pick = min(step, target - len(train), len(remaining)); s = np.sort(sc)[::-1]
        info["boundary_gap"].append(float(s[pick - 1] - s[pick])); info["exact_ties_at_boundary"] += int(s[pick - 1] == s[pick]); info["n_dup_scores"] += int(len(sc) - len(np.unique(sc)))
        train += [remaining[i] for i in np.argsort(sc)[::-1][:pick]]
    return train


def job(a):
    p, mode = a; df = hy["load_dat"](1, p); info = dict(svss_fits=[], boundary_gap=[], exact_ties_at_boundary=0, n_dup_scores=0, final_fits=[])
    al = svss(df, info) if mode == "A" else hy["svts_select"](df)
    if mode == "A": assert al == hy["svts_select"](df), "instrumented copy differs"
    s = hy["adaptive_sde"](df, al); _, nf = r.sde(df, al, None); pos = {(int(x), int(y)) for x, y in df.loc[al + s, ["X", "Y"]].values}; res = dict(PATH=p, fallback=nf, n_sel=len(pos), pos=sorted(pos), anchors=sorted((int(df.X[i]), int(df.Y[i])) for i in al), info=info)
    for f in (1, 2):
        d = hy["load_dat"](f, p); ds = d.Value.max() - d.Value.min(); tr = [i for i in d.index if (d.X[i], d.Y[i]) in pos]; te = [i for i in d.index if i not in set(tr)]
        mu, _ = fit(d.loc[tr], d.loc[te], info["final_fits"] if mode == "A" else None) if mode == "A" else hy["gpr"](d.loc[tr], d.loc[te])
        res[f"delta{f}"] = float(np.mean(np.abs(mu - d.Value[te].values) / ds)); res[f"ntrain{f}"] = len(tr)
    return res


if __name__ == "__main__":
    mode = sys.argv[1]; paths = BAD if mode == "C" else list(range(1, 33))
    out = [job((p, mode)) for p in paths] if mode == "C" else Pool(8).map(job, [(p, mode) for p in paths], chunksize=1)
    json.dump(out, open(f"{H}/diag/run{mode}.json", "w"))
