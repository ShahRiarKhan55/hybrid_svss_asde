"""Audit of run.py / stages.py (read-only): shared indices, leakage, FPGA1 path1 deltas, GPR fit details."""
import os
for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS"): os.environ[v] = "2"
import numpy as np, pandas as pd, GPy
import run, stages

fpga, path = 1, 1
df = run.load_dat(fpga, path); dspec = df.Value.max() - df.Value.min()
print("unique sites:", not df.duplicated(["X", "Y"]).any(), "| n =", len(df))
masks = stages.sites(fpga, path, df)
print("hybrid n_train =", masks["hybrid"].sum(), "| random n_train =", masks["random"].sum())

# (1) every model sees identical train/test indices; split is a partition
for cond, mk in masks.items():
    seen = {}
    for name, fit in run.MODELS.items():
        def spy(tr, te, fit=fit, name=name):
            seen[name] = (tuple(tr.index), tuple(te.index)); return fit(tr, te)
        tr, te = df[mk], df[~mk]
        pred = spy(tr, te)
        seen[name] += (np.mean(np.abs(pred - te.Value.values)) / dspec,)
    ref = next(iter(seen.values()))
    assert all(v[0] == ref[0] and v[1] == ref[1] for v in seen.values()), "index mismatch"
    assert not set(ref[0]) & set(ref[1]) and len(ref[0]) + len(ref[1]) == len(df)
    print(f"[{cond}] OK: {len(seen)} models share train({len(ref[0])})/test({len(ref[1])}), disjoint, cover all sites")
    print("   delta:", {k: round(v[2], 4) for k, v in seen.items()})
saved = pd.read_csv(f"{run.ROOT}/hybrid_svss_asde/hybrid_svts_sde_delta_per_path.csv").query("FPGA==1 and PATH==1").Delta.iloc[0]
print("existing hybrid GPR FPGA1 path1 delta:", round(saved, 4))

# (2)/(4) grid-search details and GPR fit, hybrid sites
mk = masks["hybrid"]; tr, te = df[mk], df[~mk]
for name in ("svr", "rf", "mlp"):
    c = dict(zip(run.MODELS[name].__code__.co_freevars, [x.cell_contents for x in run.MODELS[name].__closure__]))
    from sklearn.compose import TransformedTargetRegressor as TTR
    from sklearn.model_selection import GridSearchCV
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    g = GridSearchCV(TTR(make_pipeline(StandardScaler(), c["est"]), transformer=StandardScaler()),
                     {f"regressor__{k}": v for k, v in c["grid"].items()}, cv=3, n_jobs=1).fit(tr[["X", "Y"]].values, tr.Value.values)
    print(name, "best:", {k.split("__", 1)[1]: v for k, v in g.best_params_.items()}, "| grid:", c["grid"], "| cv R2 =", round(g.best_score_, 3))
m = GPy.models.GPRegression(tr[["X", "Y"]].values, tr[["Value"]].values, GPy.kern.RBF(input_dim=2), normalizer=True)
print("GPR init: variance=%.2f lengthscale=%.2f noise=%.2f ARD=%s" % (m.kern.variance[0], m.kern.lengthscale[0], m.Gaussian_noise.variance[0], m.kern.ARD))
m.optimize(max_iters=200, messages=False)
print("GPR fitted:", m.kern.variance[0], m.kern.lengthscale[0], m.Gaussian_noise.variance[0])
print("optimizer:", m.optimization_runs[-1].status, "| f_eval:", m.optimization_runs[-1].funct_eval, "| restarts: none")
mu = m.predict(te[["X", "Y"]].values)[0].ravel()
print("same prediction as run.gpr:", np.allclose(mu, run.gpr(tr, te)[0], atol=1e-6))
print("X range:", df.X.min(), df.X.max(), "Y range:", df.Y.min(), df.Y.max())
# stage-1 / pilot consistency for this path
for f in ("results/stage1/raw_fpga1.csv", "pilot/raw.csv"):
    r = pd.read_csv(f); r = r[(r.FPGA == 1) & (r.PATH == 1) & (r.cond == "hybrid")]
    print(f, r.set_index("model").Delta.round(4).to_dict())
# how often does GPR stop at the iteration cap? (path 1..8, FPGA 1 hybrid)
st = []
for p in range(1, 9):
    d = run.load_dat(1, p); k = np.array([(x, y) in run.hybrid_sites(p) for x, y in zip(d.X, d.Y)])
    mm = GPy.models.GPRegression(d[k][["X", "Y"]].values, d[k][["Value"]].values, GPy.kern.RBF(2), normalizer=True)
    mm.optimize(max_iters=200); st.append((p, mm.optimization_runs[-1].status, round(mm.kern.lengthscale[0], 2)))
print("GPR status/lengthscale paths 1-8:", st)
