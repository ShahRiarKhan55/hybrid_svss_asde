"""Stage 6: predictors on the frozen hybrid sites (317 sites, 5 seeds) and on Random sites (317 sites, 5 seeds). TEST FPGAs 27-50, held-out 30% test set.
python stage6.py [run|report|timing]"""
from lib6 import *
import json, sys
from multiprocessing import Pool
from scipy.stats import wilcoxon
OUT = f"{HERE}/stage6"; os.makedirs(OUT, exist_ok=True); S1 = f"{HERE}/stage1"
PRED = ["add", "iso", "svr", "rf", "mlp", "dgp"]; NM = {"add": "Additive RBF GPR", "iso": "One-lengthscale RBF GPR", "svr": "SVR (RBF)", "rf": "Random forest", "mlp": "MLP (tanh, 2 hidden layers)", "dgp": "Deep GP [1,3,2]"}
SITES = {"hybrid": "hybrid", "random": "random"}


def timed(kind, Xtr, ytr, Xte, seed):
    best, var = "", None
    if kind in ("add", "iso"):
        t0 = time.perf_counter()
        if kind == "add": m = fit(Xtr, ytr, "p_add")
        else:
            m = GPy.models.GPRegression(Xtr, ytr.reshape(-1, 1), GPy.kern.RBF(2), normalizer=True)
            m.kern.lengthscale.constrain_bounded(*LS, warning=False); m.kern.variance.constrain_bounded(*VAR, warning=False); m.Gaussian_noise.variance.constrain_bounded(*NOISE, warning=False); m.optimize(max_iters=200, messages=False)
        t1 = time.perf_counter(); mu, var = m.predict(Xte); mu, var = mu.ravel(), var.ravel(); t2 = time.perf_counter()
    elif kind == "dgp":
        t0 = time.perf_counter(); mu = dgp_pred(Xtr, ytr, Xte, seed); t2 = time.perf_counter(); t1 = t2  # the deep GP function fits and predicts in one call
    else:
        pipe = TTR(make_pipeline(StandardScaler(), EST[kind]), transformer=StandardScaler())
        t0 = time.perf_counter(); gs = GridSearchCV(pipe, {f"regressor__{k}": v for k, v in GRIDS[kind].items()}, cv=CV(), scoring="neg_mean_absolute_error", n_jobs=1).fit(Xtr, ytr); t1 = time.perf_counter()
        mu = gs.predict(Xte); t2 = time.perf_counter(); best = json.dumps({k.split("__", 1)[1]: (list(v) if isinstance(v, tuple) else v) for k, v in gs.best_params_.items()})
    return mu, var, best, t1 - t0, t2 - t1


def job(a):
    f, p = a; out = f"{OUT}/res_f{f}_p{p}.csv"
    if os.path.exists(out): return
    df = hy.load_dat(f, p); ds = df.Value.max() - df.Value.min(); co = list(zip(df.X, df.Y)); cid = {c: i for i, c in enumerate(co)}; held = heldout(p)[1]; it = np.array([c in held for c in co])
    X, v = df[["X", "Y"]].values.astype(float), df.Value.values; rows = []
    for seed in SEEDS:
        O = pd.read_csv(f"{S1}/orders_s{seed}_p{p}.csv")
        for site, meth in SITES.items():
            g = O[O.method == meth].sort_values("rank"); ix = np.array([cid.get(c, -1) for c in zip(g.X, g.Y)][:T]); ix = ix[ix >= 0]; assert not it[ix].any()
            mk = np.zeros(len(df), bool); mk[ix] = True; Xtr, ytr, Xte, yte, tt = X[mk], v[mk], X[~mk], v[~mk], it[~mk]
            for kind in PRED:
                mu, var, best, tf, tp = timed(kind, Xtr, ytr, Xte, seed * 100000 + f * 100 + p); e = np.abs(mu - yte) / ds; r = dict(FPGA=f, PATH=p, site=site, seed=seed, predictor=kind, err=e[tt].mean(), err_old=e.mean(), n_train=int(mk.sum()), n_test=int(tt.sum()), fit_s=tf, predict_s=tp, best=best)
                if var is not None: sd = np.sqrt(var) / ds; r["cov90"] = float(np.mean(e[tt] <= 1.645 * sd[tt])); r["cov95"] = float(np.mean(e[tt] <= 1.96 * sd[tt]))
                rows.append(r)
    pd.DataFrame(rows).to_csv(out + ".tmp", index=False); os.replace(out + ".tmp", out); print(time.strftime("%H:%M:%S"), "done", a, flush=True)


def holm(p):
    p = np.asarray(p); o = np.argsort(p); m = len(p); adj = np.empty(m); run_ = 0
    for rank, i in enumerate(o): run_ = max(run_, min(1, (m - rank) * p[i])); adj[i] = run_
    return adj


def report():
    FP = [f for f in TEST if all(os.path.exists(f"{OUT}/res_f{f}_p{p}.csv") for p in range(1, 33))]  # FPGAs with all 32 paths finished
    R = pd.concat(pd.read_csv(f"{OUT}/res_f{f}_p{p}.csv") for f in FP for p in range(1, 33)); R.to_csv(f"{OUT}/res_all.csv", index=False)
    seedm = R.groupby(["site", "predictor", "seed", "FPGA"]).err.mean().groupby(["site", "predictor", "seed"]).mean(); G = seedm.groupby(["site", "predictor"]).agg(mean="mean", sd="std")
    Go = R.groupby(["site", "predictor", "seed", "FPGA"]).err_old.mean().groupby(["site", "predictor", "seed"]).mean().groupby(["site", "predictor"]).agg(mean="mean", sd="std")
    P = R.groupby(["site", "predictor", "PATH"]).err.mean(); rng = np.random.default_rng(0); rows = []

    def add(label, d, kind):
        d = np.asarray(d); lo, hi = np.quantile([d[rng.integers(0, 32, 32)].mean() for _ in range(10000)], [.025, .975])
        rows.append(dict(comparison=label, kind=kind, mean_diff=d.mean(), ci_low=lo, ci_high=hi, first_wins=int((d < 0).sum()), second_wins=int((d > 0).sum()), p_raw=wilcoxon(d).pvalue))
    for s in SITES:
        for k in PRED[1:]: add(f"{NM[k]} − additive GPR, {s} sites", P[s, k].values - P[s, "add"].values, "vs additive GPR")
    for k in PRED: add(f"{NM[k]}: hybrid sites − Random sites", P["hybrid", k].values - P["random", k].values, "site effect")
    D = pd.DataFrame(rows); D["p_holm"] = holm(D.p_raw.values); D["significant"] = D.p_holm < 0.05; D.to_csv(f"{OUT}/tests.csv", index=False)
    tab = lambda s, g=G: ["| predictor | error (mean ± SD over 5 seeds) | rank |", "|---|---|---|"] + [f"| {NM[k]} | {g.loc[(s, k), 'mean']:.4f} ± {g.loc[(s, k), 'sd']:.4f} | {int(g.loc[s]['mean'].rank().loc[k])} |" for k in PRED]
    dt = lambda kind: [f"| comparison (negative = first lower error) | mean diff | 95% CI | first / second wins (of 32) | Wilcoxon p | Holm p ({len(D)} tests) | significant |", "|---|---|---|---|---|---|---|"] + [f"| {r.comparison} | {r.mean_diff:+.5f} | [{r.ci_low:+.5f}, {r.ci_high:+.5f}] | {r.first_wins} / {r.second_wins} | {r.p_raw:.2g} | {r.p_holm:.2g} | {'yes' if r.significant else 'no'} |" for _, r in D[D.kind == kind].iterrows()]
    tm = R.groupby(["site", "predictor"])[["fit_s", "predict_s"]].mean(); tmt = ["| predictor | fit incl. tuning, s | predict, s |", "|---|---|---|"] + [f"| {NM[k]} | {tm.loc[('hybrid', k), 'fit_s']:.2f} | {tm.loc[('hybrid', k), 'predict_s']:.3f} |" for k in PRED]
    cv = R[R.predictor.isin(["add", "iso"])].groupby(["site", "predictor"])[["cov90", "cov95"]].mean(); ct = ["| GPR | sites | 90% interval coverage | 95% interval coverage |", "|---|---|---|---|"] + [f"| {NM[k]} | {s} | {100 * cv.loc[(s, k), 'cov90']:.1f}% | {100 * cv.loc[(s, k), 'cov95']:.1f}% |" for s in SITES for k in ["add", "iso"]]
    # grid-edge frequency of tuned hyperparameters
    edge = ["| predictor | hyperparameter | share of fits with the best value at a grid edge |", "|---|---|---|"]
    for k in ["svr", "rf", "mlp"]:
        B = R[R.predictor == k].best.apply(json.loads)
        for hp, vals in GRIDS[k].items():
            name = hp; vals = [list(x) if isinstance(x, tuple) else x for x in vals]; lo, hi = vals[0], vals[-1]
            edge.append(f"| {NM[k]} | {hp.split("__", 1)[1]} (grid {vals[0]} … {vals[-1]}) | {100 * B.apply(lambda b: b[name] in (lo, hi)).mean():.0f}% |")
    nt = R.groupby("site").agg(n_train=("n_train", "mean"), n_test=("n_test", "mean"))
    r_h, r_r = G.loc["hybrid"]["mean"].rank(), G.loc["random"]["mean"].rank(); same = bool((r_h == r_r).all())
    md = f"""# Stage 6 – Predictor comparison on the frozen sites (TEST FPGAs: {len(FP)} of the 24)

## Protocol
- Sites: the saved FPGA-1 selection orders, first 317 sites, 5 seeds: **hybrid** (the frozen current hybrid) and **Random**. Same sites for every predictor. Mean training sites per FPGA {nt.n_train.round(1).to_dict()}.
- Test: the held-out 30% sites (about {nt.n_test.mean():.0f} per path), {len(FP)} of the 24 TEST FPGAs ({', '.join(map(str, FP))}); error = |mu − v| / d_spec averaged test sites → path → FPGA → mean over the {len(FP)} FPGAs; mean ± SD over the 5 seeds.
- Predictors: additive RBF GPR (the bounded GPR used in all stages); one-lengthscale RBF GPR (same normaliser, bounds and 200 iterations); SVR (RBF; C ∈ {{0.1, 1, 10, 100}}, gamma ∈ {{0.01, 0.1, 1, 10}}, epsilon ∈ {{0.01, 0.1, 0.3}}; scaled inputs, standardised target); random forest (200 trees; min_samples_leaf ∈ {{1, 3, 5, 10}}, max_features ∈ {{1.0, 0.5}}); MLP (tanh, two hidden layers of width w, w ∈ {{16, 32, 64}}, alpha ∈ {{0.001, 0.1, 10}}, Adam, up to 1000 iterations); deep GP (layers [1,3,2], 50 inducing points, 200 iterations, isotropic RBF kernels, as in `predictor_comparison`). **Assumptions:** the grids above are my choice; SVR, RF and MLP are tuned by 5-fold CV (shuffled, fixed seed) with mean absolute error as the score, on the training sites only, separately for every FPGA, path, seed and site type, and the final model is refit on all training sites.
- Paired differences: path-level values (mean over seeds and the {len(FP)} FPGAs), 32 paths, 95% bootstrap CI (10,000 resamples, seed 0), Wilcoxon signed-rank (two-sided); Holm correction over all {len(D)} comparisons below ({len(D[D.kind == 'vs additive GPR'])} predictor-vs-additive-GPR comparisons plus {len(D[D.kind == 'site effect'])} hybrid-vs-Random site contrasts).

## Error on the hybrid sites
{chr(10).join(tab("hybrid"))}

## Error on Random sites
{chr(10).join(tab("random"))}

Ranking by site type: {"the ranking of the six predictors is identical on hybrid and Random sites." if same else "the ranking differs between hybrid and Random sites (see the rank columns)."}

## Paired differences to the additive GPR
{chr(10).join(dt("vs additive GPR"))}

## Does the site choice matter for each predictor? (hybrid sites − Random sites)
{chr(10).join(dt("site effect"))}

## Interval coverage of the GPRs (held-out test sites, predictive variance including noise)
{chr(10).join(ct)}

## Fit and predict time per path (means over all fits of the main run on hybrid sites; 8 jobs ran in parallel on a laptop, so these are indicative; see the sequential timing below if present)
{chr(10).join(tmt)}

## Tuned hyperparameters at a grid edge
{chr(10).join(edge)}

## Secondary: paper's old test set (all non-training sites)
| predictor | hybrid sites | Random sites |
|---|---|---|
""" + "\n".join(f"| {NM[k]} | {Go.loc[('hybrid', k), 'mean']:.4f} ± {Go.loc[('hybrid', k), 'sd']:.4f} | {Go.loc[('random', k), 'mean']:.4f} ± {Go.loc[('random', k), 'sd']:.4f} |" for k in PRED) + "\n"
    open(f"{HERE}/stage6_report.md", "w", encoding="utf8").write(md); print("report written")


def timing():
    """sequential, one thread, nothing else started by me: FPGA 27, hybrid seed-1 sites, paths 1,5,...,29, 3 repeats, median; appended to the report"""
    rows = []; f = 27
    for rep in range(3):
        for p in range(1, 33, 4):
            df = hy.load_dat(f, p); co = list(zip(df.X, df.Y)); cid = {c: i for i, c in enumerate(co)}; O = pd.read_csv(f"{S1}/orders_s1_p{p}.csv"); g = O[O.method == "hybrid"].sort_values("rank")
            ix = np.array([cid.get(c, -1) for c in zip(g.X, g.Y)]); ix = ix[ix >= 0]; mk = np.zeros(len(df), bool); mk[ix] = True; X, v = df[["X", "Y"]].values.astype(float), df.Value.values
            for kind in PRED:
                _, _, _, tf, tp = timed(kind, X[mk], v[mk], X[~mk], 7); rows.append(dict(rep=rep, PATH=p, predictor=kind, fit_s=tf, predict_s=tp))
    T_ = pd.DataFrame(rows); T_.to_csv(f"{OUT}/sequential_timing.csv", index=False); M = T_.groupby(["predictor", "PATH"]).median().groupby("predictor").mean()
    t = ["| predictor | fit incl. tuning, s | predict on ~2,850 sites, s |", "|---|---|---|"] + [f"| {NM[k]} | {M.loc[k, 'fit_s']:.2f} | {M.loc[k, 'predict_s']:.3f} |" for k in PRED]
    add = "\n## Sequential timing (one thread, nothing else started by me, FPGA 27, hybrid seed-1 sites, 8 paths, median of 3 repeats)\n" + "\n".join(t) + "\n\nThe deep GP function fits and predicts in one call; its time is reported under fit.\n"
    md = open(f"{HERE}/stage6_report.md", encoding="utf8").read().split("\n## Sequential timing")[0].rstrip() + "\n" + add
    open(f"{HERE}/stage6_report.md", "w", encoding="utf8").write(md); print(add)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    if mode == "run":
        with Pool(8) as pool: pool.map(job, [(f, p) for f in TEST[::2] + TEST[1::2] for p in range(1, 33) if not os.path.exists(f"{OUT}/res_f{f}_p{p}.csv")], chunksize=1)  # every second TEST FPGA first
        report()
    elif mode == "report": report()
    else: timing()
