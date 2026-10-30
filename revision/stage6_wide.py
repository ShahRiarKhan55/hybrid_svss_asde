"""Stage 6 wider-grid check on TEST FPGAs 27 and 29: SVR and MLP with wider hyperparameter grids, same sites (hybrid, Random; 5 seeds), same test set and error.
python stage6_wide.py [run|report]"""
from lib6 import *
import json, sys
from multiprocessing import Pool
from scipy.stats import wilcoxon
OUT = f"{HERE}/stage6_wide"; os.makedirs(OUT, exist_ok=True); S1 = f"{HERE}/stage1"; FP = [27, 29]
WG = {"svr": {"svr__C": [0.1, 1, 10, 100, 1000], "svr__gamma": [0.03, 0.1, 0.3, 1, 3, 10], "svr__epsilon": [0.001, 0.01, 0.1]},
      "mlp": {"mlpregressor__hidden_layer_sizes": [(16, 16), (32, 32), (64, 64), (128, 128)], "mlpregressor__alpha": [1e-5, 1e-4, 1e-3, 1e-2, 1e-1], "mlpregressor__activation": ["tanh", "relu"]}}
WEST = {"svr": SVR(kernel="rbf"), "mlp": MLPRegressor(max_iter=1000, random_state=0)}


def wide(kind, Xtr, ytr, Xte):
    pipe = TTR(make_pipeline(StandardScaler(), WEST[kind]), transformer=StandardScaler())
    t0 = time.perf_counter(); gs = GridSearchCV(pipe, {f"regressor__{k}": v for k, v in WG[kind].items()}, cv=CV(), scoring="neg_mean_absolute_error", n_jobs=1).fit(Xtr, ytr); t1 = time.perf_counter()
    return gs.predict(Xte), json.dumps({k.split("__", 1)[1]: (list(v) if isinstance(v, tuple) else v) for k, v in gs.best_params_.items()}), t1 - t0


def job(a):
    f, p = a; out = f"{OUT}/res_f{f}_p{p}.csv"
    if os.path.exists(out): return
    df = hy.load_dat(f, p); ds = df.Value.max() - df.Value.min(); co = list(zip(df.X, df.Y)); cid = {c: i for i, c in enumerate(co)}; held = heldout(p)[1]; it = np.array([c in held for c in co])
    X, v = df[["X", "Y"]].values.astype(float), df.Value.values; rows = []
    for seed in SEEDS:
        O = pd.read_csv(f"{S1}/orders_s{seed}_p{p}.csv")
        for site in ("hybrid", "random"):
            g = O[O.method == site].sort_values("rank"); ix = np.array([cid.get(c, -1) for c in zip(g.X, g.Y)][:T]); ix = ix[ix >= 0]; mk = np.zeros(len(df), bool); mk[ix] = True
            for kind in ("svr", "mlp"):
                mu, best, tf = wide(kind, X[mk], v[mk], X[~mk]); e = np.abs(mu - v[~mk]) / ds
                rows.append(dict(FPGA=f, PATH=p, site=site, seed=seed, predictor=kind + "_wide", err=e[it[~mk]].mean(), fit_s=tf, best=best))
    pd.DataFrame(rows).to_csv(out + ".tmp", index=False); os.replace(out + ".tmp", out); print(time.strftime("%H:%M:%S"), "done", a, flush=True)


def holm(p):
    p = np.asarray(p); o = np.argsort(p); m = len(p); adj = np.empty(m); r = 0
    for rank, i in enumerate(o): r = max(r, min(1, (m - rank) * p[i])); adj[i] = r
    return adj


def report():
    W = pd.concat(pd.read_csv(f"{OUT}/res_f{f}_p{p}.csv") for f in FP for p in range(1, 33)); N = pd.read_csv(f"{HERE}/stage6/res_all.csv"); N = N[N.FPGA.isin(FP) & N.predictor.isin(["add", "svr", "mlp"])].assign(predictor=lambda d: d.predictor.replace({"svr": "svr_narrow", "mlp": "mlp_narrow"}))
    R = pd.concat([W[["FPGA", "PATH", "site", "seed", "predictor", "err", "fit_s", "best"]], N[["FPGA", "PATH", "site", "seed", "predictor", "err", "fit_s", "best"]]]); R.to_csv(f"{OUT}/res_all.csv", index=False)
    nm = {"add": "Additive RBF GPR", "svr_narrow": "SVR, Stage 6 grid", "svr_wide": "SVR, wider grid", "mlp_narrow": "MLP, Stage 6 grid (tanh)", "mlp_wide": "MLP, wider grid (tanh or relu)"}; order = list(nm)
    G = R.groupby(["site", "predictor", "seed", "FPGA"]).err.mean().groupby(["site", "predictor", "seed"]).mean().groupby(["site", "predictor"]).agg(mean="mean", sd="std")
    P = R.groupby(["site", "predictor", "PATH"]).err.mean(); rng = np.random.default_rng(0); rows = []

    def add(site, a, b, label):
        d = P[site, a].values - P[site, b].values; lo, hi = np.quantile([d[rng.integers(0, 32, 32)].mean() for _ in range(10000)], [.025, .975])
        rows.append(dict(site=site, comparison=label, mean_diff=d.mean(), ci_low=lo, ci_high=hi, first_wins=int((d < 0).sum()), second_wins=int((d > 0).sum()), p_raw=wilcoxon(d).pvalue))
    for s in ("hybrid", "random"):
        add(s, "svr_wide", "add", "SVR (wider grid) − additive GPR"); add(s, "mlp_wide", "add", "MLP (wider grid) − additive GPR"); add(s, "svr_wide", "svr_narrow", "SVR wider grid − SVR Stage 6 grid")
        add(s, "mlp_wide", "mlp_narrow", "MLP wider grid − MLP Stage 6 grid"); add(s, "svr_wide", "mlp_wide", "SVR (wider grid) − MLP (wider grid)")
    D = pd.DataFrame(rows); D["p_holm"] = holm(D.p_raw.values); D["significant"] = D.p_holm < 0.05; D.to_csv(f"{OUT}/tests.csv", index=False)
    tab = ["| predictor | hybrid sites | rank among the 3 families* | Random sites | rank* |", "|---|---|---|---|---|"]
    def rk(site, wide_):  # rank of GPR / SVR / MLP using the narrow or wide grids
        v = {"Additive GPR": G.loc[(site, "add"), "mean"], "SVR": G.loc[(site, "svr_wide" if wide_ else "svr_narrow"), "mean"], "MLP": G.loc[(site, "mlp_wide" if wide_ else "mlp_narrow"), "mean"]}
        return pd.Series(v).rank().to_dict(), pd.Series(v).sort_values().index.tolist()
    for k in order:
        fam = "Additive GPR" if k == "add" else ("SVR" if k.startswith("svr") else "MLP"); w = k.endswith("wide")
        tab.append(f"| {nm[k]} | {G.loc[('hybrid', k), 'mean']:.4f} ± {G.loc[('hybrid', k), 'sd']:.4f} | {int(rk('hybrid', w)[0][fam])} | {G.loc[('random', k), 'mean']:.4f} ± {G.loc[('random', k), 'sd']:.4f} | {int(rk('random', w)[0][fam])} |")
    dt = [f"| site | comparison (negative = first lower error) | mean diff | 95% CI | first / second wins (of 32) | Wilcoxon p | Holm p ({len(D)} tests) | significant |", "|---|---|---|---|---|---|---|---|"]
    dt += [f"| {r.site} | {r.comparison} | {r.mean_diff:+.5f} | [{r.ci_low:+.5f}, {r.ci_high:+.5f}] | {r.first_wins} / {r.second_wins} | {r.p_raw:.2g} | {r.p_holm:.2g} | {'yes' if r.significant else 'no'} |" for _, r in D.iterrows()]
    edge = ["| predictor | hyperparameter | share at a grid edge | most frequent best value |", "|---|---|---|---|"]
    for kind in ("svr", "mlp"):
        B = W[W.predictor == kind + "_wide"].best.apply(json.loads)
        for hp, vals in WG[kind].items():
            vals = [list(x) if isinstance(x, tuple) else x for x in vals]; col = B.apply(lambda b: json.dumps(b[hp])); edge.append(f"| {kind.upper()} | {hp.split('__', 1)[1]} (grid {vals[0]} … {vals[-1]}) | {100 * B.apply(lambda b: b[hp] in (vals[0], vals[-1])).mean():.0f}% | {col.value_counts().index[0]} ({100 * col.value_counts(normalize=True).iloc[0]:.0f}%) |")
    act = W[W.predictor == "mlp_wide"].best.apply(lambda b: json.loads(b)["mlpregressor__activation"]).value_counts(normalize=True)
    chg = {s: (rk(s, False)[1], rk(s, True)[1]) for s in ("hybrid", "random")}
    md = f"""# Stage 6 – wider-grid check (TEST FPGAs 27 and 29 only)

## Protocol
Same sites (hybrid and Random, 5 seeds, 317 sites), held-out test set, error definition, 5-fold CV on the training sites (mean absolute error), standardised inputs/target, FPGAs 27 and 29, all 32 paths. Wider grids: SVR with C ∈ {{0.1, 1, 10, 100, 1000}}, gamma ∈ {{0.03, 0.1, 0.3, 1, 3, 10}}, epsilon ∈ {{0.001, 0.01, 0.1}}; MLP (Adam, up to 1000 iterations) with hidden widths (16,16), (32,32), (64,64), (128,128), alpha ∈ {{1e-5, 1e-4, 1e-3, 1e-2, 1e-1}} and activation ∈ {{tanh, relu}} chosen by CV. The "Stage 6 grid" rows and the additive GPR are the Stage 6 results for the same two FPGAs. Errors are the mean over the 2 FPGAs; SD over the 5 seeds; paired differences are path-level (32 paths, mean over seeds and the 2 FPGAs) with a 95% bootstrap CI (10,000 resamples) and Wilcoxon, Holm over the {len(D)} comparisons below. Only 2 devices were used, so the numbers differ from the 24-FPGA table.

## Error
{chr(10).join(tab)}

\\* rank among the three families (additive GPR, SVR, MLP), using the Stage 6 grids for the SVR and MLP rows marked "Stage 6 grid" and the wider grids for the rows marked "wider grid".

Order of the three families (lowest error first): hybrid sites — Stage 6 grids: {' < '.join(chg['hybrid'][0])}; wider grids: {' < '.join(chg['hybrid'][1])}. Random sites — Stage 6 grids: {' < '.join(chg['random'][0])}; wider grids: {' < '.join(chg['random'][1])}.

## Paired differences
{chr(10).join(dt)}

## Wider-grid hyperparameters
{chr(10).join(edge)}

MLP activation chosen: {', '.join(f'{k} {100 * v:.0f}%' for k, v in act.items())}.

Fit time including the wider grid search (mean over fits, 8 jobs in parallel, indicative): SVR {W[W.predictor == 'svr_wide'].fit_s.mean():.1f} s, MLP {W[W.predictor == 'mlp_wide'].fit_s.mean():.1f} s.
"""
    open(f"{HERE}/stage6_wide_report.md", "w", encoding="utf8").write(md); print("report written")


if __name__ == "__main__":
    if sys.argv[1:] == ["report"]: report()
    else:
        with Pool(8) as pool: pool.map(job, [(f, p) for f in FP for p in range(1, 33) if not os.path.exists(f"{OUT}/res_f{f}_p{p}.csv")], chunksize=1)
        report()
