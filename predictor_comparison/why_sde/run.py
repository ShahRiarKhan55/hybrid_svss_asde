"""Why do hybrid SVSS+ASDE and SDE tie at 10 %?  Code audit, site geometry, ablations, coverage-vs-error.
Selection is on FPGA 1 per path (as in all earlier runs); hybrid and SDE orders are read from budget/results/sel_path{p}.csv,
ASDE-only (random 1 % start) and SVSS-only (SVSS continued to 10 %) are generated here with the existing functions.
Predictor: additive RBF GPR (GPy, raw coordinates). Test set per (FPGA, path, budget): sites outside ALL training sets of that budget
(hybrid, SDE, ASDE-only, SVSS-only and the 10 random sets), identical for every strategy.  Usage: python run.py (8 workers, resumable)"""
import os
for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[v] = "1"
import sys
HERE = os.path.dirname(os.path.abspath(__file__)); PC = os.path.dirname(HERE); ROOT = os.path.dirname(PC)
sys.path.insert(0, PC)  # `import run` is the earlier predictor_comparison/run.py (same file name as this script)
import time
from multiprocessing import Pool
import GPy, numpy as np, pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr, pearsonr, wilcoxon
import run

R = f"{HERE}/results"; BR = f"{PC}/budget/results"
hy = run.hy  # hybrid_svss_asde.ipynb cell 0: load_dat, svts_select (SVSS), adaptive_sde (ASDE)
BUDGETS = [3, 5, 7, 10]; STRATS = ["hybrid", "sde", "asde_only", "svss_only"]; SEEDS = range(1, 11)


def asde_split(df, train, target):
    """Copy of hy.adaptive_sde that also reports how many additions passed the exclusion rule (the rest are random leftovers) and the windows."""
    K, TH, SEED = hy.K_NEIGHBORS, hy.ROUGHNESS_THRESH, hy.SEED
    ab = {}
    for j in train:
        xj, yj = df.loc[j, "X"], df.loc[j, "Y"]
        nn = sorted([p for p in train if p != j], key=lambda p: (df.loc[p, "X"] - xj) ** 2 + (df.loc[p, "Y"] - yj) ** 2)[:K]
        if not nn:
            ab[j] = (1, 1); continue
        dx = np.mean([abs(df.loc[p, "X"] - xj) for p in nn]); dy = np.mean([abs(df.loc[p, "Y"] - yj) for p in nn])
        bump = 1 if np.std([df.loc[p, "Value"] for p in nn]) > TH else 0
        ab[j] = (int(np.ceil(dx)) + bump, int(np.ceil(dy)) + bump)
    default_ab = (int(np.mean([v[0] for v in ab.values()])), int(np.mean([v[1] for v in ab.values()])))
    anchor_ab = list(ab.values()); selected, sde = list(train), []
    for i in np.random.RandomState(SEED).permutation(df.index):
        if i in train or len(sde) >= target: continue
        x, y = df.loc[i, "X"], df.loc[i, "Y"]
        if all(abs(df.loc[j, "X"] - x) > ab.get(j, default_ab)[0] or abs(df.loc[j, "Y"] - y) > ab.get(j, default_ab)[1] for j in selected):
            sde.append(i); selected.append(i); ab[i] = default_ab
    n_acc = len(sde)
    if len(sde) < target:
        sde += [i for i in np.random.RandomState(SEED).permutation(df.index) if i not in train and i not in sde][:target - len(sde)]
    return sde, n_acc, default_ab, anchor_ab


def geom(df, idx, tree_all=None):
    """Geometry of a training set (indices into df) on the active-site grid, distances in cells (Euclidean)."""
    T = df.loc[idx, ["X", "Y"]].values.astype(float); A = df[["X", "Y"]].values.astype(float); t = cKDTree(T)
    nn = t.query(T, k=2)[0][:, 1]; d = t.query(A)[0]  # nearest other training site; nearest training site of every active site
    return dict(nn_mean=nn.mean(), nn_p5=np.quantile(nn, .05), max_gap=d.max(), mean_fill=d.mean(), share2=(d <= 2).mean())


def select(p):
    f = f"{R}/sel_path{p}.csv"
    if os.path.exists(f):
        return
    df = hy.load_dat(1, p); n = len(df); S = pd.read_csv(f"{BR}/sel_path{p}.csv"); nb = {b: int(S[f"n_{b}"].iloc[0]) for b in BUDGETS}
    cid = {(x, y): i for i, (x, y) in enumerate(zip(df.X, df.Y))}; ids = lambda g: [cid[c] for c in zip(g.X, g.Y)]
    hyb = S[S.strategy == "hybrid"]; anchors, asde_h = ids(hyb[hyb.phase == "svss"]), ids(hyb[hyb.phase == "asde"])
    # hybrid audit: windows, accepted vs leftover, and a check that the copied function reproduces the saved ASDE order
    h_adds, h_acc, h_ab, h_anchor_ab = asde_split(df, anchors, len(asde_h))
    # SVSS-only: the existing svts_select continued to the full budget; ASDE-only: ASDE from the same random 1 % start
    hy.AL_RATIO = nb[10] / n; svss_only = [int(i) for i in hy.svts_select(df)]
    step = max(1, round(n * 0.01)); a_adds, a_acc, a_ab, a_anchor_ab = asde_split(df, svss_only[:step], nb[10] - step); asde_only = svss_only[:step] + a_adds
    order = {"hybrid": anchors + asde_h, "sde": ids(S[S.strategy == "sde"]), "asde_only": asde_only, "svss_only": svss_only}
    phase = {"hybrid": ["svss"] * len(anchors) + ["asde_accepted"] * h_acc + ["asde_leftover"] * (len(asde_h) - h_acc), "sde": ["sde"] * len(order["sde"]),
             "asde_only": ["random_start"] * step + ["asde_accepted"] * a_acc + ["asde_leftover"] * (len(a_adds) - a_acc), "svss_only": ["svss"] * len(svss_only)}
    o = pd.concat([df.loc[order[s], ["X", "Y"]].assign(strategy=s, phase=phase[s], rank=range(1, len(order[s]) + 1)) for s in STRATS])
    rows = []
    for b in BUDGETS:
        sets = {s: order[s][:nb[b]] for s in STRATS}
        for s in STRATS: rows.append(dict(strategy=s, budget=b, **geom(df, sets[s])))
        g = [geom(df, list(np.random.default_rng([i, 1, p]).permutation(n)[:nb[b]])) for i in SEEDS]
        rows.append(dict(strategy="random", budget=b, **{k: np.mean([x[k] for x in g]) for k in g[0]}))
    G = pd.DataFrame(rows); G["PATH"] = p
    ov = {}  # overlap at 10 %: share of the row set that is also in the column set
    sets = {s: set(order[s][:nb[10]]) for s in STRATS}; rnd = [set(np.random.default_rng([i, 1, p]).permutation(n)[:nb[10]]) for i in SEEDS]
    for a in STRATS + ["random"]:
        for c in STRATS + ["random"]:
            A = [sets[a]] if a != "random" else rnd; C = [sets[c]] if c != "random" else rnd
            ov[f"{a}>{c}"] = np.mean([len(x & y) / len(x) for x in A for y in C if x is not y])
    xs, ys = df.X.values[order["sde"]], df.Y.values[order["sde"]]  # SDE: sites that pass the |dx|>2 or |dy|>2 test against all earlier sites (the rest are leftovers)
    sde_acc = sum(k == 0 or bool(np.all((np.abs(xs[:k] - xs[k]) > 2) | (np.abs(ys[:k] - ys[k]) > 2))) for k in range(len(xs)))
    meta = dict(PATH=p, sde_n=len(xs), sde_acc=sde_acc, n_active=n, step=step, n_anchor=len(anchors), n_asde_hyb=len(asde_h), hyb_acc=h_acc, hyb_default_win=h_ab, hyb_anchor_win_mean=np.mean(h_anchor_ab, axis=0).tolist(),
                hyb_reproduced=bool(h_adds == asde_h), asde_only_acc=a_acc, asde_only_n_add=len(a_adds), asde_only_default_win=a_ab, asde_only_anchor_win_mean=np.mean(a_anchor_ab, axis=0).tolist(),
                svss_prefix_same_as_hybrid_anchors=len(set(svss_only[:len(anchors)]) & set(anchors)) / len(anchors), **{f"ov_{k}": v for k, v in ov.items()})
    G.to_csv(f"{R}/geom_path{p}.csv", index=False); pd.DataFrame([meta]).to_csv(f"{R}/meta_path{p}.csv", index=False)
    o.to_csv(f + ".tmp", index=False); os.replace(f + ".tmp", f)
    print(time.strftime("%H:%M:%S"), "selected path", p, flush=True)


def gp_err(tr, te, ds):  # additive RBF GPR, paper recipe, raw coordinates
    m = GPy.models.GPRegression(tr[["X", "Y"]].values, tr[["Value"]].values, GPy.kern.RBF(1, active_dims=[0]) + GPy.kern.RBF(1, active_dims=[1]), normalizer=True)
    m.optimize(max_iters=200, messages=False)
    return np.abs(m.predict(te[["X", "Y"]].values)[0].ravel() - te.Value.values).mean() / ds


def evaluate(a):
    f, p = a
    out = f"{R}/path_f{f}_p{p}.csv"
    if os.path.exists(out):
        return
    df = run.load_dat(f, p); ds = df.Value.max() - df.Value.min(); coords = list(zip(df.X, df.Y))
    S = pd.read_csv(f"{R}/sel_path{p}.csv")  # all four strategies (hybrid and SDE orders copied from budget/results)
    nbs = pd.read_csv(f"{BR}/sel_path{p}.csv").iloc[0]; perm = {i: np.random.default_rng([i, f, p]).permutation(len(df)) for i in SEEDS}
    v = df.set_index(["X", "Y"]).Value  # neighbour-difference noise-floor estimate: copy-the-neighbour error x 0.707 (see analysis)
    dx = (v.reindex(pd.MultiIndex.from_arrays([df.X + 1, df.Y])).values - df.Value.values); dy = (v.reindex(pd.MultiIndex.from_arrays([df.X, df.Y + 1])).values - df.Value.values)
    row = dict(FPGA=f, PATH=p, adj_mae=np.nanmean(np.abs(np.concatenate([dx, dy]))) / ds)
    for b in BUDGETS:
        nb = int(nbs[f"n_{b}"]); mask = {}
        for s in STRATS:
            g = S[S.strategy == s].sort_values("rank").iloc[:nb]; gs = set(zip(g.X, g.Y)); mask[s] = np.array([c in gs for c in coords])
        for i in SEEDS:
            m = np.zeros(len(df), bool); m[perm[i][:mask["hybrid"].sum()]] = True; mask[f"rand{i}"] = m
        test = ~np.any(list(mask.values()), axis=0); te = df[test]; row[f"{b}_n_test"] = len(te)
        for s, mk in mask.items():
            np.random.seed(f * 100 + p); row[f"{b}_{s}_n"] = int(mk.sum()); row[f"{b}_{s}_delta"] = gp_err(df[mk], te, ds)
        row[f"{b}_const_delta"] = np.abs(te.Value.values - df[mask["hybrid"]].Value.mean()).mean() / ds
    pd.DataFrame([row]).to_csv(out + ".tmp", index=False); os.replace(out + ".tmp", out)
    print(time.strftime("%H:%M:%S"), a, "done", flush=True)


def boot(d, rng): return [d[rng.integers(0, len(d), len(d))].mean() for _ in range(10000)]


def analyse():
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
    M = pd.concat(pd.read_csv(f"{R}/meta_path{p}.csv") for p in range(1, 33)); G = pd.concat(pd.read_csv(f"{R}/geom_path{p}.csv") for p in range(1, 33))
    P = pd.concat(pd.read_csv(f"{R}/path_f{f}_p{p}.csv") for f in range(1, 51) for p in range(1, 33)); P.to_csv(f"{R}/per_path_all.csv", index=False)
    print("== 1. how many final hybrid sites come from SVSS / ASDE (mean over 32 paths, 10 % budget)")
    h = M[["n_anchor", "n_asde_hyb", "hyb_acc"]].mean(); tot = h.n_anchor + h.n_asde_hyb
    print(f"SVSS anchors {h.n_anchor:.1f} ({h.n_anchor / tot:.1%}); ASDE additions {h.n_asde_hyb:.1f} ({h.n_asde_hyb / tot:.1%}), of which accepted by the exclusion rule {h.hyb_acc:.1f} "
          f"({h.hyb_acc / h.n_asde_hyb:.1%} of additions), random leftovers {h.n_asde_hyb - h.hyb_acc:.1f}; copied ASDE reproduces the saved order on {int(M.hyb_reproduced.sum())}/32 paths")
    print(f"SDE: {M.sde_n.mean():.1f} sites, {M.sde_acc.mean():.1f} pass the fixed 5x5 exclusion test ({(M.sde_acc / M.sde_n).mean():.1%}), the rest are random leftovers")
    print(f"ASDE-only (random 1 % start): additions {M.asde_only_n_add.mean():.1f}, accepted by rule {M.asde_only_acc.mean():.1f} ({(M.asde_only_acc / M.asde_only_n_add).mean():.1%}), leftovers random")
    win = lambda c: np.array([eval(x) for x in M[c]]).mean(axis=0)
    print("ASDE default window (dx,dy; a candidate is excluded when |dx|<=a and |dy|<=b): hybrid", win("hyb_default_win").round(2), "| anchor mean", win("hyb_anchor_win_mean").round(2),
          "| ASDE-only", win("asde_only_default_win").round(2), "| anchor mean", win("asde_only_anchor_win_mean").round(2), "| SDE fixed (2,2)")
    print("SVSS-only first 3 % equals hybrid anchors (share):", M.svss_prefix_same_as_hybrid_anchors.mean().round(3))
    print("\n== 2. geometry on FPGA 1 (mean over 32 paths; distances in cells), budget 10 %")
    g10 = G[G.budget == 10].groupby("strategy")[["nn_mean", "nn_p5", "max_gap", "mean_fill", "share2"]].mean(); g10.to_csv(f"{R}/geometry_10pct.csv"); print(g10.round(3).to_string())
    ov = M[[c for c in M if c.startswith("ov_")]].mean(); names = STRATS + ["random"]
    O = pd.DataFrame([[ov[f"ov_{a}>{c}"] for c in names] for a in names], index=names, columns=names); O.to_csv(f"{R}/overlap_10pct.csv"); print("overlap (share of ROW set also in COLUMN set):\n", O.round(3).to_string())
    print("\n== 3. ablations")
    SS = STRATS + ["random"]; rng, rows, tests = np.random.default_rng(0), [], []
    for b in BUDGETS:
        for s in ("hybrid", "sde", "asde_only", "svss_only"):
            rows.append(dict(budget=f"{b}%", strategy=s, delta_bar=P.groupby("FPGA")[f"{b}_{s}_delta"].mean().mean(), n_train=P[f"{b}_{s}_n"].mean(), n_test=P[f"{b}_n_test"].mean()))
        sd = [P.groupby("FPGA")[f"{b}_rand{i}_delta"].mean().mean() for i in SEEDS]
        rows.append(dict(budget=f"{b}%", strategy="random (10 seeds)", delta_bar=np.mean(sd), seed_sd=np.std(sd, ddof=1), seed_min=min(sd), seed_max=max(sd), n_train=P[f"{b}_rand1_n"].mean(), n_test=P[f"{b}_n_test"].mean()))
        rows.append(dict(budget=f"{b}%", strategy="constant (train mean)", delta_bar=P.groupby("FPGA")[f"{b}_const_delta"].mean().mean()))
        P[f"{b}_random_delta"] = P[[f"{b}_rand{i}_delta" for i in SEEDS]].mean(axis=1)
        V = {s: P.groupby("PATH")[f"{b}_{s}_delta"].mean() for s in ("hybrid", "sde", "asde_only", "svss_only", "random")}
        for o in ("sde", "asde_only", "svss_only", "random"):
            d = (V["hybrid"] - V[o]).values; bs = boot(d, rng)
            tests.append(dict(budget=f"{b}%", comparison=f"hybrid - {o}", mean_diff=d.mean(), ci95_low=np.quantile(bs, .025), ci95_high=np.quantile(bs, .975), paths_hybrid_better=int((d < 0).sum()), of=32, wilcoxon_p=wilcoxon(d).pvalue))
    A, T = pd.DataFrame(rows), pd.DataFrame(tests); T["bonf_16"] = T.wilcoxon_p < 0.05 / len(T); T["bonf_4_budgets"] = T.wilcoxon_p < 0.05 / 4
    A.to_csv(f"{R}/ablation_table.csv", index=False); T.to_csv(f"{R}/ablation_tests.csv", index=False); print(A.round(4).to_string(index=False)); print(); print(T.round(5).to_string(index=False))
    print("\n== 4. coverage vs error (rows = strategy x budget x path, path-level delta averaged over 50 FPGAs)")
    D = []
    for b in BUDGETS:
        for s in SS:
            col = f"{b}_{s}_delta"; dd = P.groupby("PATH")[col].mean()
            g = G[(G.budget == b) & (G.strategy == s)].set_index("PATH"); D.append(g.assign(delta=dd, budget=b, strategy=s)[["delta", "budget", "strategy", "nn_mean", "nn_p5", "max_gap", "mean_fill", "share2"]].reset_index())
    D = pd.concat(D); D.to_csv(f"{R}/coverage_vs_delta.csv", index=False)
    X = D[D.strategy.isin(SS)].copy(); X["cell"] = X.PATH.astype(str) + "_" + X.budget.astype(str)
    cor = []
    for c in ("nn_mean", "nn_p5", "max_gap", "mean_fill", "share2"):
        within = X.assign(dc=X[c] - X.groupby("cell")[c].transform("mean"), dd=X.delta - X.groupby("cell").delta.transform("mean"))
        bud = X.groupby("budget")[[c, "delta"]].mean(); strat = X.groupby("strategy")[[c, "delta"]].mean()
        cor.append(dict(measure=c, spearman_pooled=spearmanr(X[c], X.delta)[0], spearman_within_path_budget=spearmanr(within.dc, within.dd)[0], pearson_within=pearsonr(within.dc, within.dd)[0],
                        corr_across_budgets=pearsonr(bud[c], bud.delta)[0], corr_across_strategies_at_10pct=pearsonr(X[X.budget == 10].groupby("strategy")[c].mean(), X[X.budget == 10].groupby("strategy").delta.mean())[0]))
    C = pd.DataFrame(cor); C.to_csv(f"{R}/coverage_correlations.csv", index=False); print(C.round(3).to_string(index=False))
    print("\n== 5. noise-floor evidence")
    fl = 0.7071 * P.groupby("FPGA").adj_mae.mean().mean()  # E|n| / E|n1-n2| = 1/sqrt(2) for Gaussian noise
    print(f"copy-the-neighbour error {P.groupby('FPGA').adj_mae.mean().mean():.4f} -> upper bound on the noise floor (Gaussian noise, perfect signal) ~ {fl:.4f}")
    print("delta vs budget (hybrid):", {f"{b}%": round(P.groupby('FPGA')[f'{b}_hybrid_delta'].mean().mean(), 4) for b in BUDGETS},
          "| constant predictor 10%:", round(P.groupby("FPGA")["10_const_delta"].mean().mean(), 4))
    print("non-finite deltas:", int((~np.isfinite(P.filter(like="_delta").values)).sum()))


if __name__ == "__main__":
    os.makedirs(R, exist_ok=True)
    with Pool(8) as pool:
        pool.map(select, [p for p in range(1, 33) if not os.path.exists(f"{R}/sel_path{p}.csv")], chunksize=1)
        pool.map(evaluate, [(f, p) for f in range(1, 51) for p in range(1, 33) if not os.path.exists(f"{R}/path_f{f}_p{p}.csv")], chunksize=1)
    analyse()
