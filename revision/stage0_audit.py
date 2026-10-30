import os
for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[v] = "1"
import sys, json, re, time, subprocess, importlib.util
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); PC = f"{ROOT}/predictor_comparison"; OUT = f"{HERE}/stage0"; os.makedirs(OUT, exist_ok=True)
import numpy as np, pandas as pd, GPy
sys.path.insert(0, PC)
spec = importlib.util.spec_from_file_location("pc_cost", f"{PC}/practical_cost/run.py"); pc = importlib.util.module_from_spec(spec); spec.loader.exec_module(pc)
hy = pc.hy


def src(nb):
    return "".join(json.load(open(f"{ROOT}/{nb}.ipynb", encoding="utf8"))["cells"][0]["source"]).split("\n")


def lines(nb, pats):
    rows = []
    for i, l in enumerate(src(nb), 1):
        for p in pats:
            if re.search(p, l): rows.append((nb, i, l.strip()))
    return rows


print("== code lines"); L = []
for nb, pats in {"raw_al": [r"rng.choice\(idx_all\)", r"argmin\(var", r"n_train = ", r"\.optimize\("], "svts": [r"ALPHA, BETA", r"score = ", r"argmax\(score", r"rng.choice\(idx_all\)", r"n_train = ", r"centroid_c", r"feat_r_c ="],
                 "sde": [r"ALPHA, BETA", r"abs\(xs\[idx\]", r"def sde_select", r"RATIO = ", r"rng.shuffle", r"GPy"], "hybrid_al_sde": [r"AL_RATIO", r"argsort\(var\)", r"step ="],
                 "hybrid_svss_asde": [r"AL_RATIO =", r"SDE_RATIO =", r"SVTS_ALPHA|SVTS_BETA", r"K_NEIGHBORS =|ROUGHNESS_THRESH =", r"default_ab = ", r"int\(np\.ceil", r"step = ", r"target = ", r"score = ", r"\.optimize\(", r"argsort\(scores\)", r"leftover", r"CALIB_FPGA"]}.items():
    L += lines(nb, pats)
pd.DataFrame(L, columns=["notebook(cell 0)", "line", "code"]).to_csv(f"{OUT}/code_lines.csv", index=False)
for r in L: print(*r)

df = hy.load_dat(1, 1); X = df[["X", "Y"]].values.astype(float); Y = df.Value.values.astype(float); n = len(df)
print("\n== svts similarity term (path 1, additive GPR on a random 1% set)")
rng = np.random.default_rng(42); tr = list(rng.choice(n, size=32, replace=False)); rem = np.setdiff1d(np.arange(n), tr)
m = pc.fit(X[tr], Y[tr]); mr, vr = m.predict(X[rem]); ms, vs = m.predict(X[tr])
am, av = np.concatenate([mr.ravel(), ms.ravel()]), np.concatenate([vr.ravel(), vs.ravel()]); mz, vz = (am - am.mean()) / am.std(), (av - av.mean()) / av.std(); k = len(mr)
fr, fs = np.column_stack([mz[:k], vz[:k]]), np.column_stack([mz[k:], vz[k:]]); c = fs.mean(0); cc = c - c.mean(); frc = fr - fr.mean(1, keepdims=True)
sim = (frc @ cc) / (np.linalg.norm(frc, axis=1) * np.linalg.norm(cc) + 1e-9); sn = (sim - sim.min()) / (sim.max() - sim.min() + 1e-9)
print(pd.Series(dict(min_sim=sim.min(), max_sim=sim.max(), share_abs_sim_gt_0_999=(np.abs(sim) > 0.999).mean(), share_sim_n_gt_0_999=(sn > 0.999).mean(), share_sim_n_lt_0_001=(sn < 0.001).mean())).round(5).to_string())
print("rank of feat_r_c:", np.linalg.matrix_rank(frc), "(rows are multiples of (1,-1))")

print("\n== predictive variance includes the noise term"); a = m.predict(X[rem][:3])[1].ravel(); b = m.predict(X[rem][:3], include_likelihood=False)[1].ravel(); print("var with noise - var without:", (a - b).round(4), "noise var (original units):", round(float(m.Gaussian_noise.variance[0]) * float(m.normalizer.std ** 2), 4))

print("\n== budgets: orders of practical_cost (FPGA 1)")
O = pd.concat(pd.read_csv(f"{PC}/practical_cost/results/orders_path{p}.csv").assign(PATH=p) for p in range(1, 33)); B = []
for p in range(1, 33):
    d = hy.load_dat(1, p); act = set(zip(d.X, d.Y)); g = O[O.PATH == p]
    for mth, h in g.groupby("method"):
        c_ = list(zip(h.X, h.Y)); B.append(dict(PATH=p, method=mth, n_active=len(d), n=len(c_), unique=len(set(c_)), in_active=all(c in act for c in c_), in_range=bool(h.X.between(0, 32).all() and h.Y.between(0, 119).all()),
                                                  floor_10pct=int(np.floor(0.10 * len(d))), round_10pct=int(round(0.10 * len(d)))))
B = pd.DataFrame(B); B.to_csv(f"{OUT}/budget_check.csv", index=False)
print(B.groupby("method").agg(n_min=("n", "min"), n_max=("n", "max"), all_unique=("unique", lambda s: bool((s == B.loc[s.index, "n"]).all())), all_active=("in_active", "all"), all_in_range=("in_range", "all")).to_string())
print("n_active per path (min/mean/max):", B.n_active.min(), round(B.n_active.mean(), 1), B.n_active.max(), "| paths with n_active = 3173:", int((B[B.method == "sde"].n_active == 3173).sum()), "of 32")
P = pd.read_csv(f"{PC}/practical_cost/results/per_path_all.csv"); sh = {m: (317 - P[f"full_{m}_n"]).mean() for m in ["al", "svts", "sde", "hybrid"]}
print("mean number of selected sites lost on other FPGAs to the 3-sigma filter (vs 317):", {k: round(v, 2) for k, v in sh.items()})
print("same-count test: sizes by method:", B.groupby("method").n.agg(["min", "max"]).to_dict("index"))

print("\n== ASDE windows (why_sde meta, 32 paths)")
M = pd.concat(pd.read_csv(f"{PC}/why_sde/results/meta_path{p}.csv") for p in range(1, 33)); aw = np.array([eval(x) for x in M.hyb_anchor_win_mean]); dw = np.array([eval(x) for x in M.hyb_default_win])
print(pd.DataFrame({"anchor_mean_window_x": aw[:, 0], "anchor_mean_window_y": aw[:, 1], "code_default_floor_x": dw[:, 0], "code_default_floor_y": dw[:, 1], "round_x": np.round(aw[:, 0]), "round_y": np.round(aw[:, 1]), "ceil_x": np.ceil(aw[:, 0]), "ceil_y": np.ceil(aw[:, 1])}).agg(["mean", "min", "max"]).round(3).to_string())
print("share of final hybrid sites: anchors %.3f, ASDE accepted %.3f, random fill %.3f" % (M.n_anchor.mean() / 317, M.hyb_acc.mean() / 317, (M.n_asde_hyb - M.hyb_acc).mean() / 317))
print("SDE: sites passing the 5x5 test %.3f of sites, random fill %.3f" % ((M.sde_acc / M.sde_n).mean(), 1 - (M.sde_acc / M.sde_n).mean()))


def asde(df, train, target, mode):
    K, TH, SEED = hy.K_NEIGHBORS, hy.ROUGHNESS_THRESH, hy.SEED; ab = {}
    for j in train:
        xj, yj = df.loc[j, "X"], df.loc[j, "Y"]; nn = sorted([p for p in train if p != j], key=lambda p: (df.loc[p, "X"] - xj) ** 2 + (df.loc[p, "Y"] - yj) ** 2)[:K]
        dx = np.mean([abs(df.loc[p, "X"] - xj) for p in nn]); dy = np.mean([abs(df.loc[p, "Y"] - yj) for p in nn]); bump = 1 if np.std([df.loc[p, "Value"] for p in nn]) > TH else 0
        ab[j] = (int(np.ceil(dx)) + bump, int(np.ceil(dy)) + bump)
    mx, my = np.mean([v[0] for v in ab.values()]), np.mean([v[1] for v in ab.values()])
    d = {"floor": (int(mx), int(my)), "float": (mx, my), "round": (round(mx), round(my)), "ceil": (int(np.ceil(mx)), int(np.ceil(my)))}[mode]
    sel, add = list(train), []
    for i in np.random.RandomState(SEED).permutation(df.index):
        if i in train or len(add) >= target: continue
        x, y = df.loc[i, "X"], df.loc[i, "Y"]
        if all(abs(df.loc[j, "X"] - x) > ab.get(j, d)[0] or abs(df.loc[j, "Y"] - y) > ab.get(j, d)[1] for j in sel): add.append(i); sel.append(i); ab[i] = d
    return add, d


print("\n== ASDE default window: floor (code) vs real-valued mean (paper text) vs round vs ceil, paths 1-4")
rows = []
for p in range(1, 5):
    d_ = hy.load_dat(1, p); anchors = pc.grow(d_, round(len(d_) * 0.03), 0.5, 0.5, True); tgt = round(len(d_) * 0.07); res = {}
    for md in ("floor", "float", "round", "ceil"): res[md] = asde(d_, anchors, tgt, md)
    for md in ("float", "round", "ceil"):
        rows.append(dict(PATH=p, mode=md, window=res[md][1], window_floor=res["floor"][1], accepted_additions=len(res[md][0]), same_accepted_as_floor=res[md][0] == res["floor"][0]))
    rows.append(dict(PATH=p, mode="floor", window=res["floor"][1], window_floor=res["floor"][1], accepted_additions=len(res["floor"][0]), same_accepted_as_floor=True))
A = pd.DataFrame(rows); A.to_csv(f"{OUT}/asde_window_modes.csv", index=False); print(A.to_string(index=False))

print("\n== leakage test: perturb the values of every NON-selected site of FPGA 1 (path 1) and re-run the selection; the order must not change")
rng = np.random.default_rng(7); res = []
for mth in ("hybrid", "svts", "al", "sde"):
    idx = pc.select(mth, df); mask = np.zeros(n, bool); mask[idx] = True; d2 = df.copy(); d2.loc[~mask, "Value"] = d2.loc[~mask, "Value"].values + rng.normal(0, 5 * df.Value.std(), (~mask).sum())
    idx2 = pc.select(mth, d2); res.append(dict(method=mth, n_selected=len(idx), order_identical=list(idx) == list(idx2)))
    print(res[-1], flush=True)
pd.DataFrame(res).to_csv(f"{OUT}/leakage_test.csv", index=False)

print("\n== GPR settings: selection loops use the default optimiser (max 1000 iterations), the predictor max_iters=200; same data, compare")
S = pc.select("hybrid", df); Xs, ys = X[S], Y[S]; pcs = []
for mi in (None, 200, 200):
    m = pc.fit(Xs, ys, "add", mi); pcs.append(dict(max_iters=mi, LL=float(m.log_likelihood()), ls=[round(float(p.lengthscale[0]), 4) for p in m.kern.parts], var=[round(float(p.variance[0]), 4) for p in m.kern.parts], noise=round(float(m.Gaussian_noise.variance[0]), 5), status=m.optimization_runs[-1].status, iters=int(m.optimization_runs[-1].funct_eval)))
print(pd.DataFrame(pcs).to_string(index=False))
m = GPy.models.GPRegression(Xs, ys.reshape(-1, 1), pc.kern(), normalizer=True)
print("initial values: lengthscales", [float(p.lengthscale[0]) for p in m.kern.parts], "variances", [float(p.variance[0]) for p in m.kern.parts], "noise", float(m.Gaussian_noise.variance[0]), "| constraints:", m.kern.parts[0].lengthscale.constraints.properties()[0] if m.kern.parts[0].lengthscale.constraints.properties() else "none")

print("\n== timing environment")
from threadpoolctl import threadpool_info
print("BLAS/OpenMP threads in this process:", [(i["internal_api"], i["num_threads"]) for i in threadpool_info()])
ps = lambda q: subprocess.run(["powershell", "-NoProfile", "-Command", q], capture_output=True, text=True).stdout.strip()
print("power plan:", ps("powercfg /getactivescheme"), "| current clock MHz:", ps("(Get-CimInstance Win32_Processor | Select -First 1).CurrentClockSpeed"), "| max MHz:", ps("(Get-CimInstance Win32_Processor | Select -First 1).MaxClockSpeed"))
print("processes using >20 CPU s:", ps("Get-Process | Where-Object {$_.CPU -gt 20} | Sort CPU -desc | Select -First 8 | ForEach-Object { $_.Name + ' ' + [int]$_.CPU }").replace("\n", "; "))
T = []
for r in range(3):
    t0 = time.perf_counter(); pc.select("hybrid", df); T.append(time.perf_counter() - t0)
print("hybrid selection path 1, 3 repeats (s):", np.round(T, 2), "-> spread %.0f%%" % (100 * (max(T) - min(T)) / np.median(T)))
