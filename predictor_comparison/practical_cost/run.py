"""Practical cost of the site-selection strategies, with the additive RBF GPR inside every selection loop.
  python run.py time   -> single thread, sequential, nothing else running: selects on FPGA 1 for all 32 paths (AL, SVTS, SDE, hybrid),
                          saves the selection ORDERS, times each method, counts GPR fits; then times the final GPR fit + predict on 3 other FPGAs x 32 paths.
  python run.py eval   -> 8 workers, resumable: delta-bar vs number of training points on the 49 other FPGAs (common test set), thresholds, hardware info.
Selection rules are the existing ones (raw_al.ipynb: add the minimum-variance site; svts.ipynb: alpha=0.7, beta=0.3 one site per step;
hybrid: SVSS 3 % batches alpha=beta=0.5 then the existing ASDE; sde.ipynb: unchanged, it has no GPR), only the GPR kernel is the additive RBF."""
import os
for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[v] = "1"
import sys
HERE = os.path.dirname(os.path.abspath(__file__)); PC = os.path.dirname(HERE); ROOT = os.path.dirname(PC)
sys.path.insert(0, PC)  # `import run` is the earlier predictor_comparison/run.py (same file name as this script)
import json, platform, random, subprocess, time, types
from multiprocessing import Pool
import GPy, numpy as np, pandas as pd
import run

R = f"{HERE}/results"; os.makedirs(R, exist_ok=True)
hy = run.hy  # hybrid_svss_asde.ipynb cell 0: load_dat, adaptive_sde (ASDE, reused unchanged)
sde = types.ModuleType("sde"); exec("".join(json.load(open(f"{ROOT}/sde.ipynb", encoding="utf8"))["cells"][0]["source"]), sde.__dict__)  # sde_select, SEED
METHODS = ["al", "svts", "sde", "hybrid"]; N_FITS = [0]


def kern(kind="add"):
    return GPy.kern.RBF(1, active_dims=[0]) + GPy.kern.RBF(1, active_dims=[1]) if kind == "add" else GPy.kern.RBF(2)


def fit(X, y, kind="add", max_iters=None):  # one GPR fit as used inside the selection loops (default optimiser, like the notebooks)
    N_FITS[0] += 1
    m = GPy.models.GPRegression(X, y.reshape(-1, 1), kern(kind), normalizer=True)
    m.optimize(messages=False) if max_iters is None else m.optimize(max_iters=max_iters, messages=False)
    return m


def al(df, kind="add", mode="min"):  # raw_al.ipynb: add the remaining site with the LOWEST predictive variance, one per step (mode="max": highest, extra)
    X, Y = df[["X", "Y"]].values.astype(float), df.Value.values.astype(float); n = len(Y); idx = np.arange(n)
    tr = [np.random.default_rng(hy.SEED).choice(idx)]
    while len(tr) < int(n * 0.10):
        rem = np.setdiff1d(idx, tr); _, var = fit(X[tr], Y[tr], kind).predict(X[rem]); tr.append(rem[np.argmin(var.ravel()) if mode == "min" else np.argmax(var.ravel())])
    return [int(i) for i in tr]


def grow(df, target, alpha, beta, batch, kind="add"):  # svts.ipynb (one site per step) / SVSS of the hybrid (batches of 1 %): similarity + variance score
    X, Y = df[["X", "Y"]].values.astype(float), df.Value.values.astype(float); n = len(Y); rng = np.random.default_rng(hy.SEED)
    step = max(1, round(n * 0.01)); tr = [int(i) for i in rng.choice(n, size=step, replace=False)] if batch else [int(rng.choice(np.arange(n)))]
    while len(tr) < target:
        rem = np.setdiff1d(np.arange(n), tr); m = fit(X[tr], Y[tr], kind)
        mr, vr = m.predict(X[rem]); ms, vs = m.predict(X[tr])
        am, av = np.concatenate([mr.ravel(), ms.ravel()]), np.concatenate([vr.ravel(), vs.ravel()])
        mz, vz = (am - am.mean()) / (am.std() + 1e-9), (av - av.mean()) / (av.std() + 1e-9); nr = len(mr)
        fr, fs = np.column_stack([mz[:nr], vz[:nr]]), np.column_stack([mz[nr:], vz[nr:]])
        c = fs.mean(axis=0); cc = c - c.mean(); frc = fr - fr.mean(axis=1, keepdims=True)
        sim = (frc @ cc) / (np.linalg.norm(frc, axis=1) * np.linalg.norm(cc) + 1e-9)
        sn = (sim - sim.min()) / (sim.max() - sim.min() + 1e-9); vn = (vr.ravel() - vr.min()) / (vr.max() - vr.min() + 1e-9)
        score = alpha * sn + beta * vn; pick = min(step if batch else 1, target - len(tr), len(rem))
        tr += [int(rem[np.argmax(score)])] if not batch else [int(rem[i]) for i in np.argsort(score)[::-1][:pick]]  # argmax tie-break (first index) as in svts.ipynb
    return tr


def select(method, df, kind="add"):
    n = len(df)
    if method == "al": return al(df, kind)
    if method == "almax": return al(df, kind, "max")
    if method == "almax": return al(df, kind, "max")
    if method == "svts": return grow(df, int(n * 0.10), 0.7, 0.3, False, kind)
    if method == "sde": return sde.sde_select(df.X.values, df.Y.values, int(round(0.10 * n)), random.Random(sde.SEED))
    anchors = grow(df, max(1, round(n * 0.03)), 0.5, 0.5, True, kind)  # hybrid = SVSS 3 % (GPR) + existing ASDE to 10 % (no GPR)
    hy.SDE_RATIO = 0.07; return anchors + [int(i) for i in hy.adaptive_sde(df, anchors)]


def hardware():
    ps = lambda q: subprocess.run(["powershell", "-NoProfile", "-Command", q], capture_output=True, text=True).stdout.strip()
    return dict(cpu=ps("(Get-CimInstance Win32_Processor | Select -First 1).Name"), cores_threads=ps("$p=Get-CimInstance Win32_Processor; \"$($p.NumberOfCores) cores / $($p.NumberOfLogicalProcessors) threads\""),
                ram_gb=ps("[math]::Round((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory/1GB,1)"), os=platform.platform(), python=sys.version.split()[0],
                numpy=np.__version__, gpy=GPy.__version__)


def time_stage():
    print("hardware:", hardware(), flush=True)
    part = f"{R}/selection_timing_partial.csv"; rows = pd.read_csv(part).to_dict("records") if os.path.exists(part) else []
    for p in range(1, 33):
        f = f"{R}/orders_path{p}.csv"
        if os.path.exists(f): continue
        df = hy.load_dat(1, p); o = []
        for m in METHODS:
            N_FITS[0] = 0; t0, c0 = time.perf_counter(), time.process_time(); idx = select(m, df)
            rows.append(dict(PATH=p, method=m, wall_s=time.perf_counter() - t0, cpu_s=time.process_time() - c0, gpr_fits=N_FITS[0], n_sites=len(idx)))
            o.append(df.loc[idx, ["X", "Y"]].assign(method=m, rank=range(1, len(idx) + 1)))
            pd.DataFrame(rows).to_csv(f"{R}/selection_timing_partial.csv", index=False)
        pd.concat(o).to_csv(f, index=False); print(time.strftime("%H:%M:%S"), "path", p, {r["method"]: round(r["wall_s"], 1) for r in rows[-4:]}, flush=True)
    pd.DataFrame(rows).to_csv(f"{R}/selection_timing.csv", index=False)
    # final GPR (additive RBF, max_iters=200, as in the evaluations): fit + predict on 3 other FPGAs x 32 paths, training set = the hybrid 10 % sites
    rows = []
    for fp in (2, 26, 50):
        for p in range(1, 33):
            df = run.load_dat(fp, p); O = pd.read_csv(f"{R}/orders_path{p}.csv"); h = O[O.method == "hybrid"]; pos = set(zip(h.X, h.Y))
            mk = np.array([c in pos for c in zip(df.X, df.Y)]); tr, te = df[mk], df[~mk]
            t0 = time.perf_counter(); m = fit(tr[["X", "Y"]].values, tr.Value.values, "add", 200); t1 = time.perf_counter(); m.predict(te[["X", "Y"]].values); t2 = time.perf_counter()
            rows.append(dict(FPGA=fp, PATH=p, n_train=len(tr), n_test=len(te), fit_s=t1 - t0, predict_s=t2 - t1))
    pd.DataFrame(rows).to_csv(f"{R}/final_gpr_timing.csv", index=False); print("final GPR timing done", flush=True)


def time_max():  # extra: AL with the HIGHEST predictive variance, same single-thread sequential protocol
    part = f"{R}/selection_timing_almax.csv"; rows = pd.read_csv(part).to_dict("records") if os.path.exists(part) else []
    for p in range(1, 33):
        f = f"{R}/orders_max_path{p}.csv"
        if os.path.exists(f): continue
        df = hy.load_dat(1, p); N_FITS[0] = 0; t0, c0 = time.perf_counter(), time.process_time(); idx = select("almax", df)
        rows.append(dict(PATH=p, method="almax", wall_s=time.perf_counter() - t0, cpu_s=time.process_time() - c0, gpr_fits=N_FITS[0], n_sites=len(idx)))
        pd.DataFrame(rows).to_csv(part, index=False); df.loc[idx, ["X", "Y"]].assign(method="almax", rank=range(1, len(idx) + 1)).to_csv(f, index=False)
        print(time.strftime("%H:%M:%S"), "almax path", p, round(rows[-1]["wall_s"], 1), flush=True)


GRID = list(range(16, 305, 16)) + ["full"]; SETS = METHODS + ["almax", "rand1", "rand2", "rand3"]; SHOW = [32, 64, 96, 128, 192, 256, "full"]


def orders(p):  # {method: (X, Y) list in selection order}; SDE = the SAVED order (budget/results), AL / SVTS / hybrid = this run, almax = extra run
    O = pd.read_csv(f"{R}/orders_path{p}.csv"); o = {m: list(zip(g.sort_values("rank").X, g.sort_values("rank").Y)) for m, g in O.groupby("method")}
    b = pd.read_csv(f"{PC}/budget/results/sel_path{p}.csv"); b = b[b.strategy == "sde"].sort_values("rank"); o["sde"] = list(zip(b.X, b.Y))
    a = pd.read_csv(f"{R}/orders_max_path{p}.csv"); o["almax"] = list(zip(a.X, a.Y)); return o


def evaluate(a):
    f, p = a
    out = f"{R}/path_f{f}_p{p}.csv"
    if os.path.exists(out): return
    df = run.load_dat(f, p); ds = df.Value.max() - df.Value.min(); cid = {c: i for i, c in enumerate(zip(df.X, df.Y))}
    pos = {m: np.array([cid.get(c, -1) for c in v]) for m, v in orders(p).items()}
    nfull = {m: len(v) for m, v in pos.items()}; perm = {f"rand{i}": np.random.default_rng([i, f, p]).permutation(len(df)) for i in (1, 2, 3)}
    hfull = (pos["hybrid"] >= 0).sum(); union = np.zeros(len(df), bool)
    for m in pos: union[pos[m][pos[m] >= 0]] = True
    for k, v in perm.items(): union[v[:hfull]] = True
    te = df[~union]; Xt, yt = te[["X", "Y"]].values, te.Value.values; row = dict(FPGA=f, PATH=p, n_test=len(te))
    for g in GRID:
        for s in SETS:
            if s in pos: idx = pos[s][:(nfull[s] if g == "full" else g)]; idx = idx[idx >= 0]
            else: idx = perm[s][:(hfull if g == "full" else g)]
            np.random.seed(f * 100 + p); m = fit(df[["X", "Y"]].values[idx], df.Value.values[idx], "add", 200)
            row[f"{g}_{s}_delta"] = np.abs(m.predict(Xt)[0].ravel() - yt).mean() / ds; row[f"{g}_{s}_n"] = len(idx)
    pd.DataFrame([row]).to_csv(out + ".tmp", index=False); os.replace(out + ".tmp", out); print(time.strftime("%H:%M:%S"), a, "done", flush=True)


def analyse():
    from scipy.stats import wilcoxon
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40)
    P = pd.concat(pd.read_csv(f"{R}/path_f{f}_p{p}.csv") for f in range(1, 51) for p in range(1, 33)); P.to_csv(f"{R}/per_path_all.csv", index=False)
    ALL = ["al", "almax", "svts", "sde", "hybrid", "random"]; ST = [str(x) for x in SHOW]
    col = lambda g, s: [f"{g}_rand{i}_delta" for i in (1, 2, 3)] if s == "random" else [f"{g}_{s}_delta"]
    val = lambda g, s, Q=P: Q[col(g, s)].mean(axis=1)
    rows = []
    for tag, Q in (("all50", P), ("no_fpga1", P[P.FPGA >= 2])):
        for g in GRID:
            for s in ALL:
                fp = Q.assign(v=val(g, s, Q)).groupby("FPGA").v.mean(); n = Q[f"{g}_{'rand1' if s == 'random' else s}_n"].mean()
                rows.append(dict(fpgas=tag, grid=str(g), n=round(n, 1), method=s, delta_bar=fp.mean(), device_min=fp.min(), device_max=fp.max(), device_sd=fp.std()))
    C = pd.DataFrame(rows); C.to_csv(f"{R}/delta_vs_n.csv", index=False)
    for tag in ("all50", "no_fpga1"):
        c = C[C.fpgas == tag]
        print(f"\n== delta-bar vs n ({tag}); full = end of each method's own order (315-317 sites)")
        print(c.pivot(index="grid", columns="method", values="delta_bar")[ALL].reindex(ST).round(4).to_string())
        print("cross-device range of the per-FPGA mean delta-bar (min - max over the 50 FPGAs):")
        r = c.assign(r=c.device_min.round(4).astype(str) + " - " + c.device_max.round(4).astype(str)).pivot(index="grid", columns="method", values="r")[ALL].reindex(ST); print(r.to_string())
        th = {}
        for s in ALL:
            for t in (0.10, 0.09):
                hit = [str(g) for g in GRID if c[(c.grid == str(g)) & (c.method == s)].delta_bar.iloc[0] < t]
                th[(s, t)] = "--" if not hit else str(int(c[(c.grid == hit[0]) & (c.method == s)].n.iloc[0])) + (" (full)" if hit[0] == "full" else "")
        T = pd.DataFrame({"first n < 0.10": {s: th[(s, .10)] for s in ALL}, "first n < 0.09": {s: th[(s, .09)] for s in ALL}}); T.to_csv(f"{R}/thresholds_{tag}.csv"); print("first n below the thresholds (grid step 16):\n", T.to_string())
    rng, rows = np.random.default_rng(0), []
    for g in (64, 128, "full"):
        V = {s: val(g, s).groupby(P.PATH).mean() for s in ALL}
        for o in ("al", "almax", "svts", "sde", "random"):
            d = (V["hybrid"] - V[o]).values; bs = [d[rng.integers(0, 32, 32)].mean() for _ in range(10000)]
            rows.append(dict(n=str(g), comparison=f"hybrid - {o}", mean_diff=d.mean(), ci95_low=np.quantile(bs, .025), ci95_high=np.quantile(bs, .975), paths_hybrid_better=int((d < 0).sum()), of=32, wilcoxon_p=wilcoxon(d).pvalue))
    X = pd.DataFrame(rows); X["bonf_5_methods"] = X.wilcoxon_p < 0.05 / 5; X["bonf_15_all"] = X.wilcoxon_p < 0.05 / len(X); X.to_csv(f"{R}/path_level_tests.csv", index=False)
    print("\n== path-level tests (32 path numbers, mean over 50 FPGAs; negative = hybrid better)"); print(X.round(5).to_string(index=False))
    S = pd.concat([pd.read_csv(f"{R}/selection_timing.csv"), pd.read_csv(f"{R}/selection_timing_almax.csv")]); F = pd.read_csv(f"{R}/final_gpr_timing.csv")
    print("\n== selection timing on FPGA 1 (single thread, sequential), per path"); print(S.groupby("method")[["wall_s", "cpu_s", "gpr_fits", "n_sites"]].agg(["mean", "min", "max"]).round(2).to_string())
    print("\n== final additive GPR (200 iterations) on FPGAs 2, 26, 50 x 32 paths"); print(F[["n_train", "n_test", "fit_s", "predict_s"]].agg(["mean", "std", "max"]).round(3).to_string())
    print("\ncommon test set: mean", round(P.n_test.mean()), "sites per path; hardware:", hardware())


if __name__ == "__main__":
    if sys.argv[1] == "time": time_stage()
    elif sys.argv[1] == "timemax": time_max()
    else:
        with Pool(8) as pool: pool.map(evaluate, [(f, p) for f in range(1, 51) for p in range(1, 33) if not os.path.exists(f"{R}/path_f{f}_p{p}.csv")], chunksize=1)
        analyse()
