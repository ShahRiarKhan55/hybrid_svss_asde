"""python timing.py [repeats=5] [outdir=Timing dir] [methods=AL,SVTS,Hybrid,SDE,Final] [paths=1,10,20]"""
import os, sys
for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"): os.environ[v] = "1"
import json, time, random, platform, subprocess, traceback, tracemalloc, numpy as np, pandas as pd
import GPy, sklearn
from sklearn.svm import SVR
from sklearn.neural_network import MLPRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GridSearchCV, KFold

H = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); A = sys.argv
REP = int(A[1]) if len(A) > 1 else 5
OUT = A[2] if len(A) > 2 else f"{H}/Timing"
METH = (A[3] if len(A) > 3 else "AL,SVTS,Hybrid,SDE,Final").split(",")
PATHS = [int(p) for p in (A[4] if len(A) > 4 else "1,10,20").split(",")]
os.makedirs(OUT, exist_ok=True); RAW = f"{OUT}/timing_raw.csv"; CALIB, TEST_FPGA, SEED = 1, 2, 42


def ns(nb, cut=None):
    s = "".join(json.load(open(f"{H}/{nb}.ipynb", encoding="utf8"))["cells"][0]["source"])
    d = {"__name__": "nb"}; exec(s[:s.index(cut)] if cut else s, d); return d


al, sv, sd = ns("raw_al", "calib_positions = {}"), ns("svts", "calib_positions = {}"), ns("sde", "calib_positions = {}")
hy, dg = ns("hybrid_svss_asde"), ns("hybrid_svss_asde_dgpr")
NFIT = [0]; _opt = GPy.core.GP.optimize


def _count(self, *a, **k): NFIT[0] += 1; return _opt(self, *a, **k)


GPy.core.GP.optimize = _count


def grid(g, X, y, te):
    t = time.perf_counter(); g.fit(X, y); t1 = time.perf_counter(); g.predict(te[["X", "Y"]].values)
    return dict(cv=t1 - t - g.refit_time_, refit=g.refit_time_, n_fits=len(g.cv_results_["params"]) * g.n_splits_ + 1)


def svr(tr, te):
    X, y = tr[["X", "Y"]].values, tr["Value"].values; ym, ys = y.mean(), y.std() + 1e-9
    g = GridSearchCV(make_pipeline(StandardScaler(), SVR(kernel="rbf")), {"svr__C": [1, 10, 100], "svr__gamma": [0.1, 0.3, 1, 3], "svr__epsilon": [0.01, 0.1]}, cv=KFold(5, shuffle=True, random_state=SEED))
    return grid(g, X, (y - ym) / ys, te)


def mlp(tr, te):
    X, y = tr[["X", "Y"]].values, tr["Value"].values; ym, ys = y.mean(), y.std() + 1e-9
    g = GridSearchCV(make_pipeline(StandardScaler(), MLPRegressor(activation="tanh", max_iter=2000, random_state=SEED)), {"mlpregressor__hidden_layer_sizes": [(32, 32), (64, 64)], "mlpregressor__alpha": [1e-3, 1e-1]}, cv=KFold(3, shuffle=True, random_state=SEED), n_jobs=1)
    return grid(g, X, (y - ym) / ys, te)


REG = {"GPR": (lambda tr, te: hy["gpr"](tr, te), 1), "SVR": (svr, None), "MLP": (mlp, None), "DGPR": (lambda tr, te: dg["dgpr"](tr, te), 1)}
if not os.path.exists(RAW): open(RAW, "w").write("method,regressor,path,repeat,seconds,n_fits,peak_MB\n")
done = set(map(tuple, pd.read_csv(RAW)[["method", "regressor", "path", "repeat"]].astype(str).values))


def log(m, r, p, k, s, n, mb=""):
    open(RAW, "a").write(f"{m},{r},{p},{k},{s:.6f},{n},{mb}\n"); done.add((m, r, str(p), str(k)))


def todo(m, r, p, k): return (m, r, str(p), str(k)) not in done


def sel_fn(m, df):
    n = len(df)
    if m == "AL": return lambda: al["select_train_coords"](df)
    if m == "SVTS": return lambda: sv["select_train_coords"](df)
    if m == "SDE": return lambda: sd["sde_select"](df.X.values, df.Y.values, int(round(0.10 * n)), random.Random(SEED))
    if m == "Hybrid": return lambda: hy["adaptive_sde"](df, hy["svts_select"](df))


def hybrid_pos(df):
    a = hy["svts_select"](df); return set(map(tuple, df.loc[a + hy["adaptive_sde"](df, a), ["X", "Y"]].values))


def timed(f):
    NFIT[0] = 0; t = time.perf_counter(); r = f(); return time.perf_counter() - t, NFIT[0], r


def guard(tag, f):
    try: f()
    except Exception:
        print("ERROR", tag, flush=True); open(f"{OUT}/errors.log", "a").write(f"{tag}\n{traceback.format_exc()}\n")


if __name__ == "__main__":
    P = {p: hy["load_dat"](CALIB, p) for p in PATHS}; T = {p: hy["load_dat"](TEST_FPGA, p) for p in PATHS}; warmed = set()
    for m in [x for x in METH if x != "Final"]:
        for p in PATHS:
            def run():
                f = sel_fn(m, P[p])
                if m not in warmed: f(); warmed.add(m)
                for k in range(1, REP + 1):
                    if todo(m, "selection", p, k):
                        s, n, _ = timed(f); log(m, "selection", p, k, s, n); print(m, p, k, f"{s:.2f}s fits={n}", flush=True)
            guard(f"{m} {p}", run)
    if "Final" in METH:
        for p in PATHS:
            pos = hybrid_pos(P[p]); d = T[p]; tr = d[[(x, y) in pos for x, y in zip(d.X, d.Y)]]; te = d.drop(tr.index)
            print("path", p, "train", len(tr), "test", len(te), flush=True)
            for r, (f, nf) in REG.items():
                def run():
                    if r not in warmed: f(tr, te); warmed.add(r)
                    for k in range(1, REP + 1):
                        if not todo("Final", r, p, k): continue
                        s, n, out = timed(lambda: f(tr, te)); nfit = out["n_fits"] if isinstance(out, dict) else nf
                        log("Final", r, p, k, s, nfit)
                        if isinstance(out, dict): log("Final_cv", r, p, k, out["cv"], nfit); log("Final_refit", r, p, k, out["refit"], 1)
                        print(r, p, k, f"{s:.2f}s", flush=True)
                        if s > 600: print(r, "run > 10 min, remaining repeats skipped", flush=True); break
                    if todo("Final_mem", r, p, 1) and not todo("Final", r, p, 1):
                        tracemalloc.start(); f(tr, te); mb = tracemalloc.get_traced_memory()[1] / 2**20; tracemalloc.stop(); log("Final_mem", r, p, 1, 0, 0, f"{mb:.1f}")
                guard(f"Final {r} {p}", run)
    D = pd.read_csv(RAW); S = D[D.method != "Final_mem"].groupby(["method", "regressor", "path"]).agg(med=("seconds", "median"), lo=("seconds", "min"), hi=("seconds", "max"), n_fits=("n_fits", "median")).reset_index()
    S.groupby(["method", "regressor"]).agg(median_s_per_path=("med", "mean"), min_s=("lo", "min"), max_s=("hi", "max"), n_fits_per_path=("n_fits", "mean")).reset_index().to_csv(f"{OUT}/timing_summary.csv", index=False)
    M = D[D.method == "Final_mem"].groupby("regressor").peak_MB.mean()
    try: cpu = subprocess.check_output(["powershell", "-NoProfile", "-Command", "(Get-CimInstance Win32_Processor).Name"], text=True).strip()
    except Exception: cpu = platform.processor()
    try: ram = f"{int(subprocess.check_output(['powershell', '-NoProfile', '-Command', '(Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory'], text=True)) / 2**30:.1f} GB"
    except Exception: ram = "?"
    print(f"\nMachine: {cpu}, {ram} RAM, Python {platform.python_version()}, numpy {np.__version__}, GPy {GPy.__version__}, scikit-learn {sklearn.__version__}, deepgp (git install), 1 thread, {platform.platform()}")
    SS = pd.read_csv(f"{OUT}/timing_summary.csv").set_index(["method", "regressor"])
    print("\n| Method | GPR fits per path | Selection s/path (median; min-max) |\n|---|---|---|")
    for m in ("AL", "SVTS", "Hybrid", "SDE"):
        if (m, "selection") in SS.index: x = SS.loc[(m, "selection")]; print(f"| {m} | {x.n_fits_per_path:.0f} | {x.median_s_per_path:.1f} ({x.min_s:.1f}-{x.max_s:.1f}) |")
    print("\n| Regressor | fit+predict s/path (median; min-max) | CV tuning s | refit s | fits | peak MB (tracemalloc) |\n|---|---|---|---|---|---|")
    for r in REG:
        if ("Final", r) in SS.index:
            x = SS.loc[("Final", r)]
            cv = SS.loc[("Final_cv", r)].median_s_per_path if ("Final_cv", r) in SS.index else float("nan")
            rf = SS.loc[("Final_refit", r)].median_s_per_path if ("Final_refit", r) in SS.index else float("nan")
            print(f"| {r} | {x.median_s_per_path:.2f} ({x.min_s:.2f}-{x.max_s:.2f}) | {cv:.2f} | {rf:.3f} | {x.n_fits_per_path:.0f} | {M.get(r, float('nan')):.1f} |")
