"""Stage 5: cost of the final designs. One process, one thread (OMP/MKL/OpenBLAS = 1, set in lib), sequential, 3 repeats, median.
Selection on FPGA 1 for all 32 paths (seed 1), data loading excluded; final GPR fit and prediction on FPGAs 2, 14, 26 (TUNE; timing only) x 32 paths.
CPU frequency counters, utilisation and the power plan are recorded at every path block. python stage5.py [run|report]"""
from lib import *
import time, subprocess, sys
OUT = f"{HERE}/stage5"; os.makedirs(OUT, exist_ok=True)
METHODS = ["al", "almax", "svts", "sde", "random", "maximin", "hybrid"]; REPS = 3


def ps(q): return subprocess.run(["powershell", "-NoProfile", "-Command", q], capture_output=True, text=True).stdout.strip()


def sysinfo():
    c = ps("$c=Get-Counter '\\Processor Information(_Total)\\Processor Frequency','\\Processor Information(_Total)\\% Processor Utility','\\Processor Information(_Total)\\% of Maximum Frequency'; ($c.CounterSamples | % { [math]::Round($_.CookedValue,1) }) -join ','").split(",")
    plan = ps("powercfg /getactivescheme")
    try: f, u, pm = (float(x) for x in c)
    except Exception: f = u = pm = np.nan
    return dict(freq_mhz=f, cpu_utility_pct=u, pct_of_max_freq=pm, power_plan=plan.split("(")[-1].rstrip(")") if "(" in plan else plan)


def methods(X, Y, n, p, S0):
    return {"al": lambda: al(X, Y, S0, "min", "al"), "almax": lambda: al(X, Y, S0, "max", "almax"), "svts": lambda: svts(X, Y, S0, 0.3, 0.7, "svts"), "sde": lambda: sde_extend(X, T, 1),
            "random": lambda: rand_order(n, 1, p), "maximin": lambda: maximin(X, S0, T), "hybrid": lambda: (lambda a: a + asde(X, Y, a, T - A, 1)[0])(svss(X, Y, S0, A, "hybrid"))}


def run_selection():
    f = f"{OUT}/selection_timing.csv"; done = set(); rows = []
    if os.path.exists(f): rows = pd.read_csv(f).to_dict("records"); done = {(r["rep"], r["PATH"], r["method"]) for r in rows}
    for rep in range(1, REPS + 1):
        for p in range(1, 33):
            if all((rep, p, m) in done for m in METHODS): continue
            pool, _ = heldout(p); X, Y = pool[["X", "Y"]].values.astype(float), pool.Value.values; n = len(pool); S0 = start_set(n, 1, p); si = sysinfo()
            for m, fn in methods(X, Y, n, p, S0).items():
                if (rep, p, m) in done: continue
                STATS.pop(m, None); t0, c0 = time.perf_counter(), time.process_time(); order = fn(); w, c = time.perf_counter() - t0, time.process_time() - c0
                assert len(order) == T == len(set(order)); rows.append(dict(rep=rep, PATH=p, method=m, wall_s=w, cpu_s=c, gpr_fits=STATS[m]["fits"], **si))
            pd.DataFrame(rows).to_csv(f, index=False); print(time.strftime("%H:%M:%S"), "rep", rep, "path", p, {r["method"]: round(r["wall_s"], 1) for r in rows[-7:]}, flush=True)


def run_final():
    f = f"{OUT}/final_timing.csv"; rows = pd.read_csv(f).to_dict("records") if os.path.exists(f) else []; done = {(r["rep"], r["FPGA"], r["PATH"]) for r in rows}
    for rep in range(1, REPS + 1):
        for fp in (2, 14, 26):
            for p in range(1, 33):
                if (rep, fp, p) in done: continue
                df = hy.load_dat(fp, p); O = pd.read_csv(f"{HERE}/stage1/orders_s1_p{p}.csv"); h = O[O.method == "hybrid"].sort_values("rank"); pos = set(zip(h.X, h.Y))
                mk = np.array([c in pos for c in zip(df.X, df.Y)]); X, v = df[["X", "Y"]].values.astype(float), df.Value.values
                t0 = time.perf_counter(); m = fit(X[mk], v[mk], "final"); t1 = time.perf_counter(); m.predict(X[~mk]); t2 = time.perf_counter()
                rows.append(dict(rep=rep, FPGA=fp, PATH=p, n_train=int(mk.sum()), n_predict=int((~mk).sum()), fit_s=t1 - t0, predict_s=t2 - t1, **sysinfo()))
        pd.DataFrame(rows).to_csv(f, index=False); print(time.strftime("%H:%M:%S"), "final timing rep", rep, flush=True)


def report():
    S = pd.read_csv(f"{OUT}/selection_timing.csv"); F = pd.read_csv(f"{OUT}/final_timing.csv")
    M = S.groupby(["method", "PATH"]).agg(wall=("wall_s", "median"), cpu=("cpu_s", "median"), fits=("gpr_fits", "first"), rep_spread=("wall_s", lambda s: (s.max() - s.min()) / s.median() if s.median() > 0 else 0)).reset_index()
    names = {"al": "AL", "almax": "AL-max", "svts": "SVTS (0.3 sim, 0.7 var)", "sde": "SDE", "random": "Random", "maximin": "Maximin", "hybrid": "Hybrid (SVSS 95 + ASDE 222)"}
    ct = ["| method | GPR fits per path | selection time per path, s (median of 3 repeats; mean over 32 paths) | median over paths | min – max over paths | CPU time, s | relative to hybrid | mean repeat-to-repeat spread |", "|---|---|---|---|---|---|---|---|"]
    hb = M[M.method == "hybrid"].wall.mean()
    for m in METHODS:
        g = M[M.method == m]; ct.append(f"| {names[m]} | {g.fits.mean():.1f} | {g.wall.mean():.3f} | {g.wall.median():.3f} | {g.wall.min():.3f} – {g.wall.max():.3f} | {g.cpu.mean():.3f} | {g.wall.mean() / hb:.1f}x | {100 * g.rep_spread.mean():.0f}% |")
    FM = F.groupby(["FPGA", "PATH"]).agg(fit=("fit_s", "median"), pred=("predict_s", "median"), ntr=("n_train", "first"), npr=("n_predict", "first"))
    sysm = S.drop_duplicates(["rep", "PATH"])
    md = f"""# Stage 5 – Cost (final designs)

## Protocol
- One Python process, BLAS/OpenMP threads = 1 (set before importing NumPy; checked in Stage 0), sequential, nothing else started by me during the run. Selection on FPGA 1, all 32 paths, seed 1, held-out pool as in Stages 1–4; data loading and the held-out set are outside the timer. Each (path, method) was timed 3 times (repeats interleaved: repeat 1 over all paths, then repeat 2, then repeat 3); the table uses the median of the 3 repeats per path, then the mean over the 32 paths.
- Same additive bounded RBF GPR (`lib.fit`, up to 200 optimiser iterations) inside every selection loop and for the final fit. Final fit and prediction: hybrid training set (317 sites) on FPGAs 2, 14 and 26 (TUNE FPGAs, timing only) × 32 paths, 3 repeats, median; prediction on all non-training sites of the FPGA.
- Machine: i7-12700H (14 cores / 20 threads), 31.7 GB RAM, Windows 11, Python 3.14.7, NumPy 2.3.5, GPy 1.14.2.
- Power plan during the run: {', '.join(sorted(S.power_plan.dropna().unique()))}. I did not change it (the machine was not switched to High performance). Windows reported processor frequency {sysm.freq_mhz.min():.0f}–{sysm.freq_mhz.max():.0f} MHz (mean {sysm.freq_mhz.mean():.0f}), {sysm.pct_of_max_freq.min():.0f}–{sysm.pct_of_max_freq.max():.0f}% of maximum frequency, whole-machine CPU utilisation {sysm.cpu_utility_pct.min():.0f}–{sysm.cpu_utility_pct.max():.0f}% (mean {sysm.cpu_utility_pct.mean():.0f}%) at the start of each path block. Background services (vendor utilities, the Claude app) were running, so the machine was not idle.

## Selection cost
{chr(10).join(ct)}

## Final fit and prediction (317 training sites)
| quantity | mean over 96 (FPGA, path) | min – max |
|---|---|---|
| final GPR fit, s | {FM.fit.mean():.3f} | {FM.fit.min():.3f} – {FM.fit.max():.3f} |
| prediction on {FM.npr.mean():.0f} sites, s | {FM.pred.mean():.3f} | {FM.pred.min():.3f} – {FM.pred.max():.3f} |

The final fit and prediction are the same for every method (same GPR, same number of training sites).

## Notes
- Times are wall-clock seconds on a laptop; the repeat-to-repeat spread column shows how reproducible they were. Absolute values depend on the machine; the ratios and the fit counts are the portable numbers.
- The hybrid's 2 fits are on 32 and 64 sites; the AL, AL-max and SVTS loops refit after every added site (285 fits, on 32 to 316 sites).
- SDE, Random and Maximin use no GPR; their times are for the Python selection loops only.
"""
    open(f"{HERE}/stage5_report.md", "w", encoding="utf8").write(md); print("report written")


if __name__ == "__main__":
    if sys.argv[1:] == ["report"]: report()
    else: run_selection(); run_final(); report()
