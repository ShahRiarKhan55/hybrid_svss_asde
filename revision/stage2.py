"""Stage 2: factorial ablation at the 10% budget (317 sites, 95 stage-1 + 222 stage-2), TUNE FPGAs only, same protocol as Stage 1.
Stage-1 choice {random, SVSS} x stage-2 choice {random fill, fixed SDE, adaptive ASDE}, plus SVSS-only, AL-max, ASDE-only (random seed / maximin seed).
Reuses the Stage 1 orders (random, hybrid = SVSS+ASDE, AL-max, SDE). python stage2.py"""
from lib import *
import time
from multiprocessing import Pool
S1 = f"{HERE}/stage1"; OUT = f"{HERE}/stage2"; os.makedirs(OUT, exist_ok=True)
NEW = ["rand+SDE", "rand+ASDE", "SVSS+rand", "SVSS+SDE", "SVSS-only", "ASDE-only(rand seed)", "ASDE-only(maximin seed)"]
NAMES = {"random": "rand+rand", "hybrid": "SVSS+ASDE (hybrid)", "almax": "AL-max (1% seed)", "sde": "SDE (Stage 1 reference)"}


def sel_task(a):
    seed, p = a; f = f"{OUT}/orders_s{seed}_p{p}.csv"
    if os.path.exists(f): return
    pool, _ = heldout(p); X, Y = pool[["X", "Y"]].values.astype(float), pool.Value.values; n = len(pool); S0 = start_set(n, seed, p)
    O = pd.read_csv(f"{S1}/orders_s{seed}_p{p}.csv"); cid = {c: i for i, c in enumerate(zip(pool.X, pool.Y))}
    o1 = lambda m: [cid[c] for c in zip(O[O.method == m].sort_values("rank").X, O[O.method == m].sort_values("rank").Y)]
    anc_r, anc_s = o1("random")[:A], o1("hybrid")[:A]; perm = [int(i) for i in np.random.default_rng([seed, p, 77]).permutation(n)]; orders, info = {}, []

    def go(name, fn):
        STATS.pop(name, None); t = time.perf_counter(); r = fn(); w = time.perf_counter() - t
        order, acc = (r if isinstance(r, tuple) else (r, None)); assert len(order) == T == len(set(order)); orders[name] = order
        info.append(dict(seed=seed, PATH=p, cell=name, wall_s=w, gpr_fits=STATS[name]["fits"], asde_accepted=acc, **{k: v for k, v in STATS[name].items() if k != "fits"}))
    go("rand+SDE", lambda: sde_extend(X, T, seed, start=anc_r))
    go("rand+ASDE", lambda: (lambda r: (anc_r + r[0], r[1]))(asde(X, Y, anc_r, T - A, seed)))
    go("SVSS+rand", lambda: anc_s + [i for i in perm if i not in set(anc_s)][:T - A])
    go("SVSS+SDE", lambda: sde_extend(X, T, seed, start=anc_s))
    go("SVSS-only", lambda: svss(X, Y, S0, T, "SVSS-only"))
    go("ASDE-only(rand seed)", lambda: (lambda r: (S0 + r[0], r[1]))(asde(X, Y, S0, T - STEP, seed)))
    mm = maximin(X, S0[:1], STEP); go("ASDE-only(maximin seed)", lambda: (lambda r: (mm + r[0], r[1]))(asde(X, Y, mm, T - STEP, seed)))
    pd.DataFrame(info).to_csv(f"{OUT}/cost_s{seed}_p{p}.csv", index=False)
    pd.concat(pd.DataFrame(dict(method=m, rank=range(1, T + 1), X=pool.X.values[o], Y=pool.Y.values[o])) for m, o in orders.items()).to_csv(f + ".tmp", index=False); os.replace(f + ".tmp", f)


def eval_task(a):
    f, p = a; out = f"{OUT}/eval_f{f}_p{p}.csv"
    if os.path.exists(out): return
    orders = {}
    for s in SEEDS:
        O = pd.read_csv(f"{OUT}/orders_s{s}_p{p}.csv")
        for m, g in O.groupby("method"): orders[(m, s)] = list(zip(g.sort_values("rank").X, g.sort_values("rank").Y))
    pd.DataFrame(evaluate(f, p, orders, [T], "eval2")).to_csv(out + ".tmp", index=False); os.replace(out + ".tmp", out); print(time.strftime("%H:%M:%S"), "eval", a, flush=True)


def boot(d, rng): return np.quantile([d[rng.integers(0, len(d), len(d))].mean() for _ in range(10000)], [.025, .975])


def report():
    pd.set_option("display.width", 250)
    R1 = pd.read_csv(f"{S1}/eval_all.csv"); R1 = R1[(R1.budget == T) & R1.method.isin(NAMES)]
    R = pd.concat([R1, pd.concat(pd.read_csv(f"{OUT}/eval_f{f}_p{p}.csv") for f in [1] + TUNE for p in range(1, 33))]); R["method"] = R.method.replace(NAMES); R.to_csv(f"{OUT}/eval_all.csv", index=False)
    order = ["rand+rand", "rand+SDE", "rand+ASDE", "SVSS+rand", "SVSS+SDE", "SVSS+ASDE (hybrid)", "SVSS-only", "AL-max (1% seed)", "ASDE-only(rand seed)", "ASDE-only(maximin seed)", "SDE (Stage 1 reference)"]
    rows = []
    for tag, fp in (("TUNE 2-26", TUNE), ("TUNE + FPGA 1", [1] + TUNE)):
        for col, tn in (("err", "held-out 30%"), ("err_old", "old (all non-selected)")): rows.append(agg(R, fp, col).reset_index().assign(fpgas=tag, test=tn))
    G = pd.concat(rows); G.to_csv(f"{OUT}/results.csv", index=False)
    for (fp, tn), g in G.groupby(["fpgas", "test"]): print(f"\n== {fp} | test: {tn} | budget 317 | mean error ± SD over 5 seeds"); print((g.set_index("method")["mean"].round(4).astype(str) + " ± " + g.set_index("method")["sd"].round(4).astype(str)).reindex(order).to_string())
    P = R[R.FPGA.isin(TUNE)].groupby(["method", "PATH"]).err.mean().unstack(0)  # path-level: mean over seeds and TUNE FPGAs
    rng, rows = np.random.default_rng(0), []

    def con(name, d):
        d = np.asarray(d); lo, hi = boot(d, rng); rows.append(dict(contrast=name, mean_diff=d.mean(), ci95_low=lo, ci95_high=hi, paths_first_better=int((d < 0).sum()), of=32, ci_excludes_0=bool(lo > 0 or hi < 0)))
    H = "SVSS+ASDE (hybrid)"
    for o in order:
        if o != H: con(f"hybrid - {o}", P[H] - P[o])
    con("AL-max - SVSS-only", P["AL-max (1% seed)"] - P["SVSS-only"]); con("ASDE-only(maximin seed) - ASDE-only(rand seed)", P["ASDE-only(maximin seed)"] - P["ASDE-only(rand seed)"])
    con("main effect stage 1: SVSS - random (mean of 3 stage-2 choices)", ((P["SVSS+rand"] - P["rand+rand"]) + (P["SVSS+SDE"] - P["rand+SDE"]) + (P[H] - P["rand+ASDE"])) / 3)
    con("main effect stage 2: ASDE - SDE-fixed (mean of 2 stage-1 choices)", ((P["rand+ASDE"] - P["rand+SDE"]) + (P[H] - P["SVSS+SDE"])) / 2)
    con("main effect stage 2: ASDE - random fill", ((P["rand+ASDE"] - P["rand+rand"]) + (P[H] - P["SVSS+rand"])) / 2)
    con("main effect stage 2: SDE-fixed - random fill", ((P["rand+SDE"] - P["rand+rand"]) + (P["SVSS+SDE"] - P["SVSS+rand"])) / 2)
    con("interaction: (SVSS - random) under ASDE minus under SDE-fixed", (P[H] - P["rand+ASDE"]) - (P["SVSS+SDE"] - P["rand+SDE"]))
    Cn = pd.DataFrame(rows); Cn.to_csv(f"{OUT}/contrasts.csv", index=False); print("\n== paired differences (negative = first is better), path-level (32 paths), TUNE, primary test, 95% bootstrap CI"); print(Cn.round(5).to_string(index=False))
    C = pd.concat(pd.read_csv(f"{OUT}/cost_s{s}_p{p}.csv") for s in SEEDS for p in range(1, 33)).fillna(0)
    print("\n== selection cells: GPR fits per path, share of sites passing the ASDE rule"); print(C.groupby("cell").agg(gpr_fits=("gpr_fits", "mean"), asde_accepted=("asde_accepted", "mean")).round(1).to_string())
    prm = [f"{q}_{e}" for q in ("ls0", "ls1", "var0", "var1", "noise") for e in ("lo", "hi")]; be = R.groupby("method")[prm].mean()
    print("\n== share of predictor fits with a hyperparameter at its bound (%)"); print((100 * be).round(2).to_string()); be.to_csv(f"{OUT}/bounds_eval.csv")
    bs = C.groupby("cell")[[c for c in prm if c in C]].sum().div(C.groupby("cell").gpr_fits.sum().replace(0, np.nan), axis=0); print("\n== share of selection-loop fits (cells with fits) with a hyperparameter at its bound (%)"); print((100 * bs).round(2).dropna(how="all").to_string())


if __name__ == "__main__":
    with Pool(8) as pool:
        pool.map(sel_task, [(s, p) for s in SEEDS for p in range(1, 33) if not os.path.exists(f"{OUT}/orders_s{s}_p{p}.csv")], chunksize=1)
        pool.map(eval_task, [(f, p) for f in [1] + TUNE for p in range(1, 33) if not os.path.exists(f"{OUT}/eval_f{f}_p{p}.csv")], chunksize=1)
    report()
