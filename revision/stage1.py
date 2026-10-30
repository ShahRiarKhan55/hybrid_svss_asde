"""Stage 1: fair baseline table on TUNE (FPGAs 2-26, plus FPGA 1 as secondary). TEST FPGAs are not touched (Stage 4).
python stage1.py   -> selection on FPGA 1 (160 tasks), evaluation (26 FPGAs x 32 paths), tables."""
from lib import *
import time
from multiprocessing import Pool
OUT = f"{HERE}/stage1"; os.makedirs(OUT, exist_ok=True)
METHODS = ["al", "almax", "svts_a", "svts_b", "sde", "random", "maximin", "hybrid"]


def sel_task(a):
    seed, p = a; f = f"{OUT}/orders_s{seed}_p{p}.csv"
    if os.path.exists(f): return
    pool, _ = heldout(p); X, Y = pool[["X", "Y"]].values.astype(float), pool.Value.values; n = len(pool); S0 = start_set(n, seed, p); orders, info = {}, []

    def go(name, fn):
        STATS.pop(name, None); t = time.perf_counter(); r = fn(); w = time.perf_counter() - t
        order, extra = (r if isinstance(r, tuple) else (r, None)); assert len(order) == T == len(set(order))
        orders[name] = order; c = STATS[name]; info.append(dict(seed=seed, PATH=p, method=name, wall_s=w, gpr_fits=c["fits"], asde_accepted=extra, **{k: v for k, v in c.items() if k != "fits"}))
    go("random", lambda: rand_order(n, seed, p)); go("maximin", lambda: maximin(X, S0, T)); go("sde", lambda: sde_extend(X, T, seed))
    go("al", lambda: al(X, Y, S0, "min", "al")); go("almax", lambda: al(X, Y, S0, "max", "almax"))
    go("svts_a", lambda: svts(X, Y, S0, 0.7, 0.3, "svts_a")); go("svts_b", lambda: svts(X, Y, S0, 0.3, 0.7, "svts_b"))
    def hyb():
        anc = svss(X, Y, S0, A, "hybrid"); add, acc = asde(X, Y, anc, T - A, seed); return anc + add, acc
    go("hybrid", hyb)
    pd.DataFrame(info).to_csv(f"{OUT}/cost_s{seed}_p{p}.csv", index=False)
    pd.concat(pd.DataFrame(dict(method=m, rank=range(1, T + 1), X=pool.X.values[o], Y=pool.Y.values[o])) for m, o in orders.items()).to_csv(f + ".tmp", index=False); os.replace(f + ".tmp", f)
    print(time.strftime("%H:%M:%S"), "selected", a, flush=True)


def eval_task(a):
    f, p = a; out = f"{OUT}/eval_f{f}_p{p}.csv"
    if os.path.exists(out): return
    orders = {}
    for s in SEEDS:
        O = pd.read_csv(f"{OUT}/orders_s{s}_p{p}.csv")
        for m, g in O.groupby("method"): orders[(m, s)] = list(zip(g.sort_values("rank").X, g.sort_values("rank").Y))
    pd.DataFrame(evaluate(f, p, orders, BUDGETS, "eval")).to_csv(out + ".tmp", index=False); os.replace(out + ".tmp", out); print(time.strftime("%H:%M:%S"), "eval", a, flush=True)


def report():
    pd.set_option("display.width", 250)
    R = pd.concat(pd.read_csv(f"{OUT}/eval_f{f}_p{p}.csv") for f in [1] + TUNE for p in range(1, 33)); R.to_csv(f"{OUT}/eval_all.csv", index=False)
    a, b = agg(R, TUNE).loc["svts_a"]["mean"].mean(), agg(R, TUNE).loc["svts_b"]["mean"].mean()
    pick = "svts_a" if a <= b else "svts_b"; pd.DataFrame([dict(svts_a_alpha0_7_beta0_3=a, svts_b_alpha0_3_beta0_7=b, chosen_on_TUNE_only=pick)]).to_csv(f"{OUT}/svts_choice.csv", index=False)
    print("SVTS weights, mean TUNE error over budgets: (0.7 sim, 0.3 var) %.4f | (0.3 sim, 0.7 var) %.4f -> chosen %s" % (a, b, pick))
    rows = []
    for tag, fp in (("TUNE 2-26", TUNE), ("TUNE + FPGA 1", [1] + TUNE)):
        for col, tn in (("err", "held-out 30%"), ("err_old", "old (all non-selected)")):
            g = agg(R, fp, col).reset_index(); rows.append(g.assign(fpgas=tag, test=tn))
    G = pd.concat(rows); G["method"] = G.method.replace({pick: "svts", "svts_a" if pick == "svts_b" else "svts_b": "svts_other"}); G.to_csv(f"{OUT}/results.csv", index=False)
    for (fp, tn), g in G.groupby(["fpgas", "test"]):
        t = g.pivot(index="method", columns="budget", values=["mean", "sd"]); print(f"\n== {fp} | test: {tn}  (mean error, SD over 5 seeds)")
        print(pd.DataFrame({b: t["mean"][b].round(4).astype(str) + " ± " + t["sd"][b].round(4).astype(str) for b in BUDGETS}).reindex(["al", "almax", "svts", "svts_other", "sde", "random", "maximin", "hybrid"]).to_string())
    C = pd.concat(pd.read_csv(f"{OUT}/cost_s{s}_p{p}.csv") for s in SEEDS for p in range(1, 33)).fillna(0); C.to_csv(f"{OUT}/cost_all.csv", index=False)
    cc = C.groupby("method").agg(wall_s_per_path=("wall_s", "mean"), gpr_fits_per_path=("gpr_fits", "mean")); cc["asde_accepted"] = C[C.method == "hybrid"].asde_accepted.mean() if "hybrid" in cc.index else 0
    print("\n== selection cost (8 workers in parallel, laptop: indicative only; Stage 5 times it properly)"); print(cc.round(2).to_string())
    prm = ["ls0", "ls1", "var0", "var1", "noise"]; bs = C.groupby("method")[[f"{q}_{e}" for q in prm for e in ("lo", "hi") if f"{q}_{e}" in C]].sum().div(C.groupby("method").gpr_fits.sum().replace(0, np.nan), axis=0)
    print("\n== share of selection-loop GPR fits with a hyperparameter at its bound (lo = lower, hi = upper)"); print((100 * bs).round(2).dropna(how="all").to_string()); bs.to_csv(f"{OUT}/bounds_selection.csv")
    be = R.groupby("method")[[f"{q}_{e}" for q in prm for e in ("lo", "hi")]].mean(); print("\n== share of predictor GPR fits (evaluation) with a hyperparameter at its bound"); print((100 * be).round(2).to_string()); be.to_csv(f"{OUT}/bounds_eval.csv")
    print("\nn_train per budget (mean):", R.groupby("budget").n_train.mean().round(2).to_dict(), "| n_test per path (mean):", round(R.n_test.mean(), 1))


if __name__ == "__main__":
    with Pool(8) as pool:
        pool.map(sel_task, [(s, p) for s in SEEDS for p in range(1, 33) if not os.path.exists(f"{OUT}/orders_s{s}_p{p}.csv")], chunksize=1)
        pool.map(eval_task, [(f, p) for f in [1] + TUNE for p in range(1, 33) if not os.path.exists(f"{OUT}/eval_f{f}_p{p}.csv")], chunksize=1)
    report()
