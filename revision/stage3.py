"""Stage 3: trimmed candidate list for the hybrid, TUNE FPGAs only (FPGA 1 as secondary). Rule and definitions: stage3_rule_section.md (written before running).
python stage3.py -> selection (160 tasks), evaluation (4 budgets x 5 seeds x 5 candidates), tables, verdict, stage3_report.md"""
from lib import *
import lib, time
from multiprocessing import Pool
S1 = f"{HERE}/stage1"; OUT = f"{HERE}/stage3"; os.makedirs(OUT, exist_ok=True)
CANDS = ["SVSS-cos", "anchors 1%", "anchors 2%", "maximin anchors", "ASDE relax"]


def sel_task(a):
    seed, p = a; f = f"{OUT}/orders_s{seed}_p{p}.csv"
    if os.path.exists(f): return
    pool, _ = heldout(p); X, Y = pool[["X", "Y"]].values.astype(float), pool.Value.values; n = len(pool); S0 = start_set(n, seed, p)
    O = pd.read_csv(f"{S1}/orders_s{seed}_p{p}.csv"); cid = {c: i for i, c in enumerate(zip(pool.X, pool.Y))}; h = O[O.method == "hybrid"].sort_values("rank")
    anc_h = [cid[c] for c in zip(h.X, h.Y)][:A]; orders, info = {}, []

    def go(name, anchors_fn, kw=None):
        STATS.pop(name, None); lib.LAST.clear(); t = time.perf_counter(); anc = anchors_fn(); add, acc = asde(X, Y, anc, T - len(anc), seed, **(kw or {})); w = time.perf_counter() - t
        order = anc + add; assert len(order) == T == len(set(order)); orders[name] = order
        info.append(dict(seed=seed, PATH=p, cand=name, n_anchor=len(anc), wall_s=w, gpr_fits=STATS[name]["fits"], asde_accepted=acc, n_add=T - len(anc), relax_level=lib.LAST.get("relax_level", 0), **{k: v for k, v in STATS[name].items() if k != "fits"}))
    go("SVSS-cos", lambda: svss(X, Y, S0, A, "SVSS-cos", center=False)); go("anchors 1%", lambda: list(S0)); go("anchors 2%", lambda: svss(X, Y, S0, 63, "anchors 2%"))
    go("maximin anchors", lambda: maximin(X, S0, A)); go("ASDE relax", lambda: list(anc_h), dict(relax=True))
    pd.DataFrame(info).to_csv(f"{OUT}/cost_s{seed}_p{p}.csv", index=False)
    pd.concat(pd.DataFrame(dict(method=m, rank=range(1, T + 1), X=pool.X.values[o], Y=pool.Y.values[o])) for m, o in orders.items()).to_csv(f + ".tmp", index=False); os.replace(f + ".tmp", f)


def eval_task(a):
    f, p = a; out = f"{OUT}/eval_f{f}_p{p}.csv"
    if os.path.exists(out): return
    orders = {}
    for s in SEEDS:
        O = pd.read_csv(f"{OUT}/orders_s{s}_p{p}.csv")
        for m, g in O.groupby("method"): orders[(m, s)] = list(zip(g.sort_values("rank").X, g.sort_values("rank").Y))
    pd.DataFrame(evaluate(f, p, orders, BUDGETS, "eval3")).to_csv(out + ".tmp", index=False); os.replace(out + ".tmp", out); print(time.strftime("%H:%M:%S"), "eval", a, flush=True)


def report():
    R1 = pd.read_csv(f"{S1}/eval_all.csv"); R1 = R1[R1.method.isin(["random", "hybrid", "sde"])].replace({"random": "Random", "hybrid": "Hybrid (current)", "sde": "SDE"})
    R = pd.concat([R1, pd.concat(pd.read_csv(f"{OUT}/eval_f{f}_p{p}.csv") for f in [1] + TUNE for p in range(1, 33))]); R.to_csv(f"{OUT}/eval_all.csv", index=False)
    ROWS = ["Random", "SDE", "Hybrid (current)"] + CANDS; G = {}
    for tag, fp in (("TUNE 2-26", TUNE), ("TUNE + FPGA 1", [1] + TUNE)):
        for col, tn in (("err", "held-out"), ("err_old", "old")): G[(tag, tn)] = agg(R, fp, col)

    def tab(tag, tn):
        g = G[(tag, tn)]; r = ["| method | 95 | 158 | 222 | 317 |", "|---|---|---|---|---|"]
        for m in ROWS: r.append(f"| {m} | " + " | ".join(f"{g.loc[(m, b), 'mean']:.4f} ± {g.loc[(m, b), 'sd']:.4f}" for b in BUDGETS) + " |")
        return "\n".join(r)
    P = R[R.FPGA.isin(TUNE)].groupby(["method", "budget", "PATH"]).err.mean(); rng = np.random.default_rng(0); rows = []
    for m in ["Hybrid (current)"] + CANDS:
        for b in BUDGETS:
            for ref in ("Random", "Hybrid (current)"):
                if ref == m: continue
                d = (P[m, b] - P[ref, b]).values; lo, hi = np.quantile([d[rng.integers(0, 32, 32)].mean() for _ in range(10000)], [.025, .975])
                rows.append(dict(candidate=m, budget=b, vs=ref, mean_diff=d.mean(), ci_low=lo, ci_high=hi, paths_candidate_better=int((d < 0).sum()), ci_below_0=bool(hi < 0)))
    D = pd.DataFrame(rows); D.to_csv(f"{OUT}/contrasts.csv", index=False)
    V = D[D.vs == "Random"].groupby("candidate").agg(budgets_ci_below_0=("ci_below_0", "sum")); V["mean_over_budgets"] = [np.mean([G[("TUNE 2-26", "held-out")].loc[(m, b), "mean"] for b in BUDGETS]) for m in V.index]
    V["improvement"] = (V.index != "Hybrid (current)") & (V.budgets_ci_below_0 >= 2); V.to_csv(f"{OUT}/verdict.csv")
    win = V[V.improvement]; frozen = win.mean_over_budgets.idxmin() if len(win) else "Hybrid (current)"
    ct = ["| candidate | budget | candidate − random | 95% CI | paths better (of 32) | CI below 0 |", "|---|---|---|---|---|---|"]
    for _, r in D[D.vs == "Random"].iterrows(): ct.append(f"| {r.candidate} | {r.budget} | {r.mean_diff:+.5f} | [{r.ci_low:+.5f}, {r.ci_high:+.5f}] | {r.paths_candidate_better} | {'yes' if r.ci_below_0 else 'no'} |")
    ch = ["| candidate | budget | candidate − current hybrid | 95% CI | paths better (of 32) |", "|---|---|---|---|---|"]
    for _, r in D[(D.vs == "Hybrid (current)") & (D.candidate != "Hybrid (current)")].iterrows(): ch.append(f"| {r.candidate} | {r.budget} | {r.mean_diff:+.5f} | [{r.ci_low:+.5f}, {r.ci_high:+.5f}] | {r.paths_candidate_better} |")
    vt = ["| candidate | budgets with CI below 0 (of 4) | mean TUNE error over budgets | improvement under the rule |", "|---|---|---|---|"]
    for m, r in V.iterrows(): vt.append(f"| {m} | {int(r.budgets_ci_below_0)} | {r.mean_over_budgets:.4f} | {'YES' if r.improvement else 'no'} |")
    C = pd.concat(pd.read_csv(f"{OUT}/cost_s{s}_p{p}.csv") for s in SEEDS for p in range(1, 33)); C1 = pd.read_csv(f"{S1}/cost_all.csv"); hh = C1[C1.method == "hybrid"]
    cs = ["| design | anchors | GPR fits per path | additions passing the ASDE rule | window relaxation level (mean) |", "|---|---|---|---|---|", f"| Hybrid (current) | 95 | 2.0 | {hh.asde_accepted.mean():.1f} of 222 ({100 * hh.asde_accepted.mean() / 222:.1f}%) | – |"]
    for m, g in C.groupby("cand", sort=False): cs.append(f"| {m} | {int(g.n_anchor.mean())} | {g.gpr_fits.mean():.1f} | {g.asde_accepted.mean():.1f} of {int(g.n_add.mean())} ({100 * g.asde_accepted.mean() / g.n_add.mean():.1f}%) | " + (f"{g.relax_level.mean():.2f}" if m == "ASDE relax" else "–") + " |")
    prm = ["ls0", "ls1", "var0", "var1", "noise"]; nm = ["column lengthscale", "row lengthscale", "column variance", "row variance", "noise"]
    be = R.groupby("method")[[f"{q}_{e}" for q in prm for e in ("lo", "hi")]].mean(); bt = ["| design | " + " | ".join(f"{x} lower / upper" for x in nm) + " |", "|---|" + "---|" * 5]
    for m in ROWS: bt.append(f"| {m} | " + " | ".join(f"{100 * be.loc[m, q + '_lo']:.1f}% / {100 * be.loc[m, q + '_hi']:.1f}%" for q in prm) + " |")
    sb = C[C.gpr_fits > 0].groupby("cand")[[f"{q}_{e}" for q in prm for e in ("lo", "hi") if f"{q}_{e}" in C]].sum().div(C[C.gpr_fits > 0].groupby("cand").gpr_fits.sum(), axis=0)
    st = ["| design (selection-loop fits) | " + " | ".join(f"{x} lower" for x in nm) + " |", "|---|" + "---|" * 5]
    for m in sb.index: st.append(f"| {m} | " + " | ".join(f"{100 * sb.loc[m, q + '_lo']:.1f}%" for q in prm) + " |")
    head = "# Stage 3 – Improving the hybrid (TUNE FPGAs only)\n\n" + open(f"{HERE}/stage3_rule_section.md", encoding="utf8").read().rstrip()
    md = f"""{head}

## Protocol
Same as Stages 1–2: held-out 30% test set, additive bounded RBF GPR (same function inside selection and for prediction), 5 seeds, selection on FPGA 1, TUNE FPGAs 2–26 (25 devices), FPGA 1 only in a secondary table, error = |mu − v| / d_spec averaged test sites → path → FPGA. TEST FPGAs (27–50) were not used. Random, SDE and the current hybrid are the Stage 1 results (same seeds, same test sets).

## Verdict (applying the rule written above)
{chr(10).join(vt)}

**Frozen design: {frozen}.** {"No candidate met the rule, so the current hybrid is frozen unchanged." if frozen == "Hybrid (current)" else "Selected by the pre-written tie-break among the candidates that met the rule."}

## Error by budget (TUNE 2–26, held-out test set; mean ± SD over 5 seeds)
{tab("TUNE 2-26", "held-out")}

## Candidate − random, path-level, 95% bootstrap CI over 32 paths (TUNE; negative = candidate better)
{chr(10).join(ct)}

For reference, the current hybrid against random in the same format is in `stage3/contrasts.csv` (rows with candidate = "Hybrid (current)").

## Candidate − current hybrid (information only, not part of the rule)
{chr(10).join(ch)}

## Secondary tables
TUNE + FPGA 1, held-out test set:

{tab("TUNE + FPGA 1", "held-out")}

TUNE 2–26, paper's old test set (all non-training sites):

{tab("TUNE 2-26", "old")}

## Selection facts
{chr(10).join(cs)}

Hyperparameters at a bound (within 5%), predictor fits, all budgets and FPGAs:

{chr(10).join(bt)}

Selection-loop fits (only designs with SVSS fits):

{chr(10).join(st)}
"""
    open(f"{HERE}/stage3_report.md", "w", encoding="utf8").write(md); print(md[md.index("## Verdict"):])


if __name__ == "__main__":
    with Pool(8) as pool:
        pool.map(sel_task, [(s, p) for s in SEEDS for p in range(1, 33) if not os.path.exists(f"{OUT}/orders_s{s}_p{p}.csv")], chunksize=1)
        pool.map(eval_task, [(f, p) for f in [1] + TUNE for p in range(1, 33) if not os.path.exists(f"{OUT}/eval_f{f}_p{p}.csv")], chunksize=1)
    report()
