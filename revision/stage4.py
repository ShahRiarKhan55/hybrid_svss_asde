"""Stage 4: the one TEST run (FPGAs 27-50). Frozen design = current hybrid. Selection orders are the saved FPGA-1 orders of Stages 1-3 (selection never saw TEST).
python stage4.py -> evaluation (24 FPGAs x 32 paths), tables, tests, stage4_report.md (claims section is added by hand afterwards)."""
from lib import *
import time
from multiprocessing import Pool
from scipy.stats import wilcoxon
S1, S2, S3 = f"{HERE}/stage1", f"{HERE}/stage2", f"{HERE}/stage3"; OUT = f"{HERE}/stage4"; os.makedirs(OUT, exist_ok=True)
H1 = ["al", "almax", "svts_b", "sde", "random", "maximin", "hybrid"]; EXPL = ["anchors 1%", "maximin anchors"]
NEW2 = ["rand+SDE", "rand+ASDE", "SVSS+rand", "SVSS+SDE", "SVSS-only", "ASDE-only(rand seed)", "ASDE-only(maximin seed)"]
NM = {"hybrid": "Hybrid (frozen)", "al": "AL", "almax": "AL-max", "svts_b": "SVTS", "sde": "SDE", "random": "Random", "maximin": "Maximin", "anchors 1%": "anchors 1% (exploratory)", "maximin anchors": "maximin anchors (exploratory)"}
HEAD = ["Hybrid (frozen)", "AL", "AL-max", "SVTS", "SDE", "Random", "Maximin"]; EX = [NM[e] for e in EXPL]; K50 = 1587


def read(root, p, methods):
    o = {}
    for s in SEEDS:
        O = pd.read_csv(f"{root}/orders_s{s}_p{p}.csv")
        for m in methods: g = O[O.method == m].sort_values("rank"); o[(m, s)] = list(zip(g.X, g.Y))
    return o


def eval_task(a):
    f, p = a; out = f"{OUT}/eval_f{f}_p{p}.csv"
    if os.path.exists(out): return
    rows = evaluate(f, p, read(S1, p, H1), BUDGETS, "t4a") + evaluate(f, p, read(S3, p, EXPL), BUDGETS, "t4b") + evaluate(f, p, read(S2, p, NEW2), [T], "t4c")
    pool, held = heldout(p); n = len(pool)  # references: GPR on a random 50% of all sites (1587 of the 70% pool), constant = mean of a random 317-site training set
    o50 = {("gpr50", s): list(zip(pool.X.values[q], pool.Y.values[q])) for s in [1] for q in [np.random.default_rng([s, p, 88]).permutation(n)[:K50]]}
    rows += evaluate(f, p, o50, [K50], "t4d")
    df = hy.load_dat(f, p); ds = df.Value.max() - df.Value.min(); co = list(zip(df.X, df.Y)); cid = {c: i for i, c in enumerate(co)}; it = np.array([c in held for c in co]); v = df.Value.values
    for (m, s), order in read(S1, p, ["random"]).items():
        ix = np.array([cid.get(c, -1) for c in order[:T]]); ix = ix[ix >= 0]; mk = np.zeros(len(df), bool); mk[ix] = True; e = np.abs(v[~mk] - v[mk].mean()) / ds
        rows.append(dict(FPGA=f, PATH=p, method="constant", seed=s, budget=T, n_train=int(mk.sum()), n_test=int(it[~mk].sum()), err=e[it[~mk]].mean(), err_old=e.mean()))
    pd.DataFrame(rows).to_csv(out + ".tmp", index=False); os.replace(out + ".tmp", out); print(time.strftime("%H:%M:%S"), "eval", a, flush=True)


def holm(p):
    p = np.asarray(p); o = np.argsort(p); m = len(p); adj = np.empty(m); run = 0
    for rank, i in enumerate(o): run = max(run, min(1, (m - rank) * p[i])); adj[i] = run
    return adj


def report():
    R = pd.concat(pd.read_csv(f"{OUT}/eval_f{f}_p{p}.csv") for f in TEST for p in range(1, 33)); R["method"] = R.method.replace(NM); R.to_csv(f"{OUT}/eval_all.csv", index=False)
    alias = {"rand+rand": "Random", "SVSS+ASDE": "Hybrid (frozen)", "AL-max (1% seed)": "AL-max", "SDE ref": "SDE"}
    g1 = agg(R, TEST); g1.to_csv(f"{OUT}/results_primary.csv")
    # earlier-stage evaluations of FPGA 1 for the secondary "TEST + FPGA 1" table (headline and exploratory rows only)
    F1 = pd.concat([pd.read_csv(f"{S1}/eval_all.csv").query("FPGA == 1 and method in @H1").replace(NM), pd.read_csv(f"{S3}/eval_all.csv").query("FPGA == 1 and method in @EXPL").replace(NM)])
    RF = pd.concat([R[R.method.isin(HEAD + EX)], F1])
    def tab(g, rows, budgets=BUDGETS):
        t = ["| method | " + " | ".join(str(b) for b in budgets) + " |", "|---|" + "---|" * len(budgets)]
        for m in rows: t.append(f"| {m} | " + " | ".join(f"{g.loc[(m, b), 'mean']:.4f} ± {g.loc[(m, b), 'sd']:.4f}" for b in budgets) + " |")
        return "\n".join(t)
    gF, gO = agg(RF, TEST + [1]), agg(R, TEST, "err_old")
    const, g50 = g1.loc[("constant", T)], g1.loc[("gpr50", K50)]; head = const["mean"] - g50["mean"]
    # path-level values: mean over seeds and the 24 TEST FPGAs
    P = R.groupby(["method", "budget", "PATH"]).err.mean(); rng = np.random.default_rng(0); tests = {}
    def val(m, b): return P[alias.get(m, m), b].values
    def add(key, label, d, first, second):
        tests[key] = dict(comparison=label, d=np.asarray(d), first=first, second=second)
    for b in BUDGETS:
        for m in HEAD[1:]: add(("H", m, b), f"Hybrid − {m} @ {b}", val("Hybrid (frozen)", b) - val(m, b), "hybrid", m)
        for m in EX: add(("E", m, b), f"Hybrid − {m} @ {b}", val("Hybrid (frozen)", b) - val(m, b), "hybrid", m)
    c = lambda m: P[m, T].values
    for m in NEW2:
        if m != "SVSS-only" or True: add(("A", m), f"Hybrid − {m} @ 317", val("Hybrid (frozen)", T) - c(m), "hybrid", m)
    add(("A", "almax-svss"), "AL-max − SVSS-only @ 317", val("AL-max", T) - c("SVSS-only"), "AL-max", "SVSS-only"); add(("A", "mm-rs"), "ASDE-only(maximin seed) − ASDE-only(rand seed) @ 317", c("ASDE-only(maximin seed)") - c("ASDE-only(rand seed)"), "ASDE-only(maximin)", "ASDE-only(rand)")
    H_ = val("Hybrid (frozen)", T); Rd = val("Random", T)
    add(("A", "me1"), "Stage-1 main effect: SVSS − random (mean of 3 stage-2 choices) @ 317", ((c("SVSS+rand") - Rd) + (c("SVSS+SDE") - c("rand+SDE")) + (H_ - c("rand+ASDE"))) / 3, "SVSS", "random")
    add(("A", "me2"), "Stage-2 main effect: ASDE − fixed SDE @ 317", ((c("rand+ASDE") - c("rand+SDE")) + (H_ - c("SVSS+SDE"))) / 2, "ASDE", "fixed SDE"); add(("A", "me3"), "Stage-2 main effect: ASDE − random fill @ 317", ((c("rand+ASDE") - Rd) + (H_ - c("SVSS+rand"))) / 2, "ASDE", "random fill")
    add(("A", "me4"), "Stage-2 main effect: fixed SDE − random fill @ 317", ((c("rand+SDE") - Rd) + (c("SVSS+SDE") - c("SVSS+rand"))) / 2, "fixed SDE", "random fill"); add(("A", "int"), "Interaction: (SVSS − random) under ASDE minus under fixed SDE @ 317", (H_ - c("rand+ASDE")) - (c("SVSS+SDE") - c("rand+SDE")), "interaction", "–")
    rows = []
    for k, t in tests.items():
        d = t["d"]; lo, hi = np.quantile([d[rng.integers(0, 32, 32)].mean() for _ in range(10000)], [.025, .975])
        rows.append(dict(key=str(k), family=k[0], comparison=t["comparison"], mean_diff=d.mean(), ci_low=lo, ci_high=hi, first_wins=int((d < 0).sum()), second_wins=int((d > 0).sum()), ties=int((d == 0).sum()), p_raw=wilcoxon(d).pvalue))
    D = pd.DataFrame(rows); D["p_holm_all"] = holm(D.p_raw.values); hm = D.family == "H"; D.loc[hm, "p_holm_headline_only"] = holm(D[hm].p_raw.values); D["sig_holm_all"] = D.p_holm_all < 0.05; D.to_csv(f"{OUT}/tests.csv", index=False)
    def ttab(fam):
        t = [f"| comparison (negative = first lower error) | mean diff | 95% CI | first wins / second wins (of 32) | Wilcoxon p | Holm p (all {len(D)} tests) | significant after Holm |", "|---|---|---|---|---|---|---|"]
        for _, r in D[D.family.isin(fam)].iterrows(): t.append(f"| {r.comparison} | {r.mean_diff:+.5f} | [{r.ci_low:+.5f}, {r.ci_high:+.5f}] | {r.first_wins} / {r.second_wins} | {r.p_raw:.2g} | {r.p_holm_all:.2g} | {'yes' if r.sig_holm_all else 'no'} |")
        return "\n".join(t)
    # thresholds
    th = ["| method | first budget below 0.09 | first budget below 0.085 |", "|---|---|---|"]
    for m in HEAD + EX:
        f_ = lambda t: next((str(b) for b in BUDGETS if g1.loc[(m, b), "mean"] < t), "--"); th.append(f"| {m} | {f_(0.09)} | {f_(0.085)} |")
    # best method per path among headline
    pw = pd.DataFrame({b: pd.Series({m: 0 for m in HEAD}) for b in BUDGETS})
    for b in BUDGETS:
        M = pd.DataFrame({m: P[m, b] for m in HEAD}); best = M.idxmin(axis=1).value_counts()
        for m in best.index: pw.loc[m, b] = best[m]
    pwt = ["| method | " + " | ".join(str(b) for b in BUDGETS) + " |", "|---|" + "---|" * 4] + [f"| {m} | " + " | ".join(str(int(pw.loc[m, b])) for b in BUDGETS) + " |" for m in HEAD]
    # noise floor
    nf = ["| method | " + " | ".join(f"share of headroom closed @ {b}" for b in BUDGETS) + " |", "|---|" + "---|" * 4]
    for m in HEAD + EX: nf.append(f"| {m} | " + " | ".join(f"{100 * (const['mean'] - g1.loc[(m, b), 'mean']) / head:.1f}%" for b in BUDGETS) + " |")
    hr = D[(D.family == "H") & D.comparison.str.contains("Random")]; hr_txt = ", ".join(f"{r.mean_diff / head * 100:+.1f}% @ {b}" for (_, r), b in zip(hr.iterrows(), BUDGETS))
    # ablation tables
    gA = g1.loc[[(m, T) for m in NEW2]].droplevel(1) if False else None
    A = lambda m: f"{g1.loc[(alias.get(m, m), T), 'mean']:.4f} ± {g1.loc[(alias.get(m, m), T), 'sd']:.4f}"
    fac = ["| stage 1 (95 sites) \\ stage 2 (222 sites) | random fill | fixed SDE (2,2) | adaptive ASDE |", "|---|---|---|---|", f"| random | {A('rand+rand')} | {A('rand+SDE')} | {A('rand+ASDE')} |", f"| SVSS | {A('SVSS+rand')} | {A('SVSS+SDE')} | {A('SVSS+ASDE')} (hybrid) |"]
    oth = ["| cell @ 317 sites | error (mean ± SD over seeds) |", "|---|---|"] + [f"| {m} | {A(m)} |" for m in ["SVSS-only", "AL-max (1% seed)", "ASDE-only(rand seed)", "ASDE-only(maximin seed)", "SDE ref"]]
    prm = ["ls0", "ls1", "var0", "var1", "noise"]; nmz = ["column lengthscale", "row lengthscale", "column variance", "row variance", "noise"]
    be = R[R.method.isin(HEAD + EX + NEW2 + ["gpr50"])].groupby("method")[[f"{q}_{e}" for q in prm for e in ("lo", "hi")]].mean()
    bt = ["| method | " + " | ".join(f"{x} lower / upper" for x in nmz) + " |", "|---|" + "---|" * 5] + [f"| {m} | " + " | ".join(f"{100 * be.loc[m, q + '_lo']:.1f}% / {100 * be.loc[m, q + '_hi']:.1f}%" for q in prm) + " |" for m in be.index]
    nt = R.groupby("budget").agg(n_train=("n_train", "mean"), n_test=("n_test", "mean"))
    md = f"""# Stage 4 – Final TEST run (FPGAs 27–50, run once)

## Protocol
- Frozen design: the current hybrid (SVSS 95 anchors, alpha = beta = 0.5, then ASDE with random fill), unchanged; chosen by the Stage 3 rule (no candidate qualified). SVTS weights (0.3 similarity, 0.7 variance) were chosen on TUNE in Stage 1. The two exploratory rows (anchors 1%, maximin anchors) were defined in `stage3_rule_section.md` before any result.
- Selection orders are the saved FPGA-1 orders of Stages 1–3 (5 seeds, same sites on all FPGAs); selection never used TEST data. Evaluation: additive bounded RBF GPR (same function as in selection), held-out 30% test set (about {R.n_test.mean():.0f} test sites per path; never selectable, never fitted), error = |mu − v| / d_spec, test sites → path → FPGA → mean over the 24 TEST FPGAs. Mean ± SD over 5 seeds of that FPGA-averaged error. FPGA 1 is excluded from the primary tables.
- Paired differences: path-level values (mean over the 5 seeds and 24 TEST FPGAs), 32 paths, 95% percentile bootstrap CI (10,000 resamples, seed 0), Wilcoxon signed-rank (two-sided, zero differences dropped) on the 32 path-level differences. Holm correction over all {len(D)} comparisons reported in this stage (24 headline, 8 exploratory, 14 ablation) as one family; the headline-only Holm p-values are in `stage4/tests.csv`.
- Budgets 95 / 158 / 222 / 317 sites including seed points. Noise-floor references: constant predictor = mean of a random 317-site training set; "GPR 50%" = the same GPR trained on a random {K50} sites (50% of the 3,173 sites, taken from the selectable 70% pool), both on the same held-out test sets. The GPR 50% reference uses one seed only (each fit takes about 20 s).
- Training sites actually available on TEST FPGAs: mean {nt.n_train.round(1).to_dict()} (the 3-sigma filter removes a few selected sites).

## Headline table (TEST 27–50, held-out test set, mean error ± SD over 5 seeds)
{tab(g1, HEAD)}

Exploratory rows (not part of the headline comparison):

{tab(g1, EX)}

References: constant predictor {const['mean']:.4f} ± {const['sd']:.4f}; GPR trained on 50% of the sites {g50['mean']:.4f} (seed 1 only, one fit of {K50} sites per path); headroom between them {head:.4f}.

## First budget below error 0.09 and 0.085
{chr(10).join(th)}

## Paired differences: hybrid against each baseline (headline)
{ttab(["H"])}

## Paired differences: hybrid against the exploratory designs
{ttab(["E"])}

## Paths where each method is the best of the seven headline methods (of 32 paths)
{chr(10).join(pwt)}

## Noise floor
Share of the headroom (constant − GPR 50%) closed by each method:

{chr(10).join(nf)}

The hybrid − random differences relative to the headroom: {hr_txt}.

## Ablation on TEST (Stage 2 repeated, 317 sites)
{chr(10).join(fac)}

{chr(10).join(oth)}

{ttab(["A"])}

## Secondary: TEST + FPGA 1 (headline and exploratory rows)
{tab(gF, HEAD + EX)}

## Secondary: paper's old test set (all non-training sites), TEST 27–50
{tab(gO, HEAD + EX)}

## Hyperparameters at a bound (share of predictor fits, within 5% of the bound)
{chr(10).join(bt)}
"""
    open(f"{HERE}/stage4_report.md", "w", encoding="utf8").write(md); print("report written", len(md))


if __name__ == "__main__":
    with Pool(8) as pool:
        pool.map(eval_task, [(f, p) for f in TEST for p in range(1, 33) if not os.path.exists(f"{OUT}/eval_f{f}_p{p}.csv")], chunksize=1)
    report()
