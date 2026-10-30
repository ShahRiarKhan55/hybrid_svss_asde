from lib import *
OUT = f"{HERE}/stage1"
G = pd.read_csv(f"{OUT}/results.csv"); C = pd.read_csv(f"{OUT}/cost_all.csv"); be = pd.read_csv(f"{OUT}/bounds_eval.csv", index_col=0); bs = pd.read_csv(f"{OUT}/bounds_selection.csv", index_col=0)
ch = pd.read_csv(f"{OUT}/svts_choice.csv").iloc[0]
NM = {"al": "AL (lowest variance)", "almax": "AL-max (highest variance)", "svts": "SVTS (0.3 sim, 0.7 var; chosen on TUNE)", "svts_other": "SVTS (0.7 sim, 0.3 var)", "sde": "SDE (fixed 5x5)", "random": "Random",
      "maximin": "Maximin (farthest point)", "hybrid": "Hybrid SVSS+ASDE (current)"}


def tab(fp, tn):
    g = G[(G.fpgas == fp) & (G.test == tn)].pivot(index="method", columns="budget", values=["mean", "sd"]); rows = ["| method | 95 (3%) | 158 (5%) | 222 (7%) | 317 (10%) |", "|---|---|---|---|---|"]
    for m in NM: rows.append(f"| {NM[m]} | " + " | ".join(f"{g['mean'][b][m]:.4f} ± {g['sd'][b][m]:.4f}" for b in BUDGETS) + " |")
    return "\n".join(rows)


a_ = G[(G.fpgas=="TUNE 2-26")&(G.test=="held-out 30%")].set_index(["method","budget"])["mean"]; b_ = G[(G.fpgas=="TUNE + FPGA 1")&(G.test=="held-out 30%")].set_index(["method","budget"])["mean"]; dmax = (b_-a_).abs().max(); same_order = all(list(a_.xs(b, level="budget").sort_values().index)==list(b_.xs(b, level="budget").sort_values().index) for b in BUDGETS)
h = C[C.method == "hybrid"]; acc = h.asde_accepted.mean()
cost = ["| method | GPR fits per path | seconds per path (8 workers in parallel) |", "|---|---|---|"]
for m, n in [("al", "AL"), ("almax", "AL-max"), ("svts_b", "SVTS (chosen)"), ("svts_a", "SVTS (0.7/0.3)"), ("hybrid", "Hybrid"), ("sde", "SDE"), ("maximin", "Maximin"), ("random", "Random")]:
    g = C[C.method == m]; cost.append(f"| {n} | {g.gpr_fits.mean():.1f} | {g.wall_s.mean():.2f} |")
P = ["ls0", "ls1", "var0", "var1", "noise"]; names = {"ls0": "column lengthscale", "ls1": "row lengthscale", "var0": "column variance", "var1": "row variance", "noise": "noise"}


def btab(df, ms):
    r = ["| method | " + " | ".join(f"{names[q]} lower / upper" for q in P) + " |", "|---|" + "---|" * len(P)]
    for m in ms: r.append(f"| {m} | " + " | ".join(f"{100 * df.loc[m, q + '_lo']:.1f}% / {100 * df.loc[m, q + '_hi']:.1f}%" for q in P) + " |")
    return "\n".join(r)


md = f"""# Stage 1 – Fair baseline table (TUNE FPGAs)

## Protocol as run (assumptions in bold)
- Selection once per path on FPGA 1, same coordinates on all FPGAs, selection uses only FPGA 1 values of already selected sites. Predictor and the GPR inside selection: one function (`lib.fit`): additive RBF (column + row), raw coordinates, normalizer, initial lengthscales/variances/noise 1.0, bounds lengthscale 0.5–300 cells, kernel variance 1e-4–1e3, noise 1e-4–10, no restarts, at most 200 optimiser iterations, a fresh model for every fit.
- Held-out test set: a fixed random 30% of the active sites of each path (taken on FPGA 1, 946 sites on path 1, 951.8 test sites per path on average over FPGAs), never selectable and never fitted. The selectable pool is the other 70%. Secondary: the paper's old test set (all non-training sites).
- Budgets 95 / 158 / 222 / 317 sites including seed points, fixed for every method and path. 5 seeds (1–5): the random 1% seed set (32 sites) shared by AL, AL-max, SVTS, SVSS, maximin; the SDE shuffle; the ASDE scan order and fill; the random baseline.
- AL, AL-max, SVTS: 32-site random seed, then one site per refit (285 fits). SVSS: batches of 32 (32 / 32 / 31), 2 fits. Maximin starts from the same 32-site seed set (**assumption**: Euclidean distance on raw coordinates).
- Error = mean |mu − v| / d_spec: test sites → path → FPGA → mean over FPGAs. Mean ± SD over the 5 seeds of that FPGA-averaged error. Primary evaluation excludes FPGA 1; tables with FPGA 1 included are below.
- **Deviation from the Stage 1 text: TEST FPGAs (27–50) are not run here.** The protocol says TEST is touched once, at the end (Stage 4). Stage 1 and 2 use TUNE (FPGAs 2–26, 25 devices) only.
- **SVTS weights were chosen on TUNE only.** Mean TUNE error over the four budgets: (0.7 similarity, 0.3 variance) {ch.svts_a_alpha0_7_beta0_3:.4f}; (0.3 similarity, 0.7 variance) {ch.svts_b_alpha0_3_beta0_7:.4f}. The second is used as "SVTS" below; both are shown.
- Checks done before running: the copied ASDE gives the same sites as `hybrid_svss_asde.ipynb`'s `adaptive_sde` (same seed, same anchors); `sde_extend` gives the same sites as `sde.ipynb`'s `sde_select`; no held-out site is in any training set (asserted in every evaluation).

## Main table – TUNE FPGAs 2–26, held-out 30% test set (mean error ± SD over 5 seeds)
{tab("TUNE 2-26", "held-out 30%")}

## Same, FPGA 1 included (26 devices)
{tab("TUNE + FPGA 1", "held-out 30%")}

## Secondary: paper's old test set (all non-training sites), TUNE 2–26
{tab("TUNE 2-26", "old (all non-selected)")}

## Selection cost (indicative only)
Measured while 8 selection tasks ran in parallel on a laptop that throttles; Stage 5 repeats it single-threaded. Hybrid: {acc:.1f} of the 222 ASDE additions ({100 * acc / 222:.1f}%) passed the exclusion rule, the rest is random fill.

{chr(10).join(cost)}

## Hyperparameters at a bound (share of GPR fits; lower = at the lower bound, upper = at the upper bound; within 5% of the bound)
Selection-loop fits (all fits of the method's selection run; hybrid has only 2 fits per path, on 32 and 64 sites):

{btab(bs, ["al", "almax", "svts_b", "svts_a", "hybrid"])}

Predictor fits in the evaluation (all budgets, all FPGAs incl. FPGA 1):

{btab(be, ["al", "almax", "svts_b", "svts_a", "sde", "random", "maximin", "hybrid"])}

## Reading the table (plain words)
- AL (lowest variance) is the worst method at every budget (0.0846 at 317 sites; random 0.0744).
- At 95 sites (the SVSS anchors alone) the hybrid is worse than random, SDE, maximin, AL-max and SVTS (0.0859 against 0.0815 for random).
- At 317 sites the order is random 0.0744, SDE 0.0746, maximin 0.0746, AL-max 0.0747, hybrid 0.0751, SVTS 0.0752, AL 0.0846. The spread between random and the hybrid is 0.0007, about two seed SDs (0.0003); the spread between the best six methods is 0.0008.
- The hybrid is not the lowest-error method at any budget in this table. Random, SDE and maximin are at or below the hybrid at all four budgets. Significance tests are deferred to Stage 4.
- The seed SD is 0.0002–0.0011; it is largest for SVTS at 95 sites and for AL.
- The column lengthscale reaches its lower bound (0.5 cells) in 22% of AL predictor fits and in 3–7% of the predictor fits of the other methods; in the selection loops it is 24% for AL and for the hybrid's two small fits. Upper bounds are almost never reached.
- With FPGA 1 included, every mean moves by at most {dmax:.4f}; the order of the methods is {"unchanged at all four budgets" if same_order else "not identical at every budget (see the two tables)"}.

## Differences from the old results
Different held-out test set (about 950 test sites per path; the old protocol tested on all non-selected sites, about 2,850), selection restricted to the 70% pool, fixed budgets for all methods, 1% random seed for AL / AL-max / SVTS (old: one site), bounded GPR, 5 seeds instead of one, TUNE devices instead of all 50.
"""
open(f"{HERE}/stage1_report.md", "w", encoding="utf8").write(md); print(md[:400])
