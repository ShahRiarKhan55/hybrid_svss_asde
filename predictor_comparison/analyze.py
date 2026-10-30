"""Build delta CSVs (all FPGAs and without FPGA 1), Wilcoxon tests and summary table from results/stage*/raw_fpga*.csv."""
import glob, os
import pandas as pd
from scipy.stats import wilcoxon

os.makedirs("results/final", exist_ok=True)
raw = pd.concat(pd.read_csv(f) for f in glob.glob("results/stage[12]/raw_fpga*.csv"))
raw.to_csv("results/final/raw_all.csv", index=False)
SETS = {"all": raw, "noFPGA1": raw[raw.FPGA != 1]}  # FPGA 1 is the calibration FPGA

summary = []
for tag, d in SETS.items():
    for (cond, m), g in d.groupby(["cond", "model"]):
        pre = f"results/final/{m}_{cond}_delta" + ("" if tag == "all" else "_noFPGA1")
        g[["FPGA", "PATH", "Delta"]].to_csv(f"{pre}_per_path.csv", index=False)  # same formats as the existing delta CSVs
        pf = g.groupby("FPGA").Delta.mean().reset_index()
        pf.to_csv(f"{pre}_per_fpga.csv", index=False)
        pd.DataFrame({"Overall Delta": [pf.Delta.mean()]}).to_csv(f"{pre}_overall.csv", index=False)
        if tag == "all":
            g[["FPGA", "PATH", "sec"]].to_csv(f"results/final/{m}_{cond}_timing.csv", index=False)
        summary.append((tag, cond, m, pf.Delta.mean(), g.Delta.std(), g.sec.mean(), len(g)))
s = pd.DataFrame(summary, columns=["set", "cond", "model", "mean_delta", "std_per_path", "sec_per_path", "n_paths"])
s.to_csv("results/final/summary.csv", index=False)

tests = []  # paired Wilcoxon over (FPGA, PATH) pairs; diff = first - second
for tag, d in SETS.items():
    w = d.pivot_table(index=["FPGA", "PATH"], columns=["cond", "model"], values="Delta")
    pairs = [((c, "gpr"), (c, m), f"gpr vs {m} [{c}]") for c in ("hybrid", "random") for m in ("svr", "rf", "mlp", "deepgp")]
    pairs += [(("hybrid", m), ("random", m), f"hybrid vs random [{m}]") for m in ("gpr", "svr", "rf", "mlp")]
    for a, b, name in pairs:
        if a in w and b in w:
            x = w[[a, b]].dropna()
            diff = x[a] - x[b]
            tests.append((tag, name, len(x), diff.mean(), diff.median(), wilcoxon(x[a], x[b]).pvalue))
pd.DataFrame(tests, columns=["set", "comparison", "n_pairs", "mean_diff", "median_diff", "p_value"]).to_csv("results/final/wilcoxon.csv", index=False)
print(s.round(4).to_string(index=False)); print(pd.read_csv("results/final/wilcoxon.csv").round(4).to_string(index=False))
