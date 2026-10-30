"""Deep GP, RBF in every layer: iso vs ARD, 50 inducing points, 200 iterations, hybrid sites only.
Usage: python dgpr.py   (parallel over FPGAs, resumable, then writes the CSVs and the comparison table)"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); PC = os.path.dirname(HERE)
sys.path.insert(0, PC)
import time
from multiprocessing import Pool
import GPy, deepgp, numpy as np, pandas as pd
from scipy.stats import wilcoxon
import run, stages

R = f"{HERE}/results"


def dgp_fit(tr, te, ard):  # same recipe as run.deepgp_fit (dims [1,3,2], 50 inducing, 200 iters); returns pred and training ELBO
    X = tr[["X", "Y"]].values.astype(float)
    mx, sx, my, sy = X.mean(0), X.std(0), tr.Value.mean(), tr.Value.std()
    m = deepgp.DeepGP([1, 3, 2], ((tr.Value.values - my) / sy)[:, None], X=(X - mx) / sx,
                      kernels=[GPy.kern.RBF(3, ARD=ard), GPy.kern.RBF(2, ARD=ard)], num_inducing=50, back_constraint=False)
    m.optimize(max_iters=200, messages=False)
    return m.predict((te[["X", "Y"]].values - mx) / sx)[0].ravel() * sy + my, m.log_likelihood()


def job(fpga):
    f = f"{R}/raw_fpga{fpga}.csv"
    if os.path.exists(f):
        return
    rows = []
    for path in range(1, 33):
        df = run.load_dat(fpga, path); dspec = df.Value.max() - df.Value.min()
        mk = stages.sites(fpga, path, df)["hybrid"]; tr, te = df[mk], df[~mk]
        for name, ard in (("dgp_iso", False), ("dgp_ard", True)):
            t = time.time(); p, elbo = dgp_fit(tr, te, ard)
            rows.append((fpga, path, name, np.mean(np.abs(p - te.Value.values)) / dspec, elbo, time.time() - t))
    pd.DataFrame(rows, columns=["FPGA", "PATH", "model", "Delta", "TrainLL", "sec"]).to_csv(f + ".tmp", index=False)
    os.replace(f + ".tmp", f); print(time.strftime("%H:%M:%S"), f"FPGA-{fpga} done", flush=True)


def build_table():
    ld = lambda pat, cond=None: pd.concat(pd.read_csv(f) for f in (f"{pat}{i}.csv" for i in range(1, 51)))
    dg = ld(f"{R}/raw_fpga"); un = ld(f"{PC}/untuned/results/raw_fpga")
    t1 = ld(f"{PC}/results/stage1/raw_fpga"); t1 = t1[t1.cond == "hybrid"]; t1.model += "_tuned"
    t2 = ld(f"{PC}/results/stage2/raw_fpga"); t2 = t2[t2.cond == "hybrid"]; t2.model = "dgp_iso_earlier"
    raw = pd.concat([dg, un, t1, t2]); raw.to_csv(f"{R}/raw_all.csv", index=False)
    w = raw.pivot(index=["FPGA", "PATH"], columns="model", values="Delta"); ll = raw.pivot(index=["FPGA", "PATH"], columns="model", values="TrainLL")
    fam = {"GPR": ("rbf_iso", "rbf_ard"), "DGP": ("dgp_iso", "dgp_ard")}
    sel = {}
    for k, (a, b) in fam.items():  # pick iso vs ARD per path from the training criterion only
        pick_ard = ll[b] > ll[a]
        sel[k] = (w[b].where(pick_ard, w[a]), f"ARD {pick_ard.sum()}/{len(w)}")
    rows = [("GPR  RBF iso (paper)", "gpr", "GPR"), ("GPR  RBF ARD", "rbf_ard", "GPR"),
            ("SVR  untuned", "svr", None), ("SVR  CV-tuned", "svr_tuned", None), ("RF   untuned", "rf", None), ("RF   CV-tuned", "rf_tuned", None),
            ("MLP  untuned", "mlp", None), ("MLP  CV-tuned", "mlp_tuned", None),
            ("Deep GP RBF iso", "dgp_iso", "DGP"), ("Deep GP RBF ARD", "dgp_ard", "DGP"), ("(Deep GP iso, earlier run)", "dgp_iso_earlier", None)]
    out = []
    for label, m, f in rows:
        d = w[m].rename("Delta").reset_index()
        out.append((label, d.groupby("FPGA").Delta.mean().mean(), w[m].std(), raw[raw.model == m].sec.mean(),
                    wilcoxon(w[m], w["gpr"]).pvalue if m != "gpr" else np.nan,
                    sel[f][0].groupby("FPGA").mean().mean() if f else np.nan, sel[f][1] if f else ""))
    t = pd.DataFrame(out, columns=["model", "mean_delta", "std_per_path", "sec_per_path", "p_vs_paper_GPR", "delta_by_train_criterion", "picked"])
    t.to_csv(f"{R}/table.csv", index=False); print(t.round(4).to_string(index=False))
    print("iso rerun vs earlier run, max |diff| per path:", (w.dgp_iso - w.dgp_iso_earlier).abs().max())


if __name__ == "__main__":
    os.makedirs(R, exist_ok=True)
    with Pool(16) as p:
        p.map(job, [f for f in range(1, 51) if not os.path.exists(f"{R}/raw_fpga{f}.csv")], chunksize=1)
    build_table()
