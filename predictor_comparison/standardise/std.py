"""Exact RBF GPR, isotropic/ARD x raw/standardised coordinates, hybrid sites, 50 FPGAs x 32 paths.
Same recipe as the paper's gpr(): GPRegression(normalizer=True), optimize(max_iters=200). Usage: python std.py (parallel, resumable)"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); PC = os.path.dirname(HERE)
sys.path.insert(0, PC)
import time
from multiprocessing import Pool
import GPy, numpy as np, pandas as pd
from scipy.stats import wilcoxon
import run, stages

R = f"{HERE}/results"
FITS = {"iso_raw": (False, False), "ard_raw": (True, False), "iso_std": (False, True), "ard_std": (True, True)}  # (ARD, standardise)
LO, HI = 0.01, 100.0  # "extreme" = lengthscale below 0.01 or above 100 training-std units of its coordinate


def fit(ard, std, tr, te):
    X, Xt = tr[["X", "Y"]].values, te[["X", "Y"]].values
    mx, sx = (X.mean(0), X.std(0)) if std else (0, 1)
    sc = X.std(0)  # only used to express lengthscales in training-std units
    m = GPy.models.GPRegression((X - mx) / sx, tr[["Value"]].values, GPy.kern.RBF(2, ARD=ard), normalizer=True)
    m.optimize(max_iters=200, messages=False)
    ls = np.broadcast_to(m.kern.lengthscale.values, 2) * sx / sc  # in training-std units
    ok = m.optimization_runs[-1].status == "Converged"
    return m.predict((Xt - mx) / sx)[0].ravel(), ls, ok, float(m.log_likelihood()), float(m.Gaussian_noise.variance[0])


def job(f):
    out = f"{R}/raw_fpga{f}.csv"
    if os.path.exists(out):
        return
    rows = []
    for p in range(1, 33):
        df = run.load_dat(f, p); dspec = df.Value.max() - df.Value.min()
        mk = stages.sites(f, p, df)["hybrid"]; tr, te = df[mk], df[~mk]
        for n, (ard, std) in FITS.items():
            np.random.seed(f * 100 + p); t = time.time()
            pr, ls, ok, ll, noise = fit(ard, std, tr, te)
            rows.append((f, p, n, np.mean(np.abs(pr - te.Value.values)) / dspec, ls[0], ls[1], ok, ll, noise, time.time() - t))
    pd.DataFrame(rows, columns=["FPGA", "PATH", "model", "Delta", "ls_x_std", "ls_y_std", "converged", "LL", "noise", "sec"]).to_csv(out + ".tmp", index=False)
    os.replace(out + ".tmp", out); print(time.strftime("%H:%M:%S"), f"FPGA-{f} done", flush=True)


if __name__ == "__main__":
    os.makedirs(R, exist_ok=True)
    with Pool(16) as pool:
        pool.map(job, [f for f in range(1, 51) if not os.path.exists(f"{R}/raw_fpga{f}.csv")], chunksize=1)
    r = pd.concat(pd.read_csv(f"{R}/raw_fpga{f}.csv") for f in range(1, 51)); r.to_csv(f"{R}/raw_all.csv", index=False)
    for m, g in r.groupby("model"):  # existing delta-CSV formats
        g[["FPGA", "PATH", "Delta"]].to_csv(f"{R}/{m}_hybrid_delta_per_path.csv", index=False)
        pf = g.groupby("FPGA").Delta.mean().reset_index(); pf.to_csv(f"{R}/{m}_hybrid_delta_per_fpga.csv", index=False)
        pd.DataFrame({"Overall Delta": [pf.Delta.mean()]}).to_csv(f"{R}/{m}_hybrid_delta_overall.csv", index=False)
    ref = pd.read_csv(f"{PC}/../hybrid_svss_asde/hybrid_svts_sde_delta_per_path.csv").set_index(["FPGA", "PATH"]).Delta  # the paper's GPR
    r["ext"] = ((r[["ls_x_std", "ls_y_std"]] < LO) | (r[["ls_x_std", "ls_y_std"]] > HI)).any(axis=1)
    r["ext_lo"] = (r[["ls_x_std", "ls_y_std"]] < LO).any(axis=1); r["ext_hi"] = (r[["ls_x_std", "ls_y_std"]] > HI).any(axis=1)
    out = []
    for m in FITS:
        g = r[r.model == m].set_index(["FPGA", "PATH"]); d = g.Delta
        out.append((m, g.groupby("FPGA").Delta.mean().mean(), d.std(), (d - ref.loc[d.index]).abs().max(),
                    wilcoxon(d, ref.loc[d.index]).pvalue if m != "iso_raw" else np.nan,
                    g.ext.mean(), g.ext_lo.mean(), g.ext_hi.mean(), 1 - g.converged.mean(), g.LL.mean(), g.sec.mean()))
    t = pd.DataFrame(out, columns=["fit", "mean_delta", "std", "max|diff vs paper|", "p_vs_paper_GPR", "frac_extreme", "frac_ls<0.01", "frac_ls>100", "frac_not_converged", "mean_train_LL", "sec"])
    t.to_csv(f"{R}/table.csv", index=False); print(t.round(4).to_string(index=False))
