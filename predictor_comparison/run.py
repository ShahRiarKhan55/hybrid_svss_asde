"""Compare predictors on the same ~10% training sites (hybrid-selected vs random).
Usage: python run.py OUT_DIR FPGA_START FPGA_END
Reuses load_dat / gpr from hybrid_svss_asde.ipynb (cell 0) without modifying it."""
import json, os, sys, time, types, warnings
import numpy as np, pandas as pd
from sklearn.compose import TransformedTargetRegressor as TTR
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GridSearchCV
from sklearn.neural_network import MLPRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

# --- import existing code (notebook cell 0, __main__ guard keeps it from running)
hy = types.ModuleType("hy")
nb = json.load(open(f"{ROOT}/hybrid_svss_asde.ipynb", encoding="utf8"))
exec("".join(nb["cells"][0]["source"]), hy.__dict__)
load_dat, gpr, NUM_PATHS = hy.load_dat, hy.gpr, hy.NUM_PATHS


def hybrid_sites(path):  # saved FPGA1 training sites (selection is calibrated on FPGA1 only)
    d = pd.read_csv(f"{ROOT}/hybrid_svss_asde/FPGA1T/fpga1_path{path}_train_hybrid_svts_sde.csv")
    return set(zip(d.X, d.Y))


NJ = -1  # GridSearchCV workers; the parallel driver sets this to 1


# --- models: each takes (train_df, test_df) -> predictions
def sk(est, grid):
    def f(tr, te):
        m = GridSearchCV(TTR(make_pipeline(StandardScaler(), est), transformer=StandardScaler()),
                         {f"regressor__{k}": v for k, v in grid.items()}, cv=3, n_jobs=NJ)
        m.fit(tr[["X", "Y"]].values, tr.Value.values)
        return m.predict(te[["X", "Y"]].values)
    return f


def deepgp_fit(tr, te, ind=50, iters=200):
    import GPy, deepgp
    X = tr[["X", "Y"]].values.astype(float)
    mx, sx, my, sy = X.mean(0), X.std(0), tr.Value.mean(), tr.Value.std()
    m = deepgp.DeepGP([1, 3, 2], ((tr.Value.values - my) / sy)[:, None], X=(X - mx) / sx,
                      kernels=[GPy.kern.RBF(3), GPy.kern.RBF(2)], num_inducing=ind, back_constraint=False)
    m.optimize(max_iters=iters, messages=False)
    return m.predict((te[["X", "Y"]].values - mx) / sx)[0].ravel() * sy + my


MODELS = {
    "gpr": lambda tr, te: gpr(tr, te)[0],
    "svr": sk(SVR(), {"svr__C": [1, 10, 100], "svr__gamma": [0.1, 1, 10]}),
    "rf": sk(RandomForestRegressor(n_estimators=200, random_state=0),
             {"randomforestregressor__min_samples_leaf": [1, 3, 5], "randomforestregressor__max_features": [1.0, 0.5]}),
    "mlp": sk(MLPRegressor(max_iter=2000, random_state=0),
              {"mlpregressor__hidden_layer_sizes": [(64, 64), (128, 128)], "mlpregressor__alpha": [1e-4, 1e-2]}),
    "deepgp": deepgp_fit,
}

if __name__ == "__main__":
    out, f0, f1 = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    os.makedirs(out, exist_ok=True)
    rows = []
    for fpga in range(f0, f1 + 1):
        for path in range(1, NUM_PATHS + 1):
            df = load_dat(fpga, path)
            dspec = df.Value.max() - df.Value.min()
            pos = hybrid_sites(path)
            is_h = np.array([(x, y) in pos for x, y in zip(df.X, df.Y)])
            is_r = np.zeros(len(df), bool)  # random sites, same count as hybrid
            is_r[np.random.default_rng(fpga * 100 + path).choice(len(df), is_h.sum(), replace=False)] = True
            for cond, mask in [("hybrid", is_h), ("random", is_r)]:
                tr, te = df[mask], df[~mask]
                for name, fit in MODELS.items():
                    t = time.time()
                    pred = fit(tr, te)
                    rows.append((cond, name, fpga, path, np.mean(np.abs(pred - te.Value.values)) / dspec,
                                 time.time() - t))
            print(f"FPGA-{fpga} Path-{path} done", flush=True)
            pd.DataFrame(rows, columns=["cond", "model", "FPGA", "PATH", "Delta", "sec"]).to_csv(f"{out}/raw.csv", index=False)

    r = pd.read_csv(f"{out}/raw.csv")
    for (cond, name), g in r.groupby(["cond", "model"]):  # same formats as the existing delta CSVs
        pre = f"{out}/{name}_{cond}_delta"
        g[["FPGA", "PATH", "Delta"]].to_csv(f"{pre}_per_path.csv", index=False)
        pf = g.groupby("FPGA").Delta.mean().reset_index()
        pf.to_csv(f"{pre}_per_fpga.csv", index=False)
        pd.DataFrame({"Overall Delta": [pf.Delta.mean()]}).to_csv(f"{pre}_overall.csv", index=False)
    print(r.groupby(["cond", "model"]).agg(Delta=("Delta", "mean"), sec_per_path=("sec", "mean"), sec_total=("sec", "sum")))
