"""Items 1-5 from results/pred_*.csv and meta_*.csv. All errors are on the delta scale: |pred - value| / d_spec of the path."""
import glob, os
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import kurtosis, spearmanr, wilcoxon

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
M = ["gpr_iso", "gpr_ard", "gpr_add", "svr", "rf", "mlp"]
pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
meta = pd.concat(pd.read_csv(f) for f in glob.glob(f"{R}/meta_*.csv"))
D = []
for f in glob.glob(f"{R}/pred_f*_p*.csv"):
    b = os.path.basename(f)[6:-4].split("_p"); d = pd.read_csv(f).assign(FPGA=int(b[0]), PATH=int(b[1])); D.append(d)
D = pd.concat(D).merge(meta[["FPGA", "PATH", "dspec"]].drop_duplicates(), on=["FPGA", "PATH"])
for m in M: D["e_" + m] = (D["pred_" + m] - D.Value) / D.dspec  # signed residual, delta scale
E = ["e_" + m for m in M]; T = D[D.train == 0]; TR = D[D.train == 1]
print(f"{D.groupby(['FPGA', 'PATH']).ngroups} paths, {len(T)} test sites\n")
paper = lambda s: s.groupby(level=0).mean().mean()  # per-path -> per-FPGA mean -> mean, as in the paper
pp = lambda df, fn: df.groupby(["FPGA", "PATH"])[E].apply(fn).set_axis(M, axis=1)  # per-path statistic of each model


def save(name, t):
    t.to_csv(f"{R}/{name}.csv"); print(f"== {name}\n{t.round(4).to_string()}\n")


# 1. loss ----------------------------------------------------------------------------------------
mae, rmse = pp(T, lambda g: g.abs().mean()), pp(T, lambda g: np.sqrt((g ** 2).mean()))
i1 = pd.DataFrame({"delta_bar": mae.apply(paper), "rmse": rmse.apply(paper), "median_abs": T[E].abs().median().set_axis(M),
                   "p90_abs": T[E].abs().quantile(.9).set_axis(M), "p99_abs": T[E].abs().quantile(.99).set_axis(M),
                   "kurtosis": T[E].apply(kurtosis).set_axis(M)})
i1["rank_delta"], i1["rank_rmse"] = i1.delta_bar.rank().astype(int), i1.rmse.rank().astype(int)
save("item1_loss", i1); print("ranking changes between delta_bar and RMSE:", not (i1.rank_delta == i1.rank_rmse).all(), "\n")

# 2. structure -----------------------------------------------------------------------------------
prof = {}
for k in ("X", "Y"):  # MAE per column (X) and per row (Y), pooled over paths
    prof[k] = T.assign(**{m: T["e_" + m].abs() for m in M}).groupby(k)[M].mean(); prof[k].to_csv(f"{R}/profile_by_{k}.csv")
i2a = pd.DataFrame({f"{k}_{s}": getattr(prof[k], s)() for k in prof for s in ("std", "max", "min")}); i2a["Xcol_max/min"] = i2a.X_max / i2a.X_min; i2a["Yrow_max/min"] = i2a.Y_max / i2a.Y_min
save("item2a_mae_by_column_row_summary", i2a[["X_std", "X_max", "X_min", "Xcol_max/min", "Y_std", "Y_max", "Y_min", "Yrow_max/min"]])


def effects(tr, ev):  # column-only, row-only, additive (backfitted) on training sites; predictions at ev
    mu, y = tr.Value.mean(), tr.Value.values; a = pd.Series(0., index=tr.X.unique()); b = pd.Series(0., index=tr.Y.unique())
    for _ in range(30):
        a = pd.Series(y - mu - b.reindex(tr.Y).values, index=tr.index).groupby(tr.X).mean()
        b = pd.Series(y - mu - a.reindex(tr.X).values, index=tr.index).groupby(tr.Y).mean()
    cm, rm = tr.groupby("X").Value.mean(), tr.groupby("Y").Value.mean()
    return {"col_only": ev.X.map(cm).fillna(mu).values, "row_only": ev.Y.map(rm).fillna(mu).values,
            "additive_col+row": mu + a.reindex(ev.X).fillna(0).values + b.reindex(ev.Y).fillna(0).values}


rows = []
for (f, p), g in D.groupby(["FPGA", "PATH"]):
    tr, te = g[g.train == 1], g[g.train == 0]; ds = g.dspec.iloc[0]
    for n in effects(tr, tr):
        r2 = lambda d, pr: 1 - ((d.Value - pr) ** 2).sum() / ((d.Value - d.Value.mean()) ** 2).sum()
        rows.append((f, p, n, r2(tr, effects(tr, tr)[n]), r2(te, effects(tr, te)[n]), np.abs(te.Value - effects(tr, te)[n]).mean() / ds))
S = pd.DataFrame(rows, columns=["FPGA", "PATH", "model", "R2_train", "R2_test", "delta"])
i2b = S.groupby("model")[["R2_train", "R2_test"]].mean().assign(delta_bar=S.set_index(["FPGA", "PATH"]).groupby("model").delta.apply(paper))
i2b.loc["gpr_add (RBF)"] = [np.nan, np.nan, i1.delta_bar["gpr_add"]]; i2b.loc["gpr_ard"] = [np.nan, np.nan, i1.delta_bar["gpr_ard"]]
save("item2b_variance_explained", i2b)

# 3. edges and gaps ------------------------------------------------------------------------------
T = T.copy(); parts = []
for (f, p), g in D.groupby(["FPGA", "PATH"]):
    te = g[g.train == 0]; d = cKDTree(g[g.train == 1][["X", "Y"]].values).query(te[["X", "Y"]].values)[0]
    edge = np.minimum.reduce([te.X - g.X.min(), g.X.max() - te.X, te.Y - g.Y.min(), g.Y.max() - te.Y])
    parts.append(pd.DataFrame({"dist": d, "edge": np.asarray(edge)}, index=te.index.astype(str) + f"_{f}_{p}"))
T.index = T.index.astype(str) + "_" + T.FPGA.astype(str) + "_" + T.PATH.astype(str); T = T.join(pd.concat(parts))
T["bin"] = pd.cut(T.dist, [0, 1, 2, 3, 5, 10, 1e9], labels=["<=1", "(1,2]", "(2,3]", "(3,5]", "(5,10]", ">10"])
A = T[E].abs().set_axis(M, axis=1)
i3 = A.groupby(T.bin, observed=True).mean().assign(share=T.bin.value_counts(normalize=True).sort_index())
bd = (T.edge <= 2).map({True: "border(<=2)", False: "interior"})
i3 = pd.concat([i3, A.groupby(bd).mean().assign(share=bd.value_counts(normalize=True))]); save("item3_distance_border", i3)
print("median distance to nearest training site: border %.2f vs interior %.2f\n" % (T.dist[T.edge <= 2].median(), T.dist[T.edge > 2].median()))

# 4. paths ---------------------------------------------------------------------------------------
diff = mae.sub(mae.gpr_ard, axis=0).drop(columns="gpr_ard")
ard = meta[meta.model == "gpr_ard"].set_index(["FPGA", "PATH"]).sort_index()
cov = pd.DataFrame({"value_std": ard.vstd, "d_spec": ard.dspec, "std/d_spec": ard.vstd / ard.dspec, "noise_frac": ard.noise_frac,
                    "ls_x": ard.ls1, "ls_y": ard.ls2, "ls_y/ls_x": ard.ls2 / ard.ls1})
i4a = pd.DataFrame({"mean_diff": diff.mean(), "median_diff": diff.median(), "frac_paths_better": (diff < 0).mean(), "wilcoxon_p": diff.apply(lambda s: wilcoxon(s).pvalue)})
save("item4a_paired_vs_ard", i4a)
i4b = pd.DataFrame({m: {c: spearmanr(diff[m], cov.loc[diff.index, c])[0] for c in cov} for m in diff}).T; save("item4b_spearman_rho_diff_vs_covariate", i4b)
top = []
for m in diff:
    t = diff[m].nsmallest(10).reset_index(); t.columns = ["FPGA", "PATH", "diff"]; top.append(t.assign(model=m, rank=range(1, 11)))
t = (mae.gpr_ard - mae.drop(columns="gpr_ard").min(axis=1)).nsmallest(10).reset_index(); t.columns = ["FPGA", "PATH", "diff"]  # ARD vs best other model
top.append(t.assign(model="gpr_ard (vs best other)", rank=range(1, 11)))
top = pd.concat(top); top.to_csv(f"{R}/item4c_top10_wins.csv", index=False)
print("== item4c top-10 wins (diff = model - ARD, delta scale; negative = model better)")
for m, g in top.groupby("model", sort=False): print(f"{m:26s}", "  ".join(f"{a}/{b}:{c:+.3f}" for a, b, c in zip(g.FPGA, g.PATH, g["diff"])))

# 5. fit -----------------------------------------------------------------------------------------
tr_mae = pp(TR, lambda g: g.abs().mean())
i5 = pd.DataFrame({"in_sample": tr_mae.apply(paper), "test": mae.apply(paper)}); i5["ratio_test/in_sample"] = i5.test / i5.in_sample
i5["median_path_ratio"] = (mae / tr_mae).median(); print(); save("item5_fit", i5)
