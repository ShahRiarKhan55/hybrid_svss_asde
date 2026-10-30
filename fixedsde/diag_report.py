import json, numpy as np, pandas as pd
BAD = [2, 15, 16, 19, 23, 30, 31]; L = lambda k: {d["PATH"]: d for d in json.load(open(f"diag/run{k}.json"))}
A, B, C = L("A"), L("B"), L("C"); o = pd.read_csv("../hybrid_svss_asde/hybrid_svts_sde_delta_per_path.csv").set_index(["FPGA", "PATH"]).Delta
first = {(f, p): json.load(open(f"hybrid/f{f}_p{p}.json"))["Delta"] for f in (1, 2) for p in range(1, 33)}
print("== 1. per path (old, rerun=run A, diff)")
rows = []
for f in (1, 2):
    for p in BAD: rows.append((f, p, o[f, p], A[p][f"delta{f}"], A[p][f"delta{f}"] - o[f, p]))
print(pd.DataFrame(rows, columns=["FPGA", "PATH", "old", "rerun", "diff"]).to_string(index=False))
oldall = np.mean([o[f, p] for f in (1, 2) for p in range(1, 33)]); newall = np.mean([A[p][f"delta{f}"] for f in (1, 2) for p in range(1, 33)])
oldr = np.mean([o[f, p] for f in (1, 2) for p in range(1, 33) if p not in BAD]); print(f"mean of 2 FPGAs (64 values): old {oldall:.6f} rerun {newall:.6f} diff {newall-oldall:+.6f}; change due to 7 paths = all of it (other 25 paths max diff {max(abs(A[p][f'delta{f}']-o[f,p]) for f in (1,2) for p in range(1,33) if p not in BAD):.1e})")
print("per-FPGA mean old/new:", [(f, round(np.mean([o[f, p] for p in range(1, 33)]), 6), round(np.mean([A[p][f'delta{f}'] for p in range(1, 33)]), 6)) for f in (1, 2)])
print("== 2. run-to-run")
for nm, X in (("A vs B", B), ("A vs C(7 paths)", C), ("B vs C", None)):
    Y = A if nm.startswith("A") else B; X = B if nm == "A vs B" else C
    ps = [p for p in X]; d = max(abs(Y[p][f"delta{f}"] - X[p][f"delta{f}"]) for p in ps for f in (1, 2)); same = all(Y[p]["pos"] == X[p]["pos"] for p in ps); print(nm, "max |delta diff|", d, "same selected sites on all paths:", same)
print("first run (cache run.py) vs A max diff", max(abs(first[f, p] - A[p][f"delta{f}"]) for f in (1, 2) for p in range(1, 33)))
print("A,B vs old: paths with |diff|>1e-9:", sorted({p for p in range(1, 33) for f in (1, 2) if abs(A[p][f"delta{f}"] - o[f, p]) > 1e-9}), "/", sorted({p for p in range(1, 33) for f in (1, 2) if abs(B[p][f"delta{f}"] - o[f, p]) > 1e-9}))
print("== 3. diagnostics, bad paths vs others")
tr1 = {p: pd.read_csv(f"../hybrid_svss_asde/FPGA1T/fpga1_path{p}_train_hybrid_svts_sde.csv") for p in range(1, 33)}
rows = []
for p in range(1, 33):
    i = A[p]["info"]; old = {(int(a), int(b)) for a, b in zip(tr1[p].X, tr1[p].Y)}; mine = set(map(tuple, A[p]["pos"])); anc = set(map(tuple, A[p]["anchors"]))
    fits = i["svss_fits"] + i["final_fits"]; f = lambda x: x
    rows.append(dict(path=p, bad=p in BAD, fallback=A[p]["fallback"], n_sel=A[p]["n_sel"], anchors_not_in_old=len(anc - old), sites_not_in_old=len(mine - old), old_not_in_mine=len(old - mine), min_boundary_gap=min(i["boundary_gap"]), exact_ties=i["exact_ties_at_boundary"], dup_scores=i["n_dup_scores"],
        not_converged=sum(("onverged" not in x["status"]) for x in fits), maxiter=sum("maxiter" in x["status"].lower() for x in fits), ls_min=min(x["ls"] for x in fits), ls_max=max(x["ls"] for x in fits), noise_min=min(x["noise"] for x in fits), var_max=max(x["var"] for x in fits)))
D = pd.DataFrame(rows); pd.set_option("display.width", 250); print(D[D.bad].to_string(index=False)); print("-- others (summary)"); print(D[~D.bad].describe().loc[["min", "50%", "max"]].to_string())
print(D.groupby("bad")[["fallback", "anchors_not_in_old", "sites_not_in_old", "not_converged", "maxiter"]].sum())
print("statuses:", pd.Series([x["status"] for p in A for x in A[p]["info"]["svss_fits"] + A[p]["info"]["final_fits"]]).value_counts().to_dict())
