"""SVSS 3% + ASDE 7% (hybrid) vs SVSS 3% + fixed (2,2) SDE 7%. Reuses cell 0 of hybrid_svss_asde.ipynb unchanged.
python run.py check | run | report"""
import os, sys, json, numpy as np, pandas as pd
from multiprocessing import Pool
R = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); H = os.path.dirname(os.path.abspath(__file__))
hy = {"__name__": "hy"}; exec("".join(json.load(open(f"{R}/hybrid_svss_asde.ipynb", encoding="utf8"))["cells"][0]["source"]), hy)
SEED = hy["SEED"]; NP = 32


def sde(df, train, fixed):
    """copy of hy.adaptive_sde; fixed=(a,b) replaces per-anchor windows. Returns (picked, n_fallback)."""
    target = max(1, round(len(df) * hy["SDE_RATIO"])); ab = {}
    for j in train:
        xj, yj = df.loc[j, "X"], df.loc[j, "Y"]
        nn = sorted([p for p in train if p != j], key=lambda p: (df.loc[p, "X"] - xj) ** 2 + (df.loc[p, "Y"] - yj) ** 2)[:hy["K_NEIGHBORS"]]
        dx = np.mean([abs(df.loc[p, "X"] - xj) for p in nn]); dy = np.mean([abs(df.loc[p, "Y"] - yj) for p in nn])
        bump = 1 if np.std([df.loc[p, "Value"] for p in nn]) > hy["ROUGHNESS_THRESH"] else 0
        ab[j] = (int(np.ceil(dx)) + bump, int(np.ceil(dy)) + bump)
    d = fixed or (int(np.mean([v[0] for v in ab.values()])), int(np.mean([v[1] for v in ab.values()])))
    if fixed: ab = {j: fixed for j in train}
    sel, out = list(train), []
    for i in np.random.RandomState(SEED).permutation(df.index):
        if i in train or len(out) >= target: continue
        x, y = df.loc[i, "X"], df.loc[i, "Y"]
        if all(abs(df.loc[j, "X"] - x) > ab.get(j, d)[0] or abs(df.loc[j, "Y"] - y) > ab.get(j, d)[1] for j in sel):
            out.append(i); sel.append(i); ab[i] = d
    nf = target - len(out)
    if nf > 0:
        out += [i for i in np.random.RandomState(SEED).permutation(df.index) if i not in train and i not in out][:nf]
    return out, nf


def select(p):
    f = f"{H}/sel/p{p}.json"
    if os.path.exists(f): return json.load(open(f))
    df = hy["load_dat"](1, p); al = hy["svts_select"](df); r = {}
    for name, fx in (("hybrid", None), ("fixed", (2, 2))):
        s, nf = sde(df, al, fx); r[name] = dict(pos=[[int(a), int(b)] for a, b in df.loc[al + s, ["X", "Y"]].values], fallback=nf, n=len(al + s))
        if name == "hybrid": assert sorted(s) == sorted(hy["adaptive_sde"](df, al)), "replica differs from original adaptive_sde"
    r["n_svss"] = len(al); os.makedirs(f"{H}/sel", exist_ok=True); json.dump(r, open(f, "w")); return r


def job(a):
    name, fpga, p = a; out = f"{H}/{name}/f{fpga}_p{p}.json"
    if os.path.exists(out): return
    df = hy["load_dat"](fpga, p); ds = df.Value.max() - df.Value.min(); pos = {tuple(x) for x in select(p)[name]["pos"]}
    tr = [i for i in df.index if (df.X[i], df.Y[i]) in pos]; te = [i for i in df.index if i not in set(tr)]
    mu, _ = hy["gpr"](df.loc[tr], df.loc[te]); d = float(np.mean(np.abs(mu - df.Value[te].values) / ds))
    os.makedirs(f"{H}/{name}", exist_ok=True); json.dump(dict(FPGA=fpga, PATH=p, Delta=d, n_train=len(tr)), open(out, "w"))


def load(name, fp):
    return pd.DataFrame([json.load(open(f"{H}/{name}/f{f}_p{p}.json")) for f in fp for p in range(1, NP + 1)])


if __name__ == "__main__":
    mode = sys.argv[1]
    with Pool(8) as pool:
        pool.map(select, range(1, NP + 1))
        if mode == "check":
            pool.map(job, [("hybrid", f, p) for f in (1, 2) for p in range(1, NP + 1)])
            o = pd.read_csv(f"{R}/hybrid_svss_asde/hybrid_svts_sde_delta_per_path.csv"); n = load("hybrid", (1, 2)).merge(o, on=["FPGA", "PATH"], suffixes=("_new", "_old"))
            n["abs_diff"] = (n.Delta_new - n.Delta_old).abs(); print(n.abs_diff.describe()); print("max abs diff", n.abs_diff.max(), "| mean new", n.Delta_new.mean(), "old", n.Delta_old.mean())
        elif mode == "run":
            pool.map(job, [(nm, f, p) for nm in ("hybrid", "fixed") for f in range(1, 51) for p in range(1, NP + 1)], chunksize=1)
