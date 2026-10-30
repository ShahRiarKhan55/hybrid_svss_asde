import os
for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[v] = "1"
import sys, random, warnings
from collections import Counter, defaultdict
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); PC = f"{ROOT}/predictor_comparison"
sys.path.insert(0, PC)
import numpy as np, pandas as pd, GPy
from scipy.spatial.distance import cdist
warnings.filterwarnings("ignore")
import run  # predictor_comparison/run.py: hy.load_dat (3-sigma filter)
hy = run.hy

TUNE, TEST = list(range(2, 27)), list(range(27, 51))
BUDGETS = [95, 158, 222, 317]; T, A, SEEDS = 317, 95, [1, 2, 3, 4, 5]; STEP = 32
LS, VAR, NOISE = (0.5, 300.0), (1e-4, 1e3), (1e-4, 10.0)
STATS = defaultdict(Counter)


def fit(X, y, tag="x"):
    m = GPy.models.GPRegression(X, y.reshape(-1, 1), GPy.kern.RBF(1, active_dims=[0]) + GPy.kern.RBF(1, active_dims=[1]), normalizer=True)
    for k in m.kern.parts:
        k.lengthscale.constrain_bounded(*LS, warning=False); k.variance.constrain_bounded(*VAR, warning=False)
    m.Gaussian_noise.variance.constrain_bounded(*NOISE, warning=False)
    m.optimize(max_iters=200, messages=False)
    p, c = m.kern.parts, STATS[tag]; c["fits"] += 1; m.hits = {}
    for nm, val, (lo, hi) in [("ls0", p[0].lengthscale[0], LS), ("ls1", p[1].lengthscale[0], LS), ("var0", p[0].variance[0], VAR), ("var1", p[1].variance[0], VAR), ("noise", m.Gaussian_noise.variance[0], NOISE)]:
        m.hits[nm + "_lo"], m.hits[nm + "_hi"] = int(val <= lo * 1.05), int(val >= hi / 1.05)
        c[nm + "_lo"] += m.hits[nm + "_lo"]; c[nm + "_hi"] += m.hits[nm + "_hi"]
    return m


def heldout(p):  # fixed random 30% of the active sites of the calibration FPGA, never selectable, same for every method, budget and seed
    df = hy.load_dat(1, p); n = len(df); h = np.zeros(n, bool); h[np.random.default_rng([12345, p]).choice(n, round(0.3 * n), replace=False)] = True
    return df[~h].reset_index(drop=True), set(zip(df.X[h], df.Y[h]))


def start_set(n, seed, p): return [int(i) for i in np.random.default_rng([seed, p, 1]).choice(n, STEP, replace=False)]


def rand_order(n, seed, p): return [int(i) for i in np.random.default_rng([seed, p, 77]).permutation(n)[:T]]


def score(m, X, tr, rem, a, b, center=True):
    mr, vr = m.predict(X[rem]); ms, vs = m.predict(X[tr])
    am, av = np.concatenate([mr.ravel(), ms.ravel()]), np.concatenate([vr.ravel(), vs.ravel()])
    mz, vz = (am - am.mean()) / (am.std() + 1e-9), (av - av.mean()) / (av.std() + 1e-9); k = len(mr)
    fr, fs = np.column_stack([mz[:k], vz[:k]]), np.column_stack([mz[k:], vz[k:]]); c = fs.mean(0)
    if center: cc, f = c - c.mean(), fr - fr.mean(1, keepdims=True)
    else: cc, f = c, fr
    sim = (f @ cc) / (np.linalg.norm(f, axis=1) * np.linalg.norm(cc) + 1e-9)
    return a * (sim - sim.min()) / (sim.max() - sim.min() + 1e-9) + b * (vr.ravel() - vr.min()) / (vr.max() - vr.min() + 1e-9)


def al(X, Y, S0, mode, tag):  # random 1% seed, then one site per refit: lowest (AL) or highest (AL-max) predictive variance
    tr = list(S0)
    while len(tr) < T:
        rem = np.setdiff1d(np.arange(len(X)), tr); v = fit(X[tr], Y[tr], tag).predict(X[rem])[1].ravel(); tr.append(int(rem[np.argmin(v) if mode == "min" else np.argmax(v)]))
    return tr


def svts(X, Y, S0, a, b, tag, center=True):  # random 1% seed, then one site per refit: argmax of alpha*sim + beta*var
    tr = list(S0)
    while len(tr) < T:
        rem = np.setdiff1d(np.arange(len(X)), tr); tr.append(int(rem[np.argmax(score(fit(X[tr], Y[tr], tag), X, tr, rem, a, b, center))]))
    return tr


def svss(X, Y, S0, target, tag, a=0.5, b=0.5, center=True):  # batches of 1% (32), refit per batch, top-score sites
    tr = list(S0)
    while len(tr) < target:
        rem = np.setdiff1d(np.arange(len(X)), tr); s = score(fit(X[tr], Y[tr], tag), X, tr, rem, a, b, center)
        tr += [int(rem[i]) for i in np.argsort(s)[::-1][:min(STEP, target - len(tr))]]
    return tr


LAST = {}


def asde(X, Y, anchors, target, seed, K=2, tau=2.0, window="floor", relax=False):  # the existing ASDE (hybrid_svss_asde.ipynb), seed-parametrised; returns additions and how many passed the rule
    an = np.array(anchors); D = cdist(X[an], X[an], "sqeuclidean"); np.fill_diagonal(D, np.inf); ab = []
    for r, j in enumerate(an):
        nn = an[np.argsort(D[r], kind="stable")[:K]]
        bump = int(np.std(Y[nn]) > tau); ab.append((int(np.ceil(np.abs(X[nn, 0] - X[j, 0]).mean())) + bump, int(np.ceil(np.abs(X[nn, 1] - X[j, 1]).mean())) + bump))
    ab = np.array(ab, float); dm = ab.mean(0); d = {"floor": np.floor(dm), "round": np.round(dm), "ceil": np.ceil(dm)}[window]
    sx, sy, sa, sb = list(X[an, 0]), list(X[an, 1]), list(ab[:, 0]), list(ab[:, 1]); add, aset = [], set(anchors)
    perm = np.random.RandomState(seed).permutation(len(X))
    for i in perm:
        if len(add) >= target: break
        if i in aset: continue
        if np.all((np.abs(np.array(sx) - X[i, 0]) > np.array(sa)) | (np.abs(np.array(sy) - X[i, 1]) > np.array(sb))):
            add.append(int(i)); sx.append(X[i, 0]); sy.append(X[i, 1]); sa.append(d[0]); sb.append(d[1])
    n_acc, got = len(add), set(add)
    if relax:  # shrink every window by 1 cell per pass (not below 0) and rescan in the same order until the budget is met; no random fill
        r = 0
        while len(add) < target:
            r += 1; start = len(add); got = set(add)
            ax, ay = list(sx), list(sy); wa, wb = [max(a - r, 0) for a in sa], [max(b - r, 0) for b in sb]
            for i in perm:
                if len(add) >= target: break
                if i in aset or int(i) in got: continue
                if np.all((np.abs(np.array(ax) - X[i, 0]) > np.array(wa)) | (np.abs(np.array(ay) - X[i, 1]) > np.array(wb))):
                    add.append(int(i)); got.add(int(i)); ax.append(X[i, 0]); ay.append(X[i, 1]); wa.append(max(d[0] - r, 0)); wb.append(max(d[1] - r, 0))
            for i in add[start:]: sx.append(X[i, 0]); sy.append(X[i, 1]); sa.append(d[0]); sb.append(d[1])
        LAST["relax_level"] = r
    elif n_acc < target: add += [int(i) for i in perm if i not in aset and int(i) not in got][:target - n_acc]
    return add, n_acc


def sde_extend(X, target, seed, start=None, a=2, b=2):  # the existing SDE rule (fixed window), optionally continuing from given anchors; start=None = sde.ipynb's sde_select
    order = list(range(len(X))); random.Random(seed).shuffle(order); sel = list(start) if start is not None else [order[0]]; ss, left = set(sel), []
    sx, sy = list(X[sel, 0]), list(X[sel, 1])
    for i in order:
        if len(sel) >= target: break
        if i in ss: continue
        if np.all((np.abs(np.array(sx) - X[i, 0]) > a) | (np.abs(np.array(sy) - X[i, 1]) > b)): sel.append(i); ss.add(i); sx.append(X[i, 0]); sy.append(X[i, 1])
        else: left.append(i)
    return sel + left[:target - len(sel)] if len(sel) < target else sel


def maximin(X, S0, target):  # greedy farthest-point from the same random 1% seed
    sel = list(S0); d = cdist(X, X[sel]).min(1)
    while len(sel) < target:
        i = int(np.argmax(d)); sel.append(i); d = np.minimum(d, np.hypot(X[:, 0] - X[i, 0], X[:, 1] - X[i, 1]))
    return sel


def evaluate(f, p, orders, budgets, tag):
    """orders {(method, seed): [(x, y), ...]}. Fit on the first b sites of each order (those present on FPGA f), predict all other sites of FPGA f.
    primary error = held-out 30% sites; secondary = all non-training sites (the paper's old test set). Error = |mu - v| / d_spec."""
    df = hy.load_dat(f, p); ds = df.Value.max() - df.Value.min(); co = list(zip(df.X, df.Y)); cid = {c: i for i, c in enumerate(co)}
    held = heldout(p)[1]; is_t = np.array([c in held for c in co]); X, v = df[["X", "Y"]].values.astype(float), df.Value.values; rows = []
    for (m, s), order in orders.items():
        idx = np.array([cid.get(c, -1) for c in order])
        for b in budgets:
            ix = idx[:b]; ix = ix[ix >= 0]; assert not is_t[ix].any(), "held-out site in a training set"
            mk = np.zeros(len(df), bool); mk[ix] = True; md = fit(X[mk], v[mk], tag)
            e = np.abs(md.predict(X[~mk])[0].ravel() - v[~mk]) / ds; tt = is_t[~mk]
            rows.append(dict(FPGA=f, PATH=p, method=m, seed=s, budget=b, n_train=int(mk.sum()), n_test=int(tt.sum()), err=e[tt].mean(), err_old=e.mean(), **md.hits))
    return rows


def agg(R, fpgas, col="err"):  # test sites -> path (already) -> FPGA mean -> mean over FPGAs, per (method, budget, seed)
    Q = R[R.FPGA.isin(fpgas)]; g = Q.groupby(["method", "budget", "seed", "FPGA"])[col].mean().groupby(["method", "budget", "seed"]).mean()
    return g.groupby(["method", "budget"]).agg(mean="mean", sd="std")
