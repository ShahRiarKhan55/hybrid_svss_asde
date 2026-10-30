import json, random, time, types
from lib import *
p = 1; pool, held = heldout(p); X, Y = pool[["X", "Y"]].values.astype(float), pool.Value.values; n = len(pool)
print("pool", n, "held-out", len(held), "total", n + len(held))
S0 = start_set(n, 1, p); anchors = rand_order(n, 1, p)[:A]
hy.SDE_RATIO = (T - A) / n; old = [int(i) for i in hy.adaptive_sde(pool, anchors)]; new, acc = asde(X, Y, anchors, T - A, 42)
print("ASDE copy == hybrid_svss_asde.adaptive_sde (seed 42, 95 random anchors):", old == new, "| accepted by rule", acc)
sde = types.ModuleType("sde"); exec("".join(json.load(open(f"{ROOT}/sde.ipynb", encoding="utf8"))["cells"][0]["source"]), sde.__dict__)
print("sde_extend == sde.ipynb sde_select (seed 42):", sde.sde_select(X[:, 0], X[:, 1], T, random.Random(42)) == sde_extend(X, T, 42))
t = time.time(); m = fit(X[S0], Y[S0], "chk"); print("one bounded fit on 32 sites: %.2fs" % (time.time() - t), m.hits)
S = svss(X, Y, S0, A, "chk"); t = time.time(); m = fit(X[S], Y[S], "chk"); print("fit on 95 sites %.2fs" % (time.time() - t), {k: v for k, v in m.hits.items() if v})
print("budget sizes unique:", len(set(S)) == A, "| held-out sites in pool:", len(set(zip(pool.X, pool.Y)) & held))
print("maximin 95:", len(set(maximin(X, S0, A))), "| sde 317 unique:", len(set(sde_extend(X, T, 3))))
print(dict(STATS["chk"]))
