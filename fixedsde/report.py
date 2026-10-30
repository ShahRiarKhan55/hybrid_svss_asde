import json, numpy as np, pandas as pd
from scipy.stats import binomtest, wilcoxon
import run as r
R = {}
for nm, pre in (("hybrid", "hybrid_rerun"), ("fixed", "hybrid_svss_fixedsde")):
    d = r.load(nm, range(1, 51)); R[nm] = d
    d[["FPGA", "PATH", "Delta"]].to_csv(f"{pre}_delta_per_path.csv", index=False)
    pf = d.groupby("FPGA").Delta.mean().reset_index(); pf.to_csv(f"{pre}_delta_per_fpga.csv", index=False)
    pd.DataFrame({"Overall Delta": [pf.Delta.mean()]}).to_csv(f"{pre}_delta_overall.csv", index=False)
a, b = R["hybrid"], R["fixed"]; fa, fb = [x.groupby("FPGA").Delta.mean() for x in (a, b)]; pa, pb = [x.groupby("PATH").Delta.mean() for x in (a, b)]
print(f"overall a (adaptive, rerun) {fa.mean():.6f}  b (fixed 2x2) {fb.mean():.6f}  b-a {fb.mean()-fa.mean():+.6f}  old 0.079824 -> a-old {fa.mean()-0.079824498:+.6f}")
na, nb = int((fa < fb).sum()), int((fb < fa).sum()); print(f"FPGAs lower: a {na}, b {nb}, ties {50-na-nb}; sign test p = {binomtest(na, na+nb, .5).pvalue:.3g}")
d = pb - pa; print(f"path-level (32): mean b-a {d.mean():+.6f}, b lower on {(d<0).sum()}, a lower on {(d>0).sum()}; Wilcoxon p = {wilcoxon(d).pvalue:.3g}")
S = {nm: [json.load(open(f"sel/p{p}.json"))[nm] for p in range(1, 33)] for nm in ("hybrid", "fixed")}
for nm in S:
    fb_ = np.array([s["fallback"] for s in S[nm]]); print(f"{nm}: calibration set size mean {np.mean([s['n'] for s in S[nm]]):.2f} (min {min(s['n'] for s in S[nm])}); mean train sites actually present per FPGA-path {R[nm].n_train.mean():.2f} (min {R[nm].n_train.min()}, share==317 {(R[nm].n_train==317).mean():.3f}); fallback fired on {int((fb_>0).sum())}/32 paths, mean {fb_.mean():.1f} of 222 sites ({fb_.mean()/222:.1%}), total {fb_.sum()}, accepted by window rule {222-fb_.mean():.1f}")
print("per-FPGA delta a-b by FPGA (first 5):", (fa-fb).head().round(5).to_dict())
