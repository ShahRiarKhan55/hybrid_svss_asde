## Success rule and candidates – written 2026-10-21 08:53 before any Stage 3 result exists

**Rule (as instructed).** A candidate is an improvement only if its 95% bootstrap CI of (candidate − random) lies entirely below 0 at 2 or more of the 4 budgets (95 / 158 / 222 / 317 sites) on the TUNE FPGAs (2–26).
- Bootstrap: path-level values (mean over the 5 seeds and the 25 TUNE FPGAs), 32 paths resampled, 10,000 resamples, seed 0, percentile CI [2.5%, 97.5%]; "below 0" means the upper CI bound is < 0. Primary held-out test set. No multiplicity correction.
- If no candidate qualifies, the current hybrid is frozen and this is stated plainly.
- Tie-break, fixed now: if more than one candidate qualifies, the one with the lowest TUNE error averaged over the four budgets is frozen.
- Candidate − current hybrid is reported for information only; it does not enter the rule.
- Stage 4 (TEST) is not started until you confirm.

**Candidates (trimmed list, nothing else changes; budgets fixed at 95 / 158 / 222 / 317 including anchors).** All use the same 32-site random seed set per seed and path, the current ASDE (K = 2, tau = 2.0, per-anchor windows, window of new sites = floor of the mean anchor window) and random fill unless stated.
1. **SVSS-cos:** SVSS anchors (95 sites, alpha = beta = 0.5, 2 fits) with the similarity term replaced by the cosine without mean-centring (definition in stage3_variant_defined.md), then ASDE.
2. **Anchor share (assumption: percent of the 3,173 sites).** (a) 1% = the 32 random seed sites only, no SVSS step, then ASDE to 317 (this is the Stage 2 "ASDE-only, random seed" design, evaluated here at all four budgets); (b) 2% = 63 anchors (32 random seed + one SVSS batch of 31, 1 fit), then ASDE to 317. The current hybrid is 3% = 95 anchors.
3. **Maximin anchors:** the 95 anchors are chosen by greedy farthest-point selection starting from the 32 random seed sites (Euclidean distance, raw coordinates), then ASDE.
4. **ASDE with window relaxation:** current SVSS anchors, current ASDE first pass; if fewer additions than needed pass the rule, all windows are reduced by 1 cell (not below 0) and the remaining candidates are rescanned in the same random order, repeating until the budget is met. No random fill.

Dropped as instructed: the K and tau grids and the alpha = beta grid.
