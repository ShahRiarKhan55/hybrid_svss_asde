# Stage 3 – Improving the hybrid (TUNE FPGAs only)

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

## Protocol
Same as Stages 1–2: held-out 30% test set, additive bounded RBF GPR (same function inside selection and for prediction), 5 seeds, selection on FPGA 1, TUNE FPGAs 2–26 (25 devices), FPGA 1 only in a secondary table, error = |mu − v| / d_spec averaged test sites → path → FPGA. TEST FPGAs (27–50) were not used. Random, SDE and the current hybrid are the Stage 1 results (same seeds, same test sets).

## Verdict (applying the rule written above)
| candidate | budgets with CI below 0 (of 4) | mean TUNE error over budgets | improvement under the rule |
|---|---|---|---|
| ASDE relax | 0 | 0.0791 | no |
| Hybrid (current) | 0 | 0.0791 | no |
| SVSS-cos | 0 | 0.0789 | no |
| anchors 1% | 0 | 0.0776 | no |
| anchors 2% | 0 | 0.0778 | no |
| maximin anchors | 0 | 0.0777 | no |

**Frozen design: Hybrid (current).** No candidate met the rule, so the current hybrid is frozen unchanged.

## Error by budget (TUNE 2–26, held-out test set; mean ± SD over 5 seeds)
| method | 95 | 158 | 222 | 317 |
|---|---|---|---|---|
| Random | 0.0815 ± 0.0007 | 0.0780 ± 0.0005 | 0.0761 ± 0.0005 | 0.0744 ± 0.0003 |
| SDE | 0.0812 ± 0.0004 | 0.0774 ± 0.0003 | 0.0759 ± 0.0003 | 0.0746 ± 0.0002 |
| Hybrid (current) | 0.0859 ± 0.0004 | 0.0790 ± 0.0007 | 0.0766 ± 0.0003 | 0.0751 ± 0.0003 |
| SVSS-cos | 0.0855 ± 0.0012 | 0.0788 ± 0.0008 | 0.0764 ± 0.0005 | 0.0749 ± 0.0005 |
| anchors 1% | 0.0811 ± 0.0010 | 0.0782 ± 0.0006 | 0.0764 ± 0.0006 | 0.0749 ± 0.0003 |
| anchors 2% | 0.0822 ± 0.0010 | 0.0780 ± 0.0008 | 0.0763 ± 0.0004 | 0.0748 ± 0.0003 |
| maximin anchors | 0.0815 ± 0.0005 | 0.0782 ± 0.0004 | 0.0764 ± 0.0004 | 0.0748 ± 0.0005 |
| ASDE relax | 0.0859 ± 0.0004 | 0.0790 ± 0.0007 | 0.0765 ± 0.0002 | 0.0749 ± 0.0002 |

## Candidate − random, path-level, 95% bootstrap CI over 32 paths (TUNE; negative = candidate better)
| candidate | budget | candidate − random | 95% CI | paths better (of 32) | CI below 0 |
|---|---|---|---|---|---|
| Hybrid (current) | 95 | +0.00436 | [+0.00313, +0.00561] | 3 | no |
| Hybrid (current) | 158 | +0.00102 | [+0.00046, +0.00159] | 8 | no |
| Hybrid (current) | 222 | +0.00051 | [+0.00023, +0.00079] | 9 | no |
| Hybrid (current) | 317 | +0.00063 | [+0.00031, +0.00097] | 9 | no |
| SVSS-cos | 95 | +0.00392 | [+0.00276, +0.00510] | 4 | no |
| SVSS-cos | 158 | +0.00080 | [+0.00009, +0.00152] | 8 | no |
| SVSS-cos | 222 | +0.00037 | [-0.00009, +0.00084] | 12 | no |
| SVSS-cos | 317 | +0.00044 | [+0.00002, +0.00089] | 13 | no |
| anchors 1% | 95 | -0.00045 | [-0.00114, +0.00028] | 19 | no |
| anchors 1% | 158 | +0.00019 | [-0.00047, +0.00092] | 16 | no |
| anchors 1% | 222 | +0.00032 | [-0.00020, +0.00088] | 16 | no |
| anchors 1% | 317 | +0.00042 | [+0.00002, +0.00083] | 11 | no |
| anchors 2% | 95 | +0.00063 | [-0.00031, +0.00156] | 13 | no |
| anchors 2% | 158 | +0.00000 | [-0.00050, +0.00051] | 15 | no |
| anchors 2% | 222 | +0.00027 | [-0.00007, +0.00063] | 13 | no |
| anchors 2% | 317 | +0.00034 | [+0.00001, +0.00070] | 10 | no |
| maximin anchors | 95 | -0.00000 | [-0.00080, +0.00080] | 17 | no |
| maximin anchors | 158 | +0.00016 | [-0.00048, +0.00082] | 16 | no |
| maximin anchors | 222 | +0.00033 | [-0.00016, +0.00084] | 17 | no |
| maximin anchors | 317 | +0.00034 | [-0.00007, +0.00077] | 14 | no |
| ASDE relax | 95 | +0.00436 | [+0.00311, +0.00565] | 3 | no |
| ASDE relax | 158 | +0.00102 | [+0.00048, +0.00159] | 8 | no |
| ASDE relax | 222 | +0.00046 | [+0.00019, +0.00071] | 8 | no |
| ASDE relax | 317 | +0.00049 | [+0.00019, +0.00080] | 10 | no |

For reference, the current hybrid against random in the same format is in `stage3/contrasts.csv` (rows with candidate = "Hybrid (current)").

## Candidate − current hybrid (information only, not part of the rule)
| candidate | budget | candidate − current hybrid | 95% CI | paths better (of 32) |
|---|---|---|---|---|
| SVSS-cos | 95 | -0.00044 | [-0.00142, +0.00054] | 17 |
| SVSS-cos | 158 | -0.00023 | [-0.00083, +0.00042] | 17 |
| SVSS-cos | 222 | -0.00015 | [-0.00056, +0.00028] | 17 |
| SVSS-cos | 317 | -0.00019 | [-0.00049, +0.00016] | 21 |
| anchors 1% | 95 | -0.00481 | [-0.00582, -0.00380] | 30 |
| anchors 1% | 158 | -0.00083 | [-0.00148, -0.00015] | 21 |
| anchors 1% | 222 | -0.00019 | [-0.00070, +0.00037] | 19 |
| anchors 1% | 317 | -0.00022 | [-0.00055, +0.00013] | 22 |
| anchors 2% | 95 | -0.00373 | [-0.00461, -0.00284] | 30 |
| anchors 2% | 158 | -0.00102 | [-0.00146, -0.00060] | 26 |
| anchors 2% | 222 | -0.00024 | [-0.00061, +0.00013] | 19 |
| anchors 2% | 317 | -0.00029 | [-0.00053, -0.00004] | 22 |
| maximin anchors | 95 | -0.00436 | [-0.00549, -0.00322] | 28 |
| maximin anchors | 158 | -0.00086 | [-0.00147, -0.00023] | 22 |
| maximin anchors | 222 | -0.00019 | [-0.00062, +0.00029] | 19 |
| maximin anchors | 317 | -0.00029 | [-0.00060, +0.00002] | 23 |
| ASDE relax | 95 | +0.00000 | [+0.00000, +0.00000] | 0 |
| ASDE relax | 158 | +0.00000 | [+0.00000, +0.00000] | 0 |
| ASDE relax | 222 | -0.00006 | [-0.00013, +0.00000] | 13 |
| ASDE relax | 317 | -0.00014 | [-0.00031, +0.00002] | 20 |

## Secondary tables
TUNE + FPGA 1, held-out test set:

| method | 95 | 158 | 222 | 317 |
|---|---|---|---|---|
| Random | 0.0820 ± 0.0007 | 0.0784 ± 0.0005 | 0.0765 ± 0.0005 | 0.0749 ± 0.0003 |
| SDE | 0.0816 ± 0.0004 | 0.0779 ± 0.0003 | 0.0763 ± 0.0003 | 0.0750 ± 0.0002 |
| Hybrid (current) | 0.0864 ± 0.0004 | 0.0795 ± 0.0007 | 0.0770 ± 0.0003 | 0.0755 ± 0.0003 |
| SVSS-cos | 0.0859 ± 0.0012 | 0.0793 ± 0.0008 | 0.0769 ± 0.0005 | 0.0753 ± 0.0005 |
| anchors 1% | 0.0815 ± 0.0010 | 0.0786 ± 0.0006 | 0.0768 ± 0.0005 | 0.0753 ± 0.0003 |
| anchors 2% | 0.0826 ± 0.0010 | 0.0785 ± 0.0007 | 0.0768 ± 0.0004 | 0.0752 ± 0.0003 |
| maximin anchors | 0.0820 ± 0.0005 | 0.0786 ± 0.0004 | 0.0768 ± 0.0004 | 0.0752 ± 0.0005 |
| ASDE relax | 0.0864 ± 0.0004 | 0.0795 ± 0.0007 | 0.0770 ± 0.0002 | 0.0754 ± 0.0002 |

TUNE 2–26, paper's old test set (all non-training sites):

| method | 95 | 158 | 222 | 317 |
|---|---|---|---|---|
| Random | 0.0820 ± 0.0007 | 0.0785 ± 0.0005 | 0.0765 ± 0.0005 | 0.0749 ± 0.0003 |
| SDE | 0.0816 ± 0.0005 | 0.0779 ± 0.0003 | 0.0763 ± 0.0003 | 0.0749 ± 0.0003 |
| Hybrid (current) | 0.0867 ± 0.0004 | 0.0797 ± 0.0005 | 0.0773 ± 0.0002 | 0.0757 ± 0.0002 |
| SVSS-cos | 0.0860 ± 0.0011 | 0.0794 ± 0.0006 | 0.0771 ± 0.0004 | 0.0755 ± 0.0003 |
| anchors 1% | 0.0816 ± 0.0009 | 0.0786 ± 0.0005 | 0.0769 ± 0.0005 | 0.0753 ± 0.0002 |
| anchors 2% | 0.0828 ± 0.0009 | 0.0787 ± 0.0006 | 0.0769 ± 0.0004 | 0.0753 ± 0.0002 |
| maximin anchors | 0.0820 ± 0.0004 | 0.0786 ± 0.0003 | 0.0769 ± 0.0004 | 0.0753 ± 0.0005 |
| ASDE relax | 0.0867 ± 0.0004 | 0.0797 ± 0.0005 | 0.0773 ± 0.0001 | 0.0756 ± 0.0001 |

## Selection facts
| design | anchors | GPR fits per path | additions passing the ASDE rule | window relaxation level (mean) |
|---|---|---|---|---|
| Hybrid (current) | 95 | 2.0 | 168.1 of 222 (75.7%) | – |
| SVSS-cos | 95 | 2.0 | 176.2 of 222 (79.4%) | – |
| anchors 1% | 32 | 0.0 | 44.8 of 285 (15.7%) | – |
| anchors 2% | 63 | 1.0 | 107.2 of 254 (42.2%) | – |
| maximin anchors | 95 | 0.0 | 18.5 of 222 (8.3%) | – |
| ASDE relax | 95 | 2.0 (the current hybrid's SVSS anchors are reused, no new fits) | 168.1 of 222 (75.7%) | 0.93 |

Hyperparameters at a bound (within 5%), predictor fits, all budgets and FPGAs:

| design | column lengthscale lower / upper | row lengthscale lower / upper | column variance lower / upper | row variance lower / upper | noise lower / upper |
|---|---|---|---|---|---|
| Random | 5.4% / 0.0% | 6.0% / 0.0% | 0.4% / 0.0% | 3.5% / 0.0% | 0.0% / 0.0% |
| SDE | 4.7% / 0.0% | 7.6% / 0.0% | 0.5% / 0.0% | 3.3% / 0.0% | 0.0% / 0.0% |
| Hybrid (current) | 6.5% / 0.1% | 5.5% / 0.0% | 1.1% / 0.0% | 2.1% / 0.0% | 0.0% / 0.0% |
| SVSS-cos | 6.3% / 0.0% | 6.1% / 0.0% | 1.1% / 0.0% | 2.6% / 0.0% | 0.0% / 0.0% |
| anchors 1% | 4.4% / 0.0% | 6.2% / 0.0% | 0.4% / 0.0% | 2.3% / 0.0% | 0.0% / 0.0% |
| anchors 2% | 4.4% / 0.1% | 5.5% / 0.0% | 0.7% / 0.0% | 2.1% / 0.0% | 0.0% / 0.0% |
| maximin anchors | 3.4% / 0.0% | 7.2% / 0.0% | 0.1% / 0.0% | 2.6% / 0.0% | 0.0% / 0.0% |
| ASDE relax | 6.7% / 0.1% | 5.5% / 0.0% | 1.1% / 0.0% | 2.2% / 0.0% | 0.0% / 0.0% |

Selection-loop fits (only designs with SVSS fits):

| design (selection-loop fits) | column lengthscale lower | row lengthscale lower | column variance lower | row variance lower | noise lower |
|---|---|---|---|---|---|
| SVSS-cos | 27.8% | 21.6% | 15.0% | 10.9% | 6.9% |
| anchors 2% | 32.5% | 26.2% | 21.9% | 13.1% | 13.8% |
