#!/usr/bin/env python3
"""LEAD-chi validation_local -- POST HOC check (2026-10-07) of the two CSVs added in the revision.

1. D_putamen_shortfall_by_Dr.csv (written by the R script) vs the verifier's independent recomputation
   verification/mine/D_putamen_shortfall_by_Dr.csv (written by verify_independent.py).
2. A_X2_BH_family_sensitivity.csv vs an own Benjamini-Hochberg implementation applied to the verifier's recomputed
   X2 p values (verification/mine/A_X2_chisep_group_effects.csv) and the 5 registered Tier-4 p values
   (results_exploratory/tier4_exploratory/fdr_table.csv, read-only).
Writes ONLY into verification/rerun_2026-10-07/.

    python3 check_new_outputs.py
"""
import os

import numpy as np
import pandas as pd

VAL = "/path/to/Penn_LEAD/LEADchi/results_posthoc/validation_local"
OUT = f"{VAL}/verification/rerun_2026-10-07"
T4 = "/path/to/Penn_LEAD/LEADchi/results_exploratory/tier4_exploratory/fdr_table.csv"
assert os.path.realpath(OUT).startswith(os.path.realpath(VAL) + "/verification/")
ROWS = []


def bh(p):
    p = np.asarray(p, float)
    m = len(p)
    o = np.argsort(p)
    q = p[o] * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(m)
    out[o] = np.minimum(q, 1.0)
    return out


def cmp(what, theirs, mine, rtol):
    ad = abs(float(theirs) - float(mine))
    rel = ad / max(abs(float(theirs)), 1e-300)
    ROWS.append(dict(check=what, theirs=float(theirs), mine=float(mine), abs_diff=ad, rel_diff=rel, rtol=rtol,
                     status="ok" if (ad <= 1e-9 or rel <= rtol) else "OUTSIDE"))


# 1. shortfall by Dr
th = pd.read_csv(f"{VAL}/D_putamen_shortfall_by_Dr.csv")
mi = pd.read_csv(f"{VAL}/verification/mine/D_putamen_shortfall_by_Dr.csv")
mt = th.merge(mi, on=["Dr", "sample", "covariates"], suffixes=("_R", "_py"), how="outer", indicator=True)
assert (mt["_merge"] == "both").all() and len(mt) == 24, "row mismatch"
for _, r in mt.iterrows():
    lme = not r["sample"].startswith("first")
    for c in ("slope", "ci_lo", "ci_hi", "p_two"):
        # OLS rows must agree to 1e-6; lme rows to 1e-4 (nlminb default tolerance, see README "Verification after the fixes")
        cmp(f"shortfall Dr {r.Dr:g} | {r['sample']} | {r.covariates} | {c}", r[c + "_R"], r[c + "_py"], 1e-4 if lme else 1e-6)

# 2. BH family sensitivity
A = pd.read_csv(f"{VAL}/verification/mine/A_X2_chisep_group_effects.csv")
prim = A[(A.block == "A_X2") & (A.Dr == 274.0) & (A.R2 == 0.0)].reset_index(drop=True)
t4 = pd.read_csv(T4)
regx = t4[t4.family == "registered X"].p_two.to_numpy(float)
assert len(regx) == 5
p6 = prim.p_two.to_numpy(float)
p18 = np.r_[p6, prim.PSonly_minus_TD_p.to_numpy(float), prim.PSADHD_minus_ADHD_p.to_numpy(float)]
q6, q18 = bh(p6), bh(p18)
q11, q23 = bh(np.r_[p6, regx])[:6], bh(np.r_[p18, regx])[:18]
S = pd.read_csv(f"{VAL}/A_X2_BH_family_sensitivity.csv")
for i, r in S.iterrows():
    j = prim.index[(prim.roi == r.roi) & (prim.component == r.component)][0]
    key = f"{r.roi} {r.component}"
    cmp(f"BH {key} q_X2_family6", r.q_X2_family6, q6[j], 1e-6)
    cmp(f"BH {key} q_X2_plus_regX_11", r.q_X2_plus_regX_11, q11[j], 1e-6)
    cmp(f"BH {key} q_X2_ext18", r.q_X2_ext18, q18[j], 1e-6)
    cmp(f"BH {key} q_X2_ext18_plus_regX_23", r.q_X2_ext18_plus_regX_23, q23[j], 1e-6)
    cmp(f"BH {key} q_PSonly_ext18", r.q_PSonly_ext18, q18[6 + j], 1e-6)
    cmp(f"BH {key} q_PSonly_ext18_plus_regX_23", r.q_PSonly_ext18_plus_regX_23, q23[6 + j], 1e-6)
    cmp(f"BH {key} q_PSADHD_ext18", r.q_PSADHD_ext18, q18[12 + j], 1e-6)
    cmp(f"BH {key} q_PSADHD_ext18_plus_regX_23", r.q_PSADHD_ext18_plus_regX_23, q23[12 + j], 1e-6)
R = pd.DataFrame(ROWS)
R.to_csv(f"{OUT}/compare_new_outputs.csv", index=False)
bad = R[R.status != "ok"]
ols = R[R.check.str.contains("first-usable")]
lme = R[R.check.str.contains("lme")]
bhr = R[R.check.str.startswith("BH")]
lines = [f"checks: {len(R)}; outside tolerance: {len(bad)}",
         f"shortfall OLS rows: max rel diff {ols.rel_diff.max():.2e}",
         f"shortfall lme rows: max rel diff {lme.rel_diff.max():.2e}; max abs diff {lme.abs_diff.max():.2e}",
         f"BH q values: max abs diff {bhr.abs_diff.max():.2e}",
         f"smallest q: X2-6 {q6.min():.4f}; X2+regX-11 {q11.min():.4f}; X2-18 {q18.min():.4f}; X2+regX-23 {q23.min():.4f}"]
with open(f"{OUT}/check_new_outputs.log", "w") as fh:
    fh.write("\n".join(lines) + "\n")
print("\n".join(lines))
if len(bad):
    print(bad.to_string())
