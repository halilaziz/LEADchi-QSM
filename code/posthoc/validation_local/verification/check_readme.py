#!/usr/bin/env python3
"""LEAD-chi -- check every number printed in validation_local/README.md against the INDEPENDENT recomputation
(verification/mine/*.csv and mine/scalars.json written by verify_independent.py). POST HOC check.

A printed value passes when it equals the recomputed value after rounding to the printed number of decimals
(tolerance 0.5 * 10^-decimals; for e-notation, half a unit in the last printed significant digit).
Writes ONLY verification/compare_readme.csv and verification/compare_readme_failures.csv.

    python3 check_readme.py
"""
import json
import math
import os
import re

import numpy as np
import pandas as pd

VAL = "/path/to/Penn_LEAD/LEADchi/results_posthoc/validation_local"
OUT = f"{VAL}/verification"
MINE = f"{OUT}/mine"
R = open(f"{VAL}/README.md", encoding="utf-8").read().split("\n")
L = {i + 1: s for i, s in enumerate(R)}          # 1-based line numbers as in the README

A = pd.read_csv(f"{MINE}/A_X2_chisep_group_effects.csv")
AD = pd.read_csv(f"{MINE}/A_X2_component_descriptives.csv")
B = pd.read_csv(f"{MINE}/B_nuisance_models.csv")
JC = pd.read_csv(f"{MINE}/B_joint_model_coefficients.csv")
R2T = pd.read_csv(f"{MINE}/B_r2s_models.csv")
TT = pd.read_csv(f"{MINE}/B_dice_tertile_cells.csv")
TR = pd.read_csv(f"{MINE}/B_dice_tertile_models.csv")
CS = pd.read_csv(f"{MINE}/C_low_TD_tail_sessions.csv", dtype={"code": str})
QC = pd.read_csv(f"{MINE}/C_low_TD_tail_qc_comparison.csv")
CR = pd.read_csv(f"{MINE}/C_low_TD_tail_refits.csv")
DR = pd.read_csv(f"{MINE}/D_age_slopes.csv")
IM = pd.read_csv(f"{MINE}/D_implied_chi_slope_from_R2s.csv")
SC = json.load(open(f"{MINE}/scalars.json"))
FROZEN = json.load(open("/path/to/Penn_LEAD/LEADchi/results_confirmatory/results.json"))

ROWS = []


def clean(s):
    s = s.replace("−", "-").replace("**", "").replace("s^-1", "s_inv").replace("mm^3", "mm3")
    s = re.sub(r"R2\*?", "Rtwo", s).replace("M1", "Mone").replace("HC3", "HCthree")
    s = re.sub(r"(?<=th)-(?=\d)", " ~ ", s)
    s = re.sub(r"(?<=\d)-(?=[\d.])", " ~ ", s)        # a hyphen between two numbers is a range separator
    return s


NUM = re.compile(r"-?(?:\d+\.\d+|\d+|\.\d+)(?:e-?\d+)?")


def nums(s):
    return NUM.findall(clean(s))


def ok(printed, val):
    if isinstance(val, str):
        return printed == val
    if val is None or (isinstance(val, float) and math.isnan(val)):
        return False
    x = float(printed)
    if "e" in printed:
        mant, ex = printed.split("e")
        sig = len(mant.replace("-", "").replace(".", "").lstrip("0"))
        tol = 0.5 * 10 ** (int(ex) - (sig - 1))
    else:
        dec = len(printed.split(".")[1]) if "." in printed else 0
        tol = 0.5 * 10 ** (-dec)
    return abs(val - x) <= tol * (1 + 1e-9) + 1e-12


def check(where, printed_list, values, what=""):
    if len(printed_list) != len(values):
        ROWS.append(dict(where=where, what=what, printed=" ".join(printed_list), recomputed=str(values),
                         status="PARSE-MISMATCH"))
        return
    for p, v in zip(printed_list, values):
        if v is None:
            continue
        st = "ok" if ok(p, v) else "MISMATCH"
        ROWS.append(dict(where=where, what=what, printed=p, recomputed=v if isinstance(v, str) else float(v), status=st))


def cells(line):
    return [c.strip() for c in L[line].split("|")]


def row(df, **kw):
    m = np.ones(len(df), bool)
    for k, v in kw.items():
        m &= (df[k] == v).to_numpy()
    assert m.sum() == 1, (kw, m.sum())
    return df[m].iloc[0]


def bci(r, pre=""):
    if pre:
        return [r[pre], r[pre + "_lo"], r[pre + "_hi"]]
    return [r["beta_PS"], r["ci_lo"], r["ci_hi"]]


def s1(r):
    return [r.PSonly_minus_TD, r.PSonly_minus_TD_lo, r.PSonly_minus_TD_hi, r.PSonly_minus_TD_p]


def s2(r):
    return [r.PSADHD_minus_ADHD, r.PSADHD_minus_ADHD_lo, r.PSADHD_minus_ADHD_hi, r.PSADHD_minus_ADHD_p]


def text(line, snippet, values, what=""):
    """prose check: `snippet` must occur in README line `line`; every number in it is compared in order."""
    assert snippet in L[line], (line, snippet)
    check(f"L{line}", nums(snippet), values, what or snippet)


M1 = row(B, model="M1 (frozen primary)")
fd = FROZEN["H1"]["decision"]
# ---------------------------------------------------------------- header and reproduction table
text(5, "beta_PS = +4.963 ppb (HC3 SE 3.125, p = .1169); outcome-B bound = -6.1318 ppb, t_B = 3.550, p_B = .0004; first-usable n = 80 (TD 22 / ADHD-only 16 / PS-only 17 / PS+ADHD 25)",
     [M1.beta_PS, M1.se_hc3, M1.p_two, SC["bound"], M1.tB_fixed_bound, M1.pB_fixed_bound, 80, 22, 16, 17, 25],
     "frozen restatement")
check("L27", nums(cells(27)[2]) + nums(cells(27)[3]), [M1.beta_PS, M1.se_hc3, fd["beta_ppb"], fd["se_hc3"]])
check("L28", nums(cells(28)[2]) + nums(cells(28)[3]), [SC["bound"], M1.tB_fixed_bound, fd["B_registered"]["bound_ppb"], fd["B_registered"]["tB"]])
sg = {v: row(B, model=f"M1 + {v} (alone)") for v in ("dice_min", "motion_gsr", "ref2_minus_wb_ppb", "icv_seg_ml")}
fs = {x["covariate"]: x for x in SC["frozen_sensitivity_compare"]}
check("L29", nums(cells(29)[2]) + nums(cells(29)[3]), [sg["dice_min"].beta_PS, sg["motion_gsr"].beta_PS, sg["icv_seg_ml"].beta_PS,
                                                       fs["dice_min"]["beta_frozen"], fs["motion_gsr"]["beta_frozen"], fs["icv_seg_ml"]["beta_frozen"]])
check("L30", nums(cells(30)[2]) + nums(cells(30)[3]), [sg["dice_min"].pB_frozen_style, sg["motion_gsr"].pB_frozen_style, sg["icv_seg_ml"].pB_frozen_style,
                                                       fs["dice_min"]["pB_frozen"], fs["motion_gsr"]["pB_frozen"], fs["icv_seg_ml"]["pB_frozen"]])
r2m = row(R2T, model="R2* M1")
check("L31", nums(cells(31)[2]) + nums(cells(31)[3]), [r2m.beta_PS, SC["K2_frozen"]["est"]])
put0 = row(DR, outcome="chi_PUT", sample="first-usable (OLS, HC3)", covariates="age only")
gp0 = row(DR, outcome="chi_GP", sample="first-usable (OLS, HC3)", covariates="age only")
check("L32", nums(cells(32)[2]) + nums(cells(32)[3]), [put0.slope_per_y, gp0.slope_per_y, SC["pc2x_frozen"]["chi_PUT"], SC["pc2x_frozen"]["chi_GP"]])

# ---------------------------------------------------------------- A.2 table
A274 = lambda roi, comp: row(A, roi=roi, component=comp, Dr=274.0, R2=0.0)  # noqa: E731
inp = lambda roi, what: row(A, model=f"{roi}: {what}")  # noqa: E731
A2 = {65: inp("SNVTA", "total chi (input)"), 66: inp("SNVTA", "R2* (input, s^-1)"), 67: A274("SNVTA", "chi_para"), 68: A274("SNVTA", "chi_dia"),
      69: inp("SNonly", "total chi (input)"), 70: inp("SNonly", "R2* (input, s^-1)"), 71: A274("SNonly", "chi_para"), 72: A274("SNonly", "chi_dia"),
      73: inp("GP", "total chi (input)"), 74: inp("GP", "R2* (input, s^-1)"), 75: A274("GP", "chi_para"), 76: A274("GP", "chi_dia")}
for ln, r in A2.items():
    c = cells(ln)
    q = [] if c[6] == "-" else [r.q_BH_family6]
    check(f"L{ln}", nums(c[3]) + nums(c[4]) + nums(c[5]) + nums(c[6]) + nums(c[7]) + nums(c[8]),
          bci(r) + [r.se_hc3, r.p_two] + q + s1(r) + s2(r), c[1] + " " + c[2])
# cell means
for ln, roi, comp in ((81, "SNVTA", "chi_para"), (82, "SNVTA", "chi_dia"), (83, "SNonly", "chi_para"), (84, "SNonly", "chi_dia"),
                      (85, "GP", "chi_para"), (86, "GP", "chi_dia")):
    r = A274(roi, comp)
    check(f"L{ln}", nums(L[ln].split(":")[1]), [r.mean_TD, r.mean_ADHDonly, r.mean_PSonly, r.mean_PSADHD], f"cell means {roi} {comp}")
prim = A[(A.block == "A_X2") & (A.Dr == 274.0) & (A.R2 == 0.0)]
text(78, "The smallest q is .096 (SN-only chi_para); SN-VTA chi_para has q .118", [prim.q_BH_family6.min(), A274("SNVTA", "chi_para").q_BH_family6])
ext_min = np.nanmin(np.r_[prim.q_BH_ext18_beta, prim.q_BH_ext18_PSonly, prim.q_BH_ext18_PSADHD])
text(78, "An extended family of 18 tests (the same 6 beta_PS plus their 12 simple effects) gives a smallest q of .137",
     [18, 6, 12, ext_min])
# A.1 declarations
assert "The maximum numerical deviation across all rows is 3.4e-14" in L[56]
ROWS.append(dict(where="L56", what="X2 identity max error (floating-point noise)", printed="3.4e-14",
                 recomputed=SC["identity_max_abs_err"], status="ok (noise; both < 1e-13)" if SC["identity_max_abs_err"] < 1e-13 else "MISMATCH"))
assert "beta_PS differs by at most 1.1e-14 and SE_HC3 by at most 3.6e-15" in L[57]
ROWS.append(dict(where="L57", what="R2 cancellation ranges of beta / SE (floating-point noise)", printed="1.1e-14 / 3.6e-15",
                 recomputed=f"{SC['R2_cancellation_max_beta_range']:.2e} / {SC['R2_cancellation_max_se_range']:.2e}",
                 status="ok (noise; both < 1e-13)" if max(SC['R2_cancellation_max_beta_range'], SC['R2_cancellation_max_se_range']) < 1e-13 else "MISMATCH"))
text(57, "At R2 = 10 s_inv, 61 of 80 SN-VTA sessions (76%)".replace("s_inv", "s^-1"), [10, SC["n_dia_pos_R2_10_SNVTA"], 80, SC["n_dia_pos_R2_10_SNVTA"] / 80 * 100])
for roi, comp_p, comp_d, ln in (("SNVTA", "chi_para_unref", "chi_dia_unref", 58),):
    rp = row(A, roi=roi, component=comp_p)
    rd = row(A, roi=roi, component=comp_d)
    text(ln, "SN-VTA chi_para 4.485 [0.214, 8.757], p .0399; chi_dia 0.477 [-2.373, 3.328], p .7394",
         bci(rp) + [rp.p_two] + bci(rd) + [rd.p_two])
text(58, "SN-only 5.708 / 0.166; GP 1.575 / 0.699",
     [row(A, roi="SNonly", component="chi_para_unref").beta_PS, row(A, roi="SNonly", component="chi_dia_unref").beta_PS,
      row(A, roi="GP", component="chi_para_unref").beta_PS, row(A, roi="GP", component="chi_dia_unref").beta_PS])
d274 = {r.roi: r for r in AD[AD.Dr == 274.0].itertuples()}
text(59, "At Dr = 274 they are met by every SN-VTA and SN-only session and by 78 of 80 GP sessions",
     [274, 80 - d274["GP"].n_dia_gt0 - d274["GP"].n_para_lt0, 80])
ROWS.append(dict(where="L59", what="positivity met by every SN-VTA / SN-only session at Dr 274",
                 printed="all", recomputed=f"SNVTA viol {d274['SNVTA'].n_para_lt0 + d274['SNVTA'].n_dia_gt0}, SNonly viol {d274['SNonly'].n_para_lt0 + d274['SNonly'].n_dia_gt0}",
                 status="ok" if (d274['SNVTA'].n_para_lt0 + d274['SNVTA'].n_dia_gt0 + d274['SNonly'].n_para_lt0 + d274['SNonly'].n_dia_gt0) == 0 else "MISMATCH"))

# ---------------------------------------------------------------- A.3 table
DRMAP = {"137": 137.0, "274": 274.0, "321": 321.0, "548": 548.0}
ROIMAP = {"SN-VTA": "SNVTA", "SN-only": "SNonly", "GP": "GP"}
for ln in range(92, 104):
    c = cells(ln)
    roi, Dr = ROIMAP[c[1]], DRMAP[c[2].replace("*", "")]
    rp, rd = row(A, roi=roi, component="chi_para", Dr=Dr, R2=0.0), row(A, roi=roi, component="chi_dia", Dr=Dr, R2=0.0)
    ndesc = row(AD, roi=roi, Dr=Dr)
    check(f"L{ln}", nums(c[3]) + nums(c[4]) + nums(c[5]) + nums(c[6]) + nums(c[7]),
          bci(rp) + [rp.p_two] + bci(rd) + [rd.p_two] + [rp.PSonly_minus_TD, rp.PSonly_minus_TD_p, rd.PSonly_minus_TD, rd.PSonly_minus_TD_p,
                                                          int(ndesc.n_dia_gt0), 80], f"{roi} Dr {Dr:g}")
text(105, "Mean R2* is 22.92 s^-1 for SN-VTA, 23.82 for SN-only and 30.53 for GP".replace("R2*", "R2*"),
     [d274["SNVTA"].mean_R2s, d274["SNonly"].mean_R2s, d274["GP"].mean_R2s])
text(106, "At Dr = 274: SN-VTA chi_para 68.32 (9.13), chi_dia -15.33 (5.14); SN-only 74.08 (10.52) and -12.87 (5.37); GP 102.12 (10.77) and -9.29 (4.81)",
     [274] + sum([[d274[k].mean_para, d274[k].sd_para, d274[k].mean_dia, d274[k].sd_dia] for k in ("SNVTA", "SNonly", "GP")], []))
d548 = {r.roi: r for r in AD[AD.Dr == 548.0].itertuples()}
d137 = {r.roi: r for r in AD[AD.Dr == 137.0].itertuples()}
text(107, "At Dr = 548, chi_dia is positive in most sessions (SN-VTA 68/80, GP 80/80)", [548, d548["SNVTA"].n_dia_gt0, 80, d548["GP"].n_dia_gt0, 80])
text(108, "At Dr = 137, chi_para and chi_dia are negatively correlated (SN-VTA r = -0.44). At 274 the correlation is +0.33",
     [137, d137["SNVTA"].r_para_dia, 274, d274["SNVTA"].r_para_dia])
text(112, "chi_para +4.49 ppb [0.23, 8.74] and chi_dia +0.48 [-2.39, 3.35]", bci(A274("SNVTA", "chi_para")) + bci(A274("SNVTA", "chi_dia")))
text(112, "(+5.71 vs +0.17)", [A274("SNonly", "chi_para").beta_PS, A274("SNonly", "chi_dia").beta_PS])
text(112, "smallest q .096", [prim.q_BH_family6.min()])

# ---------------------------------------------------------------- B.1
BMAP = {122: "M1 (frozen primary)", 123: "M1 + dice_min (alone)", 124: "M1 + motion_gsr (alone)", 125: "M1 + ref2_minus_wb_ppb (alone)",
        126: "M1 + icv_seg_ml (alone)", 127: "M1 + Dice + motion + CSF offset + ICV (joint)"}
for ln, mdl in BMAP.items():
    c = cells(ln)
    r = row(B, model=mdl)
    check(f"L{ln}", sum((nums(c[i]) for i in range(2, 12)), []),
          [r.n] + bci(r) + [r.se_hc3, r.t, r.df, r.p_two, r.PSonly_minus_TD, r.PSonly_minus_TD_p, r.PSADHD_minus_ADHD, r.PSADHD_minus_ADHD_p,
                            r.tB_fixed_bound, r.pB_fixed_bound, r.pB_frozen_style], mdl)
JMAP = {133: "scale(dice_min)", 134: "scale(motion_gsr)", 135: "scale(ref2_minus_wb_ppb)", 136: "scale(icv_seg_ml)", 137: "ADHD_e"}
for ln, term in JMAP.items():
    c = cells(ln)
    r = row(JC, term=term)
    check(f"L{ln}", nums(c[2]) + nums(c[3]) + nums(c[4]), [r.est, r.ci_lo, r.ci_hi, r.se_hc3, r.p_two], term)
text(139, "The variance inflation of PS_e is 1.18 in M1 and 1.27 in the joint model", [SC["vif_PS_M1"], SC["vif_PS_joint"]])
text(139, "Spearman rho = .358 (p .0011)", [SC["spearman_dice_chi"]["rho"], SC["spearman_dice_chi"]["p"]])
for ln, mdl in ((145, "R2* M1"), (146, "R2* M1 + Dice"), (147, "R2* M1 + joint nuisance")):
    c = cells(ln)
    r = row(R2T, model=mdl)
    check(f"L{ln}", nums(c[2]) + nums(c[3]) + nums(c[4]) + nums(c[5]) + nums(c[6]),
          bci(r) + [r.se_hc3, r.p_two, r.PSonly_minus_TD, r.PSonly_minus_TD_p, r.PSADHD_minus_ADHD, r.PSADHD_minus_ADHD_p], mdl)
rd = row(R2T, model="R2* M1 + Dice")
text(149, "Dice is not related to R2* (+0.272 s^-1 per SD, p .3624)", [rd.dice_coef, rd.dice_coef_p])
text(149, "The R2* estimate barely moves (1.10 to 1.01)", [r2m.beta_PS, rd.beta_PS])
text(149, "one-sided for *lower* R2* (p .9825)", [FROZEN["H1"]["K2"]["p_one_sided_lower"]])
text(151, "label-free cut points 0.8054 and 0.8233 over the 80", SC["dice_tertile_cuts"] + [80])
for ln, tl in ((155, "T1_low"), (156, "T2_mid"), (157, "T3_high")):
    c = cells(ln)
    r = row(TT, tertile=tl)
    lo, hi = (float(x) for x in r.dice_range.split("-"))
    check(f"L{ln}", nums(c[2]) + nums(c[3]) + nums(c[4]) + nums(c[5]),
          [lo, hi, r.TD + r["ADHD-only"] + r["PS-only"] + r["PS+ADHD"], r.TD, r["ADHD-only"], r["PS-only"], r["PS+ADHD"], r.mean_chi], tl)
text(159, "Tertile by cell: Fisher p = .32. Tertile by PS: Fisher p = .36", [SC["fisher_p_cell"], SC["fisher_p_PS"]])
TMAP = {163: "within T1_low : full M1", 164: "within T1_low : reduced (PS*ADHD + age + sex)", 165: "within T2_mid : full M1",
        166: "within T2_mid : reduced (PS*ADHD + age + sex)", 167: "within T3_high : full M1", 168: "within T3_high : reduced (PS*ADHD + age + sex)",
        169: "pooled M1 + Dice-tertile fixed effects", 170: "pooled M1 + Dice tertile FE + linear Dice",
        171: "quality-matched: drop lowest Dice tertile (M1)", 172: "quality-matched: drop lowest Dice tertile (M1 + Dice)"}
for ln, mdl in TMAP.items():
    c = cells(ln)
    r = row(TR, model=mdl)
    ncell = [r.n] + ([int(x) for x in r.cells.split("/")] if "(" in c[2] else [])
    check(f"L{ln}", sum((nums(c[i]) for i in range(2, 10)), []),
          ncell + bci(r) + [r.se_hc3, r.df, r.p_two, r.PSonly_minus_TD, r.PSonly_minus_TD_p, r.PSADHD_minus_ADHD, r.PSADHD_minus_ADHD_p,
                            r.pB_fixed_bound], mdl)
h = SC["het"]
text(174, "F(2, 64) = 2.44, p = .095", [2, h["df2"], h["F"], h["p_F"]])
text(174, "With 26-27 sessions per tertile and 12 M1 parameters", [26, 27, row(TR, model="within T1_low : full M1").n - row(TR, model="within T1_low : full M1").df])
J = row(B, model="M1 + Dice + motion + CSF offset + ICV (joint)")
text(178, "lowers beta_PS from 4.96 to 4.16 ppb (p .175). The outcome-B bound is still excluded (p_B .0006)",
     [M1.beta_PS, J.beta_PS, J.p_two, J.pB_fixed_bound])
text(178, "(+4.3 ppb per SD)", [row(JC, term="scale(dice_min)").est])
t1r = row(TT, tertile="T1_low")
text(178, "holds half of the TD sessions (11/22) but only 2 of the 17 PS-only sessions", [t1r.TD, 22, t1r["PS-only"], 17])
text(178, "(Fisher p .32)", [SC["fisher_p_cell"]])
qmr = row(TR, model="quality-matched: drop lowest Dice tertile (M1)")
text(178, "(+10.4, p .042)", [qmr.beta_PS, qmr.p_two])
text(178, "(-0.95)", [row(TR, model="within T1_low : full M1").beta_PS])
text(178, "(p .095)", [h["p_F"]])
pooled = [M1.beta_PS] + [sg[k].beta_PS for k in sg] + [J.beta_PS, row(TR, model="pooled M1 + Dice-tertile fixed effects").beta_PS,
                                                       row(TR, model="pooled M1 + Dice tertile FE + linear Dice").beta_PS]
text(178, "between +3.7 and +5.5 ppb", [min(pooled), max(pooled)])
dice_alone_drop = M1.beta_PS - sg["dice_min"].beta_PS
dice_joint_drop = row(B, model="M1 + Dice + motion + CSF offset + ICV (joint)").beta_PS
ROWS.append(dict(where="L178", what="'Dice accounts for about 1 ppb of the estimate, alone or jointly'", printed="about 1",
                 recomputed=f"alone {dice_alone_drop:.3f}; joint model total drop {M1.beta_PS - dice_joint_drop:.3f}",
                 status="ok" if 0.6 < dice_alone_drop < 1.5 else "MISMATCH"))

# ---------------------------------------------------------------- C
for ln in range(188, 192):
    c = cells(ln)
    r = row(CS, code=c[1])
    extra = [0, 80] if c[1] == "a516f3" else []
    check(f"L{ln}", nums(c[3]) + nums(c[4]) + nums(c[5]) + nums(c[6]) + nums(c[7]) + nums(c[8]) + nums(c[9]) + nums(c[10]) + nums(c[11]),
          [r.chi_SNVTA, r.qc_grade, r.dice_min, r.pctTDrest_dice_min] + extra[:1] * 0 + ([] if not extra else [80]) +
          [r.motion_gsr, r.pctTDrest_motion_gsr, r.fimp_SNVTA, r.inmask_SNVTA, r.Q_phase, r.r2s_SNVTA, r.pctTDrest_r2s_SNVTA, r.rstud_M1],
          c[1])
QMAP = {197: "dice_min", 198: "motion_gsr", 200: "inmask_SNVTA", 201: "Q_phase", 202: "roi_vol", 203: "r2s_SNVTA", 204: "ref2_minus_wb_ppb",
        205: "icv_seg_ml", 206: "age"}
for ln, v in QMAP.items():
    c = cells(ln)
    r = row(QC, metric=v)
    check(f"L{ln}", nums(c[2]) + nums(c[3]) + nums(c[4]) + nums(c[5]),
          [r.tail_mean, r.tail_min, r.tail_max, r.TDrest_median, r.TDrest_q25, r.TDrest_q75, r.TDrest_min, r.TDrest_max, r.all80_median,
           r.wilcox_p_tail_vs_TDrest], v)
qg = SC["qc_grade"]
text(208, "all 4 tail sessions are grade 1", [4, 1] if qg["tail"].get("1.0", 0) == 4 else [4, -1])
text(208, "The other 18 TD sessions are 2 grade 0 and 16 grade 1; all 80 are 4 grade 0 and 76 grade 1",
     [18, qg["TD_rest"].get("0.0", 0), 0, qg["TD_rest"].get("1.0", 0), 1, 80, qg["all80"].get("0.0", 0), 0, qg["all80"].get("1.0", 0), 1])
for ln, idx in ((212, 0), (213, 1), (214, 2)):
    c = cells(ln)
    r = CR.iloc[idx]
    pre = [80] if ln == 212 else []
    if ln == 214:
        pre = [35, 4, 1, float(r.dropped_min_chi), 80]
    elif ln == 213:
        pre = [4]
    check(f"L{ln}", nums(c[1]) + sum((nums(c[i]) for i in range(2, 11)), []),
          pre + [r.n] + [int(x) for x in r.cells.split("/")] + bci(r) + [r.se_hc3, r.p_two, r.PSonly_minus_TD, r.PSonly_minus_TD_p,
                                                                        r.PSADHD_minus_ADHD, r.PSADHD_minus_ADHD_p, r.tB_fixed_bound,
                                                                        r.pB_fixed_bound, r.huber_beta], r.model)
text(216, "(|studentised residual| <= 1.96)".replace("<=", "≤"), [CS.rstud_M1.abs().max()])
text(216, "their mean Dice is 0.771 versus 0.811 (p .019)", [row(QC, metric="dice_min").tail_mean, row(QC, metric="dice_min").TDrest_median,
                                                             row(QC, metric="dice_min").wilcox_p_tail_vs_TDrest],
     "0.811 = mean (0.8112) and also median (0.8107) of the other 18 TD sessions")
text(216, "(89th-100th TD percentile)", [sorted(CS.pctTDrest_motion_gsr)[1], CS.pctTDrest_motion_gsr.max()])
text(216, "(0th and 6th TD percentile)", [sorted(CS.pctTDrest_r2s_SNVTA)[0], sorted(CS.pctTDrest_r2s_SNVTA)[1]])
text(216, "lowers beta_PS from 4.96 to 3.33 ppb, and PS-only - TD from 8.16 (p .044) to 4.36 (p .276). The label-free cut lowers beta_PS to 1.45",
     [CR.iloc[0].beta_PS, CR.iloc[1].beta_PS, CR.iloc[0].PSonly_minus_TD, CR.iloc[0].PSonly_minus_TD_p, CR.iloc[1].PSonly_minus_TD,
      CR.iloc[1].PSonly_minus_TD_p, CR.iloc[2].beta_PS])
text(216, "(p_B ≤ .0022)", [CR.pB_fixed_bound.max()])
ROWS.append(dict(where="L216", what="'about a third ... (two-thirds with the label-free cut)'", printed="1/3 ; 2/3",
                 recomputed=f"{1 - CR.iloc[1].beta_PS / CR.iloc[0].beta_PS:.3f} ; {1 - CR.iloc[2].beta_PS / CR.iloc[0].beta_PS:.3f}",
                 status="ok (2/3 is a slight understatement of 0.71)"))

# ---------------------------------------------------------------- D.1
SAMP = {"first-usable": "first-usable (OLS, HC3)", "all usable, lme": "all usable sessions (lme REML, random intercept)",
        "lme, within-person": "all usable (lme REML), WITHIN-person age", "lme, between-person": "all usable (lme REML), BETWEEN-person mean age"}
cur = None
for ln in range(233, 265):
    c = cells(ln)
    o = c[1].replace("**", "").split(" ")[0]
    r = row(DR, outcome=o, sample=SAMP[c[2]], covariates=c[3])
    pr = nums(c[4]) + nums(c[5]) + nums(c[6])
    exp = [r.slope_per_y, r.ci_lo, r.ci_hi, r.p_two, r.slope_per_y_in_SD]
    check(f"L{ln}", pr, exp, f"{o} | {c[2]} | {c[3]}")
    if c[7]:
        flag = c[7].split(" ")[0]
        mine_flag = "yes" if r.ci_excludes_expected_range == 1 else "no"
        ROWS.append(dict(where=f"L{ln}", what=f"{o} | {c[2]} | {c[3]} | CI excludes registered range", printed=flag, recomputed=mine_flag,
                         status="ok" if flag == mine_flag else "MISMATCH"))
text(226, "putamen chi **1.33-1.74 ppb/y**; GP chi **2.26-3.54 ppb/y**", [2.6 / 1.95, 3.4 / 1.95, 4.4 / 1.95, 6.9 / 1.95])
text(226, "2.6-3.4 ppb putamen and 4.4-6.9 ppb GP over 1.95 y", [2.6, 3.4, 4.4, 6.9, 1.95])
text(224, "26 (107 - 80 - 1)", [26, 107, 80, 1])
ar = SC["correlations"]
text(224, "age 8.8-19.5 y", ar["age_range_usable"])
text(224, "(n = 107, 80 subjects", [107, 80])
text(223, "(n = 80; age 8.8-17.8 y)", [80] + ar["age_range_first"])
text(16, "26 subjects have 2 or more usable sessions", [ar["n_subjects_with_2plus_usable"], 2])
# ---------------------------------------------------------------- D.2
IMAP = {270: ("PUT", "first-usable (OLS, HC3)", "age only"), 271: ("PUT", "all usable sessions (lme REML, random intercept)", "+ sex + cell"),
        272: ("PUT", "all usable (lme REML), WITHIN-person age", "+ sex + cell"), 273: ("GP", "first-usable (OLS, HC3)", "age only"),
        274: ("GP", "all usable sessions (lme REML, random intercept)", "+ sex + cell"), 275: ("GP", "all usable (lme REML), WITHIN-person age", "+ sex + cell")}
for ln, (roi, s, cv) in IMAP.items():
    c = cells(ln)
    r = row(IM, roi=roi, sample=s, covariates=cv)
    check(f"L{ln}", nums(c[3]) + nums(c[4]) + nums(c[5]),
          [r.slope_per_y_r2s, r.ci_lo_r2s, r.ci_hi_r2s, r.implied_chi_slope_ppb_per_y, r.implied_lo, r.implied_hi, r.slope_per_y_chi,
           r.ci_lo_chi, r.ci_hi_chi], f"{roi} {s} {cv}")
dp0 = row(DR, outcome="chi_dia_PUT", sample=SAMP["first-usable"], covariates="age only")
dpl = row(DR, outcome="chi_dia_PUT", sample=SAMP["all usable, lme"], covariates="age only")
dplc = row(DR, outcome="chi_dia_PUT", sample=SAMP["all usable, lme"], covariates="+ sex + cell")
text(278, "**-0.358 [-0.715, -0.001] ppb/y, p .0492** (first-usable); **-0.501 [-0.840, -0.161], p .0055** (all usable, lme)",
     [dp0.slope_per_y, dp0.ci_lo, dp0.ci_hi, dp0.p_two, dpl.slope_per_y, dpl.ci_lo, dpl.ci_hi, dpl.p_two],
     "chi_dia_PUT (the lme value is the AGE-ONLY lme; the D.2 table above uses the '+ sex + cell' lme)")
dg0 = row(DR, outcome="chi_dia_GP", sample=SAMP["first-usable"], covariates="age only")
dgl = row(DR, outcome="chi_dia_GP", sample=SAMP["all usable, lme"], covariates="age only")
text(279, "-0.056 [-0.561, 0.449], p .83 (first-usable); -0.111 [-0.555, 0.333], p .61", [dg0.slope_per_y, dg0.ci_lo, dg0.ci_hi, dg0.p_two,
                                                                                          dgl.slope_per_y, dgl.ci_lo, dgl.ci_hi, dgl.p_two])
pp0 = row(DR, outcome="chi_para_PUT", sample=SAMP["first-usable"], covariates="age only")
pg0 = row(DR, outcome="chi_para_GP", sample=SAMP["first-usable"], covariates="age only")
text(280, "+0.740 [0.335, 1.145], p .0005", [pp0.slope_per_y, pp0.ci_lo, pp0.ci_hi, pp0.p_two])
text(280, "+2.453 [1.446, 3.459], p < .0001", [pg0.slope_per_y, pg0.ci_lo, pg0.ci_hi, None])
ROWS.append(dict(where="L280", what="GP chi_para p < .0001", printed="<.0001", recomputed=float(pg0.p_two),
                 status="ok" if pg0.p_two < 1e-4 else "MISMATCH"))
text(283, "putamen r = .21 (age-partialled .16); GP r = .68", [ar["r_chi_r2s_first"], ar["r_chi_r2s_first_age_partial"], ar["r_chi_r2s_GP_first"]])
# ---------------------------------------------------------------- D.3
D3 = {289: [("ref_ppb", "age only")], 290: [("ref2_minus_wb_ppb", "age only")], 291: [("chi_PUT_unref", "age only")],
      292: [("chi_PUT_csfref", "age only")], 293: [("chi_GP_unref", "age only"), ("chi_GP_csfref", "age only")],
      294: [("chi_SNVTA", "age only"), ("r2s_SNVTA", "age only")]}
for ln, lst in D3.items():
    c = cells(ln)
    rs = [row(DR, outcome=o, sample=SAMP["first-usable"], covariates=cv) for o, cv in lst]
    exp = sum(([r.slope_per_y, r.ci_lo, r.ci_hi] for r in rs), []) + [r.p_two for r in rs]
    check(f"L{ln}", nums(c[2]) + nums(c[3]), exp, " / ".join(o for o, _ in lst))
# ---------------------------------------------------------------- D.4-5 prose
rp0 = row(DR, outcome="r2s_PUT", sample=SAMP["first-usable"], covariates="age only")
rpw = row(DR, outcome="r2s_PUT", sample=SAMP["lme, within-person"], covariates="+ sex + cell")
rg0 = row(DR, outcome="r2s_GP", sample=SAMP["first-usable"], covariates="age only")
text(302, "Our putamen R2* rise of 0.30 s^-1/y (0.25 SD/y) over 8.8-17.8 y".replace("R2*", "R2*"), [rp0.slope_per_y, rp0.slope_per_y_in_SD] + ar["age_range_first"])
text(306, "Putamen R2* rises by 0.30 s^-1/y (p ≈ 1e-6), including within persons (0.30, p .007)", [rp0.slope_per_y, rp0.p_two, rpw.slope_per_y, rpw.p_two])
text(306, "Its standardised rate (0.25 SD/y) matches GP R2* (0.24 SD/y)", [rp0.slope_per_y_in_SD, rg0.slope_per_y_in_SD])
text(306, "(8.8-17.8 / 19.5 y)", ar["age_range_first"] + [ar["age_range_usable"][1]])
pl = [row(DR, outcome="chi_PUT", sample=SAMP["all usable, lme"], covariates=cv).slope_per_y for cv in ("age only", "+ sex", "+ sex + cell")]
text(307, "Its slope is 0.38 ppb/y (first-usable) and 0.08-0.15 ppb/y (all usable)", [put0.slope_per_y, min(pl), max(pl)])
pi0 = row(IM, roi="PUT", sample=SAMP["first-usable"], covariates="age only")
text(307, "below the 1.10 ppb/y implied by its own R2* slope (shortfall test p .049 first-usable, p .0055 all usable)",
     [pi0.implied_chi_slope_ppb_per_y, dp0.p_two, dpl.p_two])
gi0 = row(IM, roi="GP", sample=SAMP["first-usable"], covariates="age only")
text(308, "observed 2.40 vs implied 2.51 ppb/y", [gi0.slope_per_y_chi, gi0.implied_chi_slope_ppb_per_y])
text(308, "(r = .68 vs .21 in putamen)", [ar["r_chi_r2s_GP_first"], ar["r_chi_r2s_first"]])
text(312, "(chi_para up 0.74, chi_dia down 0.36 ppb/y)", [pp0.slope_per_y, -dp0.slope_per_y])
snr = [row(DR, outcome="r2s_SNVTA", sample=s, covariates=cv).slope_per_y for s, cv in
       ((SAMP["first-usable"], "age only"), (SAMP["all usable, lme"], "age only"))]
snc = [row(DR, outcome="chi_SNVTA", sample=s, covariates=cv).slope_per_y for s, cv in
       ((SAMP["first-usable"], "age only"), (SAMP["all usable, lme"], "age only"))]
text(313, "age-related R2* (0.31-0.50 s^-1/y) and chi (about 1.2-1.3 ppb/y)".replace("R2*", "R2*"), [snr[0], snr[1], min(snc), max(snc)])

RES = pd.DataFrame(ROWS)
RES.to_csv(f"{OUT}/compare_readme.csv", index=False)
bad = RES[~RES.status.str.startswith("ok")]
bad.to_csv(f"{OUT}/compare_readme_failures.csv", index=False)
print(f"README numbers checked: {len(RES)}; not ok: {len(bad)}")
pd.set_option("display.width", 250)
pd.set_option("display.max_colwidth", 90)
print(bad.to_string())
