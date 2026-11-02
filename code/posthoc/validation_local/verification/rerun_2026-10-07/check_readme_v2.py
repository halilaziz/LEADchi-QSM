#!/usr/bin/env python3
"""LEAD-chi -- check every number printed in the REVISED validation_local/README.md (2026-10-07) against the
INDEPENDENT recomputation (verification/mine/*.csv and mine/scalars.json written by verify_independent.py). POST HOC.

ADAPTED COPY of the reviewer's verification/check_readme.py (unchanged there). Changes, all for the revised README:
  * README lines are addressed by their PRE-FIX line numbers (archive_v1_prefix_2026-10-07/README.md); each one is
    mapped to the revised README by identical text or, for edited lines, by a unique anchor (ANCHOR below). The
    output reports both numbers;
  * checks of reworded sentences use the new wording; checks for new sentences and tables are added at the end;
  * the A.2 simple-effect cells now also carry BH q values (family of 18), which are checked.
Round 2 (2026-10-07; the round-1 copy is in archive_v2_round1_2026-10-07/): pre-fix line 58 (declaration 5) is mapped by
an anchor because it now also labels its p values nominal and gives q; checks were added for that sentence, for the
rewritten closing sentence of Reading (C) (Dice of the five label-free-cut sessions, recomputed from the raw table) and
for the numbers in the new "Open decisions" section.
Writes ONLY verification/rerun_2026-10-07/compare_readme_v2.csv and compare_readme_v2_failures.csv.

A printed value passes when it equals the recomputed value after rounding to the printed number of decimals
(tolerance 0.5 * 10^-decimals; for e-notation, half a unit in the last printed significant digit).

    python3 check_readme_v2.py
"""
import json
import math
import os
import re

import numpy as np
import pandas as pd

VAL = "/path/to/Penn_LEAD/LEADchi/results_posthoc/validation_local"
OUT = f"{VAL}/verification/rerun_2026-10-07"
MINE = f"{VAL}/verification/mine"
OLD = open(f"{VAL}/archive_v1_prefix_2026-10-07/README.md", encoding="utf-8").read().split("\n")
NEW = open(f"{VAL}/README.md", encoding="utf-8").read().split("\n")
ANCHOR = {  # pre-fix line number -> unique substring of the revised line (edited lines only)
    56: "3. Because the design matrix is shared", 57: "4. **R2 cancels in group contrasts.**",
    58: "5. **The whole-brain chi reference shifts both components equally**",
    67: "| SN-VTA | chi_para | 4.49", 68: "| SN-VTA | chi_dia | 0.48", 71: "| SN-only | chi_para | 5.71",
    72: "| SN-only | chi_dia | 0.17", 75: "| GP | chi_para | 1.58", 76: "| GP | chi_dia | 0.70",
    78: "**BH-FDR within the X2 block.**", 112: "Under the model's assumptions, the point estimates",
    178: "Adding all four nuisance covariates together", 216: "**Reading (C).**", 224: "- All usable sessions: `nlme::lme`",
    245: "| chi_PUT | all usable, lme | + sex |", 246: "| chi_PUT | all usable, lme | + sex + cell |",
    249: "| **r2s_GP** (s^-1) | first-usable | age only |", 253: "| r2s_GP | all usable, lme | + sex |",
    278: "- Putamen chi_dia at Dr = 274, first-usable", 279: "- GP chi_dia (age only)",
    280: "- Putamen chi_para (first-usable, age only)", 302: "Our putamen R2* rise of",
    306: "- **The magnitude signal develops in the same sessions.**", 307: "- **Putamen chi does not follow.**",
    308: "- **In GP, chi tracks R2* as expected**", 149: "Dice is not related to R2* (+0.272 s^-1 per SD, p .3624)", 309: "- **We found no sign that referencing drives it.**",
    312: "  - (b) A real concurrent diamagnetic increase", 313: "- SN-VTA co-report:",
}
NEWNO = {}


def _map(old_no):
    if old_no in ANCHOR:
        hits = [i for i, t in enumerate(NEW) if ANCHOR[old_no] in t]
    else:
        hits = [i for i, t in enumerate(NEW) if t == OLD[old_no - 1]]
    assert len(hits) == 1, (old_no, len(hits))
    return hits[0] + 1


class _Lines(dict):
    def __missing__(self, k):
        n = _map(k)
        NEWNO[k] = n
        self[k] = NEW[n - 1]
        return self[k]


L = _Lines()                                       # L[pre-fix line number] -> text of the matching revised line


def find(anchor):
    hits = [t for t in NEW if anchor in t]
    assert len(hits) == 1, (anchor, len(hits))
    return hits[0]

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
    q1 = [r.q_BH_ext18_PSonly] if "(q" in c[7] else []
    q2 = [r.q_BH_ext18_PSADHD] if "(q" in c[8] else []
    assert bool(q) == bool(q1) == bool(q2), ln
    check(f"L{ln}", nums(c[3]) + nums(c[4]) + nums(c[5]) + nums(c[6]) + nums(c[7]) + nums(c[8]),
          bci(r) + [r.se_hc3, r.p_two] + q + s1(r) + q1 + s2(r) + q2, c[1] + " " + c[2])
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
assert "The numerical deviation is below 1e-13 in every row (floating-point noise)" in L[56]
ROWS.append(dict(where="L56", what="X2 identity max error (floating-point noise)", printed="< 1e-13",
                 recomputed=SC["identity_max_abs_err"], status="ok (noise; both < 1e-13)" if SC["identity_max_abs_err"] < 1e-13 else "MISMATCH"))
assert "beta_PS and SE_HC3 differ by less than 1e-13 (floating-point noise)" in L[57]
ROWS.append(dict(where="L57", what="R2 cancellation ranges of beta / SE (floating-point noise)", printed="< 1e-13",
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
_sp, _sd, _op, _od = A274("SNVTA", "chi_para"), A274("SNVTA", "chi_dia"), A274("SNonly", "chi_para"), A274("SNonly", "chi_dia")
text(112, "chi_para +4.49 ppb [0.23, 8.74] (p .039, q .118) and chi_dia +0.48 [-2.39, 3.35] (p .74, q .889)",
     bci(_sp) + [_sp.p_two, _sp.q_BH_family6] + bci(_sd) + [_sd.p_two, _sd.q_BH_family6])
text(112, "SN-only gives +5.71 (p .016, q .096) and +0.17 (p .91, q .913)",
     [_op.beta_PS, _op.p_two, _op.q_BH_family6, _od.beta_PS, _od.p_two, _od.q_BH_family6])
text(112, "(smallest q .096)", [prim.q_BH_family6.min()])
text(112, "(chi_para +1.58 [-3.09, 6.24], chi_dia +0.70 [-1.97, 3.37])", bci(A274("GP", "chi_para")) + bci(A274("GP", "chi_dia")))

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
text(149, "The R2* estimate moves little (1.10 to 1.01, 8 % lower) and its SE grows by 7 %; together these move p across .05 (.0350 to .0684)",
     [r2m.beta_PS, rd.beta_PS, 100 * (1 - rd.beta_PS / r2m.beta_PS), 100 * (rd.se_hc3 / r2m.se_hc3 - 1), 0.05, r2m.p_two, rd.p_two], "E1 reworded")
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
text(178, "(+10.4, p .042, n 53", [qmr.beta_PS, qmr.p_two, qmr.n])
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
text(216, "Their mean Dice is 0.771 versus 0.811 (p .019)", [row(QC, metric="dice_min").tail_mean, row(QC, metric="dice_min").TDrest_median,
                                                             row(QC, metric="dice_min").wilcox_p_tail_vs_TDrest],
     "0.811 = mean (0.8112) and also median (0.8107) of the other 18 TD sessions")
text(216, "(89th-100th TD percentile)", [sorted(CS.pctTDrest_motion_gsr)[1], CS.pctTDrest_motion_gsr.max()])
text(216, "(0th and 6th TD percentile)", [sorted(CS.pctTDrest_r2s_SNVTA)[0], sorted(CS.pctTDrest_r2s_SNVTA)[1]])
text(216, "lowers beta_PS from 4.96 to 3.33 ppb, and PS-only - TD from 8.16 (p .044) to 4.36 (p .276). The label-free cut lowers beta_PS to 1.45",
     [CR.iloc[0].beta_PS, CR.iloc[1].beta_PS, CR.iloc[0].PSonly_minus_TD, CR.iloc[0].PSonly_minus_TD_p, CR.iloc[1].PSonly_minus_TD,
      CR.iloc[1].PSonly_minus_TD_p, CR.iloc[2].beta_PS])
text(216, "(p_B ≤ .0022)", [CR.pB_fixed_bound.max()])
_f1, _f2 = 1 - CR.iloc[1].beta_PS / CR.iloc[0].beta_PS, 1 - CR.iloc[2].beta_PS / CR.iloc[0].beta_PS
assert "about a third of the positive point estimate (about 70 % with the label-free cut)" in L[216]
ROWS.append(dict(where="L216", what="'about a third ... (about 70 % with the label-free cut)'", printed="about 1/3 ; about 70 %",
                 recomputed=f"{_f1:.3f} ; {_f2:.3f}", status="ok" if (0.28 < _f1 < 0.38 and 0.65 < _f2 < 0.75) else "MISMATCH"))

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
dps = row(DR, outcome="chi_dia_PUT", sample=SAMP["first-usable"], covariates="+ sex")
dpc = row(DR, outcome="chi_dia_PUT", sample=SAMP["first-usable"], covariates="+ sex + cell")
dpls = row(DR, outcome="chi_dia_PUT", sample=SAMP["all usable, lme"], covariates="+ sex")
dpw = row(DR, outcome="chi_dia_PUT", sample=SAMP["lme, within-person"], covariates="+ sex + cell")
text(278, "-0.358 [-0.715, -0.001] ppb/y, p .0492 (age only); -0.352, p .0512 (+ sex); -0.346, p .0551 (+ sex + cell)",
     [dp0.slope_per_y, dp0.ci_lo, dp0.ci_hi, dp0.p_two, dps.slope_per_y, dps.p_two, dpc.slope_per_y, dpc.p_two], "chi_dia_PUT first-usable")
ln_lme = find("- Putamen chi_dia at Dr = 274, all usable (lme)")
for snip, vals in (("-0.501 [-0.840, -0.161], p .0055 (age only); -0.481, p .0076 (+ sex); -0.515 [-0.869, -0.160], p .0061",
                    [dpl.slope_per_y, dpl.ci_lo, dpl.ci_hi, dpl.p_two, dpls.slope_per_y, dpls.p_two, dplc.slope_per_y, dplc.ci_lo, dplc.ci_hi, dplc.p_two]),
                   ("-0.589 [-1.570, 0.393], p .2290", [dpw.slope_per_y, dpw.ci_lo, dpw.ci_hi, dpw.p_two])):
    assert snip in ln_lme, snip
    check("D.2 lme shortfall line", nums(snip), vals, snip)
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
text(302, "Our putamen R2* rise of 0.30 s^-1/y (0.24 SD/y) over 8.8-17.8 y", [rp0.slope_per_y, rp0.slope_per_y_in_SD] + ar["age_range_first"])
text(306, "Putamen R2* rises by 0.30 s^-1/y (p ≈ 1e-6), including within persons (0.30, p .007)", [rp0.slope_per_y, rp0.p_two, rpw.slope_per_y, rpw.p_two])
text(306, "Its standardised rate (0.24 SD/y) is similar to that of GP R2* (0.24 SD/y)", [rp0.slope_per_y_in_SD, rg0.slope_per_y_in_SD])
text(306, "(8.8-17.8 / 19.5 y)", ar["age_range_first"] + [ar["age_range_usable"][1]])
pl = [row(DR, outcome="chi_PUT", sample=SAMP["all usable, lme"], covariates=cv).slope_per_y for cv in ("age only", "+ sex", "+ sex + cell")]
text(307, "Its slope is 0.38 ppb/y (first-usable) and 0.08-0.14 ppb/y (all usable)", [put0.slope_per_y, min(pl), max(pl)])
pi0 = row(IM, roi="PUT", sample=SAMP["first-usable"], covariates="age only")
text(307, "below the 1.10 ppb/y implied by its own R2* slope", [pi0.implied_chi_slope_ppb_per_y])
SH = pd.read_csv(f"{MINE}/D_putamen_shortfall_by_Dr.csv")
_fu = SH[SH["sample"].str.startswith("first")]
_lm = SH[~SH["sample"].str.startswith("first")]
_p = lambda d, Dr: d[d.Dr == Dr].p_two  # noqa: E731
text(307, "(p .049-.055 at Dr = 274, depending on covariates; p ≥ .11 at Dr = 321 and 548)",
     [_p(_fu, 274.0).min(), _p(_fu, 274.0).max(), 274, min(_p(_fu, 321.0).min(), _p(_fu, 548.0).min()), 321, 548])
text(307, "(p .0055-.0076 at Dr = 274; .014-.020 at Dr = 321; .14-.19 at Dr = 548)",
     [_p(_lm, 274.0).min(), _p(_lm, 274.0).max(), 274, _p(_lm, 321.0).min(), _p(_lm, 321.0).max(), 321, _p(_lm, 548.0).min(), _p(_lm, 548.0).max(), 548])
gi0 = row(IM, roi="GP", sample=SAMP["first-usable"], covariates="age only")
text(308, "observed 2.40 vs implied 2.51 ppb/y", [gi0.slope_per_y_chi, gi0.implied_chi_slope_ppb_per_y])
text(308, "(r = .68 vs .21 in putamen)", [ar["r_chi_r2s_GP_first"], ar["r_chi_r2s_first"]])
text(312, "(chi_para up 0.74, chi_dia down 0.36 ppb/y)", [pp0.slope_per_y, -dp0.slope_per_y])
_r = lambda o, s, cv: row(DR, outcome=o, sample=SAMP[s], covariates=cv)  # noqa: E731
r0, rl, rc = _r("r2s_SNVTA", "first-usable", "age only"), _r("r2s_SNVTA", "all usable, lme", "age only"), _r("r2s_SNVTA", "first-usable", "+ sex + cell")
text(313, "R2* rises with age (0.31 [0.05, 0.57] s^-1/y first-usable, p .019; 0.50 [0.28, 0.71] lme; first-usable with + sex + cell 0.25, p .053)",
     [r0.slope_per_y, r0.ci_lo, r0.ci_hi, r0.p_two, rl.slope_per_y, rl.ci_lo, rl.ci_hi, rc.slope_per_y, rc.p_two])
c0, cl, cc, cw = (_r("chi_SNVTA", "first-usable", "age only"), _r("chi_SNVTA", "all usable, lme", "age only"),
                  _r("chi_SNVTA", "first-usable", "+ sex + cell"), _r("chi_SNVTA", "lme, within-person", "+ sex + cell"))
text(313, "1.31 [-0.07, 2.69] ppb/y first-usable (p .063); 1.24 [0.21, 2.27] lme; first-usable with + sex + cell 1.02 [-0.70, 2.75] (p .24); within-person 0.33 [-1.69, 2.34]",
     [c0.slope_per_y, c0.ci_lo, c0.ci_hi, c0.p_two, cl.slope_per_y, cl.ci_lo, cl.ci_hi, cc.slope_per_y, cc.ci_lo, cc.ci_hi, cc.p_two,
      cw.slope_per_y, cw.ci_lo, cw.ci_hi])
text(309, "(p .55 and .99;", [_r("ref_ppb", "first-usable", "age only").p_two, _r("ref2_minus_wb_ppb", "first-usable", "age only").p_two])

# ================================================================ NEW sentences and tables of the revised README
def textn(anchor, snippet, values, what=""):
    t = find(anchor)
    assert snippet in t, (anchor, snippet)
    check("new: " + anchor[:40], nums(snippet), values, what or snippet)


# A.2 BH notes
_p18 = A274("SNVTA", "chi_para")
textn("**BH-FDR within the X2 block.**", "in it, SN-VTA chi_para beta_PS has q .177", [_p18.q_BH_ext18_beta])
T4 = pd.read_csv("/path/to/Penn_LEAD/LEADchi/results_exploratory/tier4_exploratory/fdr_table.csv")
_rx = T4[T4.family == "registered X"].p_two.to_numpy(float)


def _bh(p):
    p = np.asarray(p, float); m = len(p); o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    out = np.empty(m); out[o] = np.minimum(q, 1); return out


_p6 = prim.p_two.to_numpy(float)
_p18all = np.r_[_p6, prim.PSonly_minus_TD_p.to_numpy(float), prim.PSADHD_minus_ADHD_p.to_numpy(float)]
textn("**The BH family is a choice (post hoc note).**", "the 6 X2 tests join the 5 registered tests", [len(_p6), 2, len(_rx)])
textn("**The BH family is a choice (post hoc note).**", "all p ≥ .576). That family of 11 gives a smallest q of .177, and the corresponding family of 23 (18 + 5) gives .175",
      [_rx.min(), len(_p6) + len(_rx), _bh(np.r_[_p6, _rx])[:6].min(), len(_p18all) + len(_rx), len(_p18all), len(_rx), _bh(np.r_[_p18all, _rx])[:18].min()])
# A.3 note
_a548 = row(A, roi="SNVTA", component="chi_para", Dr=548.0, R2=0.0)
textn("The p values in this table are nominal", "(3.48 [-0.13, 7.10], p .0586)", bci(_a548) + [_a548.p_two])
_lo = {Dr: row(A, roi="SNVTA", component="chi_para", Dr=Dr, R2=0.0).ci_lo for Dr in (137.0, 274.0, 321.0, 548.0)}
ROWS.append(dict(where="A.3 note / A.4", what="SN-VTA chi_para CI excludes 0 only for Dr <= 321", printed="137/274/321 yes; 548 no",
                 recomputed=str({k: round(v, 3) for k, v in _lo.items()}),
                 status="ok" if (_lo[137.0] > 0 and _lo[274.0] > 0 and _lo[321.0] > 0 and _lo[548.0] < 0) else "MISMATCH"))
_ord = all(row(A, roi=ro, component="chi_para", Dr=Dr, R2=0.0).beta_PS > row(A, roi=ro, component="chi_dia", Dr=Dr, R2=0.0).beta_PS
           for ro in ("SNVTA", "SNonly") for Dr in (137.0, 274.0, 321.0, 548.0))
ROWS.append(dict(where="A.3 note / A.4", what="ordering chi_para > chi_dia at all four Dr (SN-VTA, SN-only)", printed="yes",
                 recomputed=str(_ord), status="ok" if _ord else "MISMATCH"))
# B.4 new sentence
_t1f, _t1r, _t3r = (row(TR, model="within T1_low : full M1"), row(TR, model="within T1_low : reduced (PS*ADHD + age + sex)"),
                    row(TR, model="within T3_high : reduced (PS*ADHD + age + sex)"))
textn("Adding all four nuisance covariates together", "(lowest tertile p_B .2157 full M1 and .2993 reduced; reduced high-tertile model .0965)",
      [_t1f.pB_fixed_bound, _t1r.pB_fixed_bound, _t3r.pB_fixed_bound])
_inc = all(row(TR, model=f"within {t} : full M1").ci_lo < 0 < row(TR, model=f"within {t} : full M1").ci_hi for t in ("T1_low", "T2_mid", "T3_high"))
ROWS.append(dict(where="B.4", what="all three within-tertile CIs include 0 (full M1)", printed="yes", recomputed=str(_inc), status="ok" if _inc else "MISMATCH"))
ROWS.append(dict(where="B.4", what="'adding it lowers the estimate by about 1 ppb, alone or jointly'", printed="about 1",
                 recomputed=f"alone {M1.beta_PS - sg['dice_min'].beta_PS:.3f}", status="ok" if 0.6 < M1.beta_PS - sg['dice_min'].beta_PS < 1.5 else "MISMATCH"))
# C new sentence
_pc = {c: float(CS[CS.code == c].pctTDrest_dice_min.iloc[0]) for c in CS.code}
textn("**Reading (C).**", "The fourth (10a248) is at the TD median (50th percentile)", [10, 248, _pc["10a248"]])
ROWS.append(dict(where="C", what="three of four have a lower Dice than every other TD session (pct 0)", printed="3",
                 recomputed=str(sum(v == 0 for v in _pc.values())), status="ok" if sum(v == 0 for v in _pc.values()) == 3 else "MISMATCH"))
# D.2 shortfall table
for Dr in (137.0, 274.0, 321.0, 548.0):
    _h = [x for x in NEW if x.startswith(f"| {Dr:g} | ")]
    assert len(_h) == 1, Dr
    t = _h[0]
    cs = [c.strip() for c in t.split("|")][2:8]
    vals = []
    for smp, cvs in (("first-usable (OLS, HC3)", ("age only", "+ sex", "+ sex + cell")),
                     ("all usable sessions (lme REML, random intercept)", ("age only", "+ sex", "+ sex + cell"))):
        for cv in cvs:
            rr = SH[(SH.Dr == Dr) & (SH["sample"] == smp) & (SH.covariates == cv)].iloc[0]
            vals += [rr.slope, rr.p_two]
    check(f"D.2 shortfall table Dr {Dr:g}", sum((nums(c) for c in cs), []), vals, f"shortfall Dr {Dr:g}")
textn("In the first-usable sample the estimate is borderline", "(p .049-.055 across covariate sets)", [_p(_fu, 274.0).min(), _p(_fu, 274.0).max()])
textn("In the first-usable sample the estimate is borderline", "(p ≤ .020) but not at Dr = 548 (p .14-.19)",
      [max(_p(_lm, 137.0).max(), _p(_lm, 274.0).max(), _p(_lm, 321.0).max()), 548, _p(_lm, 548.0).min(), _p(_lm, 548.0).max()])
# Verification section
CSUM = pd.read_csv(f"{VAL}/verification/compare_csv_summary.csv")
CDIF = pd.read_csv(f"{VAL}/verification/compare_csv_discrepancies.csv")
CDIF["absd"] = (CDIF.theirs.astype(float) - CDIF.mine.astype(float)).abs()
OLDSUM = pd.read_csv(f"{VAL}/archive_v1_prefix_2026-10-07/verification_outputs/compare_csv_summary.csv")
_d = CSUM.set_index("file")
_isp = CDIF.column.isin(["p_two", "p_one_sided_pos"])
textn("- D mixed-model rows: the largest difference", "fell from 6.2e-3 relative",
      [OLDSUM.set_index("file").loc["D_age_slopes", "max_rel_diff"]])
textn("- D mixed-model rows: the largest difference", "to 3.3e-5 relative; no slope or CI bound differs by more than 1.2e-5 and no p value by more than 5.1e-6",
      [_d.loc[["D_age_slopes", "D_implied_chi_slope_from_R2s"], "max_rel_diff"].max(), CDIF[~_isp].absd.max(), CDIF[_isp].absd.max()])
textn("- D mixed-model rows: the largest difference", "174 cells (144 in",
      [int(_d.loc["D_age_slopes", "cells_outside_tol"] + _d.loc["D_implied_chi_slope_from_R2s", "cells_outside_tol"]), int(_d.loc["D_age_slopes", "cells_outside_tol"])])
textn("- D mixed-model rows: the largest difference", "30 derived ones", [int(_d.loc["D_implied_chi_slope_from_R2s", "cells_outside_tol"])])
_ac = _d.drop(index=["D_age_slopes", "D_implied_chi_slope_from_R2s"]).cells_outside_tol.sum()
ROWS.append(dict(where="Verification", what="all A-C cells agree (cells outside 1e-6)", printed="0", recomputed=str(int(_ac)), status="ok" if _ac == 0 else "MISMATCH"))

# Revision paragraph: size of the change in the lme rows (revised CSV vs the archived pre-fix CSV; not the verifier)
_new = pd.read_csv(f"{VAL}/D_age_slopes.csv")
_old = pd.read_csv(f"{VAL}/archive_v1_prefix_2026-10-07/D_age_slopes.csv")
_m = _new.merge(_old, on=["outcome", "sample", "covariates"], suffixes=("", "_old"))
_l = _m[~_m["sample"].str.startswith("first")]
_f = _m[_m["sample"].str.startswith("first")]
textn("**Revision 2026-10-07", "the 90 mixed-model rows were refitted with the nlme default optimiser: lme slopes change by at most 2.2e-4 ppb/y, SEs by at most 0.014 % and p values by at most 0.0005",
      [len(_l), (_l.slope_per_y - _l.slope_per_y_old).abs().max(), 100 * ((_l.se - _l.se_old).abs() / _l.se_old).max(), (_l.p_two - _l.p_two_old).abs().max()])
_fz = max(float((_f[c] - _f[c + "_old"]).abs().max()) for c in ("slope_per_y", "se", "ci_lo", "ci_hi", "p_two"))
ROWS.append(dict(where="Revision paragraph", what="first-usable OLS rows unchanged", printed="unchanged", recomputed=f"max abs diff {_fz:.1e}",
                 status="ok" if _fz == 0 else "MISMATCH"))

# ================================================================ Round 2 (2026-10-07) additions
# A.1 declaration 5: unreferenced p values labelled nominal (no q in any family) + q of the matching referenced tests
assert "These p values are nominal: the unreferenced rows are outside every BH family." in L[58]
_unref = A[A.component.astype(str).str.endswith("_unref")]
_noq = bool(_unref.filter(like="q_BH").isna().all().all()) and len(_unref) == 6
ROWS.append(dict(where="L58", what="unreferenced rows carry no BH q (outside every family)", printed="nominal; no q",
                 recomputed=f"{len(_unref)} rows, all q NaN: {_noq}", status="ok" if _noq else "MISMATCH"))
text(58, "q .118 (chi_para) and .889 (chi_dia) in the family of 6",
     [A274("SNVTA", "chi_para").q_BH_family6, A274("SNVTA", "chi_dia").q_BH_family6, len(_p6)])
# Reading (C), closing sentence: Dice of the sessions removed by the label-free cut, recomputed here from the raw table
_raw = pd.read_csv("/path/to/Penn_LEAD/LEADchi/mirror/unlocked/merged_qsmxt_primary_SNVTA.csv")
_raw["cell"] = np.where(_raw.study_group == "PRO/CHR", np.where(_raw.dx_adhd == 1, "PS+ADHD", "PS-only"),
                        np.where(_raw.dx_adhd == 1, "ADHD-only", "TD"))
_fu80 = _raw[_raw.usable == 1].sort_values(["participant_id", "ses_order"]).groupby("participant_id", sort=True).head(1)
assert len(_fu80) == 80 and set(_fu80.participant_id + "_" + _fu80.session_id) == set(
    (_raw.participant_id + "_" + _raw.session_id)[_raw.first_usable == 1])
_lf = _fu80[_fu80.chi_SNVTA < 35]
assert sorted(_lf.cell) == ["ADHD-only", "TD", "TD", "TD", "TD"]
_crlf = row(CR, model="M1 without every first-usable session with chi < 35 ppb (label-free; 5 sessions)")
assert _crlf.dropped_cells == "4/1/0/0" and abs(_crlf.dropped_min_chi - _lf.chi_SNVTA.min()) < 1e-6
_lf_adhd = _lf[_lf.cell == "ADHD-only"].iloc[0]
_td3 = set(_fu80[_fu80.cell == "TD"].nsmallest(3, "dice_min").chi_SNVTA.round(6))
_tail_low = set(_lf[(_lf.cell == "TD")].nsmallest(3, "dice_min").chi_SNVTA.round(6))
ROWS.append(dict(where="L216", what="three of the low-chi sessions have the three lowest TD Dice values", printed="3",
                 recomputed=str(_td3 == _tail_low), status="ok" if _td3 == _tail_low else "MISMATCH"))
text(216, "but 10a248 (0.814) and the ADHD-only session removed by the label-free cut (0.810) are close to the median Dice of all 80 sessions (0.814)",
     [10, 248, float(CS[CS.code == "10a248"].dice_min.iloc[0]), float(_lf_adhd.dice_min), len(_fu80), row(QC, metric="dice_min").all80_median])
assert "low-Dice sessions," not in L[216] and "low-chi, low-Dice" not in L[216]
# Open decisions section
textn("1. **nlme convergence tolerance (D).**", "by up to 3.3e-5 relative (174 cells above the strict 1e-6 test)",
      [_d.loc[["D_age_slopes", "D_implied_chi_slope_from_R2s"], "max_rel_diff"].max(),
       int(_d.loc["D_age_slopes", "cells_outside_tol"] + _d.loc["D_implied_chi_slope_from_R2s", "cells_outside_tol"]), 1e-6])
textn("3. **BH family for X2 (A.2).**", "(smallest q .096 / .137 / .177 / .175)",
      [prim.q_BH_family6.min(), ext_min, _bh(np.r_[_p6, _rx])[:6].min(), _bh(np.r_[_p18all, _rx])[:18].min()])

RES = pd.DataFrame(ROWS)
RES["readme_line_prefix"] = RES["where"]
RES["readme_line_revised"] = [f"L{NEWNO[int(w[1:])]}" if (w.startswith("L") and w[1:].isdigit() and int(w[1:]) in NEWNO) else "" for w in RES["where"]]
RES.to_csv(f"{OUT}/compare_readme_v2.csv", index=False)
bad = RES[~RES.status.str.startswith("ok")]
bad.to_csv(f"{OUT}/compare_readme_v2_failures.csv", index=False)
print(f"README numbers checked: {len(RES)}; not ok: {len(bad)}")
pd.set_option("display.width", 250)
pd.set_option("display.max_colwidth", 90)
print(bad.to_string())
