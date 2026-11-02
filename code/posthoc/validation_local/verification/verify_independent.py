#!/usr/bin/env python3
"""LEAD-chi -- INDEPENDENT VERIFICATION (POST HOC) of results_posthoc/validation_local (analyses A-D).

Language: Python (the code under check is R). Written from scratch; it does not source the frozen R helpers or the
validation R script. The model, samples and coding are re-implemented from their written definitions:
  sample   first-usable = earliest usable session per participant (checked against the frozen first_usable flag);
           all usable = usable == 1; plausibility rules re-applied (envelopes; none fires in these rows, checked);
  coding   PS_e = +0.5 PRO/CHR else -0.5; ADHD_e = +0.5 if dx_adhd == 1; sex_e = +0.5 male;
           race3 = Black / White / Other (all other categories);
  M1       y ~ PS_e*ADHD_e + ns(age, df = 3) + sex_e + race3 + z(roi_vol) + z(Q_phase), OLS, HC3, t(df = n - p).
Reads (read-only): the merged table, gates_preunlock.json, results_confirmatory/results.json and the CSV / JSON
outputs under check. Writes ONLY into results_posthoc/validation_local/verification/.

    python3 verify_independent.py
"""
import json
import math
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import vlib as V  # noqa: E402

ROOT = "/path/to/Penn_LEAD/LEADchi"
INPUT = f"{ROOT}/mirror/unlocked/merged_qsmxt_primary_SNVTA.csv"
GATES = f"{ROOT}/mirror/gates/gates_preunlock.json"
FROZEN = f"{ROOT}/results_confirmatory/results.json"
VAL = f"{ROOT}/results_posthoc/validation_local"
OUT = f"{VAL}/verification"
assert os.path.realpath(OUT).startswith(os.path.realpath(VAL) + "/verification")

LOG = []


def log(*a):
    s = " ".join(str(x) for x in a)
    LOG.append(s)
    print(s)


# ================================================================================================================
# 0. data
# ================================================================================================================
raw = pd.read_csv(INPUT)
d = raw.copy()
# plausibility re-application (sec. 12.4 rules that touch any column used below)
ENV = dict(chi_GP=(-50, 300), chi_PUT=(-80, 200), chi_SNVTA=(-50, 300))
plaus_changes = {}
for c, (lo, hi) in ENV.items():
    bad = d[c].notna() & ((d[c] < lo) | (d[c] > hi))
    if c == "chi_SNVTA":
        assert not (bad & (d.usable == 1)).any(), "SN-VTA envelope failure in a usable row"
    else:
        plaus_changes[c] = int((bad & (d.usable == 1)).sum())
        d.loc[bad, c] = np.nan
bad = d.r2s_SNVTA.notna() & ((d.r2s_SNVTA < 5) | (d.r2s_SNVTA > 100))
plaus_changes["r2s_SNVTA"] = int((bad & (d.usable == 1)).sum())
d.loc[bad, "r2s_SNVTA"] = np.nan
assert not ((d.fimp_SNVTA > 0.05) & (d.usable == 1)).any()
assert not ((d.inmask_SNVTA < 0.95) & (d.usable == 1)).any()
log("plausibility edits in usable rows:", plaus_changes)

d["PS_e"] = np.where(d.study_group == "PRO/CHR", 0.5, -0.5)
d["ADHD_e"] = np.where(d.dx_adhd == 1, 0.5, -0.5)
d["sex_e"] = np.where(d.sex == "M", 0.5, -0.5)
d["race3"] = np.where(d.race == "Black or African American", "Black", np.where(d.race == "White", "White", "Other"))
d["cell"] = np.where(d.PS_e > 0, np.where(d.ADHD_e > 0, "PS+ADHD", "PS-only"), np.where(d.ADHD_e > 0, "ADHD-only", "TD"))
d["sess_key"] = d.participant_id + "_" + d.session_id
CELLS = ["TD", "ADHD-only", "PS-only", "PS+ADHD"]

usable = d[d.usable == 1].sort_values(["participant_id", "ses_order"]).reset_index(drop=True)
first = usable.groupby("participant_id", sort=True).head(1).reset_index(drop=True)
assert set(first.sess_key) == set(d.sess_key[d.first_usable == 1]), "first-usable flag mismatch"
cells_first = first.cell.value_counts().reindex(CELLS).tolist()
log("first-usable n =", len(first), "cells", cells_first, "; usable sessions", len(usable), "subjects",
    usable.participant_id.nunique(), "; subjects with >= 2 usable:", int((usable.participant_id.value_counts() >= 2).sum()))

gates = json.load(open(GATES))
rL = float(gates["rL_primary"])
frozen = json.load(open(FROZEN))


def cellstr(dd):
    return "/".join(str(int((dd.cell == c).sum())) for c in CELLS)


def fit_m1(dd, y, extra_z=(), extra_raw_fn=None, need=()):
    keep = dd[y].notna()
    for v in list(extra_z) + list(need):
        keep &= dd[v].notna()
    dd = dd[keep].reset_index(drop=True)
    extra_raw = extra_raw_fn(dd) if extra_raw_fn else None
    X, names, yy = V.m1_design(dd, y, extra_z, extra_raw)
    m = V.OLS(X, yy, names)
    m.data, m.yname = dd, y
    return m


SIMPLE1 = {"PS_e": 1, "PS_e:ADHD_e": -0.5}     # PS-only - TD
SIMPLE2 = {"PS_e": 1, "PS_e:ADHD_e": 0.5}      # PS+ADHD - ADHD-only


def huber_scale(m):
    return V.rlm_huber(m.X, m.y)


def ps_row(m, model, bound=None):
    t = m.coef("PS_e")
    s1, s2 = m.lin(SIMPLE1), m.lin(SIMPLE2)
    r = dict(model=model, outcome=m.yname, n=m.n, cells=cellstr(m.data), beta_PS=t["est"], se_hc3=t["se"],
             ci_lo=t["lo"], ci_hi=t["hi"], t=t["t"], df=t["df"], p_two=t["p"],
             PSonly_minus_TD=s1["est"], PSonly_minus_TD_se=s1["se"], PSonly_minus_TD_lo=s1["lo"],
             PSonly_minus_TD_hi=s1["hi"], PSonly_minus_TD_p=s1["p"],
             PSADHD_minus_ADHD=s2["est"], PSADHD_minus_ADHD_se=s2["se"], PSADHD_minus_ADHD_lo=s2["lo"],
             PSADHD_minus_ADHD_hi=s2["hi"], PSADHD_minus_ADHD_p=s2["p"], sigma_resid=m.sigma)
    if bound is not None:
        tB, pB = V.bound_test(t["est"], t["se"], t["df"], bound)
        r.update(tB_fixed_bound=tB, pB_fixed_bound=pB)
    return r


# ================================================================================================================
# reproduction of the frozen M1 and outcome-B bound
# ================================================================================================================
M1 = fit_m1(first, "chi_SNVTA")
hb = huber_scale(M1)
sigma_B = min(M1.sigma, hb["scale"])
dB = 0.6 * math.sqrt(min(1.0, rL))
BOUND = -dB * sigma_B
t1 = M1.coef("PS_e")
tB0, pB0 = V.bound_test(t1["est"], t1["se"], t1["df"], BOUND)
fd = frozen["H1"]["decision"]
REPRO = dict(beta=(t1["est"], fd["beta_ppb"]), se=(t1["se"], fd["se_hc3"]), p=(t1["p"], fd["p_two"]),
             sigma_resid=(M1.sigma, fd["sigma_resid"]), sigma_rob=(hb["scale"], fd["sigma_rob"]),
             bound=(BOUND, fd["B_registered"]["bound_ppb"]), tB=(tB0, fd["B_registered"]["tB"]),
             pB=(pB0, fd["B_registered"]["pB"]), huber_beta=(float(hb["coef"][M1.idx("PS_e")]), frozen["H1"]["huber"]["est"]))
for k, (a, b) in REPRO.items():
    log(f"repro {k:12s} mine {a:.12g} frozen {b:.12g} rel {abs(a - b) / max(abs(b), 1e-12):.2e}")


def frozen_style_pB(m):
    h = huber_scale(m)
    sB = min(m.sigma, h["scale"])
    t = m.coef("PS_e")
    return V.bound_test(t["est"], t["se"], t["df"], -dB * sB)[1]


# ================================================================================================================
# A. X2
# ================================================================================================================
DR0 = 274.0
DR_SET = {"Dr x 0.5 (137)": 137.0, "Dr x 1 (274, paper)": 274.0, "theory (321)": 321.0, "Dr x 2 (548)": 548.0}
ROIS = {"SNVTA": ("chi_SNVTA", "r2s_SNVTA"), "SNonly": ("chi_SNonly", "r2s_SNonly"), "GP": ("chi_GP", "r2s_GP")}


def sep(dd, chi, r2s, Dr, R2=0.0):
    q = 1000.0 * (dd[r2s] - R2) / Dr
    return (dd[chi] + q) / 2.0, (dd[chi] - q) / 2.0


A_rows, A_desc = [], []
for roi, (cv, rv) in ROIS.items():
    base = first[first[cv].notna() & first[rv].notna()].reset_index(drop=True)
    m_chi, m_r2 = fit_m1(base, cv), fit_m1(base, rv)
    for m, lab in ((m_chi, "total chi (input)"), (m_r2, "R2* (input, s^-1)")):
        r = ps_row(m, f"{roi}: {lab}")
        r.update(block="A_X2_inputs")
        A_rows.append(r)
    for nm, Dr in DR_SET.items():
        for R2 in (0.0, 10.0, 15.0):
            if R2 != 0 and Dr != DR0:
                continue
            dd = base.copy()
            dd["chi_para"], dd["chi_dia"] = sep(dd, cv, rv, Dr, R2)
            for comp in ("chi_para", "chi_dia"):
                m = fit_m1(dd, comp)
                r = ps_row(m, f"{roi}: {comp}, {nm}, R2 = {R2:g}")
                sg = 1 if comp == "chi_para" else -1
                r.update(block="A_X2", roi=roi, component=comp, Dr=Dr, R2=R2,
                         beta_from_identity=(m_chi.coef("PS_e")["est"] + sg * 1000 * m_r2.coef("PS_e")["est"] / Dr) / 2,
                         mean_all=dd[comp].mean(),
                         share_constraint_violated=float((dd[comp] < 0).mean() if comp == "chi_para" else (dd[comp] > 0).mean()))
                for c in CELLS:
                    r["mean_" + "".join(ch for ch in c if ch.isalpha())] = dd.loc[dd.cell == c, comp].mean()
                A_rows.append(r)
    dd = base.copy()
    dd["chi_unref"] = dd[cv] + dd["ref_ppb"]
    dd["chi_para_unref"], dd["chi_dia_unref"] = sep(dd, "chi_unref", rv, DR0)
    for comp in ("chi_para_unref", "chi_dia_unref"):
        r = ps_row(fit_m1(dd, comp), f"{roi}: {comp} (chi + ref_ppb), Dr 274")
        r.update(block="A_X2_reference", roi=roi, component=comp, Dr=DR0, R2=0.0)
        A_rows.append(r)
    for nm, Dr in DR_SET.items():
        pa, di = sep(base, cv, rv, Dr)
        A_desc.append(dict(roi=roi, Dr_label=nm, Dr=Dr, mean_chi=base[cv].mean(), mean_R2s=base[rv].mean(),
                           mean_R2s_over_Dr_ppb=(1000 * base[rv] / Dr).mean(), mean_para=pa.mean(), sd_para=pa.std(ddof=1),
                           mean_dia=di.mean(), sd_dia=di.std(ddof=1), n_para_lt0=int((pa < 0).sum()),
                           n_dia_gt0=int((di > 0).sum()), r_para_dia=float(np.corrcoef(pa, di)[0, 1])))
A = pd.DataFrame(A_rows)
A_desc = pd.DataFrame(A_desc)
ident_err = float(np.nanmax(np.abs(A.beta_PS - A.beta_from_identity)))
canc = A[(A.block == "A_X2") & (A.Dr == DR0)]
canc_b = float(canc.groupby(["roi", "component"]).beta_PS.agg(lambda b: b.max() - b.min()).max())
canc_se = float(canc.groupby(["roi", "component"]).se_hc3.agg(lambda b: b.max() - b.min()).max())
prim = (A.block == "A_X2") & (A.Dr == DR0) & (A.R2 == 0)
A["q_BH_family6"] = np.nan
A.loc[prim, "q_BH_family6"] = V.bh(A.loc[prim, "p_two"])
ext = np.concatenate([A.loc[prim, "p_two"], A.loc[prim, "PSonly_minus_TD_p"], A.loc[prim, "PSADHD_minus_ADHD_p"]])
eq = V.bh(ext)
k = int(prim.sum())
A["q_BH_ext18_beta"] = np.nan
A["q_BH_ext18_PSonly"] = np.nan
A["q_BH_ext18_PSADHD"] = np.nan
A.loc[prim, "q_BH_ext18_beta"] = eq[:k]
A.loc[prim, "q_BH_ext18_PSonly"] = eq[k:2 * k]
A.loc[prim, "q_BH_ext18_PSADHD"] = eq[2 * k:]
# counts of chi_dia > 0 at R2 = 10 (README A.1 decl. 4)
dd = first.copy()
_, dia10 = sep(dd, "chi_SNVTA", "r2s_SNVTA", DR0, 10.0)
n_dia_pos_R2_10 = int((dia10 > 0).sum())
log(f"A: identity max err {ident_err:.2e}; R2 cancellation beta range {canc_b:.2e}, SE range {canc_se:.2e}; "
    f"min q6 {A.q_BH_family6.min():.4f}; min q18 {np.nanmin(eq):.4f}; SN-VTA dia>0 at R2=10: {n_dia_pos_R2_10}/80")

# ================================================================================================================
# B. nuisance
# ================================================================================================================
NUIS = ["dice_min", "motion_gsr", "ref2_minus_wb_ppb", "icv_seg_ml"]
B_rows = []
r = ps_row(M1, "M1 (frozen primary)", bound=BOUND)
r["pB_frozen_style"] = frozen_style_pB(M1)
B_rows.append(r)
single = {}
for v in NUIS:
    m = fit_m1(first, "chi_SNVTA", extra_z=(v,))
    single[v] = m
    r = ps_row(m, f"M1 + {v} (alone)", bound=BOUND)
    r["pB_frozen_style"] = frozen_style_pB(m)
    B_rows.append(r)
joint = fit_m1(first, "chi_SNVTA", extra_z=tuple(NUIS))
r = ps_row(joint, "M1 + Dice + motion + CSF offset + ICV (joint)", bound=BOUND)
r["pB_frozen_style"] = frozen_style_pB(joint)
B_rows.append(r)
B = pd.DataFrame(B_rows)
JC = []
for term, nm in [("PS_e", "PS_e"), ("ADHD_e", "ADHD_e")] + [(f"z_{v}", f"scale({v})") for v in NUIS]:
    c = joint.coef(term)
    JC.append(dict(term=nm, est=c["est"], se_hc3=c["se"], t=c["t"], p_two=c["p"], ci_lo=c["lo"], ci_hi=c["hi"]))
JC = pd.DataFrame(JC)


def vif(m, term="PS_e"):
    j = m.idx(term)
    others = [i for i in range(m.p) if i != j]
    Xo = m.X[:, others]
    b = np.linalg.lstsq(Xo, m.X[:, j], rcond=None)[0]
    res = m.X[:, j] - Xo @ b
    tss = np.sum((m.X[:, j] - m.X[:, j].mean()) ** 2)
    return 1.0 / (np.sum(res ** 2) / tss)


vif_joint, vif_m1 = vif(joint), vif(M1)
# frozen sensitivity rows
fs = {row["sensitivity"]: row for row in frozen["sensitivity"]["table"]}
FROZEN_SENS = {"dice_min": fs["+ registration Dice (minimum of the G5 Dice values)"],
               "motion_gsr": fs["+ magnitude ghost-to-signal motion index"],
               "icv_seg_ml": fs["+ intracranial volume (SynthSeg, coded side)"]}
sens_cmp = []
for v, row in FROZEN_SENS.items():
    m = single[v]
    t = m.coef("PS_e")
    sens_cmp.append(dict(covariate=v, beta_mine=t["est"], beta_frozen=row["beta_ppb"], se_mine=t["se"],
                         se_frozen=row["se_hc3"], p_mine=t["p"], p_frozen=row["p_two"],
                         pB_frozen_style_mine=frozen_style_pB(m), pB_frozen=row["pB"]))
sens_cmp = pd.DataFrame(sens_cmp)
# R2*
R2rows = []
for lab, ez in (("R2* M1", ()), ("R2* M1 + Dice", ("dice_min",)), ("R2* M1 + joint nuisance", tuple(NUIS))):
    m = fit_m1(first, "r2s_SNVTA", extra_z=ez)
    r = ps_row(m, lab)
    if lab == "R2* M1 + Dice":
        c = m.coef("z_dice_min")
        r.update(dice_coef=c["est"], dice_coef_p=c["p"])
    R2rows.append(r)
R2T = pd.DataFrame(R2rows)
K2 = frozen["H1"]["K2"]
# Dice tertiles
cuts = V.q7(first.dice_min, [1 / 3, 2 / 3])
tert = np.where(first.dice_min <= cuts[0], "T1_low", np.where(first.dice_min <= cuts[1], "T2_mid", "T3_high"))
first["dice_t"] = tert
TT = []
for tl in ("T1_low", "T2_mid", "T3_high"):
    sub = first[first.dice_t == tl]
    TT.append(dict(tertile=tl, **{c: int((sub.cell == c).sum()) for c in CELLS},
                   dice_range=f"{sub.dice_min.min():.4f}-{sub.dice_min.max():.4f}", mean_chi=sub.chi_SNVTA.mean()))
TT = pd.DataFrame(TT)
tab_cell = TT[CELLS].to_numpy()
fisher_cell = V.fisher_rxc(tab_cell)
tab_ps = np.array([[int(((first.dice_t == tl) & (first.PS_e == s)).sum()) for s in (-0.5, 0.5)]
                   for tl in ("T1_low", "T2_mid", "T3_high")])
fisher_ps = V.fisher_rxc(tab_ps)
TR = []
for tl in ("T1_low", "T2_mid", "T3_high"):
    sub = first[first.dice_t == tl].reset_index(drop=True)
    m = fit_m1(sub, "chi_SNVTA")
    TR.append(ps_row(m, f"within {tl} : full M1", bound=BOUND))
    X = np.column_stack([np.ones(len(sub)), sub.PS_e, sub.ADHD_e, sub.age, sub.sex_e, sub.PS_e * sub.ADHD_e])
    mr = V.OLS(X, sub.chi_SNVTA.to_numpy(float), ["(Intercept)", "PS_e", "ADHD_e", "age", "sex_e", "PS_e:ADHD_e"])
    mr.data, mr.yname = sub, "chi_SNVTA"
    TR.append(ps_row(mr, f"within {tl} : reduced (PS*ADHD + age + sex)", bound=BOUND))


def tdum(dd):
    return {"t_T2": (dd.dice_t == "T2_mid").astype(float).to_numpy(), "t_T3": (dd.dice_t == "T3_high").astype(float).to_numpy()}


m = fit_m1(first, "chi_SNVTA", extra_raw_fn=tdum)
TR.append(ps_row(m, "pooled M1 + Dice-tertile fixed effects", bound=BOUND))
m = fit_m1(first, "chi_SNVTA", extra_z=("dice_min",), extra_raw_fn=tdum)
TR.append(ps_row(m, "pooled M1 + Dice tertile FE + linear Dice", bound=BOUND))
qm = first[first.dice_t != "T1_low"].reset_index(drop=True)
TR.append(ps_row(fit_m1(qm, "chi_SNVTA"), "quality-matched: drop lowest Dice tertile (M1)", bound=BOUND))
TR.append(ps_row(fit_m1(qm, "chi_SNVTA", extra_z=("dice_min",)), "quality-matched: drop lowest Dice tertile (M1 + Dice)", bound=BOUND))
TR = pd.DataFrame(TR)


def het_fn(dd):
    t2 = (dd.dice_t == "T2_mid").astype(float).to_numpy()
    t3 = (dd.dice_t == "T3_high").astype(float).to_numpy()
    return {"t_T2": t2, "t_T3": t3, "PSxT2": dd.PS_e.to_numpy() * t2, "PSxT3": dd.PS_e.to_numpy() * t3}


mh = fit_m1(first, "chi_SNVTA", extra_raw_fn=het_fn)
jj = [mh.idx("PSxT2"), mh.idx("PSxT3")]
bj = mh.b[jj]
W = float(bj @ np.linalg.solve(mh.V[np.ix_(jj, jj)], bj))
HET = dict(wald_chi2=W, df1=2, F=W / 2, df2=mh.df, p_F=float(stats.f.sf(W / 2, 2, mh.df)))
rho, rho_p = V.spearman(first.dice_min, first.chi_SNVTA)
log(f"B: joint beta {B.beta_PS.iloc[-1]:.6f}; het F {HET['F']:.4f} p {HET['p_F']:.4f}; fisher cell {fisher_cell:.4f} PS {fisher_ps:.4f}; "
    f"spearman {rho:.4f} p {rho_p:.5f}; VIF M1 {vif_m1:.4f} joint {vif_joint:.4f}")

# ================================================================================================================
# C. low TD tail
# ================================================================================================================
first["rstud_M1"] = M1.rstudent()     # M1.data row order == first row order (no NA dropped)
assert (M1.data.sess_key.to_numpy() == first.sess_key.to_numpy()).all()
tail = first[(first.cell == "TD") & (first.chi_SNVTA < 35)].sort_values("chi_SNVTA")
td_rest = first[(first.cell == "TD") & ~first.sess_key.isin(tail.sess_key)]
QCV = ["dice_min", "motion_gsr", "fimp_SNVTA", "inmask_SNVTA", "Q_phase", "roi_vol", "r2s_SNVTA", "ref2_minus_wb_ppb", "icv_seg_ml", "age"]
CS = tail[["code", "session_id", "ses_order", "chi_SNVTA", "qc_grade", "dice_min", "motion_gsr", "fimp_SNVTA", "inmask_SNVTA",
           "Q_phase", "r2s_SNVTA", "roi_vol", "chi_SNVTA_median", "chi_SNonly", "chi_VTAPBP", "ref2_minus_wb_ppb", "age",
           "rstud_M1"]].copy().reset_index(drop=True)
for v in QCV:
    CS["pctTDrest_" + v] = [100 * float(np.mean(td_rest[v].to_numpy() <= x)) for x in tail[v]]
    CS["pctAll80_" + v] = [100 * float(np.mean(first[v].to_numpy() <= x)) for x in tail[v]]
QC = []
for v in QCV:
    a, b = tail[v].to_numpy(float), td_rest[v].to_numpy(float)
    p, _, how = V.wilcox_rs(a, b)
    QC.append(dict(metric=v, tail_mean=a.mean(), tail_min=a.min(), tail_max=a.max(), TDrest_n=len(b), TDrest_median=np.median(b),
                   TDrest_q25=V.q7(b, 0.25), TDrest_q75=V.q7(b, 0.75), TDrest_min=b.min(), TDrest_max=b.max(),
                   all80_median=np.median(first[v]), wilcox_p_tail_vs_TDrest=p, wilcox_method=how))
QC = pd.DataFrame(QC)
qc_grade = dict(tail=tail.qc_grade.value_counts().to_dict(), TD_rest=td_rest.qc_grade.value_counts().to_dict(),
                all80=first.qc_grade.value_counts().to_dict())
CR = [ps_row(M1, "M1 (all 80)", bound=BOUND)]
CR[0]["huber_beta"] = float(hb["coef"][M1.idx("PS_e")])
m_nt = fit_m1(first[~first.sess_key.isin(tail.sess_key)].reset_index(drop=True), "chi_SNVTA")
r = ps_row(m_nt, "M1 without the 4 TD sessions with chi < 35 ppb", bound=BOUND)
r["huber_beta"] = float(V.rlm_huber(m_nt.X, m_nt.y)["coef"][m_nt.idx("PS_e")])
CR.append(r)
lf = first[first.chi_SNVTA < 35]
m_lf = fit_m1(first[first.chi_SNVTA >= 35].reset_index(drop=True), "chi_SNVTA")
r = ps_row(m_lf, f"M1 without every first-usable session with chi < 35 ppb (label-free; {len(lf)} sessions)", bound=BOUND)
r["huber_beta"] = float(V.rlm_huber(m_lf.X, m_lf.y)["coef"][m_lf.idx("PS_e")])
r["dropped_cells"] = cellstr(lf)
r["dropped_min_chi"] = float(lf.chi_SNVTA.min())
CR.append(r)
CR = pd.DataFrame(CR)
log(f"C: tail n {len(tail)}; refit beta {CR.beta_PS.tolist()}")

# ================================================================================================================
# D. putamen / GP age slopes
# ================================================================================================================


def add_derived(dd):
    dd = dd.copy()
    dd["chi_PUT_unref"] = dd.chi_PUT + dd.ref_ppb
    dd["chi_PUT_csfref"] = dd.chi_PUT - dd.ref2_minus_wb_ppb
    dd["chi_GP_unref"] = dd.chi_GP + dd.ref_ppb
    dd["chi_GP_csfref"] = dd.chi_GP - dd.ref2_minus_wb_ppb
    dd["chi_para_PUT"], dd["chi_dia_PUT"] = sep(dd, "chi_PUT", "r2s_PUT", DR0)
    dd["chi_para_GP"], dd["chi_dia_GP"] = sep(dd, "chi_GP", "r2s_GP", DR0)
    return dd


firstD, usableD = add_derived(first), add_derived(usable)
D_OUT = ["r2s_PUT", "chi_PUT", "chi_PUT_raw", "r2s_GP", "chi_GP", "chi_GP_raw", "chi_PUT_unref", "chi_PUT_csfref",
         "chi_GP_unref", "chi_GP_csfref", "ref_ppb", "ref2_minus_wb_ppb", "chi_para_PUT", "chi_dia_PUT", "chi_para_GP",
         "chi_dia_GP", "r2s_SNVTA", "chi_SNVTA"]
COVS = {"age only": [], "+ sex": ["sex_e"], "+ sex + cell": ["sex_e", "PS_e", "ADHD_e", "PS_e:ADHD_e"]}


def design_cols(dd, terms):
    cols = {"(Intercept)": np.ones(len(dd))}
    for t in terms:
        if t == "PS_e:ADHD_e":
            cols[t] = (dd.PS_e * dd.ADHD_e).to_numpy(float)
        else:
            cols[t] = dd[t].to_numpy(float)
    return np.column_stack(list(cols.values())), list(cols)


DR = []
for y in D_OUT:
    sdy = firstD[y].std(ddof=1)
    for cn, cv in COVS.items():
        dd = firstD[firstD[y].notna()]
        X, nm = design_cols(dd, ["age"] + cv)
        m = V.OLS(X, dd[y].to_numpy(float), nm)
        c = m.coef("age")
        DR.append(dict(outcome=y, sample="first-usable (OLS, HC3)", covariates=cn, n=m.n, n_subjects=dd.participant_id.nunique(),
                       slope_per_y=c["est"], se=c["se"], ci_lo=c["lo"], ci_hi=c["hi"], df=c["df"], p_two=c["p"],
                       p_one_sided_pos=c["p_upper"], slope_per_y_in_SD=c["est"] / sdy))
        du = usableD[usableD[y].notna()]
        X, nm = design_cols(du, ["age"] + cv)
        mm = V.LMMRI(X, du[y].to_numpy(float), du.participant_id.to_numpy(), nm)
        c = mm.coef("age")
        DR.append(dict(outcome=y, sample="all usable sessions (lme REML, random intercept)", covariates=cn, n=mm.N,
                       n_subjects=mm.G, slope_per_y=c["est"], se=c["se"], ci_lo=c["lo"], ci_hi=c["hi"], df=c["df"],
                       p_two=c["p"], p_one_sided_pos=c["p_upper"], slope_per_y_in_SD=c["est"] / sdy, theta=mm.theta))
    du = usableD[usableD[y].notna()].copy()
    du["age_b"] = du.groupby("participant_id").age.transform("mean")
    du["age_w"] = du.age - du.age_b
    X, nm = design_cols(du, ["age_w", "age_b", "sex_e", "PS_e", "ADHD_e", "PS_e:ADHD_e"])
    mm = V.LMMRI(X, du[y].to_numpy(float), du.participant_id.to_numpy(), nm)
    for tm in ("age_w", "age_b"):
        c = mm.coef(tm)
        DR.append(dict(outcome=y, sample="all usable (lme REML), " + ("WITHIN-person age" if tm == "age_w" else "BETWEEN-person mean age"),
                       covariates="+ sex + cell", n=mm.N, n_subjects=mm.G, slope_per_y=c["est"], se=c["se"], ci_lo=c["lo"],
                       ci_hi=c["hi"], df=c["df"], p_two=c["p"], p_one_sided_pos=c["p_upper"], slope_per_y_in_SD=c["est"] / sdy,
                       theta=mm.theta))
DRT = pd.DataFrame(DR)
EXPECT = {"chi_PUT": (2.6 / 1.95, 3.4 / 1.95), "chi_GP": (4.4 / 1.95, 6.9 / 1.95)}
DRT["expected_lo"] = np.nan
DRT["expected_hi"] = np.nan
for y, (lo, hi) in EXPECT.items():
    for yy in (y, y + "_raw", y + "_unref", y + "_csfref"):
        i = DRT.outcome == yy
        DRT.loc[i, "expected_lo"] = lo
        DRT.loc[i, "expected_hi"] = hi
DRT["ci_excludes_expected_range"] = np.where(DRT.expected_lo.isna(), np.nan,
                                             ((DRT.ci_hi < DRT.expected_lo) | (DRT.ci_lo > DRT.expected_hi)).astype(float))
IMPL = []
for roi in ("PUT", "GP"):
    a = DRT[DRT.outcome == "r2s_" + roi].set_index(["sample", "covariates"])
    b = DRT[DRT.outcome == "chi_" + roi].set_index(["sample", "covariates"])
    for key in a.index:
        IMPL.append(dict(roi=roi, sample=key[0], covariates=key[1], slope_per_y_r2s=a.loc[key, "slope_per_y"],
                         ci_lo_r2s=a.loc[key, "ci_lo"], ci_hi_r2s=a.loc[key, "ci_hi"], slope_per_y_chi=b.loc[key, "slope_per_y"],
                         ci_lo_chi=b.loc[key, "ci_lo"], ci_hi_chi=b.loc[key, "ci_hi"],
                         implied_chi_slope_ppb_per_y=1000 * a.loc[key, "slope_per_y"] / DR0,
                         implied_lo=1000 * a.loc[key, "ci_lo"] / DR0, implied_hi=1000 * a.loc[key, "ci_hi"] / DR0))
IMPL = pd.DataFrame(IMPL)
rp = lambda a, b: float(np.corrcoef(a, b)[0, 1])  # noqa: E731
res_age = lambda v: v - np.polyval(np.polyfit(firstD.age, v, 1), firstD.age)  # noqa: E731
CORR = dict(r_chi_r2s_first=rp(firstD.chi_PUT, firstD.r2s_PUT),
            r_chi_r2s_first_age_partial=rp(res_age(firstD.chi_PUT), res_age(firstD.r2s_PUT)),
            r_chi_r2s_GP_first=rp(firstD.chi_GP, firstD.r2s_GP),
            age_range_first=[float(firstD.age.min()), float(firstD.age.max())],
            age_range_usable=[float(usableD.age.min()), float(usableD.age.max())],
            n_subjects_with_2plus_usable=int((usableD.participant_id.value_counts() >= 2).sum()),
            raw_equal_PUT=bool(np.array_equal(usableD.chi_PUT.to_numpy(), usableD.chi_PUT_raw.to_numpy())),
            raw_equal_GP=bool(np.array_equal(usableD.chi_GP.to_numpy(), usableD.chi_GP_raw.to_numpy())))
log("D: corr", CORR)

# Dr-dependence of the putamen "shortfall" test (chi_dia_PUT age slope), first-usable OLS-HC3 and lme, age only
SHORT = []
for Dr in (137.0, 274.0, 321.0, 548.0):
    for lab, dd in (("first-usable (OLS, HC3)", firstD), ("all usable sessions (lme REML, random intercept)", usableD)):
        yv = ((dd.chi_PUT - 1000 * dd.r2s_PUT / Dr) / 2).to_numpy(float)
        for cn, cv in COVS.items():
            X, nm = design_cols(dd, ["age"] + cv)
            mdl = V.OLS(X, yv, nm) if lab.startswith("first") else V.LMMRI(X, yv, dd.participant_id.to_numpy(), nm)
            c = mdl.coef("age")
            SHORT.append(dict(Dr=Dr, sample=lab, covariates=cn, slope=c["est"], ci_lo=c["lo"], ci_hi=c["hi"], p_two=c["p"]))
SHORT = pd.DataFrame(SHORT)
log("putamen shortfall (chi_dia_PUT age slope) by Dr:\n" + SHORT.to_string())

# optional second check of the mixed models with statsmodels MixedLM (REML)
try:
    import warnings
    import statsmodels.formula.api as smf
    warnings.filterwarnings("ignore")
    sm_chk = []
    for y in ("r2s_PUT", "chi_PUT", "chi_GP", "chi_dia_PUT"):
        du = usableD[usableD[y].notna()].copy()
        for cn, f in (("age only", f"{y} ~ age"), ("+ sex + cell", f"{y} ~ age + sex_e + PS_e*ADHD_e")):
            fit = smf.mixedlm(f, du, groups=du["participant_id"]).fit(reml=True, method="powell")
            mine = DRT[(DRT.outcome == y) & (DRT["sample"] == "all usable sessions (lme REML, random intercept)") & (DRT.covariates == cn)].iloc[0]
            sm_chk.append(dict(outcome=y, covariates=cn, sm_est=float(fit.fe_params["age"]), mine_est=float(mine.slope_per_y),
                               sm_se=float(fit.bse_fe["age"]), mine_se=float(mine.se)))
    SMCHK = pd.DataFrame(sm_chk)
    SMCHK["rel_est"] = (SMCHK.sm_est - SMCHK.mine_est).abs() / SMCHK.mine_est.abs()
    SMCHK["rel_se"] = (SMCHK.sm_se - SMCHK.mine_se).abs() / SMCHK.mine_se.abs()
    log("statsmodels MixedLM (powell) vs own REML: max rel diff of the age estimate %.2e "
        "(statsmodels SEs use a different covariance route and are listed only, not compared)" % SMCHK.rel_est.max())
except Exception as e:  # pragma: no cover
    SMCHK = pd.DataFrame()
    log("statsmodels MixedLM check failed:", e)

# ================================================================================================================
# write own outputs
# ================================================================================================================
mine = dict(A_X2_chisep_group_effects=A, A_X2_component_descriptives=A_desc, B_nuisance_models=B, B_joint_model_coefficients=JC,
            B_r2s_models=R2T, B_dice_tertile_cells=TT, B_dice_tertile_models=TR, C_low_TD_tail_sessions=CS,
            C_low_TD_tail_qc_comparison=QC, C_low_TD_tail_refits=CR, D_age_slopes=DRT, D_implied_chi_slope_from_R2s=IMPL)
os.makedirs(f"{OUT}/mine", exist_ok=True)
for k, v in mine.items():
    v.to_csv(f"{OUT}/mine/{k}.csv", index=False)
SHORT.to_csv(f"{OUT}/mine/D_putamen_shortfall_by_Dr.csv", index=False)
SC = dict(bound=BOUND, sigma_B=sigma_B, delta_B=dB, rL=rL, repro={k: list(v) for k, v in REPRO.items()},
          identity_max_abs_err=ident_err, R2_cancellation_max_beta_range=canc_b, R2_cancellation_max_se_range=canc_se,
          n_dia_pos_R2_10_SNVTA=n_dia_pos_R2_10, vif_PS_joint=vif_joint, vif_PS_M1=vif_m1, dice_tertile_cuts=list(cuts),
          fisher_p_cell=fisher_cell, fisher_p_PS=fisher_ps, het=HET, spearman_dice_chi=dict(rho=rho, p=rho_p),
          frozen_sensitivity_compare=sens_cmp.to_dict("records"), K2_frozen=dict(est=K2["est"], se=K2["se"]),
          qc_grade=qc_grade, correlations=CORR, plausibility_edits_usable=plaus_changes,
          statsmodels_mixedlm_check=SMCHK.to_dict("records"),
          pc2x_frozen=dict(chi_PUT=frozen["PC2x"]["tests"]["chi_PUT"]["slope_ppb_per_y"], chi_GP=frozen["PC2x"]["tests"]["chi_GP"]["slope_ppb_per_y"]))
with open(f"{OUT}/mine/scalars.json", "w") as fh:
    json.dump(SC, fh, indent=2, default=float)

# ================================================================================================================
# compare with the outputs under check (every numeric cell; relative tolerance 1e-6, absolute floor 1e-9)
# ================================================================================================================
KEYS = dict(A_X2_chisep_group_effects=["model"], A_X2_component_descriptives=["roi", "Dr_label"], B_nuisance_models=["model"],
            B_joint_model_coefficients=["term"], B_r2s_models=["model"], B_dice_tertile_cells=["tertile"],
            B_dice_tertile_models=["model"], C_low_TD_tail_sessions=["code"], C_low_TD_tail_qc_comparison=["metric"],
            C_low_TD_tail_refits=["model"], D_age_slopes=["outcome", "sample", "covariates"],
            D_implied_chi_slope_from_R2s=["roi", "sample", "covariates"])
RTOL, AFLOOR = 1e-6, 1e-9
diffs, summary = [], []
for name, keys in KEYS.items():
    theirs = pd.read_csv(f"{VAL}/{name}.csv")
    my = mine[name].copy()
    if name == "C_low_TD_tail_sessions":
        theirs["code"] = theirs["code"].astype(str)
        my["code"] = my["code"].astype(str)
    mt = theirs.merge(my, on=keys, how="outer", suffixes=("_their", "_mine"), indicator=True)
    unmatched = int((mt["_merge"] != "both").sum())
    ncmp, nbad, maxrel = 0, 0, 0.0
    common = [c for c in theirs.columns if c in my.columns and c not in keys]
    for c in common:
        a, b = mt[c + "_their"], mt[c + "_mine"]
        for i in range(len(mt)):
            ta, tb = a.iloc[i], b.iloc[i]
            if isinstance(ta, str) or isinstance(tb, str):
                if pd.isna(ta) and pd.isna(tb):
                    continue
                ncmp += 1
                if str(ta) != str(tb):
                    nbad += 1
                    diffs.append(dict(file=name, key=" | ".join(str(mt[k].iloc[i]) for k in keys), column=c, theirs=ta, mine=tb, rel=np.nan))
                continue
            if pd.isna(ta) and pd.isna(tb):
                continue
            ncmp += 1
            if pd.isna(ta) != pd.isna(tb):
                nbad += 1
                diffs.append(dict(file=name, key=" | ".join(str(mt[k].iloc[i]) for k in keys), column=c, theirs=ta, mine=tb, rel=np.nan))
                continue
            ad = abs(float(ta) - float(tb))
            rel = ad / max(abs(float(ta)), 1e-300)
            if ad > AFLOOR and rel > RTOL:
                nbad += 1
                diffs.append(dict(file=name, key=" | ".join(str(mt[k].iloc[i]) for k in keys), column=c, theirs=ta, mine=tb, rel=rel))
            if ad > AFLOOR:
                maxrel = max(maxrel, rel)
    summary.append(dict(file=name, rows_theirs=len(theirs), rows_mine=len(my), unmatched_rows=unmatched, cells_compared=ncmp,
                        cells_outside_tol=nbad, max_rel_diff=maxrel, their_cols_not_recomputed=";".join(c for c in theirs.columns if c not in my.columns and c not in keys)))
SUM = pd.DataFrame(summary)
DIFF = pd.DataFrame(diffs)
SUM.to_csv(f"{OUT}/compare_csv_summary.csv", index=False)
DIFF.to_csv(f"{OUT}/compare_csv_discrepancies.csv", index=False)
log(SUM.to_string())
if len(DIFF):
    log(DIFF.to_string())

# results.json scalars
rj = json.load(open(f"{VAL}/results.json"))
SCMP = []


def sc(name, theirs, mine_v):
    ad = abs(theirs - mine_v)
    SCMP.append(dict(quantity=name, theirs=theirs, mine=mine_v, abs_diff=ad, rel_diff=ad / max(abs(theirs), 1e-300)))


sc("bound_ppb", rj["meta"]["bound_ppb"], BOUND)
sc("vif_PS_joint", rj["B_nuisance"]["vif_PS_joint"], vif_joint)
sc("vif_PS_M1", rj["B_nuisance"]["vif_PS_M1"], vif_m1)
sc("dice_cut_1", rj["B_nuisance"]["dice_tertile_cuts"][0], cuts[0])
sc("dice_cut_2", rj["B_nuisance"]["dice_tertile_cuts"][1], cuts[1])
sc("fisher_p_cell", rj["B_nuisance"]["dice_tertile_fisher_p_cell"], fisher_cell)
sc("fisher_p_PS", rj["B_nuisance"]["dice_tertile_fisher_p_PS"], fisher_ps)
for kk in ("wald_chi2", "F", "df2", "p_F"):
    sc("het_" + kk, rj["B_nuisance"]["dice_heterogeneity"][kk], HET[kk])
sc("spearman_rho", rj["B_nuisance"]["spearman_dice_chi"]["rho"], rho)
sc("spearman_p", rj["B_nuisance"]["spearman_dice_chi"]["p"], rho_p)
for kk in ("r_chi_r2s_first", "r_chi_r2s_first_age_partial", "r_chi_r2s_GP_first", "n_subjects_with_2plus_usable"):
    sc("D_" + kk, rj["D_positive_control"]["correlations"][kk], CORR[kk])
sc("A_identity_max_abs_err(<1e-12 expected)", rj["A_X2"]["identity_max_abs_err"], ident_err)
sc("A_R2_cancel_beta_range(<1e-12 expected)", rj["A_X2"]["R2_cancellation_max_beta_range"], canc_b)
SCMP = pd.DataFrame(SCMP)
SCMP.to_csv(f"{OUT}/compare_results_json_scalars.csv", index=False)
log(SCMP.to_string())
with open(f"{OUT}/verify_independent.log", "w") as fh:
    fh.write("\n".join(LOG) + "\n")
