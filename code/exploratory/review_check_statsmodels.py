#!/usr/bin/env python3
"""LEAD-chi TASK R: second-engine check (statsmodels, HC3, use_t) of the OLS decision quantities.
EXPLORATORY / REVIEW ONLY. Own sample rules and own natural-spline basis (truncated-power construction with knots at the
age terciles and boundary knots at the range; it spans the same space as splines::ns(df = 3), so the non-spline
coefficients and their HC3 SEs are invariant to the basis choice)."""
import json
import numpy as np
import pandas as pd
import statsmodels.api as sm

IN = "/path/to/Penn_LEAD/LEADchi/mirror/unlocked/merged_qsmxt_primary_SNVTA.csv"
GATE = "/path/to/Penn_LEAD/LEADchi/mirror/gates/gates_preunlock.json"
RES = "/path/to/Penn_LEAD/LEADchi/mirror/results/confirmatory/results.json"
OUT = "/path/to/Penn_LEAD/LEADchi/results_exploratory/review/review_check_statsmodels.json"

d = pd.read_csv(IN)
g = json.load(open(GATE))
R = json.load(open(RES))
d["PS_e"] = np.where(d.study_group == "PRO/CHR", 0.5, -0.5)
d["ADHD_e"] = np.where(d.dx_adhd == 1, 0.5, -0.5)
d["sex_e"] = np.where(d.sex == "M", 0.5, -0.5)
d["race3"] = np.where(d.race == "Black or African American", "Black", np.where(d.race == "White", "White", "Other"))
d = d.sort_values(["participant_id", "ses_order"])
u = d[d.usable == 1]
first = u.drop_duplicates("participant_id").copy()
assert len(first) == 80


def ns3(x):
    x = np.asarray(x, float)
    k = [x.min(), *np.quantile(x, [1 / 3, 2 / 3]), x.max()]
    def dk(j):
        return (np.clip(x - k[j], 0, None) ** 3 - np.clip(x - k[3], 0, None) ** 3) / (k[3] - k[j])
    return np.column_stack([x, dk(0) - dk(2), dk(1) - dk(2)])


def zs(v):
    v = np.asarray(v, float)
    return (v - v.mean()) / v.std(ddof=1)


def design(df, lead, cont):
    cols, names = [np.ones(len(df))], ["const"]
    for nm, v in lead:
        cols.append(np.asarray(v, float)); names.append(nm)
    B = ns3(df.age)
    for j in range(3):
        cols.append(B[:, j]); names.append(f"ns{j + 1}")
    cols.append(df.sex_e.to_numpy(float)); names.append("sex_e")
    lev = [l for l in ("Black", "White", "Other") if l in set(df.race3)]
    for l in lev[1:]:   # treatment dummies: a reparameterisation of contr.sum, same column space
        cols.append((df.race3 == l).to_numpy(float)); names.append("race_" + l)
    for c in cont:
        cols.append(zs(df[c])); names.append("z_" + c)
    return np.column_stack(cols), names


def fit(y, X, names):
    return sm.OLS(np.asarray(y, float), X).fit(cov_type="HC3", use_t=True), names


def cvec(names, w):
    c = np.zeros(len(names))
    for k, v in w.items():
        c[names.index(k)] = v
    return c


checks = []
def add(q, py, r):
    checks.append({"quantity": q, "python": float(py), "R_frozen": float(r), "rel_diff": float(abs(py - r) / max(abs(r), 1e-12))})

# M1
X, nm = design(first, [("PS_e", first.PS_e), ("ADHD_e", first.ADHD_e), ("PSxADHD", first.PS_e * first.ADHD_e)], ["roi_vol", "Q_phase"])
f, nm = fit(first.chi_SNVTA, X, nm)
i = nm.index("PS_e")
add("M1 beta_PS", f.params[i], R["H1"]["decision"]["beta_ppb"])
add("M1 SE_HC3", f.bse[i], R["H1"]["decision"]["se_hc3"])
add("M1 p_two", f.pvalues[i], R["H1"]["decision"]["p_two"])
add("M1 df", f.df_resid, R["H1"]["decision"]["df"])
add("M1 sigma_resid", np.sqrt(f.scale), R["H1"]["decision"]["sigma_resid"])
for lab, w, key in [("K1 PS-only - TD", {"PS_e": 1, "PSxADHD": -0.5}, "simple_PSonly_minus_TD"),
                    ("K1 PS+ADHD - ADHD-only", {"PS_e": 1, "PSxADHD": 0.5}, "simple_PSADHD_minus_ADHDonly")]:
    t = f.t_test(cvec(nm, w))
    add(lab + " est", float(np.squeeze(t.effect)), R["H1"]["K1"][key]["est"])
    add(lab + " SE", float(np.squeeze(t.sd)), R["H1"]["K1"][key]["se"])
    add(lab + " p", float(np.squeeze(t.pvalue)), R["H1"]["K1"][key]["p_two"])
# bound test (registered quantities from the frozen scale; the Huber scale itself is an R-only quantity)
dB = 0.6 * np.sqrt(min(1, g["rL_primary"]))
sigB = R["H1"]["decision"]["sigma_B"]
from scipy import stats
tB = (f.params[i] + dB * sigB) / f.bse[i]
add("M1 p_B (frozen sigma_B)", stats.t.sf(tB, f.df_resid), R["H1"]["decision"]["B_registered"]["pB"])
# K2, K3
f2, nm2 = fit(first.r2s_SNVTA, X, nm); add("K2 R2* beta_PS", f2.params[i], R["H1"]["K2"]["est"])
f3, nm3 = fit(first.chi_SNVTA - first.chi_SNVTA_ref2, X, nm); add("K3 beta_dref", f3.params[i], R["H1"]["K3"]["beta_dref"])
add("K3 p", f3.pvalues[i], R["H1"]["K3"]["p_two"])
# S1
s1 = first.copy()
s1["D"] = (s1.chi_SNVTA - g["pooled"]["mean_SNVTA"]) / g["pooled"]["sd_SNVTA"] - (s1.chi_STR - g["pooled"]["mean_STR"]) / g["pooled"]["sd_STR"]
X1, nm1 = design(s1, [("PS_e", s1.PS_e), ("ADHD_e", s1.ADHD_e), ("PSxADHD", s1.PS_e * s1.ADHD_e)], ["roi_vol", "str_vol", "Q_phase"])
fs1, nm1 = fit(s1.D, X1, nm1)
t = fs1.t_test(cvec(nm1, {"PS_e": 1, "ADHD_e": -1}))
add("S1 L", float(np.squeeze(t.effect)), R["S1"]["L"]["est"]); add("S1 SE(L)", float(np.squeeze(t.sd)), R["S1"]["L"]["se"]); add("S1 p", float(np.squeeze(t.pvalue)), R["S1"]["L"]["p_two"])
# S2
s2 = first[(first.age >= 11) & first.PRIME_total.notna() & first.BDI_total.notna() & first.ARI_total.notna()].copy()
X2, nm2 = design(s2, [("z_PRIME", zs(s2.PRIME_total)), ("z_BDI", zs(s2.BDI_total)), ("z_ARI", zs(s2.ARI_total)), ("ADHD_e", s2.ADHD_e)], ["roi_vol", "Q_phase"])
fs2, nm2 = fit(s2.chi_SNVTA, X2, nm2); j = nm2.index("z_PRIME")
add("S2 n", len(s2), R["S2"]["n"]); add("S2 beta_PRIME", fs2.params[j], R["S2"]["beta_PRIME"]["est"]); add("S2 SE", fs2.bse[j], R["S2"]["beta_PRIME"]["se"]); add("S2 p", fs2.pvalues[j], R["S2"]["beta_PRIME"]["p_two"])
# Holm
from statsmodels.stats.multitest import multipletests
ph = multipletests([R["S1"]["L"]["p_two"], R["S2"]["beta_PRIME"]["p_two"]], method="holm")[1]
add("Holm p S1", ph[0], R["holm"]["p_holm"]["S1"]); add("Holm p S2", ph[1], R["holm"]["p_holm"]["S2"])

ok = all(c["rel_diff"] < 1e-6 for c in checks)
json.dump({"all_within_1e-6": ok, "n_checks": len(checks), "checks": checks,
           "note": "second-engine check written for the adversarial review; sample membership and spline basis re-derived independently"},
          open(OUT, "w"), indent=1)
for c in checks:
    print(f"{c['quantity']:<28} py {c['python']:.10g}  R {c['R_frozen']:.10g}  rel {c['rel_diff']:.1e}")
print("ALL WITHIN 1e-6:", ok)
