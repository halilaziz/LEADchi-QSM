#!/usr/bin/env python3
"""T3a independent re-implementation (statsmodels / numpy) of the OLS-based Tier-3 quantities:
   (i)  E-int: the four adjusted cell means of M1 (first-usable, HC3) and the interaction;
   (ii) E-z : TD-referenced z for every first-usable subject (TD leave-one-subject-out; fixed TD knots).
Compares against the R outputs in results_exploratory/tier3_estimation/. Tolerance 1e-6 relative (OLS).
EXPLORATORY / TIER 3; a verification aid only. Declared: both implementations come from the same AI system."""
import json, sys
import numpy as np, pandas as pd
import statsmodels.api as sm
from scipy.stats import norm

ROOT = "/path/to/Penn_LEAD/LEADchi"
OUT = f"{ROOT}/results_exploratory/tier3_estimation"
d = pd.read_csv(f"{ROOT}/mirror/unlocked/merged_qsmxt_primary_SNVTA.csv")
d["PS_e"] = np.where(d.study_group == "PRO/CHR", 0.5, -0.5)
d["ADHD_e"] = np.where(d.dx_adhd == 1, 0.5, -0.5)
d["sex_e"] = np.where(d.sex == "M", 0.5, -0.5)
d["race3"] = np.where(d.race == "Black or African American", "Black", np.where(d.race == "White", "White", "Other"))
d["cell"] = np.where(d.PS_e > 0, np.where(d.ADHD_e > 0, "PS+ADHD", "PS-only"), np.where(d.ADHD_e > 0, "ADHD-only", "TD"))
fu = d[d.first_usable == 1].copy().reset_index(drop=True)
assert len(fu) == 80

def ns_basis(x, knots_from):
    x = np.asarray(x, float); ref = np.asarray(knots_from, float)
    b0, b1 = ref.min(), ref.max(); k1, k2 = np.quantile(ref, [1 / 3, 2 / 3]); kn = [b0, k1, k2, b1]
    def dd(k):
        return (np.clip(x - kn[k], 0, None) ** 3 - np.clip(x - kn[3], 0, None) ** 3) / (kn[3] - kn[k])
    return np.column_stack([x, dd(0) - dd(2), dd(1) - dd(2)])

def zsc(v):
    v = np.asarray(v, float); return (v - v.mean()) / v.std(ddof=1)

def race_cols(r):
    lev = [l for l in ["Black", "White", "Other"] if l in set(r)]
    return np.column_stack([np.where(r == lev[j], 1.0, np.where(r == lev[-1], -1.0, 0.0)) for j in range(len(lev) - 1)])

checks = []
def chk(name, a, b, tol=1e-6):
    a = np.asarray(a, float); b = np.asarray(b, float)
    rel = np.max(np.abs(a - b) / np.maximum(1e-12, np.abs(b)))
    checks.append({"check": name, "max_rel_diff": float(rel), "pass": bool(rel < tol)})

# ---------------- (i) E-int ----------------
X = np.column_stack([np.ones(80), fu.PS_e, fu.ADHD_e, ns_basis(fu.age, fu.age), fu.sex_e, race_cols(fu.race3.values),
                     zsc(fu.roi_vol), zsc(fu.Q_phase), fu.PS_e * fu.ADHD_e])
ols = sm.OLS(fu.chi_SNVTA.values, X).fit(cov_type="HC3")
b = ols.params; V = ols.cov_params()
frozen = json.load(open(f"{ROOT}/mirror/results/confirmatory/results.json"))
chk("E-int: beta_PS vs frozen", b[1], frozen["H1"]["decision"]["beta_ppb"])
chk("E-int: interaction vs frozen", b[-1], frozen["H1"]["K1"]["interaction"]["est"])
chk("E-int: interaction HC3 SE vs frozen", np.sqrt(V[-1, -1]), frozen["H1"]["K1"]["interaction"]["se"])
cm_r = pd.read_csv(f"{OUT}/eint_cell_means.csv")
codes = {"TD": (-0.5, -0.5), "ADHD-only": (-0.5, 0.5), "PS-only": (0.5, -0.5), "PS+ADHD": (0.5, 0.5)}
est_py, se_py = [], []
for cell in cm_r.cell:
    x = X.mean(axis=0).copy(); x[1], x[2] = codes[cell]; x[-1] = codes[cell][0] * codes[cell][1]
    est_py.append(float(x @ b)); se_py.append(float(np.sqrt(x @ V @ x)))
chk("E-int: adjusted cell means (4)", est_py, cm_r.adjusted_mean_ppb.values)
chk("E-int: adjusted cell-mean HC3 SEs (4)", se_py, cm_r.se_hc3.values)

# ---------------- (ii) E-z ----------------
ez_r = pd.read_csv(f"{OUT}/ez_subject_z.csv")
td = fu[fu.cell == "TD"].reset_index(drop=True); ref_age = td.age.values
def fit_td(sub):
    X = np.column_stack([np.ones(len(sub)), ns_basis(sub.age, ref_age), sub.sex_e]); f = sm.OLS(sub.chi_SNVTA.values, X).fit()
    sigma = np.sqrt(f.ssr / (len(sub) - X.shape[1])); return f.params, sigma
def pred(bb, row):
    return float(np.r_[1.0, ns_basis([row.age], ref_age)[0], row.sex_e] @ bb)
b_full, s_full = fit_td(td)
mu_py, sig_py, z_py = [], [], []
for _, r in ez_r.iterrows():
    row = fu[fu.participant_id == r["sub"]].iloc[0]
    if r.cell == "TD":
        bb, ss = fit_td(td[td.participant_id != r["sub"]])
    else:
        bb, ss = b_full, s_full
    mu = pred(bb, row); mu_py.append(mu); sig_py.append(ss); z_py.append((row.chi_SNVTA - mu) / ss)
chk("E-z: mu_TD (80; TD LOSO)", mu_py, ez_r.mu_TD.values)
chk("E-z: sigma_TD (80)", sig_py, ez_r.sigma_TD.values)
chk("E-z: z (80)", z_py, ez_r.z.values)
chk("E-z: centile (80)", 100 * norm.cdf(z_py), ez_r.centile.values)
sumr = pd.read_csv(f"{OUT}/ez_cell_summary.csv")
zz = pd.Series(z_py, index=ez_r.index); cells = ez_r.cell
chk("E-z: mean z per cell", [zz[cells == c].mean() for c in sumr.cell], sumr.mean_z.values)
chk("E-z: share z <= -1 per cell", [(zz[cells == c] <= -1).mean() for c in sumr.cell], sumr.share_z_le_m1.values)

res = pd.DataFrame(checks); print(res.to_string(index=False))
res.to_csv(f"{OUT}/t3a_verify_statsmodels.csv", index=False)
n_fail = int((~res["pass"]).sum())
print(f"\n{len(res) - n_fail}/{len(res)} checks within tolerance; max relative difference {res.max_rel_diff.max():.3g}")
sys.exit(1 if n_fail else 0)
