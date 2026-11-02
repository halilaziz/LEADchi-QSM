#!/usr/bin/env python3
"""verify_statsmodels.py - independent re-implementation check of the LEAD-chi confirmatory estimates (Q19 verification).

  python3 verify_statsmodels.py --input <merged table> --gates <gates JSON> --results <run_all.R results.json> --out <json>

Re-implements, in Python (numpy / scipy / statsmodels, not R), the decision quantities of M1 (beta_PS, HC3 SE, two-sided p;
the two simple effects), S1 (L with the pooled gate constants), S2 (beta_PRIME), S3a (beta of scale(chi_DN_R)) and T1
(beta and SE of PS_e:occ at the REML maximum, by a direct numpy/scipy REML profile), on the same samples as run_all.R (the sample membership is re-derived
from the table with the registered first-usable / S2 / S3 / T1 rules, after the same plausibility re-verification is
applied in R: rows that R set to missing are read back from results.json's plausibility log).
Tolerances: 1e-6 (relative) for OLS quantities, 1e-4 for the mixed model. The natural-spline basis is built here
(knots at the 1/3 and 2/3 quantiles, boundary knots at the range, as splines::ns(df = 3)); it spans the same space
as R's basis, so every non-spline coefficient and its HC3 SE are identical. Declared limitation: both implementations
were written by the same AI system (not independent analysts), and the parametric bootstrap and permutation p values
are Monte Carlo quantities that are not re-computed here. A difference beyond tolerance is logged in the report and,
if it concerns a decision, the analysis is halted until the coding error is found (Annex C.4).
"""
from __future__ import annotations

import argparse
import json
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

TOL_OLS, TOL_MIXED = 1e-6, 1e-4


def ns_basis(x, knots_from=None):
    x = np.asarray(x, float)
    ref = x if knots_from is None else np.asarray(knots_from, float)
    b0, b1 = ref.min(), ref.max()
    k1, k2 = np.quantile(ref, [1 / 3, 2 / 3])
    kn = [b0, k1, k2, b1]
    def d(k):
        return (np.clip(x - kn[k], 0, None) ** 3 - np.clip(x - kn[3], 0, None) ** 3) / (kn[3] - kn[k])
    return np.column_stack([x, d(0) - d(2), d(1) - d(2)])


def z(v):
    v = np.asarray(v, float)
    return (v - v.mean()) / v.std(ddof=1)


def code(d):
    d = d.copy()
    d["PS_e"] = np.where(d.study_group == "PRO/CHR", 0.5, -0.5)
    d["ADHD_e"] = np.where(d.dx_adhd == 1, 0.5, -0.5)
    d["sex_e"] = np.where(d.sex == "M", 0.5, -0.5)
    d["race3"] = np.where(d.race == "Black or African American", "Black", np.where(d.race == "White", "White", "Other"))
    d["sess_key"] = d.participant_id + "_" + d.session_id
    return d


def design(d, extra=None, group=("PS_e", "ADHD_e"), with_int=True):
    cols = [np.ones(len(d))]
    names = ["(Intercept)"]
    for g in group:
        cols.append(d[g].to_numpy(float)); names.append(g)
    if with_int and len(group) == 2:
        cols.append((d[group[0]] * d[group[1]]).to_numpy(float)); names.append("PS_e:ADHD_e")
    B = ns_basis(d["age"])
    for j in range(3):
        cols.append(B[:, j]); names.append(f"ns{j + 1}")
    cols.append(d["sex_e"].to_numpy(float)); names.append("sex_e")
    lev = [l for l in ("Black", "White", "Other") if l in set(d.race3)]
    for l in lev[1:]:
        cols.append((d.race3 == l).to_numpy(float)); names.append(f"race_{l}")
    for e in (extra or []):
        cols.append(z(d[e])); names.append(f"scale({e})")
    return np.column_stack(cols), names


def ols(y, X, names):
    f = sm.OLS(np.asarray(y, float), X).fit(cov_type="HC3", use_t=True)
    return {n: (float(f.params[i]), float(f.bse[i]), float(f.pvalues[i])) for i, n in enumerate(names)}, f


def rel(a, b):
    return abs(a - b) / max(1.0, abs(b))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True); ap.add_argument("--gates", required=True)
    ap.add_argument("--results", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    R = json.load(open(a.results)); G = json.load(open(a.gates))
    d = code(pd.read_csv(a.input, na_values=["NA", ""]))
    # apply the plausibility actions R logged (values set to missing; CSF reference invalidated)
    for row in (R.get("plausibility", {}).get("log") or []):
        if isinstance(row, dict) and row.get("action", "").startswith(("set to missing", "CSF reference invalid")):
            d.loc[d.sess_key == row["sess_key"], row["column"]] = np.nan
    d = d.sort_values(["participant_id", "ses_order"])
    u = d[d.usable == 1].copy()
    u["k"] = u.groupby("participant_id").cumcount() + 1
    first = u[u.k == 1].copy(); occ1 = u[u.k == 2].copy()
    checks = []
    def add(name, mine, theirs, tol):
        checks.append({"quantity": name, "python": mine, "R": theirs, "rel_diff": rel(mine, theirs), "pass": rel(mine, theirs) <= tol})
    # M1
    X, nm = design(first, ["roi_vol", "Q_phase"])
    res, f = ols(first.chi_SNVTA, X, nm)
    dec = R["H1"]["decision"]
    add("M1 beta_PS", res["PS_e"][0], dec["beta_ppb"], TOL_OLS); add("M1 SE_HC3", res["PS_e"][1], dec["se_hc3"], TOL_OLS)
    add("M1 p_two", res["PS_e"][2], dec["p_two"], TOL_OLS)
    V = f.cov_params(); ip, ii = nm.index("PS_e"), nm.index("PS_e:ADHD_e")
    for w, key in ((-0.5, "simple_PSonly_minus_TD"), (0.5, "simple_PSADHD_minus_ADHDonly")):
        c = np.zeros(len(nm)); c[ip] = 1; c[ii] = w
        add(f"K1 {key}", float(c @ f.params), R["H1"]["K1"][key]["est"], TOL_OLS)
    # S1
    s1 = first.dropna(subset=["chi_STR", "str_vol"]).copy()
    P = G["pooled"]
    s1["D"] = (s1.chi_SNVTA - P["mean_SNVTA"]) / P["sd_SNVTA"] - (s1.chi_STR - P["mean_STR"]) / P["sd_STR"]
    X, nm = design(s1, ["roi_vol", "str_vol", "Q_phase"])
    _, f1 = ols(s1.D, X, nm)
    c = np.zeros(len(nm)); c[nm.index("PS_e")] = 1; c[nm.index("ADHD_e")] = -1
    add("S1 L", float(c @ f1.params), R["S1"]["L"]["est"], TOL_OLS)
    add("S1 SE(L)", float(np.sqrt(c @ f1.cov_params() @ c)), R["S1"]["L"]["se"], TOL_OLS)
    # S2
    s2 = first[(first.age >= 11) & first.PRIME_total.notna() & first.BDI_total.notna() & first.ARI_total.notna()].copy()
    X, nm = design(s2, ["roi_vol", "Q_phase"], group=("ADHD_e",), with_int=False)
    X = np.column_stack([X, z(s2.PRIME_total), z(s2.BDI_total), z(s2.ARI_total)]); nm += ["PRIME", "BDI", "ARI"]
    r2, _ = ols(s2.chi_SNVTA, X, nm)
    add("S2 beta_PRIME", r2["PRIME"][0], R["S2"]["beta_PRIME"]["est"], TOL_OLS); add("S2 SE", r2["PRIME"][1], R["S2"]["beta_PRIME"]["se"], TOL_OLS)
    # S3
    s3v = ["twICVF_pSCP_R", "b0loc_pSCP_R", "tract_len", "nSL", "FD", "NDC", "chi_DN_R", "chi_DGM", "DN_R_vol", "abs_pitch"]
    s3 = first[first.s3_eligible == 1].dropna(subset=s3v).copy()
    s3["lognSL"] = np.log(s3.nSL)
    X, nm = design(s3, ["chi_DGM", "b0loc_pSCP_R", "abs_pitch", "tract_len", "FD", "NDC", "lognSL", "DN_R_vol", "Q_phase"])
    X = np.column_stack([X, z(s3.chi_DN_R)]); nm += ["DN"]
    r3, _ = ols(z(s3.twICVF_pSCP_R), X, nm)
    add("S3a beta", r3["DN"][0], R["S3"]["beta_DN"]["est"], TOL_OLS); add("S3a SE", r3["DN"][1], R["S3"]["beta_DN"]["se"], TOL_OLS)
    # T1 (REML, random intercept)
    first["occ"] = 0; occ1["occ"] = 1
    two = pd.concat([first, occ1])
    two["y"] = two.chi_SNVTA / P["sd_total_SNVTA"]; two["PSxocc"] = two.PS_e * two.occ
    X, nm = design(two, ["roi_vol", "Q_phase"])
    X = np.column_stack([X, two.occ, two.PSxocc]); nm += ["occ", "PSxocc"]
    # REML for the random-intercept model, computed here directly (numpy/scipy; independent of lme4 and of the
    # statsmodels MixedLM optimiser, which can fail numerically): profile the REML log-likelihood over
    # theta = sigma_subject / sigma_resid (grid 0..5 step 0.01, then bounded Brent), GLS fixed effects and
    # cov(beta) = sigma^2 (X' V^-1 X)^-1 at the maximum (as lme4's vcov).
    from scipy.optimize import minimize_scalar
    yv = two.y.to_numpy(float); grp = two.participant_id.to_numpy(); groups = [np.where(grp == g)[0] for g in np.unique(grp)]
    n_, p_ = X.shape
    def reml_parts(theta):
        A = np.zeros((p_, p_)); bx = np.zeros(p_); logdet = 0.0
        for ix in groups:
            ng = len(ix); Vg = np.eye(ng) + theta ** 2 * np.ones((ng, ng)); Vi = np.linalg.inv(Vg)
            A += X[ix].T @ Vi @ X[ix]; bx += X[ix].T @ Vi @ yv[ix]; logdet += np.linalg.slogdet(Vg)[1]
        beta = np.linalg.solve(A, bx); q = 0.0
        for ix in groups:
            ng = len(ix); r = yv[ix] - X[ix] @ beta
            q += r @ np.linalg.solve(np.eye(ng) + theta ** 2 * np.ones((ng, ng)), r)
        s2 = q / (n_ - p_)
        ll = -0.5 * ((n_ - p_) * np.log(s2) + logdet + np.linalg.slogdet(A)[1])
        return ll, beta, s2, A
    grid = np.arange(0, 5.0001, 0.01); lls = [reml_parts(t)[0] for t in grid]; t0 = grid[int(np.argmax(lls))]
    opt = minimize_scalar(lambda t: -reml_parts(t)[0], bounds=(max(0, t0 - 0.02), t0 + 0.02), method="bounded", options={"xatol": 1e-8})
    ll, beta, s2, A = reml_parts(opt.x)
    j = nm.index("PSxocc")
    se_gls = float(np.sqrt(s2 * np.linalg.inv(A)[j, j]))
    class _MM: pass
    mm = _MM(); mm.fe_params = beta
    add("T1 beta(PS x occ)", float(mm.fe_params[j]), R["T1"]["beta_PSxocc"], TOL_MIXED)
    add("T1 SE", se_gls, R["T1"]["se"], TOL_MIXED)
    out = {"n_checks": len(checks), "all_pass": all(c["pass"] for c in checks), "checks": checks,
           "note": "same-AI re-implementation; resampling p values not re-computed; tolerance 1e-6 (OLS) / 1e-4 (mixed)"}
    json.dump(out, open(a.out, "w"), indent=1)
    for c in checks:
        print(f"{c['quantity']:34s} python {c['python']: .8f}  R {c['R']: .8f}  rel {c['rel_diff']:.1e}  {'PASS' if c['pass'] else 'FAIL'}")
    print(f"{sum(c['pass'] for c in checks)} of {len(checks)} within tolerance")
    sys.exit(0 if out["all_pass"] else 1)


if __name__ == "__main__":
    main()
