#!/usr/bin/env python3
"""Independent re-computation of LEADchi/results_posthoc/engine2 (POST HOC engine-2 comparison).

Written from the task specification only (no R code run or copied). Read-only on every input; writes only
verification/recomputed.json next to this script.
"""
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
E2 = os.path.dirname(HERE)
LX = os.path.abspath(os.path.join(E2, "..", ".."))
KEY = os.path.join(LX, "mirror", "unlocked", "key_plain.csv")
MERGED = os.path.join(LX, "mirror", "unlocked", "merged_qsmxt_primary_SNVTA.csv")
REC_PRIMARY_DIR = os.path.join(LX, "results_posthoc", "recovery", "raw")

RL = 0.938278557806515
DELTA_B = 0.6 * np.sqrt(RL)
ADULT = -14.24
HUBER_K = 1.345


# ----------------------------------------------------------------------------------------------- design
def code_vars(d):
    d = d.copy()
    d["PS_e"] = np.where(d["study_group"] == "PRO/CHR", 0.5, -0.5)
    d["ADHD_e"] = np.where(d["dx_adhd"] == 1, 0.5, -0.5)
    d["sex_e"] = np.where(d["sex"] == "M", 0.5, -0.5)
    d["race3"] = np.where(d["race"] == "Black or African American", "Black",
                          np.where(d["race"] == "White", "White", "Other"))
    d["sess_key"] = d["participant_id"] + "_" + d["session_id"]
    return d


def ns_basis(x, df=3):
    """Natural cubic spline basis (ESL eq. 5.4-5.5), same column space (with intercept) as R splines::ns(x, df)."""
    x = np.asarray(x, float)
    inner = np.quantile(x, [1 / 3, 2 / 3])           # type-7 quantiles, as R
    knots = np.r_[x.min(), inner, x.max()]
    K = len(knots)

    def dk(k):
        return (np.clip(x - knots[k], 0, None) ** 3 - np.clip(x - knots[K - 1], 0, None) ** 3) / (knots[K - 1] - knots[k])

    cols = [x] + [dk(k) - dk(K - 2) for k in range(K - 2)]
    return np.column_stack(cols)


def zscore(v):
    v = np.asarray(v, float)
    return (v - v.mean()) / v.std(ddof=1)


def design(d):
    """Columns: Intercept, PS, ADHD, ns1-3, sex, race3 (sum-to-zero over levels present), z(roi_vol), z(Q_phase), PS:ADHD."""
    n = len(d)
    cols = {"(Intercept)": np.ones(n), "PS_e": d["PS_e"].values, "ADHD_e": d["ADHD_e"].values}
    B = ns_basis(d["age"].values)
    for j in range(B.shape[1]):
        cols[f"ns{j + 1}"] = B[:, j]
    cols["sex_e"] = d["sex_e"].values
    lev = [l for l in ["Black", "White", "Other"] if l in set(d["race3"])]
    for j, l in enumerate(lev[:-1]):
        cols[f"race3_{j + 1}"] = np.where(d["race3"] == l, 1.0, np.where(d["race3"] == lev[-1], -1.0, 0.0))
    cols["z_roi_vol"] = zscore(d["roi_vol"])
    cols["z_Q_phase"] = zscore(d["Q_phase"])
    cols["PS_e:ADHD_e"] = d["PS_e"].values * d["ADHD_e"].values
    names = list(cols)
    return np.column_stack([cols[c] for c in names]), names


# ----------------------------------------------------------------------------------------------- estimators
def ols_hc3(X, y):
    n, p = X.shape
    XtXi = np.linalg.inv(X.T @ X)
    b = XtXi @ X.T @ y
    e = y - X @ b
    h = np.einsum("ij,jk,ik->i", X, XtXi, X)
    meat = (X * (e / (1 - h))[:, None] ** 2).T @ X
    V = XtXi @ meat @ XtXi
    df = n - p
    sigma = np.sqrt(e @ e / df)
    return dict(b=b, V=V, e=e, df=df, sigma=sigma, n=n, p=p)


def huber_rlm(X, y, k=HUBER_K, k2=HUBER_K, maxit=200, acc=1e-4):
    """Huber M-estimate with Huber proposal-2 scale; IRLS as described for MASS::rlm (LS start, MAD start scale,
    scale updated before each reweighting, convergence on the relative change of the residual vector)."""
    n, p = X.shape
    b = np.linalg.lstsq(X, y, rcond=None)[0]
    r = y - X @ b
    n1 = n - p
    theta = 2 * stats.norm.cdf(k2) - 1
    gamma = theta + k2 ** 2 * (1 - theta) - 2 * k2 * stats.norm.pdf(k2)
    scale = 1.4826 * np.median(np.abs(r))
    converged = False
    it = 0
    for it in range(1, maxit + 1):
        r_old = r.copy()
        scale = np.sqrt(np.sum(np.minimum(r ** 2, (k2 * scale) ** 2)) / (n1 * gamma))
        u = r / scale
        w = np.minimum(1.0, k / np.maximum(np.abs(u), 1e-300))
        sw = np.sqrt(w)
        b = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)[0]
        r = y - X @ b
        conv = np.sqrt(np.sum((r_old - r) ** 2) / max(1e-20, np.sum(r_old ** 2)))
        if conv <= acc:
            converged = True
            break
    return dict(b=b, scale=scale, converged=converged, iters=it)


def huber_fixed_point(X, y, k=HUBER_K, tol=1e-14, maxit=10000):
    """Same estimator iterated to full numerical convergence (to show how much the acc = 1e-4 stop matters)."""
    n, p = X.shape
    b = np.linalg.lstsq(X, y, rcond=None)[0]
    r = y - X @ b
    theta = 2 * stats.norm.cdf(k) - 1
    gamma = theta + k ** 2 * (1 - theta) - 2 * k * stats.norm.pdf(k)
    s = 1.4826 * np.median(np.abs(r))
    for _ in range(maxit):
        s_new = np.sqrt(np.sum(np.minimum(r ** 2, (k * s) ** 2)) / ((n - p) * gamma))
        w = np.minimum(1.0, k / np.maximum(np.abs(r / s_new), 1e-300))
        sw = np.sqrt(w)
        b_new = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)[0]
        r = y - X @ b_new
        if abs(s_new - s) < tol * s and np.max(np.abs(b_new - b)) < tol * (1 + np.max(np.abs(b))):
            b, s = b_new, s_new
            break
        b, s = b_new, s_new
    return dict(b=b, scale=s)


def lin(m, c):
    est = float(c @ m["b"])
    se = float(np.sqrt(c @ m["V"] @ c))
    t = est / se
    q = stats.t.ppf(0.975, m["df"])
    return dict(est=est, se=se, t=t, p_two=float(2 * stats.t.sf(abs(t), m["df"])), ci95=[est - q * se, est + q * se])


def run_m1(d):
    X, names = design(d)
    y = d["chi"].values.astype(float)
    assert np.linalg.matrix_rank(X) == X.shape[1] == 12
    m = ols_hc3(X, y)
    iPS, iINT = names.index("PS_e"), names.index("PS_e:ADHD_e")
    cPS = np.zeros(len(names)); cPS[iPS] = 1
    tt = lin(m, cPS)
    hub = huber_rlm(X, y)
    hub_fp = huber_fixed_point(X, y)
    sigma_B = min(m["sigma"], hub["scale"])
    tB = (tt["est"] + DELTA_B * sigma_B) / tt["se"]
    pB = float(stats.t.sf(tB, m["df"]))
    A = tt["est"] < 0 and tt["p_two"] < 0.05
    B = pB < 0.05
    outcome = "A+B" if (A and B) else "A" if A else "B" if B else "C"
    c1 = cPS.copy(); c1[iINT] = -0.5
    c2 = cPS.copy(); c2[iINT] = 0.5
    s1, s2 = lin(m, c1), lin(m, c2)
    for s in (s1, s2):
        s["tB"] = (s["est"] + DELTA_B * sigma_B) / s["se"]
        s["pB_one_sided"] = float(stats.t.sf(s["tB"], m["df"]))
    cI = np.zeros(len(names)); cI[iINT] = 1
    return dict(n=m["n"], df=m["df"], beta=tt["est"], se=tt["se"], p_two=tt["p_two"], ci95=tt["ci95"],
                sigma_resid=m["sigma"], huber_scale=hub["scale"], huber_converged=hub["converged"],
                huber_iters=hub["iters"], huber_beta=float(hub["b"][iPS]),
                huber_scale_fixed_point=hub_fp["scale"], huber_beta_fixed_point=float(hub_fp["b"][iPS]),
                sigma_B=sigma_B, delta=tt["est"] / m["sigma"], delta_B=DELTA_B, bound=-DELTA_B * sigma_B,
                tB=tB, pB=pB, outcome=outcome, disconfirming=bool(tt["est"] > 0 and tt["p_two"] < 0.05),
                PSonly_minus_TD=s1, PSADHD_minus_ADHD=s2, interaction=lin(m, cI))


def statsmodels_crosscheck(d):
    import statsmodels.api as sm
    X, names = design(d)
    y = d["chi"].values.astype(float)
    ols = sm.OLS(y, X).fit(cov_type="HC3")
    rlm = sm.RLM(y, X, M=sm.robust.norms.HuberT(t=HUBER_K)).fit(scale_est=sm.robust.scale.HuberScale(d=HUBER_K),
                                                                   maxiter=500, tol=1e-12, conv="coefs")
    i = names.index("PS_e")
    return dict(beta=float(ols.params[i]), se_hc3=float(ols.bse[i]), rlm_beta=float(rlm.params[i]),
                rlm_scale=float(rlm.scale))


# ----------------------------------------------------------------------------------------------- agreement
def icc31(a, b):
    Y = np.column_stack([a, b]).astype(float)
    n, k = Y.shape
    gm = Y.mean()
    ssr = k * np.sum((Y.mean(1) - gm) ** 2)
    ssc = n * np.sum((Y.mean(0) - gm) ** 2)
    sse = np.sum((Y - gm) ** 2) - ssr - ssc
    msr = ssr / (n - 1)
    mse = sse / ((n - 1) * (k - 1))
    return (msr - mse) / (msr + (k - 1) * mse)


def read_variant(path, key):
    x = pd.read_csv(path, dtype={"code": str})
    x["ses_order"] = x["ses_order"].astype(int)
    x = x.merge(key[["code", "ses_order", "participant_id", "session_id"]], on=["code", "ses_order"], how="left",
                validate="one_to_one")
    assert x["participant_id"].notna().all()
    x["sess_key"] = x["participant_id"] + "_" + x["session_id"]
    return x


def recovery_k(base, rec, first_keys):
    kk = base[["sess_key", "status", "chi_SNVTA_raw"]].merge(rec[["sess_key", "status", "chi_SNVTA_raw"]],
                                                             on="sess_key", suffixes=("", "_rec"), validate="one_to_one")
    kk = kk[(kk["status"] == "ok") & (kk["status_rec"] == "ok")].copy()
    kk["k"] = (kk["chi_SNVTA_raw_rec"] - kk["chi_SNVTA_raw"]) / ADULT
    kf = kk[kk["sess_key"].isin(first_keys)]
    f = lambda z: [int(len(z)), float(z["k"].mean()), float(z["k"].std(ddof=1)), float(z["k"].median())]
    excl = sorted(set(first_keys) - set(kf["sess_key"]))
    return dict(all=f(kk), first_usable=f(kf), first_usable_excluded=excl)


def adult_test(beta, se, df, k):
    expected = ADULT * k
    t = (beta - expected) / se
    return dict(k=k, expected_ppb=expected, t=t, p_one_sided=float(stats.t.sf(t, df)))


# ----------------------------------------------------------------------------------------------- main
def main():
    key = pd.read_csv(KEY, dtype={"code": str})
    key["ses_order"] = key["ses_order"].astype(int)
    mer = code_vars(pd.read_csv(MERGED, dtype={"code": str}))
    first = mer[mer["first_usable"] == 1].copy()
    assert len(first) == 80 and first["participant_id"].is_unique
    first_keys = set(first["sess_key"])
    out = {}

    p = first.copy(); p["chi"] = p["chi_SNVTA"]
    assert p["chi"].notna().all()
    out["primary"] = run_m1(p)
    out["primary"]["statsmodels"] = statsmodels_crosscheck(p)
    out["primary"]["adult_test_given"] = adult_test(4.9628, 3.1251, 68, 0.4728)
    out["primary"]["adult_test_recomputed"] = None  # filled after primary k

    # primary recovery k (recovery/raw, read-only) for the README's "0.47 (0.05)"
    pb = read_variant(os.path.join(REC_PRIMARY_DIR, "roi_qsmxt_primary.csv"), key)
    pr = read_variant(os.path.join(REC_PRIMARY_DIR, "roi_qsmxt_rec_m14.csv"), key)
    out["primary"]["recovery_k_L14"] = recovery_k(pb, pr, first_keys)
    kp = out["primary"]["recovery_k_L14"]["first_usable"][1]
    out["primary"]["adult_test_recomputed"] = adult_test(out["primary"]["beta"], out["primary"]["se"],
                                                         out["primary"]["df"], kp)

    roi_ratio_cols = ["chi_SNVTA", "chi_GP", "chi_PUT", "chi_CAU", "chi_RN", "chi_DGM", "chi_STR", "chi_THAL"]
    for inv in ("rts", "tv"):
        e2 = read_variant(os.path.join(E2, "raw", f"roi_qsmxt_e2_lap_vsharp_{inv}.csv"), key)
        rec = read_variant(os.path.join(E2, "raw", f"roi_qsmxt_e2rec_m14_{inv}.csv"), key)
        keep = ["sess_key", "status", "chi_SNVTA", "chi_SNVTA_raw", "chi_GP", "fail_SNVTA", "fimp_SNVTA",
                "inmask_SNVTA", "roi_vol", "Q_phase"] + [c for c in roi_ratio_cols if c not in ("chi_SNVTA", "chi_GP")]
        m = first.merge(e2[keep].add_prefix("e2_").rename(columns={"e2_sess_key": "sess_key"}), on="sess_key",
                        how="left", validate="one_to_one")
        ok = np.isfinite(m["e2_chi_SNVTA"].astype(float))
        a, b = m.loc[ok, "chi_SNVTA"].values, m.loc[ok, "e2_chi_SNVTA"].values
        gp = m[["chi_GP", "e2_chi_GP"]].dropna()
        ag = dict(n_first_usable=int(len(m)), n_valid_e2=int(ok.sum()), icc31=icc31(a, b),
                  pearson=float(stats.pearsonr(a, b)[0]), spearman=float(stats.spearmanr(a, b)[0]),
                  mean_primary=float(a.mean()), mean_e2=float(b.mean()), sd_primary=float(a.std(ddof=1)),
                  sd_e2=float(b.std(ddof=1)), bland_altman_mean_diff=float((b - a).mean()),
                  bland_altman_sd_diff=float((b - a).std(ddof=1)), gp_pearson=float(stats.pearsonr(gp.iloc[:, 0], gp.iloc[:, 1])[0]),
                  gp_n=int(len(gp)))
        dropped = m.loc[~ok, ["sess_key", "e2_status", "e2_chi_SNVTA", "e2_fail_SNVTA"]].to_dict("records")
        kept_fail = m.loc[ok & (m["e2_fail_SNVTA"] == 1), ["sess_key", "e2_fimp_SNVTA", "e2_inmask_SNVTA"]]
        # contrast ratios e2/primary per ROI (first-usable, both finite)
        ratios = {}
        for c in roi_ratio_cols:
            pc, ec = m[c] if c in m else None, m["e2_" + c]
            if pc is None:
                continue
            z = pd.concat([pc, ec], axis=1).dropna()
            if len(z):
                ratios[c] = dict(n=int(len(z)), mean_primary=float(z.iloc[:, 0].mean()), mean_e2=float(z.iloc[:, 1].mean()),
                                 ratio=float(z.iloc[:, 1].mean() / z.iloc[:, 0].mean()))
        # M1 on engine-2 chi; covariates from the merged (primary) table, as the task specifies
        d2 = m[ok].copy(); d2["chi"] = d2["e2_chi_SNVTA"].astype(float)
        r = run_m1(d2)
        r["statsmodels"] = statsmodels_crosscheck(d2)
        # sensitivity 1 (informational): engine-2's own roi_vol / Q_phase columns as covariates
        d3 = d2.copy(); d3["roi_vol"] = d3["e2_roi_vol"]; d3["Q_phase"] = d3["e2_Q_phase"]
        r_owncov = run_m1(d3)
        # sensitivity 2 (informational): drop sessions whose engine-2 SN-VTA fails the registered ROI plausibility rules
        d4 = d2[d2["e2_fail_SNVTA"] != 1].copy()
        r_plaus = run_m1(d4)
        kres = recovery_k(e2, rec, first_keys)
        out[f"e2_{inv}"] = dict(
            agreement=ag, dropped_rows=dropped,
            n_kept_with_e2_fail_SNVTA=int(len(kept_fail)), kept_fail_keys=kept_fail["sess_key"].tolist(),
            kept_fail_fimp_range=[float(kept_fail["e2_fimp_SNVTA"].min()), float(kept_fail["e2_fimp_SNVTA"].max())] if len(kept_fail) else None,
            kept_fail_inmask_range=[float(kept_fail["e2_inmask_SNVTA"].min()), float(kept_fail["e2_inmask_SNVTA"].max())] if len(kept_fail) else None,
            covariate_identity=dict(
                roi_vol_max_abs_diff=float(np.nanmax(np.abs(d2["roi_vol"] - d2["e2_roi_vol"]))),
                Q_phase_max_abs_diff=float(np.nanmax(np.abs(d2["Q_phase"] - d2["e2_Q_phase"])))),
            contrast_ratio=ratios, M1=r, recovery_k_L14=kres,
            adult_test=adult_test(r["beta"], r["se"], r["df"], kres["first_usable"][1]),
            sens_e2_own_covariates=dict(beta=r_owncov["beta"], se=r_owncov["se"], p_two=r_owncov["p_two"],
                                        pB=r_owncov["pB"], outcome=r_owncov["outcome"]),
            sens_drop_e2_plausibility_fail=dict(n=r_plaus["n"], beta=r_plaus["beta"], se=r_plaus["se"],
                                                p_two=r_plaus["p_two"], ci95=r_plaus["ci95"], delta=r_plaus["delta"],
                                                sigma_B=r_plaus["sigma_B"], pB=r_plaus["pB"], outcome=r_plaus["outcome"]))
        # agreement CSV check
        csv = pd.read_csv(os.path.join(E2, f"agreement_first_usable_{inv}.csv"))
        cm = csv.merge(m[["sess_key", "chi_SNVTA", "e2_chi_SNVTA", "e2_status"]], on="sess_key", suffixes=("_csv", ""))
        out[f"e2_{inv}"]["agreement_csv_check"] = dict(
            n_rows=int(len(csv)), n_matched=int(len(cm)),
            max_abs_diff_primary=float(np.nanmax(np.abs(cm["chi_SNVTA_csv"] - cm["chi_SNVTA"]))),
            max_abs_diff_e2=float(np.nanmax(np.abs(cm["e2_chi_SNVTA_csv"] - cm["e2_chi_SNVTA"]))),
            na_pattern_equal=bool((cm["e2_chi_SNVTA_csv"].isna() == cm["e2_chi_SNVTA"].isna()).all()),
            status_equal=bool((cm["e2_status_csv"].astype(str) == cm["e2_status"].astype(str)).all()))

    def conv(o):
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.bool_,)):
            return bool(o)
        raise TypeError(type(o))

    with open(os.path.join(HERE, "recomputed.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=conv)
    return out


if __name__ == "__main__":
    main()
