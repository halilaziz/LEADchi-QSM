#!/usr/bin/env python3
"""Independent re-computation of the LEAD-chi POST HOC SN-VTA recovery analysis.

Written from scratch (numpy / pandas / scipy / statsmodels); the R script in ../code/ is not run or copied.
Read-only on all inputs; writes only into this verification/ folder.
"""
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
REC = HERE.parent
UNL = Path("/path/to/Penn_LEAD/LEADchi/mirror/unlocked")
RL = 0.938278557806515
DELTA_B = 0.6 * math.sqrt(RL)
HUBER_K = 1.345
LEVELS = [("rec_m05", -5.0), ("rec_m10", -10.0), ("rec_m14", -14.24), ("rec_m20", -20.0), ("rec_p10", 10.0)]


# ----------------------------------------------------------------------------------------------- design helpers
def nat_spline_basis(x, df=3):
    """Natural cubic spline basis without intercept (truncated-power form, ESL 5.4-5.5).
    Interior knots at quantiles 1/3, 2/3 (type-7), boundary knots at min / max: same column space as R splines::ns."""
    x = np.asarray(x, float)
    probs = np.linspace(0, 1, df + 1)[1:-1]
    inner = np.quantile(x, probs)  # numpy 'linear' == R type 7
    knots = np.concatenate([[x.min()], inner, [x.max()]])
    K = len(knots)

    def dk(k):
        return (np.clip(x - knots[k], 0, None) ** 3 - np.clip(x - knots[K - 1], 0, None) ** 3) / (knots[K - 1] - knots[k])

    cols = [x]
    for k in range(K - 2):
        cols.append(dk(k) - dk(K - 2))
    B = np.column_stack(cols)
    return B, knots


def design(d, with_qc=True):
    """M1 design: 1, PS_e, ADHD_e, PS_e:ADHD_e, ns(age,3), sex_e, race3 (contr.sum over levels present), z(roi_vol), z(Q_phase)."""
    n = len(d)
    cols, names = [np.ones(n)], ["(Intercept)"]
    cols += [d.PS_e.values, d.ADHD_e.values, (d.PS_e * d.ADHD_e).values]
    names += ["PS_e", "ADHD_e", "PS_e:ADHD_e"]
    B, knots = nat_spline_basis(d.age.values, 3)
    cols += [B[:, j] for j in range(B.shape[1])]
    names += [f"ns{j+1}" for j in range(B.shape[1])]
    cols.append(d.sex_e.values); names.append("sex_e")
    lev = [l for l in ["Black", "White", "Other"] if l in set(d.race3_chr)]
    if len(lev) >= 2:
        for j, l in enumerate(lev[:-1]):
            c = np.where(d.race3_chr == l, 1.0, np.where(d.race3_chr == lev[-1], -1.0, 0.0))
            cols.append(c); names.append(f"race3_{j+1}")
    if with_qc:
        for v in ["roi_vol", "Q_phase"]:
            z = (d[v] - d[v].mean()) / d[v].std(ddof=1)
            cols.append(z.values); names.append(f"z_{v}")
    X = np.column_stack(cols)
    return X, names, knots


def ols_hc3(X, y):
    n, p = X.shape
    XtXi = np.linalg.inv(X.T @ X)
    b = XtXi @ X.T @ y
    e = y - X @ b
    h = np.einsum("ij,jk,ik->i", X, XtXi, X)
    meat = (X * (e / (1 - h))[:, None] ** 2).T @ X
    V = XtXi @ meat @ XtXi
    df = n - p
    sigma = math.sqrt(e @ e / df)
    return dict(b=b, se=np.sqrt(np.diag(V)), e=e, df=df, sigma=sigma, n=n, p=p)


def term(m, j):
    b, se, df = m["b"][j], m["se"][j], m["df"]
    t = b / se
    q975, q95 = stats.t.ppf(0.975, df), stats.t.ppf(0.95, df)
    return dict(est=b, se=se, t=t, df=df, p_two=2 * stats.t.sf(abs(t), df), p_lower=stats.t.cdf(t, df),
                p_upper=stats.t.sf(t, df), ci95=[b - q975 * se, b + q975 * se], ci90=[b - q95 * se, b + q95 * se])


def huber_prop2_scale(X, y, k=HUBER_K, maxit=200, acc=1e-4):
    """IRLS M-estimation with Huber psi and Huber's proposal-2 scale, written to the published MASS::rlm algorithm
    (init = LS, initial scale = 1.4826*median|r|, scale updated at the start of each iteration, convergence on the
    relative change of the residual vector <= 1e-4). Returns (scale, converged, iterations)."""
    n, p = X.shape
    b = np.linalg.lstsq(X, y, rcond=None)[0]
    r = y - X @ b
    s = 1.4826 * np.median(np.abs(r))
    theta = 2 * stats.norm.cdf(k) - 1
    gamma = theta + k ** 2 * (1 - theta) - 2 * k * stats.norm.pdf(k)
    n1 = n - p
    for it in range(1, maxit + 1):
        r_old = r.copy()
        s = math.sqrt(np.sum(np.minimum(r ** 2, (k * s) ** 2)) / (n1 * gamma))
        u = np.abs(r / s)
        w = np.where(u <= k, 1.0, k / np.maximum(u, 1e-300))
        sw = np.sqrt(w)
        b = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)[0]
        r = y - X @ b
        conv = math.sqrt(np.sum((r_old - r) ** 2) / max(1e-20, np.sum(r_old ** 2)))
        if conv <= acc:
            return s, True, it
    return s, False, maxit


def m1_decide(d, y="chi_SNVTA"):
    X, names, _ = design(d, with_qc=True)
    yy = d[y].values.astype(float)
    m = ols_hc3(X, yy)
    tt = term(m, names.index("PS_e"))
    s_rob, conv, it = huber_prop2_scale(X, yy)
    if not conv:
        e = m["e"]
        s_rob = 1.4826 * np.median(np.abs(e - np.median(e)))
    sigma_B = min(m["sigma"], s_rob)
    tB = (tt["est"] + DELTA_B * sigma_B) / tt["se"]
    pB = stats.t.sf(tB, m["df"])
    A = tt["est"] < 0 and tt["p_two"] < 0.05
    B = pB < 0.05
    outcome = "A+B" if (A and B) else ("A" if A else ("B" if B else "C"))
    return dict(beta=tt["est"], se=tt["se"], t=tt["t"], df=m["df"], p_two=tt["p_two"], ci95=tt["ci95"], ci90=tt["ci90"],
                sigma_resid=m["sigma"], sigma_rob=s_rob, huber_converged=conv, huber_iter=it, sigma_B=sigma_B,
                bound=-DELTA_B * sigma_B, tB=tB, pB=pB, outcome=outcome,
                disconfirming=bool(tt["est"] > 0 and tt["p_two"] < 0.05), n=m["n"], p=m["p"])


# ----------------------------------------------------------------------------------------------- Spearman (R cor.test)
def prho_upper_edgeworth(n, s):
    """AS 89 Edgeworth approximation of Pr[S >= s] (n > 9), as used by R cor.test(method='spearman', exact=TRUE)."""
    c = [0.2274, 0.2531, 0.1745, 0.0758, 0.1033, 0.3932, 0.0879, 0.0151, 0.0072, 0.0831, 0.0131, 4.6e-4]
    b = 1.0 / n
    x = (6 * (s - 1) * b / (n * n - 1) - 1) * math.sqrt(1 / b - 1)
    y = x * x
    u = x * b * (c[0] + b * (c[1] + c[2] * b) + y * (-c[3] + b * (c[4] + c[5] * b)
                 - y * b * (c[6] + c[7] * b - y * (c[8] - c[9] * b + y * b * (c[10] - c[11] * y)))))
    return min(1.0, max(0.0, u / math.exp(y / 2) + stats.norm.sf(x)))


def spearman_r(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    r = np.corrcoef(rx, ry)[0, 1]
    ties = min(len(np.unique(x)), len(np.unique(y))) < n
    if not ties and n < 1290:
        q = (n ** 3 - n) * (1 - r) / 6
        n3 = (n ** 3 - n) / 3
        if q > (n ** 3 - n) / 6:
            p = prho_upper_edgeworth(n, round(q)) if round(q) <= n3 else 0.0
        else:
            p = 1 - prho_upper_edgeworth(n, round(q) + 2)
        p = min(2 * p, 1.0)
        method = "AS89 Edgeworth (R exact branch)"
    else:
        t = r / math.sqrt((1 - r ** 2) / (n - 2))
        p = 2 * min(stats.t.cdf(t, n - 2), stats.t.sf(t, n - 2))
        method = "t approximation (ties)"
    p_tapprox = 2 * stats.t.sf(abs(r / math.sqrt((1 - r ** 2) / (n - 2))), n - 2)
    return r, p, method, ties, p_tapprox


def summ(x):
    x = np.asarray(x, float)
    return [len(x), x.mean(), x.std(ddof=1), np.median(x), np.quantile(x, .25), np.quantile(x, .75), x.min(), x.max()]


# ----------------------------------------------------------------------------------------------- data
def load():
    key = pd.read_csv(UNL / "key_plain.csv")
    mer = pd.read_csv(UNL / "merged_qsmxt_primary_SNVTA.csv")
    mer["PS_e"] = np.where(mer.study_group == "PRO/CHR", 0.5, -0.5)
    mer["ADHD_e"] = np.where(mer.dx_adhd == 1, 0.5, -0.5)
    mer["sex_e"] = np.where(mer.sex == "M", 0.5, -0.5)
    mer["race3_chr"] = np.where(mer.race == "Black or African American", "Black", np.where(mer.race == "White", "White", "Other"))
    mer["cell"] = np.where(mer.PS_e > 0, np.where(mer.ADHD_e > 0, "PS+ADHD", "PS-only"), np.where(mer.ADHD_e > 0, "ADHD-only", "TD"))
    mer["sess_key"] = mer.participant_id + "_" + mer.session_id
    cols = ["code", "ses_order", "status", "chi_SNVTA_raw", "ref_ppb", "chi_GP_raw", "chi_THAL_raw"]
    prim = pd.read_csv(REC / "raw" / "roi_qsmxt_primary.csv")[cols]
    prim = prim.merge(key[["code", "ses_order", "participant_id", "session_id"]], on=["code", "ses_order"], how="inner")
    rows = []
    for v, L in LEVELS:
        r = pd.read_csv(REC / "raw" / f"roi_qsmxt_{v}.csv")[cols]
        m = prim.merge(r, on=["code", "ses_order"], suffixes=("", "_rec"))
        m["L"], m["variant"] = L, v
        m["d"] = m.chi_SNVTA_raw_rec - m.chi_SNVTA_raw
        m["k"] = m.d / L
        m["d_ref"] = m.ref_ppb_rec - m.ref_ppb
        m["d_GP"] = m.chi_GP_raw_rec - m.chi_GP_raw
        m["d_THAL"] = m.chi_THAL_raw_rec - m.chi_THAL_raw
        rows.append(m)
    long = pd.concat(rows, ignore_index=True)
    long = long[(long.status == "ok") & (long.status_rec == "ok") & np.isfinite(long.d)].copy()
    long["sess_key"] = long.participant_id + "_" + long.session_id
    long = long.merge(mer[["sess_key", "first_usable", "PS_e", "ADHD_e", "cell", "age", "dice_min", "motion_gsr",
                           "roi_vol", "Q_phase", "sex_e", "race3_chr"]], on="sess_key", how="inner")
    return key, mer, long


def main():
    key, mer, long = load()
    out = {}
    lv = [L for _, L in LEVELS]
    out["n_long"] = len(long)
    out["k_all"] = {L: summ(long.k[long.L == L]) for L in lv}
    fu = long[long.first_usable == 1]
    out["k_first_usable"] = {L: summ(fu.k[fu.L == L]) for L in lv}
    out["d_ref_ppb"] = {L: summ(long.d_ref[long.L == L]) for L in lv}
    out["d_GP_ppb"] = {L: summ(long.d_GP[long.L == L]) for L in lv}
    out["d_THAL_ppb"] = {L: summ(long.d_THAL[long.L == L]) for L in lv}
    # per-session slopes through origin
    sl = []
    for _, s in long.groupby("sess_key"):
        if len(s) >= 3:
            sl.append(np.sum(s.L * s.d) / np.sum(s.L ** 2))
    out["k_linear_per_session"] = summ(sl)
    # pooled
    res = stats.linregress(long.L, long.d)
    out["pooled_linear"] = dict(slope=res.slope, intercept=res.intercept, r2=res.rvalue ** 2, n=len(long))
    # cross-check the pooled fit with statsmodels
    import statsmodels.api as sm
    f = sm.OLS(long.d.values, sm.add_constant(long.L.values)).fit()
    out["pooled_linear_sm"] = dict(slope=f.params[1], intercept=f.params[0], r2=f.rsquared)
    # k by group at L = -14.24 (first usable)
    k14 = fu[fu.L == -14.24].copy()
    X, names, knots = design(k14, with_qc=False)
    mk = ols_hc3(X, k14.k.values)
    out["k_by_group_L14"] = dict(PS=term(mk, names.index("PS_e")), ADHD=term(mk, names.index("ADHD_e")),
                                 n=len(k14), knots=list(knots),
                                 cell_means=[k14.k[k14.cell == c].mean() for c in sorted(k14.cell.unique())],
                                 cell_order=sorted(k14.cell.unique()))
    # cross-check HC3 with statsmodels
    f2 = sm.OLS(k14.k.values, X).fit(cov_type="HC3")
    out["k_by_group_L14_sm"] = dict(PS_est=f2.params[1], PS_se=f2.bse[1], ADHD_est=f2.params[2], ADHD_se=f2.bse[2])
    sp = {}
    for v in ["dice_min", "motion_gsr", "roi_vol", "Q_phase", "age"]:
        r, p, meth, ties, pt = spearman_r(k14.k, k14[v])
        sp[v] = dict(rho=r, p=p, method=meth, ties=ties, p_tapprox=pt)
    out["k_spearman_L14"] = sp
    # observed M1 + counterfactuals
    first = mer[mer.first_usable == 1].copy()
    obs = m1_decide(first)
    out["observed"] = obs
    cf = []
    for v, L in LEVELS:
        dd = long[long.L == L][["sess_key", "d"]]
        f_ = first.merge(dd, on="sess_key", how="left")
        miss = int(f_.d.isna().sum())
        f_["d"] = f_.d.fillna(dd.d.mean())
        f_["chi_SNVTA"] = f_.chi_SNVTA + np.where(f_.PS_e > 0, f_.d, 0.0)
        r = m1_decide(f_)
        cf.append(dict(L_true_ppb=L, n_missing_rec=miss, mean_recovered_in_PS=f_.d[f_.PS_e > 0].mean(), **r))
    out["counterfactual"] = cf
    kbar = out["k_first_usable"][-14.24][1]
    out["true_scale"] = dict(k_mean_L14=kbar, bound_true_ppb=obs["bound"] / kbar,
                             lower95_one_sided_true_ppb=obs["ci90"][0] / kbar,
                             lower95_two_sided_true_ppb=obs["ci95"][0] / kbar)
    # B threshold on the observed scale: smallest beta that still gives pB < .05 at the observed SE / sigma_B
    out["B_threshold_beta"] = obs["bound"] + stats.t.ppf(0.95, obs["df"]) * obs["se"]
    # D table (README): t = (4.9628 - k*D)/3.1251, df 68, one-sided upper p
    tab = []
    for D in [-5, -10, -14.24, -20]:
        shift = kbar * D
        t_r = (4.9628 - shift) / 3.1251
        t_f = (obs["beta"] - shift) / obs["se"]
        tab.append(dict(D=D, shift=shift, p_rounded_inputs=stats.t.sf(t_r, 68), p_full=stats.t.sf(t_f, obs["df"])))
    out["D_table"] = tab
    # n-check of sessions
    out["n_sessions_per_level"] = {L: int((long.L == L).sum()) for L in lv}
    out["n_first_usable_per_level"] = {L: int((fu.L == L).sum()) for L in lv}
    out["n_first_usable_merged"] = int(len(first))
    out["n_PS_first"] = int((first.PS_e > 0).sum())
    return out, long


if __name__ == "__main__":
    out, long = main()

    def conv(o):
        if isinstance(o, dict):
            return {str(k): conv(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [conv(v) for v in o]
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, np.bool_):
            return bool(o)
        return o

    (HERE / "recomputed.json").write_text(json.dumps(conv(out), indent=1))
    print(json.dumps(conv(out), indent=1)[:200])
