"""LEAD-chi -- independent verification library (POST HOC check of results_posthoc/validation_local).

Written from scratch in Python/numpy/scipy. It does NOT source or call the R code under check, and it does not use
statsmodels for the main fits (statsmodels is used only as an optional second check of the mixed models).
Every routine below re-implements a published algorithm:
  * natural cubic spline basis: truncated-power form (Hastie, Tibshirani & Friedman, ESL 2nd ed., eq. 5.4-5.5) with
    boundary knots at range(x) and interior knots at the 1/3 and 2/3 type-7 quantiles. Together with the intercept
    it spans exactly the column space of R's [1, splines::ns(x, df = 3)], so every PS_e coefficient and its HC3 SE
    are identical even though the basis columns differ;
  * OLS via QR, HC3 sandwich (MacKinnon & White 1985), t reference with df = n - p;
  * Huber M-estimator: a line-by-line port of the IRLS algorithm of MASS::rlm.default (psi.huber, k = 1.345,
    "proposal 2" scale, init = LS, acc = 1e-4 on the residual change, maxit = 200);
  * random-intercept linear mixed model by profiled REML (closed-form per-group inverses), nlme-style
    "containment" denominator df;
  * Fisher exact test for r x c tables by full enumeration with the R tolerance (1 + 1e-7);
  * Wilcoxon rank-sum (exact null distribution without ties; normal approximation with continuity and tie
    correction otherwise, as in R's wilcox.test);
  * Benjamini-Hochberg q values.
Nothing here writes files.
"""
import itertools
import math

import numpy as np
import pandas as pd
from scipy import optimize, stats


# ----------------------------------------------------------------------------------------------------------------
# design pieces
# ----------------------------------------------------------------------------------------------------------------
def q7(x, p):
    """Hyndman-Fan type-7 quantile (R default; numpy 'linear')."""
    return np.quantile(np.asarray(x, float), p, method="linear")


def ns3_basis(x):
    """3 columns spanning (with an intercept) the natural cubic splines with knots
    {min(x), q7(1/3), q7(2/3), max(x)}, i.e. the same space as R's splines::ns(x, df = 3)."""
    x = np.asarray(x, float)
    K = [x.min(), *q7(x, [1 / 3, 2 / 3]), x.max()]
    xK = K[-1]

    def dk(k):
        return (np.maximum(x - K[k], 0.0) ** 3 - np.maximum(x - xK, 0.0) ** 3) / (xK - K[k])

    cols = [x] + [dk(k) - dk(len(K) - 2) for k in range(len(K) - 2)]
    return np.column_stack(cols)


def zscore(v):
    v = np.asarray(v, float)
    return (v - v.mean()) / v.std(ddof=1)


def race_dummies(race3):
    """Indicator columns for the race3 levels present (first present level of Black/White/Other omitted).
    Spans the same space (with the intercept) as R's contr.sum coding over the present levels."""
    lev = [l for l in ("Black", "White", "Other") if l in set(race3)]
    return {f"race_{l}": (np.asarray(race3) == l).astype(float) for l in lev[1:]}


def m1_design(d, y, extra_z=(), extra_raw=None, with_interaction=True):
    """Model matrix of M1: PS_e*ADHD_e + ns(age,3) + sex_e + race3 + z(roi_vol) + z(Q_phase) [+ z(extra_z)] [+ extra_raw].
    `d` must already be the complete-case analysis sample. Returns (X, names, yvec)."""
    cols = {"(Intercept)": np.ones(len(d)), "PS_e": d["PS_e"].to_numpy(float), "ADHD_e": d["ADHD_e"].to_numpy(float)}
    B = ns3_basis(d["age"].to_numpy(float))
    for j in range(B.shape[1]):
        cols[f"ns_age_{j + 1}"] = B[:, j]
    cols["sex_e"] = d["sex_e"].to_numpy(float)
    cols.update(race_dummies(d["race3"].to_numpy()))
    cols["z_roi_vol"] = zscore(d["roi_vol"])
    cols["z_Q_phase"] = zscore(d["Q_phase"])
    for v in extra_z:
        cols[f"z_{v}"] = zscore(d[v])
    if extra_raw is not None:
        for k, v in extra_raw.items():
            cols[k] = np.asarray(v, float)
    if with_interaction:
        cols["PS_e:ADHD_e"] = cols["PS_e"] * cols["ADHD_e"]
    names = list(cols)
    X = np.column_stack([cols[n] for n in names])
    return X, names, d[y].to_numpy(float)


# ----------------------------------------------------------------------------------------------------------------
# OLS + HC3
# ----------------------------------------------------------------------------------------------------------------
class OLS:
    def __init__(self, X, y, names):
        X = np.asarray(X, float)
        y = np.asarray(y, float)
        n, p = X.shape
        if np.linalg.matrix_rank(X) < p:
            raise ValueError("rank-deficient design")
        Q, R = np.linalg.qr(X)
        self.b = np.linalg.solve(R, Q.T @ y)
        Ri = np.linalg.inv(R)
        self.XtXi = Ri @ Ri.T
        self.e = y - X @ self.b
        self.h = np.sum(Q ** 2, axis=1)
        u = self.e / (1.0 - self.h)
        meat = (X * (u ** 2)[:, None]).T @ X
        self.V = self.XtXi @ meat @ self.XtXi
        self.se = np.sqrt(np.diag(self.V))
        self.n, self.p, self.df = n, p, n - p
        self.sigma = math.sqrt(float(self.e @ self.e) / self.df)
        self.names = list(names)
        self.X, self.y = X, y

    def idx(self, term):
        return self.names.index(term)

    def coef(self, term):
        j = self.idx(term)
        b, se = float(self.b[j]), float(self.se[j])
        return self._summ(b, se)

    def lin(self, cdict):
        c = np.zeros(self.p)
        for k, v in cdict.items():
            c[self.idx(k)] = v
        est = float(c @ self.b)
        se = float(math.sqrt(c @ self.V @ c))
        return self._summ(est, se)

    def _summ(self, est, se):
        t = est / se
        q = stats.t.ppf(0.975, self.df)
        return dict(est=est, se=se, t=t, df=self.df, p=2 * stats.t.sf(abs(t), self.df),
                    p_upper=stats.t.sf(t, self.df), p_lower=stats.t.cdf(t, self.df), lo=est - q * se, hi=est + q * se)

    def rstudent(self):
        r, h = self.e, self.h
        s2i = ((self.df) * self.sigma ** 2 - r ** 2 / (1 - h)) / (self.df - 1)
        return r / np.sqrt(s2i * (1 - h))


def bound_test(est, se, df, bound):
    """one-sided test of beta <= bound: t_B = (est - bound)/se, upper tail."""
    tB = (est - bound) / se
    return tB, stats.t.sf(tB, df)


# ----------------------------------------------------------------------------------------------------------------
# Huber M-estimator: port of MASS::rlm.default (method "M", psi.huber, scale.est "proposal 2", init "ls")
# ----------------------------------------------------------------------------------------------------------------
def rlm_huber(X, y, k=1.345, k2=1.345, maxit=200, acc=1e-4):
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    n, p = X.shape
    coef = np.linalg.lstsq(X, y, rcond=None)[0]
    resid = y - X @ coef
    n1 = n - p
    theta = 2 * stats.norm.cdf(k2) - 1
    gamma = theta + k2 ** 2 * (1 - theta) - 2 * k2 * stats.norm.pdf(k2)
    scale = 1.4826 * np.median(np.abs(resid))          # mad(resid, center = 0)
    done = False
    for _ in range(maxit):
        testpv = resid.copy()
        scale = math.sqrt(np.sum(np.minimum(resid ** 2, (k2 * scale) ** 2)) / (n1 * gamma))
        if scale == 0:
            done = True
            break
        u = np.abs(resid / scale)
        with np.errstate(divide="ignore"):
            w = np.minimum(1.0, k / u)
        sw = np.sqrt(w)
        coef = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)[0]
        resid = y - X @ coef
        conv = math.sqrt(np.sum((testpv - resid) ** 2) / max(1e-20, np.sum(testpv ** 2)))
        if conv <= acc:
            done = True
            break
    return dict(coef=coef, scale=scale, converged=done)


# ----------------------------------------------------------------------------------------------------------------
# random-intercept LMM by profiled REML
# ----------------------------------------------------------------------------------------------------------------
class LMMRI:
    """y = X b + u_g + e, u_g ~ N(0, s2*theta), e ~ N(0, s2). REML, theta profiled numerically (bounded scalar search
    on log theta, refined), s2 and b in closed form. Cov(b) = s2 * (X' H^-1 X)^-1 (as nlme's varFix)."""

    def __init__(self, X, y, groups, names):
        X = np.asarray(X, float)
        y = np.asarray(y, float)
        g = pd.factorize(pd.Series(groups))[0]
        self.X, self.y, self.g, self.names = X, y, g, list(names)
        self.N, self.p = X.shape
        self.G = g.max() + 1
        self.ng = np.bincount(g)
        # group sums
        self.Sx = np.zeros((self.G, self.p))
        np.add.at(self.Sx, g, X)
        self.Sy = np.bincount(g, weights=y)
        self.XtX, self.Xty, self.yty = X.T @ X, X.T @ y, y @ y
        res = optimize.minimize_scalar(lambda lt: self._nll(math.exp(lt)), bounds=(-25, 10), method="bounded",
                                       options=dict(xatol=1e-12, maxiter=2000))
        lt = res.x
        # polish with a fine golden search around the optimum
        res2 = optimize.minimize_scalar(lambda lt: self._nll(math.exp(lt)), bracket=(lt - 0.1, lt, lt + 0.1),
                                        method="brent", tol=1e-14) if lt > -24 else res
        self.theta = math.exp(res2.x) if res2.fun <= res.fun else math.exp(lt)
        self.theta0_nll = self._nll(0.0)
        if self.theta0_nll < self._nll(self.theta):
            self.theta = 0.0
        self._fit(self.theta)

    def _parts(self, theta):
        c = theta / (1.0 + self.ng * theta)
        XHX = self.XtX - (self.Sx * c[:, None]).T @ self.Sx
        XHy = self.Xty - (self.Sx * (c * self.Sy)[:, None]).sum(0)
        yHy = self.yty - np.sum(c * self.Sy ** 2)
        return XHX, XHy, yHy, np.sum(np.log1p(self.ng * theta))

    def _nll(self, theta):
        XHX, XHy, yHy, ldH = self._parts(theta)
        b = np.linalg.solve(XHX, XHy)
        rHr = yHy - XHy @ b
        s2 = rHr / (self.N - self.p)
        sign, ld = np.linalg.slogdet(XHX)
        return 0.5 * ((self.N - self.p) * math.log(s2) + ldH + ld)

    def _fit(self, theta):
        XHX, XHy, yHy, _ = self._parts(theta)
        self.b = np.linalg.solve(XHX, XHy)
        self.s2 = (yHy - XHy @ self.b) / (self.N - self.p)
        self.V = self.s2 * np.linalg.inv(XHX)
        self.se = np.sqrt(np.diag(self.V))

    def within_varying(self, j):
        col = self.X[:, j]
        df = pd.DataFrame(dict(g=self.g, v=col))
        return bool((df.groupby("g")["v"].nunique() > 1).any())

    def containment_df(self):
        """nlme rule for one grouping level: terms that vary within groups get N - (G + p_inner);
        group-level terms get G - (1 + p_outer); the intercept gets the inner df."""
        inner = [j for j in range(self.p) if self.names[j] != "(Intercept)" and self.within_varying(j)]
        outer = [j for j in range(self.p) if self.names[j] != "(Intercept)" and j not in inner]
        dfi = self.N - (self.G + len(inner))
        dfo = self.G - (1 + len(outer))
        return {self.names[j]: (dfi if (j in inner or self.names[j] == "(Intercept)") else dfo) for j in range(self.p)}

    def coef(self, term):
        j = self.names.index(term)
        df = self.containment_df()[term]
        est, se = float(self.b[j]), float(self.se[j])
        t = est / se
        q = stats.t.ppf(0.975, df)
        return dict(est=est, se=se, t=t, df=df, p=2 * stats.t.sf(abs(t), df), p_upper=stats.t.sf(t, df),
                    lo=est - q * se, hi=est + q * se)


# ----------------------------------------------------------------------------------------------------------------
# tests
# ----------------------------------------------------------------------------------------------------------------
def bh(p):
    p = np.asarray(p, float)
    m = len(p)
    o = np.argsort(p)[::-1]
    q = np.minimum.accumulate(p[o] * m / np.arange(m, 0, -1))
    out = np.empty(m)
    out[o] = np.minimum(q, 1.0)
    return out


def _log_table_prob(tab, rs, cs, n):
    return (sum(math.lgamma(r + 1) for r in rs) + sum(math.lgamma(c + 1) for c in cs) - math.lgamma(n + 1)
            - sum(math.lgamma(x + 1) for x in np.ravel(tab)))


def fisher_rxc(tab):
    """Exact Fisher test for an r x c table: sum of probabilities of all tables with the observed margins whose
    probability is <= the observed one * (1 + 1e-7) (R's fisher.test convention)."""
    tab = np.asarray(tab, int)
    rs, cs, n = tab.sum(1), tab.sum(0), int(tab.sum())
    r, c = tab.shape
    lp_obs = _log_table_prob(tab, rs, cs, n)
    const = sum(math.lgamma(x + 1) for x in rs) + sum(math.lgamma(x + 1) for x in cs) - math.lgamma(n + 1)
    lf = np.array([math.lgamma(i + 1) for i in range(n + 1)])

    # enumerate column by column; vectorise over the final two columns
    def col_splits(total, caps):
        # all vectors of length r, non-negative ints, sum = total, element-wise <= caps
        if len(caps) == 1:
            if total <= caps[0]:
                yield (total,)
            return
        for a in range(min(total, caps[0]) + 1):
            for rest in col_splits(total - a, caps[1:]):
                yield (a,) + rest

    thr = lp_obs + math.log1p(1e-7)
    tot = 0.0
    split_cache = {}

    def splits_arr(total, caps):
        key = (total, caps)
        if key not in split_cache:
            a = np.array(list(col_splits(total, caps)), dtype=int).reshape(-1, r)
            split_cache[key] = (a, lf[a].sum(1) if len(a) else np.zeros(0))
        return split_cache[key]

    def rec(ci, rem, acc_lf):
        nonlocal tot
        if ci == c - 2:
            # second-last column enumerated as an array, last column forced
            S, Sl = splits_arr(int(cs[ci]), rem)
            if len(S) == 0:
                return
            last = np.asarray(rem)[None, :] - S
            ok = np.all(last >= 0, axis=1)
            lp = const - acc_lf - Sl[ok] - lf[last[ok]].sum(1)
            tot += float(np.exp(lp[lp <= thr]).sum())
            return
        for s in col_splits(int(cs[ci]), rem):
            rec(ci + 1, tuple(rr - ss for rr, ss in zip(rem, s)), acc_lf + sum(lf[x] for x in s))

    if c == 1:
        return 1.0
    rec(0, tuple(int(x) for x in rs), 0.0)
    return min(1.0, tot)


def wilcox_rs(x, y):
    """Two-sided Wilcoxon rank-sum p as R's wilcox.test(x, y, exact = TRUE, correct = TRUE)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    nx, ny = len(x), len(y)
    allv = np.concatenate([x, y])
    rk = stats.rankdata(allv)
    W = rk[:nx].sum() - nx * (nx + 1) / 2
    ties = len(np.unique(allv)) < len(allv)
    if not ties and nx < 50 and ny < 50:
        # exact null distribution of W by counting subsets (DP over rank sums)
        N = nx + ny
        maxs = nx * ny
        # Mann-Whitney counting recursion f(m, n, w) = f(m-1, n, w-n) + f(m, n-1, w)
        f = {}

        def F(m, nn, w):
            if w < 0 or w > m * nn:
                return 0.0
            if m == 0 or nn == 0:
                return 1.0 if w == 0 else 0.0
            key = (m, nn, w)
            if key not in f:
                f[key] = F(m - 1, nn, w - nn) + F(m, nn - 1, w)
            return f[key]

        tot = math.comb(N, nx)
        cdf = sum(F(nx, ny, w) for w in range(int(W) + 1)) / tot
        sf = sum(F(nx, ny, w) for w in range(int(W), maxs + 1)) / tot
        return min(1.0, 2 * min(cdf, sf)), W, "exact"
    _, tc = np.unique(allv, return_counts=True)
    z = W - nx * ny / 2
    sigma = math.sqrt((nx * ny / 12) * ((nx + ny + 1) - np.sum(tc ** 3 - tc) / ((nx + ny) * (nx + ny - 1))))
    if sigma == 0:
        return float("nan"), W, "all tied"
    corr = 0.5 * np.sign(z)
    zz = (z - corr) / sigma
    return 2 * min(stats.norm.cdf(zz), stats.norm.sf(zz)), W, "normal+cc"


def spearman(x, y):
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    r = np.corrcoef(rx, ry)[0, 1]
    n = len(x)
    t = r / math.sqrt((1 - r ** 2) / (n - 2))
    return r, 2 * stats.t.sf(abs(t), n - 2)
