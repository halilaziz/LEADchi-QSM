#!/usr/bin/env python3
"""gates_preunlock.py - every pre-unlock gate, branch and decision constant of LEAD-chi (registration Q24, A.6).

Label-free and age-free: reads only code-keyed ROI tables (hidden), the frozen usable lists, the visible flags,
the QC kappa and the DWI-side PC5/n_min JSON. It must be run through audit_run.py, which rejects any output that
contains a per-code value. Output JSON is aggregate-only with a fixed schema (the analysis code validates it).

Stage `g2` (amendment 1, before the full run): engine agreement on the 20 coded sessions of blind/gate20.txt.
  ICC(3,1) (two-way mixed, consistency; McGraw & Wong 1996, F-based 95% CI) of SN-VTA chi between QSMxT and SEPIA:
  G2 passes if the point estimate >= 0.90 AND the lower 95% bound >= 0.75; the same for chi_DN_R gives G2-DN.
  Co-reported: ICC(2,1) (absolute agreement, McGraw & Wong CI), Pearson r and the Deming (orthogonal) slope.
  If A2 was declined (--a2 no): within-QSMxT algorithm agreement over V-SHARP+MEDI, V-SHARP+RTS and PDF+iLSQR,
  ICC(3,1) point >= 0.75; the engine stays QSMxT whatever the outcome (a failure is a registered caveat).
Stage `final` (amendment 2, after build_usable.py --stage usable): branch order of Q24.
  1. engine from gates/g2.json; 2. PC1, PC2w, PC3 on that engine's primary pipeline (SN-VTA usable list), and, if any
  fails, the same on the PDF+iLSQR fallback (fallback_PC_pass); 3. G3 on the 9-node common-support multiverse table
  (ICC(3,1) point >= 0.75 over the first-usable sessions) selects the ROI (SN-VTA or SN-only); 4. PC4-DN, G2-DN and
  PC5 decide S3's status. Then, for the selected branch: PC4 and r_L (split-echo, Spearman-Brown, Fisher-z lower
  95% limit of r), the pooled label-free constants (mean/SD of chi_primary and chi_STR over the first-usable
  sessions; sd_total = the same SD), the cross-session reliability (Pearson r and ICC(3,1) over the returners'
  first two usable sessions, not age-adjusted; T1 informativeness), SD of dref, rho(SN-VTA, STR),
  r(chi_DN_R, chi_DGM), G4 (qc/g4_final.json), G5 (flags_coded/g5_summary.json) and counts.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lx_common as L  # noqa: E402

BRANCH_VARIANT = {"primary": "primary", "fallback_PDF_iLSQR": "mv_pdf_ilsqr"}


# ------------------------------------------------------------------------------------------------ statistics

def anova2(Y):
    """Two-way ANOVA mean squares for an n x k matrix (rows = sessions, columns = methods)."""
    import numpy as np
    Y = np.asarray(Y, float)
    n, k = Y.shape
    gm = Y.mean()
    ssr = k * ((Y.mean(1) - gm) ** 2).sum()
    ssc = n * ((Y.mean(0) - gm) ** 2).sum()
    sst = ((Y - gm) ** 2).sum()
    sse = sst - ssr - ssc
    return n, k, ssr / (n - 1), ssc / (k - 1), sse / ((n - 1) * (k - 1))


def icc3_1(Y, alpha=0.05):
    from scipy import stats
    n, k, msr, msc, mse = anova2(Y)
    icc = (msr - mse) / (msr + (k - 1) * mse)
    F = msr / mse
    fl = F / stats.f.ppf(1 - alpha / 2, n - 1, (n - 1) * (k - 1))
    fu = F * stats.f.ppf(1 - alpha / 2, (n - 1) * (k - 1), n - 1)
    return {"icc": float(icc), "lo": float((fl - 1) / (fl + k - 1)), "hi": float((fu - 1) / (fu + k - 1)), "n": int(n), "k": int(k)}


def icc2_1(Y, alpha=0.05):
    from scipy import stats
    n, k, msr, msc, mse = anova2(Y)
    icc = (msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n)
    a = k * icc / (n * (1 - icc))
    b = 1 + k * icc * (n - 1) / (n * (1 - icc))
    v = (a * msc + b * mse) ** 2 / ((a * msc) ** 2 / (k - 1) + (b * mse) ** 2 / ((n - 1) * (k - 1)))
    fs, fss = stats.f.ppf(1 - alpha / 2, n - 1, v), stats.f.ppf(1 - alpha / 2, v, n - 1)
    lo = n * (msr - fs * mse) / (fs * (k * msc + (k * n - k - n) * mse) + n * msr)
    hi = n * (fss * msr - mse) / (k * msc + (k * n - k - n) * mse + n * fss * msr)
    return {"icc": float(icc), "lo": float(lo), "hi": float(hi), "n": int(n), "k": int(k)}


def pearson_ci(x, y, alpha=0.05):
    import numpy as np
    from scipy import stats
    x, y = np.asarray(x, float), np.asarray(y, float)
    n = len(x)
    r = float(np.corrcoef(x, y)[0, 1])
    z, q = math.atanh(max(min(r, 0.999999), -0.999999)), stats.norm.ppf(1 - alpha / 2)
    return {"r": r, "lo": math.tanh(z - q / math.sqrt(n - 3)), "hi": math.tanh(z + q / math.sqrt(n - 3)), "n": n}


def deming_slope(x, y):
    import numpy as np
    x, y = np.asarray(x, float), np.asarray(y, float)
    sxx, syy, sxy = np.var(x, ddof=1), np.var(y, ddof=1), np.cov(x, y)[0, 1]
    return float((syy - sxx + math.sqrt((syy - sxx) ** 2 + 4 * sxy ** 2)) / (2 * sxy)) if sxy != 0 else float("nan")


def sb(r):
    return 2 * r / (1 + r) if r > -1 else float("nan")


# ------------------------------------------------------------------------------------------------ tables

def key(r):
    return (r["code"], int(r["ses_order"]))


def fnum(v):
    try:
        x = float(v)
        return x if x == x else None
    except (TypeError, ValueError):
        return None


def load_roi(lx, eng, var):
    p = os.path.join(lx, "tables_coded", f"roi_{eng}_{var}.csv")
    return {key(r): r for r in L.read_csv(p)} if os.path.exists(p) else {}


def load_usable(lx, eng, pipe, rd):
    p = os.path.join(lx, "freeze", f"usable_{eng}_{pipe}_{rd}.csv")
    if not os.path.exists(p):
        sys.exit(f"{p} missing: run build_usable.py --stage usable first")
    return [r for r in L.read_csv(p)]


def first_and_pairs(ulist):
    use = sorted((r["code"], int(r["ses_order"])) for r in ulist if r["usable"] == "1")
    first, second = {}, {}
    for c, s in use:
        if c not in first:
            first[c] = s
        elif c not in second:
            second[c] = s
    return [(c, s) for c, s in first.items()], [((c, first[c]), (c, second[c])) for c in second]


def vals(tab, keys, col):
    out = []
    for k in keys:
        v = fnum(tab.get(k, {}).get(col))
        out.append(v)
    return out


# ------------------------------------------------------------------------------------------------ stage g2

def stage_g2(args):
    import numpy as np
    lx = L.lx_root()
    gate = [ln.split("\t") for ln in open(os.path.join(lx, "blind", "gate20.txt")).read().split("\n") if ln.strip()]
    gk = [(c, int(s)) for c, s in gate]
    out = {"utc": L.utcnow(), "A2_approved": args.a2 == "yes", "n_gate_sessions": len(gk)}
    if args.a2 == "yes":
        tabs = [load_roi(lx, "qsmxt", "g2"), load_roi(lx, "sepia", "g2")]
        labels = ["qsmxt", "sepia"]
    else:
        tabs = [load_roi(lx, "qsmxt", "g2"), load_roi(lx, "qsmxt", "g2_mv_vsharp_rts"), load_roi(lx, "qsmxt", "g2_mv_pdf_ilsqr")]
        labels = ["vsharp_medi", "vsharp_rts", "pdf_ilsqr"]
    for roi, name in (("chi_SNVTA", "G2"), ("chi_DN_R", "G2_DN")):
        M = [vals(t, gk, roi) for t in tabs]
        ok = [i for i in range(len(gk)) if all(m[i] is not None for m in M)]
        Y = np.array([[m[i] for m in M] for i in ok])
        if len(ok) < 10:
            out[name] = False
            out[f"{name}_detail"] = {"n": len(ok), "note": "fewer than 10 complete gate sessions: gate fails"}
            continue
        c3, c2 = icc3_1(Y), icc2_1(Y)
        det = {"n": len(ok), "methods": labels, "ICC3_1": c3, "ICC2_1": c2}
        if args.a2 == "yes":
            det["pearson"] = pearson_ci(Y[:, 0], Y[:, 1])
            det["deming_slope_sepia_on_qsmxt"] = deming_slope(Y[:, 0], Y[:, 1])
            det["mean_difference_ppb"] = float(np.mean(Y[:, 1] - Y[:, 0]))
            passed = c3["icc"] >= 0.90 and c3["lo"] >= 0.75
        else:
            passed = c3["icc"] >= 0.75
        out[name] = bool(passed)
        out[f"{name}_detail"] = det
    out["engine"] = "sepia" if (args.a2 == "yes" and not out["G2"]) else "qsmxt"
    if args.a2 != "yes" and not out["G2"]:
        out["caveat"] = "inversion-algorithm agreement (within-QSMxT G2) not established; engine stays QSMxT"
    L.write_json(args.out, out)
    print(json.dumps({k: out[k] for k in ("G2", "G2_DN", "engine")}, indent=1))


# ------------------------------------------------------------------------------------------------ stage final

def pc_checks(lx, eng, pipe):
    import numpy as np
    from scipy import stats
    tab = load_roi(lx, eng, BRANCH_VARIANT[pipe])
    ul = load_usable(lx, eng, pipe, "SNVTA")
    usable = [(r["code"], int(r["ses_order"])) for r in ul if r["usable"] == "1"]
    first, pairs = first_and_pairs(ul)
    rois = ["chi_GP", "chi_SNVTA", "chi_RN", "chi_PUT", "chi_CAU", "chi_THAL"]
    M = [[fnum(tab.get(k, {}).get(c)) for c in rois] for k in usable]
    M = np.array([m for m in M if all(v is not None for v in m)], float)
    pc1 = {"n": int(len(M))}
    if len(M):
        share = float(np.mean((M[:, 0] > 0) & (M[:, 0] > M[:, 5])))
        med = np.median(M, axis=0)
        pc1.update({"share_gp_pos_and_gt_thal": share, "gp_median_highest": bool(np.argmax(med) == 0)})
        pc1["pass"] = share >= 0.95 and pc1["gp_median_highest"]
    else:
        pc1["pass"] = False
    # PC2w: first two usable sessions of returners
    swc = {key(r): r.get("sw_class") for r in L.read_csv(os.path.join(lx, "blind", "session_flags.csv"))}
    d = []
    for (k0, k1) in pairs:
        a0, a1 = tab.get(k0, {}), tab.get(k1, {})
        v = [fnum(a0.get(c)) for c in ("chi_PUT", "chi_GP", "chi_THAL")] + [fnum(a1.get(c)) for c in ("chi_PUT", "chi_GP", "chi_THAL")]
        if all(x is not None for x in v):
            d.append((v[3] - v[0], v[4] - v[1], v[5] - v[2], swc.get(k0), swc.get(k1)))
    pc2 = {"n_pairs": len(d), "share_pairs_software_class_differs": (sum(x[3] != x[4] for x in d) / len(d)) if d else None}
    if d and pc2["share_pairs_software_class_differs"] > 0.20:
        same = [x for x in d if x[3] == x[4]]
        pc2["same_class_pairs"] = {"n": len(same), "mean_dbar_ppb": float(np.mean([(x[0] + x[1]) / 2 for x in same])) if same else None}
    if len(d) >= 5:
        dp = np.array([x[0] for x in d]); dg = np.array([x[1] for x in d]); dt = np.array([x[2] for x in d])
        dbar = (dp + dg) / 2
        t = float(dbar.mean() / (dbar.std(ddof=1) / math.sqrt(len(dbar))))
        p = float(stats.t.sf(t, len(dbar) - 1))
        state = "PASS" if (p < 0.05 and dp.mean() >= 0 and dg.mean() >= 0) else ("FAIL" if dbar.mean() <= 0 else "INDETERMINATE")
        pc2.update({"mean_dbar_ppb": float(dbar.mean()), "p_one_sided": p, "mean_dPUT": float(dp.mean()), "mean_dGP": float(dg.mean()),
                    "mean_dTHAL": float(dt.mean()), "mean_dbar_minus_dTHAL": float((dbar - dt).mean()), "state": state})
    else:
        pc2["state"] = "INDETERMINATE"
    # PC3: r(chi_GP, R2*_GP) over first-usable sessions
    xs = [(fnum(tab.get(k, {}).get("chi_GP")), fnum(tab.get(k, {}).get("r2s_GP"))) for k in first]
    xs = [x for x in xs if None not in x]
    pc3 = {"n": len(xs)}
    if len(xs) >= 10:
        r = float(np.corrcoef([x[0] for x in xs], [x[1] for x in xs])[0, 1])
        pc3.update({"r": r, "pass": r >= 0.30})
    else:
        pc3["pass"] = False
    return pc1, pc2, pc3


def stage_final(args):
    import numpy as np
    lx = L.lx_root()
    g2 = json.load(open(os.path.join(lx, "gates", "g2.json")))
    eng = g2["engine"]
    out = {"data_origin": "real", "utc": L.utcnow(), "engine": eng, "A2_approved": bool(g2["A2_approved"]),
           "G2": bool(g2["G2"]), "G2_DN": bool(g2["G2_DN"])}
    pc1, pc2, pc3 = pc_checks(lx, eng, "primary")
    out.update({"PC1": bool(pc1["pass"]), "PC2w": pc2["state"], "PC3": bool(pc3["pass"]),
                "PC_detail": {"PC1": pc1, "PC2w": pc2, "PC3": pc3}})
    pipe = "primary"
    if not out["PC1"] or out["PC2w"] == "FAIL" or not out["PC3"]:
        eng = "qsmxt"                         # the registered fallback is QSMxT PDF + iLSQR for either engine
        f1, f2, f3 = pc_checks(lx, eng, "fallback_PDF_iLSQR")
        out["fallback_PC_pass"] = bool(f1["pass"] and f2["state"] != "FAIL" and f3["pass"])
        out["PC_detail_fallback"] = {"PC1": f1, "PC2w": f2, "PC3": f3}
        pipe = "fallback_PDF_iLSQR"
    out["pipeline"] = pipe
    out["engine_analysed"] = eng
    # G3 on the common-support multiverse (SN-VTA), first-usable sessions of the pipeline's SN-VTA list
    ul = load_usable(lx, eng, pipe, "SNVTA")
    first, _ = first_and_pairs(ul)
    mvp = os.path.join(lx, "tables_coded", f"multiverse_{eng}.csv")
    def g3_for(colname):
        res = {"n": 0}
        if os.path.exists(mvp):
            mv = {}
            for r in L.read_csv(mvp):
                mv.setdefault(key(r), {})[(r["bf"], r["inv"])] = fnum(r[colname])
            nodes = [(b, i) for b in ("vsharp", "pdf", "lbv") for i in ("medi", "ilsqr", "rts")]
            Y = [[mv.get(k, {}).get(n) for n in nodes] for k in first]
            Y = np.array([y for y in Y if all(v is not None for v in y)], float)
            if len(Y) >= 10:
                res = {"n": int(len(Y)), "ICC3_1": icc3_1(Y), "ICC2_1": icc2_1(Y)}
        return res
    g3 = g3_for("chi_SNVTA_wb_mean")
    out["G3"] = bool(g3.get("ICC3_1", {}).get("icc", -1) >= 0.75)
    out["G3_detail"] = g3
    if not out["G3"]:                          # the SN-only fallback ROI must pass the same check (else a caveat)
        g3s = g3_for("chi_SNonly_wb_mean")
        out["G3_SNonly"] = bool(g3s.get("ICC3_1", {}).get("icc", -1) >= 0.75)
        out["G3_SNonly_detail"] = g3s
    rd = "SNVTA" if out["G3"] else "SN_only"
    col = "chi_SNVTA" if rd == "SNVTA" else "chi_SNonly"
    out["constants_roi"] = rd
    ul = load_usable(lx, eng, pipe, rd)
    first, pairs = first_and_pairs(ul)
    tab = load_roi(lx, eng, BRANCH_VARIANT[pipe])
    s13, s24 = load_roi(lx, eng, "split13"), load_roi(lx, eng, "split24")

    def split_rel(c, fail_col):
        xs = [(fnum(s13.get(k, {}).get(c)), fnum(s24.get(k, {}).get(c))) for k in first
              if s13.get(k, {}).get(fail_col, "0") != "1" and s24.get(k, {}).get(fail_col, "0") != "1"]
        xs = [x for x in xs if None not in x]
        if len(xs) < 10:
            return {"n": len(xs), "R_SB": None, "rL": None, "pass": False}
        pc = pearson_ci([x[0] for x in xs], [x[1] for x in xs])
        return {"n": len(xs), "r": pc["r"], "R_SB": sb(pc["r"]), "r_lower": pc["lo"], "rL": sb(pc["lo"]),
                "pass": sb(pc["r"]) >= 0.70, "label": "split-echo (noise) reliability, an upper bound on measurement reliability"}
    pc4 = split_rel(col, "fail_SNVTA" if rd == "SNVTA" else "fail_SNonly")
    pc4dn = split_rel("chi_DN_R", "fail_DN_R")
    # an undefined split-echo reliability is recorded as -1 (PC4 fails; outcome B is then not evaluable)
    nz = lambda v: -1.0 if v is None or not (v == v) else float(max(-1.0, min(1.0, v)))
    out.update({"PC4": bool(pc4["pass"]), "rL_primary": nz(pc4["rL"]), "r_point_primary": nz(pc4["R_SB"]), "PC4_detail": pc4,
                "PC4_DN": bool(pc4dn["pass"]), "rL_DN": pc4dn["rL"], "r_point_DN": pc4dn["R_SB"], "PC4_DN_detail": pc4dn})
    # pooled label-free constants over the first-usable sessions of the selected branch
    x = np.array([v for v in vals(tab, first, col) if v is not None])
    st = np.array([v for v in vals(tab, first, "chi_STR") if v is not None])
    out["pooled"] = {"mean_SNVTA": float(x.mean()), "sd_SNVTA": float(x.std(ddof=1)), "mean_STR": float(st.mean()),
                     "sd_STR": float(st.std(ddof=1)), "sd_total_SNVTA": float(x.std(ddof=1)), "n_first_usable": int(len(x)),
                     "roi": rd}
    # cross-session reliability (label- and age-free; returners' first two usable sessions)
    pr = [(fnum(tab.get(a, {}).get(col)), fnum(tab.get(b, {}).get(col))) for a, b in pairs]
    pr = [p for p in pr if None not in p]
    if len(pr) >= 10:
        pc = pearson_ci([p[0] for p in pr], [p[1] for p in pr])
        ic = icc3_1(np.array(pr))
        out.update({"r_cross": pc["r"], "rL_cross": sb(pc["lo"]) if pc["lo"] > -1 else None, "icc_cross_session": ic["icc"],
                    "cross_detail": {"n_pairs": len(pr), "pearson": pc, "ICC3_1": ic,
                                     "T1_informative": ic["icc"] >= 0.40,
                                     "note": "not age-adjusted; heterogeneous developmental change lowers it (a lower bound)"}})
    else:
        out.update({"r_cross": None, "rL_cross": None, "icc_cross_session": None})
    # SD of dref (valid CSF only), rho(SN-VTA, STR), r(DN_R, DGM)
    dref = [fnum(tab.get(k, {}).get(col)) - fnum(tab.get(k, {}).get(f"{col}_ref2")) for k in first
            if fnum(tab.get(k, {}).get(col)) is not None and fnum(tab.get(k, {}).get(f"{col}_ref2")) is not None]
    out["sd_dref_ppb"] = float(np.std(dref, ddof=1)) if len(dref) > 2 else None
    out["sd_dref_over_sd_primary"] = out["sd_dref_ppb"] / out["pooled"]["sd_SNVTA"] if out["sd_dref_ppb"] else None
    pairs2 = [(fnum(tab.get(k, {}).get(col)), fnum(tab.get(k, {}).get("chi_STR"))) for k in first]
    pairs2 = [p for p in pairs2 if None not in p]
    out["rho_SNVTA_STR"] = float(np.corrcoef([p[0] for p in pairs2], [p[1] for p in pairs2])[0, 1]) if len(pairs2) > 10 else None
    pairs3 = [(fnum(tab.get(k, {}).get("chi_DN_R")), fnum(tab.get(k, {}).get("chi_DGM"))) for k in first]
    pairs3 = [p for p in pairs3 if None not in p]
    out["r_DN_R_DGM"] = float(np.corrcoef([p[0] for p in pairs3], [p[1] for p in pairs3])[0, 1]) if len(pairs3) > 10 else None
    out["ref2_type_counts_first_usable"] = {t: sum(1 for k in first if tab.get(k, {}).get("ref2_type") == t) for t in ("CSF", "CSF_erode1", "none")}
    # G4, G5, PC5, n_min
    g4 = json.load(open(os.path.join(lx, "qc", "g4_final.json")))
    out["G4"] = bool(g4["G4_pass"])
    g5 = json.load(open(os.path.join(lx, "flags_coded", "g5_summary.json")))
    # amendment 1b: G5 is a per-session exclusion only (G5-H1 / G5-DN / IQ fail); the registered ">10% -> methods report"
    # chain trigger was removed after the menu revision (1a) was shown to worsen the T1w->MNI step. The gate field is
    # therefore TRUE whenever the per-session flags exist; the shares are reported, not used as a trigger.
    out["G5"] = True
    out["G5_detail"] = {"fail_share_all": g5["g5_fail_share"], "fail_share_nonIQ": g5.get("g5_fail_share_nonIQ"),
                        "n_iq_fail": g5.get("n_iq_fail"), "g5_pass": g5.get("g5_pass"), "g5_dn_pass": g5.get("g5_dn_pass"),
                        "vdc_threshold": g5.get("vdc_threshold"), "revised_once": True,
                        "rule": "per-session exclusion only (amendment 1b); chain trigger removed; registered rule would have read "
                                + ("pass" if g5["g5_fail_share"] <= 0.10 else "fail")}
    p5 = os.path.join(lx, "dwi", "pc5.json")
    out["PC5"] = bool(json.load(open(p5))["PC5"]) if os.path.exists(p5) else False
    nm = os.path.join(lx, "dwi", "n_min.json")
    out["n_min"] = json.load(open(nm))["n_min"] if os.path.exists(nm) else None
    geo = {key(r): fnum(r.get("msd_bs_vdc_mm")) for r in L.read_csv(os.path.join(lx, "tables_coded", "session_geom.csv"))}
    msd = [geo.get(k) for k in first if geo.get(k) is not None]
    out["boundary_distance_mm"] = {"n": len(msd), "median": float(np.median(msd)) if msd else None,
                                   "p95": float(np.percentile(msd, 95)) if msd else None,
                                   "note": "brainstem + ventral-DC symmetric mean surface distance, echo-1 SynthSeg vs T1w SynthSeg"}
    out["counts"] = {"n_first_usable": len(first), "n_returner_pairs": len(pairs)}
    out["D10_dR2s"] = None
    L.write_json(args.out, out)
    print(json.dumps({k: out[k] for k in ("engine", "pipeline", "constants_roi", "PC1", "PC2w", "PC3", "PC4", "rL_primary",
                                          "G3", "G4", "G5", "PC4_DN", "G2_DN", "PC5")}, indent=1))


def main():
    a = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    a.add_argument("--stage", choices=["g2", "final"], required=True)
    a.add_argument("--a2", choices=["yes", "no"], default="no", help="was approval A2 (SEPIA + Cornell MEDI) given?")
    a.add_argument("--out", default=None)
    args = a.parse_args()
    lx = L.lx_root()
    args.out = args.out or os.path.join(lx, "gates", "g2.json" if args.stage == "g2" else "gates_preunlock.json")
    L.assert_inside(args.out)
    if os.path.exists(args.out):
        sys.exit(f"{args.out} exists (computed once)")
    (stage_g2 if args.stage == "g2" else stage_final)(args)


if __name__ == "__main__":
    main()
