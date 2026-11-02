#!/usr/bin/env python3
"""analyze_kit2.py - LEAD-chi human visual QC, kit 2: complete-human-QC refit of M1 and test-retest (POST HOC).

Implements ../PLAN_kit2.md (fixed before any kit-2 rating). Nothing here changes a frozen file or re-decides H1.

Usage
  python3 analyze_kit2.py --kit2-ratings <returned kit-2 CSV> --out <new folder>
         [--kit1-ratings ../../returned/LEADchi_humanQC_rater_20261028T1724Z.csv] [--simulated] [--rscript Rscript]

Per session (plan section 2)
  * 101 kit-1 sessions: kit-1 rating, worse grade of the kit-1 original and its kit-1 repeat, placement fail if either
    fails (identical to kit-1 analyze_human_qc.py part_b, i.e. the n = 89 run);
  * 79 new sessions: the kit-2 rating; the 21 kit-2 repeats are used only for test-retest;
  * a session without a valid human grade (placement) keeps the AI grade (placement), as in H_replace; counted.
  * usable = no frozen automated fail_reason component, human grade <= 1, human placement pass, and (if it becomes
    usable) the frozen SN-VTA re-check; first usable = earliest usable ses_order per youth. The first-usable logic is
    kit 1's reconstruct(), imported unchanged, after the same reproduction check of the frozen table.
Model (plan section 3)
  * The unchanged kit-1 script human_qc_kit/analysis/refit_m1_human_qc.R sources the frozen R helpers read-only and
    fits M1 for scenarios AI_frozen (must give n = 80, beta = 4.96280154409223), K1_H_replace (kit-1 grades only; with
    the real kit-1 CSV it must give n = 89, beta = 6.29013622956428) and H_all180 (primary of the plan).
Test-retest (plan section 4): kit-2 repeat vs the kit-1 primary item of the same session; exact agreement (Wilson CI),
  linear-weighted kappa, 3x3 table, placement agreement; also without the example-montage session if it is drawn.
Secondary (plan section 3): census agreement human vs AI over the 180 sessions (primary without the example-montage
  session q15499e9a, also with it) and the human fail rate among the 119 AI-pass sessions (Wilson CI).

Kit-1 functions (read_ratings, wkappa, wilson, reconstruct, clean) are imported from
human_qc_kit/analysis/analyze_human_qc.py without modification (no bytecode is written there).
"""
from __future__ import annotations

import argparse
import datetime
import hashlib
import json
import os
import subprocess
import sys

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

K1 = "/path/to/Penn_LEAD/LEADchi/human_qc_kit"
K2 = os.path.join(K1, "kit2_2026-10-28")
sys.path.insert(0, os.path.join(K1, "analysis"))
import analyze_human_qc as A1  # noqa: E402

REFIT_R = os.path.join(K1, "analysis", "refit_m1_human_qc.R")
REAL_K1_CSV = os.path.join(K1, "returned", "LEADchi_humanQC_rater_20261028T1724Z.csv")
REAL_K1_SHA = "1975118e501066fd87cb7e18a096f537a034caf3b38897590a927649cc2118c3"
K1_EXPECT = {"n": 89, "beta": 6.29013622956428}           # kit-1 H_replace, real run 2026-10-28
EXAMPLE_QC = "q15499e9a"                                    # kit-1 item hfaec73, example montage seen by the rater
SCEN = ["AI_frozen", "K1_H_replace", "H_all180"]


def sha256(p):
    return A1.sha256(p)


def session_ratings(key, r):
    """Per-session human rating from a key + ratings: worse grade, placement fail if any (kit-1 part_b rule)."""
    al = key.merge(r, on="kit_id", how="left")
    hs = {}
    for (c, s), g in al.groupby(["code", "ses_order"]):
        gg, pp = g[g.h_grade.notna()], g[g.h_place.notna()]
        hs[(c, int(s))] = {"h_grade": float(gg.h_grade.max()) if len(gg) else float("nan"),
                           "h_place": ("fail" if (pp.h_place == "fail").any() else "pass") if len(pp) else None,
                           "kit_ids": ";".join(g.kit_id)}
    return hs, al


def kappa_block(h, a):
    h, a = np.asarray(h, int), np.asarray(a, int)
    T = np.zeros((3, 3), int)
    np.add.at(T, (h, a), 1)
    k = A1.wkappa(h, a, 3, "linear")[0]
    return {"n": int(len(h)), "exact_agreement": float((h == a).mean()) if len(h) else float("nan"),
            "exact_wilson95": A1.wilson(int((h == a).sum()), len(h)),
            "kappa_linear": k, "kappa_quadratic": A1.wkappa(h, a, 3, "quadratic")[0],
            "table_rows_first_cols_second": T.tolist()}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--kit2-ratings", dest="kit2", required=True)
    ap.add_argument("--kit1-ratings", dest="kit1", default=REAL_K1_CSV)
    ap.add_argument("--out", required=True)
    ap.add_argument("--simulated", action="store_true", help="label every output as a SIMULATED test run")
    ap.add_argument("--key1", default=os.path.join(K1, "_key_DO_NOT_SEND", "kit_key.csv"))
    ap.add_argument("--key2", default=os.path.join(K2, "_key_DO_NOT_SEND", "kit2_key.csv"))
    for k in ("merged", "grades", "gates", "rdir"):
        ap.add_argument(f"--{k}", default=A1.DEF[k])
    ap.add_argument("--rscript", default="Rscript")
    args = ap.parse_args()
    label = "SIMULATED_" if args.simulated else ""

    key1 = pd.read_csv(args.key1, dtype={"code": str, "is_duplicate_of": str}, keep_default_na=False)
    key2 = pd.read_csv(args.key2, dtype={"code": str, "kit1_id": str}, keep_default_na=False).rename(
        columns={"kit2_id": "kit_id"})
    for k in (key1, key2):
        k["ai_grade"] = k.ai_grade.astype(int)
        k["ses_order"] = k.ses_order.astype(int)
    real_k1 = sha256(args.kit1) == REAL_K1_SHA         # checked before any rating file is parsed
    if args.simulated and (real_k1 or os.path.realpath(args.kit1) == os.path.realpath(REAL_K1_CSV)):
        sys.exit("--simulated run given the real kit-1 ratings: test runs use simulated ratings only; no rating file "
                 "was parsed and nothing was written")
    # ---- read and validate both rating files before anything is written
    r1, iss1 = A1.read_ratings(args.kit1, key1)
    r2, iss2 = A1.read_ratings(args.kit2, key2)
    n2 = int(r2.h_grade.notna().sum())
    if n2 == 0:
        sys.exit(f"no kit-2 item has a valid grade in {args.kit2} (is it the kit-2 CSV?); nothing was written")
    if len(iss2["unknown_kit_ids"]) > len(r2):
        sys.exit(f"{args.kit2} holds more unknown kit IDs ({len(iss2['unknown_kit_ids'])}) than kit-2 IDs: wrong file? "
                 "nothing was written")
    if not r1.h_grade.notna().any():
        sys.exit(f"no kit-1 item has a valid grade in {args.kit1}; nothing was written")
    os.makedirs(args.out, exist_ok=True)
    if any(f.endswith("kit2_results.json") for f in os.listdir(args.out)):
        sys.exit(f"{args.out} already holds results: use a new output folder (nothing is overwritten)")

    # ---- merged table and frozen reconstruction check (identical to kit-1 part_b)
    mg = pd.read_csv(args.merged, dtype={"code": str}, keep_default_na=True)
    gf = pd.read_csv(args.grades, dtype={"code": str})
    mg = mg.merge(gf[["code", "ses_order", "qc_grade", "placement"]].rename(
        columns={"qc_grade": "ai_grade_gf", "placement": "ai_place"}), on=["code", "ses_order"], how="left")
    if not ((mg.qc_grade == mg.ai_grade_gf) | (mg.qc_grade.isna() & mg.ai_grade_gf.isna())).all():
        sys.exit("merged-table qc_grade differs from grades_final.csv: stop")
    mg["reason_set"] = [set(x.split(";")) if isinstance(x, str) and x else set() for x in mg.fail_reason]
    if set().union(*mg.reason_set) - set(A1.REASONS):
        sys.exit("unknown fail_reason components: stop")
    mg["_posthoc_plaus_fail"] = ((mg.chi_SNVTA < -50) | (mg.chi_SNVTA > 300) | (mg.fimp_SNVTA > 0.05) |
                                 (mg.inmask_SNVTA < 0.95) | mg.chi_SNVTA.isna() | mg.roi_vol.isna() |
                                 mg.Q_phase.isna()).fillna(False)
    ai_g = mg.ai_grade_gf.to_numpy(float)
    ai_p = mg.ai_place.fillna("missing").to_numpy(object)
    u0, f0, r0 = A1.reconstruct(mg, ai_g, ai_p)
    repro = {"usable_equal": bool((u0 == mg.usable.to_numpy()).all()),
             "first_usable_equal": bool((f0 == mg.first_usable.to_numpy()).all()),
             "fail_reason_equal": bool(all(set(a.split(";")) - {""} == b for a, b in zip(r0, mg.reason_set)))}
    if not all(repro.values()):
        sys.exit(f"the first-usable reconstruction does not reproduce the frozen table: {repro}")

    # ---- per-session human ratings
    hs1, al1 = session_ratings(key1, r1)                                # kit 1: original + kit-1 repeat
    new2 = key2[key2.stratum == "AI_pass_new"]
    hs2, al2 = session_ratings(new2, r2)                                # kit 2: the 79 new sessions only
    if set(hs1) & set(hs2):
        sys.exit("kit-1 and kit-2 new sessions overlap: stop")
    keys = list(zip(mg.code, mg.ses_order.astype(int)))
    get = lambda hs, f: np.array([hs[k][f] if k in hs else (np.nan if f == "h_grade" else None) for k in keys], object)
    g1, p1 = get(hs1, "h_grade").astype(float), get(hs1, "h_place")
    hsa = {**hs1, **hs2}
    ga, pa = get(hsa, "h_grade").astype(float), get(hsa, "h_place")
    graded = ~np.isnan(ai_g)
    cov = {"AI_graded_sessions": int(graded.sum()),
           "sessions_with_human_grade": int((graded & ~np.isnan(ga)).sum()),
           "sessions_with_human_placement": int((graded & np.array([x is not None for x in pa])).sum()),
           "kit1_sessions": len(hs1), "kit2_new_sessions": len(hs2),
           "kit2_new_without_grade": int(sum(np.isnan(v["h_grade"]) for v in hs2.values())),
           "kit2_new_without_placement": int(sum(v["h_place"] is None for v in hs2.values())),
           "kit1_sessions_without_grade": int(sum(np.isnan(v["h_grade"]) for v in hs1.values()))}
    cov["complete_human_qc"] = bool(cov["sessions_with_human_grade"] == cov["AI_graded_sessions"] ==
                                    cov["sessions_with_human_placement"])
    sel = lambda h, a: np.where(np.array([x is not None and x == x for x in h], bool) if h.dtype == object else ~np.isnan(h), h, a)
    sc = {"AI_frozen": (ai_g, ai_p),
          "K1_H_replace": (sel(g1, ai_g).astype(float), sel(p1, ai_p)),
          "H_all180": (sel(ga, ai_g).astype(float), sel(pa, ai_p))}
    tab = pd.DataFrame({"sess_key": mg.participant_id + "_" + mg.session_id, "code": mg.code, "ses_order": mg.ses_order,
                        "ai_grade": ai_g, "ai_placement": ai_p, "human_grade_kit1": g1, "human_placement_kit1": p1,
                        "human_grade_all": ga, "human_placement_all": pa,
                        "source": ["kit1" if k in hs1 else "kit2" if k in hs2 else "" for k in keys]})
    summ, firsts = {}, {}
    for name, (g, p) in sc.items():
        u, f, rr = A1.reconstruct(mg, g, p)
        tab[f"usable__{name}"], tab[f"first_usable__{name}"], tab[f"fail_reason__{name}"] = u, f, rr
        firsts[name] = set(tab.sess_key[f == 1])
        summ[name] = {"n_usable_sessions": int(u.sum()), "n_first_usable": int(f.sum()),
                      "posthoc_plausibility_recheck_excluded": tab.sess_key[[
                          "posthoc_snvta_plausibility_recheck" in x for x in rr]].tolist()}
    subj = lambda s: {k.rsplit("_ses-", 1)[0] for k in s}
    for ref in ("AI_frozen", "K1_H_replace"):
        a0, a1 = firsts[ref], firsts["H_all180"]
        summ["H_all180"][f"vs_{ref}"] = {"youths_added": sorted(subj(a1) - subj(a0)),
                                         "youths_removed": sorted(subj(a0) - subj(a1)),
                                         "youths_with_changed_first_session": sorted(subj(a0 ^ a1) & subj(a0) & subj(a1))}
    scen_csv = os.path.join(args.out, f"{label}kit2_scenario_usable.csv")
    tab.to_csv(scen_csv, index=False)

    # ---- M1 refit with the unchanged kit-1 R script (frozen helpers, read-only)
    rj = os.path.join(args.out, f"{label}kit2_m1_refit.json")
    cmd = [args.rscript, REFIT_R, "--merged", args.merged, "--scen", scen_csv, "--gates", args.gates,
           "--rdir", args.rdir, "--out", rj, "--bound", repr(A1.FROZEN_BOUND)]
    pr = subprocess.run(cmd, capture_output=True, text=True)
    fit = json.load(open(rj))["scenarios"] if pr.returncode == 0 else None
    m1 = {}
    for s in SCEN:
        x = (fit or {}).get(s, {})
        if not x or "error" in x:
            m1[s] = {"error": x.get("error", "R refit failed")}
            continue
        m1[s] = {"n": x["n"], "beta_PS_ppb": x["beta_PS_ppb"], "se_hc3": x["se_hc3"], "ci95_ppb": x["ci95_ppb"],
                 "p_two": x["p_two"], "pB_fixed_bound": x["bound_fixed"]["pB_one_sided"],
                 "bound_fixed_ppb": x["bound_fixed"]["bound_ppb"], "pB_rederived_bound": x["bound_rederived"]["pB"],
                 "A_beta_negative_p05": bool(x["A_beta_negative_p05"]), "label_descriptive": x["outcome_label_descriptive"],
                 "cells": x["cells"]}
    checks = {"frozen_reconstruction": repro,
              "AI_frozen_reproduces_confirmatory": bool(m1["AI_frozen"].get("n") == 80 and
                                                        abs(m1["AI_frozen"].get("beta_PS_ppb", np.nan) - A1.FROZEN_BETA) < 1e-8),
              "kit1_ratings_are_real_returned_csv": real_k1}
    if real_k1:
        checks["K1_H_replace_reproduces_n89_run"] = bool(m1["K1_H_replace"].get("n") == K1_EXPECT["n"] and abs(
            m1["K1_H_replace"].get("beta_PS_ppb", np.nan) - K1_EXPECT["beta"]) < 1e-8)
    else:
        checks["K1_H_replace_reproduces_n89_run"] = "not applicable (kit-1 ratings are not the real returned CSV)"

    # ---- test-retest on the 21 repeats (kit-2 repeat vs the kit-1 primary item of the same session)
    rep = key2[key2.is_kit1_repeat.astype(int) == 1].merge(r2, on="kit_id", how="left")
    k1p = al1[["kit_id", "h_grade", "h_place"]].rename(columns={"kit_id": "kit1_id", "h_grade": "k1_grade",
                                                                "h_place": "k1_place"})
    rep = rep.merge(k1p, on="kit1_id", how="left")
    both = rep[rep.h_grade.notna() & rep.k1_grade.notna()]
    tr = {"n_repeats": int(len(rep)), "n_both_graded": int(len(both)),
          "all": kappa_block(both.h_grade.astype(int), both.k1_grade.astype(int)),
          "table_note": "rows = kit-2 grade, columns = kit-1 grade",
          "by_stratum_counts": {s: {"n": int(len(g)), "identical": int((g.h_grade == g.k1_grade).sum())}
                                for s, g in both.groupby("stratum")}}
    bp = rep[rep.h_place.notna() & rep.k1_place.notna()]
    tr["placement_agreement"] = {"n": int(len(bp)), "identical": int((bp.h_place == bp.k1_place).sum())}
    tr["example_montage_session_among_repeats"] = bool((rep.qc_id == EXAMPLE_QC).any())
    if tr["example_montage_session_among_repeats"]:
        b2 = both[both.qc_id != EXAMPLE_QC]
        tr["without_example_montage_session"] = kappa_block(b2.h_grade.astype(int), b2.k1_grade.astype(int))
    rep.drop(columns=["note"], errors="ignore").to_csv(os.path.join(args.out, f"{label}kit2_test_retest.csv"), index=False)

    # ---- secondary: census agreement (human session grade vs AI) and human fail rate among AI-pass sessions
    qid = dict(zip(zip(pd.concat([key1.code, key2.code]), pd.concat([key1.ses_order, key2.ses_order])),
                   pd.concat([key1.qc_id, key2.qc_id])))
    cen = pd.DataFrame({"qc_id": [qid.get(k, "") for k in keys], "h": ga, "a": ai_g})
    cen = cen[cen.h.notna() & cen.a.notna()]
    secondary = {"census_with_example_session": kappa_block(cen.h.astype(int), cen.a.astype(int)),
                 "census_primary_without_example_session": kappa_block(cen[cen.qc_id != EXAMPLE_QC].h.astype(int),
                                                                       cen[cen.qc_id != EXAMPLE_QC].a.astype(int)),
                 "table_note": "rows = human grade, columns = AI grade"}
    apass = cen[cen.a < 2]
    m_fail = int((apass.h == 2).sum())
    secondary["human_fail_rate_among_AI_pass"] = {"n": int(len(apass)), "human_grade_2": m_fail,
                                                  "rate": m_fail / len(apass) if len(apass) else float("nan"),
                                                  "wilson95": A1.wilson(m_fail, len(apass))}

    # ---- interpretation (plan section 5)
    h = m1["H_all180"]
    if "error" in h:
        interp = "H_all180 refit failed: no interpretation"
    elif not checks["AI_frozen_reproduces_confirmatory"] or checks["K1_H_replace_reproduces_n89_run"] is False:
        interp = "reproduction check failed: do not interpret"
    elif h["pB_fixed_bound"] < 0.05 and not h["A_beta_negative_p05"]:
        interp = ("the result is robust to replacing the AI visual-QC grade with a first-author human grade for all 180 "
                  "sessions" + ("" if cov["complete_human_qc"] else " (NOTE: not every session has a human grade; see coverage)"))
    else:
        interp = "the result is sensitive to the visual-QC rater"
    nominal = "error" not in h and h["p_two"] < 0.05 and h["beta_PS_ppb"] > 0
    meta = {"simulated": args.simulated, "utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "kit1_ratings": os.path.abspath(args.kit1), "kit1_ratings_sha256": sha256(args.kit1),
            "kit2_ratings": os.path.abspath(args.kit2), "kit2_ratings_sha256": sha256(args.kit2),
            "key1_sha256": sha256(args.key1), "key2_sha256": sha256(args.key2),
            "plan_sha256": sha256(os.path.join(K2, "PLAN_kit2.md")), "script_sha256": sha256(os.path.abspath(__file__)),
            "kit1_module_sha256": sha256(A1.__file__), "refit_R_sha256": sha256(REFIT_R),
            "rscript_returncode": pr.returncode, "issues_kit1": {k: (len(v) if isinstance(v, (list, dict)) else v) for k, v in iss1.items()},
            "issues_kit2": {k: (len(v) if isinstance(v, (list, dict)) else v) for k, v in iss2.items()}}
    if pr.returncode != 0:
        meta["rscript_stderr"] = pr.stderr[-4000:]
    res = {"SIMULATED": args.simulated, "label": "POST HOC: complete human visual QC (kit 1 + kit 2); not confirmatory",
           "meta": meta, "coverage": cov, "checks": checks, "m1": m1, "sample_changes": summ, "test_retest": tr,
           "secondary": secondary, "interpretation": interp,
           "nominally_significant_positive_beta": bool(nominal),
           "nominal_note": "a nominally significant positive beta is NOT interpreted as a finding (plan section 5)"}
    with open(os.path.join(args.out, f"{label}kit2_results.json"), "w") as fo:
        json.dump(A1.clean(res), fo, indent=1, allow_nan=False)
    report(res, os.path.join(args.out, f"{label}kit2_report.md"))
    print(json.dumps({"simulated": args.simulated, "complete": cov["complete_human_qc"], "checks": checks,
                      "H_all180": {k: h.get(k) for k in ("n", "beta_PS_ppb", "ci95_ppb", "p_two", "pB_fixed_bound")},
                      "retest_exact": tr["all"]["exact_agreement"], "retest_kappa": tr["all"]["kappa_linear"],
                      "interpretation": interp}, indent=1, default=str))


def f(x, d=3):
    if x is None or (isinstance(x, float) and x != x):
        return "NA"
    return f"{x:.{d}f}" if isinstance(x, float) else str(x)


def report(res, path):
    sim = res["SIMULATED"]
    L = []
    if sim:
        L += ["# SIMULATED TEST RUN: NOT REAL RATINGS. DO NOT REPORT.", ""]
    L += ["# LEAD-χ human visual QC, kit 2: complete-human-QC result (POST HOC)", "",
          "Plan: `PLAN_kit2.md` (SHA-256 `" + res["meta"]["plan_sha256"] + "`). Rater: first author (not independent).", "",
          "## Coverage", ""]
    L += [f"- {k}: {v}" for k, v in res["coverage"].items()]
    L += ["", "## Checks", ""] + [f"- {k}: {v}" for k, v in res["checks"].items()]
    L += ["", "## M1 (frozen helpers; bound test at the frozen bound −6.1318 ppb)", "",
          "| Scenario | n | β_PS (ppb) | 95 % CI | p (two-sided) | p_B | label |", "|---|---|---|---|---|---|---|"]
    for s, x in res["m1"].items():
        if "error" in x:
            L.append(f"| {s} | error: {x['error']} | | | | | |")
        else:
            L.append(f"| {s} | {x['n']} | {f(x['beta_PS_ppb'])} | [{f(x['ci95_ppb'][0])}, {f(x['ci95_ppb'][1])}] | "
                     f"{f(x['p_two'], 4)} | {x['pB_fixed_bound']:.2e} | {x['label_descriptive']} |")
    L += ["", "`AI_frozen` = registered primary (n = 80); `K1_H_replace` = kit-1 grades only (the n = 89 run when the real kit-1 CSV is used); **`H_all180` = this plan's "
          "primary**. Reported next to the registered result; the registered outcome stays B.", "",
          f"**Interpretation (plan section 5):** {res['interpretation']}.", ""]
    if res["nominally_significant_positive_beta"]:
        L += ["Nominally significant positive β_PS: **not interpreted as a finding** (plan section 5).", ""]
    ch = res["sample_changes"]["H_all180"]
    L += ["### Sample changes of H_all180", ""]
    for ref in ("vs_AI_frozen", "vs_K1_H_replace"):
        x = ch[ref]
        L.append(f"- {ref}: added {len(x['youths_added'])}, removed {len(x['youths_removed'])}, changed first session "
                 f"{len(x['youths_with_changed_first_session'])}")
    L.append(f"- SN–VTA re-check excluded: {len(ch['posthoc_plausibility_recheck_excluded'])}")
    t = res["test_retest"]
    L += ["", "## Test-retest (21 hidden repeats vs kit-1 grade)", "",
          f"- pairs with both grades: {t['n_both_graded']} / {t['n_repeats']}",
          f"- exact agreement {f(t['all']['exact_agreement'])} (Wilson 95 % CI {f(t['all']['exact_wilson95'][0])}–"
          f"{f(t['all']['exact_wilson95'][1])}); linear-weighted κ {f(t['all']['kappa_linear'])}",
          f"- table (rows kit 2, columns kit 1, grades 0/1/2): {t['all']['table_rows_first_cols_second']}",
          f"- placement identical: {t['placement_agreement']['identical']} / {t['placement_agreement']['n']}",
          f"- by item type: {t['by_stratum_counts']}"]
    if "without_example_montage_session" in t:
        w = t["without_example_montage_session"]
        L.append(f"- without the example-montage session: exact {f(w['exact_agreement'])}, κ {f(w['kappa_linear'])}")
    L.append("- Descriptive only; recall is possible, so agreement may be optimistic.")
    s = res["secondary"]
    c0, c1 = s["census_primary_without_example_session"], s["census_with_example_session"]
    fr = s["human_fail_rate_among_AI_pass"]
    L += ["", "## Secondary (descriptive): census agreement human vs AI", "",
          f"- without the example-montage session (primary): n {c0['n']}, linear κ {f(c0['kappa_linear'])}, exact "
          f"{f(c0['exact_agreement'])}; table (rows human, cols AI) {c0['table_rows_first_cols_second']}",
          f"- with it: n {c1['n']}, linear κ {f(c1['kappa_linear'])}, exact {f(c1['exact_agreement'])}",
          f"- human grade 2 among AI-pass sessions: {fr['human_grade_2']} / {fr['n']} = {f(fr['rate'])} "
          f"(Wilson 95 % CI {f(fr['wilson95'][0])}–{f(fr['wilson95'][1])})", ""]
    open(path, "w", encoding="utf-8").write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
