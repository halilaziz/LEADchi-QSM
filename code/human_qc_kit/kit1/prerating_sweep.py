#!/usr/bin/env python3
"""prerating_sweep.py - POST HOC, run BEFORE any human rating exists (2026-10-07; review item M3).

How much can the human grade of the kit sessions move the H1 sensitivity label at all? This script reproduces the
checker's sweep (REVIEW_by_checker.md, M3; same logic and seed as _checker_tmp/sweep/sweep.py) and adds the two
extremes. It uses NO human rating: every scenario is synthetic.

  * 300 random human gradings of the kit sessions: P(human passes an AI fail) ~ U(0, 1) and P(human fails an AI pass)
    ~ U(0, 0.4), drawn per scenario (numpy default_rng(1)); placement kept from the AI; sessions outside the kit keep
    the AI grade.
  * extremes: every kit session passed (grade 1 for AI grade 2) and every kit session failed (grade 2).
The first-usable rule is re-applied with analyze_human_qc.reconstruct() and M1 is refitted with refit_m1_human_qc.R
(frozen R helpers), exactly as in the analysis.

Usage: python3 prerating_sweep.py [--out prerating_sweep_2026-10-07] [--n 300]   (refuses to overwrite)
"""
import argparse
import json
import os
import subprocess
import sys

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import analyze_human_qc as A  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=os.path.join(HERE, "prerating_sweep_2026-10-07"))
    ap.add_argument("--n", type=int, default=300)
    a = ap.parse_args()
    if os.path.exists(os.path.join(a.out, "sweep_results.csv")):
        sys.exit(f"{a.out} already holds results (nothing is overwritten)")
    os.makedirs(a.out, exist_ok=True)
    mg = pd.read_csv(A.DEF["merged"], dtype={"code": str})
    gf = pd.read_csv(A.DEF["grades"], dtype={"code": str})
    mg = mg.merge(gf[["code", "ses_order", "qc_grade", "placement"]].rename(
        columns={"qc_grade": "ai_grade_gf", "placement": "ai_place"}), on=["code", "ses_order"], how="left")
    mg["reason_set"] = [set(x.split(";")) if isinstance(x, str) and x else set() for x in mg.fail_reason]
    mg["_posthoc_plaus_fail"] = ((mg.chi_SNVTA < -50) | (mg.chi_SNVTA > 300) | (mg.fimp_SNVTA > 0.05) |
                                 (mg.inmask_SNVTA < 0.95) | mg.chi_SNVTA.isna() | mg.roi_vol.isna() |
                                 mg.Q_phase.isna()).fillna(False)
    key = pd.read_csv(A.DEF["key"], dtype={"code": str}, keep_default_na=False)
    kit = set(zip(key.code, key.ses_order.astype(int)))
    inkit = np.array([(c, int(s)) in kit for c, s in zip(mg.code, mg.ses_order)])
    ai_g = mg.ai_grade_gf.to_numpy(float)
    ai_p = mg.ai_place.fillna("missing").to_numpy(object)
    u0, f0, _ = A.reconstruct(mg, ai_g, ai_p)
    assert (u0 == mg.usable.to_numpy()).all() and (f0 == mg.first_usable.to_numpy()).all(), "reconstruction check"
    rng = np.random.default_rng(1)
    cols, meta = {"sess_key": (mg.participant_id + "_" + mg.session_id).to_numpy()}, []
    for s in range(a.n):
        pf, pm = rng.uniform(0, 1), rng.uniform(0, 0.4)
        g = ai_g.copy()
        for i in np.where(inkit)[0]:
            if ai_g[i] == 2 and rng.random() < pf:
                g[i] = 1
            elif ai_g[i] < 2 and rng.random() < pm:
                g[i] = 2
        u, f, _ = A.reconstruct(mg, g, ai_p)
        cols[f"usable__s{s}"], cols[f"first_usable__s{s}"] = u, f
        meta.append({"s": f"s{s}", "pf": pf, "pm": pm})
    for name, g in (("all_kit_pass", np.where(inkit & (ai_g == 2), 1.0, ai_g)),
                    ("all_kit_fail", np.where(inkit, 2.0, ai_g))):
        u, f, _ = A.reconstruct(mg, g, ai_p)
        cols[f"usable__{name}"], cols[f"first_usable__{name}"] = u, f
        meta.append({"s": name, "pf": float("nan"), "pm": float("nan")})
    scen = os.path.join(a.out, "scenarios.csv")
    pd.DataFrame(cols).to_csv(scen, index=False)
    fj = os.path.join(a.out, "m1_refit_sweep.json")
    cmd = ["Rscript", os.path.join(HERE, "refit_m1_human_qc.R"), "--merged", A.DEF["merged"], "--scen", scen,
           "--gates", A.DEF["gates"], "--rdir", A.DEF["rdir"], "--out", fj, "--bound", repr(A.FROZEN_BOUND)]
    pr = subprocess.run(cmd, capture_output=True, text=True)
    if pr.returncode != 0:
        sys.exit(pr.stderr[-2000:])
    fit = json.load(open(fj))["scenarios"]
    rows = []
    for m in meta:
        x = fit[m["s"]]
        if "error" in x:
            rows.append(dict(m, error=x["error"]))
            continue
        rows.append(dict(m, n=x["n"], beta=x["beta_PS_ppb"], se=x["se_hc3"], pB_fixed=x["bound_fixed"]["pB_one_sided"],
                         pB_rederived=x["bound_rederived"]["pB"], A=x["A_beta_negative_p05"],
                         label=x["outcome_label_descriptive"]))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(a.out, "sweep_results.csv"), index=False)
    rnd = df[df.s.str.match(r"s\d+$")]
    summ = {"label": "POST HOC pre-rating sweep (no human rating used); review item M3",
            "n_random_scenarios": int(len(rnd)), "label_counts_random": rnd.label.value_counts().to_dict(),
            "A_any_random": bool(rnd.A.any()), "n_range_random": [int(rnd.n.min()), int(rnd.n.max())],
            "beta_range_random": [float(rnd.beta.min()), float(rnd.beta.max())],
            "pB_fixed_max_random": float(rnd.pB_fixed.max()),
            "extremes": df[~df.s.str.match(r"s\d+$")].drop(columns=["pf", "pm"]).to_dict(orient="records")}
    json.dump(summ, open(os.path.join(a.out, "sweep_summary.json"), "w"), indent=1)
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
