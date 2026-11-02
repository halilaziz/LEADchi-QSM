#!/usr/bin/env python3
"""LEAD-chi post hoc validation: independent statsmodels re-fit of selected R models (POST HOC check).

Reads the exact design matrices exported by validation_local.R (verify/design_<tag>.csv: column y, then the model
matrix including the intercept) and the R estimates (verify/R_reference_estimates.csv); re-fits OLS with HC3 in
statsmodels and writes verify/statsmodels_check.json. Checks the estimation and HC3 algebra, not the design coding.

    python3 verify_statsmodels.py
"""
import json
import os

import numpy as np
import pandas as pd
import statsmodels
import statsmodels.api as sm

OUT = "/path/to/Penn_LEAD/LEADchi/results_posthoc/validation_local"
VER = os.path.join(OUT, "verify")

ref = pd.read_csv(os.path.join(VER, "R_reference_estimates.csv"))
rows = []
for _, r in ref.iterrows():
    d = pd.read_csv(os.path.join(VER, f"design_{r['tag']}.csv"))
    y = d.pop("y").to_numpy()
    fit = sm.OLS(y, d.to_numpy()).fit(cov_type="HC3")
    j = list(d.columns).index(r["term"])
    est, se = float(fit.params[j]), float(fit.bse[j])
    rows.append({"tag": r["tag"], "term": r["term"], "R_est": float(r["est"]), "sm_est": est,
                 "R_se_hc3": float(r["se_hc3"]), "sm_se_hc3": se, "n": int(fit.nobs), "df_resid": float(fit.df_resid),
                 "abs_diff_est": abs(est - r["est"]), "abs_diff_se": abs(se - r["se_hc3"]),
                 "pass": bool(abs(est - r["est"]) < 1e-8 and abs(se - r["se_hc3"]) < 1e-8)})
res = {"statsmodels": statsmodels.__version__, "numpy": np.__version__, "checks": rows,
       "all_pass": all(x["pass"] for x in rows)}
with open(os.path.join(VER, "statsmodels_check.json"), "w") as fh:
    json.dump(res, fh, indent=2)
for x in rows:
    print(f"{x['tag']:<22} est R {x['R_est']:.6f} sm {x['sm_est']:.6f} | SE R {x['R_se_hc3']:.6f} sm {x['sm_se_hc3']:.6f} | "
          f"{'OK' if x['pass'] else 'MISMATCH'}")
print("all pass:", res["all_pass"])
