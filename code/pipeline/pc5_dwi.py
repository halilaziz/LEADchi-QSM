#!/usr/bin/env python3
"""pc5_dwi.py - gate PC5 (DWI positive control for S3; registration Q24) on the real-ID DWI tables.

Runs after dwi_geometry.py (stage geometry) on the DWI side (non-blind by construction; no chi, no group label).
  PC5(a) on each subject's first PennLINC scalar-pass DWI session with status ok: pSCP_R tw-ICVF and whole-tract
         MCP_R tw-ICVF each rise with age, one-sided p < .05 for the age slope in OLS y ~ age + sex + FD + NDC.
  PC5(b) over the retained pSCP_R node observations of the same sessions: |r(ISOVF, tw-ICVF)| < |r(ISOVF, ICVF)|.
  PC5 passes if (a) and (b) both hold. Output: dwi/pc5.json (aggregate only).
Covariates: --covariates is a CSV from which ONLY the columns (participant_id|sub, session_id|ses, age, sex) are read;
the script refuses to proceed if it would need any other column, and it never reads a group or diagnosis column.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lx_common as L  # noqa: E402

WHITELIST = {"participant_id", "sub", "session_id", "ses", "age", "sex"}


def read_covariates(path):
    out = {}
    with open(path, newline="", encoding="utf-8-sig") as f:
        rd = csv.reader(f)
        head = next(rd)
        keep = [i for i, h in enumerate(head) if h in WHITELIST]
        names = [head[i] for i in keep]
        for row in rd:
            d = {names[j]: row[i] for j, i in enumerate(keep)}      # only whitelisted columns are ever held
            sub = d.get("participant_id") or d.get("sub")
            ses = d.get("session_id") or d.get("ses")
            ses = ses if str(ses).startswith("ses-") else f"ses-{ses}"
            sub = sub if str(sub).startswith("sub-") else f"sub-{sub}"
            out[(sub, ses)] = d
    return out


def ols_slope_p(y, X, j):
    import numpy as np
    from scipy import stats
    X = np.column_stack([np.ones(len(y))] + X)
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ b
    df = len(y) - X.shape[1]
    s2 = float(r @ r) / df
    se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
    t = b[j + 1] / se[j + 1]
    return float(b[j + 1]), float(t), float(stats.t.sf(t, df)), int(len(y))


def main():
    import numpy as np
    a = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    a.add_argument("--covariates", required=True, help="CSV with participant_id/sub, session_id/ses, age, sex (other columns ignored)")
    a.add_argument("--dry-run", action="store_true")
    args = a.parse_args()
    lx = L.lx_root()
    d = os.path.join(lx, "dwi")
    out = os.path.join(d, "pc5.json")
    L.say(f"PC5: {d}/dwi_sessions.csv + dwi_nodes.csv + covariates -> {out}", args.dry_run)
    if args.dry_run:
        return
    cov = read_covariates(args.covariates)
    sess = [r for r in L.read_csv(os.path.join(d, "dwi_sessions.csv")) if r["status"] == "ok" and r["scalar_pass"] == "1"]
    sess.sort(key=lambda r: (r["participant_id"], int(r["session_id"].split("-")[1])))
    first = {}
    for r in sess:
        first.setdefault(r["participant_id"], r)
    rows = [r for r in first.values() if (r["participant_id"], r["session_id"]) in cov]
    res = {}
    for y in ("twICVF_pSCP_R", "twICVF_MCP_R_whole"):
        ok = [r for r in rows if r.get(y) not in ("", "NA") and r.get("FD") not in ("", "NA") and r.get("NDC") not in ("", "NA")]
        c = [cov[(r["participant_id"], r["session_id"])] for r in ok]
        yy = np.array([float(r[y]) for r in ok])
        X = [np.array([float(x["age"]) for x in c]), np.array([1.0 if x["sex"] == "M" else 0.0 for x in c]),
             np.array([float(r["FD"]) for r in ok]), np.array([float(r["NDC"]) for r in ok])]
        b, t, p, n = ols_slope_p(yy, X, 0)
        res[y] = {"n": n, "age_slope_per_y": b, "t": t, "p_one_sided": p, "pass": p < 0.05 and b > 0}
    keys = {(r["participant_id"], r["session_id"]) for r in rows}
    nd = [n for n in L.read_csv(os.path.join(d, "dwi_nodes.csv"))
          if n["bundle"] == "SCP_R" and n.get("retained") == "1" and (n["participant_id"], n["session_id"]) in keys]
    iso = np.array([float(n["isovf"]) for n in nd]); tw = np.array([float(n["twicvf"]) for n in nd])
    ic = np.array([float(n["icvf"]) for n in nd])
    r_tw, r_ic = float(np.corrcoef(iso, tw)[0, 1]), float(np.corrcoef(iso, ic)[0, 1])
    b_ok = abs(r_tw) < abs(r_ic)
    summ = {"utc": L.utcnow(), "PC5a": res, "PC5b": {"n_node_obs": len(nd), "r_isovf_twicvf": r_tw, "r_isovf_icvf": r_ic, "pass": b_ok},
            "PC5": bool(all(v["pass"] for v in res.values()) and b_ok)}
    L.write_json(out, summ)
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
