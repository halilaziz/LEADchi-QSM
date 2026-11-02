#!/usr/bin/env python3
"""Compare results.json (and recovery_long.csv) against the independent recomputation in recomputed.json."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REC = HERE.parent
rj = json.loads((REC / "results.json").read_text())
me = json.loads((HERE / "recomputed.json").read_text())
rows = []


def ok(a, b):
    # results.json is written with 8 decimal digits: accept rel <= 1e-6 OR abs <= half a unit of the 8th decimal
    if isinstance(b, str) or isinstance(a, str):
        return a == b
    a, b = float(a), float(b)
    return abs(a - b) <= max(1e-6 * abs(b), 5.01e-9)


def add(path, rep, mine):
    rows.append((path, rep, mine, "OK" if ok(rep, mine) else "MISMATCH"))


lvmap = {"-20": "-20.0", "-14.24": "-14.24", "-10": "-10.0", "-5": "-5.0", "10": "10.0"}
fields = ["n", "mean", "sd", "median", "q25", "q75", "min", "max"]
for blk in ["k_all", "k_first_usable", "d_ref_ppb", "d_GP_ppb", "d_THAL_ppb"]:
    for L, v in rj[blk].items():
        for f, a, b in zip(fields, v, me[blk][lvmap[L]]):
            add(f"{blk}[{L}].{f}", a, b)
for f, a, b in zip(fields, rj["k_linear_per_session"], me["k_linear_per_session"]):
    add(f"k_linear_per_session.{f}", a, b)
for f in ["slope", "intercept", "r2", "n"]:
    add(f"pooled_linear.{f}", rj["pooled_linear"][f], me["pooled_linear"][f])
for g in ["PS", "ADHD"]:
    for f in ["est", "se", "t", "df", "p_two", "p_lower", "p_upper"]:
        add(f"k_by_group_L14.{g}.{f}", rj["k_by_group_L14"][g][f], me["k_by_group_L14"][g][f])
    for f in ["ci95", "ci90"]:
        for i in range(2):
            add(f"k_by_group_L14.{g}.{f}[{i}]", rj["k_by_group_L14"][g][f][i], me["k_by_group_L14"][g][f][i])
for i, c in enumerate(me["k_by_group_L14"]["cell_order"]):
    add(f"k_by_group_L14.cell_means[{c}]", rj["k_by_group_L14"]["cell_means"][i], me["k_by_group_L14"]["cell_means"][i])
for j, v in enumerate(["dice_min", "motion_gsr", "roi_vol", "Q_phase", "age"]):
    add(f"k_spearman_L14.rho[{v}]", rj["k_spearman_L14"][0][j], me["k_spearman_L14"][v]["rho"])
    add(f"k_spearman_L14.p[{v}]", rj["k_spearman_L14"][1][j], me["k_spearman_L14"][v]["p"])
o, mo = rj["observed"], me["observed"]
add("observed.beta_ppb", o["beta_ppb"], mo["beta"])
add("observed.se_hc3", o["se_hc3"], mo["se"])
add("observed.p_two", o["p_two"], mo["p_two"])
add("observed.ci95_ppb[0]", o["ci95_ppb"][0], mo["ci95"][0])
add("observed.ci95_ppb[1]", o["ci95_ppb"][1], mo["ci95"][1])
add("observed.sigma_B", o["sigma_B"], mo["sigma_B"])
add("observed.outcome", o["outcome"], mo["outcome"])
add("observed.bound_ppb", o["bound_ppb"], mo["bound"])
add("observed.lower95_one_sided", o["lower95_one_sided"], mo["ci90"][0])
for c, m in zip(rj["counterfactual"], me["counterfactual"]):
    L = c["L_true_ppb"]
    for f, g in [("L_true_ppb", "L_true_ppb"), ("n_missing_rec", "n_missing_rec"), ("mean_recovered_in_PS", "mean_recovered_in_PS"),
                 ("beta_ppb", "beta"), ("se", "se"), ("p_two", "p_two"), ("pB", "pB"), ("bound_ppb", "bound"), ("outcome", "outcome")]:
        add(f"counterfactual[{L}].{f}", c[f], m[g])
for f in ["k_mean_L14", "bound_true_ppb", "lower95_one_sided_true_ppb", "lower95_two_sided_true_ppb"]:
    add(f"true_scale.{f}", rj["true_scale"][f], me["true_scale"][f])
add("meta.rL", rj["meta"]["rL"], 0.938278557806515)
for i, (v, L) in enumerate([("rec_m05", -5), ("rec_m10", -10), ("rec_m14", -14.24), ("rec_m20", -20), ("rec_p10", 10)]):
    add(f"meta.levels[{i}].variant", rj["meta"]["levels"][i]["variant"], v)
    add(f"meta.levels[{i}].L", rj["meta"]["levels"][i]["L"], L)

# recovery_long.csv vs recomputed long table
import importlib.util
spec = importlib.util.spec_from_file_location("v", HERE / "verify_recovery.py")
v = importlib.util.module_from_spec(spec); spec.loader.exec_module(v)
_, _, mylong = v.load()
rl = pd.read_csv(REC / "recovery_long.csv")
m = rl.merge(mylong, on=["sess_key", "L"], suffixes=("_f", "_m"))
long_check = {"rows_file": len(rl), "rows_mine": len(mylong), "rows_matched": len(m)}
for c in ["d", "k", "d_ref", "d_GP", "d_THAL", "first_usable", "PS_e", "ADHD_e", "age"]:
    long_check[f"maxabs_{c}"] = float(np.nanmax(np.abs(m[c + "_f"].astype(float) - m[c + "_m"].astype(float))))
long_check["cell_equal"] = bool((m.cell_f == m.cell_m).all())

n_ok = sum(r[3] == "OK" for r in rows)
lines = [f"{p}\t{a}\t{b}\t{s}" for p, a, b, s in rows]
(HERE / "results_json_check.tsv").write_text("path\treported\trecomputed\tstatus\n" + "\n".join(lines) + "\n")
print(f"{n_ok}/{len(rows)} OK")
for r in rows:
    if r[3] != "OK":
        print("MISMATCH", r)
print(json.dumps(long_check, indent=1))
(HERE / "long_check.json").write_text(json.dumps(long_check, indent=1))
