#!/usr/bin/env python3
"""Compare results.json and the README's printed numbers with verification/recomputed.json.
Writes verification/number_checks.md (a table pasted into VERIFY.md)."""
import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
E2 = os.path.dirname(HERE)
R = json.load(open(os.path.join(E2, "results.json")))
V = json.load(open(os.path.join(HERE, "recomputed.json")))
rows = []


def rel_ok(a, b, tol=1e-6):
    """OK if relative difference <= 1e-6, or within the 8-decimal rounding used when results.json was written."""
    return abs(a - b) <= tol * abs(b) or abs(a - b) <= 0.5e-8 + 1e-15


def add_json(path, rep, mine):
    if isinstance(rep, str) or isinstance(rep, bool):
        ok = rep == mine
    else:
        ok = rel_ok(float(rep), float(mine))
    rows.append(("results.json", path, rep, mine, "OK" if ok else "MISMATCH"))


def add_readme(label, printed, mine, dec=None, sig=None):
    if isinstance(printed, str) and not any(ch.isdigit() for ch in printed):
        ok = printed == mine
    elif sig is not None:
        ok = float(f"{mine:.{sig}g}") == float(printed)
    else:
        ok = abs(round(float(mine), dec) - float(printed)) < 10 ** (-dec) / 10
    rows.append(("README", label, printed, mine, "OK" if ok else "MISMATCH"))


# ---------------- results.json
p, vp = R["primary"], V["primary"]
for k_rep, k_me in [("beta_ppb", "beta"), ("se_hc3", "se"), ("p_two", "p_two"), ("sigma_resid", "sigma_resid"),
                    ("sigma_B", "sigma_B"), ("delta", "delta"), ("outcome", "outcome")]:
    add_json(f"primary.{k_rep}", p[k_rep], vp[k_me])
for i in range(2):
    add_json(f"primary.ci95_ppb[{i}]", p["ci95_ppb"][i], vp["ci95"][i])
for inv in ("e2_rts", "e2_tv"):
    r, v = R[inv], V[inv]
    for k in r["agreement"]:
        add_json(f"{inv}.agreement.{k}", r["agreement"][k], v["agreement"][k])
    for k_rep, k_me in [("beta_ppb", "beta"), ("se_hc3", "se"), ("p_two", "p_two"), ("sigma_resid", "sigma_resid"),
                        ("sigma_B", "sigma_B"), ("delta", "delta"), ("outcome", "outcome"), ("disconfirming", "disconfirming")]:
        add_json(f"{inv}.M1.{k_rep}", r["M1"][k_rep], v["M1"][k_me])
    for i in range(2):
        add_json(f"{inv}.M1.ci95_ppb[{i}]", r["M1"]["ci95_ppb"][i], v["M1"]["ci95"][i])
    add_json(f"{inv}.bound_ppb_registered_deltaB", r["bound_ppb_registered_deltaB"], v["M1"]["bound"])
    add_json(f"{inv}.pB", r["pB"], v["M1"]["pB"])
    add_json(f"{inv}.huber_beta", r["huber_beta"], v["M1"]["huber_beta"])
    for c in ("PSonly_minus_TD", "PSADHD_minus_ADHD"):
        for k in ("est", "se", "p_two"):
            add_json(f"{inv}.K1.{c}.{k}", r["K1"][c][k], v["M1"][c][k])
        for i in range(2):
            add_json(f"{inv}.K1.{c}.ci95[{i}]", r["K1"][c]["ci95"][i], v["M1"][c]["ci95"][i])
    add_json(f"{inv}.n_dropped_invalid", r["n_dropped_invalid"], 80 - v["agreement"]["n_valid_e2"])
    for s in ("all", "first_usable"):
        for i, nm in enumerate(["n", "mean", "sd", "median"]):
            add_json(f"{inv}.recovery_k_L14.{s}.{nm}", r["recovery_k_L14"][s][i], v["recovery_k_L14"][s][i])

# ---------------- README
rt, tv = V["e2_rts"], V["e2_tv"]
A = lambda x: x["agreement"]
add_readme("Agreement n", "79", rt["agreement"]["n_valid_e2"], dec=0)
for lab, x in (("RTS", rt), ("TV", tv)):
    add_readme(f"ICC(3,1) {lab}", {"RTS": "0.68", "TV": "0.69"}[lab], A(x)["icc31"], dec=2)
    add_readme(f"Pearson {lab}", "0.73", A(x)["pearson"], dec=2)
    add_readme(f"Spearman {lab}", "0.71", A(x)["spearman"], dec=2)
    add_readme(f"Mean chi {lab}", {"RTS": "27.1", "TV": "24.0"}[lab], A(x)["mean_e2"], dec=1)
    add_readme(f"SD {lab}", {"RTS": "8.1", "TV": "8.7"}[lab], A(x)["sd_e2"], dec=1)
    add_readme(f"GP Pearson {lab}", {"RTS": "0.89", "TV": "0.90"}[lab], A(x)["gp_pearson"], dec=2)
add_readme("Mean chi primary", "53.0", A(rt)["mean_primary"], dec=1)
add_readme("SD primary", "12.0", A(rt)["sd_primary"], dec=1)
cols = {"Primary": V["primary"], "RTS": rt["M1"], "TV": tv["M1"]}
printed = {
    "Primary": dict(b="4.96", lo="-1.27", hi="11.20", p=".117", d="0.44", bd="-6.13", pB=".0004", hub="3.94",
                    s1="8.16", s1p=".044", s2="1.77", s2p=".71", s2pB=".052"),
    "RTS": dict(b="3.81", lo="-0.51", hi="8.14", p=".083", d="0.50", bd="-3.92", pB=".0003", hub="2.97",
                s1="3.56", s1p=".23", s2="4.07", s2p=".20", s2pB=".006"),
    "TV": dict(b="4.37", lo="-0.25", hi="8.99", p=".063", d="0.53", bd="-4.48", pB=".0001", hub="3.34",
               s1="4.52", s1p=".15", s2="4.22", s2p=".21", s2pB=".006")}
for lab, x in cols.items():
    pr = printed[lab]
    add_readme(f"beta_PS {lab}", pr["b"], x["beta"], dec=2)
    add_readme(f"CI lower {lab}", pr["lo"], x["ci95"][0], dec=2)
    add_readme(f"CI upper {lab}", pr["hi"], x["ci95"][1], dec=2)
    add_readme(f"p two-sided {lab}", pr["p"], x["p_two"], dec=len(pr["p"]) - 1)
    add_readme(f"delta {lab}", pr["d"], x["delta"], dec=2)
    add_readme(f"bound ppb {lab}", pr["bd"], x["bound"], dec=2)
    add_readme(f"p_B {lab}", pr["pB"], x["pB"], dec=4)
    add_readme(f"Outcome {lab}", "B", x["outcome"])
    add_readme(f"Huber beta {lab}", pr["hub"], x["huber_beta"], dec=2)
    add_readme(f"PS-only - TD est {lab}", pr["s1"], x["PSonly_minus_TD"]["est"], dec=2)
    add_readme(f"PS-only - TD p {lab}", pr["s1p"], x["PSonly_minus_TD"]["p_two"], dec=len(pr["s1p"]) - 1)
    add_readme(f"PS+ADHD - ADHD-only est {lab}", pr["s2"], x["PSADHD_minus_ADHD"]["est"], dec=2)
    add_readme(f"PS+ADHD - ADHD-only p {lab}", pr["s2p"], x["PSADHD_minus_ADHD"]["p_two"], dec=len(pr["s2p"]) - 1)
    add_readme(f"comorbid-stratum p_B {lab}", pr["s2pB"], x["PSADHD_minus_ADHD"]["pB_one_sided"], dec=3)
    add_readme(f"comorbid-stratum bound rejected {lab}", {"Primary": "no", "RTS": "yes", "TV": "yes"}[lab],
               "yes" if x["PSADHD_minus_ADHD"]["pB_one_sided"] < 0.05 else "no")
add_readme("delta_B", "0.581", V["primary"]["delta_B"], dec=3)
kp = V["primary"]["recovery_k_L14"]["first_usable"]
add_readme("k primary mean", "0.47", kp[1], dec=2)
add_readme("k primary SD", "0.05", kp[2], dec=2)
for lab, x, m_, s_ in (("RTS", rt, "0.36", "0.04"), ("TV", tv, "0.25", "0.06")):
    add_readme(f"k {lab} mean", m_, x["recovery_k_L14"]["first_usable"][1], dec=2)
    add_readme(f"k {lab} SD", s_, x["recovery_k_L14"]["first_usable"][2], dec=2)
add_readme("expected adult deficit primary", "-6.7", V["primary"]["adult_test_recomputed"]["expected_ppb"], dec=1)
add_readme("expected adult deficit RTS", "-5.2", rt["adult_test"]["expected_ppb"], dec=1)
add_readme("expected adult deficit TV", "-3.5", tv["adult_test"]["expected_ppb"], dec=1)
add_readme("adult-deficit p primary (given inputs)", ".0002", V["primary"]["adult_test_given"]["p_one_sided"], dec=4)
add_readme("adult-deficit p primary (recomputed inputs)", ".0002", V["primary"]["adult_test_recomputed"]["p_one_sided"], dec=4)
add_readme("adult-deficit p RTS", "5e-5", rt["adult_test"]["p_one_sided"], sig=1)
add_readme("adult-deficit p TV", "5e-4", tv["adult_test"]["p_one_sided"], sig=1)
add_readme("Reading: strata RTS +3.6", "3.6", rt["M1"]["PSonly_minus_TD"]["est"], dec=1)
add_readme("Reading: strata RTS +4.1", "4.1", rt["M1"]["PSADHD_minus_ADHD"]["est"], dec=1)
add_readme("Reading: primary PS-only - TD +8.2", "8.2", V["primary"]["PSonly_minus_TD"]["est"], dec=1)
add_readme("Reading: ICC approx 0.68 (RTS)", "0.68", rt["agreement"]["icc31"], dec=2)
add_readme("Reading: delta approx +0.5 (RTS)", "0.5", rt["M1"]["delta"], dec=1)
add_readme("Reading: delta approx +0.5 (TV)", "0.5", tv["M1"]["delta"], dec=1)
add_readme("'about half' SN-VTA contrast RTS (ratio)", "0.5", rt["contrast_ratio"]["chi_SNVTA"]["ratio"], dec=1)
add_readme("'about half' SN-VTA contrast TV (ratio)", "0.5", tv["contrast_ratio"]["chi_SNVTA"]["ratio"], dec=1)


def fmt(x):
    if isinstance(x, float):
        return f"{x:.8g}"
    return str(x)


lines = ["| source | item | reported | recomputed | check |", "|---|---|---|---|---|"]
for s, k, a, b, c in rows:
    lines.append(f"| {s} | {k} | {fmt(a)} | {fmt(b)} | {c} |")
open(os.path.join(HERE, "number_checks.md"), "w").write("\n".join(lines) + "\n")
n_bad = sum(r[4] != "OK" for r in rows)
print(f"{len(rows)} checks, {n_bad} not OK")
for r in rows:
    if r[4] != "OK":
        print(r)
