#!/usr/bin/env python3
"""make_results.py - POST HOC (2026-10-28): assemble results.json and README.md of the complete-human-QC sample.

Reads only outputs inside results_posthoc/humanqc_primary/ plus frozen reference files (read-only):
  frozen_run_H_all180/results.json       frozen run_all.R (unchanged) on input/merged_H_all180.csv
  supplement_H_all180.json               code/supplement_humanqc.R on the same table
  sanity_frozen_run_AI/results.json      frozen run_all.R on the frozen AI-QC table (reproduction check)
  sanity_supplement_AI/supplement_AI_frozen.json   supplement on the frozen AI-QC table (Table 1 check)
  results_confirmatory/results.json, kit2_m1_refit.json, engine2/results.json   (frozen references)
Writes ../results.json and ../README.md only.
"""
import json, math, os, sys
from decimal import Decimal, ROUND_HALF_UP

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.normpath(os.path.join(HERE, ".."))
LX = os.path.normpath(os.path.join(BASE, "../.."))
sys.path.insert(0, HERE); sys.dont_write_bytecode = True
import build_input  # noqa: E402

J = lambda p: json.load(open(p))
FR = J(f"{BASE}/frozen_run_H_all180/results.json")
SUP = J(f"{BASE}/supplement_H_all180.json")
FRA = J(f"{BASE}/sanity_frozen_run_AI/results.json")
SUPA = J(f"{BASE}/sanity_supplement_AI/supplement_AI_frozen.json")
CONF = J(f"{LX}/results_confirmatory/results.json")
KIT2 = J(f"{LX}/human_qc_kit/kit2_2026-10-28/analysis/real_run_20261028/kit2_m1_refit.json")
E2F = J(f"{LX}/results_posthoc/engine2/results.json")
CELLS = ["TD", "ADHD-only", "PS-only", "PS+ADHD"]
FROZEN_BOUND = -6.13180464348368


def r(x, d=1):
    """round half up (the manuscript's convention: 10.45 -> 10.5)"""
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return None
    return float(Decimal(repr(x)).quantize(Decimal(1).scaleb(-d), rounding=ROUND_HALF_UP))


def fp(p):
    if p is None: return "NA"
    if p < 0.0001: return f"{p:.1e}"
    return f"{p:.4f}"


def deep_equal(a, b, path="", out=None):
    out = [] if out is None else out
    if isinstance(a, dict) and isinstance(b, dict):
        for k in set(a) | set(b):
            if path == "" and k == "meta": continue
            if k not in a or k not in b: out.append(path + "/" + k); continue
            deep_equal(a[k], b[k], path + "/" + k, out)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b): out.append(path); return out
        for i, (u, v) in enumerate(zip(a, b)): deep_equal(u, v, f"{path}[{i}]", out)
    elif isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool):
        if not (a == b or abs(a - b) <= 1e-12 * max(1, abs(b))): out.append(path)
    elif a != b:
        out.append(path)
    return out


# ------------------------------------------------------------------------------------------------ sanity checks
checks = {}
txt, raw = build_input.build("AI_frozen")
checks["input_rebuild_AI_frozen_byte_identical_to_frozen_merged_table"] = txt == raw
d = deep_equal(FRA, CONF)
checks["frozen_run_all_on_AI_table_reproduces_results_confirmatory_json_exactly"] = {"ok": not d, "n_differing_fields": len(d)}
checks["AI_n80_beta4.9628"] = FRA["H1"]["n"] == 80 and abs(FRA["H1"]["decision"]["beta_ppb"] - 4.96280154409223) < 1e-12
T1_MS = {  # manuscript/LEADchi_manuscript_source.md Table 1 (n = 80), as printed
    "youths_with_megre": [30, 22, 24, 37], "usable_sessions": [31, 19, 22, 35], "analysed_first_usable": [22, 16, 17, 25],
    "earliest_megre_not_usable": [10, 12, 8, 19], "analysed_at_ses2_or_ses3": [3, 6, 2, 11],
    "female": [9, 9, 9, 11], "male": [13, 7, 8, 14], "race_black": [9, 10, 11, 15], "race_white": [9, 6, 4, 7],
    "race_other": [4, 0, 2, 3], "threshold_psychotic_disorder": [0, 0, 0, 3],
    "age": [(13.1, 2.2, 8.8, 17.8), (12.8, 2.4, 9.1, 16.9), (14.0, 1.7, 11.1, 16.0), (14.1, 1.7, 9.7, 17.0)],
    "PRIME_total": [(10.5, 14.5, 20), (14.1, 13.2, 15), (17.8, 14.6, 17), (21.7, 16.3, 24)],
    "BDI_total": [(4.5, 3.0, 19), (6.2, 6.7, 15), (3.9, 3.2, 17), (10.6, 10.8, 24)],
    "ARI_total": [(2.7, 2.6, 22), (3.9, 2.7, 16), (3.1, 2.7, 17), (4.8, 4.0, 25)],
    "chi_raw": [(50.1, 12.6), (50.9, 13.8), (56.8, 10.0), (54.3, 10.8)],
    "chi_adj": [(49.7, 44.1, 55.3), (51.6, 43.5, 59.8), (57.9, 52.3, 63.5), (53.4, 48.9, 57.9)]}
TA = SUPA["table1"]; bad = []
for k in ["youths_with_megre", "usable_sessions", "analysed_first_usable", "earliest_megre_not_usable", "analysed_at_ses2_or_ses3",
          "female", "male", "race_black", "race_white", "race_other", "threshold_psychotic_disorder"]:
    if [TA[k][c] for c in CELLS] != T1_MS[k]: bad.append(k)
if [(r(x["mean"]), r(x["sd"]), r(x["min"]), r(x["max"])) for x in TA["age"]] != T1_MS["age"]: bad.append("age")
for v in ["PRIME_total", "BDI_total", "ARI_total"]:
    if [(r(x["mean"]), r(x["sd"]), x["n"]) for x in TA[v]] != T1_MS[v]: bad.append(v)
if [(r(x["mean"]), r(x["sd"])) for x in TA["chi_SNVTA_raw"]] != T1_MS["chi_raw"]: bad.append("chi_raw")
if [(r(x["adj_mean"]), r(x["lo"]), r(x["hi"])) for x in TA["chi_SNVTA_adjusted"]] != T1_MS["chi_adj"]: bad.append("chi_adj")
if not (abs(TA["qc_by_group_fisher_p_not_in_H1"] - 0.9670) < 5e-5 and abs(TA["age_PS_welch_p"] - 0.018) < 5e-4
        and r(TA["age_PS_minus_nonPS_years"]) == 1.1): bad.append("Fisher/age text")
if TA["exclusion_reason_component_counts"]["c_qc_grade"] != 62 or TA["exclusion_reason_component_counts"]["e_g5"] != 45: bad.append("exclusion 62/45")
checks["AI_pipeline_reproduces_manuscript_Table1_and_Results_sample_text"] = {"ok": not bad, "mismatches": bad}
sb = {s["sensitivity"]: s for s in SUPA["sensitivity"]["battery_frozen_function"]}
mx = max(abs((sb[s["sensitivity"]][k] or 0) - (s[k] or 0)) for s in CONF["sensitivity"]["table"]
         for k in ("beta_ppb", "se_hc3", "p_two", "pB") if s["sensitivity"] in sb)
checks["AI_sensitivity_battery_with_df_wrapper_equals_frozen_battery_max_abs_diff"] = mx
e2a = SUPA["engine2"]
checks["AI_engine2_reproduces_engine2_results_json"] = all(
    abs(e2a[i]["main"]["M1"]["beta_ppb"] - E2F[f"e2_{i}"]["main"]["M1"]["beta_ppb"]) < 1e-6 and
    e2a[i]["main"]["n"] == E2F[f"e2_{i}"]["main"]["n"] and abs(e2a[i]["main"]["pB"] - E2F[f"e2_{i}"]["main"]["pB"]) < 1e-8
    for i in ("rts", "tv"))
k2 = KIT2["scenarios"]["H_all180"]
checks["H_all180_reproduces_kit2_m1_refit"] = {
    "ok": FR["H1"]["n"] == k2["n"] == 85 and abs(FR["H1"]["decision"]["beta_ppb"] - k2["beta_PS_ppb"]) < 1e-12
          and sorted(SUP["table1"]["first_usable_keys"]) == sorted(k2["first_usable_keys"])
          and abs(SUP["H1_extras"]["bound_frozen"]["pB"] - k2["bound_fixed"]["pB_one_sided"]) < 1e-15,
    "n": FR["H1"]["n"], "beta_PS": FR["H1"]["decision"]["beta_ppb"], "kit2_beta": k2["beta_PS_ppb"]}
checks["supplement_M1_equals_frozen_run_M1_H_all180"] = abs(SUP["H1_extras"]["beta_PS"] - FR["H1"]["decision"]["beta_ppb"]) < 1e-12

# ------------------------------------------------------------------------------------------------ assemble
H, X, T, S2, S1 = FR["H1"], SUP["H1_extras"], SUP["table1"], FR["S2"], FR["S1"]
dec = H["decision"]
battery = SUP["sensitivity"]["battery_frozen_function"]
pick = {b["sensitivity"]: b for b in battery}
fsens = {b["sensitivity"]: b for b in FR["sensitivity"]["table"]}


def srow(b, label, src):
    return {"label": label, "n": b["n"], "beta_ppb": b["beta_ppb"], "se_hc3": b["se_hc3"], "df": b.get("df"), "p_two": b["p_two"],
            "pB_frozen_bound": b.get("pB_frozen_bound"), "bound_recomputed_ppb": b.get("bound_recomputed_ppb"),
            "pB_recomputed_bound": b["pB"], "outcome_recomputed": b["outcome"], "source": src}


SRC_B = "supplement_H_all180.json sensitivity.battery_frozen_function (frozen R/sensitivity.R run_h1_sensitivities; identical to frozen_run_H_all180/results.json sensitivity.table + df and frozen-bound p_B)"
sens_rows = [
    srow(SUP["sensitivity"]["ses1_only"][0], "session-1-only first-usable",
         "supplement sensitivity.ses1_only (frozen run_m1_core on S$ses1_only; = frozen_run H1.ses1_only)"),
    srow(pick["without the 6 threshold-psychosis youths"], f"without threshold-psychosis youths ({SUP['sensitivity']['threshold_psychosis_removed']} removed; frozen label says 6)", SRC_B),
    srow(pick["without prescribed-stimulant endorsers (same session)"], f"without prescribed-stimulant endorsers ({SUP['sensitivity']['stimulant_endorsers_removed']} removed)", SRC_B),
    srow(pick["+ registration Dice (minimum of the G5 Dice values)"], "+ registration Dice", SRC_B),
    srow(pick["+ magnitude ghost-to-signal motion index"], "+ motion index", SRC_B),
    srow(SUP["sensitivity"]["joint_nuisance"][0], "+ Dice + motion + CSF offset + ICV (joint; post hoc)",
         "supplement sensitivity.joint_nuisance (frozen fit_ols_hc3/huber_fit/h1_decide; extra terms as validation_local.R B.1)"),
    srow(pick["SN-only complex (descriptive)"], "SN-only complex (descriptive)", SRC_B),
    srow(pick["VTA + PBP only (descriptive)"], "VTA + PBP only (descriptive)", SRC_B),
]
loo = H["loo"]
res = {
    "label": "POST HOC: LEAD-chi analyses on the complete-human-QC first-usable sample (scenario H_all180; first-author "
             "visual QC of all 180 graded sessions, kit 1 + kit 2). Not confirmatory; the registered outcome (AI-QC, n = 80) stays B.",
    "inputs": {"merged_table": "input/merged_H_all180.csv", "merged_sha256": SUP["merged_sha256"],
               "scenario_source": "human_qc_kit/kit2_2026-10-28/analysis/real_run_20261028/kit2_scenario_usable.csv",
               "gates": "mirror/gates/gates_preunlock.json (frozen; rL_primary 0.9383, PC4 TRUE; S1 pooled constants)",
               "multiverse": "mirror/unlocked/multiverse_qsmxt_SNVTA.csv (frozen)",
               "replaced_columns": "usable, first_usable, fail_reason <- H_all180; qc_grade <- human grade (AI grade kept as qc_grade_ai)"},
    "sanity_checks": checks,
    "1_sample": {
        "cells_analysed": T["analysed_first_usable"], "n": T["n_analysed"],
        "sessions_usable_of_181": T["sessions_usable"], "sessions_total": T["sessions_total"],
        "usable_sessions_by_cell": T["usable_sessions"], "youths_with_megre": T["youths_with_megre"],
        "earliest_megre_not_usable": T["earliest_megre_not_usable"],
        "exclusion_reason_component_counts_sessions": T["exclusion_reason_component_counts"],
        "exclusion_reason_combinations": T["exclusion_reason_combinations"],
        "analysed_at_ses2_or_ses3": T["analysed_at_ses2_or_ses3"], "analysed_session": T["analysed_session"],
        "age": T["age"], "female": T["female"], "male": T["male"],
        "race_black": T["race_black"], "race_white": T["race_white"], "race_other": T["race_other"],
        "threshold_psychotic_disorder": T["threshold_psychotic_disorder"],
        "chi_SNVTA_raw": T["chi_SNVTA_raw"], "chi_SNVTA_M1_adjusted": T["chi_SNVTA_adjusted"],
        "PRIME_total": T["PRIME_total"], "BDI_total": T["BDI_total"], "ARI_total": T["ARI_total"],
        "fisher_p_not_in_H1": T["qc_by_group_fisher_p_not_in_H1"], "fisher_p_earliest_failed": T["qc_by_group_fisher_p_earliest_failed"],
        "age_PS_minus_nonPS_years": T["age_PS_minus_nonPS_years"], "age_PS_welch_p": T["age_PS_welch_p"],
        "robust_z_flags": FR["plausibility"]["robust_z"]},
    "2_H1": {
        "n": H["n"], "df": dec["df"], "beta_PS_ppb": dec["beta_ppb"], "se_hc3": dec["se_hc3"], "ci95": dec["ci95_ppb"], "ci90": dec["ci90_ppb"],
        "p_two": dec["p_two"], "one_sided_95_lower_limit": X["one_sided_95_lower_limit"], "delta": dec["delta"],
        "delta_boot_ci95": H["boot_delta"]["ci95"], "delta_nct_ci95": H["nct_ci_delta"], "hedges_g_raw": H["hedges_g_raw"],
        "sigma_resid": dec["sigma_resid"], "sigma_rob": dec["sigma_rob"], "sigma_B": dec["sigma_B"],
        "freedman_lane": H["freedman_lane"], "huber": H["huber"],
        "bound_test_frozen_bound": X["bound_frozen"], "bound_test_recomputed": dec["B_registered"],
        "B_coreports_recomputed_sigma": dec["B_coreport"], "B_coreport_cross_session": H["B_coreport_cross_session"],
        "B_coreport_residualised_split": H["B_coreport_residualised_split"],
        "outcome_M1": dec["outcome"], "final_outcome_after_K3": H["final_outcome"], "disconfirming_rule": dec["disconfirming"],
        "K3": {**H["K3"], "csf_se": X["K3_csf"]["se"], "csf_ci95": X["K3_csf"]["ci95"], "csf_pB_frozen_bound": X["K3_csf"]["pB_frozen_bound"]},
        "bayes": H["bayes"], "adult_reference_ppb": H["ppb"]["adult_reference"]},
    "3_K1_K2": {"K1": X["K1"], "K2_R2s": {"one_sided_lower_p": H["K2"]["p_one_sided_lower"], "est": H["K2"]["est"], "se": H["K2"]["se"],
                                          "ci95": H["K2"]["ci95"], "p_two": SUP["fig3B_r2s"]["p_two"], "n": H["K2"]["n"]}},
    "4_sensitivity": {"rows": sens_rows,
                      "leave_one_out": {"n_refits": loo["n_refits"], "n_outcome_changes": loo["n_category_change"], "fragile": loo["fragile"]},
                      "influence": H["influence"]["diagnostics"], "exposure_sensitivity": H["exposure_sensitivity"],
                      "full_battery": battery,
                      "spec_curve": {k: FR["spec_curve"][k] for k in ("n_specifications", "n_evaluable", "assessable", "median_delta",
                                                                         "share_delta_negative", "robust", "robust_reason")}},
    "5_secondary": {"S1_L": {k: S1["L"][k] for k in ("est", "se", "t", "df", "p_two", "ci95")}, "S1_n": S1["n"],
                    "S2_PRIME_per_SD": {k: S2["beta_PRIME"][k] for k in ("est", "se", "t", "df", "p_two", "ci95")}, "S2_n": S2["n"],
                    "holm": FR["holm"], "S2_plus_PS_e": S2["plus_PS_e"] if "plus_PS_e" in S2 else None},
    "6_PC2x": SUP["PC2x"],
    "7_figures": {"fig3A_regional": SUP["fig3A_regional"], "fig3B_r2s": SUP["fig3B_r2s"]},
    "8_engine2": SUP["engine2"],
    "other_registered_outputs_not_requested": {"T1_label": FR["T1"]["label"], "T1_followup_only": FR["T1"]["followup_only"],
                                               "pooled_all_session": FR["pooled"], "qc_ipw": FR["qc"]["ipw"]},
}
json.dump(res, open(f"{BASE}/results.json", "w"), indent=1, allow_nan=False, default=lambda o: None)

# ------------------------------------------------------------------------------------------------ README
c = lambda dct, f=str: " / ".join(f(dct[k]) for k in CELLS)
mean_sd_rng = lambda L: " / ".join(f"{r(x['mean'])} ({r(x['sd'])}; {r(x['min'])}–{r(x['max'])})" for x in L)
mean_sd_n = lambda L: " / ".join(f"{r(x['mean'])} ({r(x['sd'])}; n {x['n']})" for x in L)
L = []
A = L.append
A("# LEAD-χ on the complete-human-QC sample (H_all180): manuscript numbers\n")
A("**POST HOC.** This is the complete-human-QC sample: the first author graded all 180 graded sessions (kit 1 + kit 2). "
  "It is not the registered analysis, which used the AI QC rater (n = 80) and gave outcome B. Nothing here re-decides H1. "
  "The rater was not independent of the study.\n")
A("## How it was computed\n")
A("1. `code/build_input.py` builds `input/merged_H_all180.csv` from the frozen merged table. It replaces only four columns, with the "
  "H_all180 values from `kit2_scenario_usable.csv`: `usable`, `first_usable`, `fail_reason` and `qc_grade` (the human grade; the AI "
  "grade is kept as `qc_grade_ai`). Check: rebuilding with the AI_frozen columns gives the frozen table byte for byte.")
A("2. **`frozen_run_H_all180/`**: the frozen `mirror/scripts/analysis/run_all.R` was run **unchanged** on that table, with the frozen "
  "gates and multiverse files and registered draw counts. It writes `report.md` and `results.json`.")
A("3. `code/supplement_humanqc.R` (output `supplement_H_all180.json`) computes what `run_all.R` does not print. It sources the frozen "
  "helpers read-only. The non-frozen pieces are listed in its header:")
A("   - Table 1 descriptives;")
A("   - stratum bound tests;")
A("   - the joint nuisance row;")
A("   - Figure 3 regional rows and adjusted means (helpers copied verbatim from `tier4_exploratory.R`);")
A("   - the engine-2 refit (`engine2_analysis.R` logic).")
A("4. `code/figs_humanqc.py` is a copy of `manuscript/code/figs_pcn.py` with the same plotting code. It reads only the files above and "
  "runs 26 assertions. It writes `figures/Figure{1,2,3}.tif` (600 dpi) and `.jpg` (300 dpi).")
A("5. `code/make_results.py` assembles `results.json` and this README. It also re-runs the sanity checks.\n")
A("Path abbreviations used below:")
A("- **FR** = `frozen_run_H_all180/results.json` (frozen `R/m1.R`, `secondary.R`, `sensitivity.R`, `spec_curve.R`);")
A("- **SUP** = `supplement_H_all180.json`.\n")
A("## Sanity checks\n")
A("| check | result |\n|---|---|")
for k, v in checks.items():
    A(f"| {k} | {json.dumps(v) if not isinstance(v, bool) else v} |")
A("\nTable 1 values are rounded half up, as in the manuscript (for example, TD PRIME 10.45 is shown as 10.5).\n")

A("## 1. Sample (cells: TD / ADHD-only / PS-only / PS+ADHD)\n")
A("| quantity | value | code path |\n|---|---|---|")
A(f"| analysed (first usable) | {c(T['analysed_first_usable'])} (n = {T['n_analysed']}) | FR samples; SUP table1.analysed_first_usable (frozen build_samples) |")
A(f"| usable sessions | {T['sessions_usable']} of {T['sessions_total']} ({c(T['usable_sessions'])}) | SUP table1.sessions_usable |")
A(f"| youths with ≥ 1 MEGRE session | {c(T['youths_with_megre'])} | SUP table1.youths_with_megre |")
A(f"| earliest MEGRE session not usable | {c(T['earliest_megre_not_usable'])} | frozen qc_by_group() table_earliest |")
A(f"| analysed at ses-2 or ses-3 | {c(T['analysed_at_ses2_or_ses3'])} (ses-1 {T['analysed_session'].get('ses-1')}, ses-2 {T['analysed_session'].get('ses-2')}, ses-3 {T['analysed_session'].get('ses-3', 0)}) | SUP table1 |")
ex = T["exclusion_reason_component_counts"]
A(f"| exclusion reasons (sessions; a session can have several) | registration/G5 (e_g5) {ex.get('e_g5')}; visual grade 2 or ungradable (c_qc_grade) {ex.get('c_qc_grade')}; "
  f"placement fail (d_placement) {ex.get('d_placement')}; (c') plausibility {ex.get('cprime_plausibility')}; automated grade 2 {ex.get('c_auto_grade2')}; "
  f"T1w (b_t1w) {ex.get('b_t1w')}; not registered / no QSM {ex.get('e_not_registered')} / {ex.get('no_qsm')} | SUP table1.exclusion_reason_component_counts (fail_reason__H_all180) |")
A(f"| age, y, mean (SD; range) | {mean_sd_rng(T['age'])} | SUP table1.age |")
A("| female/male | " + "; ".join(f"{T['female'][k]}/{T['male'][k]}" for k in CELLS) + " | SUP table1 |")
A("| race Black / White / Other | " + "; ".join(f"{T['race_black'][k]}/{T['race_white'][k]}/{T['race_other'][k]}" for k in CELLS) + " | SUP table1 (race3 coding of code_vars) |")
A(f"| threshold psychotic disorder | {c(T['threshold_psychotic_disorder'])} | SUP table1 (dx_psychosis) |")
A(f"| SN–VTA χ raw, mean (SD) | " + " / ".join(f"{r(x['mean'])} ({r(x['sd'])})" for x in T["chi_SNVTA_raw"]) + " | SUP table1.chi_SNVTA_raw |")
A(f"| SN–VTA χ M1-adjusted [95% CI] | " + " / ".join(f"{r(x['adj_mean'])} [{r(x['lo'])}, {r(x['hi'])}]" for x in T["chi_SNVTA_adjusted"]) + " | SUP table1.chi_SNVTA_adjusted (tier4 adj_cell_means) |")
for v in ("PRIME_total", "BDI_total", "ARI_total"):
    A(f"| {v.split('_')[0]} mean (SD; n) | {mean_sd_n(T[v])} | SUP table1.{v} |")
A(f"| Fisher, not in H1 by cell; earliest failed by cell | p = {T['qc_by_group_fisher_p_not_in_H1']:.4f}; p = {T['qc_by_group_fisher_p_earliest_failed']:.4f} | FR qc (frozen qc_by_group) |")
A(f"| age PS − non-PS (Welch) | {T['age_PS_minus_nonPS_years']:+.2f} y, p = {T['age_PS_welch_p']:.4f} | SUP table1 (tier4 4f t.test) |")
rz = {x["roi"]: x["n_flag_absz_gt5"] for x in FR["plausibility"]["robust_z"]}
A(f"| robust-z flags (abs z > 5, retained by rule) | chi_SNVTA {rz['chi_SNVTA']} (the TD youth sub-22045 ses-1, χ = −6.6 ppb), chi_DGM {rz['chi_DGM']}, chi_PUT_R {rz['chi_PUT_R']} | FR plausibility.robust_z |\n")

A("## 2. H1 (M1)\n")
A("| quantity | value | code path |\n|---|---|---|")
bf, bR = X["bound_frozen"], dec["B_registered"]
A(f"| n; df | {H['n']}; {dec['df']} | FR H1 |")
A(f"| β_PS [95% CI]; SE (HC3); p | {dec['beta_ppb']:+.4f} [{dec['ci95_ppb'][0]:.3f}, {dec['ci95_ppb'][1]:.3f}] ppb; {dec['se_hc3']:.3f}; {dec['p_two']:.4f} | FR H1.decision (run_m1_core / h1_decide) |")
A(f"| one-sided 95% lower limit | {X['one_sided_95_lower_limit']:+.3f} ppb | FR H1.decision.ci90_ppb[0] |")
A(f"| δ = β/σ_resid [bootstrap 95% CI]; noncentral-t CI; Hedges g | {dec['delta']:.3f} [{H['boot_delta']['ci95'][0]:.3f}, {H['boot_delta']['ci95'][1]:.3f}]; [{H['nct_ci_delta'][0]:.3f}, {H['nct_ci_delta'][1]:.3f}]; {H['hedges_g_raw']:.3f} | FR H1 |")
A(f"| σ_resid; σ_rob; σ_B | {dec['sigma_resid']:.3f}; {dec['sigma_rob']:.3f}; {dec['sigma_B']:.3f} | FR H1.decision |")
A(f"| Freedman–Lane p (10,000) [MC 95% CI] | {H['freedman_lane']['p']:.4f} [{H['freedman_lane']['mc_ci95'][0]:.4f}, {H['freedman_lane']['mc_ci95'][1]:.4f}] | FR H1.freedman_lane (fl_perm) |")
A(f"| Huber β (SE) | {H['huber']['est']:+.3f} ({H['huber']['se']:.3f}) | FR H1.huber |")
A(f"| bound test, frozen bound −6.1318 ppb: t_B; p_B | {bf['tB']:.3f}; {fp(bf['pB'])} | SUP H1_extras.bound_frozen (= kit2_m1_refit.json) |")
A(f"| bound test, recomputed: δ_B {bR['delta_B']:.3f} × σ_B → bound {bR['bound_ppb']:.3f} ppb; t_B; p_B | {bR['tB']:.3f}; {fp(bR['pB'])} | FR H1.decision.B_registered |")
A(f"| B co-reports (δ_B 0.6 / 0.5; cross-session r_L 0.54) | p {fp(dec['B_coreport'][0]['pB'])} / {fp(dec['B_coreport'][1]['pB'])}; {fp(H['B_coreport_cross_session'].get('pB'))} | FR H1 |")
A(f"| outcome (M1; after K3) | {dec['outcome']}; {H['final_outcome']} | FR H1 |")
k3 = H["K3"]; xc = X["K3_csf"]
A(f"| K3: β_PS(dref); p; trigger | {k3['beta_dref']:+.3f}; {k3['p_two']:.4f}; {k3['trigger']} | FR H1.K3 |")
A(f"| CSF-referenced β [95% CI]; SE; p; p_B (recomputed / frozen bound) | {xc['beta']:+.3f} [{xc['ci95'][0]:.3f}, {xc['ci95'][1]:.3f}]; {xc['se']:.3f}; {xc['p_two']:.4f}; {fp(k3['csf_pB'])} / {fp(xc['pB_frozen_bound'])} | FR H1.K3; SUP H1_extras.K3_csf |")
bb = H["bayes"]["BF"]
A(f"| BF01 attenuated N(−0.35, 0.25); BF01 sceptical N(0, 0.5) | {bb['attenuated']['BF01']:.3f}; {bb['sceptical']['BF01']:.3f} | FR H1.bayes (bayes_normal) |")
A(f"| residualised split-echo bound co-report | not available: the frozen code returns split-echo reliability r_L = {H['B_coreport_residualised_split']['detail']['rL']:.3f} (≤ 0); see the notes below | FR H1.B_coreport_residualised_split |\n")

A("## 3. K1 and K2\n")
A("| contrast | est [95% CI] (SE) | p | Holm p (2; post hoc) | p_B frozen bound | p_B recomputed bound | code path |\n|---|---|---|---|---|---|---|")
for key, lab in (("PSonly_minus_TD", "PS-only − TD"), ("PSADHD_minus_ADHDonly", "PS+ADHD − ADHD-only")):
    z = X["K1"][key]
    A(f"| {lab} | {z['est']:+.3f} [{z['ci95'][0]:.3f}, {z['ci95'][1]:.3f}] ({z['se']:.3f}) | {z['p_two']:.4f} | {z['p_holm_2']:.4f} | {fp(z['pB_frozen_bound'])} | {fp(z['pB_recomputed_bound'])} | SUP H1_extras.K1 (frozen k1_labels/lin_test) |")
for key, lab in (("ADHD_main", "ADHD main effect"), ("PSxADHD", "PS × ADHD interaction")):
    z = X["K1"][key]
    A(f"| {lab} | {z['est']:+.3f} [{z['ci95'][0]:.3f}, {z['ci95'][1]:.3f}] ({z['se']:.3f}) | {z['p_two']:.4f} | | | | FR H1.K1 |")
A(f"\nK1 label: {X['K1']['label']}. ADHD TOST ±{X['K1']['tost']['bound_ppb']:.3f} ppb: p {X['K1']['tost']['p_lower']:.4f} / {X['K1']['tost']['p_upper']:.4f}, so the TOST does not succeed.\n")
k2v = res["3_K1_K2"]["K2_R2s"]
A(f"K2 (R2\\*): β_PS = {k2v['est']:+.3f} s⁻¹ [95% CI {k2v['ci95'][0]:.3f}, {k2v['ci95'][1]:.3f}] (SE {k2v['se']:.3f}, n {k2v['n']}). One-sided p (lower) = {k2v['one_sided_lower_p']:.4f}; two-sided p = {k2v['p_two']:.4f}. Sources: FR H1.K2 and SUP fig3B_r2s.\n")

A("## 4. Sensitivity\n")
A("| row | n | β_PS (SE) | p | p_B frozen bound | recomputed bound; p_B | code path |\n|---|---|---|---|---|---|---|")
for s in sens_rows:
    A(f"| {s['label']} | {s['n']} | {s['beta_ppb']:+.3f} ({s['se_hc3']:.3f}) | {s['p_two']:.4f} | {fp(s['pB_frozen_bound'])} | {s['bound_recomputed_ppb']:.2f}; {fp(s['pB_recomputed_bound'])} | {s['source'].split(' (')[0]} |")
A(f"\nLeave-one-out: {loo['n_category_change']} of {loo['n_refits']} refits change the outcome (fragile: {loo['fragile']}); source FR H1.loo. "
  f"Dominance: triggered {H['influence']['diagnostics']['triggered']}.")
sc = FR["spec_curve"]
A(f"Specification curve: {sc['n_evaluable']} of {sc['n_specifications']} specifications are evaluable, because strict QC is now estimable. "
  f"The median δ is {sc['median_delta']:.3f} and the share with δ < 0 is {sc['share_delta_negative']:.3f}. "
  f"Robust: {sc['robust']} ({sc['robust_reason']}).")
A(f"The full registered battery (33 rows, with df and frozen-bound p_B) is in results.json → 4_sensitivity.full_battery. "
  f"Strict QC (grade 0 only, human grade) gives n {pick['strict QC (visual grade 0 only)']['n']} and β {pick['strict QC (visual grade 0 only)']['beta_ppb']:+.3f}. "
  f"Automated exclusions only gives n 90 and β +5.928, identical to the registered run because it does not depend on the rater.\n")

A("## 5. Secondary (Holm family: S1, S2)\n")
A("| test | n | est (SE) | p | Holm p | code path |\n|---|---|---|---|---|---|")
A(f"| S1 L | {S1['n']} | {S1['L']['est']:.3f} ({S1['L']['se']:.3f}) | {S1['L']['p_two']:.4f} | {FR['holm']['p_holm']['S1']:.4f} | FR S1.L (run_s1) |")
A(f"| S2 β_PRIME per SD | {S2['n']} | {S2['beta_PRIME']['est']:+.3f} ({S2['beta_PRIME']['se']:.3f}) | {S2['beta_PRIME']['p_two']:.4f} | {FR['holm']['p_holm']['S2']:.4f} | FR S2.beta_PRIME (run_s2) |\n")

A("## 6. Developmental check (PC2x, frozen pc2x(): chi ~ age, first-usable)\n")
A("| ROI | n | slope ppb/y [95% CI] | one-sided p | two-sided p |\n|---|---|---|---|---|")
for v in ("chi_PUT", "chi_GP"):
    z = SUP["PC2x"][v]
    A(f"| {v} | {z['n']} | {z['slope_ppb_per_y']:.3f} [{z['ci95'][0]:.3f}, {z['ci95'][1]:.3f}] | {fp(z['p_one_sided'])} | {fp(z['p_two_sided'])} |")
A(f"\nPC2x pass: {SUP['PC2x']['pass']}, because the putamen one-sided p is above .05. The caveat \"developmental positive control not confirmed\" therefore still applies. Sources: FR PC2x and SUP PC2x.\n")

A("## 7. Figure data\n")
A("Figure 3A shows the PS effect in raw-SD units (β/SD [95% CI]); source SUP fig3A_regional (tier4 ps_row on fit_y):\n")
A("| ROI | n | β ppb | p | β/SD [CI] |\n|---|---|---|---|---|")
for z in SUP["fig3A_regional"]:
    A(f"| {z['test']} | {z['n']} | {z['beta_PS']:+.3f} | {z['p_two']:.4f} | {z['beta_per_sd_raw']:+.3f} [{z['lo_sd']:+.3f}, {z['hi_sd']:+.3f}] |")
A("\nFigure 3B shows SN–VTA R2\\* by cell (s⁻¹). Raw values are mean (SD); M1-adjusted values are mean [95% CI]:\n")
A("| cell | raw | adjusted |\n|---|---|---|")
for x, y in zip(T["r2s_SNVTA_raw"], T["r2s_SNVTA_adjusted"]):
    A(f"| {x['cell']} (n {x['n']}) | {x['mean']:.2f} ({x['sd']:.2f}) | {y['adj_mean']:.2f} [{y['lo']:.2f}, {y['hi']:.2f}] |")
A("\nFigures (`figures/`):")
A("- **Figure 1**: χ by cell, n = 85.")
A("- **Figure 2**: forest of the section 4 rows plus primary, Huber and CSF. All rows use fitted HC3 SEs, so none is back-calculated. Bound −6.13 ppb; adult −14.24 ppb.")
A("- **Figure 3**: regional effects and R2\\*.")
A("- Assertion log: `figures/figs_humanqc_output.txt`.")
A("- Supplementary Figure S1 (TD centile chart) was not requested and was not redrawn.\n")

A("## 8. Second QSM engine (QSM.jl) on the H_all180 first-usable sessions\n")
E = SUP["engine2"]
A(f"The {E['n_first_usable']} H_all180 first-usable sessions were used. As in engine2_analysis.R:")
A("- one session has no engine-2 value (mask_mismatch);")
A(f"- under engine 2's plausibility rule (§12.4), {E['rts']['n_fail_plausibility']} more sessions fail and are set missing for engine 2 only.\n")
A("| inversion | analysis | n (cells) | β_PS [95% CI] | p | δ | bound (registered-style); p_B | strata est (p_B) |\n|---|---|---|---|---|---|---|---|")
for inv in ("rts", "tv"):
    for k, lab in (("main", "main (plausibility rule)"), ("sensitivity_failures_kept", "failures kept")):
        m = E[inv][k]; M = m["M1"]; cl = m["cells"]
        st = "; ".join(f"{v['est']:+.2f} ({fp(v['pB'])})" for v in m["strata"].values())
        A(f"| {inv.upper()} | {lab} | {m['n']} ({cl['TD']}/{cl['ADHD-only']}/{cl['PS-only']}/{cl['PS+ADHD']}) | {M['beta_ppb']:+.3f} [{M['ci95_ppb'][0]:.2f}, {M['ci95_ppb'][1]:.2f}] | {M['p_two']:.4f} | {M['delta']:.3f} | {m['bound_ppb_registered_deltaB']:.2f}; {fp(m['pB'])} | {st} |")
A("\nThe engine-2 bound is δ_B (0.581, borrowed from the primary engine's reliability) × engine 2's own σ_B. The result is a \"B pattern\", "
  "not a registered outcome. The frozen primary bound (−6.13 ppb) does not apply to engine 2's scale (its recovery k is lower). Source: SUP engine2.\n")

A("## Not computed, or computed differently, and why\n")
A("- **Residualised split-echo B co-report**: the frozen code gives split-echo reliability r_L = −0.23 in this sample, so it returns no bound. "
  "This is the co-report affected by coding deviation #52 (split-value plausibility rule ignored). The text-conforming patch was not "
  "re-run. The co-report is non-decisive and was not requested.")
A("- **Threshold-psychosis row label**: the frozen label says \"6 threshold-psychosis youths\", but 5 are in this sample (all PS+ADHD); "
  "in the n = 80 sample there were 3.")
A("- **Strict-QC and code-D rows**: strict QC now uses the human grade, because qc_grade was replaced. The dental/metal code D (`dental_D`) "
  "is still the AI rater's rubric code, because the human kits did not re-code it.")
A("- **S1 pooled constants**: these are the frozen label-free gate constants (gates JSON). They are not recomputed in this sample, "
  "as run_all.R specifies.")
A("- **T1, pooled and spec-curve**: these were also produced by the frozen run (see `frozen_run_H_all180/report.md`; T1 is not evaluated "
  "because A does not hold). A lme4 Nelder–Mead convergence warning in one of the four REML fits does not affect the selected bobyqa maximum.")
A("- **Supplementary Figure S1 and the recovery k-scaled adult test** were not redone.")
A("- **Engine 2 frozen bound**: no p_B at the frozen −6.1318 ppb bound is given, because that bound is on the primary engine's scale.\n")
A("## Files\n")
A("- `input/merged_H_all180.csv`: analysis table (sha256 in results.json)")
A("- `frozen_run_H_all180/`: frozen run_all.R output (report.md, results.json, manifest.json, sessionInfo.txt); log in `frozen_run_H_all180.log`")
A("- `supplement_H_all180.json`: supplement output")
A("- `results.json`: every number above, with sanity checks")
A("- `sanity_frozen_run_AI/` and `sanity_supplement_AI/`: the same pipeline on the frozen AI-QC table (reproduction checks)")
A("- `figures/`: Figure1–3 .tif and .jpg, plus the assertion log")
A("- `code/`: build_input.py, supplement_humanqc.R, figs_humanqc.py, make_results.py")
open(f"{BASE}/README.md", "w", encoding="utf-8").write("\n".join(L) + "\n")
print(json.dumps(checks, indent=1))
