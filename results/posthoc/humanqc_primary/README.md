# LEAD-χ on the complete-human-QC sample (H_all180): manuscript numbers

**POST HOC.** This is the complete-human-QC sample: the first author graded all 180 graded sessions (kit 1 + kit 2). It is not the registered analysis, which used the AI QC rater (n = 80) and gave outcome B. Nothing here re-decides H1. The rater was not independent of the study.

## How it was computed

1. `code/build_input.py` builds `input/merged_H_all180.csv` from the frozen merged table. It replaces only four columns, with the H_all180 values from `kit2_scenario_usable.csv`: `usable`, `first_usable`, `fail_reason` and `qc_grade` (the human grade; the AI grade is kept as `qc_grade_ai`). Check: rebuilding with the AI_frozen columns gives the frozen table byte for byte.
2. **`frozen_run_H_all180/`**: the frozen `mirror/scripts/analysis/run_all.R` was run **unchanged** on that table, with the frozen gates and multiverse files and registered draw counts. It writes `report.md` and `results.json`.
3. `code/supplement_humanqc.R` (output `supplement_H_all180.json`) computes what `run_all.R` does not print. It sources the frozen helpers read-only. The non-frozen pieces are listed in its header:
   - Table 1 descriptives;
   - stratum bound tests;
   - the joint nuisance row;
   - Figure 3 regional rows and adjusted means (helpers copied verbatim from `tier4_exploratory.R`);
   - the engine-2 refit (`engine2_analysis.R` logic).
4. `code/figs_humanqc.py` is a copy of `manuscript/code/figs_pcn.py` with the same plotting code. It reads only the files above and runs 26 assertions. It writes `figures/Figure{1,2,3}.tif` (600 dpi) and `.jpg` (300 dpi).
5. `code/make_results.py` assembles `results.json` and this README. It also re-runs the sanity checks.

Path abbreviations used below:
- **FR** = `frozen_run_H_all180/results.json` (frozen `R/m1.R`, `secondary.R`, `sensitivity.R`, `spec_curve.R`);
- **SUP** = `supplement_H_all180.json`.

## Sanity checks

| check | result |
|---|---|
| input_rebuild_AI_frozen_byte_identical_to_frozen_merged_table | True |
| frozen_run_all_on_AI_table_reproduces_results_confirmatory_json_exactly | {"ok": true, "n_differing_fields": 0} |
| AI_n80_beta4.9628 | True |
| AI_pipeline_reproduces_manuscript_Table1_and_Results_sample_text | {"ok": true, "mismatches": []} |
| AI_sensitivity_battery_with_df_wrapper_equals_frozen_battery_max_abs_diff | 0.0 |
| AI_engine2_reproduces_engine2_results_json | True |
| H_all180_reproduces_kit2_m1_refit | {"ok": true, "n": 85, "beta_PS": 6.28085215875773, "kit2_beta": 6.28085215875773} |
| supplement_M1_equals_frozen_run_M1_H_all180 | True |

Table 1 values are rounded half up, as in the manuscript (for example, TD PRIME 10.45 is shown as 10.5).

## 1. Sample (cells: TD / ADHD-only / PS-only / PS+ADHD)

| quantity | value | code path |
|---|---|---|
| analysed (first usable) | 23 / 18 / 17 / 27 (n = 85) | FR samples; SUP table1.analysed_first_usable (frozen build_samples) |
| usable sessions | 120 of 181 (37 / 23 / 22 / 38) | SUP table1.sessions_usable |
| youths with ≥ 1 MEGRE session | 30 / 22 / 24 / 37 | SUP table1.youths_with_megre |
| earliest MEGRE session not usable | 10 / 11 / 8 / 16 | frozen qc_by_group() table_earliest |
| analysed at ses-2 or ses-3 | 4 / 7 / 2 / 11 (ses-1 61, ses-2 20, ses-3 4) | SUP table1 |
| exclusion reasons (sessions; a session can have several) | registration/G5 (e_g5) 45; visual grade 2 or ungradable (c_qc_grade) 12; placement fail (d_placement) 12; (c') plausibility 8; automated grade 2 7; T1w (b_t1w) 6; not registered / no QSM 1 / 1 | SUP table1.exclusion_reason_component_counts (fail_reason__H_all180) |
| age, y, mean (SD; range) | 12.8 (2.3; 8.8–17.8) / 12.9 (2.3; 9.1–16.9) / 13.9 (1.7; 11.1–16.0) / 14.0 (1.7; 9.7–17.0) | SUP table1.age |
| female/male | 11/12; 10/8; 10/7; 12/15 | SUP table1 |
| race Black / White / Other | 7/11/5; 12/6/0; 11/4/2; 17/7/3 | SUP table1 (race3 coding of code_vars) |
| threshold psychotic disorder | 0 / 0 / 0 / 5 | SUP table1 (dx_psychosis) |
| SN–VTA χ raw, mean (SD) | 45.8 (19.1) / 50.9 (13.0) / 55.3 (11.7) / 54.8 (10.8) | SUP table1.chi_SNVTA_raw |
| SN–VTA χ M1-adjusted [95% CI] | 45.9 [38.2, 53.6] / 51.6 [44.7, 58.5] / 56.7 [50.1, 63.2] / 53.4 [48.9, 57.9] | SUP table1.chi_SNVTA_adjusted (tier4 adj_cell_means) |
| PRIME mean (SD; n) | 8.8 (9.6; n 20) / 17.6 (16.5; n 17) / 16.9 (14.6; n 17) / 23.0 (17.0; n 27) | SUP table1.PRIME_total |
| BDI mean (SD; n) | 5.3 (5.7; n 21) / 6.2 (6.4; n 17) / 3.5 (3.1; n 17) / 11.0 (11.1; n 26) | SUP table1.BDI_total |
| ARI mean (SD; n) | 2.9 (3.0; n 23) / 3.8 (2.6; n 18) / 3.5 (2.8; n 17) / 4.9 (4.0; n 27) | SUP table1.ARI_total |
| Fisher, not in H1 by cell; earliest failed by cell | p = 0.8298; p = 0.5645 | FR qc (frozen qc_by_group) |
| age PS − non-PS (Welch) | +1.18 y, p = 0.0084 | SUP table1 (tier4 4f t.test) |
| robust-z flags (abs z > 5, retained by rule) | chi_SNVTA 1 (the TD youth sub-22045 ses-1, χ = −6.6 ppb), chi_DGM 2, chi_PUT_R 1 | FR plausibility.robust_z |

## 2. H1 (M1)

| quantity | value | code path |
|---|---|---|
| n; df | 85; 73 | FR H1 |
| β_PS [95% CI]; SE (HC3); p | +6.2809 [-0.226, 12.787] ppb; 3.265; 0.0583 | FR H1.decision (run_m1_core / h1_decide) |
| one-sided 95% lower limit | +0.842 ppb | FR H1.decision.ci90_ppb[0] |
| δ = β/σ_resid [bootstrap 95% CI]; noncentral-t CI; Hedges g | 0.481 [0.062, 1.043]; [-0.017, 0.975]; 0.489 | FR H1 |
| σ_resid; σ_rob; σ_B | 13.066; 11.667; 11.667 | FR H1.decision |
| Freedman–Lane p (10,000) [MC 95% CI] | 0.0436 [0.0396, 0.0477] | FR H1.freedman_lane (fl_perm) |
| Huber β (SE) | +4.918 (2.874) | FR H1.huber |
| bound test, frozen bound −6.1318 ppb: t_B; p_B | 3.802; 0.0001 | SUP H1_extras.bound_frozen (= kit2_m1_refit.json) |
| bound test, recomputed: δ_B 0.581 × σ_B → bound -6.781 ppb; t_B; p_B | 4.001; 7.5e-05 | FR H1.decision.B_registered |
| B co-reports (δ_B 0.6 / 0.5; cross-session r_L 0.54) | p 5.9e-05 / 0.0002; 0.0004 | FR H1 |
| outcome (M1; after K3) | B; B | FR H1 |
| K3: β_PS(dref); p; trigger | +0.199; 0.9447; False | FR H1.K3 |
| CSF-referenced β [95% CI]; SE; p; p_B (recomputed / frozen bound) | +6.082 [-1.645, 13.810]; 3.877; 0.1211; 0.0002 / 0.0012 | FR H1.K3; SUP H1_extras.K3_csf |
| BF01 attenuated N(−0.35, 0.25); BF01 sceptical N(0, 0.5) | 3.518; 0.509 | FR H1.bayes (bayes_normal) |
| residualised split-echo bound co-report | not available: the frozen code returns split-echo reliability r_L = -0.226 (≤ 0); see the notes below | FR H1.B_coreport_residualised_split |

## 3. K1 and K2

| contrast | est [95% CI] (SE) | p | Holm p (2; post hoc) | p_B frozen bound | p_B recomputed bound | code path |
|---|---|---|---|---|---|---|
| PS-only − TD | +10.765 [0.420, 21.110] (5.191) | 0.0416 | 0.0832 | 0.0009 | 0.0006 | SUP H1_extras.K1 (frozen k1_labels/lin_test) |
| PS+ADHD − ADHD-only | +1.797 [-6.532, 10.126] (4.179) | 0.6685 | 0.6685 | 0.0309 | 0.0219 | SUP H1_extras.K1 (frozen k1_labels/lin_test) |
| ADHD main effect | +1.201 [-5.332, 7.733] (3.278) | 0.7151 | | | | FR H1.K1 |
| PS × ADHD interaction | -8.968 [-22.512, 4.576] (6.796) | 0.1911 | | | | FR H1.K1 |

K1 label: none. ADHD TOST ±5.833 ppb: p 0.0176 / 0.0809, so the TOST does not succeed.

K2 (R2\*): β_PS = +0.684 s⁻¹ [95% CI -0.395, 1.762] (SE 0.541, n 85). One-sided p (lower) = 0.8948; two-sided p = 0.2103. Sources: FR H1.K2 and SUP fig3B_r2s.

## 4. Sensitivity

| row | n | β_PS (SE) | p | p_B frozen bound | recomputed bound; p_B | code path |
|---|---|---|---|---|---|---|
| session-1-only first-usable | 61 | +6.117 (4.258) | 0.1572 | 0.0030 | -6.99; 0.0017 | supplement sensitivity.ses1_only |
| without threshold-psychosis youths (5 removed; frozen label says 6) | 80 | +5.829 (3.269) | 0.0790 | 0.0002 | -6.93; 0.0001 | supplement_H_all180.json sensitivity.battery_frozen_function |
| without prescribed-stimulant endorsers (4 removed) | 81 | +4.914 (3.049) | 0.1116 | 0.0003 | -6.84; 0.0001 | supplement_H_all180.json sensitivity.battery_frozen_function |
| + registration Dice | 85 | +4.710 (3.329) | 0.1614 | 0.0009 | -5.99; 0.0010 | supplement_H_all180.json sensitivity.battery_frozen_function |
| + motion index | 85 | +5.551 (3.315) | 0.0984 | 0.0004 | -6.78; 0.0002 | supplement_H_all180.json sensitivity.battery_frozen_function |
| + Dice + motion + CSF offset + ICV (joint; post hoc) | 85 | +4.410 (3.261) | 0.1807 | 0.0009 | -6.29; 0.0008 | supplement sensitivity.joint_nuisance |
| SN-only complex (descriptive) | 85 | +6.994 (3.452) | 0.0464 | 0.0001 | -7.46; 3.9e-05 | supplement_H_all180.json sensitivity.battery_frozen_function |
| VTA + PBP only (descriptive) | 85 | +3.161 (3.277) | 0.3378 | 0.0030 | -5.03; 0.0073 | supplement_H_all180.json sensitivity.battery_frozen_function |

Leave-one-out: 0 of 85 refits change the outcome (fragile: False); source FR H1.loo. Dominance: triggered False.
Specification curve: 486 of 486 specifications are evaluable, because strict QC is now estimable. The median δ is 0.341 and the share with δ < 0 is 0.109. Robust: False (final H1 outcome is not A).
The full registered battery (33 rows, with df and frozen-bound p_B) is in results.json → 4_sensitivity.full_battery. Strict QC (grade 0 only, human grade) gives n 82 and β +6.725. Automated exclusions only gives n 90 and β +5.928, identical to the registered run because it does not depend on the rater.

## 5. Secondary (Holm family: S1, S2)

| test | n | est (SE) | p | Holm p | code path |
|---|---|---|---|---|---|
| S1 L | 85 | 0.616 (0.439) | 0.1649 | 0.1649 | FR S1.L (run_s1) |
| S2 β_PRIME per SD | 70 | +4.613 (2.474) | 0.0673 | 0.1347 | FR S2.beta_PRIME (run_s2) |

## 6. Developmental check (PC2x, frozen pc2x(): chi ~ age, first-usable)

| ROI | n | slope ppb/y [95% CI] | one-sided p | two-sided p |
|---|---|---|---|---|
| chi_PUT | 85 | 0.639 [-0.148, 1.426] | 0.0550 | 0.1100 |
| chi_GP | 85 | 2.780 [1.374, 4.185] | 8.7e-05 | 0.0002 |

PC2x pass: False, because the putamen one-sided p is above .05. The caveat "developmental positive control not confirmed" therefore still applies. Sources: FR PC2x and SUP PC2x.

## 7. Figure data

Figure 3A shows the PS effect in raw-SD units (β/SD [95% CI]); source SUP fig3A_regional (tier4 ps_row on fit_y):

| ROI | n | β ppb | p | β/SD [CI] |
|---|---|---|---|---|
| SN-VTA (primary) | 85 | +6.281 | 0.0583 | +0.437 [-0.016, +0.889] |
| SN only (SNc+SNr) | 85 | +6.994 | 0.0464 | +0.451 [+0.007, +0.894] |
| VTA + PBP | 85 | +3.161 | 0.3378 | +0.241 [-0.256, +0.737] |
| globus pallidus | 85 | +2.111 | 0.5318 | +0.146 [-0.317, +0.608] |
| putamen | 85 | +0.342 | 0.8648 | +0.045 [-0.482, +0.572] |
| caudate | 85 | -3.021 | 0.2507 | -0.290 [-0.790, +0.209] |
| striatum (CAU+PUT+NAcc) | 85 | -0.693 | 0.6632 | -0.104 [-0.575, +0.368] |
| thalamus | 85 | -0.361 | 0.8258 | -0.053 [-0.529, +0.423] |
| red nucleus | 85 | +1.975 | 0.5392 | +0.161 [-0.359, +0.680] |
| subthalamic nucleus | 85 | +4.037 | 0.1652 | +0.307 [-0.129, +0.743] |
| right dentate | 59 | -0.311 | 0.9409 | -0.024 [-0.671, +0.623] |
| left dentate | 59 | +1.278 | 0.7520 | +0.095 [-0.507, +0.698] |
| deep GM composite (CAU+GP) | 85 | -1.840 | 0.4294 | -0.185 [-0.650, +0.279] |

Figure 3B shows SN–VTA R2\* by cell (s⁻¹). Raw values are mean (SD); M1-adjusted values are mean [95% CI]:

| cell | raw | adjusted |
|---|---|---|
| TD (n 23) | 22.88 (2.07) | 22.98 [22.11, 23.86] |
| ADHD-only (n 18) | 22.18 (2.10) | 22.50 [21.37, 23.63] |
| PS-only (n 17) | 24.56 (2.44) | 24.44 [23.27, 25.60] |
| PS+ADHD (n 27) | 22.63 (2.44) | 22.41 [21.41, 23.42] |

Figures (`figures/`):
- **Figure 1**: χ by cell, n = 85.
- **Figure 2**: forest of the section 4 rows plus primary, Huber and CSF. All rows use fitted HC3 SEs, so none is back-calculated. Bound −6.13 ppb; adult −14.24 ppb.
- **Figure 3**: regional effects and R2\*.
- Assertion log: `figures/figs_humanqc_output.txt`.
- Supplementary Figure S1 (TD centile chart) was not requested and was not redrawn.

## 8. Second QSM engine (QSM.jl) on the H_all180 first-usable sessions

The 85 H_all180 first-usable sessions were used. As in engine2_analysis.R:
- one session has no engine-2 value (mask_mismatch);
- under engine 2's plausibility rule (§12.4), 22 more sessions fail and are set missing for engine 2 only.

| inversion | analysis | n (cells) | β_PS [95% CI] | p | δ | bound (registered-style); p_B | strata est (p_B) |
|---|---|---|---|---|---|---|---|
| RTS | main (plausibility rule) | 62 (12/12/14/24) | +4.170 [-1.08, 9.42] | 0.1167 | 0.531 | -4.21; 0.0012 | +4.55 (0.0083); +3.79 (0.0164) |
| RTS | failures kept | 84 (23/18/17/26) | +4.665 [0.72, 8.61] | 0.0212 | 0.639 | -3.93; 2.3e-05 | +4.58 (0.0022); +4.75 (0.0009) |
| TV | main (plausibility rule) | 62 (12/12/14/24) | +5.001 [-0.71, 10.71] | 0.0849 | 0.573 | -4.76; 0.0006 | +5.80 (0.0055); +4.20 (0.0130) |
| TV | failures kept | 84 (23/18/17/26) | +5.186 [0.83, 9.55] | 0.0204 | 0.642 | -4.54; 1.5e-05 | +5.60 (0.0013); +4.77 (0.0011) |

The engine-2 bound is δ_B (0.581, borrowed from the primary engine's reliability) × engine 2's own σ_B. The result is a "B pattern", not a registered outcome. The frozen primary bound (−6.13 ppb) does not apply to engine 2's scale (its recovery k is lower). Source: SUP engine2.

## Not computed, or computed differently, and why

- **Residualised split-echo B co-report**: the frozen code gives split-echo reliability r_L = −0.23 in this sample, so it returns no bound. This is the co-report affected by coding deviation #52 (split-value plausibility rule ignored). The text-conforming patch was not re-run. The co-report is non-decisive and was not requested.
- **Threshold-psychosis row label**: the frozen label says "6 threshold-psychosis youths", but 5 are in this sample (all PS+ADHD); in the n = 80 sample there were 3.
- **Strict-QC and code-D rows**: strict QC now uses the human grade, because qc_grade was replaced. The dental/metal code D (`dental_D`) is still the AI rater's rubric code, because the human kits did not re-code it.
- **S1 pooled constants**: these are the frozen label-free gate constants (gates JSON). They are not recomputed in this sample, as run_all.R specifies.
- **T1, pooled and spec-curve**: these were also produced by the frozen run (see `frozen_run_H_all180/report.md`; T1 is not evaluated because A does not hold). A lme4 Nelder–Mead convergence warning in one of the four REML fits does not affect the selected bobyqa maximum.
- **Supplementary Figure S1 and the recovery k-scaled adult test** were not redone.
- **Engine 2 frozen bound**: no p_B at the frozen −6.1318 ppb bound is given, because that bound is on the primary engine's scale.

## Files

- `input/merged_H_all180.csv`: analysis table (sha256 in results.json)
- `frozen_run_H_all180/`: frozen run_all.R output (report.md, results.json, manifest.json, sessionInfo.txt); log in `frozen_run_H_all180.log`
- `supplement_H_all180.json`: supplement output
- `results.json`: every number above, with sanity checks
- `sanity_frozen_run_AI/` and `sanity_supplement_AI/`: the same pipeline on the frozen AI-QC table (reproduction checks)
- `figures/`: Figure1–3 .tif and .jpg, plus the assertion log
- `code/`: build_input.py, supplement_humanqc.R, figs_humanqc.py, make_results.py
