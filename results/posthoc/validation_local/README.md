# LEAD-chi: local validation analyses A-D (POST HOC; X2 registered-exploratory with a declared deviation)

> **DO NOT SEND. Internal only.** This folder holds unlocked results: the M1 results, every post hoc estimate, group labels by coded session, and the AI visual-QC grades of named coded sessions (section C, `C_low_TD_tail_sessions.csv`, `results.json`). The note covers every file here, including `CHANGES_2026-10-07.md`, `RECHECK_2026-10-07.md`, `verification/` and both `archive_*` folders, and any excerpt or summary of them.
> - **Blinded statistician:** must not receive any of it before the sealed comparison is complete. Their instructions ask them not to look at any LEAD-chi results, drafts or code.
> - **Human QC rater:** must not receive any of it before their ratings have been returned. The rater must stay blind to group, to the AI grades and to unlocked results.
> - At the re-check (`RECHECK_2026-10-07.md`, sec. 4), neither outgoing package (statistician package v1.1, QC rater kit) contained anything from this folder. The v1.1 statistician package has since been superseded without being sent. Check every new or rebuilt package the same way before it goes out.

**STATUS: POST HOC. Nothing here is a finding or a claim.** B, C and D are post hoc. A (X2) is a registered Tier-4 exploratory module (prereg Q25), **run with a declared deviation**: its constants were fixed after unlock, taken from the paper and not tuned (A.1, declaration 2). Everything here is hypothesis-generating. Nothing moves between tiers, and the confirmatory outcome is not re-decided.

The frozen confirmatory result is in `results_confirmatory/report.md` and `results.json`. It is not restated differently here: H1 final outcome **B**; beta_PS = +4.963 ppb (HC3 SE 3.125, p = .1169); outcome-B bound = -6.1318 ppb, t_B = 3.550, p_B = .0004; first-usable n = 80 (TD 22 / ADHD-only 16 / PS-only 17 / PS+ADHD 25).

**Revision 2026-10-07 (post hoc, after the independent verification in `verification/VERIFY.md`).** Every change is listed with its review item in `CHANGES_2026-10-07.md`. The pre-fix files are kept unchanged in `archive_v1_prefix_2026-10-07/`. The A-C CSVs are byte-identical to the pre-fix run. In D, the 90 mixed-model rows were refitted with the nlme default optimiser: lme slopes change by at most 2.2e-4 ppb/y, SEs by at most 0.014 % and p values by at most 0.0005; the first-usable OLS rows are unchanged. Ten printed numbers were corrected, and several sentences were reworded so that they claim no more than the analyses allow.

**Round 2 (2026-10-07, post hoc, wording only), after the independent re-check in `RECHECK_2026-10-07.md`.** The DO-NOT-SEND note above and the "Open decisions" section were added. A.1 declaration 5 now labels its p values as nominal and gives the q values of the matching referenced tests. The closing sentence of Reading (C) no longer calls all the low-chi sessions low-Dice. No analysis, CSV, `results.json` or code changed. The round-1 README is kept in `archive_v2_round1_2026-10-07/`, and every change is listed under "Round 2" in `CHANGES_2026-10-07.md`.

## Provenance

| item | value |
|---|---|
| analysis script | `code/validation_local.R` (R 4.4.1), SHA-256 `54ea89269a334a8293fb985c6dca9b62bcf29928c0915f1caab54f01c9ac8538`, run 2026-10-07T17:56:35Z (revised version; the pre-fix script, SHA-256 `280ae26f...4b8c35`, is in `archive_v1_prefix_2026-10-07/code/`) |
| verification scripts | `code/verify_statsmodels.py` (python3, statsmodels; HC3 algebra on exported designs). Independent from-raw re-implementation in Python: `verification/verify_independent.py` (review: `verification/VERIFY.md`; re-run after the fixes, see "Verification after the fixes" below) |
| input table | `mirror/unlocked/merged_qsmxt_primary_SNVTA.csv`, SHA-256 `f669bc6f41ab4e489abbd0a7c41949d3e5ccaba3e7f3c219d1bbebc3118df215` (the same file as the confirmatory input) |
| gates | `mirror/gates/gates_preunlock.json`, SHA-256 `dcf281b8...7157` (r_L = 0.938, PC4 TRUE; used only to reproduce the frozen B test) |
| frozen helpers (sourced read-only) | `mirror/scripts/analysis/R/{constants,utils,design,plausibility,m1}.R`. All five SHA-256 values equal the frozen `mirror/scripts/analysis/SHA256SUMS` |
| samples | validate -> code_vars -> apply_plausibility -> build_samples, identical to `run_all.R`. First-usable n = 80 (22/16/17/25); all usable sessions n = 107 from 80 subjects (26 subjects have 2 or more usable sessions) |
| model M1 | `y ~ PS_e * ADHD_e + ns(age, df = 3) + sex_e + race3 + scale(roi_vol) + scale(Q_phase)`, effect codes +/-0.5, race3 contr.sum, OLS with HC3, t reference with df = n - p (frozen `m1_formula` and `fit_ols_hc3`). `roi_vol` is the SN-VTA atlas volume in every M1-form model, including SN-only and GP, as in the frozen sensitivity rows |
| fixed outcome-B bound test | t_B = (beta + 6.1318) / SE_HC3, one-sided upper tail, df = n - p. The "frozen-style" p_B instead recomputes sigma_B (Huber) inside each model, as `sensitivity.R::sens_row` does. It is given only to reproduce the frozen numbers |
| other inputs (read-only) | `results_exploratory/tier4_exploratory/fdr_table.csv` (p values of the 5 registered Tier-4 tests; used only for the BH family sensitivity in A.2) |
| no server, no images | No server was contacted and no image was opened. Every value comes from the merged ROI table |

Every number below is in the CSV named in its section heading and in `results.json`. Rounding is for display only.

### Reproduction checks (all pass; `run.log`, `results.json -> checks`)

| quantity | this run | frozen |
|---|---|---|
| M1 beta_PS / SE_HC3 | 4.9628 / 3.1251 | 4.9628 / 3.1251 |
| bound (ppb) / t_B | -6.1318 / 3.5502 | -6.1318 / 3.550 |
| M1 + Dice / + motion / + ICV beta_PS | 3.7494 / 4.2268 / 5.3368 | 3.749 / 4.227 / 5.337 |
| frozen-style p_B, + Dice / + motion / + ICV | .0019 / .0004 / .0003 | .0019 / .0004 / .0003 |
| K2 R2* beta_PS (s^-1) | 1.0982 | 1.098 |
| PC2x age slope chi_PUT / chi_GP (ppb/y) | 0.3822 / 2.3966 | 0.382 / 2.397 |

Independent check (`verify/statsmodels_check.json`): statsmodels OLS with HC3 on the exported design matrices reproduces beta_PS and SE to < 1e-8 for M1, X2 SN-VTA chi_para and chi_dia, the joint nuisance model and the low-tail refit.

---

## A. X2 chi-separation (registered Tier-4 exploratory, run with a declared deviation; `A_X2_chisep_group_effects.csv`, `A_X2_component_descriptives.csv`, `A_X2_BH_family_sensitivity.csv`)

### A.1 What the paper specifies

Source: Dimov AV, Gillen KM, Nguyen TD, Kang J, Sharma R, Pitt D, Gauthier SA, Wang Y. *Tomography* 2022;8(3):1544-1551, doi:10.3390/tomography8030127 (open access; read from PMC9228115, because mdpi.com returned HTTP 403). Section 2.1, "MR Signal Modeling for R2*QSM" (the article spans pp. 1544-1551; the page of sec. 2.1 was not verified):

- **Eq. 1** models the complex-signal exponent as R2* + i2πf = r+|χ+| + r-|χ-| + iγ·d*(χ+ + χ-). The magnitude decay rate R2* is the sum of the two source magnitudes, each multiplied by a relaxometric constant. The field is the dipole convolution of their sum. Physical constraints: χ+ > 0 and χ- < 0.
- **Constant.** The paper states that r+ and r- "are relaxometric constants equal to 274 Hz/ppm". The value comes from the authors' earlier work (their ref. 22: Dimov et al. 2022, *J Neuroimaging*, doi:10.1111/jon.13014). The theoretical value of 321 Hz/ppm (their ref. 23: Yablonskiy and Haacke 1994) is cited for comparison. The acquisition was 3 T (sec. 2.2: ex vivo multiple-sclerosis tissue at 0.5 mm isotropic, 8 echoes, first TE and ΔTE 5.9 ms).
- **R2.** Eq. 1 has **no separate R2 (non-susceptibility) term**: R2* itself, not R2' = R2* - R2, carries the susceptibility sources. In the Discussion the authors say a constant component may be missing from Eq. 1. They argue it is small because ventricular CSF has a low R2*. They give no numeric R2. **So the assumed R2 is 0 s^-1.**
- **Not a closed form in the paper.** The paper solves Eq. 1 voxelwise as a regularised inverse problem (Eq. 2: weighted data terms plus edge-masked TV, conjugate gradient). The "closed form" used here is the unregularised, unconstrained point solution of Eq. 1 with r+ = r- = Dr. QSMxT calls this route "R2*-QSM closed-form (Dimov 2022)" (`prereg/server_capture/qsmxt_run_help.txt`, l. 795).

Only one short phrase of the paper is quoted above, to respect copyright. Everything else is paraphrased, with its section and equation number given.

**The model as implemented.** chi = chi_para + chi_dia; R2' = R2* - R2 = Dr·(chi_para - chi_dia). Therefore chi_para = (chi + R2'/Dr)/2 and chi_dia = (chi - R2'/Dr)/2. Units: chi in ppb, R2* in s^-1, Dr in Hz/ppm, so 1000·R2'/Dr is in ppb. Primary constants: **Dr = 274 Hz/ppm, R2 = 0**. Sensitivity: Dr x 0.5 (137), the theoretical 321, and Dr x 2 (548). The illustrative constants R2 = 10 and 15 s^-1 are arbitrary values, not taken from any paper. They are used only to show that R2 cancels. The relation is linear, so ROI-level components come from the ROI means `chi_<ROI>` (whole-brain referenced) and `r2s_<ROI>`. The ROIs are SN-VTA (primary), SN-only and GP (comparison).

**Declarations (model-driven separation).**
1. The data come from one field strength (Siemens Prisma 3 T), there is no R2 map, and Dr is assumed spatially constant. The separation is therefore **model-driven, not measured**.
2. **Deviation (declared).** Q25 required Dr and R2 to be "written into the amendment-2 module before any χ exists". `prereg/LEADchi_amendment_2_FROZEN.md` contains no X2 constant (searched for X2, Dimov, separation, relaxometric, 274 and Hz/ppm). The constants were therefore fixed now, after unlock. They are taken directly from the paper, were not tuned, and the full Dr range is reported.
3. Because the design matrix is shared, X2 adds **no information beyond chi and R2* jointly**: beta(chi_para) = [beta(chi) + 1000·beta(R2*)/Dr]/2 exactly. The numerical deviation is below 1e-13 in every row (floating-point noise).
4. **R2 cancels in group contrasts.** A constant R2 shifts every session's chi_para by -500·R2/Dr ppb and its chi_dia by +500·R2/Dr ppb. Over R2 = 0, 10 and 15 s^-1, beta_PS and SE_HC3 differ by less than 1e-13 (floating-point noise). Only the absolute partition moves. At R2 = 10 s^-1, 61 of 80 SN-VTA sessions (76%) would have chi_dia > 0, which violates the χ- < 0 constraint.
5. **The whole-brain chi reference shifts both components equally** (by ref_ppb/2). Adding the reference back (chi + ref_ppb) changes nothing material: SN-VTA chi_para 4.485 [0.214, 8.757], p .0399; chi_dia 0.477 [-2.373, 3.328], p .7394; SN-only 5.708 / 0.166; GP 1.575 / 0.699. These p values are nominal: the unreferenced rows are outside every BH family. The matching referenced tests have q .118 (chi_para) and .889 (chi_dia) in the family of 6 (A.2). The separation cannot repair or reveal a referencing error.
6. The chi and R2* ROI means do not use identical voxel sets. Chi uses support S without implausible voxels. R2* uses Mask 3 with voxels clipped to [0, 200] s^-1 (`extract_rois.py`). The ROI-level formula therefore approximates the ROI mean of voxelwise closed-form maps. Positivity constraints are not enforced. At Dr = 274 they are met by every SN-VTA and SN-only session and by 78 of 80 GP sessions.
7. **R2* voxel set: departure from the plan text (labelled here, 2026-10-07, post hoc).** The preregistration (`prereg/LEADchi_prereg_v1.0_FINAL.md`, line 253, sec. 12.2) defines R2*_ROI on the "same voxel sets as the matching χ ROI, restricted to Mask 3". The hashed extraction code (`extract_rois.py`, lines 289-295; SHA-256 `23592ebb...a693`, as frozen in amendment 1) instead averages R2* over Mask 3 ∩ finite R2* inside the ROI, without the chi support S and without the implausible-voxel exclusion (declaration 6). This is a departure from the registered text. It affects every R2* input used here (X2 in A, B.2 and D) and also the frozen K2 and PC3. It is not listed in Supplementary Table S1 of manuscript draft v3 or in the amendment-1 log (searched 2026-10-07). Its numerical size cannot be computed from the ROI table, because that needs the voxel maps; it was not computed. Whether to add it to the deviation annex (Annex C.4) is for the authors.

### A.2 Group effects at the paper constants (Dr = 274 Hz/ppm, R2 = 0; first-usable n = 80; M1)

| ROI | outcome | beta_PS ppb [95% CI] | SE | p | BH q (X2 family of 6) | PS-only - TD [95% CI], p (BH q, family of 18) | PS+ADHD - ADHD-only [95% CI], p (BH q, family of 18) |
|---|---|---|---|---|---|---|---|
| SN-VTA | total chi (input) | 4.96 [-1.27, 11.20] | 3.13 | .1169 | - | 8.16 [0.21, 16.10], .0443 | 1.77 [-7.78, 11.32], .7129 |
| SN-VTA | R2* (input; s^-1) | 1.10 [0.08, 2.12] | 0.51 | .0350 | - | 1.74 [0.13, 3.35], .0344 | 0.46 [-0.94, 1.86], .5154 |
| SN-VTA | chi_para | 4.49 [0.23, 8.74] | 2.13 | .0392 | .118 | 7.25 [1.04, 13.46], .0228 (q .137) | 1.72 [-4.49, 7.94], .5825 (q .986) |
| SN-VTA | chi_dia | 0.48 [-2.39, 3.35] | 1.44 | .7408 | .889 | 0.91 [-2.28, 4.09], .5720 (q .986) | 0.05 [-4.42, 4.52], .9831 (q .998) |
| SN-only | total chi (input) | 5.87 [-0.73, 12.48] | 3.31 | .0804 | - | 9.16 [0.60, 17.71], .0362 | 2.59 [-7.63, 12.81], .6148 |
| SN-only | R2* (input; s^-1) | 1.52 [0.38, 2.65] | 0.57 | .0096 | - | 2.19 [0.35, 4.03], .0206 | 0.85 [-0.71, 2.42], .2819 |
| SN-only | chi_para | 5.71 [1.09, 10.32] | 2.31 | .0161 | .096 | 8.57 [1.59, 15.55], .0169 (q .137) | 2.85 [-3.87, 9.56], .4005 (q .986) |
| SN-only | chi_dia | 0.17 [-2.85, 3.19] | 1.51 | .9128 | .913 | 0.59 [-2.64, 3.82], .7162 (q .998) | -0.26 [-5.10, 4.58], .9158 (q .998) |
| GP | total chi (input) | 2.27 [-3.59, 8.13] | 2.94 | .4413 | - | 3.45 [-5.19, 12.08], .4285 | 1.10 [-6.95, 9.16], .7855 |
| GP | R2* (input; s^-1) | 0.24 [-1.09, 1.57] | 0.66 | .7191 | - | 0.79 [-1.10, 2.67], .4079 | -0.31 [-2.15, 1.54], .7414 |
| GP | chi_para | 1.58 [-3.09, 6.24] | 2.34 | .5028 | .889 | 3.16 [-3.79, 10.11], .3679 (q .986) | -0.01 [-6.44, 6.42], .9983 (q .998) |
| GP | chi_dia | 0.70 [-1.97, 3.37] | 1.34 | .6026 | .889 | 0.29 [-3.26, 3.84], .8715 (q .998) | 1.11 [-2.59, 4.81], .5519 (q .986) |

The total chi and R2* input rows are not X2 tests and are not in any BH family; their p values are nominal.

**BH-FDR within the X2 block.** The primary family has 6 tests: 3 ROIs x {chi_para, chi_dia} beta_PS, at Dr = 274 and R2 = 0. The smallest q is .096 (SN-only chi_para); SN-VTA chi_para has q .118. An extended family of 18 tests (the same 6 beta_PS plus their 12 simple effects) gives a smallest q of .137; in it, SN-VTA chi_para beta_PS has q .177. **No X2 test survives BH-FDR.**

**The BH family is a choice (post hoc note).** Q25 says "BH-FDR within block" but does not define the block. The X2-only families above were chosen after unlock. If the block is the whole registered Tier-4 set as run, the 6 X2 tests join the 5 registered tests in `results_exploratory/tier4_exploratory/fdr_table.csv` (X1 right dentate, left dentate and red nucleus, X3, X6; all p ≥ .576). That family of 11 gives a smallest q of .177, and the corresponding family of 23 (18 + 5) gives .175 (`A_X2_BH_family_sensitivity.csv`). No X2 test survives under any of these definitions. X4, X5 and X7 are not in that table, so they could not be included.

Cell means at Dr = 274 (ppb; TD / ADHD-only / PS-only / PS+ADHD):
- SN-VTA chi_para: 66.50 / 65.20 / 73.42 / 68.45
- SN-VTA chi_dia: -16.41 / -14.34 / -16.57 / -14.18
- SN-only chi_para: 71.63 / 70.19 / 79.86 / 74.79
- SN-only chi_dia: -13.90 / -11.68 / -14.30 / -11.75
- GP chi_para: 100.54 / 99.58 / 106.71 / 102.02
- GP chi_dia: -10.30 / -9.09 / -9.62 / -8.30

### A.3 Sensitivity to Dr (R2 = 0)

| ROI | Dr (Hz/ppm) | chi_para beta [95% CI], p | chi_dia beta [95% CI], p | chi_para PS-only - TD, p | chi_dia PS-only - TD, p | sessions with chi_dia > 0 |
|---|---|---|---|---|---|---|
| SN-VTA | 137 | 6.49 [0.71, 12.27], .0284 | -1.53 [-5.22, 2.17], .4123 | 10.42, .0219 | -2.26, .3358 | 0/80 |
| SN-VTA | **274** | 4.49 [0.23, 8.74], .0392 | 0.48 [-2.39, 3.35], .7408 | 7.25, .0228 | 0.91, .5720 | 0/80 |
| SN-VTA | 321 | 4.19 [0.13, 8.25], .0431 | 0.77 [-2.06, 3.60], .5887 | 6.79, .0236 | 1.37, .3876 | 4/80 |
| SN-VTA | 548 | 3.48 [-0.13, 7.10], .0586 | 1.48 [-1.37, 4.33], .3036 | 5.66, .0271 | 2.49, .1350 | 68/80 |
| SN-only | 137 | 8.48 [2.13, 14.82], .0096 | -2.60 [-6.60, 1.39], .1975 | 12.56, .0154 | -3.40, .1798 | 0/80 |
| SN-only | **274** | 5.71 [1.09, 10.32], .0161 | 0.17 [-2.85, 3.19], .9128 | 8.57, .0169 | 0.59, .7162 | 0/80 |
| SN-only | 321 | 5.30 [0.92, 9.69], .0185 | 0.57 [-2.40, 3.55], .7024 | 7.98, .0176 | 1.17, .4646 | 7/80 |
| SN-only | 548 | 4.32 [0.45, 8.20], .0294 | 1.55 [-1.44, 4.54], .3041 | 6.57, .0207 | 2.58, .1334 | 75/80 |
| GP | 137 | 2.01 [-4.82, 8.85], .5585 | 0.26 [-3.91, 4.43], .9009 | 4.59, .3665 | -1.15, .6786 | 0/80 |
| GP | **274** | 1.58 [-3.09, 6.24], .5028 | 0.70 [-1.97, 3.37], .6026 | 3.16, .3679 | 0.29, .8715 | 2/80 |
| GP | 321 | 1.51 [-2.86, 5.88], .4927 | 0.76 [-1.80, 3.33], .5551 | 2.95, .3701 | 0.50, .7748 | 31/80 |
| GP | 548 | 1.36 [-2.35, 5.06], .4674 | 0.92 [-1.61, 3.44], .4708 | 2.44, .3812 | 1.01, .5745 | 80/80 |

The p values in this table are nominal: the Dr-sensitivity rows are outside every BH family (the q values at Dr = 274 are in A.2). The ordering of the point estimates (chi_para > chi_dia) is the same at all four Dr values in SN-VTA and SN-only. The SN-VTA chi_para CI excludes 0 only for Dr ≤ 321; at Dr = 548 it includes 0 (3.48 [-0.13, 7.10], p .0586).

Component values (mean (SD), ppb; first-usable). Mean R2* is 22.92 s^-1 for SN-VTA, 23.82 for SN-only and 30.53 for GP.
- At Dr = 274: SN-VTA chi_para 68.32 (9.13), chi_dia -15.33 (5.14); SN-only 74.08 (10.52) and -12.87 (5.37); GP 102.12 (10.77) and -9.29 (4.81).
- At Dr = 548, chi_dia is positive in most sessions (SN-VTA 68/80, GP 80/80), so that Dr is physically implausible for these data.
- At Dr = 137, chi_para and chi_dia are negatively correlated (SN-VTA r = -0.44). At 274 the correlation is +0.33.

### A.4 Reading (X2)

Under the model's assumptions, the point estimates fall mainly in chi_para: SN-VTA chi_para +4.49 ppb [0.23, 8.74] (p .039, q .118) and chi_dia +0.48 [-2.39, 3.35] (p .74, q .889). SN-only gives +5.71 (p .016, q .096) and +0.17 (p .91, q .913). No component test survives BH-FDR (smallest q .096), so none of these is a finding. There is no evidence of a GP difference in either component (chi_para +1.58 [-3.09, 6.24], chi_dia +0.70 [-1.97, 3.37]); these CIs are wide, so this is not evidence of absence. The ordering of the point estimates (chi_para > chi_dia) is the same at all four Dr values, but the SN-VTA chi_para CI excludes 0 only for Dr ≤ 321 (A.3). X2 is a deterministic restatement of beta(chi) and beta(R2*), so this is the X2 form of "chi and R2* move together" (tier-4 section 4d), not new evidence. No iron model was compared with a myelin model. As in the registered K2 wording, the iron-versus-myelin source remains unresolved. The positive direction remains unclaimed, and outcome B stands.

---

## B. Joint nuisance model and Dice analyses (POST HOC; `B_nuisance_models.csv`, `B_joint_model_coefficients.csv`, `B_r2s_models.csv`, `B_dice_tertile_*.csv`)

### B.1 Nuisance covariates, alone and together (first-usable n = 80; chi_SNVTA)

| model | n | beta_PS ppb [95% CI] | SE | t (df) | p | PS-only - TD, p | PS+ADHD - ADHD-only, p | t_B (bound -6.1318) | p_B | p_B frozen-style |
|---|---|---|---|---|---|---|---|---|---|---|
| M1 (frozen primary) | 80 | 4.963 [-1.273, 11.199] | 3.125 | 1.588 (68) | .1169 | 8.16, .0443 | 1.77, .7129 | 3.550 | .0004 | .0004 |
| M1 + Dice (alone) | 80 | 3.749 [-2.450, 9.949] | 3.106 | 1.207 (67) | .2316 | 5.15, .1709 | 2.35, .6293 | 3.181 | .0011 | .0019 |
| M1 + motion (alone) | 80 | 4.227 [-1.717, 10.171] | 2.978 | 1.419 (67) | .1604 | 7.53, .0632 | 0.92, .8375 | 3.478 | .0004 | .0004 |
| M1 + CSF offset (ref2_minus_wb_ppb; alone) | 80 | 5.461 [-0.855, 11.776] | 3.164 | 1.726 (67) | .0890 | 8.76, .0310 | 2.16, .6527 | 3.664 | .0002 | .0002 |
| M1 + ICV (alone) | 80 | 5.337 [-1.156, 11.830] | 3.253 | 1.641 (67) | .1056 | 9.05, .0321 | 1.63, .7287 | 3.526 | .0004 | .0003 |
| **M1 + Dice + motion + CSF offset + ICV (joint)** | 80 | **4.162 [-1.900, 10.224]** | 3.034 | 1.372 (64) | .1750 | 6.32, .1234 | 2.00, .6421 | **3.392** | **.0006** | .0007 |

Nuisance coefficients in the joint model (ppb per SD of the covariate):

| term | estimate [95% CI] | SE | p |
|---|---|---|---|
| Dice (dice_min) | **+4.29 [1.65, 6.93]** | 1.32 | .0019 |
| motion (motion_gsr) | -1.50 [-6.04, 3.04] | 2.27 | .5114 |
| CSF offset (ref2_minus_wb_ppb) | -1.40 [-4.50, 1.71] | 1.55 | .3713 |
| ICV (icv_seg_ml) | +4.19 [-1.38, 9.77] | 2.79 | .1377 |
| ADHD_e | +0.41 [-5.20, 6.01] | 2.81 | .8852 |

The variance inflation of PS_e is 1.18 in M1 and 1.27 in the joint model, so collinearity is negligible. Dice and chi_SNVTA correlate with Spearman rho = .358 (p .0011). Note: the CSF offset is itself derived from the chi map, so it is not a pure nuisance variable.

### B.2 SN-VTA R2* with and without Dice

| model | beta_PS s^-1 [95% CI] | SE | p | PS-only - TD, p | PS+ADHD - ADHD-only, p |
|---|---|---|---|---|---|
| R2* M1 (= frozen K2 estimate) | 1.098 [0.080, 2.117] | 0.510 | .0350 | 1.738, .0344 | 0.458, .5154 |
| R2* M1 + Dice | 1.011 [-0.079, 2.101] | 0.546 | .0684 | 1.522, .0733 | 0.500, .4937 |
| R2* M1 + joint nuisance | 1.026 [-0.098, 2.151] | 0.563 | .0730 | 1.598, .0627 | 0.454, .5422 |

Dice is not related to R2* (+0.272 s^-1 per SD, p .3624). The R2* estimate moves little (1.10 to 1.01, 8 % lower) and its SE grows by 7 %; together these move p across .05 (.0350 to .0684). The frozen K2 test is one-sided for *lower* R2* (p .9825) and is not re-decided here.

### B.3 Dice tertiles (label-free cut points 0.8054 and 0.8233 over the 80 first-usable sessions)

| tertile | Dice range | n | TD / ADHD-only / PS-only / PS+ADHD | mean chi_SNVTA (ppb) |
|---|---|---|---|---|
| T1 low | 0.7371-0.8052 | 27 | 11 / 5 / 2 / 9 | 50.29 |
| T2 mid | 0.8058-0.8227 | 26 | 5 / 6 / 8 / 7 | 50.04 |
| T3 high | 0.8236-0.8667 | 27 | 6 / 5 / 7 / 9 | 58.51 |

Tertile by cell: Fisher p = .32. Tertile by PS: Fisher p = .36.

| model | n (cells) | beta_PS ppb [95% CI] | SE | df | p | PS-only - TD, p | PS+ADHD - ADHD-only, p | p_B (bound -6.1318) |
|---|---|---|---|---|---|---|---|---|
| within T1 low: full M1 | 27 (11/5/2/9) | -0.95 [-14.60, 12.70] | 6.40 | 15 | .8836 | 5.06, .6043 | -6.97, .4167 | .2157 |
| within T1 low: reduced (PS*ADHD + age + sex) | 27 | -3.31 [-14.29, 7.66] | 5.28 | 21 | .5371 | 1.10, .8987 | -7.72, .2770 | .2993 |
| within T2 mid: full M1 | 26 (5/6/8/7) | 14.22 [-6.71, 35.16] | 9.76 | 14 | .1671 | 19.59, .1230 | 8.86, .4965 | .0279 |
| within T2 mid: reduced | 26 | 13.20 [-3.16, 29.57] | 7.84 | 20 | .1080 | 14.89, .0421 | 11.51, .3624 | .0115 |
| within T3 high: full M1 | 27 (6/5/7/9) | 7.37 [-5.37, 20.10] | 5.97 | 15 | .2366 | 8.95, .1613 | 5.79, .5032 | .0196 |
| within T3 high: reduced | 27 | 2.23 [-10.70, 15.15] | 6.21 | 21 | .7236 | 2.27, .7323 | 2.18, .7964 | .0965 |
| pooled M1 + Dice-tertile fixed effects | 80 | 4.05 [-2.05, 10.15] | 3.06 | 66 | .1896 | 7.13, .0619 | 0.98, .8329 | .0007 |
| pooled M1 + tertile FE + linear Dice | 80 | 4.10 [-1.97, 10.16] | 3.04 | 65 | .1822 | 5.85, .1260 | 2.34, .6102 | .0006 |
| quality-matched: drop lowest tertile (M1) | 53 (11/11/15/16) | 10.36 [0.40, 20.32] | 4.93 | 41 | .0419 | 11.76, .0292 | 8.96, .2324 | .0009 |
| quality-matched: drop lowest tertile (M1 + Dice) | 53 | 9.17 [-0.07, 18.42] | 4.57 | 40 | .0518 | 10.50, .0468 | 7.84, .2554 | .0009 |

Heterogeneity of beta_PS across tertiles (M1 + tertile + PS_e x tertile; HC3 Wald): F(2, 64) = 2.44, p = .095. With 26-27 sessions per tertile and 12 M1 parameters, the within-tertile estimates are very imprecise.

### B.4 Reading (B)

Adding all four nuisance covariates together lowers beta_PS from 4.96 to 4.16 ppb (p .175). The outcome-B bound is still excluded (p_B .0006). Dice is the only nuisance covariate related to chi (+4.3 ppb per SD), and adding it lowers the estimate by about 1 ppb, alone or jointly. The lowest-Dice tertile holds half of the TD sessions (11/22) but only 2 of the 17 PS-only sessions. This imbalance is not significant (Fisher p .32) and is consistent with the attenuation under Dice adjustment. The point estimate does not shrink in the upper two tertiles (+10.4, p .042, n 53; a post hoc subsample), and in the lowest tertile it is about zero (-0.95). Heterogeneity across tertiles is not significant (p .095), and all three within-tertile CIs include 0, so these subsamples cannot test the registration-quality explanation. Every pooled specification keeps beta_PS between +3.7 and +5.5 ppb with the B bound excluded. This holds only for the pooled models: within single tertiles the bound is not always excluded (lowest tertile p_B .2157 full M1 and .2993 reduced; reduced high-tertile model .0965). Nothing here changes outcome B, and nothing licenses a positive-direction claim.

---

## C. Low-TD-tail check (POST HOC; `C_low_TD_tail_sessions.csv`, `C_low_TD_tail_qc_comparison.csv`, `C_low_TD_tail_refits.csv`)

Four first-usable TD sessions have chi_SNVTA < 35 ppb. All are each subject's first session (ses-1). Percentiles are within the other 18 first-usable TD sessions. No image was opened.

| code | session | chi_SNVTA | qc_grade | dice_min (TD pct) | motion_gsr (TD pct) | fimp_SNVTA | inmask_SNVTA | Q_phase | r2s_SNVTA (TD pct) | M1 studentised residual |
|---|---|---|---|---|---|---|---|---|---|---|
| 2c7b34 | ses-1 | 29.02 | 1 | 0.7589 (0) | 0.0062 (89) | 0 | 1.0000 | -2.869 | 22.95 (50) | -1.96 |
| a38f14 | ses-1 | 29.72 | 1 | 0.7751 (0) | 0.0080 (89) | 0 | 0.9962 | -3.184 | 19.02 (0) | -1.02 |
| 10a248 | ses-1 | 30.98 | 1 | 0.8143 (50) | -0.0010 (28) | 0 | 1.0000 | -3.046 | 19.81 (6) | -1.86 |
| a516f3 | ses-1 | 33.37 | 1 | 0.7371 (0; lowest of all 80) | 0.0129 (100) | 0 | 1.0000 | -2.842 | 22.69 (39) | -1.03 |

Additional per-session values are in the CSV: roi_vol, chi_SNVTA_median, chi_SNonly, chi_VTAPBP, CSF offset, age, and percentiles within all 80.

| metric | tail mean (min-max) | other TD median [IQR] (min-max), n = 18 | all-80 median | Wilcoxon p (tail vs other TD) |
|---|---|---|---|---|
| dice_min | 0.7713 (0.7371-0.8143) | 0.8107 [0.7990, 0.8268] (0.7759-0.8349) | 0.8137 | **.0194** |
| motion_gsr | 0.0065 (-0.0010-0.0129) | 0.0028 [-0.0009, 0.0042] (-0.0036-0.0123) | 0.0026 | .1178 |
| fimp_SNVTA | 0 (all) | 0 (all) | 0 | NA |
| inmask_SNVTA | 0.9991 (0.9962-1.0000) | 1.0000 [1.0000, 1.0000] (0.9690-1.0000) | 1.0000 | .3061 |
| Q_phase | -2.99 (-3.18 to -2.84) | -3.06 [-3.15, -2.94] (-3.26 to -2.57) | -3.08 | .5936 |
| roi_vol (mm^3) | 627.1 (551.0-713.0) | 629.3 [595.3, 691.2] (557.3-753.6) | 624.3 | .5936 |
| r2s_SNVTA (s^-1) | 21.12 (19.02-22.95) | 22.93 [21.79, 24.49] (19.03-26.17) | 22.77 | .1178 |
| CSF offset (ppb) | 14.34 (8.94-18.85) | 19.87 [15.16, 22.98] (6.60-38.13) | 21.59 | .1408 |
| ICV (ml) | 1552 (1317-1730) | 1583 [1535, 1660] (1412-1742) | 1525 | .9023 |
| age (y) | 12.48 (8.80-15.36) | 13.05 [12.00, 14.46] (10.03-17.80) | 13.62 | .7743 |

QC grade: all 4 tail sessions are grade 1. The other 18 TD sessions are 2 grade 0 and 16 grade 1; all 80 are 4 grade 0 and 76 grade 1.

| refit | n (cells) | beta_PS ppb [95% CI] | SE | p | PS-only - TD, p | PS+ADHD - ADHD-only, p | t_B | p_B | Huber beta |
|---|---|---|---|---|---|---|---|---|---|
| M1, all 80 | 80 (22/16/17/25) | 4.96 [-1.27, 11.20] | 3.13 | .1169 | 8.16, .0443 | 1.77, .7129 | 3.550 | .0004 | 3.94 |
| M1 without the 4 low-tail TD sessions | 76 (18/16/17/25) | **3.33 [-3.09, 9.74]** | 3.21 | .3040 | 4.36, .2763 | 2.30, .6393 | 2.946 | .0022 | 2.32 |
| M1 without every session < 35 ppb (label-free; 4 TD + 1 ADHD-only at 12.16 ppb, the lowest of all 80) | 75 (18/15/17/25) | 1.45 [-3.63, 6.53] | 2.54 | .5703 | 4.04, .2847 | -1.14, .7468 | 2.982 | .0020 | 1.47 |

**Reading (C).** None of the four sessions fails any automated rule: implausible-voxel share 0, in-support share ≥ 0.996, and all are QC grade 1 like most sessions. None is a model outlier (|studentised residual| ≤ 1.96). Three of the four share poor registration: they have a lower Dice than every other TD session, and one has the lowest Dice of all 80. The fourth (10a248) is at the TD median (50th percentile). Their mean Dice is 0.771 versus 0.811 (p .019). Three of the four also show high motion (89th-100th TD percentile), and two have low SN-VTA R2* (0th and 6th TD percentile). This would be consistent with partial-volume dilution of a bright nucleus, but that mechanism was not examined (no image was opened). Removing them lowers beta_PS from 4.96 to 3.33 ppb, and PS-only - TD from 8.16 (p .044) to 4.36 (p .276). The label-free cut lowers beta_PS to 1.45. The bound test stays excluded in every variant (p_B ≤ .0022). In short, about a third of the positive point estimate (about 70 % with the label-free cut) rests on a handful of low-chi sessions, and the nominal PS-only - TD contrast depends on them. Not all of them are low-Dice sessions: three have the three lowest TD Dice values, but 10a248 (0.814) and the ADHD-only session removed by the label-free cut (0.810) are close to the median Dice of all 80 sessions (0.814). These refits select on the outcome, so they are diagnostic bounds, not corrected estimates. Whether these sessions should be down-weighted is a question for the blinded human QC (the requested co-author re-rating); these values must not decide it.

---

## D. Putamen / GP positive-control diagnostics (POST HOC; `D_age_slopes.csv`, `D_implied_chi_slope_from_R2s.csv`, `D_putamen_shortfall_by_Dr.csv`)

Models:
- First-usable sample: OLS with HC3 (n = 80; age 8.8-17.8 y).
- All usable sessions: `nlme::lme`, REML with the nlme default optimiser (nlminb; revised 2026-10-07, see CHANGES), random intercept by subject (n = 107, 80 subjects, age 8.8-19.5 y). The df for age is nlme containment: 26 (107 - 80 - 1), which is conservative.
- Covariates: age only; + sex_e; + sex_e + PS_e * ADHD_e (cell). A within/between-person split (subject-mean age vs deviation; + sex + cell) is co-reported.
- Registered expectation (prereg PC2w: 2.6-3.4 ppb putamen and 4.4-6.9 ppb GP over 1.95 y): putamen chi **1.33-1.74 ppb/y**; GP chi **2.26-3.54 ppb/y**.
- `chi_PUT_raw` and `chi_GP_raw` (the referenced ROI means before the ROI-fail rule) are **identical** to `chi_PUT` and `chi_GP` in every usable session (no putamen or GP fails), so they give the same slopes.

### D.1 Primary outcomes

| outcome | sample | covariates | slope per y [95% CI] | p (two-sided) | slope in SD/y | CI excludes registered range |
|---|---|---|---|---|---|---|
| **r2s_PUT** (s^-1) | first-usable | age only | **0.301 [0.189, 0.413]** | 9.3e-07 | 0.245 | |
| r2s_PUT | first-usable | + sex | 0.309 [0.195, 0.422] | 6.6e-07 | 0.251 | |
| r2s_PUT | first-usable | + sex + cell | 0.327 [0.211, 0.442] | 3.1e-07 | 0.265 | |
| r2s_PUT | all usable, lme | age only | 0.299 [0.203, 0.394] | 7.9e-07 | 0.243 | |
| r2s_PUT | all usable, lme | + sex | 0.305 [0.210, 0.400] | 5.2e-07 | 0.248 | |
| r2s_PUT | all usable, lme | + sex + cell | 0.314 [0.218, 0.411] | 4.2e-07 | 0.255 | |
| r2s_PUT | lme, within-person | + sex + cell | 0.297 [0.087, 0.506] | .0072 | 0.241 | |
| r2s_PUT | lme, between-person | + sex + cell | 0.319 [0.213, 0.425] | 6.6e-08 | 0.259 | |
| **chi_PUT** (ppb) (= chi_PUT_raw) | first-usable | age only | **0.382 [-0.262, 1.026]** | .2409 | 0.065 | yes |
| chi_PUT | first-usable | + sex | 0.424 [-0.203, 1.050] | .1820 | 0.072 | yes |
| chi_PUT | first-usable | + sex + cell | 0.501 [-0.168, 1.169] | .1402 | 0.085 | yes |
| chi_PUT | all usable, lme | age only | 0.080 [-0.562, 0.721] | .8005 | 0.013 | yes |
| chi_PUT | all usable, lme | + sex | 0.145 [-0.488, 0.778] | .6419 | 0.025 | yes |
| chi_PUT | all usable, lme | + sex + cell | 0.119 [-0.532, 0.770] | .7092 | 0.020 | yes |
| chi_PUT | lme, within-person | + sex + cell | -0.095 [-1.917, 1.728] | .9158 | -0.016 | no (uninformative) |
| chi_PUT | lme, between-person | + sex + cell | 0.148 [-0.530, 0.826] | .6643 | 0.025 | yes |
| **r2s_GP** (s^-1) | first-usable | age only | **0.687 [0.414, 0.961]** | 3.3e-06 | 0.235 | |
| r2s_GP | first-usable | + sex | 0.707 [0.421, 0.993] | 4.7e-06 | 0.242 | |
| r2s_GP | first-usable | + sex + cell | 0.708 [0.403, 1.013] | 1.5e-05 | 0.242 | |
| r2s_GP | all usable, lme | age only | 0.588 [0.366, 0.811] | 1.1e-05 | 0.201 | |
| r2s_GP | all usable, lme | + sex | 0.603 [0.383, 0.824] | 6.5e-06 | 0.207 | |
| r2s_GP | all usable, lme | + sex + cell | 0.612 [0.390, 0.834] | 5.7e-06 | 0.210 | |
| r2s_GP | lme, within-person | + sex + cell | 0.444 [0.051, 0.837] | .0284 | 0.152 | |
| r2s_GP | lme, between-person | + sex + cell | 0.690 [0.431, 0.950] | 1.2e-06 | 0.236 | |
| **chi_GP** (ppb) (= chi_GP_raw) | first-usable | age only | **2.397 [1.155, 3.638]** | .0002 | 0.187 | no |
| chi_GP | first-usable | + sex | 2.473 [1.188, 3.758] | .0003 | 0.193 | no |
| chi_GP | first-usable | + sex + cell | 2.337 [0.896, 3.779] | .0018 | 0.182 | no |
| chi_GP | all usable, lme | age only | 1.914 [0.789, 3.040] | .0017 | 0.149 | no |
| chi_GP | all usable, lme | + sex | 1.993 [0.876, 3.109] | .0011 | 0.155 | no |
| chi_GP | all usable, lme | + sex + cell | 1.852 [0.707, 2.997] | .0026 | 0.144 | no |
| chi_GP | lme, within-person | + sex + cell | 0.843 [-1.318, 3.004] | .4298 | 0.066 | no |
| chi_GP | lme, between-person | + sex + cell | 2.243 [0.937, 3.549] | .0010 | 0.175 | no |

### D.2 Is the chi slope what the R2* slope implies? (Dr = 274; if all of the R2* change were paramagnetic)

| ROI | sample, covariates | R2* slope [95% CI] | implied chi slope = 1000·R2* slope/274 (ppb/y) | observed chi slope (ppb/y) |
|---|---|---|---|---|
| putamen | first-usable, age only | 0.301 [0.189, 0.413] | **1.10 [0.69, 1.51]** | **0.38 [-0.26, 1.03]** |
| putamen | all usable lme, + sex + cell | 0.314 [0.218, 0.411] | 1.15 [0.79, 1.50] | 0.12 [-0.53, 0.77] |
| putamen | within-person, + sex + cell | 0.297 [0.087, 0.506] | 1.08 [0.32, 1.85] | -0.09 [-1.92, 1.73] |
| GP | first-usable, age only | 0.687 [0.414, 0.961] | **2.51 [1.51, 3.51]** | **2.40 [1.15, 3.64]** |
| GP | all usable lme, + sex + cell | 0.612 [0.390, 0.834] | 2.23 [1.43, 3.04] | 1.85 [0.71, 3.00] |
| GP | within-person, + sex + cell | 0.444 [0.051, 0.837] | 1.62 [0.18, 3.06] | 0.84 [-1.32, 3.00] |

A model-dependent shortfall estimate is the age slope of chi_dia = (chi - 1000·R2*/Dr)/2, which equals (observed - implied)/2 under the X2 model. It is a post hoc diagnostic, not a test of a registered hypothesis, and it depends on Dr and on the covariate set (`D_putamen_shortfall_by_Dr.csv`):
- Putamen chi_dia at Dr = 274, first-usable (OLS, HC3): -0.358 [-0.715, -0.001] ppb/y, p .0492 (age only); -0.352, p .0512 (+ sex); -0.346, p .0551 (+ sex + cell).
- Putamen chi_dia at Dr = 274, all usable (lme): -0.501 [-0.840, -0.161], p .0055 (age only); -0.481, p .0076 (+ sex); -0.515 [-0.869, -0.160], p .0061 (+ sex + cell, the covariate set of the lme rows in the table above). Within-person (+ sex + cell): -0.589 [-1.570, 0.393], p .2290.
- GP chi_dia (age only): -0.056 [-0.561, 0.449], p .83 (first-usable); -0.111 [-0.555, 0.333], p .61 (all usable, lme).
- Putamen chi_para (first-usable, age only): +0.740 [0.335, 1.145], p .0005. GP chi_para: +2.453 [1.446, 3.459], p < .0001.
- Rows with sex, cell and the within/between split are in `D_age_slopes.csv`.

Putamen shortfall estimate (age slope of chi_dia_PUT, ppb/y, and p) by Dr and covariate set:

| Dr (Hz/ppm) | first-usable: age only | + sex | + sex + cell | all usable lme: age only | + sex | + sex + cell |
|---|---|---|---|---|---|---|
| 137 | -0.907, .0004 | -0.915, .0004 | -0.942, .0002 | -1.033, 4.8e-05 | -1.030, 5.7e-05 | -1.080, 4.9e-05 |
| 274 | -0.358, .0492 | -0.352, .0512 | -0.346, .0551 | -0.501, .0055 | -0.481, .0076 | -0.515, .0061 |
| 321 | -0.278, .1125 | -0.269, .1192 | -0.258, .1372 | -0.422, .0144 | -0.401, .0200 | -0.432, .0160 |
| 548 | -0.083, .6093 | -0.070, .6623 | -0.048, .7724 | -0.233, .1447 | -0.206, .1944 | -0.230, .1617 |

In the first-usable sample the estimate is borderline at Dr = 274 (p .049-.055 across covariate sets) and does not hold at Dr = 321 or 548. In all usable sessions it holds for Dr ≤ 321 (p ≤ .020) but not at Dr = 548 (p .14-.19). Dr = 548 is physically implausible for these data (A.3).

Correlation of chi and R2* across the 80 first-usable sessions: putamen r = .21 (age-partialled .16); GP r = .68.

### D.3 Referencing co-reports (first-usable, age only; all variants in the CSV)

| quantity | slope per y [95% CI] | p |
|---|---|---|
| whole-brain reference value ref_ppb | 0.042 [-0.100, 0.184] ppb/y | .554 |
| CSF offset (ref2_minus_wb_ppb) | -0.003 [-0.738, 0.732] ppb/y | .994 |
| chi_PUT before whole-brain referencing (chi_PUT + ref_ppb) | 0.425 [-0.315, 1.164] | .257 |
| chi_PUT relative to CSF (chi_PUT - ref2_minus_wb_ppb) | 0.385 [-0.410, 1.179] | .338 |
| chi_GP before referencing / relative to CSF | 2.439 [1.129, 3.748] / 2.399 [1.071, 3.728] | .0004 / .0006 |
| co-report: SN-VTA chi / R2* | 1.312 [-0.070, 2.695] ppb/y / 0.311 [0.052, 0.571] s^-1/y | .063 / .019 |

### D.4 Literature

Larsen et al. 2020, *J Neurosci* 40(9):1810-1818, doi:10.1523/JNEUROSCI.2434-19.2020 (PMC7046454):
- Sample: Philadelphia Neurodevelopmental Cohort, 922 youths aged 8-26 y at the first visit, 1,236 sessions, 3 T.
- Method: R2* from a **two-echo** GRE (TE 2.69 and 5.27 ms in the main protocol; 3.75 mm in-plane), so absolute values are not comparable with ours. Their estimated marginal mean R2* is 16.96 s^-1 in putamen and 22.33 in GP; ours are 18.60 and 30.53.
- Findings: GAMM derivatives show R2* rising significantly at ages **10.2-26.9 y in the putamen**. In the GP the windows are 8.2-25.9 y (males) and 8.2-22.2 y (females). Age F statistics are 85.8 (putamen) and 76.7 (GP).
- They give **no numeric per-year slope** in the text; trajectories appear only in figures. Our putamen R2* rise of 0.30 s^-1/y (0.24 SD/y) over 8.8-17.8 y therefore agrees with their direction and age window, but its magnitude cannot be checked against them.

### D.5 Reading (D): the putamen chi slope is low relative to putamen R2*, while GP chi tracks GP R2*

- **The magnitude signal develops in the same sessions.** Putamen R2* rises by 0.30 s^-1/y (p ≈ 1e-6), including within persons (0.30, p .007). Its standardised rate (0.24 SD/y) is similar to that of GP R2* (0.24 SD/y). The sample and its age range (8.8-17.8 / 19.5 y) therefore carry a clear developmental R2* signal in the putamen (iron-sensitive, not iron-specific: myelin and water also contribute to R2*), consistent with Larsen 2020.
- **Putamen chi does not follow.** Its slope is 0.38 ppb/y (first-usable) and 0.08-0.14 ppb/y (all usable). The 95% CIs exclude the registered 1.33-1.74 ppb/y in every cross-sectional and mixed model; the within-person slope is uninformative. The point estimate is also below the 1.10 ppb/y implied by its own R2* slope. The model-dependent shortfall estimate (D.2) is borderline in the first-usable sample (p .049-.055 at Dr = 274, depending on covariates; p ≥ .11 at Dr = 321 and 548) and clearer in all usable sessions (p .0055-.0076 at Dr = 274; .014-.020 at Dr = 321; .14-.19 at Dr = 548).
- **In GP, chi tracks R2* as expected**: observed 2.40 vs implied 2.51 ppb/y, within the registered range, with chi and R2* strongly correlated (r = .68 vs .21 in putamen). We found no sign of a global QSM failure.
- **We found no sign that referencing drives it.** Neither reference value changes with age (p .55 and .99; D.3), and unreferenced or CSF-referenced putamen chi gives the same flat slope.
- **Two explanations remain, and these data cannot separate them.**
  - (a) A putamen-specific QSM measurement problem (inversion or regularisation, background-field removal, ROI definition, or partial volume) that leaves the magnitude decay intact.
  - (b) A real concurrent diamagnetic increase in the putamen (for example, myelination of striatal fibre bundles). This would raise R2* while cancelling the chi rise. The X2 model "explains" the shortfall exactly this way (chi_para up 0.74, chi_dia down 0.36 ppb/y), but that is a restatement under the model's assumptions, not evidence.
- The registered caveat ("developmental positive control not confirmed") stays on every H1 statement. The putamen chi slope is low relative to putamen R2*, while GP chi tracks GP R2*. This is consistent with a putamen-specific QSM measurement issue (a) or with a real diamagnetic change (b); these data cannot separate them, and (b) would be neither a QSM problem nor a failure. Because putamen R2* develops in the same sessions, a sample or age-range explanation is unlikely, but these diagnostics cannot show whether the issue, if it is one, also affects the SN-VTA.
- SN-VTA co-report: R2* rises with age (0.31 [0.05, 0.57] s^-1/y first-usable, p .019; 0.50 [0.28, 0.71] lme; first-usable with + sex + cell 0.25, p .053). The SN-VTA chi slope is positive but imprecise: 1.31 [-0.07, 2.69] ppb/y first-usable (p .063); 1.24 [0.21, 2.27] lme; first-usable with + sex + cell 1.02 [-0.70, 2.75] (p .24); within-person 0.33 [-1.69, 2.34].

---

## Open decisions (for the authors)

This folder does not settle these three points; each needs an author decision.

1. **nlme convergence tolerance (D).** The D mixed-model rows keep nlme's default nlminb tolerance. They differ from the independent REML implementation by up to 3.3e-5 relative (174 cells above the strict 1e-6 test), no printed digit depends on this, and tighter tolerances make nlme fail (see "Verification after the fixes"). Decide whether to accept this residual as it is.
2. **R2* voxel set (A.1, declaration 7).** The departure from the plan text is not yet in the manuscript deviation table (Supplementary Table S1 / deviation annex, Annex C.4). Adding it, and its wording, is for the authors. It also concerns the frozen K2 and PC3. Its size can only be computed from the voxel maps on the server; that was not done.
3. **BH family for X2 (A.2).** Choose the family the manuscript uses (X2 only, or pooled with the registered Tier-4 tests) and say that it was chosen after unlock. No X2 test survives under any of the four definitions tried (smallest q .096 / .137 / .177 / .175).

---

## Files

| file | content |
|---|---|
| `code/validation_local.R` | the whole analysis (A-D). Sources the frozen helpers read-only; writes only into this folder |
| `code/verify_statsmodels.py` | independent statsmodels OLS + HC3 re-fit of 5 exported designs |
| `results.json` | every estimate in this README plus metadata, input and code hashes and the reproduction checks |
| `A_X2_chisep_group_effects.csv` | 48 rows: inputs, chi_para / chi_dia for 3 ROIs x 4 Dr values, R2 = 0/10/15 cancellation rows, unreferenced rows, BH q values, cell means, constraint-violation shares, identity check |
| `A_X2_component_descriptives.csv` | component means, SDs, constraint violations and r(para, dia) by ROI and Dr |
| `A_X2_BH_family_sensitivity.csv` | X2 q values under four BH family definitions (X2 only 6 / 18; pooled with the 5 registered Tier-4 tests: 11 / 23). Added 2026-10-07 |
| `B_nuisance_models.csv`, `B_joint_model_coefficients.csv`, `B_r2s_models.csv`, `B_dice_tertile_models.csv`, `B_dice_tertile_cells.csv` | section B |
| `C_low_TD_tail_sessions.csv`, `C_low_TD_tail_qc_comparison.csv`, `C_low_TD_tail_refits.csv` | section C (coded IDs only) |
| `D_age_slopes.csv`, `D_implied_chi_slope_from_R2s.csv` | section D (18 outcomes x samples x covariate sets) |
| `D_putamen_shortfall_by_Dr.csv` | putamen shortfall estimate (chi_dia_PUT age slope) by Dr (137 / 274 / 321 / 548), sample and covariates. Added 2026-10-07 |
| `verify/` | exported design matrices, R reference estimates, `statsmodels_check.json` |
| `run.log`, `sessionInfo.txt`, `SHA256SUMS` | run log, R session, checksums of every file in this folder and its subfolders except `SHA256SUMS` itself (`archive_v1_prefix_2026-10-07/` also keeps its own original `SHA256SUMS`) |
| `CHANGES_2026-10-07.md` | every fix of the 2026-10-07 revision, with its review item id (round 1, the recheck addendum and Round 2) |
| `RECHECK_2026-10-07.md` | the independent re-check of every review item after the round-1 fixes; its line numbers refer to the round-1 README |
| `archive_v1_prefix_2026-10-07/` | the pre-fix README, CSVs, results.json, code, logs and the reviewer's pre-fix comparison outputs, unchanged |
| `archive_v2_round1_2026-10-07/` | the round-1 README, `CHANGES_2026-10-07.md` and `SHA256SUMS`, and the round-1 `check_readme_v2.py` with its outputs, as they were before Round 2, unchanged |
| `verification/` | the independent review (`VERIFY.md`), its Python re-implementation and comparison outputs (re-run after the fixes); `verification/rerun_2026-10-07/` holds the post-fix README check and the check of the two new CSVs |

## Limits

- Nothing here is a finding. Everything is post hoc except X2, and even X2's constants were fixed after unlock. This is a declared deviation from Q25: the constants are the paper's values, not tuned.
- The R2* ROI means use a different voxel set from the one in the plan text (A.1, declaration 7). This applies to every R2* number here and to the frozen K2 and PC3.
- The BH family for X2 is a post hoc choice; no X2 test survives under any of the four definitions tried (A.2).
- The C refits select on the outcome. The quality-matched subsample in B.3 and all within-tertile estimates are small and imprecise.
- X2 is a deterministic linear recombination of chi and R2*. Its numbers depend on Dr and on the model's assumptions (single field strength, no R2 map, constant Dr, no voxelwise constraints or regularisation).
- The D mixed models use nlme's default convergence tolerance. They agree with the independent REML implementation to ≤ 3.3e-5 relative, which no printed digit depends on (see below).
- No server access, no images, no installs; only the local merged table (and, read-only, the Tier-4 FDR table) was used.

## Verification after the fixes (2026-10-07)

- `verification/verify_independent.py` (the reviewer's Python re-implementation from the raw table, unchanged) was re-run on the revised outputs. All A-C cells agree (no cell differs by more than 1e-6 relative; the A-C CSVs are byte-identical to the pre-fix run).
- D mixed-model rows: the largest difference fell from 6.2e-3 relative (2.2e-4 ppb/y on a slope) to 3.3e-5 relative; no slope or CI bound differs by more than 1.2e-5 and no p value by more than 5.1e-6 in absolute terms. 174 cells (144 in `D_age_slopes.csv`, 30 derived ones in `D_implied_chi_slope_from_R2s.csv`) still exceed the verifier's strict 1e-6 relative threshold. The cause is nlminb's default convergence tolerance: at the verifier's optimum the REML criterion is up to about 1e-9 lower than at nlme's (checked on three models; `verification/rerun_2026-10-07/diag_reml_criterion.R`). Tightening nlminb made nlme stop with "singular convergence" in 35 (rel.tol = 1e-12) or 47 (rel.tol = 1e-13) of the 72 lme fits (`verification/rerun_2026-10-07/diag_nlminb_tolerance.R`), so the default was kept. No "CI excludes registered range" flag changes.
- `verification/rerun_2026-10-07/check_new_outputs.py` checks the two new CSVs against independent Python computations (the verifier's own shortfall-by-Dr table and an own BH implementation).
- `verification/rerun_2026-10-07/check_readme_v2.py` re-checks every number printed in this README against the verifier's recomputation (it is the reviewer's `check_readme.py`, adapted to the revised line layout and new sentences). Result after Round 2: 1,550 printed numbers and statements checked, 0 fail (`verification/rerun_2026-10-07/compare_readme_v2.csv`). These are the 1,532 checks of round 1, all still passing, plus 18 for the Round 2 wording; the Dice values of the label-free-cut sessions in Reading (C) are recomputed from the raw table. The reviewer's original `verification/check_readme.py` and its outputs `verification/compare_readme*.csv` refer to the pre-fix README (line numbers of `archive_v1_prefix_2026-10-07/README.md`) and were not re-run.
