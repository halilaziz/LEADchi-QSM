# LEAD-chi confirmatory analysis report

- Code: LEAD-chi analysis v1.0-final (2026-09-26); run 2026-10-28T20:21:03Z (UTC); R version 4.4.1 (2024-06-14)
- Environment pin: all pinned versions match
- Input: `/path/to/Penn_LEAD/LEADchi/results_posthoc/humanqc_primary/input/merged_H_all180.csv` (SHA-256 `55b82d06b18382085fd04ee09050199cc306a6876a0a528582648366da5c9d02`)
- Gates: `/path/to/Penn_LEAD/LEADchi/mirror/gates/gates_preunlock.json` (SHA-256 `dcf281b8002909de58566e2af6f75c7fede7319966e2b3412932c5fb33377157`)
- Branch from gates: engine **qsmxt**, pipeline **primary**, ROI **SNVTA**; primary interpretable: **TRUE**; caveats: PC2w indeterminate: PC2x evaluated after unlock; developmental positive control not confirmed (applies to every H1 statement)
- Resampling counts: n_perm = 10000, n_perm_spec = 10000, n_boot_delta = 5000, n_boot_s3b = 5000, n_pb = 2000, n_wild = 9999 (registered)

## 1. Samples (TD / ADHD-only / PS-only / PS+ADHD)

| sample | TD | ADHD-only | PS-only | PS+ADHD | total |
|---|---|---|---|---|---|
| All MEGRE sessions in table | 51 | 38 | 33 | 59 | 181 |
| Usable sessions | 37 | 23 | 22 | 38 | 120 |
| H1 / S1 first-usable | 23 | 18 | 17 | 27 |  85 |
| S2 (age >= 11, PRIME+BDI+ARI) | 15 | 13 | 17 | 25 |  70 |
| S3 (same-session DWI, >= n_min nodes) |  0 |  0 |  0 |  0 |   0 |
| T1 occasion 1 (held-out follow-ups) | 14 |  4 |  5 |  8 |  31 |
| ses-1-only first-usable | 19 | 11 | 15 | 16 |  61 |
| Exposure sensitivity, H1 (4 inspected sessions dropped) | 22 | 18 | 17 | 26 |  83 |
| Exposure sensitivity, S3 |  0 |  0 |  0 |  0 |   0 |

QC-by-group (subject level; Fisher's exact): not in H1 p = 0.8298; earliest MEGRE session failed p = 0.5645; IPW-M1 not triggered

## 2. Plausibility re-verification (sec. 12.4)

| sess_key | column | value | rule | action |
|---|---|---|---|---|
| sub-21038_ses-2 | chi_SNVTA |  -53.59083 | envelope [-50, 300] ppb | none (session already not usable) |
| sub-20953_ses-2 | chi_DN_L | -313.79333 | envelope [-50, 300] ppb | set to missing |

Label-free robust z (first-usable; |z| > 5 flagged, retained):

| roi | n | n_flag_absz_gt5 | sd_over_robust_sd |
|---|---|---|---|
| chi_SNVTA | 85 | 1 | 1.496725 |
| chi_STR | 85 | 0 | 1.190868 |
| chi_DN_R | 59 | 0 | 1.154536 |
| chi_DGM | 85 | 2 | 1.215498 |
| chi_PUT_R | 85 | 1 | 1.202187 |
| r2s_SNVTA | 85 | 0 | 1.063068 |
| chi_GP | 85 | 0 | 1.242465 |

## 3. H1 (M1; HC3; two-sided alpha .05)

`chi_SNVTA ~ PS_e * ADHD_e + ns(age, df = 3) + sex_e + race3 + scale(roi_vol) + scale(Q_phase)`, n = 85

| Quantity | Value |
|---|---|
| beta_PS (ppb) [95% CI]; [90% CI] | 6.281 [-0.226, 12.787]; [0.842, 11.720] |
| SE_HC3; t; df; two-sided p | 3.265; 1.924; 73; 0.0583 |
| sigma_resid; sigma_rob (Huber proposal-2 scale); sigma_B = min | 13.066; 11.667; 11.667 |
| delta = beta / sigma_resid [cell-stratified bootstrap 95% CI] | 0.481 [0.062, 1.043] |
| noncentral-t 95% CI of delta; raw Hedges g | [-0.017, 0.975]; 0.489 |
| B: r_L; B evaluable; delta_B; bound in ppb (delta_B x sigma_B) | 0.938; TRUE; 0.581; -6.781 |
| t_B; one-sided p_B | 4.001; <.0001 |
| B co-reports (never decisive): bound -0.6; -0.5; cross-session bound; residualised split-echo bound | <.0001; 0.0002; 0.0004; NA |
| Freedman-Lane p (10000 permutations within ADHD strata) [MC 95% CI] | 0.0436 [0.0396, 0.0477] |
| Huber M-estimate of beta_PS (SE) | 4.918 (2.874) |
| BF10 under N(-0.35, 0.25) (attenuated); BF01 | 0.284; 3.518 |
| BF10 under N(0, 0.5) (sceptical); BF01 | 1.965; 0.509 |
| sceptical-prior shrunken delta [95% CrI] | 0.385 [-0.053, 0.823] |
| adult reference (Vano 2025; descriptive; reference, resolution, field and pipeline differ) | -14.24 ppb |

**M1 outcome (whole-brain reference): B**

### K3 reference agreement (CSF-only secondary reference)

beta_PS(dref) = 0.199 ppb (p 0.9447); threshold max(0.25|beta|, 0.10 sigma) = 1.570 and p < 0.1; trigger: **FALSE**. Sessions without a valid CSF reference: 0. CSF-referenced outcome: B (beta 6.082, p 0.1211, p_B 0.0002). Rule: not triggered.

**Final H1 outcome (M1 -> K3): B.**

### K1 specificity labels (precedence K1c > K1b > K1a)

| contrast | est_ppb | se | p_two |
|---|---|---|---|
| PS-only - TD (mandatory co-report) | 10.765 | 5.191 | 0.0416 |
| PS+ADHD - ADHD-only | 1.797 | 4.179 | 0.6685 |
| ADHD main effect | 1.201 | 3.278 | 0.7151 |
| PS x ADHD interaction | -8.968 | 6.796 | 0.1911 |

ADHD TOST at +/-0.5 sigma_B (+/-5.833 ppb): p = 0.0176 / 0.0809; success: FALSE. Interaction 95% CI includes 0: TRUE (non-rejection only).

**K1 label: none**. Mandatory caveat: all 6 threshold-psychosis youths sit in PS+ADHD.

### K2 iron vs myelin

R2* beta_PS = 0.684 s^-1 (n = 85), one-sided p = 0.8948.

### K4 precision, influence, fragility and robustness

- PC4: TRUE; r_L 0.938; B evaluable: TRUE; label: none; attenuation-corrected delta = 0.490.
- Dominance: top-1 share 0.116, top-2 share 0.211, max |t_i| 3.371, max Cook's D 0.144; triggered: FALSE; **influence-dependent: FALSE**.
- Leave-one-out (M1 A/B/C outcome (whole-brain reference, before K3)): 0 of 85 deletions change the outcome; **fragile: FALSE**.
- ses-1-only (mandatory co-report): n = 61, beta 6.117 ppb, p 0.1572, outcome B.
- Exposure sensitivity (Q18(4); 2 of the 4 listed sessions in this sample: sub-20139_ses-1, sub-21432_ses-1): n = 83, beta 7.128, p 0.0270, outcome B; changes conclusion: FALSE (flagged, never overrides).
- Specification curve: 486 of 486 evaluable; median delta 0.341; share delta < 0 0.109; share significant 0.000; Stouffer Z 26.488; joint permutation p (median t, 10000) = 0.9777; **robust: FALSE** (final H1 outcome is not A).
- PC2x (age -> putamen / GP chi, one-sided): chi_PUT slope 0.639 ppb/y p 0.0550; chi_GP slope 2.780 ppb/y p <.0001; pass: FALSE.

### H1 sensitivity battery (sec. 24.4; reported separately, never decisive)

| sensitivity | n | beta_ppb | se_hc3 | p_two | pB | outcome | note |
|---|---|---|---|---|---|---|---|
| without the 6 threshold-psychosis youths | 80 | 5.829 | 3.269 | 0.0790 | 0.0001 | B |  |
| without 4 TD youths (dx-coding gap / prescribed stimulants) | 82 | 6.824 | 3.304 | 0.0426 | <.0001 | B |  |
| without prescribed-stimulant endorsers (same session) | 81 | 4.914 | 3.049 | 0.1116 | 0.0001 | B |  |
| without T1w 'artifact' sessions | 78 | 5.501 | 3.468 | 0.1174 | 0.0004 | B |  |
| without sub-21568 (no T1-QC entry) | 84 | 7.252 | 3.188 | 0.0259 | <.0001 | B |  |
| without early-epoch sessions | 77 | 7.272 | 3.222 | 0.0274 | <.0001 | B |  |
| baseline age >= 11 only | 64 | 6.393 | 4.887 | 0.1965 | 0.0044 | B |  |
| exposure-sensitivity set (4 inspected sessions dropped; Q18) | 83 | 7.128 | 3.157 | 0.0270 | <.0001 | B |  |
| without the (c') plausibility exclusions | 85 | 6.281 | 3.265 | 0.0583 | <.0001 | B |  |
| implausible voxels retained | 85 | 6.281 | 3.265 | 0.0583 | <.0001 | B |  |
| without rubric code D (dental/metal) sessions | 82 | 5.548 | 3.268 | 0.0940 | 0.0002 | B |  |
| strict QC (visual grade 0 only) | 82 | 6.725 | 3.368 | 0.0498 | <.0001 | B |  |
| automated exclusions only (no visual grade or placement rating) | 90 | 5.928 | 3.099 | 0.0595 | <.0001 | B |  |
| without race | 85 | 6.194 | 3.159 | 0.0536 | <.0001 | B | race3 removed |
| + Hispanic ethnicity | 85 | 6.386 | 3.349 | 0.0606 | <.0001 | B |  |
| + puberty (PDS), complete case | 60 | 8.005 | 3.967 | 0.0493 | 0.0002 | B |  |
| + puberty (PDS), one stochastic regression imputation (mice fallback) | 74 | 7.365 | 3.610 | 0.0456 | <.0001 | B | observed 60, imputed 14 |
| + menarche (girls; 0 for boys) | 75 | 6.432 | 3.383 | 0.0619 | <.0001 | B |  |
| + intracranial volume (SynthSeg, coded side) | 85 | 6.107 | 3.387 | 0.0756 | 0.0001 | B |  |
| + mood diagnosis (MDD or mood NOS) | 85 | 6.353 | 3.584 | 0.0806 | 0.0002 | B |  |
| + protocol epoch | 85 | 6.450 | 3.235 | 0.0500 | <.0001 | B |  |
| + f0 cluster | 85 | 5.831 | 3.291 | 0.0807 | <.0001 | B |  |
| + years since the diagnostic assessment | 85 | 6.260 | 3.350 | 0.0657 | <.0001 | B |  |
| + |pitch| | 85 | 6.237 | 3.241 | 0.0582 | <.0001 | B |  |
| + DWI FD | NA | NA | NA | NA | NA | not run | dwi_FD absent |
| + registration Dice (minimum of the G5 Dice values) | 85 | 4.710 | 3.329 | 0.1614 | 0.0010 | B |  |
| + magnitude ghost-to-signal motion index | 85 | 5.551 | 3.315 | 0.0984 | 0.0002 | B |  |
| + reference-domain WM and CSF fractions and brain-in-FOV fraction | 85 | 6.644 | 3.467 | 0.0594 | <.0001 | B |  |
| coverage-harmonised whole-brain reference | 85 | 6.281 | 3.267 | 0.0584 | <.0001 | B |  |
| SN-only complex (descriptive) | 85 | 6.994 | 3.452 | 0.0464 | <.0001 | B |  |
| VTA + PBP only (descriptive) | 85 | 3.161 | 3.277 | 0.3378 | 0.0073 | B |  |
| rank-based ANCOVA (ranks of chi_SNVTA) | 85 | NA | NA | 0.1336 | NA | rank test not significant | beta in rank units not shown |
| other engine | NA | NA | NA | NA | NA | not run | no other-engine table supplied (A2 not approved or not run) |

## 4. Key secondary (Holm; family: S1, S2; S3 status: exploratory (a dentate gate or PC5 failed))

| test | n | est | se_hc3 | p_two | p_holm | supported | influence_dependent |
|---|---|---|---|---|---|---|---|
| S1 L (SN-VTA vs STR x PS vs ADHD) | 85 | 0.616 | 0.439 | 0.1649 | 0.1649 | FALSE | FALSE |
| S2 beta_PRIME (per SD) | 70 | 4.613 | 2.474 | 0.0673 | 0.1347 | FALSE | FALSE |
| S3a beta(chi_DN_R) (standardised) |  0 | NA | NA | NA | not in family | NA | FALSE |

- **S1** pooled constants: gates JSON (label-free, computed before unlock by gates_preunlock.py). beta_PS(D) = 0.594 (one-sided p 0.9679); beta_ADHD(D) = -0.022 (one-sided p 0.5287). Wording: **L not significant under Holm (no S1 claim)**. Stacked GLS L = 0.616 (p 0.1463).
- **S2** + PS_e: beta_PRIME = 4.061 (p 0.0746); 'beyond diagnosis' allowed: FALSE. Within PS: -0.471 (p 0.8792); rank-INT: 4.898 (p 0.0279).
  - S2 sensitivity PRIME_items_1_12: beta 4.613 (p 0.0673)
  - S2 sensitivity CPT_in_place_of_ADHD: beta 4.068 (p 0.1477)
  - S2 all-session between-person PRIME: 3.342 (SE 2.377); within-person 7.275 (SE 5.526)
- **S3** confirmatory gates (PC4-DN, G2-DN, PC5): FALSE. Freedman-Lane p = NA. NC1 (DN - PUT) = NA, one-sided 95% lower bound NA; NC2 (pSCP - MCP matched) = NA, lower bound NA. Wording: **no S3a claim**.

## 5. T1 temporal stability (held-out follow-ups; Delta = 0.4 total SD = 4.752 ppb)

REML maximum: fit `lmer_bobyqa`; ICC 0.179; boundary: FALSE; profile confirms the global maximum on the grid: TRUE.

| fit | theta | reml_deviance | gap_to_best |
|---|---|---|---|
| lmer_bobyqa | 0.46740 | 367.340990 | 0.00000000 |
| lmer_NelderMead | 0.46740 | 367.340990 | 0.00000000 |
| lme_nlminb | 0.46740 | 367.340990 | 0.00000000 |
| lme_optim | 0.46783 | 367.340993 | 0.00000315 |

| Quantity | Value |
|---|---|
| rows / subjects / returners | 116 / 85 / 31 |
| beta(PS x occ) (total-SD units), SE, containment df | -0.172, 0.499, 24 |
| margin test H0: beta >= Delta: decision p (containment df); PB p (co-report, 2000) | 0.1315; 0.1384 |
| decline test H0: beta = 0: decision p (containment df); PB p (co-report) | 0.7329; 0.7321 |
| OLS cluster-robust (HC2; not CR2 for the LMM): beta; SE; wild cluster p (margin, 9999) | -0.181; 0.509; 0.1266 |
| occasion-0 beta_PS (total-SD units) vs Delta | 0.529: the occasion-0 deficit exceeds Delta |
| follow-up-only M1 beta_PS (ppb), p, n | 6.238, 0.4980, 31 |
| exposure sensitivity beta; p margin | -0.278; 0.1046 |

**T1 label: not evaluated (outcome A does not hold): follow-up-only estimate reported descriptively**. Assumption: two-occasion LMM; missing follow-ups assumed missing at random given the model covariates.

## 6. Pooled all-session sensitivity (not a test)

beta_PS = 5.370 ppb (SE 3.087, containment df 75, p 0.0860); REML fit `lmer_bobyqa`; OLS cluster-robust (HC2) SE 2.702; wild cluster p 0.0464.

## 7. Registered reading (verbatim registered sentences only)

- Outcome B: An adult-magnitude SN-VTA susceptibility deficit is not a feature of clinically ascertained, help-seeking psychosis-spectrum youth as a group. This does not exclude a deficit confined to the minority who will develop psychosis, nor a deficit that emerges later.
- If a trait deficit existed only in the minority who will convert to psychosis, the group-level effect would be about 0.09-0.18 SD; such a deficit can be neither detected nor excluded here.
- T1 was not evaluated because outcome A does not hold; follow-up estimates are descriptive.
- No topographic dissociation was established (S1 not significant under Holm); S1 was powered only for large contrasts.
- No association between SN-VTA chi and PRIME was established (S2 not significant under Holm).
- S3 is exploratory because a dentate gate or PC5 failed.
- Caveats (registered branch outcomes): PC2w indeterminate: PC2x evaluated after unlock; developmental positive control not confirmed (applies to every H1 statement).

## 8. Log

```
2026-10-28T20:21:04Z  LEAD-chi run_all.R start; input results_posthoc/humanqc_primary/input/merged_H_all180.csv
2026-10-28T20:21:04Z  branch: engine=qsmxt pipeline=primary roi=SNVTA
2026-10-28T20:21:04Z  samples: first-usable 23 / 18 / 17 / 27 (total 85); S2 15 / 13 / 17 / 25 (total 70); S3 0 / 0 / 0 / 0 (total 0); T1 occ1 14 / 4 / 5 / 8 (total 31)
2026-10-28T20:21:04Z  H1
2026-10-28T20:21:14Z  S1-S3
2026-10-28T20:21:14Z  T1
2026-10-28T20:21:34Z  warning in lmer Nelder_Mead: unable to evaluate scaled gradient
2026-10-28T20:21:34Z  warning in lmer Nelder_Mead: Model failed to converge: degenerate  Hessian with 1 negative eigenvalues
  See ?lme4::convergence and ?lme4::troubleshooting.
2026-10-28T20:21:34Z  robustness: sensitivity battery
2026-10-28T20:21:34Z  robustness: specification curve
2026-10-28T20:22:06Z  robustness: pooled all-session and S2 all-session models
```
