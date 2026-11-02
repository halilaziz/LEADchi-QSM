# LEAD-chi confirmatory analysis report

- Code: LEAD-chi analysis v1.0-final (2026-09-26); run 2026-10-05T12:46:25Z (UTC); R version 4.4.1 (2024-06-14)
- Environment pin: all pinned versions match
- Input: `/path/to/Penn_LEAD/LEADchi/mirror/unlocked/merged_qsmxt_primary_SNVTA.csv` (SHA-256 `f669bc6f41ab4e489abbd0a7c41949d3e5ccaba3e7f3c219d1bbebc3118df215`)
- Gates: `/path/to/Penn_LEAD/LEADchi/mirror/gates/gates_preunlock.json` (SHA-256 `dcf281b8002909de58566e2af6f75c7fede7319966e2b3412932c5fb33377157`)
- Branch from gates: engine **qsmxt**, pipeline **primary**, ROI **SNVTA**; primary interpretable: **TRUE**; caveats: PC2w indeterminate: PC2x evaluated after unlock; developmental positive control not confirmed (applies to every H1 statement)
- Resampling counts: n_perm = 10000, n_perm_spec = 10000, n_boot_delta = 5000, n_boot_s3b = 5000, n_pb = 2000, n_wild = 9999 (registered)

## 1. Samples (TD / ADHD-only / PS-only / PS+ADHD)

| sample | TD | ADHD-only | PS-only | PS+ADHD | total |
|---|---|---|---|---|---|
| All MEGRE sessions in table | 51 | 38 | 33 | 59 | 181 |
| Usable sessions | 31 | 19 | 22 | 35 | 107 |
| H1 / S1 first-usable | 22 | 16 | 17 | 25 |  80 |
| S2 (age >= 11, PRIME+BDI+ARI) | 15 | 11 | 17 | 22 |  65 |
| S3 (same-session DWI, >= n_min nodes) |  0 |  0 |  0 |  0 |   0 |
| T1 occasion 1 (held-out follow-ups) |  9 |  2 |  5 | 10 |  26 |
| ses-1-only first-usable | 19 | 10 | 15 | 14 |  58 |
| Exposure sensitivity, H1 (4 inspected sessions dropped) | 21 | 16 | 17 | 23 |  77 |
| Exposure sensitivity, S3 |  0 |  0 |  0 |  0 |   0 |

QC-by-group (subject level; Fisher's exact): not in H1 p = 0.9670; earliest MEGRE session failed p = 0.2439; IPW-M1 not triggered

## 2. Plausibility re-verification (sec. 12.4)

| sess_key | column | value | rule | action |
|---|---|---|---|---|
| sub-21038_ses-2 | chi_SNVTA |  -53.59083 | envelope [-50, 300] ppb | none (session already not usable) |
| sub-20953_ses-2 | chi_DN_L | -313.79333 | envelope [-50, 300] ppb | set to missing |

Label-free robust z (first-usable; |z| > 5 flagged, retained):

| roi | n | n_flag_absz_gt5 | sd_over_robust_sd |
|---|---|---|---|
| chi_SNVTA | 80 | 0 | 1.2417721 |
| chi_STR | 80 | 0 | 1.0872506 |
| chi_DN_R | 58 | 0 | 1.1017986 |
| chi_DGM | 80 | 0 | 0.9424925 |
| chi_PUT_R | 80 | 0 | 1.0408923 |
| r2s_SNVTA | 80 | 0 | 1.0922149 |
| chi_GP | 80 | 0 | 1.0896238 |

## 3. H1 (M1; HC3; two-sided alpha .05)

`chi_SNVTA ~ PS_e * ADHD_e + ns(age, df = 3) + sex_e + race3 + scale(roi_vol) + scale(Q_phase)`, n = 80

| Quantity | Value |
|---|---|
| beta_PS (ppb) [95% CI]; [90% CI] | 4.963 [-1.273, 11.199]; [-0.248, 10.174] |
| SE_HC3; t; df; two-sided p | 3.125; 1.588; 68; 0.1169 |
| sigma_resid; sigma_rob (Huber proposal-2 scale); sigma_B = min | 11.256; 10.550; 10.550 |
| delta = beta / sigma_resid [cell-stratified bootstrap 95% CI] | 0.441 [-0.041, 1.019] |
| noncentral-t 95% CI of delta; raw Hedges g | [-0.110, 0.988]; 0.415 |
| B: r_L; B evaluable; delta_B; bound in ppb (delta_B x sigma_B) | 0.938; TRUE; 0.581; -6.132 |
| t_B; one-sided p_B | 3.550; 0.0004 |
| B co-reports (never decisive): bound -0.6; -0.5; cross-session bound; residualised split-echo bound | 0.0003; 0.0008; 0.0015; 0.0263 |
| Freedman-Lane p (10000 permutations within ADHD strata) [MC 95% CI] | 0.0857 [0.0802, 0.0913] |
| Huber M-estimate of beta_PS (SE) | 3.936 (2.661) |
| BF10 under N(-0.35, 0.25) (attenuated); BF01 | 0.279; 3.585 |
| BF10 under N(0, 0.5) (sceptical); BF01 | 1.273; 0.786 |
| sceptical-prior shrunken delta [95% CrI] | 0.337 [-0.139, 0.813] |
| adult reference (Vano 2025; descriptive; reference, resolution, field and pipeline differ) | -14.24 ppb |

**M1 outcome (whole-brain reference): B**

### K3 reference agreement (CSF-only secondary reference)

beta_PS(dref) = 1.486 ppb (p 0.2742); threshold max(0.25|beta|, 0.10 sigma) = 1.241 and p < 0.1; trigger: **FALSE**. Sessions without a valid CSF reference: 0. CSF-referenced outcome: B (beta 3.477, p 0.3507, p_B 0.0021). Rule: not triggered.

**Final H1 outcome (M1 -> K3): B.**

### K1 specificity labels (precedence K1c > K1b > K1a)

| contrast | est_ppb | se | p_two |
|---|---|---|---|
| PS-only - TD (mandatory co-report) | 8.157 | 3.980 | 0.0443 |
| PS+ADHD - ADHD-only | 1.768 | 4.785 | 0.7129 |
| ADHD main effect | -1.272 | 3.106 | 0.6836 |
| PS x ADHD interaction | -6.389 | 6.199 | 0.3063 |

ADHD TOST at +/-0.5 sigma_B (+/-5.275 ppb): p = 0.1009 / 0.0194; success: FALSE. Interaction 95% CI includes 0: TRUE (non-rejection only).

**K1 label: none**. Mandatory caveat: all 6 threshold-psychosis youths sit in PS+ADHD.

### K2 iron vs myelin

R2* beta_PS = 1.098 s^-1 (n = 80), one-sided p = 0.9825.

### K4 precision, influence, fragility and robustness

- PC4: TRUE; r_L 0.938; B evaluable: TRUE; label: none; attenuation-corrected delta = 0.450.
- Dominance: top-1 share 0.160, top-2 share 0.225, max |t_i| 4.051, max Cook's D 0.253; triggered: FALSE; **influence-dependent: FALSE**.
- Leave-one-out (M1 A/B/C outcome (whole-brain reference, before K3)): 0 of 80 deletions change the outcome; **fragile: FALSE**.
- ses-1-only (mandatory co-report): n = 58, beta 6.478 ppb, p 0.1782, outcome B.
- Exposure sensitivity (Q18(4); 3 of the 4 listed sessions in this sample: sub-20139_ses-1, sub-21432_ses-1, sub-20976_ses-2): n = 77, beta 5.469, p 0.0816, outcome B; changes conclusion: FALSE (flagged, never overrides).
- Specification curve: 324 of 486 evaluable; median delta 0.300; share delta < 0 0.019; share significant 0.000; Stouffer Z 19.964; joint permutation p (median t, 10000) = 0.9243; **robust: FALSE** (not assessable (< 400 evaluable specifications)).
- PC2x (age -> putamen / GP chi, one-sided): chi_PUT slope 0.382 ppb/y p 0.1225; chi_GP slope 2.397 ppb/y p 0.0003; pass: FALSE.

### H1 sensitivity battery (sec. 24.4; reported separately, never decisive)

| sensitivity | n | beta_ppb | se_hc3 | p_two | pB | outcome | note |
|---|---|---|---|---|---|---|---|
| without the 6 threshold-psychosis youths | 77 | 5.177 | 3.210 | 0.1117 | 0.0004 | B |  |
| without 4 TD youths (dx-coding gap / prescribed stimulants) | 76 | 4.723 | 3.154 | 0.1392 | 0.0008 | B |  |
| without prescribed-stimulant endorsers (same session) | 75 | 3.036 | 2.832 | 0.2878 | 0.0012 | B |  |
| without T1w 'artifact' sessions | 74 | 4.569 | 3.375 | 0.1807 | 0.0011 | B |  |
| without sub-21568 (no T1-QC entry) | 80 | 4.963 | 3.125 | 0.1169 | 0.0004 | B |  |
| without early-epoch sessions | 71 | 5.580 | 2.989 | 0.0669 | 0.0002 | B |  |
| baseline age >= 11 only | 62 | 6.225 | 4.738 | 0.1949 | 0.0047 | B |  |
| exposure-sensitivity set (4 inspected sessions dropped; Q18) | 77 | 5.469 | 3.091 | 0.0816 | 0.0002 | B |  |
| without the (c') plausibility exclusions | 80 | 4.963 | 3.125 | 0.1169 | 0.0004 | B |  |
| implausible voxels retained | 80 | 4.963 | 3.125 | 0.1169 | 0.0004 | B |  |
| without rubric code D (dental/metal) sessions | 78 | 5.111 | 3.115 | 0.1056 | 0.0003 | B |  |
| strict QC (visual grade 0 only) |  4 | NA | NA | NA | NA | fit failed |  |
| automated exclusions only (no visual grade or placement rating) | 90 | 5.928 | 3.099 | 0.0595 | <.0001 | B |  |
| without race | 80 | 4.893 | 3.083 | 0.1170 | 0.0004 | B | race3 removed |
| + Hispanic ethnicity | 80 | 4.689 | 3.106 | 0.1359 | 0.0004 | B |  |
| + puberty (PDS), complete case | 57 | 6.915 | 3.544 | 0.0574 | 0.0002 | B |  |
| + puberty (PDS), one stochastic regression imputation (mice fallback) | 69 | 6.236 | 3.389 | 0.0711 | 0.0002 | B | observed 57, imputed 12 |
| + menarche (girls; 0 for boys) | 72 | 4.871 | 3.242 | 0.1383 | 0.0005 | B |  |
| + intracranial volume (SynthSeg, coded side) | 80 | 5.337 | 3.253 | 0.1056 | 0.0003 | B |  |
| + mood diagnosis (MDD or mood NOS) | 80 | 5.888 | 3.611 | 0.1077 | 0.0006 | B |  |
| + protocol epoch | 80 | 4.909 | 3.077 | 0.1154 | 0.0003 | B |  |
| + f0 cluster | 80 | 5.059 | 3.262 | 0.1257 | 0.0005 | B |  |
| + years since the diagnostic assessment | 80 | 5.227 | 3.194 | 0.1064 | 0.0004 | B |  |
| + |pitch| | 80 | 5.076 | 3.128 | 0.1093 | 0.0003 | B |  |
| + DWI FD | NA | NA | NA | NA | NA | not run | dwi_FD absent |
| + registration Dice (minimum of the G5 Dice values) | 80 | 3.749 | 3.106 | 0.2316 | 0.0019 | B |  |
| + magnitude ghost-to-signal motion index | 80 | 4.227 | 2.978 | 0.1604 | 0.0004 | B |  |
| + reference-domain WM and CSF fractions and brain-in-FOV fraction | 80 | 4.683 | 3.249 | 0.1542 | 0.0006 | B |  |
| coverage-harmonised whole-brain reference | 80 | 4.970 | 3.126 | 0.1165 | 0.0004 | B |  |
| SN-only complex (descriptive) | 80 | 5.874 | 3.310 | 0.0804 | 0.0001 | B |  |
| VTA + PBP only (descriptive) | 80 | 1.278 | 2.962 | 0.6676 | 0.0341 | B |  |
| rank-based ANCOVA (ranks of chi_SNVTA) | 80 | NA | NA | 0.1231 | NA | rank test not significant | beta in rank units not shown |
| other engine | NA | NA | NA | NA | NA | not run | no other-engine table supplied (A2 not approved or not run) |

## 4. Key secondary (Holm; family: S1, S2; S3 status: exploratory (a dentate gate or PC5 failed))

| test | n | est | se_hc3 | p_two | p_holm | supported | influence_dependent |
|---|---|---|---|---|---|---|---|
| S1 L (SN-VTA vs STR x PS vs ADHD) | 80 | 0.654 | 0.473 | 0.1713 | 0.1996 | FALSE | FALSE |
| S2 beta_PRIME (per SD) | 65 | 3.968 | 2.368 | 0.0998 | 0.1996 | FALSE | FALSE |
| S3a beta(chi_DN_R) (standardised) |  0 | NA | NA | NA | not in family | NA | FALSE |

- **S1** pooled constants: gates JSON (label-free, computed before unlock by gates_preunlock.py). beta_PS(D) = 0.570 (one-sided p 0.9543); beta_ADHD(D) = -0.084 (one-sided p 0.6030). Wording: **L not significant under Holm (no S1 claim)**. Stacked GLS L = 0.654 (p 0.1584).
- **S2** + PS_e: beta_PRIME = 2.818 (p 0.2289); 'beyond diagnosis' allowed: FALSE. Within PS: -1.608 (p 0.5868); rank-INT: 3.729 (p 0.0606).
  - S2 sensitivity PRIME_items_1_12: beta 3.968 (p 0.0998)
  - S2 sensitivity CPT_in_place_of_ADHD: beta 2.432 (p 0.3189)
  - S2 all-session between-person PRIME: 1.770 (SE 2.014); within-person 6.462 (SE 3.239)
- **S3** confirmatory gates (PC4-DN, G2-DN, PC5): FALSE. Freedman-Lane p = NA. NC1 (DN - PUT) = NA, one-sided 95% lower bound NA; NC2 (pSCP - MCP matched) = NA, lower bound NA. Wording: **no S3a claim**.

## 5. T1 temporal stability (held-out follow-ups; Delta = 0.4 total SD = 4.752 ppb)

REML maximum: fit `lmer_bobyqa`; ICC 0.501; boundary: FALSE; profile confirms the global maximum on the grid: TRUE.

| fit | theta | reml_deviance | gap_to_best |
|---|---|---|---|
| lmer_bobyqa | 1.00275 | 286.446963 | 0.00000000 |
| lmer_NelderMead | 1.00276 | 286.446963 | 0.00000000 |
| lme_nlminb | 1.00275 | 286.446963 | 0.00000000 |
| lme_optim | 1.00283 | 286.446963 | 0.00000006 |

| Quantity | Value |
|---|---|
| rows / subjects / returners | 106 / 80 / 26 |
| beta(PS x occ) (total-SD units), SE, containment df | 0.183, 0.358, 19 |
| margin test H0: beta >= Delta: decision p (containment df); PB p (co-report, 2000) | 0.2759; 0.2889 |
| decline test H0: beta = 0: decision p (containment df); PB p (co-report) | 0.6151; 0.6202 |
| OLS cluster-robust (HC2; not CR2 for the LMM): beta; SE; wild cluster p (margin, 9999) | -0.030; 0.399; 0.1414 |
| occasion-0 beta_PS (total-SD units) vs Delta | 0.418: the occasion-0 deficit exceeds Delta |
| follow-up-only M1 beta_PS (ppb), p, n | 5.730, 0.5175, 26 |
| exposure sensitivity beta; p margin | 0.064; 0.1935 |

**T1 label: not evaluated (outcome A does not hold): follow-up-only estimate reported descriptively**. Assumption: two-occasion LMM; missing follow-ups assumed missing at random given the model covariates.

## 6. Pooled all-session sensitivity (not a test)

beta_PS = 5.112 ppb (SE 2.624, containment df 70, p 0.0554); REML fit `lmer_bobyqa`; OLS cluster-robust (HC2) SE 2.674; wild cluster p 0.0667.

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
2026-10-05T12:46:26Z  LEAD-chi run_all.R start; input /path/to/Penn_LEAD/LEADchi/mirror/unlocked/merged_qsmxt_primary_SNVTA.csv
2026-10-05T12:46:26Z  branch: engine=qsmxt pipeline=primary roi=SNVTA
2026-10-05T12:46:26Z  samples: first-usable 22 / 16 / 17 / 25 (total 80); S2 15 / 11 / 17 / 22 (total 65); S3 0 / 0 / 0 / 0 (total 0); T1 occ1 9 / 2 / 5 / 10 (total 26)
2026-10-05T12:46:26Z  H1
2026-10-05T12:46:35Z  S1-S3
2026-10-05T12:46:36Z  T1
2026-10-05T12:46:54Z  robustness: sensitivity battery
2026-10-05T12:46:55Z  ERROR in sensitivity strict QC (visual grade 0 only): rank-deficient design: ADHD_e, ns(age, df = 3)3, sex_e, race31, race32, scale(roi_vol), scale(Q_phase), PS_e:ADHD_e
2026-10-05T12:46:55Z  robustness: specification curve
2026-10-05T12:47:07Z  robustness: pooled all-session and S2 all-session models
```
