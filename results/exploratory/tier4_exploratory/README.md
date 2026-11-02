# LEAD-chi T4: Tier-4 exploratory modules (X1, X3, X6) and direction diagnostics

**STATUS: EXPLORATORY / POST HOC (Tier 4; prereg sec. 4.5 and Q25; Annex C.4). Hypothesis-generating only. Nothing here is a claim, nothing moves between tiers, and the confirmatory outcome is not re-decided.**

Frozen confirmatory result (not restated differently anywhere in this folder; see `mirror/results/confirmatory/report.md`): H1 final outcome **B**; beta_PS = +4.963 ppb (HC3 SE 3.125, p = .1169, delta = +0.44); bound test p_B = .0004; PS-only - TD = +8.157 ppb (p = .0443, mandatory co-report, K1 label none); S1/S2 not supported under Holm; S3 exploratory; T1 not evaluated; PC2x not confirmed (caveat on every H1 statement). First-usable n = 80 (TD 22 / ADHD-only 16 / PS-only 17 / PS+ADHD 25).

Code: `results_exploratory/code/tier4_exploratory.R` (run 2026-10-05T13:01:39Z; input SHA-256 `f669bc6f...` = the confirmatory input). The frozen helpers in `mirror/scripts/analysis/R/` are sourced read-only, so the model form (`y ~ PS_e * ADHD_e + ns(age, 3) + sex_e + race3 + scale(roi_vol) + scale(Q_phase)`, effect coding +/-0.5, race3 contr.sum over the levels present, OLS with HC3 and t(n - p)), the plausibility edits and the first-usable sample are identical to `run_all.R`. Reproduction checks: M1 beta_PS 4.9628 / SE 3.1251; K3 beta_dref 1.486 (p .2742); K3 CSF-referenced beta 3.477 (p .3507); K2 R2* beta 1.098; SN-only 5.874 and VTA+PBP 1.278 (frozen descriptives) all reproduce exactly.

Estimates are given with HC3 95% CIs; BH-FDR is applied within the registered X block and, separately, within the post hoc diagnostic block (`fdr_table.csv`). Positive beta_PS = higher chi in PS youth.

## Files

| file | content |
|---|---|
| `fdr_table.csv` | 32 tests of this block: estimate, CI, p, BH q within family (registered X / post hoc), BH q over the whole module |
| `group_effects_by_outcome.csv` | M1-form PS main effect, simple effects (PS-only - TD, PS+ADHD - ADHD-only), ADHD and interaction terms for every ROI, the X1 co-reports, the referencing and R2* rows |
| `raw_by_cell.csv` | n / mean / SD / median / quartiles / range by cell for every chi, R2*, reference, QC and age variable used |
| `adjusted_cell_means.csv` | M1-adjusted cell means (model-matrix means, proportional weights) with HC3 95% CI: SN-VTA chi, dentate R/L, RN, SN-VTA R2* |
| `qc_and_reference_by_PS.csv` | covariate-adjusted PS effect on each QC / nuisance / reference variable, Kruskal-Wallis over the 4 cells, Spearman correlation with chi_SNVTA |
| `age_slopes_by_cell.csv`, `X6_rsas_by_cell.csv` | per-cell OLS age slopes; RSAS by cell |
| `results_tier4.json` | everything above plus coefficient tables, X3 simple effects, X6 co-reports, specificity models, binary QC tables |
| `fig_X1_*.png`, `fig_X3_*.png`, `fig_X6_*.png`, `fig_D1..D7_*.png` | figures (each carries the EXPLORATORY banner) |
| `run.log`, `SHA256SUMS` | run log; checksums of every output and of the script |

## 0. Sample and availability (first-usable, n = 80)

- Dentate chi is present only where G5-DN passed: n = 58 (TD 13 / ADHD-only 12 / PS-only 16 / PS+ADHD 17). G5-DN failures are unevenly distributed: TD 9/22, ADHD-only 4/16, PS-only 1/17, PS+ADHD 8/25 (Fisher p = .079).
- RSAS at the same session: n = 55 (14 / 10 / 15 / 16); missing TD 8/22, ADHD-only 6/16, PS-only 2/17, PS+ADHD 9/25.
- The forced-review flag was joined from the coded QC sheet (`mirror/qc/qc_ids.csv` -> `rating_pass1.csv`) for all 80 sessions.

## 1. X1: dentate and red-nucleus group effects (registered M1 form)

| outcome | n | beta_PS ppb [95% CI] | p | in raw-SD units | PS-only - TD [95% CI], p | ADHD main, p | PS x ADHD, p |
|---|---|---|---|---|---|---|---|
| chi_DN_R (right dentate) | 58 | **-2.14 [-9.78, 5.50]** | .576 | -0.16 | -1.93 [-14.0, 10.1], .74 | -5.72, .12 | -0.41, .96 |
| chi_DN_L (left dentate) | 58 | **+1.68 [-5.87, 9.22]** | .657 | +0.13 | -1.21, .82 | -4.32, .28 | +5.77, .49 |
| chi_RN (red nucleus) | 80 | **+1.17 [-5.05, 7.39]** | .709 | +0.11 | +4.58, .22 | -0.19, .95 | -6.83, .25 |

BH q within the registered X block (5 tests) = .81 for all three. Co-reports (not in the FDR block): right dentate with its own volume (`scale(DN_R_vol)`) in place of `scale(roi_vol)`: -2.12 [-9.72, 5.49], p .578; bilateral dentate mean: -0.23 [-7.33, 6.87], p .948; SN-VTA chi with M1 refitted in the same G5-DN-pass subsample (n 58): +5.73 [-3.76, 15.21], p .23 (descriptive; same sign and size as the frozen full-sample estimate). Dentate and RN R2* show no PS effect either (DN_R -0.43 s^-1 [-2.34, 1.47], p .65; RN -0.38 [-1.44, 0.68], p .47). Adjusted cell means (ppb): DN_R TD 37.1 / ADHD-only 31.6 / PS-only 35.2 / PS+ADHD 29.3; RN 39.1 / 42.3 / 43.6 / 40.0.

**Reading:** no dentate or red-nucleus group effect; the CIs for the dentate span about +/-0.6 SD, so only large effects are excluded. The positive SN-VTA shift does not extend to the dentate or RN.

## 2. X3: sex x PS on SN-VTA chi (M1 + PS_e:sex_e)

- Interaction (PS effect in boys minus girls): **+1.83 ppb [-9.75, 13.40], p = .754** (q .81). With sex x ADHD also in the model: +0.90 [-10.55, 12.35], p .88.
- Simple PS effects: girls +3.98 [-5.58, 13.55] (p .41; n F = 38); boys +5.81 [-1.80, 13.42] (p .13; n M = 42). Sex by cell F/M: TD 9/13, ADHD-only 9/7, PS-only 9/8, PS+ADHD 11/14.

**Reading:** no evidence that the positive PS shift differs by sex; both simple effects are positive and imprecise.

## 3. X6: RSAS social anhedonia (cross-sectional) vs SN-VTA chi

Primary (S2 form with RSAS in place of PRIME; ADHD_e, ns(age, 3), sex_e, race3, roi_vol, Q_phase; n = 55): **-0.52 ppb per SD of RSAS [-4.88, 3.85], p = .813** (q .81; SD of RSAS 3.02; raw Spearman rho = -0.02, p .88).
Co-reports: + PS_e: -1.67 [-6.20, 2.86], p .46 (beta_PS in that n = 55 model: +9.89, p .064); + BDI + ARI (n 53): -0.94 [-5.96, 4.09], p .71; within PS only (n 31): +0.52 [-3.57, 4.60], p .79; rank-inverse-normal RSAS: -0.68 [-5.57, 4.21], p .78. All-session between/within-person split (68 usable sessions, 60 subjects, REML maximum `lmer_bobyqa`, ICC .52): between-person -0.35 (SE 1.78, p .84); within-person +6.0 (SE 8.0, containment df 4, p .50; uninformative).
RSAS by cell (mean, SD): TD 2.1 (2.4), ADHD-only 3.3 (2.8), PS-only 4.1 (2.5), PS+ADHD 4.7 (3.7).

**Reading:** no association between social anhedonia and SN-VTA chi in any specification; RSAS is higher in the clinical cells but does not track chi. The RSAS-complete subsample is not representative of the cells (missingness 9% in PS-only vs 36% in TD) and shows a larger PS shift than the full sample; this is a subsample feature, not a finding.

## 4. Direction diagnostics (descriptive; post hoc)

### 4a. SN-VTA chi by cell (`fig_D1_chi_SNVTA_by_cell.png`)

| cell | n | raw mean (SD) | median | M1-adjusted mean [95% CI] |
|---|---|---|---|---|
| TD | 22 | 50.1 (12.6) | 50.7 | 49.7 [44.1, 55.3] |
| ADHD-only | 16 | 50.9 (13.8) | 50.8 | 51.6 [43.5, 59.8] |
| PS-only | 17 | 56.8 (10.0) | 53.9 | 57.9 [52.3, 63.5] |
| PS+ADHD | 25 | 54.3 (10.8) | 52.1 | 53.4 [48.9, 57.9] |

The elevation sits in PS-only (raw +6.8, adjusted +8.2 vs TD) and is weak in PS+ADHD (+1.8 vs ADHD-only, frozen K1 co-report). TD has the widest spread with a low tail (three TD sessions near 30 ppb); the frozen Huber estimate (+3.94) is smaller than OLS, i.e. part of the mean difference comes from that tail. Medians differ less than means (PS-only 53.9 vs TD 50.7).

### 4b. Is the PS elevation SN-specific or global? (`fig_D2_*.png`; `group_effects_by_outcome.csv`)

| ROI | beta_PS ppb [95% CI] | p | raw-SD units | PS-only - TD, p |
|---|---|---|---|---|
| SN-VTA (primary, frozen) | +4.96 [-1.27, 11.20] | .117 | +0.42 | +8.16, .044 |
| SN only (SNc + SNr) | +5.87 [-0.73, 12.48] | .080 | +0.44 | +9.16, .036 |
| VTA + PBP | +1.28 [-4.63, 7.19] | .668 | +0.13 | +4.64, .19 |
| globus pallidus | +2.27 [-3.59, 8.13] | .441 | +0.18 | +3.45, .43 |
| subthalamic nucleus | +2.48 [-3.39, 8.34] | .403 | +0.22 | +6.20, .12 |
| red nucleus | +1.17 [-5.05, 7.39] | .709 | +0.11 | +4.58, .22 |
| putamen | -1.33 [-4.47, 1.82] | .402 | -0.22 | +0.52, .80 |
| caudate | -2.37 [-5.84, 1.11] | .180 | -0.34 | -0.52, .82 |
| striatum (CAU+PUT+NAcc) | -1.19 [-4.16, 1.79] | .429 | -0.20 | +0.98, .61 |
| thalamus | -0.45 [-2.20, 1.30] | .611 | -0.12 | +1.86, .17 |
| deep-GM composite (CAU+GP) | -1.45 [-4.27, 1.37] | .309 | -0.23 | -0.12, .95 |
| right / left dentate (n 58) | -2.14 / +1.68 | .58 / .66 | -0.16 / +0.13 | -1.93 / -1.21 |
| 6-ROI global z composite (GP, PUT, CAU, THAL, RN, STN) | -0.03 SD [-0.33, 0.28] | .851 | -0.03 | |

Global-iron adjustment: SN-VTA PS effect with the deep-GM composite as a covariate +5.54 [-0.57, 11.65] (p .075); with GP as a covariate +3.90 [-1.68, 9.48] (p .167; GP is the ROI most correlated with SN-VTA chi, r = .52, together with RN .69 and STN .53, so GP adjustment removes shared midbrain/pallidal variance, not only a global offset). None of the regional tests survives BH (min q .53 in the post hoc family).

**Reading:** the shift is confined to the SN proper (SNc+SNr); it is near zero in the VTA+PBP part of the same ROI, small and same-signed in the pallidum and STN, and zero or opposite-signed in striatum, thalamus, dentate and the global composite. A global scaling, offset or referencing artefact would move every ROI together; it does not. Incidental, uncorrected: thalamus shows a PS x ADHD interaction (-4.6 ppb, p .009) and PS+ADHD - ADHD-only -2.8 (p .013); SN-VTA R2* shows an ADHD main effect (-1.43 s^-1, p .007). Neither is part of any registered question.

### 4c. Referencing (`fig_D3_*.png`)

- Whole-brain reference value (ref_ppb, the value subtracted from every ROI): TD -0.59 (SD 0.96), ADHD-only -0.28 (0.20), PS-only -0.27 (0.32), PS+ADHD -0.39 (1.15) ppb; PS effect +0.03 ppb [-0.23, 0.29], p .82. The reference varies by < 1 ppb and cannot generate a 5-ppb group difference; SN-VTA chi before referencing (chi + ref_ppb) gives the identical beta_PS, +4.96 [-1.28, 11.20], p .117.
- CSF reference offset (CSF - whole-brain reference): TD 18.5 (6.9), ADHD-only 22.5 (6.6), PS-only 21.4 (5.6), PS+ADHD 23.7 (4.4) ppb; Kruskal-Wallis over 4 cells p .017; M1-form PS effect +1.49 [-1.20, 4.18], p .274 (this is the frozen K3 beta_dref, reproduced; K3 trigger FALSE, not re-decided); ADHD effect +2.93 (p .033, post hoc). Under the CSF reference the frozen co-report is +3.48 [-3.91, 10.86], p .351 (outcome B unchanged).
- Reference-domain composition: WM fraction PS -0.0065 [-0.014, 0.001], p .11; CSF fraction +0.0015, p .52; brain-in-FOV fraction +0.0014, p .22; CSF reference type (CSF / eroded-1) by cell p .83.

**Reading:** not a whole-brain referencing artefact. The two references disagree on the PS effect by about 1.5 ppb because the CSF offset is itself 3-5 ppb higher in the three clinical cells than in TD (ventricle partial volume, flow or CSF composition; v4 and ICV below). That disagreement lies inside the registered K3 tolerance, but it means the CSF-referenced shift is smaller than the whole-brain-referenced one and that the clinical cells differ from TD in the CSF reference itself, a mild warning sign rather than a clean bill.

### 4d. R2* (iron vs myelin direction; `fig_D4_*.png`)

- SN-VTA R2* by cell: TD 22.7 (2.0), ADHD-only 21.8 (1.9), PS-only 24.7 (2.5), PS+ADHD 22.6 (2.6) s^-1; adjusted 22.8 / 22.0 / 24.6 / 22.5.
- M1-form PS effect **+1.10 s^-1 [0.08, 2.12], two-sided p .035** (+0.45 raw SD; PS-only - TD +1.74, p .034). The frozen K2 test is one-sided for *lower* R2* (p .9825; "lower iron" wording correctly not used); the two-sided reading here is post hoc, q = .47 after BH, and is not a claim.
- chi-R2* correlation: all .54; TD .47, ADHD-only .39, PS-only .81, PS+ADHD .54. R2* of other ROIs: see `group_effects_by_outcome.csv` (GP, putamen, SN-only: no PS effect).

**Reading:** R2* moves in the same direction and by a similar standardised amount as chi in the same ROI. A paramagnetic (iron-like) source raises both chi and R2*; a myelin loss would raise chi but lower R2*; a phase-only artefact (unwrapping, background field, dipole inversion, referencing) would leave the magnitude-derived R2* untouched. This is the one diagnostic that points toward tissue rather than artefact, with the caveats that it is post hoc, nominally significant only, and that R2* is also sensitive to partial volume and registration (4e).

### 4e. QC and nuisance metrics by cell (`fig_D5_qc_by_cell.png`; `qc_and_reference_by_PS.csv`)

| metric | TD / ADHD-only / PS-only / PS+ADHD (mean) | adjusted PS effect [95% CI], p | raw-SD units | Spearman with chi_SNVTA (p) | frozen M1 sensitivity with this covariate |
|---|---|---|---|---|---|
| G5 Dice (min of four) | .804 / .814 / .822 / .809 | +0.0077 [-0.0049, 0.0203], .22 | +0.31 | **+.36 (.001)** | beta 3.75, p .23 |
| motion index (ghost-to-signal) | .0034 / .0035 / .0026 / .0019 | -0.0009 [-0.0030, 0.0013], .42 | -0.20 | -.21 (.06) | beta 4.23, p .16 |
| Q_phase | -3.02 / -3.05 / -3.03 / -3.04 | -0.007 [-0.090, 0.075], .86 | -0.04 | -.18 (.11) | in M1 |
| SN-VTA atlas volume (mm^3) | 641 / 627 / 605 / 635 | -17.0 [-44.8, 10.9], .23 | -0.28 | +.20 (.08) | in M1 |
| B0 tilt (deg) | 8.1 / 5.7 / 8.6 / 7.4 | +0.43 [-1.88, 2.74], .71 | +0.09 | -.11 (.31) | |
| head pitch (deg) | | +1.10 [-2.38, 4.58], .53 | +0.16 | +.02 (.89) | beta 5.08, p .11 |
| intracranial volume (ml) | 1578 / 1500 / 1455 / 1522 | -44.8 [-106.4, 16.8], .15 (KW p .009) | -0.33 | +.12 (.27) | beta 5.34, p .106 |
| 4th-ventricle volume (ml) | 1.9 / 1.7 / 1.6 / 1.8 | -0.08 [-0.39, 0.23], .60 | -0.15 | -.06 (.62) | |
| in-support share / slab coverage / implausible-voxel share | ~1 / ~1 / 0 | p .26 / .92 / (all zero) | | | |

Binary: visual QC grade 0/1 (grade 2 is never usable): TD 2/20, ADHD-only 0/16, PS-only 0/17, PS+ADHD 2/23 (Fisher p .55); forced review: 2 / 0 / 0 / 1 (p .54); dental code D: 1 / 0 / 0 / 1 (p 1.0); T1w artefact: 1 / 2 / 1 / 2 (p .84); early epoch: 5 / 2 / 0 / 2 (p .14); frequency cluster p .66.

**Reading:** no QC-grade, forced-review, dental, T1w, epoch or frequency imbalance. Two continuous channels run in the direction of the effect: registration Dice is slightly higher (+0.3 SD, n.s.) and motion slightly lower (-0.2 SD, n.s.) in PS, and both correlate with chi (better registration and less motion go with higher measured SN-VTA chi, plausibly through less partial-volume dilution of a bright nucleus). The frozen sensitivities that add these covariates reduce beta_PS to 3.75 (Dice) and 4.23 (motion); they do not remove it, but Dice alone accounts for roughly a quarter of the point estimate. ICV differs across cells (KW p .009; PS-only smallest) and the ICV-adjusted frozen estimate is 5.34.

### 4f. Age, sex, race by cell

- Age (y): TD 13.1 (2.2; 8.8-17.8), ADHD-only 12.8 (2.4; 9.1-16.9), PS-only 14.0 (1.7; 11.1-16.0), PS+ADHD 14.1 (1.7; 9.7-17.0). PS youths are 1.09 y older (Welch p .018); 4-cell ANOVA p .11. Sex p .76. race3 (Black / White / Other): TD 9/9/4, ADHD-only 10/6/0, PS-only 11/4/2, PS+ADHD 15/7/3 (Fisher p .51); race3 coefficients in M1 are -0.88 and +0.33 ppb (p .71 / .89); beta_PS without race 4.89 (frozen 4.893). Hispanic ethnicity 2 / 4 / 1 / 2. Threshold psychosis: 3 (all PS+ADHD).

### 4g. SN-VTA chi vs age by cell (`fig_D7_chi_vs_age_by_cell.png`)

Per-cell OLS slopes (ppb/y): TD +2.13 [-0.39, 4.66] (p .094); ADHD-only +0.49 [-2.80, 3.79]; PS-only +1.57 [-1.58, 4.72]; PS+ADHD +0.05 [-2.64, 2.73]. Pooled linear slope +0.88 [-0.79, 2.56] (p .30); PS x age (linear) +0.01 [-3.30, 3.33] ppb/y, p .994; in M1 + PS_e:age +0.91, p .60.

**Reading:** chi rises with age in TD at about 2 ppb/y, and PS youths are a year older, so part of the *raw* PS-TD difference could be developmental. The M1 spline removes it and the adjusted contrast (PS-only - TD +8.2) is not smaller than the raw one (+6.8); the PS-only age range (11-16) is narrower than TD (9-18), so the comparison at the TD extremes rests on the spline. No PS x age interaction.

## 5. BH-FDR table (`fdr_table.csv`)

32 tests. Registered X family (5: X1 right dentate, X1 left dentate, X1 red nucleus, X3 PS x sex, X6 RSAS): smallest p = .576, **all q = .81**. Post hoc diagnostic family (27: regional PS effects, global-composite and DGM-adjusted models, referencing, R2*, Dice, motion, Q_phase, roi_vol, QC grade, forced review, age, race, PS x age): smallest p = .018 (age by PS) and .035 (R2*), **smallest q = .47**. Over all 32 tests the smallest q is .56. **No test in this block survives BH-FDR at q < .05.**

## 6. What would, and would not, be consistent with "higher iron in PS youth" versus an artefact

The frozen result is outcome B: an adult-magnitude deficit is excluded, and the point estimate runs in the opposite (positive) direction without reaching alpha. The diagnostics above ask only whether that positive direction looks like tissue or like a measurement artefact. Nothing below is a claim.

**Observations that would be expected under a real, small SN-confined paramagnetic (iron-like) elevation, and are observed:**
1. Regional pattern: the shift is in the SN proper (+0.44 SD), absent in VTA+PBP (+0.13), small in GP/STN, zero or negative in striatum, thalamus, dentate, and the 6-ROI global composite is -0.03 SD. A referencing, scaling, phase-unwrapping or background-field artefact would move the whole deep grey matter together. It survives adjustment for the deep-GM composite (+5.5).
2. R2* in the same ROI shifts the same way (+1.10 s^-1, +0.45 SD, two-sided p .035, post hoc) and chi and R2* are correlated (.54; .81 in PS-only). Iron raises both; demyelination would raise chi and lower R2*; a phase-only artefact would not touch R2*.
3. The whole-brain reference is flat across cells (+0.03 ppb) and the unreferenced estimate is identical.
4. Visual QC grade, forced review, dental/metal, T1w artefact, epoch, frequency cluster, phase-noise index, B0 tilt, pitch, in-support share and slab coverage are balanced across cells.

**Observations that would be expected under an artefact or confound, and are partly observed (none individually decisive):**
1. Registration / partial volume: G5 Dice is 0.3 SD higher in PS (n.s.) and correlates with chi (rho .36, p .001); the Dice-adjusted frozen estimate falls to 3.75 (p .23). Motion runs the same way (lower in PS, negatively related to chi; motion-adjusted 4.23). Atlas SN-VTA volume is 0.28 SD smaller in PS-only. A bright nucleus measured with better registration, less motion and a smaller warped ROI reads higher; this could plausibly account for one to two ppb of the +5. R2* is also partial-volume sensitive, so point 2 above is not immune to this channel.
2. The CSF reference offset is 3-5 ppb higher in all three clinical cells than in TD (Kruskal-Wallis p .017; mostly an ADHD-related shift, p .033 post hoc), and the CSF-referenced PS effect is smaller (+3.5). The references do not disagree beyond the registered tolerance, but the reference region itself differs by group.
3. Distribution: TD has a wider spread with a low tail; the robust (Huber) estimate (+3.94, frozen) is smaller than OLS, and medians differ by about 3 ppb versus 7 for means. The effect is partly "some TD low" rather than "PS high".
4. Selection: G5-DN failures (TD 41% vs PS-only 6%) and RSAS missingness differ by cell, i.e. the cells that reach the first-usable sample are not image-quality-matched even though the subject-level QC-by-group tests (frozen) were not significant.
5. Age: PS youths are a year older and chi rises with age in TD; the spline adjustment handles it on the observed support, but the TD age range extends beyond the PS range.

**What would be needed before "higher iron" could be said, and is not available here:** a pre-specified test of the positive direction in an independent sample (or the held-out follow-ups, which are frozen as descriptive); a convergent R2* result under a pre-specified two-sided rule; registration-quality-matched or Dice-stratified replication; and an effect that survives both references. Within this block every CI includes zero except the post hoc two-sided R2* contrast, and nothing survives BH-FDR. The honest summary is: the positive direction is SN-specific and accompanied by a same-direction R2* shift, which is the pattern iron would give, but a partial-volume/registration contribution of similar order to the estimate's uncertainty cannot be excluded, and the registered reading (outcome B; no claim about a positive effect) stands unchanged.

## 7. Provenance and limits

- Everything here uses the first-usable sample (n = 80) of the primary branch (QSMxT, primary pipeline, SN-VTA), the same plausibility edits (two sessions edited, none in the H1 sample) and the frozen gates. No server access; no new imaging data; no installs.
- X6 primary form (S2 form with RSAS, without BDI/ARI) and the X1 volume covariate (M1 `roi_vol`, with the matched `DN_R_vol` as co-report) are implementation choices made after unlock and are declared as such.
- The post hoc tests in section 4 and the FDR family split (registered X vs post hoc) were defined by this task text, after the unlock; they are labelled post hoc throughout.
- Figures use one fixed colour per cell (TD blue, ADHD-only orange, PS-only aqua, PS+ADHD yellow) with the cell name and n on the axis, so identity never depends on colour alone.
