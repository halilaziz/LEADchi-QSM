# LEAD-chi -- T3a: registered Tier-3 estimation modules (E-dev, E-z, E-int, E-anchor)

**STATUS: EXPLORATORY / TIER 3 (ESTIMATION ONLY). No claim is made from anything in this folder.** Tier-3 results carry
no error control beyond CIs and the calibrated p of E-dev (prereg v1.0 FINAL, sec. 4.3 tier table; Q25); BH-FDR is applied
within each named block as a reading aid only. Nothing here re-decides or restates the frozen confirmatory outcome
(`mirror/results/confirmatory/results.json`: H1 final outcome **B**, beta_PS = +4.96 ppb, p = .117, bound test p_B = .0004;
K1 label none; S1/S2 not supported; S3 exploratory; T1 not evaluated; PC2x not confirmed). Every H1-related statement in the
registered report carries the caveat "developmental positive control not confirmed"; that caveat applies here as well.

- Code: `../code/t3a_tier3_estimation.R` (sources the frozen helpers `mirror/scripts/analysis/R/{constants,utils,design,plausibility,m1,t1}.R` unchanged)
- Input: `mirror/unlocked/merged_qsmxt_primary_SNVTA.csv` (182 rows, 113 subjects; SHA-256 in `tier3_estimation.json` -> `meta`)
- Samples (frozen rules): 107 usable sessions / 80 subjects (TD 31 / ADHD-only 19 / PS-only 22 / PS+ADHD 35 sessions); first-usable n = 80 (22 / 16 / 17 / 25); 26 subjects with two or more usable sessions
- Registered text implemented: Q19 "E-dev" and "E-z" paragraphs, the mixed-model rule (REML maximum, containment df), Q25 Tier 3, Q22 (BH within block)
- Registered draw counts used: 2,000 bootMer draws (PB-LRT), 1,000 realised-design simulations per saturation value; seed 20260924 + 300 (never the confirmatory streams)
- Not run in this module: E-ADHD (ADHD main effect on chi_STR; separate Tier-3 item), the E-dev attrition battery of Q14 (IPW-GEE, pattern-mixture, tipping point)

## Files

| File | Content |
|---|---|
| `edev_lme_estimates.csv` | E-dev fixed effects per ROI (PS_e:age_w, PS_e:age_bar_c, ADHD_e:age_w, age_w, PS_e, ADHD_e, PS_e:ADHD_e) at the verified REML maximum; containment df; 95% CI; BH q within the 9-test age-interaction block |
| `edev_pb_lrt.csv` | parametric-bootstrap LRT of PS_e:age_w (ML fits; lme4::bootMer, 2,000 draws) with Monte Carlo 95% interval; asymptotic chi-square(1) p for comparison |
| `edev_within_vs_between.csv` | Wald test PS_e:age_w - PS_e:age_bar_c (registered gate for a developmental reading: p >= .10) |
| `edev_calibration.csv` | realised-design null calibration: type I rate at .05 and the calibrated p of the observed t, per ROI and saturation b in {0, .2, .35, .5} / y |
| `edev_cell_slopes.csv` | cell-specific age slopes (ppb / y): within-person from the LME; cross-sectional linear from the first-usable HC3 model; BH q within each 4-test block |
| `edev_sensitivities.csv` | registered age-confounding sensitivities: common support 10-18 y; overlap weights; age matching (caliper 0.5 y) |
| `edev_spline_curves_firstusable.csv` | first-usable spline-by-cell display curves (adjusted at covariate means) |
| `edev_figure.png` | per ROI: (a) all usable sessions with within-person segments and the GAMM display (never inferential); (b) first-usable spline by cell |
| `ez_subject_z.csv` | E-z per first-usable subject: TD-referenced mu, sigma, z, centile; TD subjects scored leave-one-subject-out |
| `ez_cell_summary.csv` | E-z per cell: mean z with 95% CI, shares z <= -1 (registered), z < -1.96, z > 1.96, |z| > 1.96 with Wilson 95% CIs |
| `ez_cell_vs_TD.csv` | mean-z differences vs TD (Welch) with BH q |
| `ez_centile_chart.png` | TD-referenced centile chart (2.5 / 16 / 50 / 84 / 97.5), one panel per cell; threshold-psychosis anchors marked |
| `eint_cell_means.csv` | the four adjusted cell means of M1 (first-usable n = 80; HC3 95% CI; covariates at sample means) |
| `eint_frozen_K1_restated.csv` | the frozen K1 quantities (interaction, simple effects, main effects) copied verbatim from results.json |
| `eint_cell_means.png` | cell-means figure |
| `eanchor_z.csv` | the 6 threshold-psychosis youths, every MEGRE session, z against TD; non-usable sessions flagged and excluded from all summaries |
| `tier3_estimation.json` | everything above plus model diagnostics (optimiser table, profile check, variance components, knots, equality check against the frozen M1) |
| `t3a_log.txt`, `t3a_console.txt`, `sessionInfo.txt` | run log, console output, R session |

## Methods (as implemented)

### E-dev (PS x age; between/within-person decomposition; all usable sessions)
Registered model (Q19), per ROI y in {chi_SNVTA, chi_STR, chi_DN_R}:

```
y ~ PS_e*ADHD_e + ns(age_bar_c, 3) + age_w + PS_e:age_w + ADHD_e:age_w + PS_e:age_bar_c + sex_e + race3
    + scale(vol) + scale(Q_phase) + (1 | sub) + (0 + age_w | sub)        # = (1 + age_w || sub)
```
age_bar = subject mean age over usable sessions (centred at the subject-level mean); age_w = age - age_bar. vol = the ROI's own
volume (roi_vol, str_vol, DN_R_vol; the registered S1/S3 convention); race3 contr.sum over the levels present; scale() within
the analysis sample. Rules applied as registered:
- random slope singular (lme4::isSingular on the ML full, ML null and REML full fits) -> dropped from both models;
- estimates at the verified REML maximum: lmer (bobyqa, Nelder-Mead) and lme (nlminb, optim; pdDiag for the uncorrelated
  slope) scored on one lme4 REML deviance function, a grid profile (theta_int 0-3 step .05 x theta_slope 0-1 step .02 when the
  slope is kept; 0-3 step .01 otherwise) confirms the global maximum (refit from the grid argmin if better); the reported
  fit is an lmer fit restarted at the best theta; df by nlme containment (registered mixed-model decision);
- the parametric-bootstrap LRT of PS_e:age_w uses ML fits (full and null at their ML maxima), lme4::bootMer on the null fit
  with the full model refitted to each simulated response, 2,000 draws, p = (#>= + 1)/(N + 1), Monte Carlo 95% interval;
- a developmental reading is allowed only if PS_e:age_w and PS_e:age_bar_c do not differ (Wald, containment df of the
  within term; p >= .10); only the within-person term is then interpreted (Neuhaus & Kalbfleisch 1998; Sorensen et al. 2021);
- realised-design null calibration, run before the E-dev inference of each ROI: the actual subjects, sessions, ages,
  cells and covariates are kept; chi is simulated with NO PS, ADHD or PS x age effect from a common saturating-exponential
  age trajectory mu_b(age) = A * (1 - exp(-b (age - a0))) / (1 - exp(-b (a1 - a0))) (b = 0: linear), a0/a1 = the observed
  age range, A = the first-usable linear OLS age slope of that ROI x (a1 - a0) (so every b has the same total rise and
  differs only in curvature), plus the observed REML variance components (intercept, slope if kept, residual); the
  registered model is refitted (1,000 simulations per b) and the t of PS_e:age_w and PS_e:age_bar_c recorded. Reported:
  the empirical type I rate at alpha .05 (containment df) and the calibrated p = share of null |t| >= observed |t|;
- cell-specific within-person slopes = age_w + PS_e * b(PS_e:age_w) + ADHD_e * b(ADHD_e:age_w) at the cell codes
  (the model has no PS x ADHD x age term, so they are additive); cross-sectional first-usable slopes from
  `y ~ cellf + cellf:age + sex_e + race3 + scale(vol) + scale(Q_phase)` with HC3; the spline-by-cell display from
  `y ~ cellf * ns(age, 3) + covariates` (first-usable; predicted at covariate means, race averaged over its frequencies);
- GAMM display (mgcv::gam, s(age, by = cell, k = 4) + s(sub, bs = "re"), REML): display only, never inferential;
- age-confounding sensitivities: common support 10-18 y (sessions); overlap weights (subject-level propensity of PS from
  ns(first-usable age, 3) + sex_e; w = 1 - e for PS, e otherwise; weighted OLS of the fixed part with cluster-robust HC2
  SE, because lmer weights are precision weights); age matching (greedy 1:1 nearest neighbour on first-usable age,
  caliper 0.5 y, without replacement; the LME refitted on the matched subjects' usable sessions).

### E-z (TD-referenced deviation z of SN-VTA chi)
TD-only model `chi_SNVTA ~ ns(age, 3) + sex_e` on the 22 TD first-usable sessions; the two interior knots and the boundary
knots are fixed from all 22 TD so that the leave-one-subject-out refits (one per TD subject; the left-out subject is scored
against the fit without it, with that fit's residual SD) do not shift the basis. Non-TD subjects are scored against the
full TD model. z = (chi - mu_TD(age, sex)) / sigma_TD (sigma = residual SD; the column `z_pred` also carries the prediction
SE in the denominator, as a co-report). Centile = 100 * Phi(z). Per cell: mean z with a t 95% CI, the registered share
z <= -1 with a Wilson 95% CI, and the task's shares beyond |z| > 1.96 (also split by sign). Mean-z differences vs TD by
Welch t. Registered wording: this is **not** a normative model (22 reference youths).

### E-int (PS x ADHD and the four adjusted cell means)
M1 is refitted on the first-usable sample exactly as registered and checked against the frozen results.json (beta_PS, the
interaction and its HC3 SE agree to < 1e-8; the check halts the script otherwise). The adjusted cell means set
(PS_e, ADHD_e, PS_e:ADHD_e) to the cell codes and every other model-matrix column (ns(age) basis, sex_e, race3 contrasts,
scaled roi_vol and Q_phase) to its sample mean; HC3 SEs by the linear-combination rule; df = 68. The interaction, the two
simple effects and the main effects are restated verbatim from the frozen K1 block (they are not re-decided).

### E-anchor (threshold-psychosis youths)
The 6 youths with dx_psychosis = 1 (all PS+ADHD; registered caveat) are scored at every MEGRE session against the full TD
model. Only 3 of the 6 have a usable session (the other 3 failed the frozen QC: visual grade, G5); their values are listed
for completeness, flagged, and excluded from every summary. Summaries: first-usable z of the 3 (mean, 95% CI, share
z <= -1 with Wilson CI), all 6 usable sessions of the 3, and the PS+ADHD cell mean z with the anchors removed.

## Results (registered draw counts; run 2026-10-05T13:04Z, 3.2 min; R 4.4.1, lme4 1.1-38, nlme 3.1-164)

All estimates are ppb (slopes ppb / y) with 95% CIs. **Every interval below that matters includes zero; read them as nulls.**
Independent re-implementation (`../code/t3a_verify_statsmodels.py`, statsmodels / numpy): 11/11 OLS-based quantities (E-int
cell means and HC3 SEs, E-z mu / sigma / z / centiles / cell shares) agree with R to a maximum relative difference of 1.5e-13.
The LME quantities are checked by the registered convention (four optimisers on one deviance function, all within 1.3e-6 of
the best, grid profile confirms the global maximum for all three ROIs; no boundary fit).

### E-dev (3 ROIs; 107 sessions / 80 subjects / 26 returners for SN-VTA and STR; 78 / 62 / 16 for the right dentate)

Random slope: **kept** for SN-VTA (slope SD 2.13 ppb/y; intercept SD 8.87; residual 6.89; ICC .62) and STR (1.92; 5.61;
1.75; ICC .91); **singular -> dropped from both models** for the dentate (registered rule; intercept SD 10.91; residual 3.98).

| ROI | PS x age_w (within) [95% CI], p (containment df) | PB-LRT p [MC 95%] (2,000) | calibrated p across b = 0 / .2 / .35 / .5 | PS x age_bar (between) [95% CI], p | ADHD x age_w [95% CI], p | within - between [95% CI], p -> reading |
|---|---|---|---|---|---|---|
| SN-VTA | **+3.01 [-1.49, 7.51]**, p = .179 (df 22) | .173 [.156, .190] | .183 / .189 / .168 / .171 | -0.34 [-3.01, 2.32], p = .80 | -0.14 [-4.64, 4.37], p = .95 | 3.35 [-1.92, 8.63], p = .20 -> developmental reading allowed; the within-person term is the one to read |
| STR | +1.18 [-0.96, 3.32], p = .26 | .307 [.286, .327] | .282 / .260 / .295 / .289 | +0.13 [-1.31, 1.57], p = .86 | -1.89 [-4.02, 0.24], p = .080 | 1.05 [-1.56, 3.66], p = .41 -> allowed |
| dentate R | +0.69 [-2.66, 4.04], p = .66 (df 11) | .652 [.630, .672] | .652 / .632 / .618 / .626 | +0.19 [-3.25, 3.64], p = .91 | -1.15 [-4.35, 2.06], p = .45 | 0.49 [-4.48, 5.47], p = .83 -> allowed |

BH-FDR within the 9-test age-interaction block: smallest q = .72 (STR ADHD x age_w); all other q >= .79. PB-LRT block: q >= .46.
Realised-design null calibration (1,000 simulations per b; amplitude A = 14.0 ppb for SN-VTA, 5.0 for STR, 33.7 for the
dentate, i.e. the observed linear rise over the age range): the empirical type I rate of the PS x age_w test at alpha .05 is
.040-.056 for SN-VTA, .054-.075 for STR, .030-.051 for the dentate (containment df), so the containment-df test is close to
nominal under a saturating common trajectory on this design; the mean null t is within -0.29 and +0.03 of zero. The calibrated
p never differs from the containment-df p by more than .02.

Reading (estimation only): within the same youths, SN-VTA chi changed at +2.1 ppb/y [-1.1, 5.2] in PS youth versus
-0.9 ppb/y [-4.0, 2.1] in non-PS youth; the difference (+3.0 ppb/y) is imprecise and compatible with zero and with a
difference of either sign of up to ~7 ppb/y. No ROI shows a PS x age effect, within or between persons.

Cell-specific within-person slopes (ppb/y; LME; BH q within the 4-cell block):

| ROI | TD | ADHD-only | PS-only | PS+ADHD |
|---|---|---|---|---|
| SN-VTA | -0.88 [-4.21, 2.45], q .68 | -1.02 [-5.26, 3.22], q .68 | +2.13 [-1.99, 6.25], q .39 | +1.99 [-1.67, 5.66], q .39 |
| STR | +1.55 [0.03, 3.06], q .13 | -0.34 [-2.53, 1.84], q .75 | +2.73 [0.81, 4.64], q .049 | +0.84 [-0.68, 2.35], q .39 |
| dentate R | +3.34 [0.88, 5.79], q .049 | +2.19 [-0.97, 5.35], q .31 | +4.03 [1.08, 6.97], q .049 | +2.88 [-0.04, 5.80], q .13 |

Cross-sectional first-usable linear slopes (HC3) by cell are all compatible with zero for SN-VTA and STR (q >= .32); the
dentate rises with age in every cell (TD +3.80 [0.65, 6.95], q .23; the pooled within-person dentate slope is +3.11 [1.40, 4.82]
ppb/y, descriptive). The grand within-person slope of STR is +1.19 [0.20, 2.18] ppb/y (descriptive).

Age-confounding sensitivities (PS x age_w): SN-VTA common support 10-18 y +2.80 [-3.13, 8.74] (98 sessions / 75 subjects);
overlap weights +3.03 [-2.06, 8.11] (ESS 72 subjects); age-matched (28 pairs) +1.42 [-7.31, 10.15]. STR: +0.12 [-1.59, 1.83];
+1.81 [-1.52, 5.13]; +2.08 [-0.64, 4.81]. Dentate: -0.26 [-5.16, 4.65]; -0.02 [-3.72, 3.67]; -1.80 [-7.77, 4.17] (19 pairs).
None moves the estimate out of the registered interval. The six lme4 "boundary (singular)" console messages come from
individual bootstrap / calibration / sensitivity refits of the STR slope model and are expected.

Not run here (registered E-dev extras for a later task): the Q14 attrition battery (IPW-GEE, pattern-mixture returner x PS,
delta-adjusted tipping point).

### E-z (TD-referenced z of SN-VTA chi; first-usable; TD n = 22, knots 12.11 / 14.32, boundary 8.8 / 17.8, sigma_TD = 11.78 ppb)

| cell | n | mean z [95% CI] | p vs 0 (BH q) | z <= -1 (registered), Wilson 95% | z < -1.96 | z > 1.96 | abs z > 1.96, Wilson 95% |
|---|---|---|---|---|---|---|---|
| TD (LOSO) | 22 | +0.26 [-0.56, 1.09] | .51 (ref.) | 4 (18%) [7, 39] | 2 | 2 | 4 (18%) [7, 39] |
| ADHD-only | 16 | +0.14 [-0.58, 0.86] | .69 (.69) | 1 (6%) [1, 28] | 1 | 1 | 2 (13%) [3, 36] |
| PS-only | 17 | +0.45 [0.01, 0.90] | .047 (.14) | 1 (6%) [1, 27] | 0 | 1 | 1 (6%) [1, 27] |
| PS+ADHD | 25 | +0.21 [-0.18, 0.60] | .28 (.42) | 3 (12%) [4, 30] | 0 | 0 | 0 (0%) [0, 13] |

Mean z vs TD (Welch): ADHD-only -0.13 [-1.18, 0.93]; PS-only +0.19 [-0.73, 1.10]; PS+ADHD -0.05 [-0.95, 0.85]; all q = .90.
Reading: no psychosis-spectrum cell shows an excess of youths below the TD reference (the registered z <= -1 shares, 6-12%, are
at or below the TD leave-one-out share of 18%); no first-usable PS youth lies below the TD 2.5th centile. The TD leave-one-out
z are wider than N(0, 1) (SD 1.86), as expected with 22 reference youths and a 3-df spline; this is **not** a normative model.

### E-int (adjusted cell means of M1; first-usable n = 80; HC3 95% CI; frozen K1 quantities restated verbatim)

TD 49.7 [44.1, 55.3] (n 22) | ADHD-only 51.6 [43.5, 59.8] (n 16) | PS-only 57.9 [52.3, 63.5] (n 17) | PS+ADHD 53.4 [48.9, 57.9] (n 25).
Frozen K1 (not re-decided): PS x ADHD interaction -6.39 [-18.76, 5.98], p = .306; PS-only - TD +8.16 [0.21, 16.10], p = .044
(mandatory co-report; K1 label none); PS+ADHD - ADHD-only +1.77 [-7.78, 11.32], p = .713; ADHD main effect -1.27 [-7.47, 4.93].
The refit reproduces the frozen beta_PS, interaction and HC3 SE to 5e-15.

### E-anchor (6 threshold-psychosis youths, all PS+ADHD; registered caveat)

Only **3 of 6** have a usable MEGRE session (sub-20238, sub-20976, sub-21104); sub-20860, sub-20865 and sub-21141 failed the
frozen QC (visual grade; G5) at every session and are excluded from every summary (their values are listed, flagged, in
`eanchor_z.csv`). First-usable z: -0.15 (13.3 y, centile 44), +1.15 (14.3 y, centile 87), -1.39 (16.0 y, centile 8);
mean -0.13 [-3.28, 3.02]; 1 of 3 at z <= -1 (Wilson [6%, 79%]); none beyond |z| > 1.96. All 6 usable sessions of the 3:
mean z -0.36 [-1.73, 1.01]; the later sessions are -0.90 (sub-20238, 14.7 y), +1.12 (sub-20976, 15.9 y) and -1.98 (sub-21104,
17.95 y; beyond the TD boundary knot, extrapolated and flagged). PS+ADHD without the anchors: mean z +0.26 [-0.16, 0.67].
With three youths the anchor is uninformative in either direction.
