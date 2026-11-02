# VERIFY: post hoc SN–VTA recovery test (`results_posthoc/recovery`)

Verified 2026-10-07 by an independent Python re-implementation. The script was written from scratch with numpy, pandas, scipy and statsmodels; the R script in `code/` was neither run nor copied. Inputs were read only. Outputs were written only into this `verification/` folder. There was no server access.

## Verdict

- **results.json: 317 of 317 numbers OK** (relative tolerance 1e-6, or half a unit of the 8th decimal that the JSON writer printed).
- **README.md: 79 of 80 numbers OK at the printed rounding. 1 MISMATCH:** the +10 ppb row of the k table gives SD 0.059, but the recomputed value is 0.05850, which rounds to **0.058**.
- **recovery_long.csv:** all 900 rows reproduced (d, k and the off-target deltas agree to < 3e-14; first_usable, effect codes, age and cell are identical).
- **Raw inputs:** all six `raw/*.csv` match `raw/SHA256SUMS`.
- **Claims:** 14 issues found (C1–C14 below). Most are wording that over-states post hoc, non-significant or descriptive results. One is a factual error in the logic: C11, "B requires an estimate on the positive side".

## What was recomputed (files in this folder)

- `verify_recovery.py`: the full recomputation. It writes `recomputed.json`.
  - Inputs: the raw ROI tables merged with `key_plain.csv` on (code, ses_order). For each level, rows are kept only if status is ok in both the primary and the rec table. Then merged with `merged_qsmxt_primary_SNVTA.csv` on participant_session.
  - k = (χ_SNVTA_raw,rec − χ_SNVTA_raw,primary)/L.
  - Summaries: n, mean, SD (n−1), median, Q1, Q3 (type 7), min, max.
  - Pooled OLS d ~ L.
  - Per-session slopes through the origin: Σ(L·d)/Σ(L²), for sessions with ≥ 3 levels.
  - M1 design:
    - ns(age, df = 3) implemented as a truncated-power natural cubic spline. The knots are the type-7 quantiles at 1/3 and 2/3, with boundary knots at min and max. This spans the same column space as R `splines::ns`.
    - race3 uses sum-to-zero coding over the levels present.
    - roi_vol and Q_phase are z-scored within the sample.
  - Standard errors: HC3 computed by hand, cross-checked against statsmodels HC3.
  - Huber proposal-2 scale (k = 1.345): IRLS re-implemented from the published MASS::rlm algorithm.
    - LS start; initial scale 1.4826·median|r|.
    - The scale is updated at the start of each iteration.
    - Converged when the relative change in the residuals is ≤ 1e-4 (maxit 200). The run converged in 6–8 iterations for every fit.
    - Registered fallback if it does not converge: MAD of the OLS residuals.
  - Decision rule:
    - σ_B = min(σ_resid, σ_rob).
    - δ_B = 0.6·√0.938278557806515.
    - t_B = (β + δ_B·σ_B)/SE, upper tail, df = n − 12.
    - A = β < 0 and p < .05; B = p_B < .05; otherwise C.
  - Spearman: ρ is computed on average ranks. The p-value follows R `cor.test`:
    - AS 89 Edgeworth series when there are no ties (dice_min, motion_gsr, roi_vol);
    - the t approximation when there are ties (Q_phase, age).
    - All five p-values match to < 1e-8.
- `compare.py` → `results_json_check.tsv` (317 rows) and `long_check.json`.
- `check_readme.py` → `readme_check.tsv` (80 rows).

Re-run: `python3 verify_recovery.py && python3 compare.py && python3 check_readme.py`.

## Non-OK items

### Number mismatch

| # | Where | Printed | Recomputed | Fix |
|---|---|---|---|---|
| N1 | README k table, L = +10, SD | 0.059 | 0.0584965 | print **0.058** |

### Claim issues (post hoc wording and over-statement)

| # | README text | Issue | Suggested fix |
|---|---|---|---|
| C1 | "Recovery is about half, linear and **equal across groups**." | "Equal" is inferred from a non-significant test. | "Recovery is about half, approximately linear, and showed **no detectable group difference**." |
| C2 | "k **does not differ** by group (**M1-form model** at L = −14.24)" | (a) This is absence of evidence. (b) The model is M1 **without** z(roi_vol) and z(Q_phase) (df 70, not 68). | "No detectable group difference in k (M1 covariates without roi_vol/Q_phase; HC3, df 70): PS +0.0006 (95% CI −0.021 to 0.022, i.e. within about ±5% of k) …" |
| C3 | "The pipeline's sensitivity is **therefore the same in every cell**, and **no differential bias arises**." | Over-statement. The CI bounds a PS difference only to −4.4% to +4.7% of k. The test covers only a synthetic, uniform, phase-only perturbation. k also varies with age (ρ .26). | "A PS-related difference in sensitivity larger than about ±5% of k is unlikely **for this synthetic perturbation**. Whether this holds for a real iron difference, which also changes magnitude/R2*, is untested." |
| C4 | "k is **unrelated** to registration Dice (ρ .02), motion (−.09), ROI volume (−.01) and Q_phase (.04)." | "Unrelated" is inferred from non-significant ρ at n = 80. The 95% CI half-width for ρ is about ±0.22. | "k showed no clear association with Dice, motion, ROI volume or Q_phase (\|ρ\| ≤ .09; n = 80)." |
| C5 | "It **rises** slightly with age (ρ .26, p .02)." | This is an uncorrected p, one of five post hoc correlations. After Holm/Bonferroni correction p ≈ .10. It must not read as a finding. | "(ρ .26, uncorrected p .02; one of five post hoc correlations, p ≈ .10 after Holm correction; exploratory)." |
| C6 | "**Off-target effects are negligible.** Whole-brain reference shift: +0.006 ppb **on average** … GP: +0.18 ppb **on average** (SD 1.6)." | (a) No level is given. The +0.006 and +0.18 values are the L = −14.24 means only. Across levels the reference means are 0.0007–0.0073 ppb and the GP means −0.27 to +0.57 ppb. (b) The mean shifts are small, but the **per-session GP changes are not negligible**: SD 1.6–1.9 ppb, range about ±5–6 ppb, at **every** level including +10. Because this spread does not scale with L, it points to run-to-run variability of the re-reconstruction, which the README does not mention. (c) The thalamus shift scales with L, at about −1.2% of L at every level. | "At L = −14.24 the mean off-target shifts are small: reference +0.006, thalamus −0.17 (about −1.2% of L at every level), GP +0.18 ppb. Per-session GP changes vary (SD 1.6–1.9 ppb, range about ±6 ppb at every level including +10), which indicates some run-to-run variability of the re-reconstruction." |
| C7 | "The response is **linear and sign-symmetric**." | This holds for the level means (first-usable k 0.467–0.477), but per-session k is noisy at small \|L\|. Across all sessions k ranges 0.07–1.72 at −5. One session (sub-20494_ses-2, not first-usable) has k = 1.72 at −5 and a wrong-sign change at +10 (k = −0.19). | "approximately linear and sign-symmetric **on average**; per-session k is noisy at small \|L\| (all sessions: 0.07–1.72 at −5; one non-first-usable session has a wrong-sign change at +10)." |
| C8 | "The observed β_PS = +4.96 ppb … **excludes** such changes, one-sided" | "Excludes" is too strong for a post hoc test. The test treats k as known and fixed, and assumes a real deficit is recovered like the synthetic phase-only one. | "… is **incompatible at the one-sided p values below** with such changes (post hoc; k treated as known; assumes a real deficit is recovered like the synthetic one)." |
| C9 | Column header "One-sided p **that the data arise from it**" | This mis-defines a p-value: it is not the probability that the data arose from the hypothesis. | "One-sided p for H0: β_PS = 0.473·D (alternative β_PS > 0.473·D), t = (β − 0.473·D)/SE, df 68." |
| C10 | "−14.24 (**the adult raw difference**)" in the **true deficit** column; "If PS youth had had **the adult raw deficit** on top of what was observed, the study would have reported C" | This mixes scales. −14.24 is an *observed*-scale value from another pipeline with an unknown k, as the README's own "Other pipelines" caveat says. Here it is used as a *true* deficit in this pipeline. | "−14.24 (numerically equal to the adult raw difference, itself measured on another pipeline's observed scale)"; "If PS youth had carried an additional **true** deficit of 14.24 ppb …" |
| C11 | "B was therefore **not produced by an insensitive pipeline**: it **requires an estimate on the positive side**." | (a) **Factually inaccurate.** At the observed SE (3.125) and σ_B, B holds for any β_PS > **−0.92 ppb** (−6.132 + t.95,68·3.125). The −10 counterfactual gives B with β = +0.22, close to that threshold. B does not require a positive estimate. (b) "Not … insensitive" is over-stated: the pipeline recovers about 47% of the synthetic change. | "In this post hoc counterfactual, an additional true deficit of 10 ppb would still have given B, while 14.24 ppb or more would have given C. At the observed SE, B requires β_PS > about −0.9 ppb." |
| C12 | True-scale translation: "the registered B bound … **corresponds to** a true deficit of −13.0 ppb" etc.; Caveat "k is **not applied to any registered quantity**" | (a) The translation treats mean k as exact. With the 95% CI of mean k (0.461–0.484) the bound is −12.7 to −13.3 ppb. It also assumes the synthetic recovery applies to a real difference. (b) It contradicts the caveat, because the registered bound and CIs *are* divided by k here. | "would correspond, if k = 0.473 applied to a real group difference, to about −13 ppb (−12.7 to −13.3 over the 95% CI of mean k)". Caveat: "k does not alter any registered quantity; the true-scale numbers are a descriptive post hoc translation." |
| C13 | Heading "**What this means for outcome B.**" | This is interpretive framing of a post hoc analysis. | "Post hoc interpretation (descriptive; does not change the registered outcome)." |
| C14 | Caveat: "The **likely** shared causes are SMV/V-SHARP high-pass filtering and regularised inversion …" | This cause is not tested by any analysis here. | Replace "likely" with "possible" and add "(not tested)". |

### Statements that cannot be verified from this folder (not counted as mismatches)

- "An earlier test at +200 ppb in 3 sessions gave 0.50–0.57": there are no data for it here.
- The script hashes (`c6c389fd…`, `5216e7a3…`), QSMxT v9.22.0, the settings (ROMEO, V-SHARP 12 → 1.5 mm, MEDI λ = 0.002, Mask 3), zero-padding, and "3–4 voxels thick at 1.5 mm": these are on the server side.
- "B0 direction from the affine": `injection_log_rec_m14.json` gives `bdir = [0, 0, 1]` for all 180 injected sessions. That is consistent only with affines that carry no obliquity. The ROI tables give `b0_tilt_deg` 0.3–22° relative to the MNI frame, which is head orientation, not slab obliquity. Worth confirming against the NIfTI headers.
- Engine-2 k = 0.36 (RTS) and 0.25 (TV) agree with `../engine2/README.md` (lines 52–53). They were not recomputed here.
- "all 180 registered coded sessions" is consistent with the data: there are 181 coded sessions, and one (763f3b, ses_order 1) is `not_registered` / `no_registration` in every table and in the injection log.

## Key recomputed values (for reference)

| Quantity | Recomputed |
|---|---|
| Observed M1 β_PS (HC3 SE), df | +4.96280 (3.12508), 68; p_two .11691 |
| σ_resid / σ_rob (Huber prop. 2) / σ_B | 11.2559 / 10.5505 / 10.5505 |
| B bound; p_B; outcome | −6.13180; 0.000352; B |
| One-sided 95% lower limit; two-sided 95% lower limit | −0.24850; −1.27320 |
| Smallest β_PS still giving B at the observed SE | −0.92051 |
| Mean first-usable k at −14.24 (95% CI of mean) | 0.472826 (0.4614–0.4843) |
| True-scale bound / one-sided lower / two-sided lower | −12.9684 / −0.52555 / −2.69275 |
| Counterfactual β_PS (p_B; outcome) at −5, −10, −14.24, −20, +10 | +2.614 (.0033; B), +0.224 (.0228; B), −1.798 (.0845; C), −4.458 (.292; C), +9.779 (1.4e-6; B; p_two .0026, disconfirming) |
| D table one-sided p (README inputs / full precision) | −5: .010989 / .010989; −10: .0014032 / .0014031; −14.24: .00018838 / .00018836; −20: 9.006e-6 / 9.005e-6 |


## Full check: README.md numbers (printed rounding)

| Label | Printed | Recomputed | Rounded | Status |
|---|---|---|---|---|
| table k (FU) L=-5 mean | 0.470 | 0.4698109960249999 | 0.470 | OK |
| table k (FU) L=-5 SD | 0.088 | 0.08796548650263314 | 0.088 | OK |
| table k (FU) L=-5 median | 0.474 | 0.473542213 | 0.474 | OK |
| table k (FU) L=-5 min | 0.16 | 0.16090927200000066 | 0.16 | OK |
| table k (FU) L=-5 max | 0.96 | 0.9597381580000004 | 0.96 | OK |
| table k (FU) L=-10 mean | 0.467 | 0.4671466046500001 | 0.467 | OK |
| table k (FU) L=-10 SD | 0.059 | 0.05934108487701091 | 0.059 | OK |
| table k (FU) L=-10 median | 0.472 | 0.47216186500000035 | 0.472 | OK |
| table k (FU) L=-10 min | 0.30 | 0.30187988200000004 | 0.30 | OK |
| table k (FU) L=-10 max | 0.74 | 0.7446651460000006 | 0.74 | OK |
| table k (FU) L=-14.24 mean | 0.473 | 0.4728257849016854 | 0.473 | OK |
| table k (FU) L=-14.24 SD | 0.052 | 0.05155948954984628 | 0.052 | OK |
| table k (FU) L=-14.24 median | 0.473 | 0.4733352180477528 | 0.473 | OK |
| table k (FU) L=-14.24 min | 0.37 | 0.3651568061797757 | 0.37 | OK |
| table k (FU) L=-14.24 max | 0.68 | 0.6827884845505616 | 0.68 | OK |
| table k (FU) L=-20 mean | 0.469 | 0.46919192164374995 | 0.469 | OK |
| table k (FU) L=-20 SD | 0.047 | 0.04679650337295623 | 0.047 | OK |
| table k (FU) L=-20 median | 0.470 | 0.46954460125000014 | 0.470 | OK |
| table k (FU) L=-20 min | 0.37 | 0.3685279845 | 0.37 | OK |
| table k (FU) L=-20 max | 0.63 | 0.6319047925 | 0.63 | OK |
| table k (FU) L=+10 mean | 0.477 | 0.47676044584999994 | 0.477 | OK |
| table k (FU) L=+10 SD | 0.059 | 0.058496528087032415 | 0.058 | MISMATCH |
| table k (FU) L=+10 median | 0.480 | 0.47953300500000007 | 0.480 | OK |
| table k (FU) L=+10 min | 0.33 | 0.3307529449999997 | 0.33 | OK |
| table k (FU) L=+10 max | 0.70 | 0.6970264440000001 | 0.70 | OK |
| n per level (all sessions) | 180 | 180 | 180 | OK |
| pooled slope | 0.474 | 0.4740501879425916 | 0.474 | OK |
| pooled intercept | −0.05 | -0.04796552647098551 | -0.05 | OK |
| pooled R2 | 0.96 | 0.9616407395485311 | 0.96 | OK |
| pooled n | 900 | 900 | 900 | OK |
| k~group PS est | +0.0006 | 0.0006279350287484893 | 0.0006 | OK |
| k~group PS CI lo | −0.021 | -0.021016802868061842 | -0.021 | OK |
| k~group PS CI hi | 0.022 | 0.02227267292555882 | 0.022 | OK |
| k~group PS p | .95 | 0.9540245423763292 | .95 | OK |
| k~group ADHD est | +0.003 | 0.002898126832110881 | 0.003 | OK |
| k~group ADHD p | .79 | 0.7899999269172492 | .79 | OK |
| cell means min | 0.46 | 0.4624652854634832 | 0.46 | OK |
| cell means max | 0.49 | 0.4854990981179775 | 0.49 | OK |
| Spearman rho dice_min | .02 | 0.02440225035161744 | .02 | OK |
| Spearman rho motion_gsr | −.09 | -0.09045944678856073 | -.09 | OK |
| Spearman rho roi_vol | −.01 | -0.0074777308954524155 | -.01 | OK |
| Spearman rho Q_phase | .04 | 0.04032139914256531 | .04 | OK |
| Spearman rho age | .26 | 0.25823087473497497 | .26 | OK |
| Spearman p age | .02 | 0.02073857347731463 | .02 | OK |
| ref shift (L=-14.24 mean) | +0.006 | 0.006468531111111102 | 0.006 | OK |
| THAL shift L=-14.24 | −0.17 | -0.1736834839444444 | -0.17 | OK |
| GP shift (L=-14.24 mean) | +0.18 | 0.18202855866666645 | 0.18 | OK |
| GP SD (L=-14.24) | 1.6 | 1.5646818440390156 | 1.6 | OK |
| k multiplier | 0.47 | 0.4728257849016854 | 0.47 | OK |
| beta_PS observed | +4.96 | 4.962801544096132 | 4.96 | OK |
| SE observed | 3.125 | 3.125080404308326 | 3.125 | OK |
| df observed | 68 | 68 | 68 | OK |
| D=-5 expected shift | −2.36 | -2.364128924508427 | -2.36 | OK |
| D=-5 one-sided p | .011 | 0.010988601149731532 | .011 | OK |
| D=-10 expected shift | −4.73 | -4.728257849016854 | -4.73 | OK |
| D=-10 one-sided p | .0014 | 0.0014030708993288499 | .0014 | OK |
| D=-14.24 expected shift | −6.73 | -6.733039177 | -6.73 | OK |
| D=-14.24 one-sided p | .0002 | 0.00018836372620140602 | .0002 | OK |
| D=-20 expected shift | −9.46 | -9.456515698033709 | -9.46 | OK |
| D=-20 one-sided p | 9e-6 | 9.00453049802467e-06 | 9e-6 | OK |
| true-scale k | 0.473 | 0.4728257849016854 | 0.473 | OK |
| B bound observed | −6.13 | -6.131804643483686 | -6.13 | OK |
| B bound true | −13.0 | -12.968422703001849 | -13.0 | OK |
| one-sided 95% lower observed | −0.25 | -0.24849591338940158 | -0.25 | OK |
| one-sided 95% lower true | −0.53 | -0.52555491118377 | -0.53 | OK |
| two-sided 95% lower observed | −1.27 | -1.2731993109213464 | -1.27 | OK |
| two-sided 95% lower true | −2.7 | -2.6927450904270853 | -2.7 | OK |
| counterfactual L=-5.0 beta | +2.61 | 2.614261134208549 | 2.61 | OK |
| counterfactual L=-5.0 outcome | B | B | B | OK |
| counterfactual L=-10.0 beta | +0.22 | 0.22412651139658823 | 0.22 | OK |
| counterfactual L=-10.0 outcome | B | B | B | OK |
| counterfactual L=-14.24 beta | −1.80 | -1.7977552976907534 | -1.80 | OK |
| counterfactual L=-14.24 outcome | C | C | C | OK |
| counterfactual L=-20.0 beta | −4.46 | -4.458400154921755 | -4.46 | OK |
| counterfactual L=-20.0 outcome | C | C | C | OK |
| counterfactual L=10.0 beta | +9.78 | 9.77923006604423 | 9.78 | OK |
| counterfactual L=10.0 outcome | B | B | B | OK |
| counterfactual +10 p_two | .003 | 0.0026394940974218887 | .003 | OK |
| counterfactual +10 disconfirming | yes | True | True | OK |
| MEDI k (caveat) | 0.47 | 0.4728257849016854 | 0.47 | OK |

## Full check: results.json (rel. tol. 1e-6 or ≤ half a unit of the 8th printed decimal)

| Path | Reported | Recomputed | Status |
|---|---|---|---|
| `k_all[-20].n` | 180 | 180 | OK |
| `k_all[-20].mean` | 0.4756079 | 0.4756079031472221 | OK |
| `k_all[-20].sd` | 0.06105105 | 0.06105104828631487 | OK |
| `k_all[-20].median` | 0.47101612 | 0.4710161207500001 | OK |
| `k_all[-20].q25` | 0.44191747 | 0.4419174673750001 | OK |
| `k_all[-20].q75` | 0.50451183 | 0.5045118330000001 | OK |
| `k_all[-20].min` | 0.24780865 | 0.24780864750000015 | OK |
| `k_all[-20].max` | 0.77072899 | 0.7707289935 | OK |
| `k_all[-14.24].n` | 180 | 180 | OK |
| `k_all[-14.24].mean` | 0.47862644 | 0.4786264362905743 | OK |
| `k_all[-14.24].sd` | 0.07353823 | 0.07353823339310965 | OK |
| `k_all[-14.24].median` | 0.47260978 | 0.47260978265449427 | OK |
| `k_all[-14.24].q25` | 0.4403492 | 0.44034920365168506 | OK |
| `k_all[-14.24].q75` | 0.50900625 | 0.509006253511236 | OK |
| `k_all[-14.24].min` | 0.24197396 | 0.24197396278089878 | OK |
| `k_all[-14.24].max` | 0.90613794 | 0.9061379382022472 | OK |
| `k_all[-10].n` | 180 | 180 | OK |
| `k_all[-10].mean` | 0.48121187 | 0.48121187343888894 | OK |
| `k_all[-10].sd` | 0.08911458 | 0.08911458288701309 | OK |
| `k_all[-10].median` | 0.47622566 | 0.4762256619999999 | OK |
| `k_all[-10].q25` | 0.43488331 | 0.4348833085000004 | OK |
| `k_all[-10].q75` | 0.50781965 | 0.5078196522499999 | OK |
| `k_all[-10].min` | 0.22801514 | 0.22801513699999987 | OK |
| `k_all[-10].max` | 1.10015635 | 1.100156355 | OK |
| `k_all[-5].n` | 180 | 180 | OK |
| `k_all[-5].mean` | 0.47850164 | 0.4785016404555555 | OK |
| `k_all[-5].sd` | 0.14958036 | 0.14958035998871358 | OK |
| `k_all[-5].median` | 0.47354221 | 0.473542213 | OK |
| `k_all[-5].q25` | 0.43127289 | 0.4312728889999999 | OK |
| `k_all[-5].q75` | 0.5130661 | 0.5130661010000004 | OK |
| `k_all[-5].min` | 0.07457581 | 0.07457580600000142 | OK |
| `k_all[-5].max` | 1.72420368 | 1.7242036820000002 | OK |
| `k_all[10].n` | 180 | 180 | OK |
| `k_all[10].mean` | 0.46908684 | 0.4690868445166666 | OK |
| `k_all[10].sd` | 0.08619541 | 0.08619540695211465 | OK |
| `k_all[10].median` | 0.47352409 | 0.4735240935 | OK |
| `k_all[10].q25` | 0.43893261 | 0.4389326097499999 | OK |
| `k_all[10].q75` | 0.50568852 | 0.505688524 | OK |
| `k_all[10].min` | -0.18721752 | -0.18721752199999991 | OK |
| `k_all[10].max` | 0.69702644 | 0.6970264440000001 | OK |
| `k_first_usable[-20].n` | 80 | 80 | OK |
| `k_first_usable[-20].mean` | 0.46919192 | 0.46919192164374995 | OK |
| `k_first_usable[-20].sd` | 0.0467965 | 0.04679650337295623 | OK |
| `k_first_usable[-20].median` | 0.4695446 | 0.46954460125000014 | OK |
| `k_first_usable[-20].q25` | 0.43996792 | 0.4399679183750001 | OK |
| `k_first_usable[-20].q75` | 0.49304547 | 0.4930454730000002 | OK |
| `k_first_usable[-20].min` | 0.36852798 | 0.3685279845 | OK |
| `k_first_usable[-20].max` | 0.63190479 | 0.6319047925 | OK |
| `k_first_usable[-14.24].n` | 80 | 80 | OK |
| `k_first_usable[-14.24].mean` | 0.47282578 | 0.4728257849016854 | OK |
| `k_first_usable[-14.24].sd` | 0.05155949 | 0.05155948954984628 | OK |
| `k_first_usable[-14.24].median` | 0.47333522 | 0.4733352180477528 | OK |
| `k_first_usable[-14.24].q25` | 0.44444565 | 0.44444565045646023 | OK |
| `k_first_usable[-14.24].q75` | 0.49586745 | 0.49586744803370775 | OK |
| `k_first_usable[-14.24].min` | 0.36515681 | 0.3651568061797757 | OK |
| `k_first_usable[-14.24].max` | 0.68278848 | 0.6827884845505616 | OK |
| `k_first_usable[-10].n` | 80 | 80 | OK |
| `k_first_usable[-10].mean` | 0.4671466 | 0.4671466046500001 | OK |
| `k_first_usable[-10].sd` | 0.05934108 | 0.05934108487701091 | OK |
| `k_first_usable[-10].median` | 0.47216187 | 0.47216186500000035 | OK |
| `k_first_usable[-10].q25` | 0.43157644 | 0.4315764425000001 | OK |
| `k_first_usable[-10].q75` | 0.49225693 | 0.4922569274999997 | OK |
| `k_first_usable[-10].min` | 0.30187988 | 0.30187988200000004 | OK |
| `k_first_usable[-10].max` | 0.74466515 | 0.7446651460000006 | OK |
| `k_first_usable[-5].n` | 80 | 80 | OK |
| `k_first_usable[-5].mean` | 0.469811 | 0.4698109960249999 | OK |
| `k_first_usable[-5].sd` | 0.08796549 | 0.08796548650263314 | OK |
| `k_first_usable[-5].median` | 0.47354221 | 0.473542213 | OK |
| `k_first_usable[-5].q25` | 0.43125458 | 0.43125457849999976 | OK |
| `k_first_usable[-5].q75` | 0.50124531 | 0.5012453075000007 | OK |
| `k_first_usable[-5].min` | 0.16090927 | 0.16090927200000066 | OK |
| `k_first_usable[-5].max` | 0.95973816 | 0.9597381580000004 | OK |
| `k_first_usable[10].n` | 80 | 80 | OK |
| `k_first_usable[10].mean` | 0.47676045 | 0.47676044584999994 | OK |
| `k_first_usable[10].sd` | 0.05849653 | 0.058496528087032415 | OK |
| `k_first_usable[10].median` | 0.47953301 | 0.47953300500000007 | OK |
| `k_first_usable[10].q25` | 0.44683998 | 0.4468399757499998 | OK |
| `k_first_usable[10].q75` | 0.49789858 | 0.4978985784999999 | OK |
| `k_first_usable[10].min` | 0.33075294 | 0.3307529449999997 | OK |
| `k_first_usable[10].max` | 0.69702644 | 0.6970264440000001 | OK |
| `d_ref_ppb[-20].n` | 180 | 180 | OK |
| `d_ref_ppb[-20].mean` | 0.0073349 | 0.0073348998888888945 | OK |
| `d_ref_ppb[-20].sd` | 0.04759715 | 0.04759715283638647 | OK |
| `d_ref_ppb[-20].median` | 0.00555068 | 0.005550680000000011 | OK |
| `d_ref_ppb[-20].q25` | 0.00245222 | 0.002452224999999919 | OK |
| `d_ref_ppb[-20].q75` | 0.00798545 | 0.007985447499999965 | OK |
| `d_ref_ppb[-20].min` | -0.09400444 | -0.09400443999999997 | OK |
| `d_ref_ppb[-20].max` | 0.60406094 | 0.6040609399999999 | OK |
| `d_ref_ppb[-14.24].n` | 180 | 180 | OK |
| `d_ref_ppb[-14.24].mean` | 0.00646853 | 0.006468531111111102 | OK |
| `d_ref_ppb[-14.24].sd` | 0.04927577 | 0.049275770209666274 | OK |
| `d_ref_ppb[-14.24].median` | 0.00396148 | 0.003961484999999968 | OK |
| `d_ref_ppb[-14.24].q25` | 0.00194378 | 0.0019437825000000603 | OK |
| `d_ref_ppb[-14.24].q75` | 0.00599289 | 0.005992889999999973 | OK |
| `d_ref_ppb[-14.24].min` | -0.11198409 | -0.11198408999999998 | OK |
| `d_ref_ppb[-14.24].max` | 0.62257552 | 0.62257552 | OK |
| `d_ref_ppb[-10].n` | 180 | 180 | OK |
| `d_ref_ppb[-10].mean` | 0.00466479 | 0.004664793222222211 | OK |
| `d_ref_ppb[-10].sd` | 0.04868801 | 0.04868801313418568 | OK |
| `d_ref_ppb[-10].median` | 0.00245747 | 0.002457474999999999 | OK |
| `d_ref_ppb[-10].q25` | 0.00034228 | 0.000342275000000003 | OK |
| `d_ref_ppb[-10].q75` | 0.00405381 | 0.0040538100000000014 | OK |
| `d_ref_ppb[-10].min` | -0.11116686 | -0.11116685999999998 | OK |
| `d_ref_ppb[-10].max` | 0.61106123 | 0.61106123 | OK |
| `d_ref_ppb[-5].n` | 180 | 180 | OK |
| `d_ref_ppb[-5].mean` | 0.00066054 | 0.0006605395555555671 | OK |
| `d_ref_ppb[-5].sd` | 0.07059361 | 0.07059361478349231 | OK |
| `d_ref_ppb[-5].median` | 0.00136624 | 0.0013662400000000186 | OK |
| `d_ref_ppb[-5].q25` | 0.000333 | 0.00033300250000002293 | OK |
| `d_ref_ppb[-5].q75` | 0.00275898 | 0.0027589800000000775 | OK |
| `d_ref_ppb[-5].min` | -0.69250213 | -0.6925021300000003 | OK |
| `d_ref_ppb[-5].max` | 0.61057135 | 0.6105713500000001 | OK |
| `d_ref_ppb[10].n` | 180 | 180 | OK |
| `d_ref_ppb[10].mean` | 0.00129139 | 0.0012913926111111072 | OK |
| `d_ref_ppb[10].sd` | 0.04911149 | 0.04911148952852482 | OK |
| `d_ref_ppb[10].median` | -0.00307176 | -0.0030717600000000345 | OK |
| `d_ref_ppb[10].q25` | -0.00422474 | -0.004224739999999977 | OK |
| `d_ref_ppb[10].q75` | -0.00130862 | -0.0013086224999999951 | OK |
| `d_ref_ppb[10].min` | -0.09606779 | -0.09606779000000021 | OK |
| `d_ref_ppb[10].max` | 0.58453809 | 0.5845380899999999 | OK |
| `d_GP_ppb[-20].n` | 180 | 180 | OK |
| `d_GP_ppb[-20].mean` | 0.57409489 | 0.5740948893333331 | OK |
| `d_GP_ppb[-20].sd` | 1.92207306 | 1.9220730583113599 | OK |
| `d_GP_ppb[-20].median` | 0.22581291 | 0.22581291499999878 | OK |
| `d_GP_ppb[-20].q25` | 0.06374359 | 0.06374359249999983 | OK |
| `d_GP_ppb[-20].q75` | 0.47912789 | 0.4791278850000005 | OK |
| `d_GP_ppb[-20].min` | -6.14130402 | -6.141304020000007 | OK |
| `d_GP_ppb[-20].max` | 5.96845245 | 5.968452450000001 | OK |
| `d_GP_ppb[-14.24].n` | 180 | 180 | OK |
| `d_GP_ppb[-14.24].mean` | 0.18202856 | 0.18202855866666645 | OK |
| `d_GP_ppb[-14.24].sd` | 1.56468184 | 1.5646818440390156 | OK |
| `d_GP_ppb[-14.24].median` | 0.13751792 | 0.13751792499999382 | OK |
| `d_GP_ppb[-14.24].q25` | -0.02098465 | -0.02098465499999591 | OK |
| `d_GP_ppb[-14.24].q75` | 0.35124779 | 0.3512477899999986 | OK |
| `d_GP_ppb[-14.24].min` | -4.93043518 | -4.9304351800000035 | OK |
| `d_GP_ppb[-14.24].max` | 5.34202576 | 5.3420257599999985 | OK |
| `d_GP_ppb[-10].n` | 180 | 180 | OK |
| `d_GP_ppb[-10].mean` | 0.28219184 | 0.28219184311111095 | OK |
| `d_GP_ppb[-10].sd` | 1.88384457 | 1.8838445652505564 | OK |
| `d_GP_ppb[-10].median` | 0.11764145 | 0.1176414500000007 | OK |
| `d_GP_ppb[-10].q25` | -0.04449845 | -0.044498447499996985 | OK |
| `d_GP_ppb[-10].q75` | 0.31999969 | 0.3199996924999944 | OK |
| `d_GP_ppb[-10].min` | -6.27848053 | -6.278480529999996 | OK |
| `d_GP_ppb[-10].max` | 5.65403747 | 5.654037469999992 | OK |
| `d_GP_ppb[-5].n` | 180 | 180 | OK |
| `d_GP_ppb[-5].mean` | 0.15289375 | 0.15289375027777743 | OK |
| `d_GP_ppb[-5].sd` | 1.73380272 | 1.7338027215836647 | OK |
| `d_GP_ppb[-5].median` | 0.06156921 | 0.06156921499999868 | OK |
| `d_GP_ppb[-5].q25` | -0.12460136 | -0.12460135749999424 | OK |
| `d_GP_ppb[-5].q75` | 0.19900226 | 0.19900226249999875 | OK |
| `d_GP_ppb[-5].min` | -5.10383606 | -5.103836059999992 | OK |
| `d_GP_ppb[-5].max` | 5.5868988 | 5.5868988 | OK |
| `d_GP_ppb[10].n` | 180 | 180 | OK |
| `d_GP_ppb[10].mean` | -0.26700359 | -0.26700358661111145 | OK |
| `d_GP_ppb[10].sd` | 1.60341559 | 1.6034155931224108 | OK |
| `d_GP_ppb[10].median` | -0.10534668 | -0.10534667500000694 | OK |
| `d_GP_ppb[10].q25` | -0.25891113 | -0.2589111325000033 | OK |
| `d_GP_ppb[10].q75` | 0.04974365 | 0.0497436525000019 | OK |
| `d_GP_ppb[10].min` | -6.20172119 | -6.201721190000001 | OK |
| `d_GP_ppb[10].max` | 5.55840301 | 5.558403009999992 | OK |
| `d_THAL_ppb[-20].n` | 180 | 180 | OK |
| `d_THAL_ppb[-20].mean` | -0.26385334 | -0.26385334411111117 | OK |
| `d_THAL_ppb[-20].sd` | 0.27881551 | 0.27881550900538454 | OK |
| `d_THAL_ppb[-20].median` | -0.24962282 | -0.2496228200000017 | OK |
| `d_THAL_ppb[-20].q25` | -0.32442475 | -0.32442474500000085 | OK |
| `d_THAL_ppb[-20].q75` | -0.19908667 | -0.1990866674999987 | OK |
| `d_THAL_ppb[-20].min` | -1.54936742 | -1.5493674200000003 | OK |
| `d_THAL_ppb[-20].max` | 2.4234848 | 2.423484800000004 | OK |
| `d_THAL_ppb[-14.24].n` | 180 | 180 | OK |
| `d_THAL_ppb[-14.24].mean` | -0.17368348 | -0.1736834839444444 | OK |
| `d_THAL_ppb[-14.24].sd` | 0.26302695 | 0.263026953391322 | OK |
| `d_THAL_ppb[-14.24].median` | -0.17607117 | -0.17607117000000105 | OK |
| `d_THAL_ppb[-14.24].q25` | -0.21903204 | -0.21903204499999962 | OK |
| `d_THAL_ppb[-14.24].q75` | -0.13412237 | -0.13412237249999825 | OK |
| `d_THAL_ppb[-14.24].min` | -1.92744446 | -1.9274444600000038 | OK |
| `d_THAL_ppb[-14.24].max` | 2.45333481 | 2.453334810000001 | OK |
| `d_THAL_ppb[-10].n` | 180 | 180 | OK |
| `d_THAL_ppb[-10].mean` | -0.12480731 | -0.124807309611111 | OK |
| `d_THAL_ppb[-10].sd` | 0.2598199 | 0.2598199009241888 | OK |
| `d_THAL_ppb[-10].median` | -0.12285709 | -0.12285708999999834 | OK |
| `d_THAL_ppb[-10].q25` | -0.16115236 | -0.16115235750000023 | OK |
| `d_THAL_ppb[-10].q75` | -0.0901289 | -0.09012889749999964 | OK |
| `d_THAL_ppb[-10].min` | -1.40198564 | -1.4019856400000004 | OK |
| `d_THAL_ppb[-10].max` | 2.51765061 | 2.517650610000004 | OK |
| `d_THAL_ppb[-5].n` | 180 | 180 | OK |
| `d_THAL_ppb[-5].mean` | -0.06821045 | -0.06821044927777772 | OK |
| `d_THAL_ppb[-5].sd` | 0.27750534 | 0.2775053440030038 | OK |
| `d_THAL_ppb[-5].median` | -0.06050777 | -0.0605077749999996 | OK |
| `d_THAL_ppb[-5].q25` | -0.08044934 | -0.08044934249999969 | OK |
| `d_THAL_ppb[-5].q75` | -0.041852 | -0.04185200250000065 | OK |
| `d_THAL_ppb[-5].min` | -1.87934733 | -1.8793473300000008 | OK |
| `d_THAL_ppb[-5].max` | 2.63980484 | 2.6398048400000036 | OK |
| `d_THAL_ppb[10].n` | 180 | 180 | OK |
| `d_THAL_ppb[10].mean` | 0.11822689 | 0.11822688894444441 | OK |
| `d_THAL_ppb[10].sd` | 0.26062898 | 0.2606289806260253 | OK |
| `d_THAL_ppb[10].median` | 0.12647247 | 0.12647247499999814 | OK |
| `d_THAL_ppb[10].q25` | 0.08972216 | 0.08972215750000045 | OK |
| `d_THAL_ppb[10].q75` | 0.1523447 | 0.15234469749999846 | OK |
| `d_THAL_ppb[10].min` | -0.9020307 | -0.9020307000000001 | OK |
| `d_THAL_ppb[10].max` | 2.93006516 | 2.930065160000005 | OK |
| `k_linear_per_session.n` | 180 | 180 | OK |
| `k_linear_per_session.mean` | 0.47632395 | 0.4763239475354115 | OK |
| `k_linear_per_session.sd` | 0.06171934 | 0.06171934301944306 | OK |
| `k_linear_per_session.median` | 0.47219184 | 0.47219183943440846 | OK |
| `k_linear_per_session.q25` | 0.44189391 | 0.4418939136131494 | OK |
| `k_linear_per_session.q75` | 0.5056066 | 0.5056065979681017 | OK |
| `k_linear_per_session.min` | 0.24643521 | 0.24643520546484954 | OK |
| `k_linear_per_session.max` | 0.75676734 | 0.7567673359699515 | OK |
| `pooled_linear.slope` | 0.47405019 | 0.4740501879425916 | OK |
| `pooled_linear.intercept` | -0.04796553 | -0.04796552647098551 | OK |
| `pooled_linear.r2` | 0.96164074 | 0.9616407395485311 | OK |
| `pooled_linear.n` | 900 | 900 | OK |
| `k_by_group_L14.PS.est` | 0.00062794 | 0.0006279350287484893 | OK |
| `k_by_group_L14.PS.se` | 0.01085255 | 0.010852554722865358 | OK |
| `k_by_group_L14.PS.t` | 0.05786057 | 0.05786057244433761 | OK |
| `k_by_group_L14.PS.df` | 70 | 70 | OK |
| `k_by_group_L14.PS.p_two` | 0.95402454 | 0.9540245423763292 | OK |
| `k_by_group_L14.PS.p_lower` | 0.52298773 | 0.5229877288118354 | OK |
| `k_by_group_L14.PS.p_upper` | 0.47701227 | 0.4770122711881646 | OK |
| `k_by_group_L14.PS.ci95[0]` | -0.0210168 | -0.021016802868061842 | OK |
| `k_by_group_L14.PS.ci95[1]` | 0.02227267 | 0.02227267292555882 | OK |
| `k_by_group_L14.PS.ci90[0]` | -0.01746235 | -0.017462345573542877 | OK |
| `k_by_group_L14.PS.ci90[1]` | 0.01871822 | 0.018718215631039856 | OK |
| `k_by_group_L14.ADHD.est` | 0.00289813 | 0.002898126832110881 | OK |
| `k_by_group_L14.ADHD.se` | 0.01084094 | 0.010840943765710438 | OK |
| `k_by_group_L14.ADHD.t` | 0.2673316 | 0.26733159905114207 | OK |
| `k_by_group_L14.ADHD.df` | 70 | 70 | OK |
| `k_by_group_L14.ADHD.p_two` | 0.78999993 | 0.7899999269172492 | OK |
| `k_by_group_L14.ADHD.p_lower` | 0.60500004 | 0.6050000365413754 | OK |
| `k_by_group_L14.ADHD.p_upper` | 0.39499996 | 0.3949999634586246 | OK |
| `k_by_group_L14.ADHD.ci95[0]` | -0.01872345 | -0.018723453740846493 | OK |
| `k_by_group_L14.ADHD.ci95[1]` | 0.02451971 | 0.024519707405068255 | OK |
| `k_by_group_L14.ADHD.ci90[0]` | -0.0151728 | -0.01517279929758325 | OK |
| `k_by_group_L14.ADHD.ci90[1]` | 0.02096905 | 0.020969052961805013 | OK |
| `k_by_group_L14.cell_means[ADHD-only]` | 0.46246529 | 0.4624652854634832 | OK |
| `k_by_group_L14.cell_means[PS+ADHD]` | 0.4854991 | 0.4854990981179775 | OK |
| `k_by_group_L14.cell_means[PS-only]` | 0.4705156 | 0.4705155990994712 | OK |
| `k_by_group_L14.cell_means[TD]` | 0.46774434 | 0.46774434486721134 | OK |
| `k_spearman_L14.rho[dice_min]` | 0.02440225 | 0.02440225035161744 | OK |
| `k_spearman_L14.p[dice_min]` | 0.82957592 | 0.8295759165319683 | OK |
| `k_spearman_L14.rho[motion_gsr]` | -0.09045945 | -0.09045944678856073 | OK |
| `k_spearman_L14.p[motion_gsr]` | 0.42409328 | 0.4240932757972183 | OK |
| `k_spearman_L14.rho[roi_vol]` | -0.00747773 | -0.0074777308954524155 | OK |
| `k_spearman_L14.p[roi_vol]` | 0.94747289 | 0.947472894571305 | OK |
| `k_spearman_L14.rho[Q_phase]` | 0.0403214 | 0.04032139914256531 | OK |
| `k_spearman_L14.p[Q_phase]` | 0.72250419 | 0.7225041890496089 | OK |
| `k_spearman_L14.rho[age]` | 0.25823087 | 0.25823087473497497 | OK |
| `k_spearman_L14.p[age]` | 0.02073857 | 0.02073857347731463 | OK |
| `observed.beta_ppb` | 4.96280154 | 4.962801544096132 | OK |
| `observed.se_hc3` | 3.1250804 | 3.125080404308326 | OK |
| `observed.p_two` | 0.11691305 | 0.11691305305954833 | OK |
| `observed.ci95_ppb[0]` | -1.27319931 | -1.2731993109213464 | OK |
| `observed.ci95_ppb[1]` | 11.1988024 | 11.19880239911361 | OK |
| `observed.sigma_B` | 10.55045438 | 10.55045437585713 | OK |
| `observed.outcome` | B | B | OK |
| `observed.bound_ppb` | -6.13180464 | -6.131804643483686 | OK |
| `observed.lower95_one_sided` | -0.24849591 | -0.24849591338940158 | OK |
| `counterfactual[-5].L_true_ppb` | -5 | -5.0 | OK |
| `counterfactual[-5].n_missing_rec` | 0 | 0 | OK |
| `counterfactual[-5].mean_recovered_in_PS` | -2.38826098 | -2.3882609771428567 | OK |
| `counterfactual[-5].beta_ppb` | 2.61426113 | 2.614261134208549 | OK |
| `counterfactual[-5].se` | 3.12269572 | 3.1226957207239403 | OK |
| `counterfactual[-5].p_two` | 0.4054236 | 0.4054235992096764 | OK |
| `counterfactual[-5].pB` | 0.00334001 | 0.0033400140534950605 | OK |
| `counterfactual[-5].bound_ppb` | -6.12357071 | -6.123570707988813 | OK |
| `counterfactual[-5].outcome` | B | B | OK |
| `counterfactual[-10].L_true_ppb` | -10 | -10.0 | OK |
| `counterfactual[-10].n_missing_rec` | 0 | 0 | OK |
| `counterfactual[-10].mean_recovered_in_PS` | -4.7769383 | -4.776938302142857 | OK |
| `counterfactual[-10].beta_ppb` | 0.22412651 | 0.22412651139658823 | OK |
| `counterfactual[-10].se` | 3.12657547 | 3.1265754661349026 | OK |
| `counterfactual[-10].p_two` | 0.94306354 | 0.9430635367940028 | OK |
| `counterfactual[-10].pB` | 0.02275816 | 0.022758156827920078 | OK |
| `counterfactual[-10].bound_ppb` | -6.14558401 | -6.145584013943967 | OK |
| `counterfactual[-10].outcome` | B | B | OK |
| `counterfactual[-14.24].L_true_ppb` | -14.24 | -14.24 | OK |
| `counterfactual[-14.24].n_missing_rec` | 0 | 0 | OK |
| `counterfactual[-14.24].mean_recovered_in_PS` | -6.82714512 | -6.827145122857143 | OK |
| `counterfactual[-14.24].beta_ppb` | -1.7977553 | -1.7977552976907534 | OK |
| `counterfactual[-14.24].se` | 3.12853402 | 3.1285340166782882 | OK |
| `counterfactual[-14.24].p_two` | 0.56743658 | 0.5674365838305591 | OK |
| `counterfactual[-14.24].pB` | 0.08449828 | 0.08449827696388883 | OK |
| `counterfactual[-14.24].bound_ppb` | -6.14710016 | -6.1471001600800985 | OK |
| `counterfactual[-14.24].outcome` | C | C | OK |
| `counterfactual[-20].L_true_ppb` | -20 | -20.0 | OK |
| `counterfactual[-20].n_missing_rec` | 0 | 0 | OK |
| `counterfactual[-20].mean_recovered_in_PS` | -9.50998465 | -9.509984651904764 | OK |
| `counterfactual[-20].beta_ppb` | -4.45840015 | -4.458400154921755 | OK |
| `counterfactual[-20].se` | 3.13289023 | 3.132890227985163 | OK |
| `counterfactual[-20].p_two` | 0.15928049 | 0.15928049217744147 | OK |
| `counterfactual[-20].pB` | 0.29199846 | 0.29199846389558787 | OK |
| `counterfactual[-20].bound_ppb` | -6.18206181 | -6.182061810250923 | OK |
| `counterfactual[-20].outcome` | C | C | OK |
| `counterfactual[10].L_true_ppb` | 10 | 10.0 | OK |
| `counterfactual[10].n_missing_rec` | 0 | 0 | OK |
| `counterfactual[10].mean_recovered_in_PS` | 4.85292934 | 4.852929342380952 | OK |
| `counterfactual[10].beta_ppb` | 9.77923007 | 9.77923006604423 | OK |
| `counterfactual[10].se` | 3.13261659 | 3.132616593873969 | OK |
| `counterfactual[10].p_two` | 0.00263949 | 0.0026394940974218887 | OK |
| `counterfactual[10].pB` | 1.44414364e-06 | 1.4441436448073293e-06 | OK |
| `counterfactual[10].bound_ppb` | -6.20634908 | -6.206349082848962 | OK |
| `counterfactual[10].outcome` | B | B | OK |
| `true_scale.k_mean_L14` | 0.47282578 | 0.4728257849016854 | OK |
| `true_scale.bound_true_ppb` | -12.9684227 | -12.968422703001849 | OK |
| `true_scale.lower95_one_sided_true_ppb` | -0.52555491 | -0.52555491118377 | OK |
| `true_scale.lower95_two_sided_true_ppb` | -2.69274509 | -2.6927450904270853 | OK |
| `meta.rL` | 0.93827856 | 0.938278557806515 | OK |
| `meta.levels[0].variant` | rec_m05 | rec_m05 | OK |
| `meta.levels[0].L` | -5 | -5 | OK |
| `meta.levels[1].variant` | rec_m10 | rec_m10 | OK |
| `meta.levels[1].L` | -10 | -10 | OK |
| `meta.levels[2].variant` | rec_m14 | rec_m14 | OK |
| `meta.levels[2].L` | -14.24 | -14.24 | OK |
| `meta.levels[3].variant` | rec_m20 | rec_m20 | OK |
| `meta.levels[3].L` | -20 | -20 | OK |
| `meta.levels[4].variant` | rec_p10 | rec_p10 | OK |
| `meta.levels[4].L` | 10 | 10 | OK |
