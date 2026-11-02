# VERIFY: independent check of `results_posthoc/validation_local` (A–D)

**STATUS: POST HOC verification of POST HOC / registered-exploratory analyses.** Nothing in this file is a finding or a claim. The confirmatory outcome is not re-decided. It remains H1 final outcome **B**: beta_PS = +4.963 ppb, HC3 SE 3.125, n = 80. This check covers README.md (SHA-256 `566e0dda…0d9`), results.json, the 12 CSVs and `code/`.

Verifier: an independent re-implementation in **Python** (numpy/scipy); the code under check is R. Run on 2026-10-07. Nothing was written outside `verification/`. No server was contacted. No e-mail was sent and nothing was published.

---

## 1. Verdict

| item | result |
|---|---|
| A (X2), B (nuisance, Dice), C (low-TD tail) and every OLS row of D | **Reproduced.** 2,454 A–C cells plus all 54 first-usable D rows were recomputed from the raw merged table. The largest relative difference is 6e-12; no cell differs by more than 1e-6. |
| D mixed-model (lme) rows | **354 cells differ by more than 1e-6 relative**: 325 cells in 64 of the 90 lme rows of D_age_slopes, plus 29 cells of D_implied that are derived from them. The largest slope difference is 0.00022 ppb/y (0.18 %). Cause: `lmeControl(opt = "optim")` stops slightly before the REML optimum. nlme's default optimiser agrees with this verification (sec. 3.1). No conclusion changes. |
| Frozen anchors | M1 beta/SE/p, sigma_resid, Huber sigma_rob, the outcome-B bound −6.1318, t_B and p_B, the Huber beta, K2 and the PC2x slopes were all reproduced to ≤ 6e-15 relative. The three frozen single-covariate sensitivity rows (Dice 3.749 / motion 4.227 / ICV 5.337) reproduce to ≤ 1.5e-15 relative, and the frozen p_B values .0019 / .0004 / .0003 to ≤ 2.6e-14. |
| Dimov constants | **Confirmed** against the paper (sec. 2): 274 Hz/ppm for both components, theoretical 321 Hz/ppm, no R2 term in Eq. 1. The page number "p. 1545" could not be verified. |
| README printed numbers | **1,359 numbers checked.** 1,349 pass, 10 fail: 5 come from the lme optimiser issue, 5 are rounding errors (sec. 3.2). |
| Builder summary (workflow text) | Two numeric slips: "4.485 vs 4.486" (both values round to 4.485), and a truncated "Heterogeneity across tertiles: p ." (the value is .095). Several sentences over-claim (sec. 4). |
| Wording | **13 places over-claim** or read more strongly than the evidence allows (sec. 4). None calls anything a "finding", but several make causal or mechanistic statements and use iron language. |
| Outcome B | Unaffected. Every pooled post hoc specification has p_B ≤ .0022. |

---

## 2. How it was verified

- **Fresh code.** `verify_independent.py` and `vlib.py` are written from scratch. They do not source the frozen R helpers or `validation_local.R`. They read only the merged table, `gates_preunlock.json`, the frozen `results_confirmatory/results.json`, and (for comparison only) the outputs under check.
- **Sample.** The earliest usable session per participant is used. It equals the frozen `first_usable` flag: 80 sessions, cells 22/16/17/25. The usable sample is 107 sessions from 80 subjects, of whom 26 have 2 or more. The plausibility envelopes were re-applied: they change nothing in usable rows.
- **M1 re-implemented from its definition.**
  - Natural-spline age basis in truncated-power form (ESL eq. 5.4–5.5). The knots are the range and the type-7 1/3 and 2/3 quantiles. With the intercept, this basis spans exactly the column space of R's `ns(age, df = 3)`.
  - race3 dummies span the same space as contr.sum.
  - OLS is fitted by QR, with an HC3 sandwich and t(df = n − p).
  - Because the column spaces are identical, the PS_e estimates and SEs must agree to machine precision even though the basis columns differ, and they do.
- **Huber scale.** A line-by-line Python port of `MASS::rlm.default` (Huber psi, k = 1.345, "proposal 2" scale, acc 1e-4). It reproduces the frozen sigma_rob 10.55045437586 to 1.2e-15.
- **Mixed models.** Own profiled-REML random-intercept model with nlme containment df (26 for age; 74 for between-person age). It was cross-checked with statsmodels MixedLM (Powell): the age estimates agree to ≤ 5e-6.
- **Exact tests.**
  - Fisher r×c by full enumeration (R tolerance 1 + 1e-7).
  - Wilcoxon: exact when there are no ties; otherwise a normal approximation with continuity and tie correction, as in R.
  - Spearman (t approximation), BH q values, and the HC3 Wald test of tertile heterogeneity.
- **Printed numbers.** `check_readme.py` parses every number in the README tables, plus about 70 prose statements (61 snippet checks and 9 special-case checks). A value passes if it equals the recomputed value rounded to the printed precision.
- **Diagnostic only.** `diag_nlme_optimizer.R` adjudicates the lme discrepancy (sec. 3.1).

Outputs in this folder:

| file | content |
|---|---|
| `mine/*.csv`, `mine/scalars.json` | the full recomputation, with the same layout as the CSVs under check, plus `D_putamen_shortfall_by_Dr.csv` |
| `compare_csv_summary.csv`, `compare_csv_discrepancies.csv` | cell-by-cell comparison; every cell > 1e-6 relative is listed |
| `compare_results_json_scalars.csv` | scalars in results.json (VIF, cut points, Fisher, heterogeneity, Spearman, correlations) |
| `compare_readme.csv`, `compare_readme_failures.csv` | the check of every printed number |
| `verify_independent.log`, `diag_nlme_optimizer.out` | run logs |

---

## 3. Discrepancies

### 3.1 D: mixed-model rows (all usable sessions, `nlme::lme`)

`validation_local.R` fits every lme with `control = lmeControl(opt = "optim")`. On these fairly flat REML surfaces that optimiser stops early. nlme's default (`nlminb`), nlminb with tight tolerances, statsmodels, and this verification's own REML all agree with each other, and their optimum has a higher REML log-likelihood:

| model (all usable, lme) | `opt = "optim"` (= CSV under check) | nlminb default | this verification | REML logLik optim / nlminb |
|---|---|---|---|---|
| chi_PUT ~ age + sex + cell | 0.119622, SE 0.316704, p .708713 | 0.119401, SE 0.316661, p .709186 | 0.119401, p .709186 | −351.0250869 / **−351.0250699** |
| chi_PUT ~ age + sex | 0.144974, p .641830 | 0.144962, p .641855 | 0.144962, p .641855 | −356.6427959 / **−356.6427958** |

- **Size.**
  - Largest slope difference: 2.2e-4 ppb/y, in chi_PUT + sex + cell.
  - Largest relative slope difference: 0.62 %, for ref_ppb age-only. That slope is −0.000243 vs −0.000245, a 2e-6 absolute difference on an essentially zero slope.
  - SEs differ by ≤ 0.014 % and p values by ≤ 0.0005 absolute.
  - The df match in every row.
  - All first-usable OLS rows of D are exact (≤ 6e-12).
- **Consequence for the README.** Only five printed digits change, in L245–246 (table below). No "CI excludes registered range" flag changes, and no reading changes.
- **Fix.** Drop `opt = "optim"` (use the nlme default), or tighten `msTol`/`tolerance`, then regenerate D.

### 3.2 README numbers that fail the rounding check (10 of 1,359)

| README line | printed | correct (recomputed) | cause |
|---|---|---|---|
| L245 chi_PUT, all usable lme, + sex: p | .6418 | .6419 (0.641855) | lme optimiser (3.1) |
| L246 chi_PUT, all usable lme, + sex + cell | 0.120 [-0.531, 0.771], p .7087 | **0.119 [-0.532, 0.770], p .7092** | lme optimiser (3.1) |
| L249 r2s_GP, first-usable, age only: p | 3.4e-06 | **3.3e-06** (3.348e-06) | double rounding; the CSV under check holds the same 3.348e-06 |
| L253 r2s_GP, all usable lme, + sex: p | 6.6e-06 | **6.5e-06** (6.548e-06) | double rounding |
| L302 and L306 putamen R2* standardised rate | 0.25 SD/y (twice) | **0.24** (0.2445; the D.1 table correctly prints 0.245) | double rounding 0.2445 → 0.245 → 0.25. With 0.24 vs 0.24, the "matches GP" statement only gets closer. |
| L307 putamen chi slope range, all usable | "0.08-0.15 ppb/y" | **0.08-0.14** (max 0.14496, or 0.14497 in their CSV) | double rounding |

Values within tolerance, noted for completeness:

- **L56–57 floating-point noise.** "Identity error 3.4e-14" and "beta range 1.1e-14 / SE range 3.6e-15" here are 1.5e-14 and 5.3e-15 / 1.7e-14. These are machine noise and not reproducible to the digit; better to write "< 1e-13".
- **L216 "two-thirds".** The label-free refit removes 70.8 % of the point estimate, which is slightly more than two-thirds.

### 3.3 Builder summary (workflow text)

- "Adding back the whole-brain reference … (β 4.485 vs 4.486)": the referenced chi_para β is 4.48549 and the unreferenced one 4.48530, so both round to **4.485**. README L58 is correct.
- "Heterogeneity across tertiles: p ." is truncated. The value is **F(2, 64) = 2.44, p = .095** (README L174 is correct).
- Every other number in the summary reproduces at its printed precision:
  - 4.49 [0.23, 8.74], p .039; 0.48 [−2.39, 3.35], p .74; SN-only 5.71 (p .016) / 0.17;
  - 68/80; smallest q .096;
  - the joint model 4.16 [−1.90, 10.22], SE 3.03, p .175, t_B 3.39, p_B .0006;
  - the single-covariate table; Dice +4.29, p .0019;
  - R2* 1.098 (p .035) / 1.011 (p .068); Dice and R2* p .36;
  - the tertile counts and Fisher p .32; −0.95 / 14.2 / 7.4; 4.05 (p .19); 10.36 [0.40, 20.32], p .042, n 53.

---

## 4. Wording that over-claims

The README's status line is correct: everything is post hoc or registered-exploratory, and nothing is called a finding. The following sentences still read more strongly than the analyses allow.

| # | where | text (shortened) | problem | suggested wording |
|---|---|---|---|---|
| W1 | L112 (A.4) | "an iron-like (paramagnetic) increase fits better than a loss of diamagnetic (myelin-like) material" | No model comparison was done. The SN-VTA chi_para q is .118. X2 is a deterministic restatement of β_chi and β_R2*. The registered K2 rule allows iron wording only when K2 is met, and the frozen K2 sentence says the iron-versus-myelin source is unresolved. Iron language about a positive-direction effect also sits uneasily with "the positive direction remains unclaimed". | "Under the model's assumptions the point estimates fall mainly in chi_para. No component test survives BH-FDR. The iron-versus-myelin source remains unresolved (registered K2 wording)." |
| W2 | L67–68 (bold) and the builder summary bullets | SN-VTA chi_para printed in bold with "p .039"; SN-only "p .016" with no q | Headline emphasis on nominal p values in an FDR family where nothing survives. | Remove the bold. Print q (.118 / .096) next to every X2 p. |
| W3 | L112, L88–103, builder summary | "The pattern holds over Dr 137–548" | At Dr 548 the SN-VTA chi_para CI includes 0 (p .059), and the README itself calls Dr 548 physically implausible (68/80 chi_dia > 0). | "The ordering of the point estimates (chi_para > chi_dia) is the same at all four Dr values. The chi_para CI excludes 0 only for Dr ≤ 321." |
| W4 | L112, builder summary | "Neither component differs in GP" / "GP: no effect in either component" | Non-significance is not absence. The GP CIs are ±4.7 and ±2.7 ppb. | "No detectable GP difference in either component." |
| W5 | L178 (B.4) | "not significant (Fisher p .32), but it explains why Dice adjustment attenuates the estimate" | A causal claim based on a non-significant imbalance. | "…and is consistent with the attenuation under Dice adjustment." |
| W6 | L178 (B.4) | "The registration-quality explanation predicts … It does not: within the two upper tertiles the estimate is larger (+10.4, p .042)" | A post hoc subsample (n 53) with a nominal positive-direction p. Tertile heterogeneity p is .095, and all three within-tertile CIs include 0. This cannot refute an explanation. | "The point estimate does not shrink in the upper two tertiles (+10.4 ppb, n 53, post hoc). Heterogeneity across tertiles is not significant (p .095), so these subsamples cannot test the registration-quality explanation." |
| W7 | L216 (C) | "They do share poor registration." / "This fits partial-volume dilution of a bright nucleus." / "two-thirds" | Session 10a248 has the median TD Dice (50th percentile), so only 3 of 4 have poor registration. The mechanism was not examined, since no image was opened. The removed share is 71 %. | "Three of the four share poor registration… This would be consistent with partial-volume dilution (not examined)… about 70 %." |
| W8 | L304 (D.5 heading), L313 | "the putamen failure looks QSM-specific…"; "These diagnostics narrow it to a putamen-chi-specific problem, not a sample, age-range or global-QSM failure" | Contradicts L310–312: explanation (b), a real concurrent diamagnetic increase, cannot be excluded, and it would be neither a QSM problem nor a failure. | "The putamen chi slope is low relative to putamen R2*, while GP chi tracks GP R2*. This is consistent with a putamen-specific QSM measurement issue or with a real diamagnetic change; these data cannot separate them." |
| W9 | L306 | "carry a clear developmental iron signal in the putamen" | R2* is iron-sensitive, not iron-specific (myelin and water also contribute). | "a clear developmental R2* signal (iron-sensitive)". |
| W10 | L309 | "Referencing is not the cause." | An inference from the absence of age slopes in the reference values (p .55 and .99). | "We found no sign that referencing drives it: …" |
| W11 | L277–278, L307 | "The formal test of the shortfall … p .0492"; "shortfall test p .049" | The first-usable test is borderline and depends on the model and on Dr (`mine/D_putamen_shortfall_by_Dr.csv`). It gives p .051 with + sex, .055 with + sex + cell, .11 at Dr 321 and .61 at Dr 548. The lme version holds for Dr ≤ 321 (p ≤ .02) but not at 548 (p .14–.19). Also, the text's lme figure (−0.501) is the **age-only** lme, while the D.2 table uses the + sex + cell lme (chi_dia −0.515, p .0061). | Call it a "model-dependent shortfall estimate", report its dependence on Dr and covariates, and name the covariate set. |
| W12 | L313 | "The SN-VTA itself shows age-related R2* (0.31-0.50 s^-1/y) and chi (about 1.2-1.3 ppb/y) increases." | The first-usable SN-VTA chi slope is 1.31 [−0.07, 2.69], p .063, and the within-person slope 0.33 [−1.69, 2.34]. The quoted ranges omit the + sex + cell first-usable rows (R2* 0.25, p .053; chi 1.02, p .24). | Give estimates with CIs, for example: "SN-VTA R2* rises with age (0.31 s^-1/y first-usable, p .019; 0.50 lme). The SN-VTA chi slope is positive but imprecise (1.31 [−0.07, 2.69] first-usable; 1.24 [0.21, 2.27] lme; within-person 0.33 [−1.69, 2.34])." |
| W13 | L3, builder summary | "A (X2), which is a registered Tier-4 exploratory module" | The constants were fixed after unlock. This is declared at L55, but the status line reads as if X2 ran as registered. | "X2: registered Tier-4 exploratory module, run with a declared deviation (constants fixed after unlock from the paper, not tuned)." |

---

## 5. Other checks and observations

**Dimov et al. 2022 (Tomography 8(3):1544–1551, doi:10.3390/tomography8030127).** Re-read via PMC9228115 and the Europe PMC JATS XML; mdpi.com also returned 403 to this verifier.

| README statement | check |
|---|---|
| Sec. 2.1 heading "MR Signal Modeling for R2*QSM"; Eq. 1: R2* + i2πf = r+\|χ+\| + r−\|χ−\| + iγ·d∗(χ+ + χ−); constraints χ+ > 0, χ− < 0 | confirmed |
| r+ and r− "are relaxometric constants equal to 274 Hz/ppm", cited to ref. 22 (Dimov et al., J Neuroimaging 2022, doi:10.1111/jon.13014) | confirmed: one value for both components, ref. 22 |
| 321 Hz/ppm is the theoretical value (ref. 23, Yablonskiy & Haacke 1994) | confirmed. It appears in the Discussion. |
| Eq. 1 has no separate R2 term; the Discussion allows a possibly unaccounted constant, argued to be small because ventricular CSF R2* is low; no numeric R2 is given | confirmed. Setting R2 = 0 is the verifier-endorsed reading. Because R2 cancels exactly in every group contrast, the choice cannot affect any β. |
| Paper solves a regularised, constrained voxelwise problem (Eq. 2; conjugate gradient / Gauss–Newton) | confirmed (λ1 = 1500 is also reported). The README's "closed form" is the unregularised, unconstrained point solution, as stated. |
| Acquisition: 3 T, 0.5 mm isotropic, 8 echoes, first TE = ΔTE = 5.9 ms | confirmed |
| "p. 1545" | **not verifiable**: the HTML and XML carry no per-page markers. The article spans 1544–1551. |
| Algebra: chi_para = (chi + 1000·R2'/Dr)/2, chi_dia = (chi − 1000·R2'/Dr)/2 (ppb) | correct, and units consistent. The identity β_para = (β_chi + 1000·β_R2*/Dr)/2 holds to 1.5e-14. |

**Provenance and other claims checked.**
- **Deviation from Q25.** `LEADchi_amendment_2_FROZEN.md` contains none of "X2", "Dimov", "separation", "relaxometric", "274" or "Hz/ppm". The declared deviation is real.
- **QSMxT help text.** `prereg/server_capture/qsmxt_run_help.txt` line 795 lists "R2*-QSM closed-form (Dimov 2022)".
- **Hashes.**
  - `SHA256SUMS` in validation_local verifies with no failures.
  - The input SHA-256 `f669bc6f…f215` equals the confirmatory input hash.
  - The five frozen helper hashes equal `mirror/scripts/analysis/SHA256SUMS`.
- **results.json.** Every table in results.json equals the matching CSV. The `verify/statsmodels_check.json` all_pass is TRUE, with a maximum difference of 8e-14. That check covers only the HC3 algebra on the exported designs, as the README states, which is why this from-raw recomputation was needed.
- **Write boundary.** `validation_local.R` writes only to `OUT` and `OUT/verify`, and `verify_statsmodels.py` only to `OUT/verify`.
  - Other LEADchi files changed on 2026-10-07: `statistician_package/`, `human_qc_kit/`, `time_evidence/` and `manuscript/`. They belong to other concurrent tasks and are not written by this code.
  - The "no image was opened" claim cannot be verified from files. The code reads only the CSV.
- **Larsen et al. 2020** (J Neurosci 40(9):1810–1818, PMC7046454). Confirmed:
  - 922 participants, 1,236 sessions;
  - two echoes, TE 2.69 / 5.27 ms; 3.75 mm in-plane;
  - EMM R2* 16.96 (putamen) and 22.33 (GP);
  - significant-increase windows: 10.2–26.9 y in the putamen; 8.2–25.9 y (males) and 8.2–22.2 y (females) in the GP;
  - F 85.8 / 76.7; no numeric per-year slope in the text.

  The field strength (3 T) was not stated in the fetched text. This study's means (18.60 and 30.53 s^-1) are confirmed.

**Observations outside the numbers (no change needed for this verification, but worth recording).**
1. **BH family.** Q25 says "BH-FDR within block". The README uses an X2-only family of 6 (and 18). If "block" means the whole Tier-4 block (X1–X7), the registered family is larger and q can only rise, so "nothing survives" holds either way. Still, the README should say that the family definition is a choice.
2. **R2* voxel sets.** README decl. 6 matches `extract_rois.py`: R2* is taken over Mask 3 ∩ finite, chi over support S minus implausible voxels. However, prereg line 253 registered R2*_ROI on the "same voxel sets as the matching χ ROI, restricted to Mask 3". The implemented definition departs from that text. The departure also affects the frozen K2/PC3. It is not labelled a deviation in the README. Check whether it is already listed in the deviation annex.
3. **roi_vol covariate.** The SN-VTA atlas volume is used as `roi_vol` in the SN-only and GP models. This matches the frozen SN-only sensitivity row (β 5.87448, reproduced exactly) and is declared, but for GP it is not the ROI's own volume.
4. **Within-tertile bound tests.** The B bound is not excluded inside the lowest tertile (p_B .22 / .30) or in the reduced high-tertile model (.096). The README's "every pooled specification … bound excluded" is accurate only because it is restricted to pooled models. Keep that qualifier.

---

## 6. Reproduce

```
cd /path/to/Penn_LEAD/LEADchi/results_posthoc/validation_local/verification
python3 verify_independent.py      # recomputation + cell-by-cell comparison (about 3 s)
python3 check_readme.py            # every printed README number vs the recomputation
Rscript diag_nlme_optimizer.R      # diagnostic only: nlme optim vs nlminb (prints, writes nothing)
```

Environment: python3 with numpy 2.4.6, scipy 1.16.0, pandas 2.3.1 and statsmodels 0.14.5 (statsmodels is used only for the MixedLM cross-check). R 4.4.1 with nlme 3.1-164 is used only for the diagnostic.
