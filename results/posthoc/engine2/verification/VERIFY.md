# Verification of `results_posthoc/engine2` (POST HOC engine-2 comparison)

Verified 2026-10-07. The verification was read-only: the only files written are in this `verification/` folder. No server access was used.

## How it was checked
- **Own code.** `verify_engine2.py` and `compare.py` were written from scratch in Python 3 (numpy 2.3.1, pandas 2.3.1, scipy 1.16.0, statsmodels 0.14.5) from the task specification. `code/engine2_analysis.R` was neither run nor copied.
- **Inputs.**
  - `raw/roi_qsmxt_e2_lap_vsharp_{rts,tv}.csv` and `raw/roi_qsmxt_e2rec_m14_{rts,tv}.csv`. The SHA256SUMS check passed for all 6 files.
  - `mirror/unlocked/key_plain.csv` and `mirror/unlocked/merged_qsmxt_primary_SNVTA.csv`, using `first_usable == 1`, which gives n = 80.
  - For the primary k only: `results_posthoc/recovery/raw/roi_qsmxt_{primary,rec_m14}.csv`.
- **M1.** The model is fitted by OLS on 12 columns:
  - intercept, PS_e, ADHD_e and PS_e:ADHD_e;
  - a natural cubic spline basis for age with boundary knots at the min/max and interior knots at the type-7 1/3 and 2/3 quantiles (an ESL truncated-power basis, which spans the same column space as `ns(df=3)` plus the intercept);
  - sex_e, race3 with sum-to-zero coding, z(roi_vol) and z(Q_phase).

  SEs are HC3 with a t reference on df = n − 12.
- **Huber scale.** Huber ψ (k = 1.345) with the proposal-2 scale, using the IRLS scheme of MASS::rlm: LS start, MAD start scale, a scale update before each reweighting, and a stop when the relative change in the residuals is ≤ 1e-4.
  - This reproduces the reported σ_B and Huber β to about 1e-9.
  - Iterating to the exact fixed point changes the scale by about 1e-4 relative (for example, 10.5518 against 10.5505 for the primary). No rounded number changes.
  - A statsmodels OLS-HC3 and RLM cross-check agrees.
- **Covariates.** roi_vol and Q_phase come from the merged (primary) table. They are identical to the engine-2 table's columns (maximum absolute difference 0), so the choice does not matter.
- **Tolerances.** For results.json, a value passes if its relative difference is ≤ 1e-6 or it is within the file's 8-decimal rounding. For the README, a value passes if it matches the printed rounding.
- **Output files.** `recomputed.json` holds all recomputed values, the sensitivity refits and the per-ROI contrast ratios. `number_checks.md` is the full table reproduced below.

## Summary
- **Numbers.** All 179 numeric and label checks are **OK**: 97 against results.json and 82 against the README. The agreement CSVs also match my merge exactly (80 rows, same values, same NA pattern, same status).
- **Non-OK items.** All of them are in the description and wording, not in the arithmetic: D1–D5 below (factual or description mismatches) and C1–C12 (over-statements and claim issues).

## D. Factual / description mismatches

**D1. Wrong reason given for the dropped session.** MISMATCH
- README ("Sample"): "One session's engine-2 SN–VTA value failed the plausibility rules and was dropped (n = 79)."
- Recomputed: the dropped first-usable session (`sub-20866_ses-1`) has engine-2 `status = mask_mismatch`. Its chi_SNVTA is NaN because the run failed, not because of a plausibility rule. Its `fail_SNVTA` is NA.
- Fix: "One first-usable session had no engine-2 value (run status mask_mismatch) and was dropped (n = 79)."

**D2. Engine-2 plausibility failures were kept in the sample.** MISMATCH and missing disclosure
- **What happened.** 16 of the 79 analysed sessions have engine-2 `fail_SNVTA = 1`.
  - The rule they fail is the registered SN–VTA in-support share < 0.95, with values 0.82–0.95. f_imp is 0 for all of them.
  - Under the registered criterion (c') these sessions would not be usable for engine 2.
  - In the primary, all 80 first-usable sessions pass (minimum in-support share 0.969).
- **Why the README is inaccurate.** It says the plausibility rules are "shared with the primary" and implies that they were applied. They were computed but not applied to the engine-2 SN–VTA values: rows were dropped only when the value was non-finite.
- **Sensitivity check** (dropping these 16 sessions, n = 63):

  | Inversion | β_PS [95% CI] | p | δ | p_B | Bound-test outcome |
  |---|---|---|---|---|---|
  | RTS | +2.53 [−2.94, 8.01] ppb | .36 | +0.31 | .008 | B |
  | TV | +3.34 [−2.37, 9.05] ppb | .25 | +0.38 | .003 | B |

  The bound test still rejects, but the point estimates and δ shrink by about 25–40% (β: RTS −34%, TV −24%; δ: RTS −37%, TV −29%).
- Fix: disclose that the engine-2 plausibility flags were not used for sample selection (the sample was fixed to the primary's first-usable list), state the count (16/79), and add this sensitivity row.

**D3. "results.json: all the numbers above."** MISMATCH
- results.json lacks the following README numbers:
  - for the primary: the Huber β, the bound, p_B and the simple effects;
  - all per-stratum bound p_B values;
  - the primary k;
  - the expected adult deficits (−6.7/−5.2/−3.5);
  - the adult-deficit one-sided p values.
- I recomputed all of them and they match, but they are not traceable to the file.
- Fix: add them to results.json, or change the sentence to "most of the numbers above".

**D4. The δ_B parenthetical is mis-stated.** Wording mismatch
- README: "δ_B = 0.581 (the split-echo r̂_L of the primary engine …)".
- δ_B is not r_L. δ_B = 0.6·√r_L, where r_L = 0.938 is the primary engine's split-echo lower bound.
- Fix: "δ_B = 0.6·√r_L = 0.581, with r_L = 0.938 from the primary engine's split-echo reliability".

**D5. Sign of the bound header.** Label mismatch
- The column header reads "Bound (δ_B·σ_B, ppb)", but the values printed are −δ_B·σ_B (−6.13, −3.92, −4.48).
- Fix: "Bound (−δ_B·σ_B, ppb)".

## C. Claim audit (post hoc: nothing may be called a finding)

**C1. "Main result … H1 again gives B … Neither the registered outcome nor the direction depends on the QSM engine."** Over-statement
- **What the analysis covers.**
  - It is post hoc and uses one alternative engine. RTS and TV share the same unwrapping, echo fit and V-SHARP.
  - Outcome B is borrowed: δ_B comes from the primary engine's r_L, and engine 2's own precision gates (PC4, r_L ≥ .25) were not assessed.
  - All three 95% CIs include 0, so no direction is established in any engine.
  - 16/79 sessions fail engine-2 plausibility (D2).
- **Suggested fix.**
  - Rename "Main result" to "Engine-2 check (post hoc)".
  - Wording: "With one alternative engine (two inversions), the bound test also rejected the registered bound (p_B ≤ .0003), and the point estimates were positive (δ +0.50 and +0.53; all 95% CIs include 0). This is consistent with, but does not show, independence from the QSM engine. The engine-2 B uses the primary's δ_B, and engine-2 precision gates were not assessed."
  - Label the "Outcome" row for the engine-2 columns as "bound test p < .05 (B pattern, not a registered outcome)".

**C2. "Engine 2 spreads the positive shift evenly over both strata (+3.6 and +4.1 ppb)."** Over-statement
- Both simple effects are non-significant with wide CIs: RTS [−2.3, 9.4] and [−2.2, 10.3].
- The PS×ADHD interaction CI is [−7.9, 9.0] for RTS and [−9.4, 8.8] for TV, so homogeneity is not shown.
- The quoted values are the RTS ones; for TV they are +4.5 and +4.2.
- Fix: "Under engine 2 the two stratum point estimates are similar (RTS +3.6 and +4.1; TV +4.5 and +4.2 ppb; all CIs include 0)."

**C3. "In the primary, PS-only − TD (+8.2) was concentrated in one stratum."** Over-statement, and a clash with registered terms
- The primary interaction is −6.39 ppb, 95% CI [−18.8, 6.0], p = .31. The difference between strata is not established.
- "Concentrated" is the wording of the registered K1c label ("Concentrated in comorbid youth"), which refers to the other stratum. Reusing it here invites confusion.
- Fix: "In the primary the PS-only − TD point estimate (+8.2, p = .044) was larger than PS+ADHD − ADHD-only (+1.8); the interaction was not significant (p = .31)."

**C4. "The cell-level pattern … is therefore engine-dependent and should not be interpreted."** Unsupported inference
- No test compares the simple effects across engines, so "engine-dependent" is itself an interpretation.
- Fix: "The primary's larger PS-only − TD estimate was not reproduced by engine 2 and should not be interpreted."

**C5. "The ADHD-averaged effect, which is the registered estimand, is stable across engines."** Over-statement
- Fix: "is similar across engines in standardised units (δ +0.44, +0.50, +0.53; all 95% CIs include 0)."

**C6. "Absolute ppb values are not comparable … standardised effects (δ) are the better common scale."** Over-statement
- δ is not invariant across engines either. Signal recovery k and residual SD scale differently:

  | Inversion | k ratio vs primary | σ_resid ratio vs primary | A true group difference would appear in δ at |
  |---|---|---|---|
  | RTS | 0.77 | 0.68 | about 1.13 × the primary δ |
  | TV | 0.52 | 0.73 | about 0.72 × the primary δ |

- Fix: "δ is more comparable than ppb but still depends on each engine's recovery (k) and noise; the k-scaled comparison below addresses this."

**C7. "Each engine's observed estimate excludes it at one-sided p = .0002, 5e-5 and 5e-4."** Over-statement
- The values themselves are OK.
- "Excludes" is too strong for a post hoc test that treats k (the first-usable mean) and −14.24 ppb as known constants.
- Fix: "is incompatible with it (post hoc one-sided p = .0002, 5e-5 and 5e-4; k and the adult value treated as fixed)."
- Related context worth one sentence: under TV the registered bound (−4.48 ppb) is less stringent than the adult deficit as TV would record it (−3.52 ppb). TV's "B" alone therefore does not exclude an adult-magnitude deficit; only the k-scaled test does. For the primary (−6.13 vs −6.73) and RTS (−3.92 vs −5.16), the bound is the stricter of the two.

**C8. "Bound rejected in the comorbid stratum? (post hoc): no (.052) / yes (.006) / yes (.006)."** Risk of over-interpretation
- The values are OK.
- The yes/no split around .05 invites a reading of difference in significance between engines.
- The strata were not designed for this test, there is no multiplicity control, and only the comorbid stratum is quoted. The PS-only stratum p_B values are .0003, .006 and .003.
- Fix: report the p_B values for both strata without yes/no, and mark them as descriptive.

**C9. "Engine 2 gives about half the absolute deep-grey-matter contrast."** Over-generalised
- For SN–VTA the claim holds: engine-2 mean / primary mean is 0.51 for RTS and 0.45 for TV.
- For other deep grey matter ROIs on the same 79 sessions it does not:

  | ROI | RTS ratio | TV ratio |
  |---|---|---|
  | GP | 0.75 | 0.79 |
  | Caudate | 0.97 | 0.90 |
  | DGM composite | 0.84 | 0.84 |
  | Striatum | 0.79 | 0.69 |
  | Putamen | 0.57 | 0.44 |
  | Red nucleus | 0.59 | 0.44 |

- Fix: "Engine 2 gives about half the primary's SN–VTA contrast (51% RTS, 45% TV); across other deep grey matter ROIs the ratio is 44–97%."

**C10. "It ranks youths similarly, but not identically."** Mild over-statement
- Spearman is 0.71 and ICC is 0.68–0.69.
- Fix: "It ranks youths moderately consistently (Spearman 0.71)."

**C11. "This is common across QSM pipelines."** Unsupported
- No citation is given.
- Fix: add a citation, for example a multi-pipeline or QSM-challenge comparison, or delete the sentence.

**C12. "No QSMxT/QSM.rs component is used … all independent implementations" / "an engine that shares no reconstruction code".** Accurate as stated, but incomplete
- The independence is of code, not of method or inputs:
  - both chains use V-SHARP-type background removal;
  - both share Mask 3, the ROI warps, the two-step reference and the coded inputs.
- Fix: add "(independent code; both chains use V-SHARP-type background removal and share masks, ROI warps and reference)".

Not verifiable from the given inputs (no issue raised): the QSM.jl, Julia and FFTW versions, the chain settings, the NLTV divergence, and the script hashes on server.

## Number checks (all OK)

| source | item | reported | recomputed | check |
|---|---|---|---|---|
| results.json | primary.beta_ppb | 4.9628015 | 4.9628015 | OK |
| results.json | primary.se_hc3 | 3.1250804 | 3.1250804 | OK |
| results.json | primary.p_two | 0.11691305 | 0.11691305 | OK |
| results.json | primary.sigma_resid | 11.255903 | 11.255903 | OK |
| results.json | primary.sigma_B | 10.550454 | 10.550454 | OK |
| results.json | primary.delta | 0.44090655 | 0.44090655 | OK |
| results.json | primary.outcome | B | B | OK |
| results.json | primary.ci95_ppb[0] | -1.2731993 | -1.2731993 | OK |
| results.json | primary.ci95_ppb[1] | 11.198802 | 11.198802 | OK |
| results.json | e2_rts.agreement.n_first_usable | 80 | 80 | OK |
| results.json | e2_rts.agreement.n_valid_e2 | 79 | 79 | OK |
| results.json | e2_rts.agreement.icc31 | 0.67692751 | 0.67692751 | OK |
| results.json | e2_rts.agreement.pearson | 0.72903413 | 0.72903413 | OK |
| results.json | e2_rts.agreement.spearman | 0.70783836 | 0.70783836 | OK |
| results.json | e2_rts.agreement.mean_primary | 52.993468 | 52.993468 | OK |
| results.json | e2_rts.agreement.mean_e2 | 27.084023 | 27.084023 | OK |
| results.json | e2_rts.agreement.sd_primary | 11.955348 | 11.955348 | OK |
| results.json | e2_rts.agreement.sd_e2 | 8.0953339 | 8.0953339 | OK |
| results.json | e2_rts.agreement.bland_altman_mean_diff | -25.909445 | -25.909445 | OK |
| results.json | e2_rts.agreement.bland_altman_sd_diff | 8.2066577 | 8.2066577 | OK |
| results.json | e2_rts.agreement.gp_pearson | 0.89384564 | 0.89384564 | OK |
| results.json | e2_rts.M1.beta_ppb | 3.8148968 | 3.8148968 | OK |
| results.json | e2_rts.M1.se_hc3 | 2.1656548 | 2.1656548 | OK |
| results.json | e2_rts.M1.p_two | 0.08270844 | 0.082708443 | OK |
| results.json | e2_rts.M1.sigma_resid | 7.6506521 | 7.6506521 | OK |
| results.json | e2_rts.M1.sigma_B | 6.7476417 | 6.7476417 | OK |
| results.json | e2_rts.M1.delta | 0.49863681 | 0.49863681 | OK |
| results.json | e2_rts.M1.outcome | B | B | OK |
| results.json | e2_rts.M1.disconfirming | False | False | OK |
| results.json | e2_rts.M1.ci95_ppb[0] | -0.50776833 | -0.50776833 | OK |
| results.json | e2_rts.M1.ci95_ppb[1] | 8.1375619 | 8.1375618 | OK |
| results.json | e2_rts.bound_ppb_registered_deltaB | -3.921653 | -3.921653 | OK |
| results.json | e2_rts.pB | 0.00033045 | 0.00033044948 | OK |
| results.json | e2_rts.huber_beta | 2.9682804 | 2.9682804 | OK |
| results.json | e2_rts.K1.PSonly_minus_TD.est | 3.5621054 | 3.5621054 | OK |
| results.json | e2_rts.K1.PSonly_minus_TD.se | 2.9271637 | 2.9271637 | OK |
| results.json | e2_rts.K1.PSonly_minus_TD.p_two | 0.22790633 | 0.22790633 | OK |
| results.json | e2_rts.K1.PSonly_minus_TD.ci95[0] | -2.2805377 | -2.2805377 | OK |
| results.json | e2_rts.K1.PSonly_minus_TD.ci95[1] | 9.4047485 | 9.4047485 | OK |
| results.json | e2_rts.K1.PSADHD_minus_ADHD.est | 4.0676881 | 4.0676881 | OK |
| results.json | e2_rts.K1.PSADHD_minus_ADHD.se | 3.1280485 | 3.1280485 | OK |
| results.json | e2_rts.K1.PSADHD_minus_ADHD.p_two | 0.19792217 | 0.19792217 | OK |
| results.json | e2_rts.K1.PSADHD_minus_ADHD.ci95[0] | -2.1759228 | -2.1759228 | OK |
| results.json | e2_rts.K1.PSADHD_minus_ADHD.ci95[1] | 10.311299 | 10.311299 | OK |
| results.json | e2_rts.n_dropped_invalid | 1 | 1 | OK |
| results.json | e2_rts.recovery_k_L14.all.n | 177 | 177 | OK |
| results.json | e2_rts.recovery_k_L14.all.mean | 0.35378238 | 0.35378238 | OK |
| results.json | e2_rts.recovery_k_L14.all.sd | 0.05117079 | 0.051170791 | OK |
| results.json | e2_rts.recovery_k_L14.all.median | 0.35687256 | 0.35687256 | OK |
| results.json | e2_rts.recovery_k_L14.first_usable.n | 79 | 79 | OK |
| results.json | e2_rts.recovery_k_L14.first_usable.mean | 0.36208462 | 0.36208462 | OK |
| results.json | e2_rts.recovery_k_L14.first_usable.sd | 0.04462184 | 0.044621838 | OK |
| results.json | e2_rts.recovery_k_L14.first_usable.median | 0.3649145 | 0.3649145 | OK |
| results.json | e2_tv.agreement.n_first_usable | 80 | 80 | OK |
| results.json | e2_tv.agreement.n_valid_e2 | 79 | 79 | OK |
| results.json | e2_tv.agreement.icc31 | 0.6925026 | 0.6925026 | OK |
| results.json | e2_tv.agreement.pearson | 0.72686912 | 0.72686912 | OK |
| results.json | e2_tv.agreement.spearman | 0.71377799 | 0.71377799 | OK |
| results.json | e2_tv.agreement.mean_primary | 52.993468 | 52.993468 | OK |
| results.json | e2_tv.agreement.mean_e2 | 24.01098 | 24.01098 | OK |
| results.json | e2_tv.agreement.sd_primary | 11.955348 | 11.955348 | OK |
| results.json | e2_tv.agreement.sd_e2 | 8.7357393 | 8.7357393 | OK |
| results.json | e2_tv.agreement.bland_altman_mean_diff | -28.982488 | -28.982488 | OK |
| results.json | e2_tv.agreement.bland_altman_sd_diff | 8.2107734 | 8.2107734 | OK |
| results.json | e2_tv.agreement.gp_pearson | 0.89936967 | 0.89936967 | OK |
| results.json | e2_tv.M1.beta_ppb | 4.3712893 | 4.3712893 | OK |
| results.json | e2_tv.M1.se_hc3 | 2.3137094 | 2.3137094 | OK |
| results.json | e2_tv.M1.p_two | 0.06318008 | 0.063180083 | OK |
| results.json | e2_tv.M1.sigma_resid | 8.2341308 | 8.2341308 | OK |
| results.json | e2_tv.M1.sigma_B | 7.7012034 | 7.7012034 | OK |
| results.json | e2_tv.M1.delta | 0.5308744 | 0.5308744 | OK |
| results.json | e2_tv.M1.outcome | B | B | OK |
| results.json | e2_tv.M1.disconfirming | False | False | OK |
| results.json | e2_tv.M1.ci95_ppb[0] | -0.24689406 | -0.24689406 | OK |
| results.json | e2_tv.M1.ci95_ppb[1] | 8.9894726 | 8.9894726 | OK |
| results.json | e2_tv.bound_ppb_registered_deltaB | -4.4758522 | -4.4758522 | OK |
| results.json | e2_tv.pB | 0.00014538 | 0.00014537993 | OK |
| results.json | e2_tv.huber_beta | 3.3354725 | 3.3354725 | OK |
| results.json | e2_tv.K1.PSonly_minus_TD.est | 4.5209011 | 4.5209011 | OK |
| results.json | e2_tv.K1.PSonly_minus_TD.se | 3.138419 | 3.138419 | OK |
| results.json | e2_tv.K1.PSonly_minus_TD.p_two | 0.15438214 | 0.15438214 | OK |
| results.json | e2_tv.K1.PSonly_minus_TD.ci95[0] | -1.7434095 | -1.7434095 | OK |
| results.json | e2_tv.K1.PSonly_minus_TD.ci95[1] | 10.785212 | 10.785212 | OK |
| results.json | e2_tv.K1.PSADHD_minus_ADHD.est | 4.2216775 | 4.2216775 | OK |
| results.json | e2_tv.K1.PSADHD_minus_ADHD.se | 3.3607649 | 3.3607649 | OK |
| results.json | e2_tv.K1.PSADHD_minus_ADHD.p_two | 0.21341846 | 0.21341846 | OK |
| results.json | e2_tv.K1.PSADHD_minus_ADHD.ci95[0] | -2.4864373 | -2.4864373 | OK |
| results.json | e2_tv.K1.PSADHD_minus_ADHD.ci95[1] | 10.929792 | 10.929792 | OK |
| results.json | e2_tv.n_dropped_invalid | 1 | 1 | OK |
| results.json | e2_tv.recovery_k_L14.all.n | 177 | 177 | OK |
| results.json | e2_tv.recovery_k_L14.all.mean | 0.24187828 | 0.24187828 | OK |
| results.json | e2_tv.recovery_k_L14.all.sd | 0.06419018 | 0.064190182 | OK |
| results.json | e2_tv.recovery_k_L14.all.median | 0.24990794 | 0.24990794 | OK |
| results.json | e2_tv.recovery_k_L14.first_usable.n | 79 | 79 | OK |
| results.json | e2_tv.recovery_k_L14.first_usable.mean | 0.24743985 | 0.24743985 | OK |
| results.json | e2_tv.recovery_k_L14.first_usable.sd | 0.0562766 | 0.056276604 | OK |
| results.json | e2_tv.recovery_k_L14.first_usable.median | 0.25426658 | 0.25426658 | OK |
| README | Agreement n | 79 | 79 | OK |
| README | ICC(3,1) RTS | 0.68 | 0.67692751 | OK |
| README | Pearson RTS | 0.73 | 0.72903413 | OK |
| README | Spearman RTS | 0.71 | 0.70783836 | OK |
| README | Mean chi RTS | 27.1 | 27.084023 | OK |
| README | SD RTS | 8.1 | 8.0953339 | OK |
| README | GP Pearson RTS | 0.89 | 0.89384564 | OK |
| README | ICC(3,1) TV | 0.69 | 0.6925026 | OK |
| README | Pearson TV | 0.73 | 0.72686912 | OK |
| README | Spearman TV | 0.71 | 0.71377799 | OK |
| README | Mean chi TV | 24.0 | 24.01098 | OK |
| README | SD TV | 8.7 | 8.7357393 | OK |
| README | GP Pearson TV | 0.90 | 0.89936967 | OK |
| README | Mean chi primary | 53.0 | 52.993468 | OK |
| README | SD primary | 12.0 | 11.955348 | OK |
| README | beta_PS Primary | 4.96 | 4.9628015 | OK |
| README | CI lower Primary | -1.27 | -1.2731993 | OK |
| README | CI upper Primary | 11.20 | 11.198802 | OK |
| README | p two-sided Primary | .117 | 0.11691305 | OK |
| README | delta Primary | 0.44 | 0.44090655 | OK |
| README | bound ppb Primary | -6.13 | -6.1318046 | OK |
| README | p_B Primary | .0004 | 0.00035191843 | OK |
| README | Outcome Primary | B | B | OK |
| README | Huber beta Primary | 3.94 | 3.9355014 | OK |
| README | PS-only - TD est Primary | 8.16 | 8.1574215 | OK |
| README | PS-only - TD p Primary | .044 | 0.044281483 | OK |
| README | PS+ADHD - ADHD-only est Primary | 1.77 | 1.7681816 | OK |
| README | PS+ADHD - ADHD-only p Primary | .71 | 0.71291195 | OK |
| README | comorbid-stratum p_B Primary | .052 | 0.051693273 | OK |
| README | comorbid-stratum bound rejected Primary | no | no | OK |
| README | beta_PS RTS | 3.81 | 3.8148968 | OK |
| README | CI lower RTS | -0.51 | -0.50776833 | OK |
| README | CI upper RTS | 8.14 | 8.1375618 | OK |
| README | p two-sided RTS | .083 | 0.082708443 | OK |
| README | delta RTS | 0.50 | 0.49863681 | OK |
| README | bound ppb RTS | -3.92 | -3.921653 | OK |
| README | p_B RTS | .0003 | 0.00033044948 | OK |
| README | Outcome RTS | B | B | OK |
| README | Huber beta RTS | 2.97 | 2.9682804 | OK |
| README | PS-only - TD est RTS | 3.56 | 3.5621054 | OK |
| README | PS-only - TD p RTS | .23 | 0.22790633 | OK |
| README | PS+ADHD - ADHD-only est RTS | 4.07 | 4.0676881 | OK |
| README | PS+ADHD - ADHD-only p RTS | .20 | 0.19792217 | OK |
| README | comorbid-stratum p_B RTS | .006 | 0.0064632703 | OK |
| README | comorbid-stratum bound rejected RTS | yes | yes | OK |
| README | beta_PS TV | 4.37 | 4.3712893 | OK |
| README | CI lower TV | -0.25 | -0.24689406 | OK |
| README | CI upper TV | 8.99 | 8.9894726 | OK |
| README | p two-sided TV | .063 | 0.063180083 | OK |
| README | delta TV | 0.53 | 0.5308744 | OK |
| README | bound ppb TV | -4.48 | -4.4758522 | OK |
| README | p_B TV | .0001 | 0.00014537993 | OK |
| README | Outcome TV | B | B | OK |
| README | Huber beta TV | 3.34 | 3.3354725 | OK |
| README | PS-only - TD est TV | 4.52 | 4.5209011 | OK |
| README | PS-only - TD p TV | .15 | 0.15438214 | OK |
| README | PS+ADHD - ADHD-only est TV | 4.22 | 4.2216775 | OK |
| README | PS+ADHD - ADHD-only p TV | .21 | 0.21341846 | OK |
| README | comorbid-stratum p_B TV | .006 | 0.0059140085 | OK |
| README | comorbid-stratum bound rejected TV | yes | yes | OK |
| README | delta_B | 0.581 | 0.58118868 | OK |
| README | k primary mean | 0.47 | 0.47282578 | OK |
| README | k primary SD | 0.05 | 0.05155949 | OK |
| README | k RTS mean | 0.36 | 0.36208462 | OK |
| README | k RTS SD | 0.04 | 0.044621838 | OK |
| README | k TV mean | 0.25 | 0.24743985 | OK |
| README | k TV SD | 0.06 | 0.056276604 | OK |
| README | expected adult deficit primary | -6.7 | -6.7330392 | OK |
| README | expected adult deficit RTS | -5.2 | -5.156085 | OK |
| README | expected adult deficit TV | -3.5 | -3.5235434 | OK |
| README | adult-deficit p primary (given inputs) | .0002 | 0.00018845158 | OK |
| README | adult-deficit p primary (recomputed inputs) | .0002 | 0.00018836373 | OK |
| README | adult-deficit p RTS | 5e-5 | 4.9232711e-05 | OK |
| README | adult-deficit p TV | 5e-4 | 0.0005484311 | OK |
| README | Reading: strata RTS +3.6 | 3.6 | 3.5621054 | OK |
| README | Reading: strata RTS +4.1 | 4.1 | 4.0676881 | OK |
| README | Reading: primary PS-only - TD +8.2 | 8.2 | 8.1574215 | OK |
| README | Reading: ICC approx 0.68 (RTS) | 0.68 | 0.67692751 | OK |
| README | Reading: delta approx +0.5 (RTS) | 0.5 | 0.49863681 | OK |
| README | Reading: delta approx +0.5 (TV) | 0.5 | 0.5308744 | OK |
| README | 'about half' SN-VTA contrast RTS (ratio) | 0.5 | 0.51108229 | OK |
| README | 'about half' SN-VTA contrast TV (ratio) | 0.5 | 0.4530932 | OK |
