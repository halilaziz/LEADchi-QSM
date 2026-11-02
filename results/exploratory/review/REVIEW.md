# LEAD-chi: adversarial statistical review of the frozen confirmatory report

**Status: EXPLORATORY / REVIEW (Tier 3-4). Nothing in this document re-decides a registered outcome.** The confirmatory result is frozen: H1 final outcome **B** (beta_PS = +4.96 ppb, p = .117; p_B = .0004), K1 label none, K2 not met, K3 not triggered, K4 none, S1/S2 not supported, S3 exploratory, T1 not evaluated, n = 80. This review re-derives those numbers with independent code, audits the decision logic against the registration, and lists concerns ranked by how much they matter.

Reviewed objects: `mirror/results/confirmatory/report.md` and `results.json` (manifest SHA-256 of results.json `66b47f3b...`), the registration `prereg/LEADchi_prereg_v1.0_FINAL.md`, amendments 1/1a/1b/1c/2, `prereg/amendment1_log.md` #47-51, the frozen code `mirror/scripts/analysis/`, the unlocked table `mirror/unlocked/merged_qsmxt_primary_SNVTA.csv` (SHA-256 `f669bc6f...`, 181 data rows, 113 subjects) and the gate JSON `mirror/gates/gates_preunlock.json`.

Review code (own implementation, not derived from the frozen functions): `results_exploratory/code/review_rederive.R` (manual HC3 algebra, own first-usable rule, own contrasts, permutation, bootstrap), `review_check_statsmodels.py` (second engine), `review_split_echo.R`. Outputs in this folder (see README.md). Declared limitation: the frozen code, its verification script and this review were all written by AI systems (different instances); no human analyst has re-derived the numbers.

---

## 1. Verdict

| Question | Answer |
|---|---|
| **Numbers reproduced?** | **Yes.** Every decision quantity (M1 beta/SE/t/df/p, sigma_resid, Huber scale, sigma_B, delta, delta_B, t_B, p_B; K1 simple effects, ADHD, interaction, TOST; K2; K3 and its CSF refit; S1 L/SE/p; S2 beta/SE/p; Holm; PC2x; QC-by-group Fisher; sample sizes and cells; the automated-exclusions-only sensitivity) reproduces to relative difference <= 1e-13 in two engines (R with explicit HC3 algebra; Python statsmodels). Monte Carlo quantities reproduce within Monte Carlo error (Freedman-Lane p .0860 vs .0857; bootstrap delta CI [-0.02, 1.02] vs [-0.04, 1.02]; Bayes factors exact). |
| **Decision correct?** | **Yes, as registered and as amended before unlock.** A fails (beta > 0); B is evaluable (PC4 pass, r_L = 0.938 >= 0.25) and rejected at the registered bound; disconfirmation not met (p .117); K3 not triggered so the final outcome is the M1 outcome; K1 requires A (none); K2 one-sided p .98; K4 dominance not triggered (top-1 share .16 < .20, max \|t_i\| 4.05 < 5), 0/80 leave-one-out changes, "robust" not assessable; Holm family {S1, S2} at .025/.05 (S3 left the family pre-chi, amendment 1 #10d), both p_holm = .1996; T1 conditional on A, correctly not evaluated; PC2w INDETERMINATE -> PC2x evaluated -> putamen fails -> caveat appended to every H1 statement. The registered sentences are printed verbatim. |
| **But** | (i) outcome B carries **Tier-1 confirmatory status only because amendment 1b removed the registered G5 chain trigger** before unlock; under the registration as posted, G5 would have read "fail" and H1 would be reported **descriptively**; (ii) one co-reported number in the frozen report is wrong because of a **coding error** (the residualised split-echo bound, p = .0263; text-conforming value p = .00035), which under Annex C.4 must be logged and shown both ways; (iii) the realised power is lower than amendment 2 states; (iv) several report sentences need correction; (v) the co-reported PS-only - TD p = .044 does not survive any multiplicity handling and must not be promoted. |

**Concerns ranked** (details in section 4): 1 G5 trigger removal and confirmatory status; 2 split-echo co-report coding error and the two implausible {2,4} reconstructions (both unsigned-phase sessions); 3 what B does and does not exclude (bound insensitive to scale; B driven by the positive point estimate); 4 multiplicity and direction of the PS-only - TD co-report; 5 n shortfall and realised power; 6 PC2x partial failure (putamen) and age adjustment; 7 selection effects of the AI QC and the omitted by-cell tabulations; 8 wording and reporting defects.

---

## 2. (a) Re-derivation

Own code, own sample construction (earliest `usable == 1` session per subject; agrees with the frozen `first_usable` flag for all 80), own plausibility pass (2 envelope hits: `sub-21038_ses-2` chi_SNVTA -53.6 ppb in a non-usable row; `sub-20953_ses-2` chi_DN_L -313.8 ppb -> missing; same as the report), race3 with `contr.sum`, `ns(age, 3)`, covariates z-scored in-sample, HC3 by `(X'X)^-1 X' diag(e^2/(1-h)^2) X (X'X)^-1` (max abs difference from `sandwich::vcovHC` 0), t reference with df = n - p = 80 - 12 = 68.

| Quantity | Review (R, manual HC3) | Python statsmodels | Frozen report | Match |
|---|---|---|---|---|
| beta_PS (ppb) | 4.962802 | 4.962802 | 4.963 | yes |
| SE_HC3 | 3.125080 | 3.125080 | 3.125 | yes |
| t; df; p (two-sided) | 1.5881; 68; 0.116913 | same | 1.588; 68; 0.1169 | yes |
| 95 % / 90 % CI (ppb) | [-1.273, 11.199] / [-0.248, 10.174] | - | same | yes |
| sigma_resid; sigma_Huber (converged); sigma_B | 11.2559; 10.5505; 10.5505 | 11.2559 | same | yes |
| delta = beta/sigma_resid | 0.44091 | - | 0.441 | yes |
| delta_B = 0.6 sqrt(0.93828) | 0.58119 | 0.58119 | 0.581 | yes |
| bound (ppb); t_B; p_B | -6.132; 3.5502; 0.000352 | 0.000352 | -6.132; 3.550; 0.0004 | yes |
| Huber M-estimate (SE) | 3.936 (2.661) | - | 3.936 (2.661) | yes |
| K1 PS-only - TD | 8.1574 (3.9804), p .04428 | same | 8.157 (3.980), .0443 | yes |
| K1 PS+ADHD - ADHD-only | 1.7682 (4.7855), p .7129 | same | same | yes |
| ADHD main; interaction | -1.272 (3.106) p .684; -6.389 (6.199) p .306 | - | same | yes |
| ADHD TOST at +/-5.275 ppb | p .1009 / .0194 -> fail | - | same | yes |
| K2 R2* beta_PS; one-sided p | 1.098 (0.510); .9825 | 1.098 | same | yes |
| K3 beta_dref; p; threshold; trigger | 1.486; .274; 1.241; FALSE | 1.486; .274 | same | yes |
| K3 CSF-referenced beta; p; p_B | 3.477; .351; .0021 | - | same | yes |
| S1 L; SE; p (df 67) | 0.65355; 0.47262; .17131 | same | same | yes |
| S1 beta_PS(D); beta_ADHD(D) | 0.570; -0.084 | - | same | yes |
| S2 n; beta_PRIME; SE; p (df 52) | 65; 3.9682; 2.3682; .09982 | same | same | yes |
| Holm (S1, S2) | .19963 / .19963 | same | .1996 / .1996 | yes |
| Freedman-Lane p (10,000, own code) | 0.0860 | - | 0.0857 [.080, .091] | MC |
| Bootstrap delta 95 % CI (5,000, own seed) | [-0.023, 1.021] | - | [-0.041, 1.019] | MC |
| BF10 attenuated / sceptical | 0.2790 / 1.2727 | - | same | yes |
| PC2x PUT slope, p; GP slope, p | 0.382, .1225; 2.397, .00025 | - | same | yes |
| QC-by-group Fisher (in H1; earliest) | .967; .244 | - | same | yes |
| Automated-exclusions-only M1 | n 90, 5.928, p .0595 | - | same | yes |
| Samples | 181 rows; 107 usable; 80 first-usable 22/16/17/25; S2 65 | 80; 65 | same | yes |

Full numbers: `review_rederivation.json`, `review_check_statsmodels.json`.

---

## 3. (b) Decision logic against the registration

- **Step 1 (sec. 4.1).** A requires beta < 0 and p < .05: beta = +4.96 -> not A. B requires PC4 pass and r_L >= 0.25 (0.938: evaluable) and t_B > t_{.95,68} = 1.668: t_B = 3.55 -> B. Code (`h1_decide`) matches the registered pseudo-code line by line, including `min(1, r_L)`, `min(sigma_resid, rob$s)`, the Huber settings (k 1.345, proposal 2, maxit 200, MAD fallback) and the one-sided t reference.
- **Disconfirmation** needs beta > 0 with two-sided p < .05: p = .117 -> not disconfirming. Correct that no "higher susceptibility" sentence appears.
- **Step 2 (K3).** |beta_dref| = 1.49 > threshold 1.24 but p = .27 >= .10 -> not triggered (the trigger is conjunctive) -> final = M1 outcome = B. Correct. (Had it triggered, the CSF-referenced outcome is also B, so the final would still be B.)
- **K1 precedence** K1c > K1b > K1a: all three require A -> "none". The mandatory co-reports (PS-only - TD; ses-1-only n = 58; threshold psychosis) are present. The TOST and interaction CI are reported as non-decisive. Correct.
- **K2**: beta_PS(R2*) > 0 -> "iron-versus-myelin unresolved". Correct.
- **K4**: fragile = any LOO change of the M1 A/B/C category (0/80); influence-dependent requires the dominance trigger (not met); "robust" requires final A and >= 400 evaluable specifications (324). Correct; the "fragile" and "influence" refits correctly compare M1 outcomes before K3 (sec. 4.1 Step 2 wording).
- **Holm family.** S3 is exploratory (PC5 fail, amendment 1 #39; n_min none, #10b-d, all decided before any chi), so the family is {S1, S2} at .025/.05. Holm: min p .0998 x 2 = .1996; S1 .171 raised to .1996 by monotonicity. Both not rejected; both point estimates are also in the direction opposite to the prediction (L > 0; beta_PRIME > 0), so the direction rule would have blocked "supported" anyway. Correct.
- **T1** is conditional on final A: not evaluated, follow-up estimates descriptive. Correct.
- **Caveats.** PC2w INDETERMINATE -> PC2x after unlock -> putamen fails -> "developmental positive control not confirmed (applies to every H1 statement)". Correctly appended to the registered reading.
- **Registered sentences.** Section 7 of the report prints exactly the sec. 4.6 sentences for B, the dilution caveat, T1, S1, S2, S3 and the branch caveats. No prohibited phrasing ("iron loss emerges with illness", "not a vulnerability marker", "stable", "persistent") appears.
- **Environment pin, gate validation, draws**: all registered (10,000 / 10,000 / 5,000 / 5,000 / 2,000 / 9,999; seed 20260924; R 4.4.1 and package pins match).

---

## 4. (c) Concerns, ranked

### C1. Outcome B is confirmatory only because amendment 1b removed the registered G5 chain trigger (highest consequence)

- Registered rule (Q24 G5): if > 10 % of sessions fail G5, revise once from the closed menu; if > 10 % still fail, **G5 is recorded failed and the methods-report path applies**. `run_all.R` then sets `primary_interpretable = FALSE` and `report.R` (line 172) prints *"Pipeline not validated (positive-control or registration gate failure): H1 is reported descriptively only."*
- Realised: fail share 46 % (v1), the menu revision (1a) made it worse (27 % non-IQ), reverted (1b); fail share after the IQ split 24 % overall / 11.7 % non-IQ. Both exceed 10 %. The gate JSON records: `"rule": "per-session exclusion only (amendment 1b); chain trigger removed; registered rule would have read fail"`, and `gates_preunlock.py` was patched (#46) to write `G5 = TRUE`.
- Judgement: the decision was made before unlock on label-free aggregates, is documented, and its reasoning (MEGRE<->T1w Dice invariant to the registration revision; failures track echo-1 image quality) is plausible. But it is a *relaxation of a gate's consequence*, not a technical fix of a failing step, which is at the edge of what Annex C.3 permits ("a pre-named label-free metric reveals a technical failure"). It does not change any number; it changes whether "B" is a Tier-1 claim or a descriptive one. The frozen report's branch line says "primary interpretable: TRUE" without this history. **The paper must state, next to outcome B, that under the registration as posted H1 would have been descriptive, and why the trigger was removed.**
- The 1b bias safeguards are **absent from the frozen report** (G5-H1 and IQ fail by cell with Fisher; Dice by group with Wilcoxon; IPW if p < .05). Computed here (session level, n = 181 / first-usable n = 80): G5 fail by cell 16 / 32 / 21 / 29 % (TD / ADHD-only / PS-only / PS+ADHD), Fisher p = .26; by ADHD 18 % vs 30 %, p = .08; by PS p = .45; `dice_min` by cell Kruskal p = .17, by PS p = .17, by ADHD p = .80. The IPW trigger is not met, but the tabulation belongs in the report.

### C2. Coding error in a co-report: the residualised split-echo reliability ignores the registered plausibility rule; two {2,4} reconstructions are implausible

- The frozen report prints "residualised split-echo bound p = 0.0263" from `resid_split_reliability()` (m1.R), which correlates the covariate-residualised {1,3} and {2,4} halves over all 80 first-usable rows with non-missing split values. The registration says the plausibility rules "apply to the split-echo, engine and multiverse reconstructions; a failure makes the value missing for that variant only" (sec. 12.4, line 292), and the gate script implements exactly that for PC4 ("whose split values pass their own plausibility rules", n = 78). The analysis-side function has no such filter.
- Effect: two {2,4} values are far outside the [-50, 300] ppb envelope: `sub-20866_ses-1` -236.6 ppb and `sub-20712_ses-1` -189.1 ppb (both PS+ADHD; their primary 4-echo chi is normal: 52.2 and 65.0 ppb; {1,3} halves 52.5 and 68.8; robust z of the full value 0.02 and 1.36; all 9 multiverse reconstructions 35-70 ppb). With them, raw r = 0.218 and residualised r = 0.237 (r_L 0.036, delta_B 0.114, p_B .0263); without them (n = 78), raw r = 0.9245 / R_SB 0.961 / r_L 0.938 (**identical to the gate JSON, so PC4 itself is correct**) and residualised r = 0.922 / R_SB 0.960 / r_L 0.937 -> delta_B(resid) = 0.581, p_B = .00035.
- Consequence for the outcome: none (the co-report is "never decisive"). Consequence for the record: Annex C.4 requires that the coding error be logged (suggested log entry #52), both values shown, and the text-conforming value (p_B = .00035, delta_B 0.581) treated as the co-report. The sentence in the report that the covariate-residualised reliability is low (r_L 0.036) is wrong; after the registered filter the residualised reliability is as high as the raw one.
- Pipeline flag: both failing sessions are on the amendment-1 list of six **unsigned-12-bit phase** sessions (#14; `--phase-mode auto12` patch, #15). Of the four unsigned sessions in H1, the other two (`sub-20823_ses-1`, `sub-21038_ses-1`) have normal split values. The {2,4} echo-subset copies (TE1 = 15 ms, lambda_{2,4} = 2 lambda) should be checked for these two sessions before any split-echo quantity is published. Exploratory M1 without the four unsigned-phase sessions in H1: n = 76, beta = 4.93, p = .125 (no material change).

### C3. What outcome B actually excludes; the bound test is insensitive to the robust-scale choice

- B is driven by the sign of the point estimate, not by the scale. The one-sided test of H0: delta <= 0 (delta_B = 0) gives p = .058; the one-sided 95 % lower bound of beta is -0.25 ppb, i.e. delta >= -0.024 sigma_B (-0.022 sigma_resid). **Any group-level deficit larger than about 0.02 residual SD is excluded at one-sided 5 % on the observed scale.** The registered B sentence (adult magnitude excluded) is therefore correct and conservative.
- Scale sensitivity grid (`bound_test_sensitivity_grid.csv`; 6 scales x 5 SE estimators x 10 bounds = 300 cells): B holds in 294/300; the only non-B cells are delta_B = 0 with the HC3 SE. sigma_B would have to fall below **0.43 ppb** for B to fail at the registered bound. Scales: sigma_resid 11.26, Huber 10.55, 1.4826 MAD 9.12, pooled total SD 11.88, raw SD in TD 10.9, raw SD overall 11.9. SEs: classical 2.74, HC0 2.62, HC1 2.85, HC2 2.86, HC3 3.13 (HC3 is the most conservative here, for A and for B).
- Interpretive caveat the report already carries: delta_B uses split-echo reliability (an upper bound on reliability, r_L 0.938 -> delta_B 0.581), so B is an observed-scale statement; the cross-session co-report (r_L 0.54, delta_B 0.44, p .0015) is the more defensible reliability scaling. All co-reports agree.
- The direction (delta = +0.44, bootstrap CI [-0.04, 1.02]; PS-only - TD +8.2 ppb; R2* also higher in PS, two-sided p .035) should be reported descriptively with CIs and without any "higher iron" wording, since the disconfirmation rule (p < .05 two-sided on beta_PS) is not met.

### C4. Multiplicity and direction of the co-reported "PS-only - TD = +8.2 ppb, p = .044"

- It is the only p < .05 in the H1 section and is quoted next to the outcome in log #51. It is a post hoc linear contrast of the registered model, in the direction **opposite** to H1, and registered as "never decisive, no label".
- Adjusted values: Holm over the two simple effects .089; Holm over the four K1 contrasts .177; Holm with beta_PS in the family .133; Freedman-Lane max-T permutation over the two simple effects (PS labels permuted within ADHD strata, 10,000 draws) **p = .060** (unadjusted permutation p .034). Under no handling is it significant. Raw cell means of usable sessions: TD 51.5, ADHD-only 52.2, PS-only 59.0, PS+ADHD 54.2 ppb.
- The frozen report labels it correctly. The risk is downstream: the paper and abstract must not present it as a finding.

### C5. n shortfall and realised power (the registered power statements are now wrong)

- Realised n: H1 80 (registered 112, projected 100-106); S2 65 (83); T1 returners 26 (53); ses-1-only 58 (95). Loss: visual grade 2 on 62 sessions (33 %, AI rater; registered expectation 5-10 %), G5 on 45, T1w 6, auto grade 2 on 7, plausibility 8.
- Realised SE(delta) = 0.278 (HC3, residual-SD units; 0.296 in sigma_B units) versus 0.206 registered at n = 112 and **0.24 assumed in amendment 2**. Post hoc, noncentral t with df 68: power for A at delta = -0.3 / -0.4 / -0.5 / -0.6 = 0.19 / 0.30 / 0.43 / **0.57** (amendment 2 says "approximately 0.70"); MDE80 for A = 0.79 sigma; power for B at the registered bound when the true delta = 0 / -0.1 / -0.2 / -0.3 = 0.62 / 0.49 / 0.36 / 0.24; false B at the bound 0.05; interaction MDE80 = 1.57 sigma. The manuscript's power sentence should carry these values.
- Differential loss of the earliest session: 10/30 TD, 12/22 ADHD-only, 8/24 PS-only, 19/37 PS+ADHD (Fisher by cell p .24; by ADHD 33 % vs 53 %, p = .057). Consequently the ADHD cells are analysed more often at ses-2/3 (ADHD-only 6/16, PS+ADHD 11/25 vs TD 3/22, PS-only 2/17), i.e. older and later in the protocol. The registered ses-1-only co-report (beta 6.5, p .18, B) and the epoch sensitivities address this; it should be said explicitly.

### C6. PC2x partial failure (putamen) and the age adjustment

- GP: 2.40 ppb/y (p .0003), consistent with the registered expectation (about 2.3-3.5 ppb/y). Putamen: 0.38 ppb/y (p .12) against an expected about 1.3-1.7 ppb/y; on all 107 usable sessions 0.03 ppb/y (p .46); with sex and cell 0.50 (p .07); TD only 0.69 (p .07); Spearman .11. Together with the indeterminate PC2w (returner change p .16, n = 26 pairs), the pipeline's sensitivity to known putaminal iron accrual is genuinely unconfirmed; the caveat is correctly attached, and the paper should discuss candidate causes (whole-support reference, 1.5 mm partial volume in a thin lentiform structure, SynthSeg erosion) rather than treat it as a technicality.
- Age is a real confounder here: SN-VTA chi rises with age (Spearman .23, p .02), PS youths are about 1 y older (r(age, PS) = .27; age by PS p = .023). beta_PS is 5.72 (p .036) with no age term, 4.72 (p .14) with linear age, 4.96 (p .117) with the registered ns(age, 3). The registered adjustment is doing legitimate work; an unadjusted comparison would have looked "significant" in the wrong direction.

### C7. Selection effects of the AI visual QC; omitted by-cell tabulations

- Exclusion is balanced across cells: grade-2 share 33 / 34 / 27 / 37 % (Fisher p .82); usable share by cell p .55; subject in H1 by cell p .97 (reproduced).
- Exclusion is strongly outcome-related, as it should be: mean chi_SNVTA by grade 0 / 1 / 2 = 57.5 / 53.8 / 41.0 ppb (SD 10.0 / 12.3 / 44.7; Kruskal p = 2e-5); by G5 fail p .003; by usable p 4e-5. Non-usable sessions in the ADHD cells have very low chi (ADHD-only 32.9, PS+ADHD 38.9 ppb) versus TD 51.3 and PS-only 56.8. With balanced exclusion rates the bias risk to beta_PS is modest, and the sign of beta_PS is stable across QC stringency: no QC at all (earliest plausible session, n = 113) +9.55 (p .075); automated exclusions only (n = 90) +5.93 (p .059); registered (n = 80) +4.96 (p .117); stricter sets +3.0 to +4.6. **Stricter QC moves the estimate toward zero; the AI grading did not manufacture the positive direction.**
- Caveats that stand: the AI rater is unvalidated against expert QC (declared in 1c); pass-1/pass-2 kappa .66 is adequate, not high; 7 grade-0 sessions make the "strict QC" sensitivity (n = 4) and the grade-0 specification level unusable.
- The Q13 by-cell tabulation (orientation angles, motion index, in-support fractions, slab coverage, reference-domain fractions, code D, CSF validity) is **missing from the frozen report**. Computed here (`covariates_by_cell_first_usable.csv`): nothing differs by cell (all p > .12) except the CSF-minus-whole-brain reference offset `ref2_minus_wb_ppb`, larger in the ADHD cells (ADHD-only 21.9, PS+ADHD 23.8 vs TD 18.8, PS-only 18.5 ppb; Kruskal p .017; by ADHD p .004; by PS p .063). K3 (PS, not triggered, p .27) covers H1, but the ADHD-related reference offset should be reported next to the ADHD main effect and the TOST. Code D by cell 6 / 2 / 2 / 3 sessions.

### C8. Wording and reporting defects in the frozen report (none changes a decision)

1. T1 table: "occasion-0 beta_PS 0.418: **the occasion-0 deficit exceeds Delta**". The estimate is +0.418 total SD, an excess, not a deficit; `t1.R` lines 166-168 compare `abs(occ0)` to Delta. Misleading wording; T1 was not evaluated anyway.
2. "Mandatory caveat: all 6 threshold-psychosis youths sit in PS+ADHD" is a statement about the cohort; **only 3 are in the analysed sample** (sensitivity row n = 77 = 80 - 3, yet labelled "without the 6"). Say "3 of the 6 (all PS+ADHD)".
3. Log #49 says the strict-QC row had n = 7; the report says n = 4 (first-usable grade-0 sessions). The report is right; the log should be corrected.
4. The specification-curve line prints Stouffer Z = 19.96 over 324 highly dependent specifications; as evidence it is meaningless (its own joint permutation p is .93). The K4 "robust" label was **unassessable by construction**: the grade-0 sample level (n = 4) makes 162 of 486 specifications non-evaluable, so at most 324 < 400 could ever be evaluable. A registered design flaw; moot because A does not hold, but it should be stated rather than reported as "not assessable" without reason.
5. The B co-report line prints the contaminated residualised split-echo p (.0263) (C2).
6. K2 shows only the one-sided p (.9825); the paper should show the R2* estimate with its CI (1.10 s^-1 [0.08, 2.12]; two-sided p .035) descriptively.
7. The branch line "primary interpretable: TRUE" should carry the amendment-1b history (C1).
8. `verify_statsmodels.py` is a same-AI re-implementation (declared); the present review adds a second independent implementation but still not a human analyst.
9. Minor consistency: PC4 was computed on n = 78 (gate) while the residualised co-report used n = 80 (C2); the report should give n for both. The sample table's "182-row" phrase in log #47 counts the header; the table has 181 data rows.

### Minor statistical observations (no action needed)

- HC3 SE (3.13) exceeds the classical SE (2.74): the registered choice is conservative for both A and B in this sample. The Freedman-Lane p (.086) sits between the HC3 and classical p values.
- The Huber scale (10.55) is below the OLS sigma (11.26) because of a few large residuals (max |externally studentised residual| 4.05, `sub-23698_ses-1`, chi 12.2 ppb, ADHD-only, the lowest value in the sample; top-1 residual share 16 %); the registered sigma_B = min rule makes B slightly harder, as intended; LOO shows no category change.
- SD / robust SD for chi_SNVTA is 1.24 (< 1.5), so no robust-scale co-report was due for S1.

---

## 5. (d) Post-unlock deviations (amendment1_log #47-51) and whether any registered number could have been affected

| # | What happened | Could it affect a registered number? | Judgement |
|---|---|---|---|
| 47 | Unlock by the user (`unlock_and_merge.py`); 6 registered hashes matched; 4 branch tables merged; no DWI S3 table | No | Not a deviation. |
| 48 | First `run_all.R` halted in `run_s3` (empty S3 sample, `ns()` on an empty vector). `secondary.R` patched: early return with a structured "not fitted" result when the S3 sample has < 10 rows; S3 influence refit skipped when no model exists. Registered copy kept (`78ee1dea...`), patched `086d5c2b...` | **No.** Diff verified here: only the two guards were added; `confirmatory` is computed from the same gates (PC5 FALSE), so the Holm family is {S1, S2} either way; H1, S1, S2, T1 code untouched. The halted run wrote and printed nothing beyond cell counts (log), so no outcome was seen before the fix | Legitimate Annex C.4 coding fix; correctly logged. |
| 49 | Strict-QC sensitivity (grade 0 only) rank-deficient and skipped by `try_fit` | No | Registered consequence; the log's "n = 7" should read n = 4. |
| 50 | `verify_statsmodels.py` patched to skip S3 when not fitted; 11/11 checks within tolerance | No (verification only) | Legitimate C.4 fix. Same-AI limitation declared. This review's second-engine check: 24/24 within 1e-13. |
| 51 | Result recorded | No | Consistent with the report. Two wording points: quotes the PS-only - TD p = .044 next to the outcome (C4); "6 threshold-psychosis youths" (3 in the sample, C8.2). |
| (new) | Coding error found by this review: `resid_split_reliability()` applies no split-value plausibility rule (C2) | Only a co-report (never decisive) | Should be logged as #52 under Annex C.4 with both values; the text-conforming value is p_B = .00035 (delta_B 0.581, n = 78). |

Pre-unlock items that bear on the confirmatory numbers, for completeness (all label-free and declared before unlock; none altered an H1 number after unlock): amendment 1a/1b (G5 menu revision and revert; chain trigger removed -> C1; 45 sessions excluded by G5), amendment 1c (AI rater; 62 grade-2 sessions -> C5, C7), log #46 (gate JSON patched to `G5 = TRUE`), amendment 1 #10d and #39 (S3 exploratory), #15 (`--phase-mode auto12`; the two implausible {2,4} split values are unsigned-phase sessions -> C2).

---

## 6. What the manuscript should do (actionable)

1. Report outcome B with the registered sentence and the dilution caveat, **and** state that under the registration as posted the G5 trigger would have made H1 descriptive, with the amendment-1b rationale and the by-cell G5 / Dice tabulations (C1).
2. Log the split-echo coding error (#52), show both residualised-reliability values, use the text-conforming one; check the {2,4} echo-subset copies of the unsigned-phase sessions (C2).
3. State what is excluded on the observed scale (any deficit > about 0.02 residual SD, one-sided 95 %) alongside the adult-bound wording; present the positive direction descriptively with CIs, no "higher iron" claim (C3).
4. Keep PS-only - TD as a labelled co-report with its adjusted p values (.060 max-T; .089 Holm); never in the abstract (C4).
5. Replace the amendment-2 power sentence with the realised values (power 0.57 for an adult-size deficit; MDE80 0.79 SD; B power 0.62 at true delta 0) and report the differential loss of earliest sessions in the ADHD cells (C5).
6. Discuss the putamen PC2x failure as a measurement limitation; keep the caveat on every H1 statement (C6).
7. Add the Q13 by-cell tabulation; report the ADHD-related CSF-reference offset (C7).
8. Correct the T1 wording, the "6 threshold-psychosis" count, the strict-QC n in the log, and explain why "robust" was unassessable (C8).

---

## 7. Files

- `REVIEW.md` (this document); `README.md`
- `review_rederivation.json`: all re-derived quantities, comparisons with the frozen results, bound grid summary, K1 multiplicity, power, QC-by-cell, split-echo decomposition, misc checks
- `review_check_statsmodels.json`: second-engine check (24 quantities, all within 1e-6; actual max 1e-13)
- `bound_test_sensitivity_grid.csv`: 300-cell scale x SE x bound grid for the B test
- `covariates_by_cell_first_usable.csv`: the registered Q13 by-cell tabulation (medians, Kruskal, Wilcoxon by PS and by ADHD)
- `fail_reasons_by_cell_sessions.csv`: fail-reason tokens by cell over all 181 sessions
- `split_echo_extreme_sessions.csv`: the sessions driving the split-echo discrepancy
- Code: `results_exploratory/code/review_rederive.R`, `review_check_statsmodels.py`, `review_split_echo.R`
