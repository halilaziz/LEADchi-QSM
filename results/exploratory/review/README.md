# results_exploratory/review

**EXPLORATORY / REVIEW (Tier 3-4).** Adversarial statistical review of the frozen confirmatory report
(`mirror/results/confirmatory/report.md`, `results.json`). Nothing here re-decides a registered outcome; the
confirmatory result (H1 final outcome B) is frozen and is not restated differently anywhere in this folder.

| File | Content |
|---|---|
| `REVIEW.md` | The review: verdict, re-derivation table, decision-logic audit, ranked concerns, post-unlock deviation audit, actions |
| `review_rederivation.json` | Every re-derived quantity (own R code, HC3 by explicit algebra) with the frozen value and relative difference; bound-grid summary; K1 multiplicity (Holm, max-T permutation); realised power; PC2x variants; split-echo decomposition; QC and covariate tabulations by cell; misc report checks |
| `review_check_statsmodels.json` | Second-engine check (Python statsmodels, own sample rules and spline basis): 24 OLS quantities, all within 1e-13 of the frozen values |
| `bound_test_sensitivity_grid.csv` | Outcome-B bound test over 6 residual scales x 5 SE estimators x 10 bounds (300 cells) |
| `covariates_by_cell_first_usable.csv` | The registered Q13 by-cell tabulation the frozen report omitted (Dice, motion, Q_phase, FOV, in-support, slab coverage, orientation, reference fractions, ROI volume, age, CSF-reference offset) |
| `fail_reasons_by_cell_sessions.csv` | First-usable fail-reason tokens by cell, all 181 sessions |
| `split_echo_extreme_sessions.csv` | Sessions driving the split-echo reliability discrepancy (two implausible {2,4} values) |

Code: `results_exploratory/code/review_rederive.R`, `review_check_statsmodels.py`, `review_split_echo.R`
(R 4.4.1; python3 statsmodels). Inputs: `mirror/unlocked/merged_qsmxt_primary_SNVTA.csv`,
`mirror/gates/gates_preunlock.json`, `mirror/results/confirmatory/results.json`, `mirror/unlocked/multiverse_qsmxt_SNVTA.csv`.

Headline: all decision quantities reproduce exactly and the decision follows the registered (and pre-unlock-amended)
rules. Ranked concerns: (1) B is Tier-1 only because amendment 1b removed the registered G5 chain trigger; (2) the
residualised split-echo co-report is contaminated by two implausible {2,4} reconstructions (coding error, Annex C.4;
co-report only); (3) B is driven by the positive point estimate and is insensitive to the robust-scale choice;
(4) the PS-only - TD p = .044 does not survive multiplicity (max-T p .060); (5) realised power 0.57 for an adult-size
deficit, not 0.70; (6) putamen PC2x fails; (7) AI-QC exclusion is balanced across cells and stricter QC moves the
estimate toward zero; (8) several wording defects in the report.
