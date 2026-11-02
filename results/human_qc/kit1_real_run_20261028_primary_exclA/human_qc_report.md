# LEAD-chi: blinded human visual-QC validation of the AI rater (POST HOC)

Status: **post hoc / exploratory**. Not part of the registered confirmatory analysis; H1 is not re-decided.

- ratings file: `/path/to/Penn_LEAD/LEADchi/human_qc_kit/returned/LEADchi_humanQC_rater_20261028T1724Z.csv` (sha256 `1975118e501066fd87cb7e18a096f537a034caf3b38897590a927649cc2118c3`)
- kit key sha256: `a4bb713258897af408097bfc4523d8d7deadbf68c081b51e9bb15c5e8ea1d97b`; run (UTC): 2026-10-28T17:29:33Z; bootstrap 2000 (seed 20261008)
- input issues: {"file_encoding": "utf-8-sig", "file_delimiter": "comma", "unknown_kit_ids": 0, "duplicated_rows": 0, "invalid_failure_codes": 0, "missing_or_invalid_grade": 0, "missing_or_invalid_placement": 0, "excluded_from_part_A": 1}

## A. Agreement, human vs AI (primary items)

Rated primary items: 100 ({'AI_fail_all': 60, 'AI_pass_sample': 40}).

Confusion matrix (rows human, columns AI):

| | AI 0 | AI 1 | AI 2 |
|---|---|---|---|
| human 0 | 1 | 38 | 44 |
| human 1 | 0 | 1 | 7 |
| human 2 | 0 | 0 | 9 |

| estimand | linear-weighted kappa [95% CI] | binary fail-vs-pass kappa [95% CI] | exact agreement |
|---|---|---|---|
| sample_unweighted | 0.054 [0.022, 0.093] | 0.124 [0.054, 0.211] | 0.110 |
| design_weighted_to_180 | 0.057 [0.025, 0.096] | 0.190 [0.087, 0.288] | 0.084 |
| flag0_only_sample (n 79) | 0.068 [0.027, 0.121] | - | 0.114 |
| flag0_only_design_weighted (n 79) | 0.060 [0.024, 0.103] | 0.240 | 0.083 |

Pre-specified (README_coordinator section 7, 2026-10-07): the design-weighted linear kappa (primary) and the design-weighted linear kappa on forced-review flag = 0 items are reported with equal prominence; the flag is shown to both raters and is a strong cue for the AI grade. Kappa is NA when a rater is constant.

Design-weighted (to the 180 graded sessions), human as reference for fail = grade 2: AI sensitivity 1.000 [1.000, 1.000], specificity 0.700 [0.680, 0.726], PPV 0.150, NPV 1.000; estimated human fails among 180: 9.0.

- AI-fail stratum (n = 60; kappa undefined, AI constant): human confirms grade 2 in 0.150 [Wilson 0.081, 0.261]; human grades {'0': 44, '1': 7, '2': 9}.
- AI-pass stratum (n = 40): linear kappa 0.001; exact agreement 0.050; human grade 2 (AI missed fail) 0.000 [Wilson 0.000, 0.088], i.e. about 0.0 of the 119 AI-pass sessions.
- By forced-review flag: flag 0: n 79 (40 AI grade 2), linear kappa 0.068 (design-weighted 0.060); flag 1: n 21 (20 AI grade 2), reported as a confirmation rate only: P(human = 2 | AI = 2) 0.100 [Wilson 0.028, 0.301], P(human = 2) 0.095; human grades of the flagged AI-pass items {'0': 1, '1': 0, '2': 0}
- Placement: human fail 4 of 100 items, AI fail 0; agreement 0.960; kappa NA.
- Failure-code presence (both / human only / AI only; kappa): M 14/2/27 (0.34); W 0/0/21 (NA); S 0/0/8 (NA); B 0/0/24 (NA); K 0/1/5 (-0.02); F 0/0/0 (NA); A 0/0/4 (NA); C 0/0/0 (NA); D 0/0/10 (NA)

- Secondary comparator, AI pass-1 grade (the final AI grade is the worse of pass 1 and pass 2 for the sessions the AI rated twice; the human looks once): AI pass-1 grades {'0': 1, '1': 41, '2': 58}; linear kappa 0.057 (design-weighted 0.058 [0.027, 0.098]); exact agreement 0.110.

- Excluded from Part A (pre-specified rule): ['hfaec73'].

## Intra-rater reliability (hidden duplicates)

10 of 10 pairs rated twice: exact grade agreement 0.700, linear-weighted kappa 0.318, placement agreement 0.900, identical code sets (pairs with both grades) 0.800. 10 pairs: kappa is imprecise; exact agreement is the main intra-rater figure. Repeats were 21-70 items apart (about 15-50 min), so recall is possible and the figure may be optimistic.

For comparison, the AI's own repeatability: 31 of 36 re-rated sessions got the same grade in pass 1 and pass 2 (exact 0.861, linear kappa 0.664).

## B. POST HOC sensitivity: first-usable sample re-derived with the human grade, M1 refitted

Reconstruction of the frozen first-usable list from fail_reason + AI grades: {'usable_equal': True, 'first_usable_equal': True, 'fail_reason_equal': True}.
Sessions in kit: 101; with a human grade: 101.

Refit reproduces the frozen confirmatory estimate under AI_frozen: **True**.

| scenario | n | cells TD/ADHD/PS/PS+ADHD | beta_PS ppb (HC3 SE) | 95% CI | p two-sided | bound test at -6.1318 ppb: t_B, p_B | re-derived bound (ppb): p_B | A | descriptive label |
|---|---|---|---|---|---|---|---|---|---|
| AI_frozen | 80 | 22/16/17/25 | 4.96 (3.13) | [-1.27, 11.20] | 0.1169 | 3.55, 0.00035 | -6.13: 0.00035 | False | B |
| H_replace | 89 | 25/18/19/27 | 6.29 (3.07) | [0.18, 12.40] | 0.0439 | 4.05, 0.00006 | -6.76: 0.00003 | False | B |
| H_grade_only | 89 | 25/18/19/27 | 6.29 (3.07) | [0.18, 12.40] | 0.0439 | 4.05, 0.00006 | -6.76: 0.00003 | False | B |
| Strict_union | 80 | 22/16/17/25 | 4.96 (3.13) | [-1.27, 11.20] | 0.1169 | 3.55, 0.00035 | -6.13: 0.00035 | False | B |
| Lenient | 89 | 25/18/19/27 | 6.29 (3.07) | [0.18, 12.40] | 0.0439 | 4.05, 0.00006 | -6.76: 0.00003 | False | B |

Sample changes relative to the frozen first-usable list:

- AI_frozen: usable sessions 107, first-usable 80; subjects added 0, removed 0, first session changed 0; post hoc plausibility re-check exclusions 0.
- H_replace: usable sessions 126, first-usable 89; subjects added 9, removed 0, first session changed 2; post hoc plausibility re-check exclusions 0.
- H_grade_only: usable sessions 127, first-usable 89; subjects added 9, removed 0, first session changed 2; post hoc plausibility re-check exclusions 0.
- Strict_union: usable sessions 107, first-usable 80; subjects added 0, removed 0, first session changed 0; post hoc plausibility re-check exclusions 0.
- Lenient: usable sessions 127, first-usable 89; subjects added 9, removed 0, first session changed 2; post hoc plausibility re-check exclusions 0.

Labels in the last column are the registered decision arithmetic applied descriptively to a post hoc sample; they do not replace the registered outcome (B).

**Coverage (pre-specified wording, README_coordinator section 7.3).** 52 of the 80 frozen first-usable sessions and 79 of the 119 AI-pass sessions are not in the kit and keep the AI grade in every scenario above. A B label under H_replace therefore means 'robust to replacing the AI grade by the human grade for the re-rated sessions', not 'independent of the visual-QC rater'. Before any rating, a sweep of 300 random human gradings of the kit sessions gave label B in 300 of 300, so the main value of this exercise is the validation of the AI (Part A).

**Extrapolation (POST HOC, pre-specified 2026-10-07).** Starting from H_replace, human fails were imputed among the 79 AI-pass sessions outside the kit, at a rate drawn per replicate from Beta(m + 1/2, n - m + 1/2) with m = 0 of n = 40 rated AI-pass kit sessions failed by the human; 1000 of 1000 draws fitted (seed 20261009). Label counts {'B': 1000}; share B 1.000; share A 0.000; imputed fails (min, 2.5%, median, 97.5%, max) ['0.0', '0.0', '0.0', '6.0', '11.0']; n ['84', '86', '89', '89', '89']; beta_PS ppb ['4.53', '5.80', '6.29', '6.82', '8.25']; p_B at the fixed bound ['0.00002', '0.00004', '0.00006', '0.00015', '0.00045'].

## Method notes

- Kappa: weighted Cohen kappa, linear weights (registered G4 statistic), same formula as `qc_montage.py`; reported as NA when either rater uses a single category (post hoc tooling fix, review item M1).
- Design weights: AI grade-2 sessions are a census (weight 1); the 40 AI grade-0/1 sessions represent 119 (weight 119/40). Bootstrap resamples within the two strata (percentile 95% CI).
- Human per-session rating = worse of the original and its hidden duplicate (registered conservative rule).
- Part B keeps every non-visual exclusion component of the frozen `fail_reason` unchanged, including the automatic grade 2 (`c_auto_grade2`), which no rater can override (rubric section 2).
- M1 refit sources the frozen R helpers; bound test at the frozen registered bound -6.1318 ppb (t_B = (beta + 6.1318)/SE, one-sided upper p, df = n - p) and at the bound re-derived in each sample.
