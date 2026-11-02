# LEAD-χ confirmatory analysis code (registered attachment)

**Status.** Final version for the embargoed OSF secondary-data registration (task F, 2026-09-26). It implements the
registered text `LEADchi_prereg_v1.0_FINAL.md`, which governs: a deviation of this code from the text is a coding
error handled under Annex C.4, and Annex H of the registration maps every registered rule to the function below.
It was tested **only on synthetic outcomes** simulated on the real Penn LEAD design.

**No susceptibility, R2\* or other MEGRE-derived value exists or was used.** The code never reads images. Design
information comes from metadata and phenotype copies; `megre_coverage.tsv` is not read.

`run_all.R` runs unchanged on the real merged tables after unlock. The whole folder is hashed (`SHA256SUMS`).

## 1. Quick start

```bash
# the confirmatory run (after unlock; tables written by unlock_and_merge.py)
Rscript run_all.R --input unlocked/merged_<engine>_<pipeline>_<roi_def>.csv --gates gates_preunlock.json \
                  --multiverse unlocked/multiverse_<engine>_<roi_def>.csv --out results/confirmatory \
                  [--other-engine unlocked/merged_<other engine>_<pipeline>_<roi_def>.csv]
python3 verify_statsmodels.py --input <same table> --gates gates_preunlock.json \
                  --results results/confirmatory/results.json --out results/confirmatory/verify_statsmodels.json
# rebuild the whole synthetic demonstration (about 20-30 min on the Mac)
./run_demo.sh
```

- `--n-perm`, `--n-perm-spec`, `--n-boot`, `--n-pb` and `--n-wild` exist **for testing only**; using any of them
  stamps the report "NON-REGISTERED counts".
- Outputs: `report.md`, `results.json` (every number), `manifest.json` (hashes of inputs, code and outputs),
  `sessionInfo.txt`.
- **Environment pin (Q23).** R 4.4.1 with sandwich 3.1-1, lmtest 0.9-40, lme4 1.1-38, nlme 3.1-164, MASS 7.3-61,
  boot 1.3-30, jsonlite 1.8.9, digest 0.6.37, Matrix 1.7-0 (`REG$pinned`). **On real data any version difference
  stops the run**; on synthetic data it is recorded. renv is not used; pbkrtest, clubSandwich, lmerTest, MBESS, mice,
  brms, irr, specr and geepack are not installed and are never called (the registered fallbacks are the only branch).
  The verification uses Python 3.11 with numpy, pandas, scipy and statsmodels 0.14.5.

## 2. Files

| File | Role |
|---|---|
| `run_all.R` | Single entry point: validates the tables and the gate JSON (hard stops), checks the environment pin, branches on the gates, re-verifies plausibility, builds the samples and runs everything in the registered order |
| `R/constants.R` | Every registered constant (α, bounds, floors, margins, envelopes, draw counts, seed, fixed ID lists, D12, the environment pin, the required gate fields) |
| `R/utils.R` | CLI parsing, log, SHA-256, Monte Carlo precision interval, Holm |
| `R/design.R` | Schema validation; registered coding; samples (first usable, S2, S3, T1, ses-1-only, exposure sets for H1 and S3) |
| `R/plausibility.R` | ROI envelopes, voxel-rule columns, R2\* envelope, CSF-validity rule, robust z, dominance diagnostics, Huber M-estimator with the registered fallback |
| `R/m1.R` | M1 with HC3; decision A/B/C with σ_B and the B-evaluability rule; registered sentences; Freedman–Lane; bootstrap δ CI; Bayes factors; noncentral-t CI; K1–K3; B co-reports; K4 fragility and influence; ses-1-only and exposure sensitivities; subject-level QC-by-group and IPW; PC2x |
| `R/secondary.R` | S1, S2 (+ its sensitivities and the all-session model), S3 (+ its sensitivities incl. the exposure set), Holm with the direction requirement, influence refits |
| `R/t1.R` | REML-maximum helper; T1 decision by the containment-df test (PB co-reported); OLS cluster-robust (HC2) co-report and wild cluster bootstrap; pooled all-session sensitivity |
| `R/sensitivity.R` | The registered H1 sensitivity battery (sec. 24.4), the PDS score and the single stochastic imputation (mice fallback) |
| `R/spec_curve.R` | The 486-specification QSM curve and the joint Freedman–Lane permutation test behind the K4 "robust" label |
| `R/report.R` | Markdown report; the registered reading prints registered sentences only |
| `schema/roi_table_schema.csv`, `schema/multiverse_schema.csv` | Machine-readable input schemas |
| `verify_statsmodels.py` | Python re-implementation of the decision estimates (statsmodels OLS/HC3; T1 by a direct numpy/scipy REML profile; tolerance 1e-6 OLS, 1e-4 mixed) |
| `R/synthetic.R`, `make_synthetic.R` | **Test only.** Real design from metadata plus simulated outcomes (ROI, multiverse and gate files) |
| `sim_typeI_power.R` | **Test only; authoritative source of the Q21 power figures.** Operating characteristics of the registered functions |
| `run_demo.sh` | Rebuilds `synthetic/`, `demo_output/` and `sim/`, runs the verification and the tests, writes `SHA256SUMS` |
| `tests/test_guards.R` | 16 refusal tests and 4 edge-case tests |

## 3. Inputs

**ROI table** (one row per MEGRE session, usable or not): `schema/roi_table_schema.csv` (95 columns; "yes" = must
exist; "s3" = must exist, may be missing where `s3_eligible = 0`; "optional" = used when present, otherwise the
dependent sensitivity is reported as "not run"). Key rules enforced at load: one engine/pipeline/roi_def per table;
no duplicated sessions; `ref2_type` ∈ {CSF, CSF_erode1, none} with a missing CSF value for `none` (no PLIC); usable
rows have χ, roi_vol and Q_phase; `ses_order` chronological; a frozen `first_usable` flag, if present, must equal
the earliest usable session.

**Multiverse table**: `schema/multiverse_schema.csv` (one row per session × background-removal × inversion node,
primary and "other" ROI values on the 9-node common support, whole-brain and CSF references).

**Gate JSON** (`gates_preunlock.py --stage final`, aggregates only). Required fields (hard stop if absent or
invalid): `data_origin`, `A2_approved`, `PC1`, `PC2w` (PASS / INDETERMINATE / FAIL), `PC3`, `PC4`, `PC4_DN`, `PC5`,
`G2`, `G2_DN`, `G3`, `G4`, `G5`, `rL_primary` ∈ [−1, 1], `r_point_primary`, `constants_roi` (must equal the ROI G3
selects), and `pooled` with `mean_SNVTA`, `sd_SNVTA`, `mean_STR`, `sd_STR`, `sd_total_SNVTA` (> 0). `fallback_PC_pass`
is required when a positive control failed. Optional: `rL_cross` (B co-report), `icc_cross_session`, `sd_dref_ppb`.

## 4. What is implemented

Every registered rule of Parts 1–5 and Annex A that is decided after unlock; Annex H of the registration gives the
full rule → function crosswalk. Registered implementation rules (formerly the "P5 proposals") are now part of the
registered text (Annex I): the T1 decision is the containment-df test; occasion 1 is the earliest usable session after
the first usable one; the exposure sets drop rows without substitution (H1, S3, T1); QC-by-group and IPW are
subject-level; the S3b bootstrap is cell-stratified; the T1 scale is the SD over the first-usable sessions; Huber
k = 1.345, maxit 200, MAD fallback; the stacked GLS structure of S1; the direction requirement in Holm; synthetic runs
mask the exposure sessions' cells.

**Amendment-2 modules** (to be written and hashed in amendment 2, before unlock, implementing the registered text
verbatim; the diff is posted): E-dev with its realised-design calibration, E-z, E-anchor, X1–X7, the node-level
Freedman–Lane max-T (Figure 5d), the GAMM display and the DWI multiverse (sec. 24.5). They are Tier 3–4 and carry no
claim.

## 5. Synthetic demonstration (`demo_output/`)

Real metadata (181 MEGRE sessions, 113 subjects; group, dx, sex, race, age, PRIME/BDI/ARI, Tanner, CPT, substance,
epoch from ProtocolName, DWI availability, round-1 AFQ list) with simulated outcomes: 5 % label-free QC loss, subject
ICC 0.7, correlated SN–VTA/STR/R2\* noise, age trends, simulated DWI and multiverse values; two planted ROI-level
plausibility events per scenario. All runs used the registered resampling counts.

| | `null` | `delta-0.6` | `null_dominant` |
|---|---|---|---|
| H1 n; β_PS ppb [95% CI]; δ | 108; 2.90 [−5.62, 11.41]; 0.15 | 108; −9.10 [−17.62, −0.59]; −0.47 | 108; 5.35 [−4.28, 14.98]; 0.21 |
| HC3 p; Freedman–Lane p | .50; .47 | .036; .029 | .27; .25 |
| σ̂ / σ_rob → σ_B; p_B (bound −0.50) | 19.2 / 20.6 → 19.2; .002 | 19.2 / 20.6 → 19.2; .45 | **25.9 / 21.5 → 21.5**; .0006 |
| K3 trigger; final H1 | no; **B** | no; **A** | no; **B** |
| K1 label; ADHD TOST | none; fail | **K1a PS effect**; fail | none; fail |
| dominance; influence-dependent; fragile | no; no; no | no; no; **yes** (11/108) | **yes**; no; no |
| specification curve: joint p; robust | .80; no | **.0006; yes** | .95; no |
| S3a p (Holm); NC1 / NC2 lower bounds | .23 (.46) | **.0006 (.002)**; 0.28 / 0.13 → specific | .23 (.62) |
| T1 label | not evaluated | indeterminate | not evaluated |
| statsmodels verification | 13/13 | 13/13 | 13/13 |

**Tests** (`tests/test_guards.R`, 20/20): run_all.R refuses a wrong engine or ROI for the gates, a mismatch of the
gate constants' ROI, an origin mismatch, a missing r_L, r_L > 1, a missing A2 field, a zero pooled SD, a failed
positive control without the fallback flag, an SN–VTA value outside its envelope in a usable row, a frozen
first-usable mismatch, a missing column, a duplicated session, a usable row without an outcome, a PLIC reference,
and real data without a multiverse table; r_L < 0 and a PC4 failure make B non-evaluable without stopping; rlm
non-convergence uses the MAD fallback; an HC3 leverage of 1 gives an NA SE without an error.

## 6. Monte Carlo (`sim/typeI_power.md`; authoritative for Q21)

Real pre-QC design (H1 n = 112; S2 n = 83; S3 n = 108; T1 165 rows, 53 returners); registered functions only.
Key results: HC3 type I 0.038 (normal) / 0.043 (skewed, heteroscedastic); P(A) 0.26 / 0.44 / 0.63 / 0.80 at δ =
−0.3 / −0.4 / −0.502 / −0.6; P(B) 0.740 at δ = 0 and 0.038 at the bound (0.052 with two ±6 SD scans); TOST 0.47;
"PS-specific" 0.36 at δ = −0.6; K3 trigger 0.004 / 0.07 at SD(Δref) 0.2 / 0.5 residual SD with no true drift;
Freedman–Lane type I 0.054; S1 type I 0.015 at .0167; S2 power 0.38 / 0.45 / 0.56 (r 0.25) and 0.56 / 0.62 / 0.73
(r 0.30) at .0167 / .025 / .05; S3a power 0.38 / 0.46 / 0.58 (ρ 0.4) and 0.63 / 0.70 / 0.80 (ρ 0.5) at R²_x = 0,
lower when dentate χ shares variance with χ_DGM (ρ 0.5: 0.42 / 0.49 / 0.62 at R²_x 0.3; 0.28 / 0.34 / 0.47 at 0.5);
T1 containment-df boundary rate 0.056 (MC SE 0.007), power 0.61 (0.56 given A) at ICC 0.7, 0.46 at ICC 0.5, 0.70 at
ICC 0.8.

## 7. Compliance and data provenance

- **Read** (public metadata copies kept in the Paper 0 project): `paper0/data/session_inventory_with_pennlinc_qc.csv`,
  `paper0/data/forensics/bundle_geometry.csv` (existence of an SCP_R profile only), `paper0/data/forensics/participants.tsv`,
  `paper0/data/pennlinc_qc/EF_desc-{T1,diffusion}_qc.csv` (pass/fail and mean FD), `paper0/data/protocolnames_wide.csv`
  (ProtocolName), `paper0/data/forensics/exclude_afq_qc.txt`, and the phenotype TSVs in `paper0/data/forensics/phenotype/`
  (PRIME, BDI, ARI, Tanner, CPT, RSAS, substance).
- **Not read:** any MEGRE image, `megre_coverage.tsv`, the sub-20139 local copy, the PNGs. No server access, no
  downloads, no installs.
- **Incidental exposure, declared (registration Q18 row P5).** In an early test run, before synthetic runs masked
  them, the console showed that 2 of the 4 exposure sessions (sub-20139 ses-1 and sub-21432 ses-1) are in the
  synthetic first-usable sample, from which their public cells could be inferred (TD and PS+ADHD).
- All imaging and DWI values in `synthetic/` are simulated (`data_origin = synthetic`; banner in every report).

## 8. Reproducibility

`RNGkind("Mersenne-Twister", "Inversion", "Rejection")`; base seed 20260924, one fixed offset per resampling
component. `SHA256SUMS` lists every file except itself; each run's `manifest.json` records the hashes of the code it
executed.
