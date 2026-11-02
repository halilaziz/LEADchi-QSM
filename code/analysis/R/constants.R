# LEAD-chi confirmatory analysis code -- registered constants (frozen with the registration; hashed).
# Source of every value: the registered text LEADchi_prereg_v1.0_FINAL.md (which governs; Annex H maps each rule to
# the function implementing it). Protocol v1.1, the v1.2 delta (E2), E1 and E3 are superseded where they differ.

REG <- list(
  code_version     = "LEAD-chi analysis v1.0-final (2026-09-26)",
  alpha            = 0.05,                 # D1: two-sided .05 (Q22)
  # ---- H1 decision rule (sec. 4.1) ----
  deltaB_base      = 0.6,                  # delta_B = 0.6 * sqrt(r_L)
  rL_floor         = 0.25,                 # B is evaluable only if PC4 passes and r_L >= 0.25 (else "low-precision measure")
  deltaB_coreport  = c(0.6, 0.5),          # B also reported at the unadjusted -0.6 and at -0.5 (never decisive)
  huber_k          = 1.345,                # sigma_rob: MASS::rlm Huber psi, k = 1.345, proposal-2 scale
  huber_maxit      = 200,                  #   non-convergence -> sigma_rob = 1.4826 * MAD(OLS residuals), flagged
  # ---- qualifiers K1-K4 (sec. 4.2) ----
  tost_bound       = 0.5,                  # K1b ADHD TOST at +/- 0.5 * sigma_B
  comorbid_ratio   = 0.5,                  # K1c: est(PS-only - TD) > 0.5 * est(PS+ADHD - ADHD-only), comorbid est < 0
  drift_rel        = 0.25,                 # K3 trigger: |b(dref)| > max(0.25|b_PS|, 0.10 sigma_resid) ...
  drift_abs        = 0.10,
  drift_p          = 0.10,                 #   ... AND two-sided HC3 p(dref) < .10 (inferential trigger)
  dom_share        = 0.20,                 # K4 influence: top-1 share of residual SS
  dom_t            = 5,                    #            max |externally studentised residual|
  robust_share     = 0.80,                 # K4 "robust": >= 80 % of evaluable specifications with delta < 0
  spec_min_evaluable = 400,                #   ... of 486; fewer evaluable specifications -> "robust" not assessable
  spec_min_n_frac  = 0.80,                 #   evaluable if n >= 30 and >= 0.80 x the n of its sample level
  # ---- T1 (sec. 4.4) ----
  t1_margin        = 0.4,                  # Delta, total-SD units; decision = containment-df t test (PB co-reported)
  # ---- S2 ----
  s2_age_min       = 11,
  adult_ppb        = -14.24,               # Vano 2025 Table 3 adjusted SN-VTA difference (descriptive co-report)
  # ---- resampling ----
  n_perm           = 10000,                # Freedman-Lane (M1, S3), identity counted: p = (#>= + 1)/(N + 1)
  n_perm_spec      = 10000,                # joint specification-curve permutation (decision-relevant for "robust")
  n_boot_delta     = 5000,                 # cell-stratified bootstrap CI of delta
  n_boot_s3b       = 5000,                 # S3b NC1/NC2 case bootstrap (cell-stratified)
  n_pb             = 2000,                 # T1 parametric bootstrap (co-report only)
  n_wild           = 9999,                 # wild cluster bootstrap, Rademacher (co-report only)
  mc_band          = c(0.03, 0.07),        # Monte Carlo precision rule (Q22)
  mc_recompute     = 10000,
  seed             = 20260924,
  bf_priors        = list(attenuated = c(mean = -0.35, sd = 0.25), sceptical = c(mean = 0, sd = 0.5)),
  # ---- physical plausibility (sec. 12.4); composites = component envelopes ----
  env_chi = list(chi_GP = c(-50, 300), chi_PUT = c(-80, 200), chi_CAU = c(-80, 200), chi_SNVTA = c(-50, 300),
                 chi_RN = c(-50, 300), chi_DN_R = c(-50, 300), chi_DN_L = c(-50, 300), chi_THAL = c(-100, 100),
                 chi_STR = c(-80, 200), chi_DGM = c(-80, 300), chi_PUT_R = c(-80, 200), chi_DNIP_R = c(-50, 300)),
  env_r2s          = c(5, 100),
  ref2_max_abs     = 50,                   # |mean CSF - whole-brain mean| > 50 ppb -> CSF reference invalid (missing)
  robust_ratio     = 1.5,                  # SD / (1.4826 MAD) above which robust scales are co-reported
  robz             = 5,
  s3_r2s_drop      = 0.5,
  # ---- decisions ----
  D12_S3_in_holm   = TRUE,                 # D12, confirmed 2026-09-26 with the registered-test power (Q21)
  exposure_sessions = c("sub-20139_ses-1", "sub-21432_ses-1", "sub-20976_ses-2", "sub-20188_ses-3"),
  td_coding_gap    = c("sub-20139", "sub-20252", "sub-20303", "sub-20953"),   # Q13 sensitivity exclusion
  t1qc_missing_sub = "sub-21568",
  reml_theta_grid  = seq(0, 3, by = 0.01),
  # ---- environment pin (Q23): the confirmatory run on real data stops if any version differs ----
  pinned = list(R = "4.4.1", sandwich = "3.1.1", lmtest = "0.9.40", lme4 = "1.1.38", nlme = "3.1.164",
                MASS = "7.3.61", boot = "1.3.30", jsonlite = "1.8.9", digest = "0.6.37", Matrix = "1.7.0"),
  # gate-JSON fields that must be present (any data origin); missing or invalid -> hard stop
  gate_required    = c("data_origin", "PC1", "PC2w", "PC3", "PC4", "PC4_DN", "PC5", "G2", "G2_DN", "G3", "G4", "G5",
                       "A2_approved", "rL_primary", "r_point_primary", "constants_roi", "pooled")
)

CELLS <- c("TD", "ADHD-only", "PS-only", "PS+ADHD")
