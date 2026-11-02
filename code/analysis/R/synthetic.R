# LEAD-chi: the REAL Penn LEAD design (metadata only) plus SIMULATED outcomes, for testing the analysis code.
# Not part of the confirmatory run (run_all.R never sources this file).
# Inputs are metadata and phenotype copies only: the session inventory (MEGRE availability from the file list,
# T1w/DWI presence, ages), participants.tsv, PennLINC T1 and diffusion QC, the list of sessions that have a Paper 0
# SCP_R profile (existence only; no DWI value is read), and prime/bdi/ari TSVs. No MEGRE voxel data and no
# MEGRE-derived value is read. In particular megre_coverage.tsv is NOT used (E1 sec. 6.3).

# Stable local copies of PUBLIC metadata kept in the Paper 0 project (participants, PennLINC QC, ProtocolName table,
# round-1 AFQ-QC list, phenotype TSVs). No MEGRE-derived file is read.
P0_DIR   <- "/path/to/Penn_LEAD/paper0/data"
DATA_DIR <- P0_DIR

build_real_design <- function(data_dir = DATA_DIR, p0_dir = P0_DIR) {
  inv <- utils::read.csv(file.path(p0_dir, "session_inventory_with_pennlinc_qc.csv"), stringsAsFactors = FALSE)
  p <- utils::read.delim(file.path(data_dir, "forensics/participants.tsv"), stringsAsFactors = FALSE)
  t1 <- utils::read.csv(file.path(data_dir, "pennlinc_qc/EF_desc-T1_qc.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  t1$t1qc <- trimws(t1$qc_determination)
  dq <- utils::read.csv(file.path(data_dir, "pennlinc_qc/EF_desc-diffusion_qc.csv"), stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")
  dq$sub <- paste0("sub-", dq$participant_id); dq$ses <- paste0("ses-", dq$session_id); dq$scal <- trimws(dq$qc_determination_scalar_maps)
  geo <- utils::read.csv(file.path(p0_dir, "forensics/bundle_geometry.csv"), stringsAsFactors = FALSE)
  scp_prof <- paste(geo$subject[geo$bundle == "SCPR"], geo$session[geo$bundle == "SCPR"])   # existence of a profile only
  pn <- utils::read.csv(file.path(data_dir, "protocolnames_wide.csv"), stringsAsFactors = FALSE)   # MEGRE ProtocolName (metadata)
  afq <- grep("^sub-", readLines(file.path(data_dir, "forensics/exclude_afq_qc.txt")), value = TRUE)
  ph <- file.path(data_dir, "forensics/phenotype"); na <- c("n/a", "NA", "")
  rd <- function(f) utils::read.delim(file.path(ph, f), stringsAsFactors = FALSE, na.strings = na)
  pr <- rd("prime.tsv"); bd <- rd("bdi_1a_ms_child.tsv"); ar <- rd("ari.tsv"); tb <- rd("tanner_boy.tsv"); tg <- rd("tanner_girl.tsv")
  cpt <- rd("short_penn_continuous_performance_test.tsv"); rs <- rd("rsas.tsv"); su <- rd("substance.tsv")
  ses1_age <- setNames(inv$age[inv$ses == "ses-1"], inv$sub[inv$ses == "ses-1"])
  m <- inv[inv$MEGRE_full == 1, ]
  m <- merge(m, p[, c("participant_id", "dx_adhd", "race", "ethnicity", "dx_psychosis", "dx_mdd", "dx_moodnos")],
             by.x = "sub", by.y = "participant_id", all.x = TRUE)
  m <- merge(m, t1[, c("participant_id", "session_id", "t1qc")], by.x = c("sub", "ses"), by.y = c("participant_id", "session_id"), all.x = TRUE)
  m$sesn <- as.integer(sub("ses-", "", m$ses)); m <- m[order(m$sub, m$sesn), ]
  k <- paste(m$sub, m$ses); key <- function(t) paste(t$participant_id, t$session_id)
  at <- function(t, col) t[[col]][match(k, key(t))]
  items12 <- rowSums(sapply(paste0("prime_", 1:12), function(v) suppressWarnings(as.numeric(at(pr, v)))))
  stim <- as.integer(at(su, "substance_othr_010") %in% c(1, 2) & at(su, "substance_othr_011") %in% 1)
  out <- data.frame(participant_id = m$sub, session_id = m$ses, ses_order = stats::ave(m$sesn, m$sub, FUN = seq_along),
             study_group = m$study_group, dx_adhd = as.integer(m$dx_adhd), sex = m$sex, race = m$race, ethnicity = m$ethnicity,
             age = m$age, age_ses1 = unname(ses1_age[m$sub]),
             t1_ok = m$T1w == 1 & (is.na(m$t1qc) | m$t1qc != "fail"),
             t1_reason = ifelse(m$T1w != 1, "no T1w", ifelse(!is.na(m$t1qc) & m$t1qc == "fail", "T1 QC fail", "")),
             t1_artifact = as.integer(m$t1qc %in% "artifact"),
             dwi_ok = m$DWI == 1 & dq$scal[match(k, paste(dq$sub, dq$ses))] %in% "pass" & k %in% scp_prof,
             dwi_FD = dq$mean_FD[match(k, paste(dq$sub, dq$ses))],
             dx_psychosis = as.integer(m$dx_psychosis), dx_mdd = as.integer(m$dx_mdd), dx_moodnos = as.integer(m$dx_moodnos),
             epoch = ifelse(pn$anat_MEGRE[match(k, paste(pn$sub, pn$ses))] %in% "QSM_SWI_1.5mm", "early", "late"),
             stim_prescribed = stim, afq_round1_list = as.integer(m$sub %in% afq),
             PRIME_total = at(pr, "prime_total_score"), BDI_total = at(bd, "bdi_1a_ms_child_total_score"),
             ARI_total = at(ar, "ari_total_score"), PRIME_items12 = items12, CPT_fp = at(cpt, "cnb_cptnl_total_fp"),
             RSAS_total = at(rs, "rsas_total_score"), stringsAsFactors = FALSE)
  for (v in paste0("tanner_boy_", 1:5)) out[[v]] <- at(tb, v)
  for (v in paste0("tanner_girl_", c(1, 2, 3, 5, 6))) out[[v]] <- at(tg, v)
  out
}

# Scenario effects, in units of each ROI's residual SD (effect coding: coefficient = group difference).
SCENARIOS <- list(
  null = list(SN_PS = 0, SN_ADHD = 0, STR_PS = 0, STR_ADHD = 0, R2_PS = 0, rho_S3 = 0, dominant = 0),
  `delta-0.6` = list(SN_PS = -0.6, SN_ADHD = 0, STR_PS = -0.1, STR_ADHD = -0.4, R2_PS = -0.5,
                     rho_S3 = 0.4 * sqrt(0.70 * 0.55), dominant = 0),        # latent rho 0.4, R_chi 0.70, R_DWI 0.55
  null_dominant = list(SN_PS = 0, SN_ADHD = 0, STR_PS = 0, STR_ADHD = 0, R2_PS = 0, rho_S3 = 0, dominant = 2))

# Simulate one session-level ROI table. The design, QC draws, covariates and noise depend only on `seed`, so the
# scenarios differ ONLY in the group effects / coupling / dominant scans.
simulate_table <- function(des, scenario, seed, qc_loss = 0.05, inject = TRUE) {
  sc <- SCENARIOS[[scenario]]; set.seed(seed); n <- nrow(des); a <- des$age - 13
  subs <- unique(des$participant_id); ns <- length(subs); si <- match(des$participant_id, subs)
  PS <- ifelse(des$study_group == "PRO/CHR", 0.5, -0.5); AD <- ifelse(des$dx_adhd == 1, 0.5, -0.5)
  icc <- 0.7; C <- matrix(c(1, .5, .5, .5, 1, .3, .5, .3, 1), 3)          # SN-VTA, STR, R2* correlation
  zs <- matrix(stats::rnorm(ns * 3), ns) %*% chol(C); ze <- matrix(stats::rnorm(n * 3), n) %*% chol(C)
  z3 <- sqrt(icc) * zs[si, ] + sqrt(1 - icc) * ze                           # unit-variance correlated noise
  zind <- function() sqrt(icc) * stats::rnorm(ns)[si] + sqrt(1 - icc) * stats::rnorm(n)
  roi_vol <- stats::rnorm(n, 620 + 8 * a, 60); str_vol <- stats::rnorm(n, 14500 + 150 * a, 1400); Q_phase <- stats::rnorm(n, -2.8, 0.25)
  zr <- (roi_vol - mean(roi_vol)) / stats::sd(roi_vol); zq <- (Q_phase - mean(Q_phase)) / stats::sd(Q_phase)
  sSN <- 20; sST <- 12; sR2 <- 3
  chi_SNVTA <- 55 + 2.0 * a + 1.5 * zr - 1.0 * zq + sSN * (sc$SN_PS * PS + sc$SN_ADHD * AD) + sSN * z3[, 1]
  chi_STR <- 18 + 1.5 * a + 1.0 * (str_vol - mean(str_vol)) / stats::sd(str_vol) + sST * (sc$STR_PS * PS + sc$STR_ADHD * AD) + sST * z3[, 2]
  r2s_SNVTA <- 28 + 0.4 * a + sR2 * sc$R2_PS * PS + sR2 * z3[, 3]
  chi_GP <- 85 + 3 * a + 15 * zind(); chi_PUT <- 12 + 1.5 * a + 10 * zind(); chi_CAU <- 25 + 1.0 * a + 10 * zind()
  chi_THAL <- -5 + 6 * zind(); chi_RN <- 35 + 1.5 * a + 15 * zind()
  zdn <- zind(); chi_DN_R <- 40 + 1.5 * a + 15 * zdn
  chi_DGM <- 0.5 * (chi_CAU + chi_GP) + stats::rnorm(n, 0, 3); chi_PUT_R <- chi_PUT + stats::rnorm(n, 0, 4)
  dref <- 4 + 0.6 * a + stats::rnorm(n, 0, 4)
  ref2_type <- ifelse(des$age < 10 & stats::runif(n) < 0.5, "CSF_erode1", ifelse(stats::runif(n) < 0.02, "none", "CSF"))
  # S3 (DWI side simulated; only availability is real)
  rho <- sc$rho_S3
  twICVF <- 0.55 + 0.004 * a + 0.04 * (rho * zdn + sqrt(1 - rho^2) * stats::rnorm(n))
  twMCP <- 0.50 + 0.003 * a + 0.04 * stats::rnorm(n)
  nodes_ok <- stats::runif(n) < 0.97
  tab <- data.frame(data_origin = "synthetic", engine = "qsmxt", pipeline = "primary", roi_def = "SNVTA",
                    des[, c("participant_id", "session_id", "ses_order")],
                    code = sprintf("%06x", sample.int(16^6 - 1, ns))[si],
                    des[, c("study_group", "dx_adhd", "sex", "race", "ethnicity", "age", "age_ses1")],
                    usable = 0L, fail_reason = "", stringsAsFactors = FALSE)
  qc_fail <- stats::runif(n) < qc_loss                                     # label-free simulated MEGRE QC loss
  tab$usable <- as.integer(des$t1_ok & !qc_fail)
  tab$fail_reason <- ifelse(!des$t1_ok, des$t1_reason, ifelse(qc_fail, "simulated MEGRE QC grade 2", ""))
  tab <- cbind(tab, data.frame(chi_SNVTA, chi_SNVTA_ref2 = chi_SNVTA - dref, ref2_type, r2s_SNVTA, roi_vol, Q_phase, chi_STR, str_vol,
                               chi_GP, chi_PUT, chi_CAU, chi_THAL, chi_RN, des[, c("PRIME_total", "BDI_total", "ARI_total")],
                               s3_eligible = as.integer(des$dwi_ok & nodes_ok),
                               twICVF_pSCP_R = twICVF, b0loc_pSCP_R = stats::rnorm(n, 1, 0.08), tract_len = stats::rnorm(n, 53, 5),
                               nSL = round(exp(stats::rnorm(n, 6.3, 0.5))), FD = 0.3 + abs(stats::rnorm(n, 0.3, 0.2)),
                               NDC = pmin(0.95, stats::rnorm(n, 0.85, 0.03)), chi_DN_R, chi_DGM, chi_PUT_R,
                               DN_R_vol = stats::rnorm(n, 1600, 200), abs_pitch = abs(stats::rnorm(n, 0, 7)),
                               twICVF_MCP_R_matched = twMCP, r2s_pSCP_R = stats::rnorm(n, 22, 2)))
  # further optional columns (simulated imaging-side; real metadata-side columns come from `des`)
  tab$ref2_type <- ifelse(stats::runif(n) < 0.03, "none", tab$ref2_type)
  tab$chi_SNVTA_ref2[tab$ref2_type == "none"] <- NA
  tab$chi_SNVTA_split13 <- tab$chi_SNVTA + stats::rnorm(n, 0, 0.35 * sSN)
  tab$chi_SNVTA_split24 <- tab$chi_SNVTA + stats::rnorm(n, 0, 0.35 * sSN)
  tab$chi_SNVTA_impkeep <- tab$chi_SNVTA + abs(stats::rnorm(n, 0, 2))
  tab$chi_SNVTA_harm <- tab$chi_SNVTA + stats::rnorm(n, 0, 3)
  tab$ref_frac_wm <- stats::rnorm(n, 0.40, 0.03); tab$ref_frac_gm <- stats::rnorm(n, 0.45, 0.03); tab$ref_frac_csf <- 1 - tab$ref_frac_wm - tab$ref_frac_gm
  tab$fov_brain_frac <- pmin(1, stats::rnorm(n, 0.93, 0.03)); tab$motion_gsr <- exp(stats::rnorm(n, -3.5, 0.4))
  tab$dice_min <- pmin(0.97, stats::rnorm(n, 0.86, 0.03)); tab$icv_seg_ml <- stats::rnorm(n, 1450 + 10 * a, 120)
  tab$qc_grade <- ifelse(tab$usable == 1, as.integer(stats::runif(n) < 0.3), 2L); tab$dental_D <- as.integer(stats::runif(n) < 0.06)
  tab$f0_cluster <- sample(c("1", "2", "3"), n, replace = TRUE, prob = c(.3, .4, .3)); tab$sw_class <- "A"
  tab$chi_DN_L <- 40 + 1.5 * a + 15 * (0.6 * zdn + 0.8 * stats::rnorm(n))
  tab$chi_SNonly <- tab$chi_SNVTA + stats::rnorm(n, 4, 5); tab$chi_VTAPBP <- tab$chi_SNVTA + stats::rnorm(n, -6, 8)
  tab$twICVF_pSCP_R_peak <- twICVF + stats::rnorm(n, 0, 0.01); tab$t1t2_pSCP_R <- 1.2 + 0.3 * (0.3 * zdn * rho + stats::rnorm(n))
  for (v in c("t1_artifact", "dwi_FD", "dx_psychosis", "dx_mdd", "dx_moodnos", "epoch", "stim_prescribed", "afq_round1_list",
              "PRIME_items12", "CPT_fp", "RSAS_total", paste0("tanner_boy_", 1:5), paste0("tanner_girl_", c(1, 2, 3, 5, 6)))) tab[[v]] <- des[[v]]
  tab$fail_reason[tab$usable == 1] <- ""
  tab$fail_reason[tab$usable == 0 & tab$fail_reason == "simulated MEGRE QC grade 2"] <- "c_qc_grade"
  tab$fail_reason[tab$usable == 0 & tab$fail_reason %in% c("no T1w", "T1 QC fail")] <- "b_t1w"
  # frozen first-usable flag (earliest usable session)
  u <- tab[tab$usable == 1, ]; fu <- u$participant_id[!duplicated(u$participant_id)]
  fkey <- paste(u$participant_id, u$session_id)[!duplicated(u$participant_id)]
  tab$first_usable <- as.integer(paste(tab$participant_id, tab$session_id) %in% fkey)
  fi <- which(tab$first_usable == 1)
  if (inject) {   # exercise the ROI-level plausibility rules (non-primary ROIs only), same rows in every scenario
    j <- sample(fi, 2); tab$chi_DN_R[j[1]] <- 420; tab$r2s_SNVTA[j[2]] <- 150
  }
  if (sc$dominant > 0) {   # Paper 0 lesson (iv): scans that are physically plausible but dominate the residual SS
    j <- sample(setdiff(fi, if (inject) j else integer(0)), sc$dominant)
    tab$chi_SNVTA[j] <- tab$chi_SNVTA[j] + 6 * sSN; tab$chi_SNVTA_ref2[j] <- tab$chi_SNVTA_ref2[j] + 6 * sSN
  }
  tab
}

synthetic_gates <- function(tab, scenario) {
  f <- tab[tab$first_usable == 1, ]
  list(data_origin = "synthetic", scenario = scenario,
       note = "SYNTHETIC: gate outcomes are assumed values; pooled constants are computed label-free from the synthetic first-usable sessions, as gates_preunlock.py would",
       A2_approved = TRUE, PC1 = TRUE, PC2w = "PASS", PC3 = TRUE, PC4 = TRUE, r_point_primary = 0.80, rL_primary = 0.70,
       rL_cross = 0.55, r_cross = 0.68, icc_cross_session = 0.66,
       PC4_DN = TRUE, r_point_DN = 0.78, rL_DN = 0.72, PC5 = TRUE, G2 = TRUE, G2_DN = TRUE, G3 = TRUE, G4 = TRUE, G5 = TRUE,
       n_min = 10, constants_roi = "SNVTA",
       pooled = list(mean_SNVTA = mean(f$chi_SNVTA), sd_SNVTA = stats::sd(f$chi_SNVTA),
                     mean_STR = mean(f$chi_STR, na.rm = TRUE), sd_STR = stats::sd(f$chi_STR, na.rm = TRUE),
                     sd_total_SNVTA = stats::sd(f$chi_SNVTA)))
}

# Long multiverse table (schema/multiverse_schema.csv): 9 background-removal x inversion nodes per session, summarised on
# the common support. Simulated as the primary value plus a node-specific scale/offset and noise (SYNTHETIC).
simulate_multiverse <- function(tab, seed) {
  set.seed(seed + 7)
  nodes <- expand.grid(bf = c("vsharp", "pdf", "lbv"), inv = c("medi", "ilsqr", "rts"), stringsAsFactors = FALSE)
  other <- tab$chi_SNVTA + stats::rnorm(nrow(tab), 5, 6)
  dref <- tab$chi_SNVTA - tab$chi_SNVTA_ref2
  do.call(rbind, lapply(seq_len(nrow(nodes)), function(i) {
    sc <- if (nodes$inv[i] == "medi") 1 else if (nodes$inv[i] == "ilsqr") 0.85 else 0.92
    off <- c(vsharp = 0, pdf = 3, lbv = -2)[nodes$bf[i]]
    e <- if (nodes$bf[i] == "vsharp" && nodes$inv[i] == "medi") stats::rnorm(nrow(tab), 0, 1) else stats::rnorm(nrow(tab), 0, 5)
    wb <- sc * tab$chi_SNVTA + off + e
    fail <- stats::runif(nrow(tab)) < 0.02
    wb[fail] <- NA
    data.frame(participant_id = tab$participant_id, session_id = tab$session_id, bf = nodes$bf[i], inv = nodes$inv[i],
               chi_wb_mean = wb, chi_wb_median = wb + stats::rnorm(nrow(tab), 0, 3), chi_other_wb_mean = sc * other + off + e,
               chi_csf_mean = wb - dref, chi_csf_median = wb - dref + stats::rnorm(nrow(tab), 0, 3),
               chi_other_csf_mean = sc * other + off + e - dref, stringsAsFactors = FALSE)
  }))
}
