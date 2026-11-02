# LEAD-chi analysis: key secondary tests S1-S3 and the Holm family (protocol sec. 8.1-8.3, 9; v1.2 delta sec. 4-5).

# ---------------- S1: topographic contrast L and the intersection-union "dissociation" label ----------------
s1_rhs <- function(d) paste0("PS_e * ADHD_e + ns(age, 3) + sex_e", race_term(d), " + scale(roi_vol) + scale(str_vol) + scale(Q_phase)")

s1_fit <- function(d, cst) {
  d$D <- (d$chi_SNVTA - cst$mean_SNVTA) / cst$sd_SNVTA - (d$chi_STR - cst$mean_STR) / cst$sd_STR
  m <- fit_ols_hc3(stats::as.formula(paste("D ~", s1_rhs(d))), d)
  L <- lin_test(m, c(PS_e = 1, ADHD_e = -1), "L = beta_PS(D) - beta_ADHD(D)")
  list(m = m, L = L, ps = term_test(m, "PS_e"), adhd = term_test(m, "ADHD_e"), data = d)
}

run_s1 <- function(S, gates) {
  d <- set_race3(S$first[!is.na(S$first$chi_STR) & !is.na(S$first$str_vol), ])
  # pooled label-free constants are REQUIRED from the gate JSON (validated in run_all.R; no in-table fallback)
  cst <- gates$pooled[c("mean_SNVTA", "sd_SNVTA", "mean_STR", "sd_STR")]
  pooled_src <- "gates JSON (label-free, computed before unlock by gates_preunlock.py)"
  r <- s1_fit(d, cst)
  # per-ROI effects (reported)
  per_roi <- lapply(c(SNVTA = "chi_SNVTA", STR = "chi_STR"), function(v) {
    fo <- stats::as.formula(paste(v, "~", s1_rhs(d))); mm <- fit_ols_hc3(fo, d)
    list(beta_PS = term_test(mm, "PS_e")[c("est", "se", "p_two")], beta_ADHD = term_test(mm, "ADHD_e")[c("est", "se", "p_two")])
  })
  # robust-scale co-report (v1.2 sec. 5.5.4)
  ratio <- c(SNVTA = stats::sd(d$chi_SNVTA) / stats::mad(d$chi_SNVTA), STR = stats::sd(d$chi_STR) / stats::mad(d$chi_STR))
  robust <- NULL
  if (any(ratio > REG$robust_ratio)) {
    cr <- list(mean_SNVTA = stats::median(d$chi_SNVTA), sd_SNVTA = stats::mad(d$chi_SNVTA),
               mean_STR = stats::median(d$chi_STR), sd_STR = stats::mad(d$chi_STR))
    rr <- s1_fit(d, cr); robust <- rr$L[c("est", "se", "p_two")]
  }
  # stacked two-ROI GLS sensitivity with ROI-specific variances (nlme::gls, varIdent) and within-subject correlation
  gls_res <- try_fit({
    lg <- rbind(data.frame(r$data, roi = "SNVTA", z = (d$chi_SNVTA - cst$mean_SNVTA) / cst$sd_SNVTA, rk = 1L),
                data.frame(r$data, roi = "STR", z = (d$chi_STR - cst$mean_STR) / cst$sd_STR, rk = 2L))
    lg$roi <- factor(lg$roi, levels = c("SNVTA", "STR")); lg <- set_race3(lg)
    fo <- stats::as.formula(paste0("z ~ 0 + roi + roi:(PS_e * ADHD_e + ns(age, 3) + sex_e", race_term(lg),
                                   " + scale(roi_vol) + scale(str_vol) + scale(Q_phase))"))
    g <- nlme::gls(fo, data = lg, weights = nlme::varIdent(form = ~ 1 | roi), correlation = nlme::corSymm(form = ~ rk | sub),
                   method = "REML")
    b <- stats::coef(g); V <- stats::vcov(g)
    cv <- setNames(rep(0, length(b)), names(b))
    cv[c("roiSNVTA:PS_e", "roiSTR:PS_e", "roiSNVTA:ADHD_e", "roiSTR:ADHD_e")] <- c(1, -1, -1, 1)
    est <- sum(cv * b); se <- sqrt(drop(t(cv) %*% V %*% cv)); dfg <- nrow(lg) - length(b)
    list(L = est, se = se, p_two = 2 * stats::pt(-abs(est / se), dfg))
  }, "S1 stacked gls")
  list(n = r$m$n, cells = cell_counts(d), pooled_constants = cst, pooled_source = pooled_src,
       L = r$L, beta_PS_D = r$ps, beta_ADHD_D = r$adhd, coefficients = coef_table(r$m), per_roi = per_roi,
       sd_over_robust_sd = as.list(ratio), robust_scaled_L = robust, stacked_gls = gls_res,
       model = r$m, data = r$data)
}

s1_labels <- function(s1, holm_reject) {
  sigL <- isTRUE(holm_reject) && s1$L$est < 0
  cps <- s1$beta_PS_D$p_lower < REG$alpha; cad <- s1$beta_ADHD_D$p_upper < REG$alpha
  w <- if (!sigL) "L not significant under Holm (no S1 claim)" else if (cps && cad) "dissociation" else
    if (cps) "the PS effect is relatively more nigral" else if (cad) "the ADHD effect is relatively more striatal" else
      "L significant; neither component significant (worded as the contrast only)"
  list(L_significant_holm = sigL, beta_PS_D_lt0 = cps, beta_ADHD_D_gt0 = cad, wording = w)
}

# ---------------- S2: dimensional PRIME (age >= 11; BDI, ARI and ADHD adjusted) ----------------
s2_rhs <- function(d, extra = "") paste0("scale(PRIME_total) + scale(BDI_total) + scale(ARI_total) + ADHD_e + ns(age, 3) + sex_e",
                                           race_term(d), " + scale(roi_vol) + scale(Q_phase)", extra)
run_s2 <- function(S) {
  d <- set_race3(S$s2)
  m <- fit_ols_hc3(stats::as.formula(paste("chi_SNVTA ~", s2_rhs(d))), d); tp <- term_test(m, "scale(PRIME_total)")
  mp <- fit_ols_hc3(stats::as.formula(paste("chi_SNVTA ~", s2_rhs(d, " + PS_e"))), d); tpp <- term_test(mp, "scale(PRIME_total)")
  dps <- set_race3(d[d$PS_e > 0, ])
  mw <- try_fit(fit_ols_hc3(stats::as.formula(paste("chi_SNVTA ~", s2_rhs(dps))), dps), "S2 within PS")
  d$PRIME_rin <- stats::qnorm((rank(d$PRIME_total) - 0.5) / nrow(d))
  mr <- fit_ols_hc3(stats::as.formula(paste("chi_SNVTA ~", sub("scale(PRIME_total)", "PRIME_rin", s2_rhs(d), fixed = TRUE))), d)
  beyond <- isTRUE((tpp$est < 0) && (tpp$p_two < REG$alpha))
  sens <- list()
  # PRIME items 1-12 without the duration items (sum of prime_1 ... prime_12)
  if ("PRIME_items12" %in% names(d) && sum(!is.na(d$PRIME_items12)) > 30) {
    di <- set_race3(d[!is.na(d$PRIME_items12), ])
    sens$PRIME_items_1_12 <- term_test(fit_ols_hc3(stats::as.formula(paste("chi_SNVTA ~",
      sub("scale(PRIME_total)", "scale(PRIME_items12)", s2_rhs(di), fixed = TRUE))), di), "scale(PRIME_items12)")[c("est", "se", "p_two")]
  }
  # CPT commission errors (cnb_cptnl_total_fp) in place of dx_adhd
  if ("CPT_fp" %in% names(d) && sum(!is.na(d$CPT_fp)) > 30) {
    dc <- set_race3(d[!is.na(d$CPT_fp), ])
    sens$CPT_in_place_of_ADHD <- term_test(fit_ols_hc3(stats::as.formula(paste("chi_SNVTA ~",
      sub("ADHD_e", "scale(CPT_fp)", s2_rhs(dc), fixed = TRUE))), dc), "scale(PRIME_total)")[c("est", "se", "p_two")]
  }
  list(n = m$n, cells = cell_counts(d), beta_PRIME = tp, coefficients = coef_table(m),
       plus_PS_e = tpp[c("est", "se", "p_two")], beyond_diagnosis_wording_allowed = beyond,
       within_PS = if (is.null(mw)) NULL else term_test(mw, "scale(PRIME_total)")[c("est", "se", "p_two")],
       rank_inverse_normal = term_test(mr, "PRIME_rin")[c("est", "se", "p_two")], sensitivities = sens,
       registered_wording = "SN-VTA chi covaries with PRIME, a symptom dimension that partly indexes PS status.",
       model = m, data = d)
}

# S2 between-person all-session sensitivity (run after T1): every usable session with age >= 11 and PRIME/BDI/ARI;
# PRIME split into the subject mean and the within-subject deviation; REML maximum; reported, never decisive.
run_s2_allsession <- function(S) {
  d <- S$usable[S$usable$age >= REG$s2_age_min & !is.na(S$usable$PRIME_total) & !is.na(S$usable$BDI_total) &
                  !is.na(S$usable$ARI_total), ]
  if (nrow(d) < 40 || length(unique(d$sub)) < 30) return(list(note = "too few sessions"))
  d <- set_race3(d)
  d$PRIME_bar <- stats::ave(d$PRIME_total, d$sub); d$PRIME_w <- d$PRIME_total - d$PRIME_bar
  d$age_bar <- stats::ave(d$age, d$sub); d$age_w <- d$age - d$age_bar
  sdp <- stats::sd(d$PRIME_total); d$PRIME_bar_z <- (d$PRIME_bar - mean(d$PRIME_bar)) / sdp; d$PRIME_w_z <- d$PRIME_w / sdp
  fx <- stats::as.formula(paste0("chi_SNVTA ~ PRIME_bar_z + PRIME_w_z + scale(BDI_total) + scale(ARI_total) + ADHD_e + ns(age_bar, 3) + age_w + sex_e",
                                 race_term(d), " + scale(roi_vol) + scale(Q_phase)"))
  rm <- try_fit(reml_max(fx, d, profile = FALSE), "S2 all-session")
  if (is.null(rm)) return(list(note = "fit failed"))
  list(n_sessions = nrow(d), n_subjects = length(unique(d$sub)), reml_best = rm$best,
       between = list(est = unname(rm$b["PRIME_bar_z"]), se = unname(rm$se["PRIME_bar_z"])),
       within = list(est = unname(rm$b["PRIME_w_z"]), se = unname(rm$se["PRIME_w_z"])))
}

# ---------------- S3: right dentate chi <-> pSCP_R tissue-weighted ICVF ----------------
s3_rhs <- function(d, pred = "scale(chi_DN_R)", b0 = TRUE, extra = "")
  paste0(pred, " + scale(chi_DGM)", if (b0) " + scale(b0loc_pSCP_R)" else "", " + scale(abs_pitch) + scale(tract_len) + ns(age, 3) + sex_e",
         race_term(d), " + PS_e * ADHD_e + scale(FD) + scale(NDC) + scale(log(nSL)) + scale(DN_R_vol) + scale(Q_phase)", extra)
s3_fit <- function(d, y = "twICVF_pSCP_R", pred = "scale(chi_DN_R)", b0 = TRUE, extra = "")
  fit_ols_hc3(stats::as.formula(paste0("scale(", y, ") ~ ", s3_rhs(d, pred, b0, extra))), d)

run_s3 <- function(S, gates, nperm, nboot, seed) {
  conf0 <- isTRUE(gates$PC4_DN) && isTRUE(gates$G2_DN) && isTRUE(gates$PC5)
  if (is.null(S$s3) || nrow(S$s3) < 10) {   # post-unlock coding fix (Annex C.4): S3 exploratory and no DWI S3 table merged
    na_t <- list(est = NA_real_, se = NA_real_, t = NA_real_, p_two = NA_real_, n = 0L)
    nc0 <- list(obs = NA_real_, lower95_one_sided = NA_real_, n_valid = 0L, p_boot_le0 = NA_real_)
    return(list(n = 0L, cells = cell_counts(S$s3), n_dropped_incomplete = length(S$s3_dropped_incomplete), beta_DN = na_t,
                coefficients = NULL, freedman_lane = list(p = NA_real_, n_perm = 0L), negative_controls = list(NC1_DN_minus_PUT = nc0, NC2_pSCP_minus_MCPmatched = nc0),
                sensitivities = list(), confirmatory = conf0, gates_used = list(PC4_DN = gates$PC4_DN, G2_DN = gates$G2_DN, PC5 = gates$PC5),
                model = NULL, data = S$s3, not_run = "S3 not fitted: the S3 sample is empty (no frozen DWI S3 table was merged; S3 is exploratory after PC5 fail and n_min = none). Deferred to the exploratory amendment-2 modules."))
  }
  d <- set_race3(S$s3); term <- "scale(chi_DN_R)"
  m <- s3_fit(d); tt <- term_test(m, term)
  fl <- fl_perm(m$fit, term, d$ADHD_e, nperm, seed + 11)
  # S3b negative controls (NC1 putamen predictor; NC2 matched-MCP outcome): one-sided 95% case bootstrap (cell-stratified)
  dn1 <- d[!is.na(d$chi_PUT_R) & !is.na(d$twICVF_MCP_R_matched), ]
  stat <- function(dd, i) {
    x <- set_race3(dd[i, ])
    b_dn <- stats::coef(stats::lm(stats::as.formula(paste0("scale(twICVF_pSCP_R) ~ ", s3_rhs(x))), data = x))[term]
    b_put <- stats::coef(stats::lm(stats::as.formula(paste0("scale(twICVF_pSCP_R) ~ ", s3_rhs(x, "scale(chi_PUT_R)"))), data = x))["scale(chi_PUT_R)"]
    b_mcp <- stats::coef(stats::lm(stats::as.formula(paste0("scale(twICVF_MCP_R_matched) ~ ", s3_rhs(x))), data = x))[term]
    c(b_dn - b_put, b_dn - b_mcp)
  }
  bnc <- function(R, sd) { set.seed(sd); boot::boot(dn1, stat, R = R, strata = factor(dn1$cell)) }
  bt <- bnc(nboot, seed + 12)
  obs <- bt$t0
  mk <- function(bt) lapply(1:2, function(k) { v <- bt$t[, k]; v <- v[is.finite(v)]
    p <- (sum(v <= 0) + 1) / (length(v) + 1)
    list(obs = unname(bt$t0[k]), lower95_one_sided = unname(stats::quantile(v, 0.05)), n_valid = length(v), p_boot_le0 = p) })
  nc <- mk(bt)
  if (any(vapply(nc, function(z) z$p_boot_le0 >= REG$mc_band[1] && z$p_boot_le0 <= REG$mc_band[2], logical(1))) && nboot < REG$mc_recompute) {
    log_msg("S3b bootstrap p in the Monte Carlo band: recomputing with ", REG$mc_recompute, " draws")
    bt <- bnc(REG$mc_recompute, seed + 13); nc <- mk(bt)
  }
  names(nc) <- c("NC1_DN_minus_PUT", "NC2_pSCP_minus_MCPmatched")
  # sensitivities
  sens <- list()
  sens$neither_b0_nor_R2s <- term_test(s3_fit(d, b0 = FALSE), term)[c("est", "se", "p_two")]
  if ("r2s_pSCP_R" %in% names(d) && sum(!is.na(d$r2s_pSCP_R)) > 30) {
    dr <- set_race3(d[!is.na(d$r2s_pSCP_R), ])
    b_base <- term_test(s3_fit(dr), term)$est
    tr <- term_test(s3_fit(dr, extra = " + scale(r2s_pSCP_R)"), term)
    sens$plus_nodewise_R2s <- c(tr[c("est", "se", "p_two")], list(n = nrow(dr), beta_without = b_base,
      reduction = 1 - tr$est / b_base, label = if (b_base > 0 && (1 - tr$est / b_base) > REG$s3_r2s_drop)
        "not separable from local relaxation (artefact or mediation)" else "none"))
  }
  if ("twICVF_VTAuncross_R" %in% names(d) && sum(!is.na(d$twICVF_VTAuncross_R)) > 30) {
    dv <- set_race3(d[!is.na(d$twICVF_VTAuncross_R), ])
    sens$VTAuncross_R_tract <- term_test(s3_fit(dv, y = "twICVF_VTAuncross_R"), term)[c("est", "se", "p_two")]
  }
  if ("chi_DNIP_R" %in% names(d) && sum(!is.na(d$chi_DNIP_R)) > 30) {
    dc <- set_race3(d[!is.na(d$chi_DNIP_R), ])
    sens$dentate_interposed_composite <- term_test(s3_fit(dc, pred = "scale(chi_DNIP_R)"), "scale(chi_DNIP_R)")[c("est", "se", "p_two")]
  }
  if ("chi_DN_L" %in% names(d) && sum(!is.na(d$chi_DN_L)) > 30) {   # laterality control: left dentate as predictor
    dl <- set_race3(d[!is.na(d$chi_DN_L), ])
    sens$left_dentate_predictor <- term_test(s3_fit(dl, pred = "scale(chi_DN_L)"), "scale(chi_DN_L)")[c("est", "se", "p_two")]
  }
  if ("twICVF_pSCP_R_peak" %in% names(d) && sum(!is.na(d$twICVF_pSCP_R_peak)) > 30) {   # peak alignment (shift <= 5 nodes)
    dp <- set_race3(d[!is.na(d$twICVF_pSCP_R_peak), ])
    sens$peak_aligned <- term_test(s3_fit(dp, y = "twICVF_pSCP_R_peak"), term)[c("est", "se", "p_two")]
  }
  if ("t1t2_pSCP_R" %in% names(d) && sum(!is.na(d$t1t2_pSCP_R)) > 30) {   # T1w/T2w ratio in place of tw-ICVF
    dt <- set_race3(d[!is.na(d$t1t2_pSCP_R), ])
    sens$T1wT2w_outcome <- term_test(s3_fit(dt, y = "t1t2_pSCP_R"), term)[c("est", "se", "p_two")]
  }
  if ("afq_round1_list" %in% names(d)) {   # conservative sample: drop the 20-subject round-1 AFQ-QC list
    dk <- set_race3(d[d$afq_round1_list %in% 0, ])
    if (nrow(dk) > 30) sens$conservative_sample <- c(term_test(s3_fit(dk), term)[c("est", "se", "p_two")], list(n = nrow(dk)))
  }
  # exposure sensitivity (Q18(4)): the four listed session rows dropped, no substitution
  de <- set_race3(S$s3_exposure)
  sens$exposure_sensitivity <- c(term_test(s3_fit(de), term)[c("est", "se", "p_two")],
                                 list(n = nrow(de), n_dropped = nrow(d) - nrow(de)))
  conf <- isTRUE(gates$PC4_DN) && isTRUE(gates$G2_DN) && isTRUE(gates$PC5)
  list(n = m$n, cells = cell_counts(d), n_dropped_incomplete = length(S$s3_dropped_incomplete), beta_DN = tt,
       coefficients = coef_table(m), freedman_lane = fl, negative_controls = nc, sensitivities = sens,
       confirmatory = conf, gates_used = list(PC4_DN = gates$PC4_DN, G2_DN = gates$G2_DN, PC5 = gates$PC5),
       model = m, data = d)
}

s3_labels <- function(s3, holm_reject) {
  a <- isTRUE(holm_reject) && s3$beta_DN$est > 0
  spec <- a && s3$negative_controls$NC1_DN_minus_PUT$lower95_one_sided > 0 && s3$negative_controls$NC2_pSCP_minus_MCPmatched$lower95_one_sided > 0
  list(S3a_supported = a, S3b_specific = spec,
       wording = if (!a) "no S3a claim" else if (spec) "coupling specific to the dentate->SCP pair" else "coupled (specificity not established)")
}

# ---------------- Holm family ----------------
holm_family <- function(s1, s2, s3) {
  inS3 <- isTRUE(s3$confirmatory) && isTRUE(REG$D12_S3_in_holm)
  p <- c(S1 = s1$L$p_two, S2 = s2$beta_PRIME$p_two); if (inS3) p <- c(p, S3a = s3$beta_DN$p_two)
  h <- holm(p)
  dir_ok <- c(S1 = s1$L$est < 0, S2 = s2$beta_PRIME$est < 0, S3a = s3$beta_DN$est > 0)[names(p)]
  list(members = names(p), S3_in_family = inS3,
       S3_status = if (inS3) "confirmatory (Holm)" else if (!isTRUE(s3$confirmatory)) "exploratory (a dentate gate or PC5 failed)" else
         "Tier-3 estimation (D12 = demote; not used)",
       p = as.list(h$p), p_holm = as.list(h$p_holm), reject = as.list(h$reject),
       supported = as.list(h$reject & dir_ok), opposite_direction_significant = as.list(h$reject & !dir_ok))
}

# Influence refits for S1-S3 (v1.2 sec. 5.5.2): does the Holm decision change when dominant subjects are dropped?
secondary_influence <- function(s1, s2, s3, H, gates) {
  refit_p <- function(which) {
    if (which == "S1") { m <- s1$model; d <- s1$data; dom <- dominance(m$fit, d$sess_key); if (!dom$triggered) return(list(dom = dom))
      dd <- set_race3(d[!d$sess_key %in% dom$drop_ids, ]); r <- s1_fit(dd, s1$pooled_constants); return(list(dom = dom, p = r$L$p_two, est = r$L$est)) }
    if (which == "S2") { m <- s2$model; d <- s2$data; dom <- dominance(m$fit, d$sess_key); if (!dom$triggered) return(list(dom = dom))
      dd <- set_race3(d[!d$sess_key %in% dom$drop_ids, ])
      t <- term_test(fit_ols_hc3(stats::as.formula(paste("chi_SNVTA ~", s2_rhs(dd))), dd), "scale(PRIME_total)"); return(list(dom = dom, p = t$p_two, est = t$est)) }
    if (is.null(s3$model)) return(list(dom = list(triggered = FALSE, note = "S3 not fitted")))   # C.4 fix: empty S3
    m <- s3$model; d <- s3$data; dom <- dominance(m$fit, d$sess_key); if (!dom$triggered) return(list(dom = dom))
    dd <- set_race3(d[!d$sess_key %in% dom$drop_ids, ]); t <- term_test(s3_fit(dd), "scale(chi_DN_R)"); list(dom = dom, p = t$p_two, est = t$est)
  }
  out <- lapply(c(S1 = "S1", S2 = "S2", S3a = "S3a"), refit_p)
  for (k in names(out)) {
    out[[k]]$influence_dependent <- FALSE
    if (!is.null(out[[k]]$p) && k %in% H$members) {
      p2 <- unlist(H$p); p2[k] <- out[[k]]$p; r2 <- holm(p2)$reject[k]
      out[[k]]$holm_reject_after_drop <- unname(r2); out[[k]]$influence_dependent <- unname(r2) != unlist(H$reject)[k]
    }
  }
  out
}
