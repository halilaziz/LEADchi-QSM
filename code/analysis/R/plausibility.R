# LEAD-chi analysis: QSM physical-plausibility re-verification and influence / dominance rules (sec. 12.4, Q15).
# The per-session plausibility flags are applied label-free BEFORE unlock by qsm_plausibility.py (they feed `usable`).
# Here the same ROI-level rules are re-verified on the merged table, so that nothing implausible reaches a model:
#   - primary-ROI (chi_SNVTA column) value outside its envelope, or failing the voxel rules, in a usable row
#     -> hard stop (the frozen first-usable list contradicts criterion (c'); resolved only as a documented deviation);
#   - any other ROI outside its envelope, or failing the voxel rules -> set to missing for that session;
#   - ROI R2* outside [5, 100] s^-1 -> missing;
#   - |chi_SNVTA - chi_SNVTA_ref2| > 50 ppb -> the CSF reference is invalid: ref2 set to missing, type 'none'.

apply_plausibility <- function(d) {
  logt <- data.frame(sess_key = character(0), column = character(0), value = numeric(0), rule = character(0),
                     action = character(0), stringsAsFactors = FALSE)
  add <- function(i, col, rule, action) {
    if (length(i)) logt <<- rbind(logt, data.frame(sess_key = d$sess_key[i], column = col, value = d[[col]][i],
                                                   rule = rule, action = action, stringsAsFactors = FALSE))
  }
  for (col in intersect(names(REG$env_chi), names(d))) {
    e <- REG$env_chi[[col]]
    i <- which(!is.na(d[[col]]) & (d[[col]] < e[1] | d[[col]] > e[2]))
    if (col == "chi_SNVTA") {
      bad <- i[d$usable[i] == 1]
      if (length(bad)) stop("SN-VTA ROI mean outside the registered envelope in usable session(s): ",
                            paste(d$sess_key[bad], collapse = ", "),
                            ". The frozen first-usable list contradicts criterion (c'); resolve as a documented deviation.")
      add(i, col, sprintf("envelope [%g, %g] ppb", e[1], e[2]), "none (session already not usable)")
    } else {
      add(i, col, sprintf("envelope [%g, %g] ppb", e[1], e[2]), "set to missing")
      d[[col]][i] <- NA
    }
  }
  for (roi in c("SNVTA", "DN_R")) {
    fi <- paste0("fimp_", roi); im <- paste0("inmask_", roi); col <- paste0("chi_", roi)
    thr <- if (roi == "DN_R") 0.90 else 0.95
    if (fi %in% names(d) && col %in% names(d)) {
      i <- which(!is.na(d[[fi]]) & d[[fi]] > 0.05 & !is.na(d[[col]]))
      if (roi == "SNVTA" && any(d$usable[i] == 1)) stop("SN-VTA voxel-rule failure in usable session(s): ",
                                                       paste(d$sess_key[i[d$usable[i] == 1]], collapse = ", "))
      if (roi != "SNVTA") { add(i, col, "implausible-voxel share > 0.05", "set to missing"); d[[col]][i] <- NA }
    }
    if (im %in% names(d) && col %in% names(d)) {
      i <- which(!is.na(d[[im]]) & d[[im]] < thr & !is.na(d[[col]]))
      if (roi == "SNVTA" && any(d$usable[i] == 1)) stop("SN-VTA in-support failure in usable session(s): ",
                                                       paste(d$sess_key[i[d$usable[i] == 1]], collapse = ", "))
      if (roi != "SNVTA") { add(i, col, sprintf("in-support share < %.2f", thr), "set to missing"); d[[col]][i] <- NA }
    }
  }
  for (col in intersect(c("r2s_SNVTA"), names(d))) {
    i <- which(!is.na(d[[col]]) & (d[[col]] < REG$env_r2s[1] | d[[col]] > REG$env_r2s[2]))
    add(i, col, sprintf("R2* envelope [%g, %g] s^-1", REG$env_r2s[1], REG$env_r2s[2]), "set to missing")
    d[[col]][i] <- NA
  }
  dref <- d$chi_SNVTA - d$chi_SNVTA_ref2
  i <- which(!is.na(dref) & abs(dref) > REG$ref2_max_abs)
  add(i, "chi_SNVTA_ref2", "|CSF - whole-brain mean| > 50 ppb", "CSF reference invalid: set to missing (type 'none')")
  d$chi_SNVTA_ref2[i] <- NA; d$ref2_type[i] <- "none"
  list(d = d, log = logt)
}

# Label-free robust z (median / MAD) over the first-usable sessions; |z| > 5 flagged and retained (Q15).
robust_z_counts <- function(first) {
  cols <- intersect(c("chi_SNVTA", "chi_STR", "chi_DN_R", "chi_DGM", "chi_PUT_R", "r2s_SNVTA", "chi_GP"), names(first))
  do.call(rbind, lapply(cols, function(v) {
    x <- first[[v]]; z <- (x - stats::median(x, na.rm = TRUE)) / stats::mad(x, na.rm = TRUE)
    data.frame(roi = v, n = sum(!is.na(x)), n_flag_absz_gt5 = sum(abs(z) > REG$robz, na.rm = TRUE),
               sd_over_robust_sd = stats::sd(x, na.rm = TRUE) / stats::mad(x, na.rm = TRUE))
  }))
}

# Dominance diagnostics for an lm fit (Q15). Returns the rows to drop in the influence refit.
dominance <- function(fit, ids) {
  r <- stats::residuals(fit); sh <- r^2 / sum(r^2); o <- order(sh, decreasing = TRUE)
  ts <- stats::rstudent(fit); cd <- stats::cooks.distance(fit)
  drop <- unique(c(if (sh[o[1]] > REG$dom_share) o[1], which(abs(ts) > REG$dom_t)))
  list(top1_share = unname(sh[o[1]]), top2_share = unname(sum(sh[o[1:2]])), top1_id = ids[o[1]],
       max_abs_rstudent = max(abs(ts), na.rm = TRUE), max_cooks_d = max(cd, na.rm = TRUE),
       triggered = length(drop) > 0, drop_ids = ids[drop])
}

# Huber M-estimator (MASS::rlm, Huber psi k = 1.345, proposal-2 scale, maxit 200), always reported (Q15); gives
# sigma_rob for sigma_B. Registered fallback when rlm fails or does not converge: sigma_rob = 1.4826 * MAD of the OLS
# residuals (flagged in the report); the Huber estimate itself is then reported as unavailable.
huber_fit <- function(formula, data, term, ols_resid = NULL) {
  m <- try_fit(MASS::rlm(formula, data = data, psi = MASS::psi.huber, k = REG$huber_k, scale.est = "proposal 2",
                         maxit = REG$huber_maxit), paste("rlm", deparse(formula[[2]])))
  ok <- !is.null(m) && isTRUE(m$converged)
  if (!ok) {
    sc <- if (!is.null(ols_resid)) 1.4826 * stats::median(abs(ols_resid - stats::median(ols_resid))) else NA
    return(list(est = NA, se = NA, t = NA, scale = sc, converged = FALSE,
                scale_source = "fallback: 1.4826 * MAD(OLS residuals) (rlm failed or did not converge)"))
  }
  cf <- summary(m)$coefficients
  list(est = unname(cf[term, 1]), se = unname(cf[term, 2]), t = unname(cf[term, 3]), scale = m$s,
       converged = TRUE, scale_source = "Huber proposal-2 scale")
}
