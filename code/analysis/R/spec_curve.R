# LEAD-chi analysis: QSM specification curve for H1 and the joint Freedman-Lane permutation test behind the K4
# "robust" label (sec. 24.3; Simonsohn et al. 2020). Registered definition:
#   grid        3 background removal (V-SHARP*, PDF, LBV) x 3 inversion (MEDI*, iLSQR, RTS) x 2 reference (whole
#               brain*, lateral-ventricle CSF) x 3 ROI summaries (primary-ROI weighted mean*, primary-ROI median,
#               the other CIT168 ROI's weighted mean: SN-only if the primary ROI is SN-VTA, SN-VTA under the G3
#               fallback) x 3 covariate sets (full M1*, minimal = PS_e*ADHD_e + ns(age, 3) + sex_e, full + |pitch|)
#               x 3 samples (all first-usable*, without threshold psychosis, visual QC grade 0) = 486 specifications.
#               All 9 reconstructions are summarised on the per-session INTERSECTION of their supports (pipeline
#               extract_rois.py --stage multiverse), so that nodes differ only by algorithm.
#   sessions    always the primary branch's frozen first-usable sessions; a variant value that failed its own
#               plausibility rules, or a missing CSF reference, is missing in that specification only.
#   evaluable   a specification is evaluable if its n >= 30 and >= 0.80 x the number of first-usable subjects in its
#               sample level (all / without threshold psychosis / grade 0). The "robust" label is assessable only if
#               >= 400 of the 486 specifications are evaluable; otherwise it is reported as "not assessable".
#   statistic   the median over evaluable specifications of the specification-wise HC3 t of beta_PS (signed).
#   permutation one vector of i.i.d. uniform keys per permutation is drawn over the first-usable subjects and SHARED by
#               all specifications; within each specification and ADHD stratum the reduced-model (no PS_e) residuals
#               are reassigned by the rank order of the keys (Freedman-Lane), the full model is refitted and its HC3 t
#               recomputed. 10,000 permutations (identity counted); one-sided lower-tail
#               p = (#{median t* <= median t_obs} + 1) / (N + 1).
#   "robust"    final H1 outcome A (A or A+B) AND >= 80 % of evaluable specifications with delta < 0 AND joint p < .05.
#   co-reports  share of evaluable specifications significant in the predicted direction and the Stouffer Z (sum of
#               t-to-z over sqrt(k)), each with its own joint-permutation p (descriptive).

spec_grid <- function() {
  expand.grid(bf = c("vsharp", "pdf", "lbv"), inv = c("medi", "ilsqr", "rts"), ref = c("wb", "csf"),
              summ = c("mean", "median", "other"), cov = c("full", "minimal", "pitch"), samp = c("all", "nothreshold", "grade0"),
              stringsAsFactors = FALSE)
}

spec_rhs <- function(d, cov) {
  switch(cov, full = m1_rhs(d), minimal = "PS_e * ADHD_e + ns(age, df = 3) + sex_e",
         pitch = m1_rhs(d, " + scale(abs_pitch)"))
}

run_spec_curve <- function(S, mv, h1, nperm, seed, primary_roi = "SNVTA") {
  f <- set_race3(S$first); n0 <- nrow(f)
  need <- c("participant_id", "session_id", "bf", "inv", "chi_wb_mean", "chi_wb_median", "chi_other_wb_mean",
            "chi_csf_mean", "chi_csf_median", "chi_other_csf_mean")
  miss <- setdiff(need, names(mv)); if (length(miss)) stop("multiverse table lacks columns: ", paste(miss, collapse = ", "))
  mv$sess_key <- paste0(mv$participant_id, "_", mv$session_id)
  G <- spec_grid()
  nref <- c(all = n0, nothreshold = sum(!(f$dx_psychosis %in% 1)), grade0 = sum(f$qc_grade %in% 0))
  set.seed(seed); U <- matrix(stats::runif(n0 * nperm), n0, nperm)
  idx_cache <- new.env()
  perm_index <- function(keep) {           # keep: logical over rows of f -> n_s x nperm index matrix into 1..n_s
    sig <- paste(which(keep), collapse = ",")
    if (!is.null(idx_cache[[sig]])) return(idx_cache[[sig]])
    rows <- which(keep); ns <- length(rows); I <- matrix(rep(seq_len(ns), nperm), ns, nperm)
    for (st in unique(f$ADHD_e[rows])) {
      loc <- which(f$ADHD_e[rows] == st)
      if (length(loc) > 1) I[loc, ] <- loc[apply(U[rows[loc], , drop = FALSE], 2, order)]
    }
    idx_cache[[sig]] <- I; I
  }
  res <- vector("list", nrow(G)); Tperm <- matrix(NA_real_, nrow(G), nperm)
  for (s in seq_len(nrow(G))) {
    g <- G[s, ]
    col <- paste0("chi_", if (g$summ == "other") "other_" else "", g$ref, "_", if (g$summ == "median") "median" else "mean")
    v <- mv[mv$bf == g$bf & mv$inv == g$inv, c("sess_key", col)]
    d <- f; d$y_spec <- v[[col]][match(d$sess_key, v$sess_key)]
    keep <- !is.na(d$y_spec)
    if (g$samp == "nothreshold") keep <- keep & !(d$dx_psychosis %in% 1)
    if (g$samp == "grade0") keep <- keep & (d$qc_grade %in% 0)
    if (g$cov == "pitch") keep <- keep & !is.na(d$abs_pitch)
    ns <- sum(keep)
    res[[s]] <- data.frame(g, n = ns, evaluable = ns >= 30 && ns >= REG$spec_min_n_frac * nref[[g$samp]], beta = NA, delta = NA, t = NA, df = NA, p_two = NA)
    if (!res[[s]]$evaluable) next
    dd <- set_race3(d[keep, ])
    fit <- try_fit(stats::lm(stats::as.formula(paste("y_spec ~", spec_rhs(dd, g$cov))), data = dd), paste("spec", s))
    if (is.null(fit) || any(is.na(stats::coef(fit)))) { res[[s]]$evaluable <- FALSE; next }
    X <- stats::model.matrix(fit); y <- stats::model.response(stats::model.frame(fit)); j <- match("PS_e", colnames(X))
    Xr <- X[, -j, drop = FALSE]; fr <- drop(Xr %*% solve(crossprod(Xr), crossprod(Xr, y))); er <- y - fr
    XtXi <- solve(crossprod(X)); A <- XtXi %*% t(X); a <- A[j, ]; h <- rowSums((X %*% XtXi) * X)
    tstat <- function(Y) { B <- A %*% Y; R <- Y - X %*% B; B[j, ] / sqrt(colSums(a^2 * (R / (1 - h))^2)) }
    t_obs <- tstat(matrix(y, ncol = 1)); b <- unname(stats::coef(fit)["PS_e"])
    res[[s]]$beta <- b; res[[s]]$delta <- b / stats::sigma(fit); res[[s]]$t <- t_obs
    res[[s]]$df <- stats::df.residual(fit); res[[s]]$p_two <- 2 * stats::pt(-abs(t_obs), stats::df.residual(fit))
    I <- perm_index(keep)
    Tperm[s, ] <- tstat(fr + matrix(er[I], nrow(I), ncol(I)))
  }
  tab <- do.call(rbind, res); ev <- tab$evaluable; k <- sum(ev)
  zobs <- stats::qnorm(stats::pt(tab$t[ev], tab$df[ev]))
  stat_obs <- stats::median(tab$t[ev]); share_neg <- mean(tab$delta[ev] < 0)
  share_sig <- mean(tab$p_two[ev] < REG$alpha & tab$beta[ev] < 0); stouffer <- sum(zobs) / sqrt(k)
  Tp <- Tperm[ev, , drop = FALSE]
  med_p <- apply(Tp, 2, stats::median); sig_p <- colMeans(2 * stats::pt(-abs(Tp), tab$df[ev]) < REG$alpha & Tp < 0)
  st_p <- colSums(stats::qnorm(stats::pt(Tp, tab$df[ev]))) / sqrt(k)
  p_joint <- (sum(med_p <= stat_obs + 1e-12) + 1) / (nperm + 1)
  assessable <- k >= REG$spec_min_evaluable
  A_final <- h1$final_outcome %in% c("A", "A+B")
  robust <- assessable && A_final && share_neg >= REG$robust_share && p_joint < REG$alpha
  list(n_specifications = nrow(G), n_evaluable = k, assessable = assessable, primary_roi = primary_roi,
       median_t = stat_obs, median_delta = stats::median(tab$delta[ev]), share_delta_negative = share_neg,
       share_significant_predicted = share_sig, stouffer_Z = stouffer,
       joint = c(list(statistic = "median specification-wise HC3 t", alternative = "less"), mc_info(p_joint, nperm)),
       joint_share_significant_p = (sum(sig_p >= share_sig - 1e-12) + 1) / (nperm + 1),
       joint_stouffer_p = (sum(st_p <= stouffer + 1e-12) + 1) / (nperm + 1),
       robust = robust,
       robust_reason = if (!assessable) "not assessable (< 400 evaluable specifications)" else if (!A_final) "final H1 outcome is not A" else
         if (share_neg < REG$robust_share) "< 80 % of specifications with delta < 0" else if (p_joint >= REG$alpha) "joint permutation p >= .05" else "all criteria met",
       table = tab)
}
