# LEAD-chi analysis: T1 temporal-stability test (sec. 4.4, Q19) with the REML-maximum rule, and the pooled
# all-session H1 sensitivity. Registered inference (Q19, Q23; no pbkrtest / clubSandwich, none installed):
#   decision  = Wald t of beta(PS_e:occ) at the verified REML maximum with nlme containment df;
#   co-report = parametric bootstrap p (lme4, 2,000 draws; never decisive) and an OLS working-independence fit with
#               sandwich::vcovCL(type = "HC2") cluster-robust SE plus a restricted wild cluster bootstrap (Rademacher).
#               The latter is NOT CR2 for the LMM estimate; it is labelled "OLS cluster-robust (HC2)".

# Fit a random-intercept LMM with four optimisers and report the fit at the verified REML maximum.
reml_max <- function(fixed, data, profile = TRUE) {
  f_lmer <- stats::update(fixed, . ~ . + (1 | sub))
  fits <- list(
    lmer_bobyqa = try_fit(lme4::lmer(f_lmer, data = data, REML = TRUE, control = lme4::lmerControl(optimizer = "bobyqa")), "lmer bobyqa"),
    lmer_NelderMead = try_fit(lme4::lmer(f_lmer, data = data, REML = TRUE, control = lme4::lmerControl(optimizer = "Nelder_Mead")), "lmer Nelder_Mead"),
    lme_nlminb = try_fit(nlme::lme(fixed, random = ~ 1 | sub, data = data, method = "REML",
                                   control = nlme::lmeControl(opt = "nlminb", maxIter = 200, msMaxIter = 200, niterEM = 50)), "lme nlminb"),
    lme_optim = try_fit(nlme::lme(fixed, random = ~ 1 | sub, data = data, method = "REML",
                                  control = nlme::lmeControl(opt = "optim", maxIter = 500, msMaxIter = 500)), "lme optim"))
  dev <- lme4::lmer(f_lmer, data = data, REML = TRUE, devFunOnly = TRUE)
  theta_of <- function(f) if (inherits(f, "merMod")) unname(lme4::getME(f, "theta")) else {
    vc <- nlme::VarCorr(f); sqrt(as.numeric(vc[1, "Variance"])) / f$sigma }
  ok <- !vapply(fits, is.null, logical(1))
  th <- vapply(fits[ok], theta_of, numeric(1)); dv <- vapply(th, dev, numeric(1))
  best <- names(dv)[which.min(dv + seq_along(dv) * 1e-10)]   # ties -> registered order
  tab <- data.frame(fit = names(dv), theta = th, reml_deviance = dv, gap_to_best = dv - min(dv), row.names = NULL)
  prof <- NULL
  if (profile) {
    g <- sort(unique(c(REG$reml_theta_grid, th))); pd <- vapply(g, dev, numeric(1))
    prof <- list(theta_grid_argmin = g[which.min(pd)], dev_grid_min = min(pd), dev_best = min(dv),
                 best_is_global_on_grid = min(dv) <= min(pd) + 1e-8, dev_at_theta0 = dev(0))
    if (!prof$best_is_global_on_grid) {   # registered action: refit from the grid argmin and keep the better fit
      fp <- try_fit(lme4::lmer(f_lmer, data = data, REML = TRUE, start = list(theta = prof$theta_grid_argmin),
                               control = lme4::lmerControl(optimizer = "bobyqa")), "lmer from profile argmin")
      if (!is.null(fp)) { fits$lmer_profile_start <- fp; ok <- !vapply(fits, is.null, logical(1))
        th <- vapply(fits[ok], theta_of, numeric(1)); dv <- vapply(th, dev, numeric(1))
        best <- names(dv)[which.min(dv + seq_along(dv) * 1e-10)]
        tab <- data.frame(fit = names(dv), theta = th, reml_deviance = dv, gap_to_best = dv - min(dv), row.names = NULL) }
      prof$refit_from_profile <- TRUE
      prof$best_is_global_after_refit <- min(dv) <= min(pd) + 1e-8
    }
  }
  fb <- fits[[best]]
  if (inherits(fb, "merMod")) { b <- lme4::fixef(fb); V <- as.matrix(stats::vcov(fb)); vc <- as.data.frame(lme4::VarCorr(fb)); s2u <- vc$vcov[1]; s2e <- vc$vcov[2] }
  else { b <- nlme::fixef(fb); V <- stats::vcov(fb); s2u <- as.numeric(nlme::VarCorr(fb)[1, "Variance"]); s2e <- fb$sigma^2 }
  nl <- if (!is.null(fits$lme_nlminb)) fits$lme_nlminb else fits$lme_optim
  cdf <- if (!is.null(nl)) summary(nl)$tTable[, "DF"] else setNames(rep(nrow(data) - length(b), length(b)), names(b))
  list(best = best, table = tab, profile = prof, b = b, V = V, se = sqrt(diag(V)), containment_df = cdf,
       s2_subject = s2u, s2_resid = s2e, icc = s2u / (s2u + s2e), boundary = min(th[best]) < 1e-4,
       fit_lmer = fits$lmer_bobyqa, fits = fits)
}

# CR2 (Bell-McCaffrey) cluster-robust variance, identical to sandwich::vcovCL(type = "HC2") for OLS: in sandwich the
# G/(G-1) cluster adjustment cancels against its sqrt((G-1)/G) factor for HC2/HC3, so no adjustment is applied here.
cr2_setup <- function(X, cl) {
  XtXi <- solve(crossprod(X)); gs <- split(seq_len(nrow(X)), cl)
  M <- lapply(gs, function(ix) {
    Xg <- X[ix, , drop = FALSE]; Hg <- Xg %*% XtXi %*% t(Xg); e <- eigen(diag(length(ix)) - Hg, symmetric = TRUE)
    Ag <- e$vectors %*% diag(1 / sqrt(pmax(e$values, 1e-12)), length(ix)) %*% t(e$vectors)
    XtXi %*% t(Xg) %*% Ag })
  list(gs = gs, M = M, adj = 1, XtXi = XtXi)
}
cr2_se <- function(cs, j, R) {   # R: n x B residual matrix -> vector of B standard errors for coefficient j
  s <- 0
  for (g in seq_along(cs$gs)) s <- s + (cs$M[[g]][j, , drop = FALSE] %*% R[cs$gs[[g]], , drop = FALSE])^2
  sqrt(cs$adj * drop(s))
}

# Wild cluster (restricted, Rademacher) bootstrap-t for H0: beta_j = b0 with the CR2-type SE.
wild_cluster <- function(X, y, cl, term, b0, B, seed, alternative = c("less", "two.sided")) {
  alternative <- match.arg(alternative); j <- match(term, colnames(X))
  cs <- cr2_setup(X, cl); A <- cs$XtXi %*% t(X)
  b <- drop(A %*% y); r <- y - drop(X %*% b); t_obs <- (b[j] - b0) / cr2_se(cs, j, matrix(r))
  yr <- y - b0 * X[, j]; Xr <- X[, -j, drop = FALSE]; br <- solve(crossprod(Xr), crossprod(Xr, yr))
  fr <- drop(Xr %*% br) + b0 * X[, j]; ur <- yr - drop(Xr %*% br)
  set.seed(seed); ug <- unique(cl); cnt <- 0; done <- 0
  while (done < B) {
    m <- min(2000, B - done); W <- matrix(sample(c(-1, 1), length(ug) * m, replace = TRUE), length(ug), m)
    Ys <- fr + ur * W[match(cl, ug), , drop = FALSE]
    Bs <- A %*% Ys; Rs <- Ys - X %*% Bs; ts <- (Bs[j, ] - b0) / cr2_se(cs, j, Rs)
    cnt <- cnt + if (alternative == "less") sum(ts <= t_obs + 1e-12) else sum(abs(ts) >= abs(t_obs) - 1e-12)
    done <- done + m
  }
  c(list(t_obs = unname(t_obs), b0 = b0, alternative = alternative), mc_info((cnt + 1) / (B + 1), B))
}

# Parametric bootstrap (lme4) of the Wald t for H0: beta_term = b0, simulating from the null-restricted REML fit.
pb_test <- function(fixed, data, term_col, b0, nsim, seed, alternative = c("less", "two.sided")) {
  alternative <- match.arg(alternative)
  full <- lme4::lmer(stats::update(fixed, . ~ . + (1 | sub)), data = data, REML = TRUE, control = lme4::lmerControl(optimizer = "bobyqa"))
  t_obs <- (lme4::fixef(full)[term_col] - b0) / sqrt(stats::vcov(full)[term_col, term_col])
  dn <- data; yv <- all.vars(fixed)[1]; dn$.y_r <- dn[[yv]] - b0 * dn[[term_col]]
  fn <- stats::update(fixed, stats::as.formula(paste(".y_r ~ . -", term_col)))
  null <- lme4::lmer(stats::update(fn, . ~ . + (1 | sub)), data = dn, REML = TRUE, control = lme4::lmerControl(optimizer = "bobyqa"))
  set.seed(seed); sims <- stats::simulate(null, nsim = nsim)
  ts <- vapply(seq_len(nsim), function(k) {
    f <- tryCatch(suppressWarnings(suppressMessages(lme4::refit(full, newresp = sims[[k]] + b0 * data[[term_col]]))),
                  error = function(e) NULL)
    if (is.null(f)) return(NA_real_)
    (lme4::fixef(f)[term_col] - b0) / sqrt(stats::vcov(f)[term_col, term_col]) }, numeric(1))
  ts <- ts[is.finite(ts)]
  cnt <- if (alternative == "less") sum(ts <= t_obs + 1e-12) else sum(abs(ts) >= abs(t_obs) - 1e-12)
  c(list(t_obs = unname(t_obs), b0 = b0, alternative = alternative, n_valid = length(ts)), mc_info((cnt + 1) / (length(ts) + 1), length(ts)))
}

t1_fixed <- function(d) stats::as.formula(paste0("y_tsd ~ PS_e * ADHD_e + occ + PSxocc + ns(age, 3) + sex_e", race_term(d),
                                                 " + scale(roi_vol) + scale(Q_phase)"))

t1_core <- function(two, sd_total, n_pb, seed, do_pb = TRUE, profile = TRUE) {
  two <- set_race3(two); two$y_tsd <- two$chi_SNVTA / sd_total; two$PSxocc <- two$PS_e * two$occ
  fx <- t1_fixed(two); rm <- reml_max(fx, two, profile = profile)
  b <- unname(rm$b["PSxocc"]); se <- unname(rm$se["PSxocc"]); df <- unname(rm$containment_df["PSxocc"])
  tm <- (b - REG$t1_margin) / se; t0 <- b / se
  res <- list(n_rows = nrow(two), n_subjects = length(unique(two$sub)), n_returners = sum(two$occ == 1),
              reml = rm[c("best", "table", "profile", "icc", "s2_subject", "s2_resid", "boundary")],
              beta_PSxocc = b, se = se, containment_df = df,
              margin = list(Delta = REG$t1_margin, Delta_ppb = REG$t1_margin * sd_total, t = tm, p_containment = stats::pt(tm, df)),
              decline = list(t = t0, p_containment = 2 * stats::pt(-abs(t0), df)))
  if (do_pb) {
    pbm <- pb_test(fx, two, "PSxocc", REG$t1_margin, n_pb, seed + 21, "less")
    if (pbm$in_band && n_pb < REG$mc_recompute) { log_msg("T1 PB margin p in MC band: recomputing with ", REG$mc_recompute)
      pbm <- pb_test(fx, two, "PSxocc", REG$t1_margin, REG$mc_recompute, seed + 22, "less") }
    pb0 <- pb_test(fx, two, "PSxocc", 0, n_pb, seed + 23, "two.sided")
    if (pb0$in_band && n_pb < REG$mc_recompute) pb0 <- pb_test(fx, two, "PSxocc", 0, REG$mc_recompute, seed + 24, "two.sided")
    res$margin$pb <- pbm; res$decline$pb <- pb0                      # co-reported, never decisive
    res$margin$p_decision <- res$margin$p_containment                # registered decision: containment df
    res$decline$p_decision <- res$decline$p_containment
    # OLS + CR2-type cluster-robust co-report and wild cluster bootstrap
    ols <- stats::lm(fx, data = two); X <- stats::model.matrix(ols); y <- stats::model.response(stats::model.frame(ols))
    Vc <- sandwich::vcovCL(ols, cluster = two$sub, type = "HC2")
    bo <- unname(stats::coef(ols)["PSxocc"]); so <- sqrt(Vc["PSxocc", "PSxocc"])
    wc <- wild_cluster(X, y, two$sub, "PSxocc", REG$t1_margin, REG$n_wild, seed + 25, "less")
    chk <- abs(wc$t_obs - (bo - REG$t1_margin) / so); if (chk > 1e-6) log_msg("WARNING: CR2 algebra differs from vcovCL by ", chk)
    res$cluster_robust <- list(ols_beta = bo, se_vcovCL_HC2 = so, t_margin = (bo - REG$t1_margin) / so, wild_margin = wc)
  } else {
    res$margin$p_decision <- res$margin$p_containment; res$decline$p_decision <- res$decline$p_containment
  }
  res
}

run_t1 <- function(S, h1, gates, n_pb, seed) {
  # scale: pooled label-free SD of chi over the frozen first-usable sessions (gate JSON; required)
  sd_total <- gates$pooled$sd_total_SNVTA
  core <- t1_core(S$two_occ, sd_total, n_pb, seed)
  # follow-up-only M1 (about 45-53 subjects): only the sign of beta_PS is used; an undefined HC3 SE never stops the run
  fu <- set_race3(S$occ1)
  mfu <- try_fit(fit_ols_hc3(m1_formula("chi_SNVTA", fu), fu), "T1 follow-up-only M1")
  tfu <- if (is.null(mfu)) { b <- unname(stats::coef(stats::lm(m1_formula("chi_SNVTA", fu), data = fu))["PS_e"])
                             list(est = b, se = NA, p_two = NA) } else term_test(mfu, "PS_e")
  if (is.null(mfu)) mfu <- list(n = nrow(fu))
  A_holds <- h1$final_outcome %in% c("A", "A+B")
  not_decl <- isTRUE(core$margin$p_decision < REG$alpha && tfu$est < 0)
  declined <- isTRUE(core$beta_PSxocc > 0 && core$decline$p_decision < REG$alpha && tfu$est >= 0)
  label <- if (!A_holds) "not evaluated (outcome A does not hold): follow-up-only estimate reported descriptively" else
    if (not_decl) "did not decline by >= Delta" else if (declined) "declined" else "indeterminate"
  # robust ppb margin co-report (v1.2 sec. 5.5.4)
  x <- S$first$chi_SNVTA; ratio <- stats::sd(x) / (stats::mad(x))
  robust <- if (ratio > REG$robust_ratio) {
    r <- t1_core(S$two_occ, stats::mad(x), n_pb, seed, do_pb = FALSE, profile = FALSE)
    list(ratio = ratio, beta = r$beta_PSxocc, p_margin_containment = r$margin$p_containment, Delta_ppb_robust = REG$t1_margin * stats::mad(x))
  } else list(ratio = ratio, note = "SD / (1.4826 MAD) <= 1.5: registered scale only")
  # exposure sensitivity (E1 sec. 6.4(4)): drop listed session rows; containment-df test only
  te <- S$two_occ[!S$two_occ$sess_key %in% REG$exposure_sessions, ]
  ex <- t1_core(te, sd_total, n_pb, seed, do_pb = FALSE, profile = FALSE)
  # occasion-0 PS estimate in total-SD units, reported next to Delta (sec. 4.4: "stable"/"persistent" wording prohibited)
  occ0 <- h1$decision$beta_ppb / sd_total
  c(core, list(sd_total = sd_total, occasion0_beta_tsd = occ0,
               Delta_vs_occasion0 = if (is.finite(occ0) && abs(occ0) <= REG$t1_margin)
                 "the occasion-0 deficit is <= Delta: complete disappearance lies inside the alternative" else
                 "the occasion-0 deficit exceeds Delta",
               assumption = "two-occasion LMM; missing follow-ups assumed missing at random given the model covariates",
               followup_only = list(n = mfu$n, cells = cell_counts(fu), beta_PS = tfu$est, se = tfu$se, p_two = tfu$p_two,
                                    delta_tsd = tfu$est / sd_total),
               A_holds = A_holds, label = label, robust_scale = robust,
               exposure_sensitivity = list(dropped = intersect(REG$exposure_sessions, S$two_occ$sess_key),
                                           n_dropped = length(intersect(REG$exposure_sessions, S$two_occ$sess_key)), n_rows = ex$n_rows,
                                           beta = ex$beta_PSxocc, p_margin_containment = ex$margin$p_containment)))
}

# Sec. 13.3 pooled all-session sensitivity for H1 (not a test): between/within age split, REML maximum.
run_pooled <- function(S, n_wild, seed) {
  d <- set_race3(S$usable); d$age_bar <- stats::ave(d$age, d$sub); d$age_w <- d$age - d$age_bar
  fx <- stats::as.formula(paste0("chi_SNVTA ~ PS_e * ADHD_e + ns(age_bar, 3) + age_w + sex_e", race_term(d),
                                 " + scale(roi_vol) + scale(Q_phase)"))
  rm <- reml_max(fx, d)
  ols <- stats::lm(fx, data = d); X <- stats::model.matrix(ols); y <- stats::model.response(stats::model.frame(ols))
  Vc <- sandwich::vcovCL(ols, cluster = d$sub, type = "HC2")
  wc <- wild_cluster(X, y, d$sub, "PS_e", 0, n_wild, seed + 31, "two.sided")
  chk <- abs(wc$t_obs - stats::coef(ols)["PS_e"] / sqrt(Vc["PS_e", "PS_e"])); if (chk > 1e-6) log_msg("WARNING: CR2 algebra differs from vcovCL by ", chk)
  list(n_sessions = nrow(d), n_subjects = length(unique(d$sub)), reml = rm[c("best", "table", "profile", "icc", "boundary")],
       beta_PS = unname(rm$b["PS_e"]), se = unname(rm$se["PS_e"]), containment_df = unname(rm$containment_df["PS_e"]),
       p_containment = 2 * stats::pt(-abs(rm$b["PS_e"] / rm$se["PS_e"]), rm$containment_df["PS_e"]),
       ols_cluster_robust = list(beta = unname(stats::coef(ols)["PS_e"]), se_vcovCL_HC2 = sqrt(Vc["PS_e", "PS_e"]), wild = wc))
}
