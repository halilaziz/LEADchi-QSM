# LEAD-chi analysis: primary model M1, decision rule A/B/C, qualifiers K1-K4 and the H1 co-reports
# (registered text sec. 4.1-4.2, Q13, Q15, Q19; Annex H maps each rule to the function below).

m1_rhs <- function(d, extra = "")
  paste0("PS_e * ADHD_e + ns(age, df = 3) + sex_e", race_term(d), " + scale(roi_vol) + scale(Q_phase)", extra)
m1_formula <- function(y, d, extra = "") stats::as.formula(paste(y, "~", m1_rhs(d, extra)))

# OLS with HC3 (sandwich::vcovHC type "HC3"); t reference with df = n - p. If a leverage equals 1 the HC3 SE of the
# affected terms is undefined (NaN): the estimate is kept and the SE reported as NA (never an error; Q19 T1 note).
fit_ols_hc3 <- function(formula, data, weights = NULL) {
  fit <- if (is.null(weights)) stats::lm(formula, data = data) else stats::lm(formula, data = data, weights = weights)
  if (any(is.na(stats::coef(fit)))) stop("rank-deficient design: ", paste(names(which(is.na(stats::coef(fit)))), collapse = ", "))
  V <- suppressWarnings(sandwich::vcovHC(fit, type = "HC3"))
  list(fit = fit, V = V, b = stats::coef(fit), se = sqrt(diag(V)), df = stats::df.residual(fit),
       sigma = stats::sigma(fit), n = stats::nobs(fit), formula = formula)
}

term_test <- function(m, term, level = 0.95) {
  b <- unname(m$b[term]); se <- unname(m$se[term]); t <- b / se; q <- stats::qt(1 - (1 - level) / 2, m$df)
  list(term = term, est = b, se = se, t = t, df = m$df, p_two = 2 * stats::pt(-abs(t), m$df),
       p_lower = stats::pt(t, m$df), p_upper = stats::pt(t, m$df, lower.tail = FALSE),
       ci95 = c(b - q * se, b + q * se), ci90 = c(b - stats::qt(.95, m$df) * se, b + stats::qt(.95, m$df) * se))
}

lin_test <- function(m, cvec, label) {
  cv <- setNames(rep(0, length(m$b)), names(m$b)); cv[names(cvec)] <- cvec
  est <- sum(cv * m$b); se <- sqrt(drop(t(cv) %*% m$V %*% cv)); t <- est / se; q <- stats::qt(.975, m$df)
  list(contrast = label, est = est, se = se, t = t, df = m$df, p_two = 2 * stats::pt(-abs(t), m$df),
       p_lower = stats::pt(t, m$df), p_upper = stats::pt(t, m$df, lower.tail = FALSE), ci95 = c(est - q * se, est + q * se))
}

coef_table <- function(m) data.frame(term = names(m$b), est = unname(m$b), se_hc3 = unname(m$se),
                                     t = unname(m$b / m$se), p_two = 2 * stats::pt(-abs(m$b / m$se), m$df), row.names = NULL)

# Outcome-B bound test at a given delta (in residual-SD units) with scale sigma_B.
bound_test <- function(tt, df, sigma_B, db) {
  tB <- (tt$est + db * sigma_B) / tt$se
  list(delta_B = db, bound_ppb = -db * sigma_B, tB = tB, pB = stats::pt(tB, df, lower.tail = FALSE))
}

# Decision rule (sec. 4.1). A = beta < 0 and two-sided HC3 p < .05. B = one-sided rejection of delta <= -delta_B with
# delta_B = 0.6 * sqrt(r_L), sigma_B = min(sigma_resid, sigma_rob). B is EVALUABLE only if PC4 passes and
# r_L >= 0.25; otherwise B is not declared (label "low-precision measure") and the bound test is descriptive only.
h1_decide <- function(m, sigma_rob, rL, pc4 = TRUE, term = "PS_e") {
  tt <- term_test(m, term)
  sigma_B <- min(c(m$sigma, sigma_rob), na.rm = TRUE)
  rL <- as.numeric(rL)
  B_evaluable <- isTRUE(pc4) && is.finite(rL) && rL >= REG$rL_floor && rL <= 1
  dB <- if (is.finite(rL) && rL > 0) REG$deltaB_base * sqrt(min(1, rL)) else NA_real_
  B_reg <- if (is.na(dB)) list(delta_B = NA, bound_ppb = NA, tB = NA, pB = NA) else bound_test(tt, m$df, sigma_B, dB)
  A <- isTRUE(tt$est < 0 && tt$p_two < REG$alpha)
  B <- B_evaluable && isTRUE(B_reg$pB < REG$alpha)
  outcome <- if (A && B) "A+B" else if (A) "A" else if (B) "B" else "C"
  list(beta_ppb = tt$est, se_hc3 = tt$se, t = tt$t, df = m$df, p_two = tt$p_two, ci95_ppb = tt$ci95, ci90_ppb = tt$ci90,
       sigma_resid = m$sigma, sigma_rob = sigma_rob, sigma_B = sigma_B, delta = tt$est / m$sigma,
       rL = rL, B_evaluable = B_evaluable, B_registered = B_reg,
       B_coreport = lapply(REG$deltaB_coreport, function(db) bound_test(tt, m$df, sigma_B, db)),
       A = A, B = B, outcome = outcome, disconfirming = isTRUE(tt$est > 0 && tt$p_two < REG$alpha),
       tcrit_one_sided = stats::qt(1 - REG$alpha, m$df))
}

# Registered sentences (sec. 4.6). The report prints only these.
REGISTERED_SENTENCES <- list(
  A = paste("Outcome A: a lower SN-VTA susceptibility is present in clinically ascertained, help-seeking",
            "psychosis-spectrum youth as a group (averaged over ADHD status)."),
  AB = paste("Outcomes A and B: a lower SN-VTA susceptibility is present in help-seeking psychosis-spectrum youth as a group,",
             "and it is smaller than the adult-magnitude bound."),
  B = paste("Outcome B: An adult-magnitude SN-VTA susceptibility deficit is not a feature of clinically ascertained,",
            "help-seeking psychosis-spectrum youth as a group. This does not exclude a deficit confined to the minority",
            "who will develop psychosis, nor a deficit that emerges later."),
  C = paste("Outcome C: the data neither establish nor exclude an SN-VTA susceptibility deficit of the tested size",
            "in help-seeking psychosis-spectrum youth; estimates are offered for meta-analysis."),
  disconfirming = paste("Disconfirming H1: SN-VTA susceptibility is higher, not lower, in help-seeking psychosis-spectrum",
                        "youth (beta_PS > 0, two-sided p < .05)."),
  K1a = "The lower susceptibility is seen in psychosis-spectrum youth both with and without ADHD (point estimates).",
  K1b = paste("The lower susceptibility is attributed to psychosis-spectrum status and not to ADHD: the ADHD effect is",
              "equivalent to zero within +/-0.5 SD and the interaction is not significant (non-rejection only)."),
  K1c = paste("The lower susceptibility is concentrated in psychosis-spectrum youth with comorbid ADHD; no",
              "'psychosis-spectrum signature' is claimed."),
  K2_iron = "Lower R2* accompanies the lower susceptibility; 'lower iron' wording is used.",
  K2_notiron = "Lower susceptibility; the iron-versus-myelin source is unresolved.",
  low_precision = paste("The SN-VTA measure did not reach the registered precision (PC4 or r_L below the floor); outcome B",
                        "was not evaluated and H1 is labelled 'low-precision measure'."),
  fragile = "H1 is labelled 'fragile': deleting a single subject changes the A/B/C outcome.",
  influence = "H1 is labelled 'influence-dependent': removing the flagged dominant subject(s) changes the A/B/C outcome.",
  robust = "H1 is labelled 'robust' across the registered QSM specification curve.")

outcome_sentence <- function(o, disconfirming = FALSE) {
  if (isTRUE(disconfirming)) return(REGISTERED_SENTENCES$disconfirming)
  switch(o, "A" = REGISTERED_SENTENCES$A, "A+B" = REGISTERED_SENTENCES$AB, "B" = REGISTERED_SENTENCES$B,
         "C" = REGISTERED_SENTENCES$C)
}

# Freedman-Lane permutation of a single coefficient's HC3 t (Winkler et al. 2014). Reduced-model residuals are
# permuted within strata and added back to the reduced fit; the full model is refitted (vectorised OLS algebra,
# checked against sandwich::vcovHC on the observed data). p = (#{|t*| >= |t_obs|} + 1) / (N + 1).
fl_perm <- function(fit, term, strata, nperm, seed, chunk = 2000) {
  X <- stats::model.matrix(fit); y <- stats::model.response(stats::model.frame(fit)); n <- length(y)
  j <- match(term, colnames(X)); if (is.na(j)) stop("term not in model: ", term)
  Xr <- X[, -j, drop = FALSE]
  fr <- drop(Xr %*% solve(crossprod(Xr), crossprod(Xr, y))); er <- y - fr
  XtXi <- solve(crossprod(X)); A <- XtXi %*% t(X); a <- A[j, ]; h <- rowSums((X %*% XtXi) * X)
  tstat <- function(Y) { B <- A %*% Y; R <- Y - X %*% B; B[j, ] / sqrt(colSums(a^2 * (R / (1 - h))^2)) }
  t_obs <- tstat(matrix(y, ncol = 1))
  t_chk <- unname(stats::coef(fit)[term] / sqrt(sandwich::vcovHC(fit, type = "HC3")[term, term]))
  if (abs(t_obs - t_chk) > 1e-8 * max(1, abs(t_chk))) stop("Freedman-Lane HC3 algebra check failed")
  set.seed(seed); strata <- as.character(strata); cnt <- 0; done <- 0
  while (done < nperm) {
    m <- min(chunk, nperm - done); idx <- matrix(rep(seq_len(n), m), n, m)
    for (s in unique(strata)) { r <- which(strata == s); if (length(r) > 1) idx[r, ] <- apply(idx[r, , drop = FALSE], 2, function(v) v[sample.int(length(v))]) }
    Ys <- fr + matrix(er[idx], n, m)
    cnt <- cnt + sum(abs(tstat(Ys)) >= abs(t_obs) - 1e-12); done <- done + m
  }
  p <- (cnt + 1) / (nperm + 1)
  c(list(t_obs = t_obs, strata = "ADHD_e"), mc_info(p, nperm))
}

# Cell-stratified bootstrap CI of delta = beta / sigma_resid, percentile.
boot_delta <- function(formula, data, term, nboot, seed) {
  stat <- function(dd, i) { f <- stats::lm(formula, data = dd[i, ]); b <- stats::coef(f)[term]; if (is.na(b)) NA else b / stats::sigma(f) }
  set.seed(seed)
  bt <- boot::boot(data, stat, R = nboot, strata = factor(data$cell))
  v <- bt$t[, 1]; v <- v[is.finite(v)]
  list(n_boot = nboot, n_valid = length(v), ci95 = unname(stats::quantile(v, c(.025, .975))), strata = "cell")
}

# Analytic normal-likelihood Bayes factors on delta (brms fallback; Dienes 2014), and sceptical-prior shrinkage.
bayes_normal <- function(delta_hat, se_delta) {
  bfs <- lapply(REG$bf_priors, function(pr) {
    m1 <- stats::dnorm(delta_hat, pr["mean"], sqrt(se_delta^2 + pr["sd"]^2)); m0 <- stats::dnorm(delta_hat, 0, se_delta)
    list(prior = paste0("N(", pr["mean"], ", ", pr["sd"], ")"), BF10 = unname(m1 / m0), BF01 = unname(m0 / m1))
  })
  s2 <- REG$bf_priors$sceptical["sd"]^2
  post_var <- 1 / (1 / se_delta^2 + 1 / s2); post_mean <- post_var * (delta_hat / se_delta^2 + REG$bf_priors$sceptical["mean"] / s2)
  list(BF = bfs, shrunken_sceptical = list(mean = unname(post_mean), sd = unname(sqrt(post_var)),
                                            ci95 = unname(post_mean + c(-1, 1) * 1.959964 * sqrt(post_var))))
}

# Noncentral-t CI for delta (MBESS fallback: root finding on the noncentrality parameter).
nct_ci_delta <- function(t, df, se, sigma, level = 0.95) {
  f <- function(ncp, target) suppressWarnings(stats::pt(t, df, ncp)) - target
  lo <- try_fit(stats::uniroot(f, c(t - 10, t + 10), target = 1 - (1 - level) / 2, tol = 1e-10)$root, "nct lower")
  hi <- try_fit(stats::uniroot(f, c(t - 10, t + 10), target = (1 - level) / 2, tol = 1e-10)$root, "nct upper")
  if (is.null(lo) || is.null(hi)) return(c(NA, NA))
  c(lo, hi) * se / sigma
}

hedges_g_raw <- function(y, g) {
  a <- y[g > 0]; b <- y[g < 0]; na <- length(a); nb <- length(b)
  sp <- sqrt(((na - 1) * stats::var(a) + (nb - 1) * stats::var(b)) / (na + nb - 2))
  (mean(a) - mean(b)) / sp * (1 - 3 / (4 * (na + nb) - 9))
}

# K1 labels (sec. 4.2; intersection-union, no extra alpha). Precedence: K1c > K1b > K1a.
#   K1a "PS effect": A, and both simple effects < 0.
#   K1b "PS-specific (not ADHD)": K1a, plus TOST on beta_ADHD at +/- 0.5 * sigma_B (both one-sided at .05), plus the
#       interaction 95% CI includes 0 (a non-rejection; no evidence of absent interaction is claimed).
#   K1c "Concentrated in comorbid youth": A; est(PS+ADHD - ADHD-only) < 0 with p < .05; PS-only - TD p >= .05 and
#       est(PS-only - TD) > 0.5 * est(PS+ADHD - ADHD-only).
k1_labels <- function(m, dec) {
  se1 <- lin_test(m, c(PS_e = 1, `PS_e:ADHD_e` = -0.5), "PS-only - TD")
  se2 <- lin_test(m, c(PS_e = 1, `PS_e:ADHD_e` = 0.5), "PS+ADHD - ADHD-only")
  ta <- term_test(m, "ADHD_e"); ti <- term_test(m, "PS_e:ADHD_e"); tc <- stats::qt(1 - REG$alpha, m$df)
  bnd <- REG$tost_bound * dec$sigma_B
  tost_lo <- (ta$est + bnd) / ta$se; tost_hi <- (ta$est - bnd) / ta$se
  tost <- isTRUE((tost_lo > tc) && (tost_hi < -tc))
  int_ci0 <- isTRUE(ti$ci95[1] <= 0 && ti$ci95[2] >= 0)
  ps_effect <- isTRUE(dec$A && se1$est < 0 && se2$est < 0)
  ps_specific <- ps_effect && tost && int_ci0
  comorbid <- isTRUE(dec$A && se2$est < 0 && se2$p_two < REG$alpha && se1$p_two >= REG$alpha &&
                     se1$est > REG$comorbid_ratio * se2$est)
  label <- if (comorbid) "K1c Concentrated in comorbid youth" else if (ps_specific) "K1b PS-specific (not ADHD)" else
    if (ps_effect) "K1a PS effect" else "none"
  list(simple_PSonly_minus_TD = se1, simple_PSADHD_minus_ADHDonly = se2, ADHD = ta, interaction = ti,
       tost = list(bound_sd = REG$tost_bound, scale = "sigma_B", bound_ppb = bnd, t_lower = tost_lo, t_upper = tost_hi,
                   p_lower = stats::pt(tost_lo, m$df, lower.tail = FALSE), p_upper = stats::pt(tost_hi, m$df), success = tost),
       interaction_ci_includes_0 = int_ci0, ps_effect = ps_effect, ps_specific = ps_specific, comorbid = comorbid,
       label = label, ps_signature_wording_prohibited = comorbid)
}

# M1 core (one analysed sample): fit, Huber scale (with the registered fallback) and decision.
run_m1_core <- function(first, rL, pc4 = TRUE, y = "chi_SNVTA") {
  first <- set_race3(first)
  fo <- m1_formula(y, first); m <- fit_ols_hc3(fo, first)
  hb <- huber_fit(fo, first, "PS_e", ols_resid = stats::residuals(m$fit))
  dec <- h1_decide(m, hb$scale, rL, pc4)
  list(m = m, formula = fo, huber = hb, dec = dec, data = first)
}

# Split-echo reliability of covariate-residualised chi (post-unlock co-report of the B bound; sec. 4.1, M3 of R_stats).
resid_split_reliability <- function(first) {
  if (!all(c("chi_SNVTA_split13", "chi_SNVTA_split24") %in% names(first))) return(NULL)
  d <- set_race3(first[!is.na(first$chi_SNVTA_split13) & !is.na(first$chi_SNVTA_split24), ])
  if (nrow(d) < 20) return(NULL)
  r1 <- stats::residuals(stats::lm(m1_formula("chi_SNVTA_split13", d), data = d))
  r2 <- stats::residuals(stats::lm(m1_formula("chi_SNVTA_split24", d), data = d))
  r <- stats::cor(r1, r2); z <- atanh(r); lo <- tanh(z - 1.959964 / sqrt(nrow(d) - 3))
  sb <- function(x) 2 * x / (1 + x)
  list(n = nrow(d), r = r, R_SB = sb(r), rL = sb(lo))
}

run_h1 <- function(S, gates, nperm, nboot, seed) {
  first <- set_race3(S$first); rL <- gates$rL_primary; pc4 <- isTRUE(gates$PC4)
  core <- run_m1_core(first, rL, pc4); m <- core$m; dec <- core$dec
  k1 <- k1_labels(m, dec)
  fl <- fl_perm(m$fit, "PS_e", first$ADHD_e, nperm, seed + 1)
  bd <- boot_delta(core$formula, first, "PS_e", nboot, seed + 2)
  se_delta <- dec$se_hc3 / m$sigma
  bayes <- bayes_normal(dec$delta, se_delta)
  nct <- nct_ci_delta(dec$t, m$df, dec$se_hc3, m$sigma)
  g_raw <- hedges_g_raw(first$chi_SNVTA, first$PS_e)
  # K2: R2* concordance
  d2 <- set_race3(first[!is.na(first$r2s_SNVTA), ])
  m2 <- fit_ols_hc3(m1_formula("r2s_SNVTA", d2), d2); t2 <- term_test(m2, "PS_e")
  iron <- isTRUE(t2$est < 0 && t2$p_lower < REG$alpha)
  k2 <- list(n = m2$n, est = t2$est, se = t2$se, p_one_sided_lower = t2$p_lower, ci95 = t2$ci95,
             lower_iron_wording = iron,
             wording = if (iron) REGISTERED_SENTENCES$K2_iron else REGISTERED_SENTENCES$K2_notiron)
  # K3: reference agreement (CSF-only secondary reference; sessions without a valid CSF reference are missing)
  first$dref <- first$chi_SNVTA - first$chi_SNVTA_ref2
  d3 <- set_race3(first[!is.na(first$dref), ])
  m3 <- fit_ols_hc3(m1_formula("dref", d3), d3); t3 <- term_test(m3, "PS_e")
  thr <- max(REG$drift_rel * abs(dec$beta_ppb), REG$drift_abs * m$sigma)
  trigger <- isTRUE(abs(t3$est) > thr && t3$p_two < REG$drift_p)
  core_csf <- run_m1_core(d3, rL, pc4, y = "chi_SNVTA_ref2"); dec_csf <- core_csf$dec
  finalA <- if (trigger) dec$A && dec_csf$A else dec$A
  finalB <- if (trigger) dec$B && dec_csf$B else dec$B
  final <- if (finalA && finalB) "A+B" else if (finalA) "A" else if (finalB) "B" else "C"
  k3 <- list(beta_dref = t3$est, se = t3$se, p_two = t3$p_two, threshold = thr, p_threshold = REG$drift_p,
             trigger = trigger, n_valid_csf = nrow(d3), n_no_csf = nrow(first) - nrow(d3),
             no_csf_by_cell = cell_counts(first[is.na(first$dref), ]),
             ref2_type_counts = as.list(table(first$ref2_type)), csf_outcome = dec_csf$outcome,
             csf_beta = dec_csf$beta_ppb, csf_p = dec_csf$p_two, csf_pB = dec_csf$B_registered$pB,
             rule = if (trigger) "conjunctive: A and B each require agreement under both references" else "not triggered")
  # Precision (PC4 / r_L floor)
  precision <- list(PC4_pass = pc4, r_point = gates$r_point_primary, rL = rL, B_evaluable = dec$B_evaluable,
                    label = if (dec$B_evaluable) "none" else "low-precision measure",
                    delta_attenuation_corrected = if (is.numeric(gates$r_point_primary) && gates$r_point_primary > 0)
                      dec$delta / sqrt(gates$r_point_primary) else NA)
  # B co-reports at other reliability bounds (never decisive)
  tt <- term_test(m, "PS_e")
  cross <- if (is.numeric(gates$rL_cross) && is.finite(gates$rL_cross) && gates$rL_cross > 0)
    c(list(rL_cross = gates$rL_cross), bound_test(tt, m$df, dec$sigma_B, REG$deltaB_base * sqrt(min(1, gates$rL_cross)))) else
      list(rL_cross = gates$rL_cross, note = "cross-session bound not available or <= 0")
  rsr <- resid_split_reliability(first)
  resid_b <- if (!is.null(rsr) && is.finite(rsr$rL) && rsr$rL > 0)
    c(rsr, bound_test(tt, m$df, dec$sigma_B, REG$deltaB_base * sqrt(min(1, rsr$rL)))) else
      list(note = "split-echo columns absent or reliability <= 0", detail = rsr)
  # Robustness refits use the M1 (whole-brain reference) A/B/C outcome, before K3 (sec. 4.2 K4)
  dom <- dominance(m$fit, first$sess_key)
  infl <- list(diagnostics = dom, refit = NULL, influence_dependent = FALSE)
  if (dom$triggered) {
    dd <- first[!first$sess_key %in% dom$drop_ids, ]
    rc <- run_m1_core(dd, rL, pc4)
    infl$refit <- list(n = rc$m$n, outcome = rc$dec$outcome, beta = rc$dec$beta_ppb, p = rc$dec$p_two, pB = rc$dec$B_registered$pB)
    infl$influence_dependent <- rc$dec$outcome != dec$outcome
  }
  loo <- vapply(seq_len(nrow(first)), function(i) {
    r <- try_fit(run_m1_core(first[-i, ], rL, pc4), paste("LOO", i)); if (is.null(r)) NA_character_ else r$dec$outcome
  }, character(1))
  frag <- list(n_refits = length(loo), n_category_change = sum(loo != dec$outcome, na.rm = TRUE),
               changed_to = as.list(table(loo[loo != dec$outcome])), fragile = any(loo != dec$outcome, na.rm = TRUE),
               compared = "M1 A/B/C outcome (whole-brain reference, before K3)")
  s1o <- run_m1_core(S$ses1_only, rL, pc4)
  expo <- run_m1_core(S$exposure, rL, pc4)
  list(n = m$n, cells = cell_counts(first), formula = gsub("\\s+", " ", paste(deparse(core$formula), collapse = "")),
       coefficients = coef_table(m), decision = dec, sentence = outcome_sentence(dec$outcome, dec$disconfirming),
       huber = core$huber, K1 = k1, freedman_lane = fl, boot_delta = bd, se_delta = se_delta, bayes = bayes,
       nct_ci_delta = nct, hedges_g_raw = g_raw,
       ppb = list(beta = dec$beta_ppb, ci95 = dec$ci95_ppb, ci90 = dec$ci90_ppb, bound_ppb = dec$B_registered$bound_ppb,
                  adult_reference = REG$adult_ppb),
       K2 = k2, K3 = k3, final_outcome = final, final_sentence = outcome_sentence(final, dec$disconfirming),
       precision = precision, B_coreport_cross_session = cross, B_coreport_residualised_split = resid_b,
       influence = infl, loo = frag,
       ses1_only = list(n = s1o$m$n, cells = cell_counts(s1o$data), outcome = s1o$dec$outcome, beta = s1o$dec$beta_ppb,
                        delta = s1o$dec$delta, p = s1o$dec$p_two, pB = s1o$dec$B_registered$pB),
       exposure_sensitivity = list(dropped = intersect(REG$exposure_sessions, first$sess_key),
                                   n_dropped = length(intersect(REG$exposure_sessions, first$sess_key)), n = expo$m$n,
                                   outcome = expo$dec$outcome, beta = expo$dec$beta_ppb, p = expo$dec$p_two,
                                   pB = expo$dec$B_registered$pB, changes_conclusion = expo$dec$outcome != dec$outcome))
}

# QC-by-group at the SUBJECT level (Q13): (1) the subject contributes no first-usable session (not in H1);
# (2) the subject's earliest MEGRE session is not usable. Fisher's exact test by cell. If (1) gives p < .05, an
# inverse-probability-weighted M1 is co-reported, weights 1 / P(in H1 | PS_e * ADHD_e, age at the earliest MEGRE
# session, sex_e) from a subject-level logistic model.
qc_by_group <- function(S, rL, pc4 = TRUE) {
  a <- S$all[order(S$all$participant_id, S$all$ses_order), ]
  subj <- a[!duplicated(a$participant_id), c("participant_id", "cell", "PS_e", "ADHD_e", "sex_e", "age", "usable")]
  names(subj)[names(subj) == "age"] <- "age_first_megre"; names(subj)[names(subj) == "usable"] <- "earliest_usable"
  subj$in_H1 <- as.integer(subj$participant_id %in% S$first$participant_id)
  tab1 <- table(factor(subj$cell, levels = CELLS), factor(ifelse(subj$in_H1 == 1, "in H1", "not in H1"), levels = c("in H1", "not in H1")))
  tab2 <- table(factor(subj$cell, levels = CELLS), factor(ifelse(subj$earliest_usable == 1, "earliest usable", "earliest failed"),
                                                         levels = c("earliest usable", "earliest failed")))
  f1 <- stats::fisher.test(tab1); f2 <- stats::fisher.test(tab2)
  out <- list(unit = "subject", table_in_H1 = as.data.frame.matrix(tab1), fisher_p_in_H1 = f1$p.value,
              table_earliest = as.data.frame.matrix(tab2), fisher_p_earliest = f2$p.value, ipw = NULL)
  if (f1$p.value < REG$alpha) {
    g <- stats::glm(in_H1 ~ PS_e * ADHD_e + age_first_megre + sex_e, family = stats::binomial, data = subj)
    first <- set_race3(S$first)
    first$age_first_megre <- subj$age_first_megre[match(first$participant_id, subj$participant_id)]
    w <- 1 / stats::predict(g, newdata = first, type = "response")
    m <- fit_ols_hc3(m1_formula("chi_SNVTA", first), first, weights = w); tt <- term_test(m, "PS_e")
    out$ipw <- list(est = tt$est, se = tt$se, p_two = tt$p_two,
                    weight_model = "logit P(in H1) ~ PS_e*ADHD_e + age at earliest MEGRE session + sex_e (subject level)")
  }
  out
}

# PC2x (age -> putamen and GP chi, one-sided; reported after unlock, used only if PC2w is INDETERMINATE).
pc2x <- function(first) {
  res <- lapply(intersect(c("chi_PUT", "chi_GP"), names(first)), function(v) {
    dd <- first[!is.na(first[[v]]), ]; f <- stats::lm(stats::as.formula(paste(v, "~ age")), data = dd)
    tt <- summary(f)$coefficients["age", ]
    list(roi = v, n = nrow(dd), slope_ppb_per_y = unname(tt[1]), t = unname(tt[3]),
         p_one_sided = stats::pt(tt[3], stats::df.residual(f), lower.tail = FALSE))
  })
  names(res) <- vapply(res, `[[`, "", "roi")
  pass <- length(res) == 2 && all(vapply(res, function(r) r$p_one_sided < REG$alpha, logical(1)))
  list(tests = res, pass = if (length(res) == 2) pass else NA)
}
