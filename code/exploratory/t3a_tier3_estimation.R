#!/usr/bin/env Rscript
# LEAD-chi -- TASK T3a: registered TIER-3 ESTIMATION modules (prereg v1.0 FINAL, Q19 "E-dev", "E-z"; Q25 Tier 3;
# sec. 4.3 tier table). EVERYTHING HERE IS EXPLORATORY / TIER 3 (ESTIMATES ONLY; NO CLAIM). Nothing here re-decides or
# restates the frozen confirmatory outcome (H1 outcome B; results in mirror/results/confirmatory/results.json).
#
#   Rscript t3a_tier3_estimation.R [--quick]        (--quick: 100 draws instead of the registered counts; testing only)
#
# Modules
#   (1) E-dev   : PS x age with the between/within-person decomposition on ALL usable sessions (SN-VTA, STR, right
#                 dentate chi). Registered model (Q19):
#                   chi ~ PS_e*ADHD_e + ns(age_bar_c,3) + age_w + PS_e:age_w + ADHD_e:age_w + PS_e:age_bar_c + sex_e +
#                         race3 + scale(roi_vol) + scale(Q_phase) + (1 + age_w || sub)
#                 - singular random slope -> dropped from both models (registered);
#                 - estimates at the verified REML maximum (lmer bobyqa / Nelder-Mead, lme nlminb / optim, one lme4
#                   deviance function, grid profile), df by nlme containment (registered mixed-model rule);
#                 - parametric-bootstrap LRT of PS_e:age_w (lme4::bootMer, 2,000 draws, ML fits) as registered;
#                 - realised-design null calibration (1,000 simulations, saturation b in {0, .2, .35, .5} / y) run
#                   BEFORE the E-dev inference of each ROI, as registered; gives the "calibrated p" of the tier table;
#                 - developmental reading allowed only if PS x age_w and PS x age_bar do not differ (Wald p >= .10);
#                 - cell-specific within-person age slopes; cross-sectional first-usable linear slopes by cell (HC3)
#                   and the spline-by-cell display; a GAMM display (mgcv; never inferential);
#                 - age-confounding sensitivities: common support 10-18 y, overlap weights, age matching (0.5 y).
#   (2) E-z     : TD-referenced deviation z of SN-VTA chi; TD-only model chi ~ ns(age,3) + sex_e, leave-one-subject-out
#                 for TD; z = (chi - mu_TD)/sigma_TD; centile chart; share z <= -1 (registered; Wilson CI) and shares
#                 beyond |z| > 1.96 (task); mean z per cell with 95% CI. NOT a normative model (registered wording).
#   (3) E-int   : PS x ADHD interaction and the four adjusted cell means (M1 refit on first-usable; the frozen K1
#                 numbers are re-stated verbatim after an equality check) as a cell-means figure with HC3 95% CIs.
#   (4) E-anchor: the 6 threshold-psychosis youths (dx_psychosis = 1) as z against TD.
#
# Conventions follow the frozen code (R/constants.R, design.R, m1.R, t1.R are sourced unchanged): effect coding
# PS_e / ADHD_e / sex_e = +/-0.5, race3 contr.sum over the levels present, continuous covariates z-scored within the
# analysis sample, HC3 for OLS, BH-FDR within each exploratory block. Session-level ROI volume: the ROI's own volume
# (roi_vol for SN-VTA, str_vol for STR, DN_R_vol for the dentate), as in the registered S1/S3 models.

t_start <- Sys.time()
Sys.setenv(LANGUAGE = "en"); invisible(Sys.setlocale("LC_MESSAGES", "C"))
suppressPackageStartupMessages({
  library(splines); library(sandwich); library(lmtest); library(lme4); library(nlme); library(MASS); library(boot)
  library(jsonlite); library(digest); library(mgcv)
})
args <- commandArgs(trailingOnly = TRUE)
QUICK <- "--quick" %in% args

ROOT   <- "/path/to/Penn_LEAD/LEADchi"
CODE   <- file.path(ROOT, "mirror/scripts/analysis")
INPUT  <- file.path(ROOT, "mirror/unlocked/merged_qsmxt_primary_SNVTA.csv")
GATES  <- file.path(ROOT, "mirror/gates/gates_preunlock.json")
FROZEN <- file.path(ROOT, "mirror/results/confirmatory/results.json")
OUT    <- file.path(ROOT, "results_exploratory/tier3_estimation")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
for (f in c("constants.R", "utils.R", "design.R", "plausibility.R", "m1.R", "t1.R")) source(file.path(CODE, "R", f))

N_PB   <- if (QUICK) 100 else 2000     # registered: 2,000 bootMer draws
N_CAL  <- if (QUICK) 100 else 1000     # registered: 1,000 realised-design simulations per saturation value
B_SAT  <- c(0, 0.2, 0.35, 0.5)         # registered saturation rates (1/y)
SEED   <- REG$seed + 300               # offset from the confirmatory seed (never re-uses its streams)
LABEL  <- "EXPLORATORY / TIER 3 (estimation only; no claim; BH-FDR within block)"
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

log_msg("T3a start; input ", INPUT, if (QUICK) " [QUICK MODE: NON-REGISTERED DRAW COUNTS]" else "")
d <- utils::read.csv(INPUT, stringsAsFactors = FALSE, na.strings = c("NA", ""))
schema <- read_schema(file.path(CODE, "schema/roi_table_schema.csv")); validate_table(d, schema)
d <- code_vars(d); pl <- apply_plausibility(d); d <- pl$d
S <- build_samples(d)
gates <- jsonlite::read_json(GATES, simplifyVector = TRUE)
frozen <- jsonlite::read_json(FROZEN, simplifyVector = TRUE)
log_msg("samples: usable ", cell_str(S$usable), "; first-usable ", cell_str(S$first))
# categorical palette validated with the dataviz validator (all pairs, light surface): blue / green / magenta / amber
CELL_COL <- c(TD = "#2a78d6", `ADHD-only` = "#008300", `PS-only` = "#d55181", `PS+ADHD` = "#eda100")
CELL_PCH <- c(TD = 16, `ADHD-only` = 15, `PS-only` = 17, `PS+ADHD` = 18)
CELL_CODES <- data.frame(cell = CELLS, PS_e = c(-0.5, -0.5, 0.5, 0.5), ADHD_e = c(-0.5, 0.5, -0.5, 0.5))

# ------------------------------------------------------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------------------------------------------------------
wilson <- function(k, n, z = 1.959964) {
  if (n == 0) return(c(NA, NA)); p <- k / n; den <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / den; half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(max(0, ctr - half), min(1, ctr + half))
}
mean_ci <- function(x) { x <- x[is.finite(x)]; n <- length(x); if (n < 2) return(list(n = n, mean = mean(x), sd = NA, lo = NA, hi = NA, p = NA))
  m <- mean(x); s <- stats::sd(x); q <- stats::qt(.975, n - 1) * s / sqrt(n)
  list(n = n, mean = m, sd = s, lo = m - q, hi = m + q, p = 2 * stats::pt(-abs(m / (s / sqrt(n))), n - 1)) }
ci_t <- function(b, se, df) { q <- stats::qt(.975, df); c(b - q * se, b + q * se) }
p_t  <- function(b, se, df) 2 * stats::pt(-abs(b / se), df)

# Mixed model at the verified (RE)ML maximum. slope = TRUE -> (1 | sub) + (0 + age_w | sub)  [= (1 + age_w || sub)].
# Four optimisers scored on ONE lme4 deviance function; grid profile over theta confirms the maximum; refit from the grid
# argmin if better (registered rule). Containment df from the nlme fit. The returned fit object is an lmer fit
# restarted at the best theta (so refit()/simulate() act at the maximum).
mixed_max <- function(fixed, data, slope, REML = TRUE, profile = TRUE, label = "") {
  re_lmer <- if (slope) " + (1 | sub) + (0 + age_w | sub)" else " + (1 | sub)"
  f_lmer <- stats::update(fixed, stats::as.formula(paste(". ~ .", re_lmer)))
  re_nlme <- if (slope) list(sub = nlme::pdDiag(~ age_w)) else (~ 1 | sub)
  meth <- if (REML) "REML" else "ML"
  fits <- list(
    lmer_bobyqa = try_fit(lme4::lmer(f_lmer, data = data, REML = REML, control = lme4::lmerControl(optimizer = "bobyqa")), paste(label, "lmer bobyqa")),
    lmer_NelderMead = try_fit(lme4::lmer(f_lmer, data = data, REML = REML, control = lme4::lmerControl(optimizer = "Nelder_Mead")), paste(label, "lmer Nelder_Mead")),
    lme_nlminb = try_fit(nlme::lme(fixed, random = re_nlme, data = data, method = meth,
                                   control = nlme::lmeControl(opt = "nlminb", maxIter = 200, msMaxIter = 200, niterEM = 50)), paste(label, "lme nlminb")),
    lme_optim = try_fit(nlme::lme(fixed, random = re_nlme, data = data, method = meth,
                                  control = nlme::lmeControl(opt = "optim", maxIter = 500, msMaxIter = 500)), paste(label, "lme optim")))
  dev <- lme4::lmer(f_lmer, data = data, REML = REML, devFunOnly = TRUE)
  ref <- fits$lmer_bobyqa; if (is.null(ref)) ref <- fits$lmer_NelderMead
  th_names <- names(lme4::getME(ref, "theta"))
  int_first <- grepl("Intercept", th_names[1])
  theta_of <- function(f) {
    if (inherits(f, "merMod")) return(unname(lme4::getME(f, "theta")))
    vc <- nlme::VarCorr(f); v <- suppressWarnings(as.numeric(vc[, "Variance"])); rn <- rownames(vc)
    s_int <- sqrt(v[rn == "(Intercept)"]) / f$sigma
    if (!slope) return(s_int)
    s_sl <- sqrt(v[rn == "age_w"]) / f$sigma
    if (int_first) c(s_int, s_sl) else c(s_sl, s_int)
  }
  ok <- !vapply(fits, is.null, logical(1))
  th <- lapply(fits[ok], theta_of); dv <- vapply(th, dev, numeric(1))
  best <- names(dv)[which.min(dv + seq_along(dv) * 1e-10)]
  tab <- data.frame(fit = names(dv), theta = vapply(th, function(x) paste(formatC(x, digits = 5, format = "f"), collapse = "/"), ""),
                    deviance = dv, gap_to_best = dv - min(dv), row.names = NULL)
  prof <- NULL; theta_best <- th[[best]]
  if (profile) {
    if (!slope) { g <- as.list(sort(unique(c(REG$reml_theta_grid, unlist(th))))) }
    else { gi <- seq(0, 3, by = 0.05); gs <- seq(0, 1, by = 0.02)
           g <- if (int_first) lapply(seq_len(length(gi) * length(gs)), function(k) c(gi[(k - 1) %/% length(gs) + 1], gs[(k - 1) %% length(gs) + 1]))
                else lapply(seq_len(length(gi) * length(gs)), function(k) c(gs[(k - 1) %% length(gs) + 1], gi[(k - 1) %/% length(gs) + 1])) }
    pd <- vapply(g, dev, numeric(1)); k <- which.min(pd)
    prof <- list(grid_argmin = g[[k]], dev_grid_min = pd[k], dev_best = min(dv), best_is_global_on_grid = min(dv) <= pd[k] + 1e-8,
                 dev_at_theta0 = dev(rep(0, length(theta_best))))
    if (!prof$best_is_global_on_grid) {
      fp <- try_fit(lme4::lmer(f_lmer, data = data, REML = REML, start = list(theta = g[[k]]), control = lme4::lmerControl(optimizer = "bobyqa")), paste(label, "lmer from profile argmin"))
      if (!is.null(fp)) { d2 <- dev(unname(lme4::getME(fp, "theta")))
        tab <- rbind(tab, data.frame(fit = "lmer_profile_start", theta = paste(formatC(lme4::getME(fp, "theta"), digits = 5, format = "f"), collapse = "/"), deviance = d2, gap_to_best = d2 - min(dv)))
        if (d2 < min(dv)) { best <- "lmer_profile_start"; theta_best <- unname(lme4::getME(fp, "theta")); fits$lmer_profile_start <- fp; dv <- c(dv, lmer_profile_start = d2) }
        tab$gap_to_best <- tab$deviance - min(tab$deviance) }
      prof$refit_from_profile <- TRUE; prof$best_is_global_after_refit <- min(dv) <= pd[k] + 1e-8
    }
  }
  # lmer object restarted at the best theta (reported fit)
  fb <- try_fit(lme4::lmer(f_lmer, data = data, REML = REML, start = list(theta = theta_best), control = lme4::lmerControl(optimizer = "bobyqa")), paste(label, "lmer restart at best"))
  if (is.null(fb) || dev(unname(lme4::getME(fb, "theta"))) > min(dv) + 1e-6) fb <- if (inherits(fits[[best]], "merMod")) fits[[best]] else fits$lmer_bobyqa
  b <- lme4::fixef(fb); V <- as.matrix(stats::vcov(fb)); vc <- as.data.frame(lme4::VarCorr(fb))
  nl <- if (!is.null(fits$lme_nlminb)) fits$lme_nlminb else fits$lme_optim
  cdf <- if (!is.null(nl)) summary(nl)$tTable[, "DF"] else NULL
  df_source <- if (!is.null(nl)) "nlme containment" else "fallback: n_obs - n_subjects - p (nlme fits failed)"
  if (is.null(cdf)) cdf <- setNames(rep(nrow(data) - length(unique(data$sub)) - length(b), length(b)), names(b))
  list(best = best, table = tab, profile = prof, fit = fb, b = b, V = V, se = sqrt(diag(V)), containment_df = cdf, df_source = df_source,
       sd_subject = vc$sdcor[vc$grp == "sub" & vc$var1 == "(Intercept)"][1],
       sd_slope = if (slope) vc$sdcor[vc$grp == "sub.1" | (vc$grp == "sub" & vc$var1 == "age_w")][1] else NA,
       sd_resid = vc$sdcor[vc$grp == "Residual"], singular = lme4::isSingular(fb), boundary = any(theta_best < 1e-4),
       logLik = as.numeric(stats::logLik(fb)), REML = REML, slope = slope, formula_lmer = f_lmer, n_obs = nrow(data), n_sub = length(unique(data$sub)))
}
mm_term <- function(mm, term) { b <- unname(mm$b[term]); se <- unname(mm$se[term]); df <- unname(mm$containment_df[term])
  data.frame(term = term, est = b, se = se, df = df, t = b / se, p = p_t(b, se, df), lo95 = ci_t(b, se, df)[1], hi95 = ci_t(b, se, df)[2]) }
mm_lin <- function(mm, cvec, df, label) { cv <- setNames(rep(0, length(mm$b)), names(mm$b)); cv[names(cvec)] <- cvec
  est <- sum(cv * mm$b); se <- sqrt(drop(t(cv) %*% mm$V %*% cv))
  data.frame(term = label, est = est, se = se, df = df, t = est / se, p = p_t(est, se, df), lo95 = ci_t(est, se, df)[1], hi95 = ci_t(est, se, df)[2]) }

# ------------------------------------------------------------------------------------------------------------------
# (1) E-dev
# ------------------------------------------------------------------------------------------------------------------
ROIS <- list(SNVTA = list(y = "chi_SNVTA", vol = "roi_vol", name = "SN-VTA"),
             STR   = list(y = "chi_STR",   vol = "str_vol", name = "striatum (caudate + putamen + NAcc)"),
             DN_R  = list(y = "chi_DN_R",  vol = "DN_R_vol", name = "right dentate"))

prep_dev <- function(u, y, vol) {
  u <- u[!is.na(u[[y]]) & !is.na(u[[vol]]), ]
  u <- set_race3(u); u$y <- u[[y]]; u$vol <- u[[vol]]
  u$age_bar <- stats::ave(u$age, u$sub); u$age_w <- u$age - u$age_bar
  sub_mean <- mean(tapply(u$age_bar, u$sub, function(x) x[1]))      # centred at the subject-level mean of age_bar
  u$age_bar_c <- u$age_bar - sub_mean; u$cellf <- factor(u$cell, levels = CELLS)
  attr(u, "age_bar_centre") <- sub_mean; u
}
dev_fixed <- function(u, full = TRUE) stats::as.formula(paste0("y ~ PS_e * ADHD_e + ns(age_bar_c, 3) + age_w", if (full) " + PS_e:age_w" else "",
                                                               " + ADHD_e:age_w + PS_e:age_bar_c + sex_e", race_term(u), " + scale(vol) + scale(Q_phase)"))

# realised-design null calibration: same subjects / sessions / ages / cells / covariates; chi simulated with NO PS or
# ADHD effect and NO PS x age effect, from a common saturating-exponential age trajectory mu_b(age) with rate b (1/y)
# plus the observed variance components; the registered E-dev model is refitted and the t of PS_e:age_w (and of
# PS_e:age_bar_c) recorded. Amplitude A = linear OLS age slope of chi (first-usable, this ROI) x observed age range,
# so every b shares the same total rise and differs only in its curvature (b = 0: linear).
calibrate_dev <- function(u, mm, first_y, nsim, seed, label) {
  a0 <- min(u$age); a1 <- max(u$age)
  sl <- stats::coef(stats::lm(first_y$y ~ first_y$age))[2]; A <- unname(sl) * (a1 - a0)
  g_b <- function(age, b) if (b == 0) (age - a0) / (a1 - a0) else (1 - exp(-b * (age - a0))) / (1 - exp(-b * (a1 - a0)))
  subs <- unique(u$sub); isub <- match(u$sub, subs)
  tw_obs <- unname(mm$b["PS_e:age_w"] / mm$se["PS_e:age_w"]); tb_obs <- unname(mm$b["PS_e:age_bar_c"] / mm$se["PS_e:age_bar_c"])
  df_w <- unname(mm$containment_df["PS_e:age_w"]); df_b <- unname(mm$containment_df["PS_e:age_bar_c"])
  set.seed(seed)
  rows <- lapply(B_SAT, function(b) {
    mu <- A * g_b(u$age, b)
    ts <- vapply(seq_len(nsim), function(i) {
      ui <- stats::rnorm(length(subs), 0, mm$sd_subject); si <- if (mm$slope) stats::rnorm(length(subs), 0, mm$sd_slope) else 0
      y <- mu + ui[isub] + (if (mm$slope) si[isub] * u$age_w else 0) + stats::rnorm(nrow(u), 0, mm$sd_resid)
      f <- tryCatch(suppressWarnings(suppressMessages(lme4::refit(mm$fit, newresp = y))), error = function(e) NULL)
      if (is.null(f)) return(c(NA, NA))
      bb <- lme4::fixef(f); se <- sqrt(diag(as.matrix(stats::vcov(f))))
      c(bb["PS_e:age_w"] / se["PS_e:age_w"], bb["PS_e:age_bar_c"] / se["PS_e:age_bar_c"]) }, numeric(2))
    tw <- ts[1, is.finite(ts[1, ])]; tb <- ts[2, is.finite(ts[2, ])]
    data.frame(roi = label, b_saturation = b, amplitude_ppb = A, age_range = paste0(a0, "-", a1), n_valid = length(tw),
               typeI_PSxage_w_at05 = mean(abs(tw) > stats::qt(.975, df_w)), typeI_PSxage_bar_at05 = mean(abs(tb) > stats::qt(.975, df_b)),
               mean_t_w = mean(tw), sd_t_w = stats::sd(tw), mean_t_b = mean(tb), sd_t_b = stats::sd(tb),
               t_obs_w = tw_obs, p_calibrated_PSxage_w = (sum(abs(tw) >= abs(tw_obs)) + 1) / (length(tw) + 1),
               t_obs_b = tb_obs, p_calibrated_PSxage_bar = (sum(abs(tb) >= abs(tb_obs)) + 1) / (length(tb) + 1))
  })
  do.call(rbind, rows)
}

# parametric-bootstrap LRT of PS_e:age_w (lme4::bootMer on the ML null fit; the full ML model is refitted to each
# simulated response; LRT = 2 (logLik_full - logLik_null); p = (#>= + 1)/(N + 1)).
pb_lrt <- function(full_ml, null_ml, nsim, seed) {
  lrt_obs <- as.numeric(2 * (stats::logLik(full_ml) - stats::logLik(null_ml)))
  fun <- function(fit) tryCatch({ y <- lme4::getME(fit, "y"); ff <- suppressWarnings(suppressMessages(lme4::refit(full_ml, newresp = y)))
    as.numeric(2 * (stats::logLik(ff) - stats::logLik(fit))) }, error = function(e) NA_real_)
  bb <- suppressWarnings(suppressMessages(lme4::bootMer(null_ml, fun, nsim = nsim, seed = seed, type = "parametric", use.u = FALSE)))
  v <- bb$t[, 1]; v <- v[is.finite(v)]
  c(list(lrt_obs = lrt_obs, p_chisq1_asymptotic = stats::pchisq(lrt_obs, 1, lower.tail = FALSE), n_draws = nsim, n_valid = length(v)),
    mc_info((sum(v >= lrt_obs - 1e-10) + 1) / (length(v) + 1), length(v)))
}

run_edev <- function(key) {
  R <- ROIS[[key]]; log_msg("E-dev ", key)
  u <- prep_dev(S$usable, R$y, R$vol); first_y <- data.frame(y = S$first[[R$y]], age = S$first$age); first_y <- first_y[!is.na(first_y$y), ]
  fx_full <- dev_fixed(u, TRUE); fx_null <- dev_fixed(u, FALSE)
  # registered singularity rule (random slope checked on the ML full and null fits and the REML full fit)
  chk <- function(f, REML) { m <- try_fit(lme4::lmer(stats::update(f, . ~ . + (1 | sub) + (0 + age_w | sub)), data = u, REML = REML,
                                                     control = lme4::lmerControl(optimizer = "bobyqa")), "slope check"); if (is.null(m)) TRUE else lme4::isSingular(m) }
  sing <- c(full_ML = chk(fx_full, FALSE), null_ML = chk(fx_null, FALSE), full_REML = chk(fx_full, TRUE))
  slope <- !any(sing)
  log_msg("  random slope singular: ", paste(names(sing), sing, collapse = ", "), " -> slope ", if (slope) "kept" else "DROPPED from both models")
  # REML maximum (estimation fit) -- needed first for the variance components of the calibration DGP
  mm <- mixed_max(fx_full, u, slope, REML = TRUE, profile = TRUE, label = paste(key, "REML"))
  # realised-design null calibration (registered: before the E-dev inference)
  cal <- calibrate_dev(u, mm, first_y, N_CAL, SEED + 10 + match(key, names(ROIS)), key)
  # estimates
  terms <- c("PS_e:age_w", "PS_e:age_bar_c", "ADHD_e:age_w", "age_w", "PS_e", "ADHD_e", "PS_e:ADHD_e")
  est <- do.call(rbind, lapply(terms, function(t) mm_term(mm, t))); est$roi <- key
  # within vs between PS x age (Wald; conservative df = min of the two containment dfs)
  cv <- c(`PS_e:age_w` = 1, `PS_e:age_bar_c` = -1)
  wb <- mm_lin(mm, cv, min(mm$containment_df[names(cv)]), "PS_e:age_w - PS_e:age_bar_c")
  dev_reading <- isTRUE(wb$p >= 0.10)
  # ML fits and the parametric-bootstrap LRT (registered)
  ml_full <- mixed_max(fx_full, u, slope, REML = FALSE, profile = FALSE, label = paste(key, "ML full"))
  ml_null <- mixed_max(fx_null, u, slope, REML = FALSE, profile = FALSE, label = paste(key, "ML null"))
  pb <- pb_lrt(ml_full$fit, ml_null$fit, N_PB, SEED + 20 + match(key, names(ROIS)))
  # cell-specific within-person age slopes (ppb / y): age_w + PS_e * b(PS_e:age_w) + ADHD_e * b(ADHD_e:age_w)
  df_w <- unname(mm$containment_df["age_w"])
  cs <- do.call(rbind, lapply(seq_len(nrow(CELL_CODES)), function(i) {
    r <- mm_lin(mm, c(age_w = 1, `PS_e:age_w` = CELL_CODES$PS_e[i], `ADHD_e:age_w` = CELL_CODES$ADHD_e[i]), df_w, CELL_CODES$cell[i])
    r$term <- NULL; cbind(roi = key, cell = CELL_CODES$cell[i], kind = "within-person (LME, all usable sessions)", r) }))
  # PS-averaged and ADHD-averaged within-person slopes
  ps_rows <- rbind(cbind(roi = key, cell = "PS (both)", kind = "within-person (LME, all usable sessions)", { r <- mm_lin(mm, c(age_w = 1, `PS_e:age_w` = 0.5), df_w, "PS"); r$term <- NULL; r }),
                   cbind(roi = key, cell = "non-PS (both)", kind = "within-person (LME, all usable sessions)", { r <- mm_lin(mm, c(age_w = 1, `PS_e:age_w` = -0.5), df_w, "nonPS"); r$term <- NULL; r }))
  # cross-sectional first-usable linear slope by cell (HC3)
  fu <- S$first[!is.na(S$first[[R$y]]) & !is.na(S$first[[R$vol]]), ]; fu <- set_race3(fu); fu$y <- fu[[R$y]]; fu$vol <- fu[[R$vol]]; fu$cellf <- factor(fu$cell, levels = CELLS)
  fx_cs <- stats::as.formula(paste0("y ~ cellf + cellf:age + sex_e", race_term(fu), " + scale(vol) + scale(Q_phase)"))
  m_cs <- fit_ols_hc3(fx_cs, fu)
  cs2 <- do.call(rbind, lapply(CELLS, function(cl) { tt <- term_test(m_cs, paste0("cellf", cl, ":age"))
    data.frame(roi = key, cell = cl, kind = "cross-sectional linear (first-usable, HC3)", est = tt$est, se = tt$se, df = tt$df, t = tt$t, p = tt$p_two, lo95 = tt$ci95[1], hi95 = tt$ci95[2]) }))
  # spline-by-cell display (first-usable; adjusted at covariate means) and GAMM display (all usable sessions; never inferential)
  fx_sp <- stats::as.formula(paste0("y ~ cellf * ns(age, 3) + sex_e", race_term(fu), " + scale(vol) + scale(Q_phase)"))
  m_sp <- stats::lm(fx_sp, data = fu)
  spline_curves <- do.call(rbind, lapply(CELLS, function(cl) { dd <- fu[fu$cell == cl, ]; ag <- seq(min(dd$age), max(dd$age), length.out = 40)
    nd <- fu[rep(which(fu$cell == cl)[1], 40), ]; nd$age <- ag; nd$sex_e <- 0
    nd$vol <- mean(fu$vol); nd$Q_phase <- mean(fu$Q_phase)
    # race3 at its sample mean: average predictions over the race levels weighted by their frequency
    pr <- Reduce(`+`, lapply(levels(fu$race3), function(rl) { nd$race3 <- factor(rl, levels = levels(fu$race3)); stats::predict(m_sp, newdata = nd) * mean(fu$race3 == rl) }))
    data.frame(roi = key, cell = cl, age = ag, fit = pr) }))
  gam_curves <- tryCatch({
    u$subf <- factor(u$sub)
    g <- mgcv::gam(y ~ cellf + s(age, by = cellf, k = 4) + sex_e + s(subf, bs = "re"), data = u, method = "REML")
    do.call(rbind, lapply(CELLS, function(cl) { dd <- u[u$cell == cl, ]; ag <- seq(min(dd$age), max(dd$age), length.out = 40)
      nd <- data.frame(cellf = factor(cl, levels = CELLS), age = ag, sex_e = 0, subf = factor(u$sub[1], levels = levels(u$subf)))
      pr <- stats::predict(g, newdata = nd, exclude = "s(subf)", newdata.guaranteed = TRUE)
      data.frame(roi = key, cell = cl, age = ag, fit = as.numeric(pr)) })) }, error = function(e) { log_msg("GAMM display failed: ", conditionMessage(e)); NULL })
  # age-confounding sensitivities (registered for E-dev): common support 10-18 y; overlap weights; age matching 0.5 y
  sens <- list()
  sens_row <- function(name, mmx, n_obs, n_sub, note = "", est_fun = NULL) {
    if (is.null(mmx)) return(data.frame(roi = key, sensitivity = name, n_sessions = n_obs, n_subjects = n_sub, est_PSxage_w = NA, se_w = NA, lo95_w = NA, hi95_w = NA, p_w = NA,
                                         est_PSxage_bar = NA, se_b = NA, lo95_b = NA, hi95_b = NA, p_b = NA, note = paste("fit failed", note)))
    w <- mm_term(mmx, "PS_e:age_w"); b <- mm_term(mmx, "PS_e:age_bar_c")
    data.frame(roi = key, sensitivity = name, n_sessions = n_obs, n_subjects = n_sub, est_PSxage_w = w$est, se_w = w$se, lo95_w = w$lo95, hi95_w = w$hi95, p_w = w$p,
               est_PSxage_bar = b$est, se_b = b$se, lo95_b = b$lo95, hi95_b = b$hi95, p_b = b$p, note = note) }
  sens[[1]] <- sens_row("registered model (all usable sessions)", mm, nrow(u), length(unique(u$sub)), if (slope) "random intercept + uncorrelated random slope" else "random intercept only (slope singular)")
  ucs <- prep_dev(S$usable[S$usable$age >= 10 & S$usable$age <= 18, ], R$y, R$vol)
  mcs <- try_fit(mixed_max(dev_fixed(ucs, TRUE), ucs, slope, REML = TRUE, profile = FALSE, label = paste(key, "common support")), "common support")
  sens[[2]] <- sens_row("common support 10-18 y (sessions)", mcs, nrow(ucs), length(unique(ucs$sub)), "same model; sessions with 10 <= age <= 18")
  # overlap weights (Li, Morgan & Zaslavsky 2018): subject-level propensity of PS from first-usable age (ns 3) and sex;
  # w = 1 - e for PS, e for non-PS; weighted OLS of the fixed part with cluster-robust (HC2) SE (lmer weights are
  # precision weights, not sampling weights, so the weighted analysis is the OLS working-independence fit).
  ow <- tryCatch({
    sl <- fu[, c("sub", "age", "sex_e", "PS_e")]; sl$ps <- as.integer(sl$PS_e > 0)
    pm <- stats::glm(ps ~ ns(age, 3) + sex_e, family = stats::binomial, data = sl); e <- stats::fitted(pm)
    sl$w <- ifelse(sl$ps == 1, 1 - e, e); u$w <- sl$w[match(u$sub, sl$sub)]; uw <- u[!is.na(u$w), ]
    ols <- stats::lm(fx_full, data = uw, weights = w); Vc <- sandwich::vcovCL(ols, cluster = uw$sub, type = "HC2")
    bw <- stats::coef(ols); sew <- sqrt(diag(Vc)); dfc <- length(unique(uw$sub)) - 1
    data.frame(roi = key, sensitivity = "overlap weights (age ns3 + sex propensity; weighted OLS, cluster-robust HC2)", n_sessions = nrow(uw), n_subjects = length(unique(uw$sub)),
               est_PSxage_w = bw["PS_e:age_w"], se_w = sew["PS_e:age_w"], lo95_w = ci_t(bw["PS_e:age_w"], sew["PS_e:age_w"], dfc)[1], hi95_w = ci_t(bw["PS_e:age_w"], sew["PS_e:age_w"], dfc)[2], p_w = p_t(bw["PS_e:age_w"], sew["PS_e:age_w"], dfc),
               est_PSxage_bar = bw["PS_e:age_bar_c"], se_b = sew["PS_e:age_bar_c"], lo95_b = ci_t(bw["PS_e:age_bar_c"], sew["PS_e:age_bar_c"], dfc)[1], hi95_b = ci_t(bw["PS_e:age_bar_c"], sew["PS_e:age_bar_c"], dfc)[2], p_b = p_t(bw["PS_e:age_bar_c"], sew["PS_e:age_bar_c"], dfc),
               note = sprintf("effective sample size %.1f subjects; t df = clusters - 1", sum(sl$w)^2 / sum(sl$w^2))) }, error = function(e) { log_msg("overlap weights failed: ", conditionMessage(e)); NULL })
  if (!is.null(ow)) sens[[3]] <- ow
  # age matching (caliper 0.5 y): 1:1 greedy nearest-neighbour PS <-> non-PS on first-usable age, without replacement
  set.seed(SEED + 30 + match(key, names(ROIS)))
  ps_s <- fu[fu$PS_e > 0, c("sub", "age")]; np_s <- fu[fu$PS_e < 0, c("sub", "age")]; ps_s <- ps_s[sample(nrow(ps_s)), ]
  pairs <- list(); avail <- np_s
  for (i in seq_len(nrow(ps_s))) { if (!nrow(avail)) break; dd <- abs(avail$age - ps_s$age[i]); j <- which.min(dd)
    if (dd[j] <= 0.5) { pairs[[length(pairs) + 1]] <- c(ps_s$sub[i], avail$sub[j]); avail <- avail[-j, ] } }
  msub <- unlist(pairs); um <- prep_dev(S$usable[S$usable$sub %in% msub, ], R$y, R$vol)
  mmt <- if (length(pairs) >= 10) try_fit(mixed_max(dev_fixed(um, TRUE), um, slope, REML = TRUE, profile = FALSE, label = paste(key, "matched")), "age matched") else NULL
  sens[[length(sens) + 1]] <- sens_row(sprintf("age-matched subjects (caliper 0.5 y; %d PS/non-PS pairs)", length(pairs)), mmt, nrow(um), length(unique(um$sub)), "greedy 1:1 nearest neighbour on first-usable age, without replacement")
  sens <- do.call(rbind, sens)
  list(key = key, name = R$name, y = R$y, vol = R$vol, n_sessions = nrow(u), n_subjects = length(unique(u$sub)), cells_sessions = cell_counts(u),
       n_returners = sum(table(u$sub) > 1), age_bar_centre = attr(u, "age_bar_centre"), slope_singular = as.list(sing), slope_kept = slope,
       reml = list(best = mm$best, table = mm$table, profile = mm$profile, sd_subject = mm$sd_subject, sd_slope = mm$sd_slope, sd_resid = mm$sd_resid,
                   icc_intercept = mm$sd_subject^2 / (mm$sd_subject^2 + mm$sd_resid^2), boundary = mm$boundary, df_source = mm$df_source),
       formula = paste(deparse(mm$formula_lmer, width.cutoff = 500), collapse = ""), estimates = est, within_vs_between = wb, developmental_reading_allowed = dev_reading,
       ml = list(best_full = ml_full$best, best_null = ml_null$best, logLik_full = ml_full$logLik, logLik_null = ml_null$logLik, est_PSxage_w_ML = unname(ml_full$b["PS_e:age_w"])),
       pb_lrt = pb, calibration = cal, cell_slopes = rbind(cs, ps_rows, cs2), cross_sectional_formula = paste(deparse(fx_cs, width.cutoff = 500), collapse = ""),
       spline_curves = spline_curves, gam_curves = gam_curves, sensitivities = sens, data = u, first = fu)
}
EDEV <- lapply(names(ROIS), run_edev); names(EDEV) <- names(ROIS)

# BH-FDR within the E-dev block (PS x age_w, PS x age_bar, ADHD x age_w; 3 ROIs = 9 tests); PB-LRT p (3 tests)
edev_est <- do.call(rbind, lapply(EDEV, function(e) e$estimates))
blk <- edev_est$term %in% c("PS_e:age_w", "PS_e:age_bar_c", "ADHD_e:age_w")
edev_est$q_BH_block <- NA; edev_est$q_BH_block[blk] <- stats::p.adjust(edev_est$p[blk], "BH")
edev_est$block <- ifelse(blk, "E-dev age-interaction block (9 tests)", "not in the FDR block (main effects; descriptive)")
edev_est <- edev_est[, c("roi", "term", "est", "se", "df", "t", "p", "q_BH_block", "lo95", "hi95", "block")]
pb_tab <- do.call(rbind, lapply(EDEV, function(e) data.frame(roi = e$key, slope_kept = e$slope_kept, lrt_obs = e$pb_lrt$lrt_obs, p_pb = e$pb_lrt$p, mc_lo = e$pb_lrt$mc_ci95[1], mc_hi = e$pb_lrt$mc_ci95[2],
                                                               n_valid = e$pb_lrt$n_valid, p_chisq1 = e$pb_lrt$p_chisq1_asymptotic)))
pb_tab$q_BH_block <- stats::p.adjust(pb_tab$p_pb, "BH")
wb_tab <- do.call(rbind, lapply(EDEV, function(e) cbind(roi = e$key, e$within_vs_between, developmental_reading_allowed = e$developmental_reading_allowed)))
cell_sl <- do.call(rbind, lapply(EDEV, function(e) e$cell_slopes)); rownames(cell_sl) <- NULL
for (k in unique(cell_sl$kind)) { ix <- cell_sl$kind == k & cell_sl$cell %in% CELLS; cell_sl$q_BH_block[ix] <- stats::p.adjust(cell_sl$p[ix], "BH") }
cal_tab <- do.call(rbind, lapply(EDEV, function(e) e$calibration)); rownames(cal_tab) <- NULL
sens_tab <- do.call(rbind, lapply(EDEV, function(e) e$sensitivities)); rownames(sens_tab) <- NULL
utils::write.csv(edev_est, file.path(OUT, "edev_lme_estimates.csv"), row.names = FALSE)
utils::write.csv(pb_tab, file.path(OUT, "edev_pb_lrt.csv"), row.names = FALSE)
utils::write.csv(wb_tab, file.path(OUT, "edev_within_vs_between.csv"), row.names = FALSE)
utils::write.csv(cell_sl, file.path(OUT, "edev_cell_slopes.csv"), row.names = FALSE)
utils::write.csv(cal_tab, file.path(OUT, "edev_calibration.csv"), row.names = FALSE)
utils::write.csv(sens_tab, file.path(OUT, "edev_sensitivities.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, lapply(EDEV, function(e) e$spline_curves)), file.path(OUT, "edev_spline_curves_firstusable.csv"), row.names = FALSE)

# E-dev figure: per ROI, (a) all usable sessions by cell with within-person segments + GAMM display curves; (b) first-usable
# spline-by-cell fit (adjusted at covariate means) with the first-usable points.
png(file.path(OUT, "edev_figure.png"), width = 2400, height = 2700, res = 200)
par(mfrow = c(3, 2), mar = c(4.2, 4.5, 3, 1), oma = c(2, 0, 0, 0), mgp = c(2.6, 0.7, 0), las = 1, bty = "l")
for (e in EDEV) {
  u <- e$data; yl <- paste0(e$name, " chi (ppb)")
  plot(u$age, u$y, type = "n", xlab = "age at session (y)", ylab = yl, main = sprintf("(a) %s: all usable sessions (n = %d sessions / %d subjects)", e$key, e$n_sessions, e$n_subjects), cex.main = 0.95)
  grid(col = "grey88", lty = 1)
  for (s in unique(u$sub)) { dd <- u[u$sub == s, ]; if (nrow(dd) > 1) { dd <- dd[order(dd$age), ]; lines(dd$age, dd$y, col = adjustcolor(CELL_COL[dd$cell[1]], 0.55), lwd = 1.2) } }
  points(u$age, u$y, col = adjustcolor(CELL_COL[u$cell], 0.8), pch = CELL_PCH[u$cell], cex = 0.9)
  if (!is.null(e$gam_curves)) for (cl in CELLS) { g <- e$gam_curves[e$gam_curves$cell == cl, ]; lines(g$age, g$fit, col = CELL_COL[cl], lwd = 2.5) }
  legend("topleft", legend = CELLS, col = CELL_COL, pch = CELL_PCH, lwd = 2, bty = "n", cex = 0.85, title = "cell (lines: GAMM display only)")
  fu <- e$first
  plot(fu$age, fu$y, type = "n", xlab = "age at first usable session (y)", ylab = yl, main = sprintf("(b) %s: first-usable spline by cell (n = %d; adjusted display)", e$key, nrow(fu)), cex.main = 0.95)
  grid(col = "grey88", lty = 1)
  points(fu$age, fu$y, col = adjustcolor(CELL_COL[fu$cell], 0.8), pch = CELL_PCH[fu$cell], cex = 0.9)
  for (cl in CELLS) { g <- e$spline_curves[e$spline_curves$cell == cl, ]; lines(g$age, g$fit, col = CELL_COL[cl], lwd = 2.5) }
  legend("topleft", legend = CELLS, col = CELL_COL, pch = CELL_PCH, lwd = 2, bty = "n", cex = 0.85, title = "cell: ns(age, 3) x cell, HC3 model")
}
mtext(paste0("LEAD-chi E-dev display -- ", LABEL), side = 1, line = 0.6, outer = TRUE, cex = 0.7, col = "grey30")
dev.off()

# ------------------------------------------------------------------------------------------------------------------
# (2) E-z and (4) E-anchor
# ------------------------------------------------------------------------------------------------------------------
log_msg("E-z")
fu <- S$first; td <- fu[fu$cell == "TD", ]
ns_td <- ns(td$age, df = 3); K <- attr(ns_td, "knots"); Bk <- attr(ns_td, "Boundary.knots")   # knots fixed on all TD (LOSO keeps them)
ez_fit <- function(dd) stats::lm(chi_SNVTA ~ ns(age, knots = K, Boundary.knots = Bk) + sex_e, data = dd)
m_td <- ez_fit(td); sig_td <- stats::sigma(m_td)
ez_pred <- function(m, nd) { pr <- stats::predict(m, newdata = nd, se.fit = TRUE); list(mu = as.numeric(pr$fit), se = as.numeric(pr$se.fit)) }
z_of <- function(dd) {   # dd: rows to score; TD subjects scored LOSO, everyone else against the full TD model
  out <- lapply(seq_len(nrow(dd)), function(i) {
    r <- dd[i, ]
    if (r$cell == "TD" && r$sub %in% td$sub) { m <- ez_fit(td[td$sub != r$sub, ]); loso <- TRUE } else { m <- m_td; loso <- FALSE }
    p <- ez_pred(m, r); s <- stats::sigma(m)
    data.frame(mu_TD = p$mu, se_fit = p$se, sigma_TD = s, z = (r$chi_SNVTA - p$mu) / s, z_pred = (r$chi_SNVTA - p$mu) / sqrt(s^2 + p$se^2), loso = loso,
               extrapolated = r$age < Bk[1] | r$age > Bk[2]) })
  cbind(dd[, c("sub", "session_id", "cell", "age", "sex", "dx_psychosis", "chi_SNVTA")], do.call(rbind, out))
}
ez <- z_of(fu); ez$centile <- 100 * stats::pnorm(ez$z); ez$cell <- factor(ez$cell, levels = CELLS)
ez$threshold_psychosis <- ez$dx_psychosis %in% 1
ez_out <- ez[order(ez$cell, ez$age), ]; ez_out$cell <- as.character(ez_out$cell)
utils::write.csv(ez_out, file.path(OUT, "ez_subject_z.csv"), row.names = FALSE)
ez_summary <- do.call(rbind, lapply(CELLS, function(cl) {
  z <- ez$z[ez$cell == cl]; n <- length(z); mc <- mean_ci(z)
  k1 <- sum(z <= -1); k2 <- sum(z < -1.96); k3 <- sum(z > 1.96); k4 <- sum(abs(z) > 1.96)
  w1 <- wilson(k1, n); w2 <- wilson(k2, n); w3 <- wilson(k3, n); w4 <- wilson(k4, n)
  data.frame(cell = cl, n = n, mean_z = mc$mean, sd_z = mc$sd, lo95 = mc$lo, hi95 = mc$hi, p_mean_vs_0 = mc$p, median_z = stats::median(z),
             n_z_le_m1 = k1, share_z_le_m1 = k1 / n, wilson_lo_le_m1 = w1[1], wilson_hi_le_m1 = w1[2],
             n_z_lt_m196 = k2, share_z_lt_m196 = k2 / n, wilson_lo_lt_m196 = w2[1], wilson_hi_lt_m196 = w2[2],
             n_z_gt_196 = k3, share_z_gt_196 = k3 / n, wilson_lo_gt_196 = w3[1], wilson_hi_gt_196 = w3[2],
             n_absz_gt_196 = k4, share_absz_gt_196 = k4 / n, wilson_lo_absz_gt_196 = w4[1], wilson_hi_absz_gt_196 = w4[2],
             mean_z_pred = mean(ez$z_pred[ez$cell == cl]), n_extrapolated = sum(ez$extrapolated[ez$cell == cl])) }))
ez_summary$q_BH_block <- NA; ix <- ez_summary$cell != "TD"; ez_summary$q_BH_block[ix] <- stats::p.adjust(ez_summary$p_mean_vs_0[ix], "BH")
ez_summary$reference <- ifelse(ez_summary$cell == "TD", "leave-one-subject-out (TD)", "full TD model (n = 22)")
utils::write.csv(ez_summary, file.path(OUT, "ez_cell_summary.csv"), row.names = FALSE)
# also the cell contrasts of mean z vs TD (Welch t)
ez_vs_td <- do.call(rbind, lapply(CELLS[-1], function(cl) { tt <- stats::t.test(ez$z[ez$cell == cl], ez$z[ez$cell == "TD"])
  data.frame(contrast = paste(cl, "- TD (mean z; Welch)"), diff = unname(diff(rev(tt$estimate))), lo95 = tt$conf.int[1], hi95 = tt$conf.int[2], p = tt$p.value) }))
ez_vs_td$q_BH_block <- stats::p.adjust(ez_vs_td$p, "BH")
utils::write.csv(ez_vs_td, file.path(OUT, "ez_cell_vs_TD.csv"), row.names = FALSE)

# E-anchor: the 6 threshold-psychosis youths, every MEGRE session; usable sessions scored against the full TD model
anc <- d[d$dx_psychosis %in% 1, ]; anc <- anc[order(anc$sub, anc$ses_order), ]
anc_z <- z_of(anc); anc_z$centile <- 100 * stats::pnorm(anc_z$z)
anc_z <- cbind(anc_z[, c("sub", "session_id", "cell", "age", "sex", "chi_SNVTA", "mu_TD", "sigma_TD", "z", "z_pred", "centile", "extrapolated")],
               usable = anc$usable, first_usable = anc$first_usable, fail_reason = anc$fail_reason,
               status = ifelse(anc$usable == 1, ifelse(anc$first_usable %in% 1, "usable; first-usable (in H1 sample)", "usable; later session"),
                               "NOT USABLE (failed QC; value shown for completeness only; excluded from every summary)"))
utils::write.csv(anc_z, file.path(OUT, "eanchor_z.csv"), row.names = FALSE)
anc_first <- anc_z[anc_z$first_usable %in% 1, ]; anc_us <- anc_z[anc_z$usable == 1, ]
anc_sum <- list(n_threshold_psychosis_youths = length(unique(anc_z$sub)), n_with_usable_session = length(unique(anc_us$sub)),
                n_without_usable_session = length(unique(anc_z$sub)) - length(unique(anc_us$sub)), cells = as.list(table(tapply(anc_z$cell, anc_z$sub, `[`, 1))),
                first_usable = c(mean_ci(anc_first$z), list(z = anc_first$z, subjects = anc_first$sub, n_z_le_m1 = sum(anc_first$z <= -1), wilson_le_m1 = wilson(sum(anc_first$z <= -1), nrow(anc_first)),
                                                            n_absz_gt_196 = sum(abs(anc_first$z) > 1.96), mean_centile = mean(anc_first$centile))),
                all_usable_sessions = c(mean_ci(anc_us$z), list(n_sessions = nrow(anc_us), n_subjects = length(unique(anc_us$sub)), n_z_le_m1 = sum(anc_us$z <= -1), n_absz_gt_196 = sum(abs(anc_us$z) > 1.96))),
                PSADHD_without_anchors = c(mean_ci(ez$z[ez$cell == "PS+ADHD" & !ez$threshold_psychosis]), list(label = "PS+ADHD first-usable, anchors removed")),
                note = "All 6 sit in PS+ADHD (registered caveat). z against the TD-only model chi ~ ns(age, 3) + sex_e (n = 22). Estimates only; n = 3 usable youths.")

# E-z centile chart: TD-referenced centiles (sex-averaged) over the TD age range; one panel per cell; anchors marked
png(file.path(OUT, "ez_centile_chart.png"), width = 2400, height = 2000, res = 200)
par(mfrow = c(2, 2), mar = c(4.2, 4.5, 3, 1), oma = c(2, 0, 0, 0), mgp = c(2.6, 0.7, 0), las = 1, bty = "l")
ag <- seq(Bk[1], Bk[2], length.out = 100); nd <- data.frame(age = ag, sex_e = 0); mu <- ez_pred(m_td, nd)$mu
zq <- c(-1.96, -1, 0, 1, 1.96); yl <- range(c(fu$chi_SNVTA, mu + 1.96 * sig_td, mu - 1.96 * sig_td, anc_z$chi_SNVTA[anc_z$usable == 1]))
for (cl in CELLS) {
  plot(ag, mu, type = "n", ylim = yl, xlim = range(c(fu$age, ag)), xlab = "age at first usable session (y)", ylab = "SN-VTA chi (ppb)",
       main = sprintf("%s (n = %d): TD-referenced centiles", cl, sum(ez$cell == cl)), cex.main = 0.95)
  grid(col = "grey88", lty = 1)
  polygon(c(ag, rev(ag)), c(mu - 1.96 * sig_td, rev(mu + 1.96 * sig_td)), col = adjustcolor(CELL_COL["TD"], 0.08), border = NA)
  polygon(c(ag, rev(ag)), c(mu - sig_td, rev(mu + sig_td)), col = adjustcolor(CELL_COL["TD"], 0.10), border = NA)
  for (q in zq) lines(ag, mu + q * sig_td, col = adjustcolor(CELL_COL["TD"], if (q == 0) 0.9 else 0.5), lwd = if (q == 0) 2.2 else 1.2, lty = if (q == 0) 1 else 2)
  text(Bk[2], mu[100] + zq * sig_td, labels = c("2.5", "16", "50", "84", "97.5"), pos = 4, cex = 0.65, col = "grey30", xpd = NA)
  dd <- ez[ez$cell == cl, ]
  points(dd$age, dd$chi_SNVTA, col = adjustcolor(CELL_COL[cl], 0.85), pch = CELL_PCH[cl], cex = 1.1)
  if (cl == "PS+ADHD") { aa <- dd[dd$threshold_psychosis, ]; points(aa$age, aa$chi_SNVTA, pch = 24, bg = "black", col = "black", cex = 1.3)
    legend("bottomright", legend = c("PS+ADHD", "threshold-psychosis anchor (first usable)"), pch = c(CELL_PCH[cl], 24), col = c(CELL_COL[cl], "black"), pt.bg = c(NA, "black"), bty = "n", cex = 0.8) }
  if (cl == "TD") legend("bottomright", legend = c("TD (points scored leave-one-out)", "TD centile curves (full model, sex-averaged)"), pch = c(CELL_PCH[cl], NA), lty = c(NA, 2), col = CELL_COL["TD"], bty = "n", cex = 0.8)
}
mtext(paste0("LEAD-chi E-z display -- ", LABEL, " -- not a normative model"), side = 1, line = 0.6, outer = TRUE, cex = 0.7, col = "grey30")
dev.off()

# ------------------------------------------------------------------------------------------------------------------
# (3) E-int: PS x ADHD and the four adjusted cell means (M1 refit; frozen K1 numbers re-stated after an equality check)
# ------------------------------------------------------------------------------------------------------------------
log_msg("E-int")
fu1 <- set_race3(S$first); m1 <- fit_ols_hc3(m1_formula("chi_SNVTA", fu1), fu1)
chk_beta <- abs(unname(m1$b["PS_e"]) - frozen$H1$decision$beta_ppb); chk_int <- abs(unname(m1$b["PS_e:ADHD_e"]) - frozen$H1$K1$interaction$est)
chk_se <- abs(unname(m1$se["PS_e:ADHD_e"]) - frozen$H1$K1$interaction$se)
if (max(chk_beta, chk_int, chk_se) > 1e-8) stop("E-int: the M1 refit does not reproduce the frozen results.json (differences ", chk_beta, ", ", chk_int, ", ", chk_se, ")")
log_msg("  M1 refit reproduces frozen beta_PS, interaction and its HC3 SE (max abs diff ", signif(max(chk_beta, chk_int, chk_se), 3), ")")
X <- stats::model.matrix(m1$fit); xbar <- colMeans(X)
cell_means <- do.call(rbind, lapply(seq_len(nrow(CELL_CODES)), function(i) {
  x <- xbar; x["PS_e"] <- CELL_CODES$PS_e[i]; x["ADHD_e"] <- CELL_CODES$ADHD_e[i]; x["PS_e:ADHD_e"] <- CELL_CODES$PS_e[i] * CELL_CODES$ADHD_e[i]
  est <- sum(x * m1$b); se <- sqrt(drop(t(x) %*% m1$V %*% x)); q <- stats::qt(.975, m1$df)
  data.frame(cell = CELL_CODES$cell[i], n = sum(fu1$cell == CELL_CODES$cell[i]), adjusted_mean_ppb = est, se_hc3 = se, lo95 = est - q * se, hi95 = est + q * se,
             raw_mean_ppb = mean(fu1$chi_SNVTA[fu1$cell == CELL_CODES$cell[i]]), raw_sd = stats::sd(fu1$chi_SNVTA[fu1$cell == CELL_CODES$cell[i]])) }))
cell_means$note <- "M1 (first-usable, n = 80), HC3; covariates at their sample means (ns(age) basis, sex_e, race3, scaled roi_vol and Q_phase)"
utils::write.csv(cell_means, file.path(OUT, "eint_cell_means.csv"), row.names = FALSE)
K1 <- frozen$H1$K1
k1_restated <- data.frame(contrast = c("PS main effect (beta_PS; H1 decision quantity, frozen)", "ADHD main effect", "PS x ADHD interaction", "PS-only - TD (mandatory co-report)", "PS+ADHD - ADHD-only"),
                          est_ppb = c(frozen$H1$decision$beta_ppb, K1$ADHD$est, K1$interaction$est, K1$simple_PSonly_minus_TD$est, K1$simple_PSADHD_minus_ADHDonly$est),
                          se_hc3 = c(frozen$H1$decision$se_hc3, K1$ADHD$se, K1$interaction$se, K1$simple_PSonly_minus_TD$se, K1$simple_PSADHD_minus_ADHDonly$se),
                          lo95 = c(frozen$H1$decision$ci95_ppb[1], K1$ADHD$ci95[1], K1$interaction$ci95[1], K1$simple_PSonly_minus_TD$ci95[1], K1$simple_PSADHD_minus_ADHDonly$ci95[1]),
                          hi95 = c(frozen$H1$decision$ci95_ppb[2], K1$ADHD$ci95[2], K1$interaction$ci95[2], K1$simple_PSonly_minus_TD$ci95[2], K1$simple_PSADHD_minus_ADHDonly$ci95[2]),
                          p_two = c(frozen$H1$decision$p_two, K1$ADHD$p_two, K1$interaction$p_two, K1$simple_PSonly_minus_TD$p_two, K1$simple_PSADHD_minus_ADHDonly$p_two),
                          source = "frozen results.json (verbatim; not re-decided)")
utils::write.csv(k1_restated, file.path(OUT, "eint_frozen_K1_restated.csv"), row.names = FALSE)

png(file.path(OUT, "eint_cell_means.png"), width = 1800, height = 1400, res = 200)
par(mar = c(4.5, 4.8, 3.5, 1), mgp = c(2.8, 0.7, 0), las = 1, bty = "l")
xpos <- c(1, 2, 3.4, 4.4); names(xpos) <- c("TD", "PS-only", "ADHD-only", "PS+ADHD")
set.seed(SEED + 99)
plot(NA, xlim = c(0.4, 5), ylim = range(c(fu1$chi_SNVTA, cell_means$lo95, cell_means$hi95)) + c(-2, 6), xaxt = "n", xlab = "", ylab = "SN-VTA chi (ppb)",
     main = "E-int: adjusted cell means (M1, first-usable n = 80; HC3 95% CI)", cex.main = 0.95)
axis(1, at = xpos, labels = names(xpos), cex.axis = 0.9); grid(col = "grey88", lty = 1, nx = NA, ny = NULL)
mtext("without ADHD", side = 1, line = 2.6, at = mean(xpos[1:2]), cex = 0.85, col = "grey30"); mtext("with ADHD", side = 1, line = 2.6, at = mean(xpos[3:4]), cex = 0.85, col = "grey30")
for (cl in names(xpos)) { yy <- fu1$chi_SNVTA[fu1$cell == cl]; points(jitter(rep(xpos[cl], length(yy)), amount = 0.14), yy, pch = CELL_PCH[cl], col = adjustcolor(CELL_COL[cl], 0.35), cex = 0.9) }
for (i in seq_len(nrow(cell_means))) { cl <- cell_means$cell[i]; x <- xpos[cl]
  segments(x, cell_means$lo95[i], x, cell_means$hi95[i], lwd = 2.6, col = CELL_COL[cl]); points(x, cell_means$adjusted_mean_ppb[i], pch = 21, bg = CELL_COL[cl], col = "white", cex = 1.9, lwd = 1.5)
  text(x + 0.25, cell_means$adjusted_mean_ppb[i], sprintf("%.1f", cell_means$adjusted_mean_ppb[i]), cex = 0.8, adj = 0, col = "grey20") }
lines(xpos[c("TD", "PS-only")], cell_means$adjusted_mean_ppb[match(c("TD", "PS-only"), cell_means$cell)], col = "grey50", lty = 3)
lines(xpos[c("ADHD-only", "PS+ADHD")], cell_means$adjusted_mean_ppb[match(c("ADHD-only", "PS+ADHD"), cell_means$cell)], col = "grey50", lty = 3)
legend("topleft", bty = "n", cex = 0.78, text.col = "grey20",
       legend = c(sprintf("PS x ADHD interaction = %.2f ppb [%.2f, %.2f], p = %.3f (frozen K1; not re-decided)", K1$interaction$est, K1$interaction$ci95[1], K1$interaction$ci95[2], K1$interaction$p_two),
                  sprintf("PS-only - TD = %.2f [%.2f, %.2f], p = %.3f (co-report, no label)", K1$simple_PSonly_minus_TD$est, K1$simple_PSonly_minus_TD$ci95[1], K1$simple_PSonly_minus_TD$ci95[2], K1$simple_PSonly_minus_TD$p_two),
                  sprintf("PS+ADHD - ADHD-only = %.2f [%.2f, %.2f], p = %.3f", K1$simple_PSADHD_minus_ADHDonly$est, K1$simple_PSADHD_minus_ADHDonly$ci95[1], K1$simple_PSADHD_minus_ADHDonly$ci95[2], K1$simple_PSADHD_minus_ADHDonly$p_two),
                  "faint points: raw first-usable values (unadjusted)"))
mtext(LABEL, side = 1, line = 3.6, cex = 0.7, col = "grey30")
dev.off()

# ------------------------------------------------------------------------------------------------------------------
# JSON + log
# ------------------------------------------------------------------------------------------------------------------
jsonify <- function(x) {
  if (is.data.frame(x)) return(x)
  if (is.list(x)) return(lapply(x, jsonify))
  if (is.table(x)) return(as.list(x))
  if (is.atomic(x) && !is.null(names(x)) && length(x) > 1) return(as.list(x))
  if (is.atomic(x) && !is.null(names(x)) && length(x) == 1) return(unname(x))
  x
}
inst <- c(R = paste(R.version$major, R.version$minor, sep = "."), vapply(c("lme4", "nlme", "sandwich", "mgcv", "boot", "MASS"), function(p) as.character(utils::packageVersion(p)), ""))
RES <- list(
  meta = list(task = "T3a registered Tier-3 estimation modules (E-dev, E-z, E-int, E-anchor)", label = LABEL, quick_mode_nonregistered_counts = QUICK,
              input = INPUT, input_sha256 = sha256_file(INPUT), gates = GATES, gates_sha256 = sha256_file(GATES), frozen_results = FROZEN, frozen_results_sha256 = sha256_file(FROZEN),
              frozen_code_sourced = file.path(CODE, "R", c("constants.R", "utils.R", "design.R", "plausibility.R", "m1.R", "t1.R")),
              started_utc = format(t_start, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), r_version = R.version.string, packages = as.list(inst), seed = SEED,
              draws = list(n_pb_lrt = N_PB, n_calibration_per_b = N_CAL, b_saturation = B_SAT),
              confirmatory_outcome_frozen = "H1 outcome B (results.json); not re-decided or restated here"),
  samples = list(usable = cell_counts(S$usable), first_usable = cell_counts(S$first), plausibility_log = pl$log),
  E_dev = lapply(EDEV, function(e) { e$data <- NULL; e$first <- NULL; e$spline_curves <- NULL; e$gam_curves <- NULL; e }),
  E_dev_tables = list(estimates_with_BH = edev_est, pb_lrt = pb_tab, within_vs_between = wb_tab, cell_slopes = cell_sl, calibration = cal_tab, sensitivities = sens_tab),
  E_z = list(model = "chi_SNVTA ~ ns(age, 3) + sex_e on TD first-usable (n = 22); knots fixed on all TD; TD subjects leave-one-subject-out; others: full TD model",
             knots = K, boundary_knots = Bk, sigma_TD_full = sig_td, coef_full = as.list(stats::coef(m_td)), cell_summary = ez_summary, cell_vs_TD = ez_vs_td,
             registered_wording = "not a normative model"),
  E_int = list(m1_reproduces_frozen = TRUE, max_abs_diff_vs_frozen = max(chk_beta, chk_int, chk_se), cell_means = cell_means, frozen_K1_restated = k1_restated,
               k1_label_frozen = K1$label, mandatory_caveat = "all 6 threshold-psychosis youths sit in PS+ADHD"),
  E_anchor = list(sessions = anc_z, summary = anc_sum))
RES$meta$finished_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"); RES$meta$runtime_min <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))
jsonlite::write_json(jsonify(RES), file.path(OUT, "tier3_estimation.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null", null = "null")
writeLines(.LOG$lines, file.path(OUT, "t3a_log.txt"))
writeLines(capture.output(utils::sessionInfo()), file.path(OUT, "sessionInfo.txt"))
log_msg("done; runtime ", f3(RES$meta$runtime_min, 2), " min")
