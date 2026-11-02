#!/usr/bin/env Rscript
# Monte Carlo operating characteristics of the registered analysis FUNCTIONS (the same code run_all.R calls) on the
# real Penn LEAD first-usable design (pre-MEGRE-QC: n = 112, TD 30 / ADHD-only 21 / PS-only 24 / PS+ADHD 37).
# This is the AUTHORITATIVE source of the power figures in Q21 of the registration. Outcomes are simulated;
# roi_vol and Q_phase are fresh standard-normal columns in every draw (their values do not exist yet).
#   Rscript sim_typeI_power.R [--out sim] [--n-m1 5000] [--n-fl 500] [--n-sec 3000] [--n-t1 1000] [--n-q3 2000]
t_start <- Sys.time()
Sys.setenv(LANGUAGE = "en"); invisible(Sys.setlocale("LC_MESSAGES", "C"))
all_args <- commandArgs(trailingOnly = FALSE)
CODE_DIR <- dirname(normalizePath(sub("^--file=", "", all_args[grep("^--file=", all_args)][1])))
suppressPackageStartupMessages({ library(splines); library(sandwich); library(lme4); library(nlme); library(MASS); library(digest) })
for (f in c("constants.R", "utils.R", "design.R", "plausibility.R", "m1.R", "secondary.R", "t1.R", "synthetic.R"))
  source(file.path(CODE_DIR, "R", f))
opt <- parse_args(commandArgs(trailingOnly = TRUE))
OUT <- if (is.null(opt$out)) file.path(CODE_DIR, "sim") else opt$out; dir.create(OUT, showWarnings = FALSE)
num <- function(k, d) if (is.null(opt[[k]])) d else as.integer(opt[[k]])
N_M1 <- num("n_m1", 5000); N_FL <- num("n_fl", 500); N_SEC <- num("n_sec", 3000); N_T1 <- num("n_t1", 1000); N_Q3 <- num("n_q3", 2000)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
log_msg <- function(...) invisible(NULL)   # silence per-draw warnings inside the simulation loops

des <- build_real_design()
base <- simulate_table(des, "null", REG$seed, qc_loss = 0, inject = FALSE)
base$s3_eligible <- as.integer(des$dwi_ok & base$usable == 1)   # pre-QC S3 sample (no simulated node-rule loss)
base <- code_vars(base); S0 <- build_samples(base)
F0 <- S0$first; n <- nrow(F0)
cat("design:", cell_str(F0), "\n")
rows <- list(); addrow <- function(block, scenario, quantity, est, nsim, ref = "") {
  rows[[length(rows) + 1]] <<- data.frame(block = block, scenario = scenario, quantity = quantity, estimate = est,
                                          mc_se = sqrt(est * (1 - est) / nsim), n_sim = nsim, reference = ref)
}
draw_cov <- function(d) { d$roi_vol <- stats::rnorm(nrow(d)); d$Q_phase <- stats::rnorm(nrow(d)); d }
decide <- function(d, y = "y") {
  m <- fit_ols_hc3(m1_formula(y, d), d); hb <- huber_fit(m1_formula(y, d), d, "PS_e", stats::residuals(m$fit))
  dec <- h1_decide(m, hb$scale, 0.70); list(m = m, dec = dec)
}

# ---- M1: decision rule and K1 labels (HC3; sigma_B = min(sigma, Huber scale); TOST at +/-0.5 sigma_B) ----
m1_sim <- function(dPSo, dPSc, dAD, nsim, seed, err = c("normal", "skew_het"), dominant = 0) {
  err <- match.arg(err); set.seed(seed); d <- set_race3(F0)
  mu <- ifelse(d$cell == "PS-only", dPSo, ifelse(d$cell == "PS+ADHD", dAD + dPSc, ifelse(d$cell == "ADHD-only", dAD, 0)))
  acc <- c(A = 0, B = 0, C = 0, AB = 0, rej2 = 0, ps_effect = 0, ps_specific = 0, comorbid = 0, tost = 0, disconf = 0)
  for (s in seq_len(nsim)) {
    d <- draw_cov(d)
    e <- if (err == "normal") stats::rnorm(n) else (stats::rchisq(n, 4) - 4) / sqrt(8) * ifelse(d$PS_e > 0, 1.3, 1.0)
    d$y <- mu + e
    if (dominant > 0) { j <- sample.int(n, dominant); d$y[j] <- d$y[j] + 6 * sample(c(-1, 1), dominant, replace = TRUE) }
    r <- decide(d); dec <- r$dec; q <- k1_labels(r$m, dec)
    acc <- acc + c(dec$outcome == "A", dec$outcome == "B", dec$outcome == "C", dec$outcome == "A+B", dec$p_two < .05,
                   q$ps_effect, q$ps_specific, q$comorbid, q$tost$success, dec$disconfirming)
  }
  acc / nsim
}
cat("M1 simulations ...\n")
sc <- list(list("delta 0 (null)", 0, 0, 0, "normal", 0, "P(B) E2 v1.2: 0.740"),
           list("delta 0, skewed + heteroscedastic errors", 0, 0, 0, "skew_het", 0, "HC3 type I (E2 0.039-0.042)"),
           list("delta -0.3", -0.3, -0.3, 0, "normal", 0, ""),
           list("delta -0.4", -0.4, -0.4, 0, "normal", 0, ""),
           list("delta -0.502 (= outcome-B bound, r_L 0.70)", -0.502, -0.502, 0, "normal", 0, "P(B) must be <= .05"),
           list("delta -0.6", -0.6, -0.6, 0, "normal", 0, ""),
           list("comorbid only (PS-only 0, PS+ADHD -0.8)", 0, -0.8, 0, "normal", 0, ""),
           list("delta 0, two +/-6 SD dominant scans", 0, 0, 0, "normal", 2, "robustness of sigma_B"),
           list("delta -0.502, two +/-6 SD dominant scans", -0.502, -0.502, 0, "normal", 2, "false B must stay near .05"))
for (k in seq_along(sc)) {
  z <- sc[[k]]; ns <- if (z[[5]] == "skew_het" || z[[6]] > 0) min(N_M1, 3000) else N_M1
  r <- m1_sim(z[[2]], z[[3]], z[[4]], ns, REG$seed + 100 + k, z[[5]], z[[6]])
  addrow("M1", z[[1]], "P(two-sided p < .05)", r[["rej2"]], ns, z[[7]])
  addrow("M1", z[[1]], "P(outcome A or A+B)", r[["A"]] + r[["AB"]], ns)
  addrow("M1", z[[1]], "P(outcome B or A+B)", r[["B"]] + r[["AB"]], ns)
  addrow("M1", z[[1]], "P(outcome C)", r[["C"]], ns)
  addrow("M1", z[[1]], "P(disconfirming)", r[["disconf"]], ns)
  addrow("K1", z[[1]], "P('PS effect')", r[["ps_effect"]], ns)
  addrow("K1", z[[1]], "P('PS-specific')", r[["ps_specific"]], ns)
  addrow("K1", z[[1]], "P('comorbid-concentrated')", r[["comorbid"]], ns)
  addrow("K1", z[[1]], "P(ADHD TOST +/-0.5 sigma_B success)", r[["tost"]], ns)
  cat(" ", z[[1]], "done\n")
}

# ---- K3 reference-drift trigger (inferential rule) with zero true drift: cost to the primary decision ----
cat("K3 simulations ...\n")
for (dl in c(-0.6, 0)) for (sdr in c(0.2, 0.5)) {
  set.seed(REG$seed + 700 + 10 * abs(dl) + 100 * sdr); d <- set_race3(F0); trig <- A0 <- A1 <- B0 <- B1 <- 0
  for (s in seq_len(N_Q3)) {
    d <- draw_cov(d); d$y <- dl * d$PS_e + stats::rnorm(n); d$dref <- stats::rnorm(n, 0, sdr); d$y2 <- d$y - d$dref
    r <- decide(d); m3 <- fit_ols_hc3(m1_formula("dref", d), d); t3 <- term_test(m3, "PS_e")
    tr <- abs(t3$est) > max(REG$drift_rel * abs(r$dec$beta_ppb), REG$drift_abs * r$m$sigma) && t3$p_two < REG$drift_p
    r2 <- if (tr) decide(d, "y2") else NULL
    fA <- if (tr) r$dec$A && r2$dec$A else r$dec$A; fB <- if (tr) r$dec$B && r2$dec$B else r$dec$B
    trig <- trig + tr; A0 <- A0 + r$dec$A; A1 <- A1 + fA; B0 <- B0 + r$dec$B; B1 <- B1 + fB
  }
  scn <- sprintf("delta %g, no true drift, SD(dref) = %.1f residual SD", dl, sdr)
  addrow("K3", scn, "P(trigger)", trig / N_Q3, N_Q3)
  addrow("K3", scn, if (dl < 0) "P(A) before -> after K3" else "P(B) before -> after K3",
         if (dl < 0) A1 / N_Q3 else B1 / N_Q3, N_Q3, sprintf("before K3: %.3f", if (dl < 0) A0 / N_Q3 else B0 / N_Q3))
}

# ---- Freedman-Lane (999 permutations within ADHD strata) ----
cat("Freedman-Lane simulations ...\n")
for (dl in c(0, -0.6)) {
  set.seed(REG$seed + 200 + 10 * abs(dl)); d <- set_race3(F0); rej <- 0
  for (s in seq_len(N_FL)) {
    d <- draw_cov(d); d$y <- dl * d$PS_e + stats::rnorm(n)
    fit <- stats::lm(m1_formula("y", d), data = d); fl <- fl_perm(fit, "PS_e", d$ADHD_e, 999, REG$seed + 1000 + s)
    rej <- rej + (fl$p < .05)
  }
  addrow("Freedman-Lane", sprintf("delta %g", dl), "P(FL p < .05), 999 permutations", rej / N_FL, N_FL, if (dl == 0) "should be ~ .05" else "")
}

# ---- S1: L under the null and at L = -0.9 (SN PS -0.6, STR PS -0.1, STR ADHD -0.4; rho 0.5) ----
cat("S1 / S2 / S3 simulations ...\n")
for (alt in c(FALSE, TRUE)) {
  set.seed(REG$seed + 300 + alt); d <- set_race3(F0); rej <- rej05 <- lab <- 0; ns1 <- min(N_SEC, 1000)
  for (s in seq_len(ns1)) {
    d <- draw_cov(d); d$str_vol <- stats::rnorm(n)
    z1 <- stats::rnorm(n); z2 <- 0.5 * z1 + sqrt(0.75) * stats::rnorm(n)
    d$chi_SNVTA <- z1 + if (alt) -0.6 * d$PS_e else 0
    d$chi_STR <- z2 + if (alt) (-0.1 * d$PS_e - 0.4 * d$ADHD_e) else 0
    r <- s1_fit(d, list(mean_SNVTA = mean(d$chi_SNVTA), sd_SNVTA = stats::sd(d$chi_SNVTA), mean_STR = mean(d$chi_STR), sd_STR = stats::sd(d$chi_STR)))
    rej <- rej + (r$L$p_two < .0167); rej05 <- rej05 + (r$L$p_two < .05)
    lab <- lab + (r$L$p_two < .0167 && r$L$est < 0 && r$ps$p_lower < .05 && r$adhd$p_upper < .05)
  }
  scn <- if (alt) "L = -0.9 (rho 0.5)" else "null (rho 0.5)"
  addrow("S1", scn, "P(p_L < .0167)", rej / ns1, ns1); addrow("S1", scn, "P(p_L < .05)", rej05 / ns1, ns1)
  addrow("S1", scn, "P('dissociation' label at .0167)", lab / ns1, ns1)
}
# S2 (real S2 sample, real PRIME/BDI/ARI): null and partial r 0.25 / 0.30 (PRIME residualised on the covariates)
for (pr in c(0, 0.25, 0.30)) {
  set.seed(REG$seed + 400 + 100 * pr); d <- set_race3(S0$s2); ns2 <- nrow(d); r1 <- r2 <- r3 <- 0
  for (s in seq_len(N_SEC)) {
    d <- draw_cov(d)
    xr <- stats::residuals(stats::lm(stats::as.formula(paste("scale(PRIME_total) ~", sub("scale(PRIME_total) + ", "", s2_rhs(d), fixed = TRUE))), data = d))
    xr <- xr / stats::sd(xr)
    d$chi_SNVTA <- pr / sqrt(1 - pr^2) * xr + stats::rnorm(ns2)
    p <- term_test(fit_ols_hc3(stats::as.formula(paste("chi_SNVTA ~", s2_rhs(d))), d), "scale(PRIME_total)")$p_two
    r1 <- r1 + (p < .0167); r2 <- r2 + (p < .025); r3 <- r3 + (p < .05)
  }
  scn <- sprintf("partial r = %.2f (n = %d, real PRIME/BDI/ARI)", pr, ns2)
  addrow("S2", scn, "P(p < .0167)", r1 / N_SEC, N_SEC); addrow("S2", scn, "P(p < .025)", r2 / N_SEC, N_SEC)
  addrow("S2", scn, "P(p < .05)", r3 / N_SEC, N_SEC)
}
# S3 (real S3 sample): latent rho 0 / 0.4 / 0.5 with R_chi 0.70, R_DWI 0.55 (observed r = rho * sqrt(0.385)); the
# coupling is with TOTAL dentate chi, which shares R2_x with chi_DGM (a model covariate): R2_x in {0, 0.3, 0.5}
for (rho in c(0, 0.4, 0.5)) for (r2x in c(0, 0.3, 0.5)) {
  if (rho == 0 && r2x > 0) next
  ro <- rho * sqrt(0.70 * 0.55)
  set.seed(REG$seed + 500 + round(100 * rho) + round(10 * r2x)); d <- set_race3(S0$s3); n3 <- nrow(d); r1 <- r2 <- r3 <- 0
  for (s in seq_len(N_SEC)) {
    for (v in c("chi_DGM", "b0loc_pSCP_R", "abs_pitch", "tract_len", "FD", "NDC", "DN_R_vol", "Q_phase")) d[[v]] <- stats::rnorm(n3)
    d$nSL <- exp(stats::rnorm(n3)); d$chi_DN_R <- sqrt(r2x) * d$chi_DGM + sqrt(1 - r2x) * stats::rnorm(n3)
    d$twICVF_pSCP_R <- ro * d$chi_DN_R + sqrt(1 - ro^2) * stats::rnorm(n3)
    p <- term_test(s3_fit(d), "scale(chi_DN_R)")$p_two; r1 <- r1 + (p < .0167); r2 <- r2 + (p < .025); r3 <- r3 + (p < .05)
  }
  scn <- sprintf("latent rho %.1f (observed r %.3f), R2_x %.1f (n = %d)", rho, ro, r2x, n3)
  addrow("S3", scn, "P(p < .0167)", r1 / N_SEC, N_SEC); addrow("S3", scn, "P(p < .025)", r2 / N_SEC, N_SEC)
  addrow("S3", scn, "P(p < .05)", r3 / N_SEC, N_SEC)
}

# ---- T1: containment-df margin test (registered decision) at the boundary and under persistence; given A ----
cat("T1 simulations ...\n")
two0 <- S0$two_occ; cat("T1 design:", nrow(two0), "rows,", sum(two0$occ == 1), "returners\n")
t1_scen <- list(list(REG$t1_margin, 0.7, N_T1), list(0, 0.7, N_T1), list(0, 0.5, N_T1 %/% 2), list(0, 0.8, N_T1 %/% 2))
for (z in t1_scen) {
  bt <- z[[1]]; icc <- z[[2]]; nt <- z[[3]]
  set.seed(REG$seed + 600 + 10 * bt + 100 * icc); two <- two0; subs <- unique(two$sub); si <- match(two$sub, subs)
  rej <- ok <- nA <- rejA <- 0
  for (s in seq_len(nt)) {
    two$roi_vol <- stats::rnorm(nrow(two)); two$Q_phase <- stats::rnorm(nrow(two))
    two$chi_SNVTA <- -0.5 * two$PS_e + bt * two$PS_e * two$occ + sqrt(icc) * stats::rnorm(length(subs))[si] + sqrt(1 - icc) * stats::rnorm(nrow(two))
    r <- tryCatch(t1_core(two, 1, 0, REG$seed, do_pb = FALSE, profile = FALSE), error = function(e) NULL)
    if (is.null(r)) next
    ok <- ok + 1; rj <- r$margin$p_containment < .05; rej <- rej + rj
    f0 <- set_race3(two[two$occ == 0, ]); A <- decide(f0, "chi_SNVTA")$dec$A
    nA <- nA + A; rejA <- rejA + (A && rj)
  }
  scn <- sprintf("true beta(PS x occ) = %g, PS -0.5 total SD, ICC %.1f", bt, icc)
  addrow("T1", scn, "P(reject H0: beta >= 0.4), containment df (registered)", rej / ok, ok, if (bt > 0) "should be <= .05" else "")
  addrow("T1", scn, "P(reject | outcome A at occasion 0)", if (nA > 0) rejA / nA else NA, nA, sprintf("P(A) = %.3f", nA / ok))
}

res <- do.call(rbind, rows)
utils::write.csv(res, file.path(OUT, "typeI_power.csv"), row.names = FALSE)
L <- c("# LEAD-chi analysis code: Monte Carlo operating characteristics (synthetic outcomes, real design)", "",
       sprintf("- Run %s UTC, %s; runtime %.1f min; seed base %d.", format(t_start, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), R.version.string,
               as.numeric(difftime(Sys.time(), t_start, units = "mins")), REG$seed),
       sprintf("- Design: first-usable pre-QC sample %s; S2 n = %d; S3 n = %d; T1 %d rows / %d returners.", cell_str(F0), nrow(S0$s2), nrow(S0$s3),
               nrow(two0), sum(two0$occ == 1)),
       "- Errors N(0, 1) unless stated; delta in residual-SD units; outcome-B bound uses r_L = 0.70 (delta_B = 0.502) and sigma_B = min(sigma, Huber proposal-2 scale).",
       "- Every quantity is computed by the registered functions in R/ (fit_ols_hc3, huber_fit, h1_decide, k1_labels, fl_perm, s1_fit, s2_rhs, s3_fit, t1_core).",
       "- mc_se = binomial Monte Carlo standard error. This file is the authoritative source of the Q21 power figures.", "",
       "| block | scenario | quantity | estimate | mc_se | n_sim | reference |", "|---|---|---|---|---|---|---|",
       sprintf("| %s | %s | %s | %.3f | %.3f | %d | %s |", res$block, res$scenario, res$quantity, res$estimate, res$mc_se, res$n_sim, res$reference))
writeLines(L, file.path(OUT, "typeI_power.md"))
cat(paste(L, collapse = "\n"), "\n")
