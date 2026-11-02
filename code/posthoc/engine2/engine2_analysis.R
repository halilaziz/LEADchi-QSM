#!/usr/bin/env Rscript
# LEAD-chi POST HOC (2026-10-07): second, independent QSM engine (QSM.jl 0.5.4: Laplacian unwrapping, weighted
# multi-echo fit, V-SHARP, RTS [main] / TV [sensitivity]) versus the registered primary (QSMxT: ROMEO, V-SHARP, MEDI).
# Same coded inputs, Mask 3, ROI warps, plausibility rules and two-step whole-brain reference (registered scripts run
# unchanged with LX pointed to posthoc/shadow). This script:
#   (1) agreement of SN-VTA chi between engines (ICC(3,1), Pearson, Bland-Altman) in the 80 first-usable sessions;
#   (2) re-fits the frozen M1 + decision rule with engine-2 chi on the same 80 sessions; main: rows whose engine-2 value
#       is missing or whose engine-2 SN-VTA ROI fails the registered plausibility rules are missing (sec. 12.4);
#       sensitivity: plausibility failures kept;
#   (3) engine-2 recovery factor k from the -14.24 ppb injected tree (e2rec_m14_* minus e2_lap_vsharp_*).
suppressPackageStartupMessages({ library(splines); library(sandwich); library(MASS); library(jsonlite) })
ROOT <- "/path/to/Penn_LEAD/LEADchi"
OUT <- file.path(ROOT, "results_posthoc", "engine2")
RDIR <- file.path(ROOT, "mirror", "scripts", "analysis", "R")
for (f in c("constants.R", "utils.R", "design.R", "plausibility.R", "m1.R")) source(file.path(RDIR, f))
RL <- 0.938278557806515
key <- read.csv(file.path(ROOT, "mirror", "unlocked", "key_plain.csv"), stringsAsFactors = FALSE)
mer <- code_vars(read.csv(file.path(ROOT, "mirror", "unlocked", "merged_qsmxt_primary_SNVTA.csv"), stringsAsFactors = FALSE))
S <- build_samples(mer); first <- S$first
rd <- function(v) {
  x <- read.csv(file.path(OUT, "raw", sprintf("roi_qsmxt_%s.csv", v)), stringsAsFactors = FALSE)
  x$ses_order <- as.integer(x$ses_order)
  x <- merge(x, key[, c("code", "ses_order", "participant_id", "session_id")], by = c("code", "ses_order"))
  x$sess_key <- paste0(x$participant_id, "_", x$session_id); x
}
icc31 <- function(a, b) {               # two-way mixed, consistency, single measures (Shrout & Fleiss 3,1)
  Y <- cbind(a, b); Y <- Y[stats::complete.cases(Y), , drop = FALSE]; n <- nrow(Y); k <- 2
  gm <- mean(Y); rm <- rowMeans(Y); cm <- colMeans(Y)
  ssr <- k * sum((rm - gm)^2); ssc <- n * sum((cm - gm)^2); sse <- sum((Y - gm)^2) - ssr - ssc
  msr <- ssr / (n - 1); mse <- sse / ((n - 1) * (k - 1))
  list(icc = (msr - mse) / (msr + (k - 1) * mse), n = n)
}
res <- list(meta = list(run_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), label = "POST HOC engine-2 comparison"))
obs <- run_m1_core(first, RL)
res$primary <- obs$dec[c("beta_ppb", "se_hc3", "p_two", "ci95_ppb", "sigma_resid", "sigma_B", "delta", "outcome")]
bound_tests <- function(r, df) {          # per-stratum one-sided bound tests (post hoc, descriptive)
  bnd <- r$dec$B_registered$bound_ppb; tq <- stats::qt(.95, df)
  lapply(list(PSonly_minus_TD = lin_test(r$m, c(PS_e = 1, "PS_e:ADHD_e" = -0.5), "PS-only - TD"),
              PSADHD_minus_ADHD = lin_test(r$m, c(PS_e = 1, "PS_e:ADHD_e" = 0.5), "PS+ADHD - ADHD-only")),
         function(z) list(est = z$est, se = z$se, p_two = z$p_two, ci95 = z$ci95, lower95_one_sided = z$est - tq * z$se,
                          pB = stats::pt((z$est - bnd) / z$se, df, lower.tail = FALSE)))
}
adult_test <- function(beta, se, df, k) list(k = k, expected_shift = k * REG$adult_ppb,
  p_one_sided = stats::pt((beta - k * REG$adult_ppb) / se, df, lower.tail = FALSE))
res$primary$huber_beta <- obs$huber$est
res$primary$bound_ppb <- obs$dec$B_registered$bound_ppb; res$primary$pB <- obs$dec$B_registered$pB
res$primary$strata <- bound_tests(obs, obs$m$df)
kprim <- tryCatch(fromJSON(file.path(ROOT, "results_posthoc", "recovery", "results.json"))$k_first_usable[["-14.24"]][2],
                  error = function(e) NA_real_)
res$primary$adult_k_scaled <- adult_test(obs$dec$beta_ppb, obs$dec$se_hc3, obs$m$df, kprim)
fit_e2 <- function(m, keep, label) {
  d2 <- m[keep, ]; d2$chi_SNVTA <- d2$e2_chi_SNVTA
  r2 <- run_m1_core(d2, RL)
  list(label = label, n = nrow(d2),
       M1 = r2$dec[c("beta_ppb", "se_hc3", "p_two", "ci95_ppb", "sigma_resid", "sigma_B", "delta", "outcome", "disconfirming")],
       bound_ppb_registered_deltaB = r2$dec$B_registered$bound_ppb, pB = r2$dec$B_registered$pB,
       huber_beta = r2$huber$est, interaction = term_test(r2$m, "PS_e:ADHD_e")[c("est", "se", "p_two", "ci95")],
       strata = bound_tests(r2, r2$m$df), r = r2)
}
for (inv in c("rts", "tv")) {
  e2 <- rd(paste0("e2_lap_vsharp_", inv))
  cols <- c("sess_key", "status", "chi_SNVTA", "chi_SNVTA_raw", "chi_GP", "chi_STR", "chi_SNonly", "fail_SNVTA", "inmask_SNVTA")
  m <- merge(first, setNames(e2[, cols], c("sess_key", paste0("e2_", cols[-1]))), by = "sess_key", all.x = TRUE)
  finite <- is.finite(m$e2_chi_SNVTA)
  valid <- finite & !(m$e2_fail_SNVTA %in% 1)      # registered sec. 12.4: an engine reconstruction that fails the
                                                    # plausibility rules is missing for that variant only
  ag <- list(n_first_usable = nrow(m), n_finite_e2 = sum(finite), n_valid_e2 = sum(valid),
             n_fail_SNVTA_e2 = sum(finite & m$e2_fail_SNVTA %in% 1),
             status_of_missing = as.list(table(m$e2_status[!finite], useNA = "ifany")),
             icc31 = icc31(m$chi_SNVTA[finite], m$e2_chi_SNVTA[finite])$icc,
             pearson = cor(m$chi_SNVTA[finite], m$e2_chi_SNVTA[finite]),
             spearman = cor(m$chi_SNVTA[finite], m$e2_chi_SNVTA[finite], method = "spearman"),
             mean_primary = mean(m$chi_SNVTA[finite]), mean_e2 = mean(m$e2_chi_SNVTA[finite]),
             sd_primary = sd(m$chi_SNVTA[finite]), sd_e2 = sd(m$e2_chi_SNVTA[finite]),
             bland_altman_mean_diff = mean(m$e2_chi_SNVTA[finite] - m$chi_SNVTA[finite]),
             bland_altman_sd_diff = sd(m$e2_chi_SNVTA[finite] - m$chi_SNVTA[finite]),
             gp_pearson = suppressWarnings(cor(m$chi_GP, m$e2_chi_GP, use = "complete.obs")))
  main <- fit_e2(m, valid, "main: engine-2 plausibility failures missing (registered sec. 12.4 rule)")
  sens <- fit_e2(m, finite, "sensitivity: engine-2 SN-VTA plausibility failures kept")
  rec <- rd(paste0("e2rec_m14_", inv))
  kk <- merge(e2[, c("sess_key", "status", "chi_SNVTA_raw")], rec[, c("sess_key", "status", "chi_SNVTA_raw")],
              by = "sess_key", suffixes = c("", "_rec"))
  kk <- kk[kk$status == "ok" & kk$status_rec == "ok", ]
  kk$k <- (kk$chi_SNVTA_raw_rec - kk$chi_SNVTA_raw) / -14.24
  kf <- kk[kk$sess_key %in% first$sess_key, ]
  kmean <- mean(kf$k)
  for (nm in c("main", "sens")) {
    x <- get(nm); x$adult_k_scaled <- adult_test(x$M1$beta_ppb, x$M1$se_hc3, x$r$m$df, kmean); x$r <- NULL; assign(nm, x)
  }
  res[[paste0("e2_", inv)]] <- list(agreement = ag, main = main, sensitivity_failures_kept = sens,
    recovery_k_L14 = list(all = c(n = nrow(kk), mean = mean(kk$k), sd = sd(kk$k), median = median(kk$k)),
                          first_usable = c(n = nrow(kf), mean = kmean, sd = sd(kf$k), median = median(kf$k))))
  write.csv(m[, c("sess_key", "cell", "chi_SNVTA", "e2_chi_SNVTA", "e2_status", "e2_fail_SNVTA", "e2_inmask_SNVTA")],
            file.path(OUT, sprintf("agreement_first_usable_%s.csv", inv)), row.names = FALSE)
}
write_json(res, file.path(OUT, "results.json"), auto_unbox = TRUE, digits = 8, pretty = TRUE)
cat(toJSON(res, auto_unbox = TRUE, digits = 4, pretty = TRUE), "\n")
