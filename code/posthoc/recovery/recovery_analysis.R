#!/usr/bin/env Rscript
# LEAD-chi POST HOC (2026-10-07): synthetic SN-VTA susceptibility recovery test.
# A uniform change L (ppb) was injected into the SN-VTA ROI voxel set of every coded session (field added to the
# measured phase; posthoc/recovery/scripts/inject_sn.py) and the registered primary QSMxT chain + plausibility + ROI
# scripts were re-run unchanged. This script
#   (1) estimates the per-session recovery factor k = (chi_rec - chi_primary) / L of the SN-VTA mean;
#   (2) tests linearity across L, and whether k differs by group or tracks QC/registration metrics;
#   (3) re-fits the frozen M1 + decision rule (frozen helper functions, unchanged) on the 80 first-usable youths with
#       the recovered change added to the PS youths only, i.e. "what would we have concluded if PS youth had a true
#       deficit of L ppb?";
#   (4) translates the observed-scale exclusion into true-scale units (bound / k).
# Inputs: raw/roi_qsmxt_{primary,rec_*}.csv (coded, copied from server), mirror/unlocked/key_plain.csv,
#         mirror/unlocked/merged_qsmxt_primary_SNVTA.csv. Output: results.json, tables, README numbers.
suppressPackageStartupMessages({ library(splines); library(sandwich); library(MASS); library(jsonlite) })
args <- commandArgs(trailingOnly = TRUE)
ROOT <- "/path/to/Penn_LEAD/LEADchi"
OUT <- file.path(ROOT, "results_posthoc", "recovery")
RDIR <- file.path(ROOT, "mirror", "scripts", "analysis", "R")
for (f in c("constants.R", "utils.R", "design.R", "plausibility.R", "m1.R")) source(file.path(RDIR, f))
RL <- 0.938278557806515   # frozen gate JSON r_L (PC4)

key <- read.csv(file.path(ROOT, "mirror", "unlocked", "key_plain.csv"), stringsAsFactors = FALSE)
mer <- read.csv(file.path(ROOT, "mirror", "unlocked", "merged_qsmxt_primary_SNVTA.csv"), stringsAsFactors = FALSE)
mer <- code_vars(mer)
rd <- function(v) {
  x <- read.csv(file.path(OUT, "raw", sprintf("roi_qsmxt_%s.csv", v)), stringsAsFactors = FALSE)
  x$ses_order <- as.integer(x$ses_order); x
}
prim <- rd("primary")
levels_tab <- data.frame(variant = c("rec_m05", "rec_m10", "rec_m14", "rec_m20", "rec_p10"),
                         L = c(-5, -10, -14.24, -20, 10))
levels_tab <- levels_tab[file.exists(file.path(OUT, "raw", sprintf("roi_qsmxt_%s.csv", levels_tab$variant))), ]

base <- merge(prim[, c("code", "ses_order", "status", "chi_SNVTA_raw", "chi_SNVTA", "ref_ppb", "chi_GP_raw", "chi_THAL_raw")],
              key[, c("code", "ses_order", "participant_id", "session_id")], by = c("code", "ses_order"))
long <- do.call(rbind, lapply(seq_len(nrow(levels_tab)), function(i) {
  r <- rd(levels_tab$variant[i])
  m <- merge(base, r[, c("code", "ses_order", "status", "chi_SNVTA_raw", "chi_SNVTA", "ref_ppb", "chi_GP_raw", "chi_THAL_raw")],
             by = c("code", "ses_order"), suffixes = c("", "_rec"))
  m$L <- levels_tab$L[i]; m$variant <- levels_tab$variant[i]
  m$d <- m$chi_SNVTA_raw_rec - m$chi_SNVTA_raw
  m$k <- m$d / m$L
  m$d_ref <- m$ref_ppb_rec - m$ref_ppb
  m$d_GP <- m$chi_GP_raw_rec - m$chi_GP_raw
  m$d_THAL <- m$chi_THAL_raw_rec - m$chi_THAL_raw
  m
}))
long <- long[long$status == "ok" & long$status_rec == "ok" & is.finite(long$d), ]
long$sess_key <- paste0(long$participant_id, "_", long$session_id)
cov <- mer[, c("sess_key", "first_usable", "usable", "PS_e", "ADHD_e", "cell", "age", "dice_min", "motion_gsr", "roi_vol",
               "Q_phase", "qc_grade")]
long <- merge(long, cov, by = "sess_key")
write.csv(long, file.path(OUT, "recovery_long.csv"), row.names = FALSE)

summ <- function(x) c(n = length(x), mean = mean(x), sd = sd(x), median = median(x), q25 = unname(quantile(x, .25)),
                      q75 = unname(quantile(x, .75)), min = min(x), max = max(x))
res <- list(meta = list(run_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), label = "POST HOC recovery test",
                        levels = levels_tab, rL = RL))
res$k_all <- lapply(split(long$k, long$L), summ)
fu <- long[long$first_usable %in% 1, ]
res$k_first_usable <- lapply(split(fu$k, fu$L), summ)
res$d_ref_ppb <- lapply(split(long$d_ref, long$L), summ)
res$d_GP_ppb <- lapply(split(long$d_GP, long$L), summ)
res$d_THAL_ppb <- lapply(split(long$d_THAL, long$L), summ)
# linearity: per session slope of d on L through the origin, and pooled
sl <- sapply(split(long, long$sess_key), function(s) if (nrow(s) >= 3) unname(coef(lm(d ~ 0 + L, data = s))[1]) else NA)
res$k_linear_per_session <- summ(sl[is.finite(sl)])
res$pooled_linear <- {
  f <- lm(d ~ L, data = long); list(slope = unname(coef(f)[2]), intercept = unname(coef(f)[1]),
                                    r2 = summary(f)$r.squared, n = nrow(long))
}
# does k differ by group / track QC? (L = -14.24, first-usable)
k14 <- fu[fu$L == -14.24, ]
if (nrow(k14) > 20) {
  k14 <- set_race3(merge(k14, mer[, c("sess_key", "sex_e", "race3_chr")], by = "sess_key"))
  mk <- fit_ols_hc3(stats::as.formula(paste("k ~", "PS_e * ADHD_e + ns(age, df = 3) + sex_e", race_term(k14))), k14)
  res$k_by_group_L14 <- list(PS = term_test(mk, "PS_e"), ADHD = term_test(mk, "ADHD_e"),
                             cell_means = tapply(k14$k, k14$cell, mean))
  res$k_spearman_L14 <- sapply(c("dice_min", "motion_gsr", "roi_vol", "Q_phase", "age"), function(v) {
    ct <- suppressWarnings(cor.test(k14$k, k14[[v]], method = "spearman")); c(rho = unname(ct$estimate), p = ct$p.value)
  })
}
# group-level counterfactual: add the recovered change to PS youths only, re-fit the frozen M1 + decision rule
S <- build_samples(mer)
first <- S$first
obs <- run_m1_core(first, RL)
res$observed <- obs$dec[c("beta_ppb", "se_hc3", "p_two", "ci95_ppb", "sigma_B", "outcome")]
res$observed$bound_ppb <- obs$dec$B_registered$bound_ppb
res$observed$lower95_one_sided <- obs$dec$ci90_ppb[1]
res$counterfactual <- lapply(seq_len(nrow(levels_tab)), function(i) {
  L <- levels_tab$L[i]
  dd <- long[long$L == L, c("sess_key", "d")]
  f <- merge(first, dd, by = "sess_key", all.x = TRUE)
  miss <- sum(is.na(f$d))
  f$d[is.na(f$d)] <- mean(dd$d, na.rm = TRUE)    # sessions without a valid rec value (none expected): mean change
  f$chi_SNVTA <- f$chi_SNVTA + ifelse(f$PS_e > 0, f$d, 0)
  r <- run_m1_core(f, RL)
  list(L_true_ppb = L, n_missing_rec = miss, mean_recovered_in_PS = mean(f$d[f$PS_e > 0]),
       beta_ppb = r$dec$beta_ppb, se = r$dec$se_hc3, p_two = r$dec$p_two, pB = r$dec$B_registered$pB,
       bound_ppb = r$dec$B_registered$bound_ppb, outcome = r$dec$outcome)
})
# true-scale translation of the observed-scale exclusion
kbar <- res$k_first_usable[["-14.24"]][["mean"]]
if (!is.null(kbar)) {
  res$true_scale <- list(k_mean_L14 = kbar,
                         bound_true_ppb = res$observed$bound_ppb / kbar,
                         lower95_one_sided_true_ppb = res$observed$lower95_one_sided / kbar,
                         lower95_two_sided_true_ppb = res$observed$ci95_ppb[1] / kbar)
}
write_json(res, file.path(OUT, "results.json"), auto_unbox = TRUE, digits = 8, pretty = TRUE)
cat(toJSON(res[c("k_first_usable", "pooled_linear", "k_by_group_L14", "observed", "counterfactual", "true_scale")],
           auto_unbox = TRUE, digits = 4, pretty = TRUE), "\n")
