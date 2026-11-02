#!/usr/bin/env Rscript
# LEAD-chi Annex C.4 (2026-10-07): K2 and PC3 with text-conforming R2* ROI means (same voxel set as the chi ROI)
# versus the hashed-code version. Frozen helper functions unchanged; r_L as frozen.
suppressPackageStartupMessages({ library(splines); library(sandwich); library(MASS); library(jsonlite) })
ROOT <- "/path/to/Penn_LEAD/LEADchi"; OUT <- file.path(ROOT, "results_posthoc", "r2s_textconform")
for (f in c("constants.R", "utils.R", "design.R", "plausibility.R", "m1.R")) source(file.path(ROOT, "mirror/scripts/analysis/R", f))
key <- read.csv(file.path(ROOT, "mirror/unlocked/key_plain.csv"), stringsAsFactors = FALSE)
mer <- code_vars(read.csv(file.path(ROOT, "mirror/unlocked/merged_qsmxt_primary_SNVTA.csv"), stringsAsFactors = FALSE))
r2 <- merge(read.csv(file.path(OUT, "r2s_roi_code_vs_text.csv"), stringsAsFactors = FALSE),
            key[, c("code", "ses_order", "participant_id", "session_id")], by = c("code", "ses_order"))
r2$sess_key <- paste0(r2$participant_id, "_", r2$session_id)
first <- build_samples(mer)$first
d <- merge(first, r2[, c("sess_key", grep("^r2s_", names(r2), value = TRUE))], by = "sess_key")
res <- list(n = nrow(d))
res$code_reproduces_frozen_r2s_SNVTA <- max(abs(d$r2s_SNVTA - d$r2s_SNVTA_code), na.rm = TRUE)
res$code_reproduces_frozen_r2s_GP <- max(abs(d$r2s_GP - d$r2s_GP_code), na.rm = TRUE)
res$mean_text_minus_code <- sapply(c("SNVTA", "GP", "PUT"), function(n) mean(d[[paste0("r2s_", n, "_text")]] - d[[paste0("r2s_", n, "_code")]], na.rm = TRUE))
k2 <- function(col) {
  dd <- d; dd$y <- dd[[col]]; dd <- set_race3(dd[is.finite(dd$y), ])
  m <- fit_ols_hc3(m1_formula("y", dd), dd); tt <- term_test(m, "PS_e")
  list(n = m$n, beta = tt$est, se = tt$se, ci95 = tt$ci95, p_one_sided_lower = tt$p_lower, p_two = tt$p_two,
       K2_met = isTRUE(tt$est < 0 && tt$p_lower < .05))
}
res$K2_code <- k2("r2s_SNVTA_code"); res$K2_text <- k2("r2s_SNVTA_text")
pc3 <- function(col) { ok <- is.finite(d$chi_GP) & is.finite(d[[col]]); list(n = sum(ok), r = cor(d$chi_GP[ok], d[[col]][ok])) }
res$PC3_code <- pc3("r2s_GP_code"); res$PC3_text <- pc3("r2s_GP_text")
write_json(res, file.path(OUT, "results.json"), auto_unbox = TRUE, digits = 8, pretty = TRUE)
cat(toJSON(res, auto_unbox = TRUE, digits = 5, pretty = TRUE), "\n")
