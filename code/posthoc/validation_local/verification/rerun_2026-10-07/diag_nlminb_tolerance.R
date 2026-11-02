#!/usr/bin/env Rscript
# DIAGNOSTIC ONLY (POST HOC, 2026-10-07): can tighter nlminb tolerances bring nlme closer to the verifier's REML optimum?
# Refits every D lme model (54 + 18 within/between) under four nlminb settings and reports the worst relative difference
# from verification/mine/D_age_slopes.csv, plus every fit that stops with an error. Reads only; prints to stdout.
suppressPackageStartupMessages({library(splines); library(sandwich); library(MASS); library(nlme); library(jsonlite); library(digest)})
ROOT <- "/path/to/Penn_LEAD/LEADchi"; CODE <- file.path(ROOT, "mirror/scripts/analysis")
for (f in c("constants.R", "utils.R", "design.R", "plausibility.R", "m1.R")) source(file.path(CODE, "R", f))
d0 <- utils::read.csv(file.path(ROOT, "mirror/unlocked/merged_qsmxt_primary_SNVTA.csv"), stringsAsFactors = FALSE, na.strings = c("NA", ""))
d0 <- code_vars(d0); pl <- apply_plausibility(d0); S <- build_samples(pl$d); usable <- S$usable
usable$chi_PUT_unref <- usable$chi_PUT + usable$ref_ppb; usable$chi_PUT_csfref <- usable$chi_PUT - usable$ref2_minus_wb_ppb
usable$chi_GP_unref <- usable$chi_GP + usable$ref_ppb; usable$chi_GP_csfref <- usable$chi_GP - usable$ref2_minus_wb_ppb
usable$chi_para_PUT <- (usable$chi_PUT + 1000*usable$r2s_PUT/274)/2; usable$chi_dia_PUT <- (usable$chi_PUT - 1000*usable$r2s_PUT/274)/2
usable$chi_para_GP <- (usable$chi_GP + 1000*usable$r2s_GP/274)/2; usable$chi_dia_GP <- (usable$chi_GP - 1000*usable$r2s_GP/274)/2
mine <- read.csv(file.path(ROOT, "results_posthoc/validation_local/verification/mine/D_age_slopes.csv"), stringsAsFactors = FALSE)
Y <- c("r2s_PUT","chi_PUT","chi_PUT_raw","r2s_GP","chi_GP","chi_GP_raw","chi_PUT_unref","chi_PUT_csfref","chi_GP_unref","chi_GP_csfref","ref_ppb","ref2_minus_wb_ppb","chi_para_PUT","chi_dia_PUT","chi_para_GP","chi_dia_GP","r2s_SNVTA","chi_SNVTA")
COVS <- c(`age only` = "", `+ sex` = " + sex_e", `+ sex + cell` = " + sex_e + PS_e * ADHD_e")
CTL <- list(rel12 = lmeControl(rel.tol = 1e-12), rel13 = lmeControl(rel.tol = 1e-13), rel12_x = lmeControl(rel.tol = 1e-12, x.tol = 1e-12),
            rel12_iter = lmeControl(rel.tol = 1e-12, msMaxIter = 1000, msMaxEval = 2000))
for (cn_ in names(CTL)) {
  worst <- 0; wn <- ""
  for (y in Y) {
    du <- usable[!is.na(usable[[y]]), ]
    for (cn in names(COVS)) {
      m <- tryCatch(lme(as.formula(paste(y, "~ age", COVS[[cn]])), random = ~1|sub, data = du, method = "REML", control = CTL[[cn_]]), error = function(e) {cat("ERR", y, cn, conditionMessage(e), "\n"); NULL}); if (is.null(m)) next
      tb <- summary(m)$tTable
      r <- mine[mine$outcome == y & mine$sample == "all usable sessions (lme REML, random intercept)" & mine$covariates == cn, ]
      for (col in c("Value","Std.Error","p-value")) {
        mv <- switch(col, Value = r$slope_per_y, Std.Error = r$se, `p-value` = r$p_two)
        rel <- abs(tb["age", col] - mv) / abs(mv); if (abs(tb["age", col] - mv) > 1e-9 && rel > worst) { worst <- rel; wn <- paste(y, cn, col) }
      }
    }
    du$age_b <- ave(du$age, du$sub); du$age_w <- du$age - du$age_b
    m <- tryCatch(lme(as.formula(paste(y, "~ age_w + age_b + sex_e + PS_e * ADHD_e")), random = ~1|sub, data = du, method = "REML", control = CTL[[cn_]]), error = function(e) {cat("ERR wb", y, conditionMessage(e), "\n"); NULL}); if (is.null(m)) next
    tb <- summary(m)$tTable
    for (tm in c("age_w","age_b")) {
      r <- mine[mine$outcome == y & mine$sample == paste("all usable (lme REML),", if (tm == "age_w") "WITHIN-person age" else "BETWEEN-person mean age"), ]
      for (col in c("Value","Std.Error","p-value")) {
        mv <- switch(col, Value = r$slope_per_y, Std.Error = r$se, `p-value` = r$p_two)
        rel <- abs(tb[tm, col] - mv) / abs(mv); if (abs(tb[tm, col] - mv) > 1e-9 && rel > worst) { worst <- rel; wn <- paste(y, tm, col) }
      }
    }
  }
  cat(sprintf("%-10s worst rel diff vs verifier %.3g (%s)\n", cn_, worst, wn))
}
