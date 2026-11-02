#!/usr/bin/env Rscript
# Guard and edge-case tests for the LEAD-chi analysis code (synthetic inputs only).
#   Part 1: run_all.R must REFUSE inconsistent inputs (no results.json written).
#   Part 2: registered edge-case behaviour (r_L <= 0 or PC4 failure -> B not evaluable; rlm non-convergence ->
#           registered sigma_rob fallback; HC3 with a leverage of 1 -> NA SE, no error).
#   Rscript tests/test_guards.R
Sys.setenv(LANGUAGE = "en")
all_args <- commandArgs(trailingOnly = FALSE)
CODE_DIR <- dirname(dirname(normalizePath(sub("^--file=", "", all_args[grep("^--file=", all_args)][1]))))
suppressPackageStartupMessages(library(jsonlite))
tab0 <- utils::read.csv(file.path(CODE_DIR, "synthetic/SYNTHETIC_roi_table_null.csv"), stringsAsFactors = FALSE)
g0 <- jsonlite::read_json(file.path(CODE_DIR, "synthetic/SYNTHETIC_gates_null.json"), simplifyVector = TRUE)
mv0 <- file.path(CODE_DIR, "synthetic/SYNTHETIC_multiverse_null.csv")
tmp <- tempfile("leadchi_guard_"); dir.create(tmp)
FAST <- c("--n-perm", "99", "--n-perm-spec", "49", "--n-boot", "49", "--n-pb", "19", "--n-wild", "99")
run <- function(name, tab, g, mv = mv0) {
  ti <- file.path(tmp, paste0(name, ".csv")); gi <- file.path(tmp, paste0(name, ".json"))
  utils::write.csv(tab, ti, row.names = FALSE, na = "NA"); jsonlite::write_json(g, gi, auto_unbox = TRUE, digits = NA)
  a <- c(file.path(CODE_DIR, "run_all.R"), "--input", ti, "--gates", gi, "--out", file.path(tmp, name), FAST)
  if (!is.null(mv)) a <- c(a, "--multiverse", mv)
  list(out = suppressWarnings(system2("Rscript", a, stdout = TRUE, stderr = TRUE)), res = file.path(tmp, name, "results.json"))
}
refuse <- function(name, tab, g, expect, mv = mv0) {
  r <- run(name, tab, g, mv)
  ok <- any(grepl(expect, r$out, fixed = TRUE)) && !file.exists(r$res)
  cat(sprintf("%-44s %s\n", name, if (ok) "PASS (refused as expected)" else "FAIL"))
  if (!ok) cat(tail(r$out, 5), sep = "\n")
  ok
}
drop_field <- function(g, k) { g[[k]] <- NULL; g }
cat("Part 1: refusals\n")
res <- c(
  refuse("wrong_engine_for_G2_fail", tab0, modifyList(g0, list(G2 = FALSE)), "wrong input table for this branch"),
  refuse("wrong_roi_for_G3_fail", tab0, modifyList(g0, list(G3 = FALSE)), "but G3 selects"),
  refuse("constants_roi_mismatch", tab0, modifyList(g0, list(constants_roi = "SN_only")), "gate constants were computed for"),
  refuse("origin_mismatch", tab0, modifyList(g0, list(data_origin = "real")), "gates data_origin"),
  refuse("gate_rL_missing", tab0, drop_field(g0, "rL_primary"), "gate JSON lacks required fields"),
  refuse("gate_rL_above_1", tab0, modifyList(g0, list(rL_primary = 1.2)), "must be a finite number"),
  refuse("gate_A2_missing", tab0, drop_field(g0, "A2_approved"), "gate JSON lacks required fields"),
  refuse("gate_pooled_sd_zero", tab0, modifyList(g0, list(pooled = modifyList(g0$pooled, list(sd_total_SNVTA = 0)))), "pooled.sd_total_SNVTA"),
  refuse("gate_PC_fail_without_fallback_flag", tab0, modifyList(g0, list(PC1 = FALSE)), "fallback_PC_pass"),
  { t <- tab0; t$chi_SNVTA[which(t$first_usable == 1)[1]] <- 400
    refuse("snvta_outside_envelope_usable", t, g0, "outside the registered envelope") },
  { t <- tab0; i <- which(t$first_usable == 1 & duplicated(t$participant_id, fromLast = TRUE))[1]; t$first_usable[i] <- 0
    t$first_usable[i + 1] <- 1; refuse("first_usable_flag_mismatch", t, g0, "frozen first_usable flags disagree") },
  { t <- tab0; t$Q_phase <- NULL; refuse("missing_required_column", t, g0, "missing required columns") },
  { t <- rbind(tab0, tab0[1, ]); refuse("duplicate_session_row", t, g0, "duplicate (participant_id, session_id)") },
  { t <- tab0; t$chi_SNVTA[which(t$usable == 1)[2]] <- NA; refuse("usable_row_missing_outcome", t, g0, "usable rows must have non-missing chi_SNVTA") },
  { t <- tab0; t$ref2_type[which(t$usable == 1)[1]] <- "PLIC"; refuse("plic_reference_not_allowed", t, g0, "ref2_type must be") },
  { t <- tab0; t$data_origin <- "real"; refuse("real_data_without_multiverse", t, modifyList(g0, list(data_origin = "real")),
                                               "--multiverse is required", mv = NULL) })
cat("\nPart 2: registered edge-case behaviour\n")
chk <- function(name, ok) { cat(sprintf("%-44s %s\n", name, if (isTRUE(ok)) "PASS" else "FAIL")); isTRUE(ok) }
r1 <- run("rL_negative_B_not_evaluable", tab0, modifyList(g0, list(rL_primary = -0.1)))
j1 <- if (file.exists(r1$res)) jsonlite::read_json(r1$res) else NULL
r2 <- run("PC4_fail_B_not_evaluable", tab0, modifyList(g0, list(PC4 = FALSE)))
j2 <- if (file.exists(r2$res)) jsonlite::read_json(r2$res) else NULL
res2 <- c(
  chk("r_L < 0: run completes, B not evaluable", !is.null(j1) && identical(j1$H1$decision$B_evaluable, FALSE) &&
        identical(j1$H1$precision$label, "low-precision measure") && !(j1$H1$decision$outcome %in% c("B", "A+B"))),
  chk("PC4 fail: B not evaluable", !is.null(j2) && identical(j2$H1$decision$B_evaluable, FALSE)))
# unit tests on the registered functions
for (f in c("constants.R", "utils.R", "design.R", "plausibility.R", "m1.R")) source(file.path(CODE_DIR, "R", f))
suppressPackageStartupMessages({ library(splines); library(sandwich) })
set.seed(1); n <- 60; dd <- data.frame(y = stats::rnorm(n), x = stats::rnorm(n))
REG$huber_maxit <- 1
hb <- suppressWarnings(huber_fit(y ~ x, dd, "x", ols_resid = stats::residuals(stats::lm(y ~ x, dd))))
res2 <- c(res2, chk("rlm non-convergence -> MAD fallback, flagged", !hb$converged && is.finite(hb$scale) &&
                      grepl("fallback", hb$scale_source)))
dl <- data.frame(y = stats::rnorm(20), x = c(1, rep(0, 19)), z = stats::rnorm(20))
mh <- fit_ols_hc3(y ~ x + z, dl)
res2 <- c(res2, chk("HC3 with leverage 1 -> NA SE, no error", is.na(mh$se["x"]) || !is.finite(mh$se["x"])))
all_ok <- c(res, res2)
cat(sprintf("\n%d of %d tests passed\n", sum(all_ok), length(all_ok)))
unlink(tmp, recursive = TRUE)
if (!all(all_ok)) quit(status = 1)
