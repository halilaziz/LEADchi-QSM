#!/usr/bin/env Rscript
# LEAD-chi POST HOC (2026-10-28): quantities for the manuscript that the frozen run_all.R does not print, computed on
# one first-usable sample (the complete-human-QC sample H_all180, or the frozen AI-QC table as a reproduction check).
#
#   Rscript supplement_humanqc.R --merged <merged table> --out <json> [--tag H_all180]
#
# The frozen helpers (mirror/scripts/analysis/R: constants, utils, design, plausibility, m1, secondary, sensitivity)
# are SOURCED read-only and used unchanged: validate_table -> code_vars -> apply_plausibility -> build_samples, then
# run_m1_core / term_test / lin_test / h1_decide / qc_by_group / pc2x / run_h1_sensitivities exactly as run_all.R.
# Code that is NOT a frozen helper, and where it comes from:
#   * fit_y, ps_row, adj_cell_means, raw_by_cell: copied verbatim from results_exploratory/code/tier4_exploratory.R
#     (the code behind manuscript Table 1 adjusted means and Figures 1 and 3);
#   * stratum bound tests: as results_posthoc/engine2/code/engine2_analysis.R bound_tests() (and the manuscript's
#     "computed from RR sec. 3 K1 estimates and HC3 SEs" comment), at the frozen bound and at this sample's bound;
#   * joint nuisance row: M1 + scale(dice_min) + scale(motion_gsr) + scale(ref2_minus_wb_ppb) + scale(icv_seg_ml), as
#     results_posthoc/validation_local/code/validation_local.R (B.1, NUIS);
#   * sens_row_df: the frozen sens_row() body plus two extra output columns (df and the p_B at the frozen bound); it
#     is bound to the name sens_row only inside this R session so that the frozen run_h1_sensitivities() returns them;
#   * engine 2: the logic of results_posthoc/engine2/code/engine2_analysis.R (fit_e2, bound_tests, plausibility rule),
#     with the first-usable sample of this table instead of the frozen 80.
suppressPackageStartupMessages({ library(splines); library(sandwich); library(MASS); library(jsonlite); library(digest) })
Sys.setenv(LANGUAGE = "en"); invisible(Sys.setlocale("LC_MESSAGES", "C"))
args <- commandArgs(trailingOnly = TRUE); opt <- list(); i <- 1
while (i <= length(args)) { opt[[sub("^--", "", args[i])]] <- args[i + 1]; i <- i + 2 }
ROOT <- "/path/to/Penn_LEAD/LEADchi"
CODE <- file.path(ROOT, "mirror/scripts/analysis")
GATES <- file.path(ROOT, "mirror/gates/gates_preunlock.json")
E2RAW <- file.path(ROOT, "results_posthoc/engine2/raw")
KEY <- file.path(ROOT, "mirror/unlocked/key_plain.csv")
TAG <- if (is.null(opt$tag)) "H_all180" else opt$tag
FROZEN_BOUND <- -6.13180464348368
for (f in c("constants.R", "utils.R", "design.R", "plausibility.R", "m1.R", "sensitivity.R")) source(file.path(CODE, "R", f))
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

# ---- data: identical preparation to run_all.R ------------------------------------------------------------------
d0 <- utils::read.csv(opt$merged, stringsAsFactors = FALSE, na.strings = c("NA", ""))
validate_table(d0, read_schema(file.path(CODE, "schema/roi_table_schema.csv")))
d0 <- code_vars(d0)
gates <- jsonlite::read_json(GATES, simplifyVector = TRUE); rL <- gates$rL_primary; pc4 <- isTRUE(gates$PC4)
d <- apply_plausibility(d0)$d
S <- build_samples(d)
first <- set_race3(S$first); first$cell <- factor(first$cell, levels = CELLS)
P <- 12

# ---- tier4_exploratory.R helpers (verbatim) ---------------------------------------------------------------------
m1_rhs_y <- function(dd, vol = "roi_vol", extra = "")
  paste0("PS_e * ADHD_e + ns(age, df = 3) + sex_e", race_term(dd), " + scale(", vol, ") + scale(Q_phase)", extra)
fit_y <- function(dd, y, vol = "roi_vol", extra = "") {
  dd <- dd[!is.na(dd[[y]]) & !is.na(dd[[vol]]), ]
  dd <- set_race3(dd)
  fo <- stats::as.formula(paste(y, "~", m1_rhs_y(dd, vol, extra)))
  m <- fit_ols_hc3(fo, dd); m$data <- dd; m$y <- y; m
}
ps_row <- function(m, label) {
  tt <- term_test(m, "PS_e"); sdy <- stats::sd(m$data[[m$y]])
  data.frame(test = label, outcome = m$y, n = m$n, beta_PS = tt$est, se_hc3 = tt$se, ci_lo = tt$ci95[1], ci_hi = tt$ci95[2],
             df = tt$df, p_two = tt$p_two, beta_per_sd_raw = tt$est / sdy, sd_raw = sdy, sigma_resid = m$sigma,
             lo_sd = tt$ci95[1] / sdy, hi_sd = tt$ci95[2] / sdy, stringsAsFactors = FALSE)
}
adj_cell_means <- function(m) {
  X <- stats::model.matrix(m$fit); xb <- colMeans(X)
  do.call(rbind, lapply(CELLS, function(cl) {
    ps <- if (grepl("PS", cl)) 0.5 else -0.5; ad <- if (grepl("ADHD", cl)) 0.5 else -0.5
    x <- xb; x["PS_e"] <- ps; x["ADHD_e"] <- ad; x["PS_e:ADHD_e"] <- ps * ad
    est <- sum(x * m$b); se <- sqrt(drop(t(x) %*% m$V %*% x)); q <- stats::qt(.975, m$df)
    data.frame(cell = cl, n = sum(m$data$cell == cl), adj_mean = est, se_hc3 = se, lo = est - q * se, hi = est + q * se)
  }))
}
raw_by_cell <- function(dd, v) {
  do.call(rbind, lapply(CELLS, function(cl) {
    x <- dd[[v]][dd$cell == cl]; x <- x[!is.na(x)]
    data.frame(variable = v, cell = cl, n = length(x), mean = mean(x), sd = stats::sd(x), median = stats::median(x),
               q1 = unname(stats::quantile(x, .25)), q3 = unname(stats::quantile(x, .75)), min = min(x), max = max(x))
  }))
}
fixedB <- function(est, se, df, bound = FROZEN_BOUND) list(bound_ppb = bound, tB = (est - bound) / se,
                                                          pB = stats::pt((est - bound) / se, df, lower.tail = FALSE))
R <- list(label = paste0("POST HOC supplement (", TAG, "); not confirmatory; the registered outcome (n = 80) stays B"),
          tag = TAG, merged = normalizePath(opt$merged), merged_sha256 = sha256_file(opt$merged),
          frozen_code_sha256 = vapply(c("constants.R", "utils.R", "design.R", "plausibility.R", "m1.R", "sensitivity.R"),
                                      function(f) sha256_file(file.path(CODE, "R", f)), ""),
          r_version = R.version.string, run_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))

# ---- 1. sample and Table 1 ----------------------------------------------------------------------------------------
a <- S$all
subj <- a[!duplicated(a$participant_id), ]
qc <- qc_by_group(S, rL, pc4)
cnt <- function(x) as.list(setNames(as.integer(table(factor(x, levels = CELLS))), CELLS))
msd <- function(v) do.call(rbind, lapply(CELLS, function(cl) { x <- first[[v]][first$cell == cl]; x <- x[!is.na(x)]
  data.frame(cell = cl, n = length(x), n_missing = sum(first$cell == cl) - length(x), mean = mean(x), sd = stats::sd(x),
             min = min(x), max = max(x)) }))
rs <- strsplit(ifelse(is.na(a$fail_reason), "", a$fail_reason), ";")
comp <- sort(unique(unlist(rs))); comp <- comp[comp != ""]
mSN <- fit_y(first, "chi_SNVTA"); mR2 <- fit_y(first, "r2s_SNVTA")
T1 <- list(
  sessions_total = nrow(a), sessions_usable = sum(a$usable == 1), sessions_not_usable = sum(a$usable == 0),
  youths_with_megre = cnt(subj$cell), usable_sessions = cnt(S$usable$cell), analysed_first_usable = cnt(first$cell),
  n_analysed = nrow(first),
  earliest_megre_not_usable = as.list(setNames(qc$table_earliest[CELLS, "earliest failed"], CELLS)),
  analysed_at_ses2_or_ses3 = cnt(first$cell[first$session_id != "ses-1"]),
  analysed_session = as.list(table(first$session_id)),
  age = msd("age"), female = cnt(first$cell[first$sex == "F"]), male = cnt(first$cell[first$sex == "M"]),
  race_black = cnt(first$cell[first$race3_chr == "Black"]), race_white = cnt(first$cell[first$race3_chr == "White"]),
  race_other = cnt(first$cell[first$race3_chr == "Other"]),
  threshold_psychotic_disorder = cnt(first$cell[first$dx_psychosis %in% 1]),
  PRIME_total = msd("PRIME_total"), BDI_total = msd("BDI_total"), ARI_total = msd("ARI_total"),
  chi_SNVTA_raw = raw_by_cell(first, "chi_SNVTA"), chi_SNVTA_adjusted = adj_cell_means(mSN),
  r2s_SNVTA_raw = raw_by_cell(first, "r2s_SNVTA"), r2s_SNVTA_adjusted = adj_cell_means(mR2),
  exclusion_reason_component_counts = as.list(vapply(comp, function(z) sum(vapply(rs, function(r) z %in% r, TRUE)), 1L)),
  exclusion_reason_combinations = as.list(table(ifelse(is.na(a$fail_reason), "(usable)", a$fail_reason))),
  qc_by_group_fisher_p_not_in_H1 = qc$fisher_p_in_H1, qc_by_group_fisher_p_earliest_failed = qc$fisher_p_earliest,
  qc_by_group_tables = list(in_H1 = qc$table_in_H1, earliest = qc$table_earliest),
  age_PS_minus_nonPS_years = mean(first$age[first$PS_e > 0]) - mean(first$age[first$PS_e < 0]),
  age_PS_welch_p = stats::t.test(age ~ PS_e > 0, data = first)$p.value,
  first_usable_keys = sort(first$sess_key))

# ---- 2. H1 extras ---------------------------------------------------------------------------------------------------
core <- run_m1_core(first, rL, pc4); m <- core$m; dec <- core$dec; tt <- term_test(m, "PS_e")
H1x <- list(n = m$n, df = m$df, beta_PS = tt$est, se_hc3 = tt$se, ci95 = tt$ci95, p_two = tt$p_two,
            one_sided_95_lower_limit = tt$ci90[1], delta = dec$delta, sigma_resid = m$sigma, sigma_rob = core$huber$scale,
            sigma_B = dec$sigma_B, huber = core$huber[c("est", "se", "scale", "converged")],
            bound_frozen = fixedB(tt$est, tt$se, m$df),
            bound_recomputed = c(dec$B_registered, list(note = "delta_B 0.581 x sigma_B of this sample")),
            outcome_descriptive = dec$outcome, disconfirming = dec$disconfirming)
# K1 strata with bound tests (post hoc), Holm over the two simple effects (post hoc, as in the manuscript)
k1 <- k1_labels(m, dec)
strat <- function(z) list(est = z$est, se = z$se, p_two = z$p_two, ci95 = z$ci95,
                          lower95_one_sided = z$est - stats::qt(.95, m$df) * z$se,
                          pB_frozen_bound = fixedB(z$est, z$se, m$df)$pB,
                          pB_recomputed_bound = fixedB(z$est, z$se, m$df, dec$B_registered$bound_ppb)$pB)
ph <- stats::p.adjust(c(k1$simple_PSonly_minus_TD$p_two, k1$simple_PSADHD_minus_ADHDonly$p_two), "holm")
H1x$K1 <- list(PSonly_minus_TD = c(strat(k1$simple_PSonly_minus_TD), list(p_holm_2 = ph[1])),
               PSADHD_minus_ADHDonly = c(strat(k1$simple_PSADHD_minus_ADHDonly), list(p_holm_2 = ph[2])),
               ADHD_main = k1$ADHD[c("est", "se", "p_two", "ci95")], PSxADHD = k1$interaction[c("est", "se", "p_two", "ci95")],
               tost = k1$tost, label = k1$label)
# K3 CSF-referenced M1 with its SE (for Figure 2; run_all.R stores beta and p only)
first$dref <- first$chi_SNVTA - first$chi_SNVTA_ref2
d3 <- set_race3(first[!is.na(first$dref), ]); cc <- run_m1_core(d3, rL, pc4, y = "chi_SNVTA_ref2"); t3 <- term_test(cc$m, "PS_e")
H1x$K3_csf <- list(n = cc$m$n, df = cc$m$df, beta = t3$est, se = t3$se, ci95 = t3$ci95, p_two = t3$p_two,
                   pB_recomputed = cc$dec$B_registered$pB, pB_frozen_bound = fixedB(t3$est, t3$se, cc$m$df)$pB)
R$H1_extras <- H1x

# ---- 4. sensitivities with df and the frozen-bound p_B ----------------------------------------------------------------
sens_row_df <- function(name, d, rL, pc4, y = "chi_SNVTA", extra = "", note = "") {
  r <- try_fit({
    d <- set_race3(d); fo <- m1_formula(y, d, extra); m <- fit_ols_hc3(fo, d)
    hb <- huber_fit(fo, d, "PS_e", ols_resid = stats::residuals(m$fit)); dec <- h1_decide(m, hb$scale, rL, pc4)
    data.frame(sensitivity = name, n = m$n, beta_ppb = dec$beta_ppb, se_hc3 = dec$se_hc3, p_two = dec$p_two,
               pB = dec$B_registered$pB, outcome = dec$outcome, note = note, df = m$df,
               bound_recomputed_ppb = dec$B_registered$bound_ppb,
               pB_frozen_bound = fixedB(dec$beta_ppb, dec$se_hc3, m$df)$pB, stringsAsFactors = FALSE)
  }, paste("sensitivity", name))
  if (is.null(r)) data.frame(sensitivity = name, n = nrow(d), beta_ppb = NA, se_hc3 = NA, p_two = NA, pB = NA,
                             outcome = "fit failed", note = note, df = NA, bound_recomputed_ppb = NA, pB_frozen_bound = NA,
                             stringsAsFactors = FALSE) else r
}
sens_row_frozen <- sens_row; sens_row <- sens_row_df
not_run <- function(name, why) data.frame(sensitivity = name, n = NA, beta_ppb = NA, se_hc3 = NA, p_two = NA, pB = NA,
                                          outcome = "not run", note = why, df = NA, bound_recomputed_ppb = NA,
                                          pB_frozen_bound = NA, stringsAsFactors = FALSE)
fill <- function(x) { for (k in c("df", "bound_recomputed_ppb", "pB_frozen_bound")) if (!k %in% names(x)) x[[k]] <- NA; x }
rbind_fill <- function(lst) do.call(rbind, lapply(lst, fill))
# run_h1_sensitivities builds two rows inline (without race; rank ANCOVA); rbind needs equal columns -> wrap rbind
.rbind_base <- base::rbind
SE <- local({
  rbind <- function(...) { l <- list(...); if (all(vapply(l, is.data.frame, TRUE))) .rbind_base(rbind_fill(l)) else .rbind_base(...) }
  environment(run_h1_sensitivities) <- environment()
  run_h1_sensitivities(S, gates, REG$seed, NULL)$table
})
sens_row <- sens_row_frozen
# fill df / frozen-bound p for the inline 'without race' row by an explicit refit with the same formula
ix <- which(SE$sensitivity == "without race")
if (length(ix)) { dd <- set_race3(first); fo <- stats::as.formula(sub(" + race3", "", paste("chi_SNVTA ~", m1_rhs(dd)), fixed = TRUE))
  mm <- fit_ols_hc3(fo, dd); hb <- huber_fit(fo, dd, "PS_e", stats::residuals(mm$fit)); dd2 <- h1_decide(mm, hb$scale, rL, pc4)
  SE$df[ix] <- mm$df; SE$bound_recomputed_ppb[ix] <- dd2$B_registered$bound_ppb
  SE$pB_frozen_bound[ix] <- fixedB(dd2$beta_ppb, dd2$se_hc3, mm$df)$pB }
# extra rows requested for the manuscript (post hoc): session-1-only (with SE), joint nuisance model
s1o <- sens_row_df("session-1-only first-usable (mandatory co-report)", S$ses1_only, rL, pc4)
NUIS <- c("dice_min", "motion_gsr", "ref2_minus_wb_ppb", "icv_seg_ml")
fj <- first[stats::complete.cases(first[, NUIS]), ]
joint <- sens_row_df("M1 + Dice + motion + CSF offset + ICV (joint; post hoc, validation_local B.1)", fj, rL, pc4,
                     extra = paste0(" + scale(", NUIS, ")", collapse = ""))
R$sensitivity <- list(battery_frozen_function = SE, ses1_only = s1o, joint_nuisance = joint,
                      threshold_psychosis_removed = sum(first$dx_psychosis %in% 1),
                      stimulant_endorsers_removed = sum(first$stim_prescribed %in% 1))

# ---- 6. PC2x as the frozen pc2x() + two-sided p ------------------------------------------------------------------
pc <- pc2x(first)
R$PC2x <- lapply(pc$tests, function(z) { dd <- first[!is.na(first[[z$roi]]), ]
  f <- stats::lm(stats::as.formula(paste(z$roi, "~ age")), data = dd); s <- summary(f)$coefficients["age", ]
  c(z, list(se = unname(s[2]), p_two_sided = unname(s[4]), ci95 = unname(stats::confint(f)["age", ]))) })
R$PC2x$pass <- pc$pass

# ---- 7. Figure 3: PS effect by ROI (raw-SD units) and R2* by cell -----------------------------------------------------
rois <- c(chi_SNVTA = "SN-VTA (primary)", chi_SNonly = "SN only (SNc+SNr)", chi_VTAPBP = "VTA + PBP", chi_GP = "globus pallidus",
          chi_PUT = "putamen", chi_CAU = "caudate", chi_STR = "striatum (CAU+PUT+NAcc)", chi_THAL = "thalamus", chi_RN = "red nucleus",
          chi_STN = "subthalamic nucleus", chi_DN_R = "right dentate", chi_DN_L = "left dentate", chi_DGM = "deep GM composite (CAU+GP)")
R$fig3A_regional <- do.call(rbind, lapply(names(rois), function(v) ps_row(fit_y(first, v), rois[[v]])))
tR <- term_test(mR2, "PS_e")
R$fig3B_r2s <- list(raw = T1$r2s_SNVTA_raw, adjusted = T1$r2s_SNVTA_adjusted,
                    beta_PS = tR$est, se = tR$se, ci95 = tR$ci95, p_two = tR$p_two, p_one_sided_lower = tR$p_lower, n = mR2$n)
R$fig_points <- first[, c("sub", "session_id", "cell", "age", "sex", "dx_psychosis", "chi_SNVTA", "r2s_SNVTA")]
R$table1 <- T1

# ---- 8. engine 2 (engine2_analysis.R logic on this sample's first-usable sessions) ------------------------------------
key <- utils::read.csv(KEY, stringsAsFactors = FALSE)
mer <- code_vars(utils::read.csv(opt$merged, stringsAsFactors = FALSE))   # engine2_analysis.R: no plausibility step
fe <- build_samples(mer)$first
rd <- function(v) {
  x <- utils::read.csv(file.path(E2RAW, sprintf("roi_qsmxt_%s.csv", v)), stringsAsFactors = FALSE)
  x$ses_order <- as.integer(x$ses_order)
  x <- merge(x, key[, c("code", "ses_order", "participant_id", "session_id")], by = c("code", "ses_order"))
  x$sess_key <- paste0(x$participant_id, "_", x$session_id); x
}
bound_tests <- function(r, df) {
  bnd <- r$dec$B_registered$bound_ppb; tq <- stats::qt(.95, df)
  lapply(list(PSonly_minus_TD = lin_test(r$m, c(PS_e = 1, "PS_e:ADHD_e" = -0.5), "PS-only - TD"),
              PSADHD_minus_ADHD = lin_test(r$m, c(PS_e = 1, "PS_e:ADHD_e" = 0.5), "PS+ADHD - ADHD-only")),
         function(z) list(est = z$est, se = z$se, p_two = z$p_two, ci95 = z$ci95, lower95_one_sided = z$est - tq * z$se,
                          pB = stats::pt((z$est - bnd) / z$se, df, lower.tail = FALSE)))
}
fit_e2 <- function(mm, keep, label) {
  d2 <- mm[keep, ]; d2$chi_SNVTA <- d2$e2_chi_SNVTA
  r2 <- run_m1_core(d2, rL)
  list(label = label, n = nrow(d2), cells = as.list(cell_counts(d2)),
       M1 = r2$dec[c("beta_ppb", "se_hc3", "p_two", "ci95_ppb", "sigma_resid", "sigma_B", "delta", "outcome", "disconfirming")],
       bound_ppb_registered_deltaB = r2$dec$B_registered$bound_ppb, pB = r2$dec$B_registered$pB,
       huber_beta = r2$huber$est, interaction = term_test(r2$m, "PS_e:ADHD_e")[c("est", "se", "p_two", "ci95")],
       strata = bound_tests(r2, r2$m$df))
}
E2 <- list(note = paste("engine-2 'B' entries are a B pattern, not a registered outcome; the bound uses delta_B from the",
                        "primary engine's split-echo reliability and engine 2's own sigma_B (as engine2_analysis.R)"),
           n_first_usable = nrow(fe))
for (inv in c("rts", "tv")) {
  e2 <- rd(paste0("e2_lap_vsharp_", inv))
  cols <- c("sess_key", "status", "chi_SNVTA", "fail_SNVTA", "inmask_SNVTA")
  mm <- merge(fe, setNames(e2[, cols], c("sess_key", paste0("e2_", cols[-1]))), by = "sess_key", all.x = TRUE)
  finite <- is.finite(mm$e2_chi_SNVTA); valid <- finite & !(mm$e2_fail_SNVTA %in% 1)
  E2[[inv]] <- list(n_finite = sum(finite), n_valid = sum(valid), n_fail_plausibility = sum(finite & mm$e2_fail_SNVTA %in% 1),
                    missing_status = as.list(table(mm$e2_status[!finite], useNA = "ifany")),
                    main = fit_e2(mm, valid, "main: engine-2 plausibility failures missing (registered sec. 12.4 rule)"),
                    sensitivity_failures_kept = fit_e2(mm, finite, "sensitivity: engine-2 SN-VTA plausibility failures kept"))
}
R$engine2 <- E2
jsonlite::write_json(R, opt$out, auto_unbox = TRUE, digits = NA, pretty = TRUE, na = "null")
cat("wrote", opt$out, "\n")
