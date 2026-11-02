#!/usr/bin/env Rscript
# LEAD-chi -- POST HOC / registered-exploratory local validation, tasks A-D (run after the 2026-10-05 unlock).
#   A  X2 chi-separation (registered Tier-4 exploratory, prereg Q25; Dimov et al. 2022 Tomography 8:1544-1551)
#   B  joint nuisance model and Dice-quality analyses (POST HOC)
#   C  low-TD-tail check (POST HOC)
#   D  putamen / GP developmental positive-control diagnostics (POST HOC)
# NOTHING here re-decides the frozen confirmatory outcome (H1 final outcome B; results_confirmatory/report.md).
#
# Revision 2026-10-07 (after the independent verification in ../verification/VERIFY.md; details in ../CHANGES_2026-10-07.md):
#   * every nlme::lme fit now uses the nlme default optimiser (nlminb); the earlier optimiser setting (opt = "optim")
#     stopped short of the REML optimum (VERIFY.md sec. 3.1, item V-3.1);
#   * new outputs: D_putamen_shortfall_by_Dr.csv (review item W11) and A_X2_BH_family_sensitivity.csv (observation O1);
#   * results.json -> meta records the declared deviations (X2 constants fixed after unlock; R2* voxel set differs from
#     the plan text, prereg line 253; BH family definition is a choice). Still POST HOC; nothing here is a finding.
#
# Conventions: the frozen analysis helpers are sourced READ-ONLY from mirror/scripts/analysis/R, so the M1 design
# (PS_e * ADHD_e + ns(age, df = 3) + sex_e + race3 [contr.sum over levels present] + scale(roi_vol) + scale(Q_phase);
# effect codes +/-0.5; OLS + HC3, t reference with df = n - p), the plausibility edits and the first-usable sample are
# identical to run_all.R. Local only: no server access, no image is opened, nothing is written outside OUT.
#
#   Rscript validation_local.R

suppressPackageStartupMessages({
  library(splines); library(sandwich); library(MASS); library(nlme); library(jsonlite); library(digest)
})
Sys.setenv(LANGUAGE = "en")
ROOT  <- "/path/to/Penn_LEAD/LEADchi"
CODE  <- file.path(ROOT, "mirror/scripts/analysis")
INPUT <- file.path(ROOT, "mirror/unlocked/merged_qsmxt_primary_SNVTA.csv")
GATES <- file.path(ROOT, "mirror/gates/gates_preunlock.json")
OUT   <- file.path(ROOT, "results_posthoc/validation_local")
VER   <- file.path(OUT, "verify")
dir.create(VER, recursive = TRUE, showWarnings = FALSE)
for (f in c("constants.R", "utils.R", "design.R", "plausibility.R", "m1.R")) source(file.path(CODE, "R", f))
SELF <- file.path(OUT, "code", "validation_local.R")
log_msg("validation_local start; input ", INPUT)

# ------------------------------------------------------------------------------------------------------------
# 0. Data: identical preparation to run_all.R (validate -> code -> plausibility -> samples)
# ------------------------------------------------------------------------------------------------------------
d0 <- utils::read.csv(INPUT, stringsAsFactors = FALSE, na.strings = c("NA", ""))
validate_table(d0, read_schema(file.path(CODE, "schema/roi_table_schema.csv")))
d0 <- code_vars(d0)
gates <- jsonlite::read_json(GATES, simplifyVector = TRUE)
rL <- gates$rL_primary; pc4 <- isTRUE(gates$PC4)
pl <- apply_plausibility(d0); d <- pl$d
S <- build_samples(d)
first <- set_race3(S$first); usable <- S$usable
stopifnot(nrow(first) == 80, nrow(usable) == 107)
log_msg("first-usable ", cell_str(first), "; usable sessions ", nrow(usable), " from ", length(unique(usable$sub)), " subjects")

# reproduction of the frozen M1 (results_confirmatory/results.json)
core <- run_m1_core(first, rL, pc4)
BOUND <- core$dec$B_registered$bound_ppb
checks <- list()
chk <- function(name, got, want, tol) {
  ok <- isTRUE(abs(got - want) <= tol)
  checks[[name]] <<- list(got = got, frozen = want, tol = tol, pass = ok)
  log_msg(sprintf("check %-40s got %.4f frozen %.4f %s", name, got, want, if (ok) "OK" else "MISMATCH"))
  invisible(ok)
}
chk("M1 beta_PS", core$dec$beta_ppb, 4.9628, 5e-4)
chk("M1 SE_HC3", core$dec$se_hc3, 3.1251, 5e-4)
chk("outcome-B bound (ppb)", BOUND, -6.1318, 5e-4)
chk("M1 t_B", core$dec$B_registered$tB, 3.550, 5e-4)

# ------------------------------------------------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------------------------------------------------
# M1-form fit for outcome y with optional extra RHS terms; complete cases on y and the extra variables
fit_m1y <- function(dd, y, extra = "", need = character(0)) {
  keep <- !is.na(dd[[y]]); for (v in need) keep <- keep & !is.na(dd[[v]])
  dd <- set_race3(dd[keep, ])
  fo <- m1_formula(y, dd, extra)
  m <- fit_ols_hc3(fo, dd); m$data <- dd; m$y <- y; m
}
bound_fixed <- function(tt) {
  tB <- (tt$est - BOUND) / tt$se
  list(bound_ppb = BOUND, tB = tB, df = tt$df, pB = stats::pt(tB, tt$df, lower.tail = FALSE))
}
# frozen-style B test: sigma_B and delta_B recomputed in the model (as sens_row in sensitivity.R does)
frozen_style_pB <- function(m) {
  hb <- huber_fit(m$formula, m$data, "PS_e", ols_resid = stats::residuals(m$fit))
  h1_decide(m, hb$scale, rL, pc4)$B_registered$pB
}
ps_row <- function(m, label, block, bound = FALSE, note = "") {
  tt <- term_test(m, "PS_e")
  s1 <- lin_test(m, c(PS_e = 1, `PS_e:ADHD_e` = -0.5), "PS-only - TD")
  s2 <- lin_test(m, c(PS_e = 1, `PS_e:ADHD_e` = 0.5), "PS+ADHD - ADHD-only")
  b <- if (bound) bound_fixed(tt) else list(tB = NA, pB = NA)
  data.frame(block = block, model = label, outcome = m$y, n = m$n, cells = paste(cell_counts(m$data)[CELLS], collapse = "/"),
             beta_PS = tt$est, se_hc3 = tt$se, ci_lo = tt$ci95[1], ci_hi = tt$ci95[2], t = tt$t, df = tt$df, p_two = tt$p_two,
             PSonly_minus_TD = s1$est, PSonly_minus_TD_se = s1$se, PSonly_minus_TD_lo = s1$ci95[1], PSonly_minus_TD_hi = s1$ci95[2],
             PSonly_minus_TD_p = s1$p_two,
             PSADHD_minus_ADHD = s2$est, PSADHD_minus_ADHD_se = s2$se, PSADHD_minus_ADHD_lo = s2$ci95[1], PSADHD_minus_ADHD_hi = s2$ci95[2],
             PSADHD_minus_ADHD_p = s2$p_two, sigma_resid = m$sigma,
             tB_fixed_bound = b$tB, pB_fixed_bound = b$pB, note = note, stringsAsFactors = FALSE)
}
export_design <- function(m, tag) {           # for the independent statsmodels check (verify/)
  X <- stats::model.matrix(m$fit); y <- stats::model.response(stats::model.frame(m$fit))
  utils::write.csv(data.frame(y = y, X, check.names = FALSE), file.path(VER, paste0("design_", tag, ".csv")), row.names = FALSE)
  data.frame(tag = tag, term = "PS_e", est = unname(m$b["PS_e"]), se_hc3 = unname(m$se["PS_e"]), n = m$n, df = m$df)
}
m1_ref <- fit_m1y(first, "chi_SNVTA")        # identical fit to core$m, with $data / $y attached
stopifnot(isTRUE(all.equal(m1_ref$b, core$m$b)))
ver_rows <- list()
ver_rows[["M1_chi_SNVTA"]] <- export_design(core$m, "M1_chi_SNVTA")
wcsv <- function(x, f) utils::write.csv(x, file.path(OUT, f), row.names = FALSE)

# ------------------------------------------------------------------------------------------------------------
# A. X2 chi-separation (registered Tier-4 exploratory; Dimov et al. 2022)
# ------------------------------------------------------------------------------------------------------------
# Paper (Tomography 2022;8:1544-1551, sec. 2.1, Eq. 1; page within the article not verified): R2* = r+|chi+| + r-|chi-|, f = d * (chi+ + chi-),
# r+ = r- = 274 Hz/ppm (their ref. 22, Dimov et al. 2022 J Neuroimaging); theoretical 321 Hz/ppm (Yablonskiy & Haacke
# 1994, their ref. 23); no separate R2 term in Eq. 1 (the Discussion notes a possibly unaccounted constant, argued to be
# small from the low R2* of ventricular CSF) -> assumed R2 = 0 s^-1. With chi+ > 0, chi- < 0 and equal constants, the
# noiseless point solution is chi_para = (chi + R2'/Dr)/2, chi_dia = (chi - R2'/Dr)/2 with R2' = R2* - R2.
# Units: chi in ppb, R2* in s^-1, Dr in Hz/ppm -> 1000 * R2'/Dr is in ppb.
DR0 <- 274
DR_SET <- c(`Dr x 0.5 (137)` = 0.5 * DR0, `Dr x 1 (274, paper)` = DR0, `theory (321)` = 321, `Dr x 2 (548)` = 2 * DR0)
R2_SET <- c(0, 10, 15)                       # 0 = paper; 10 and 15 s^-1 = arbitrary illustrative constants (cancellation)
X2_ROIS <- list(SNVTA = c("chi_SNVTA", "r2s_SNVTA"), SNonly = c("chi_SNonly", "r2s_SNonly"), GP = c("chi_GP", "r2s_GP"))
sep <- function(dd, chi, r2s, Dr, R2 = 0) {
  q <- 1000 * (dd[[r2s]] - R2) / Dr
  list(para = (dd[[chi]] + q) / 2, dia = (dd[[chi]] - q) / 2)
}
x2_rows <- list(); x2_desc <- list(); x2_ident <- list()
for (roi in names(X2_ROIS)) {
  cv <- X2_ROIS[[roi]][1]; rv <- X2_ROIS[[roi]][2]
  m_chi <- fit_m1y(first, cv, need = rv); m_r2 <- fit_m1y(first, rv, need = cv)
  x2_rows[[length(x2_rows) + 1]] <- ps_row(m_chi, paste0(roi, ": total chi (input)"), "A_X2_inputs", note = "Dr-free")
  x2_rows[[length(x2_rows) + 1]] <- ps_row(m_r2, paste0(roi, ": R2* (input, s^-1)"), "A_X2_inputs", note = "Dr-free")
  for (nm in names(DR_SET)) for (R2 in R2_SET) {
    Dr <- DR_SET[[nm]]
    if (R2 != 0 && Dr != DR0) next          # R2 constants only illustrated at the paper's Dr
    s <- sep(first, cv, rv, Dr, R2); dd <- first; dd$chi_para <- s$para; dd$chi_dia <- s$dia
    for (comp in c("chi_para", "chi_dia")) {
      m <- fit_m1y(dd, comp)
      r <- ps_row(m, sprintf("%s: %s, %s, R2 = %g", roi, comp, nm, R2), "A_X2",
                  note = if (Dr == DR0 && R2 == 0) "primary (paper constants)" else if (R2 != 0) "R2 cancellation check" else "Dr sensitivity")
      r$roi <- roi; r$component <- comp; r$Dr <- Dr; r$R2 <- R2
      # identity check: beta(comp) = (beta(chi) +/- 1000 beta(R2*) / Dr) / 2 (same design matrix)
      sgn <- if (comp == "chi_para") 1 else -1
      r$beta_from_identity <- (m_chi$b["PS_e"] + sgn * 1000 * m_r2$b["PS_e"] / Dr) / 2
      r$mean_all <- mean(dd[[comp]]); r$share_constraint_violated <- if (comp == "chi_para") mean(dd[[comp]] < 0) else mean(dd[[comp]] > 0)
      for (cl in CELLS) r[[paste0("mean_", gsub("[^A-Za-z]", "", cl))]] <- mean(dd[[comp]][dd$cell == cl])
      x2_rows[[length(x2_rows) + 1]] <- r
      if (roi == "SNVTA" && Dr == DR0 && R2 == 0) ver_rows[[paste0("X2_SNVTA_", comp)]] <- export_design(m, paste0("X2_SNVTA_", comp))
    }
  }
  # referencing: chi before whole-brain referencing (chi + ref_ppb) -> both components shift by ref_ppb / 2
  dd <- first; dd$chi_unref <- dd[[cv]] + dd$ref_ppb
  s <- sep(dd, "chi_unref", rv, DR0, 0); dd$chi_para_unref <- s$para; dd$chi_dia_unref <- s$dia
  for (comp in c("chi_para_unref", "chi_dia_unref")) {
    r <- ps_row(fit_m1y(dd, comp), sprintf("%s: %s (chi + ref_ppb), Dr 274", roi, comp), "A_X2_reference",
                note = "whole-brain reference added back; shifts both components by ref_ppb/2")
    r$roi <- roi; r$component <- comp; r$Dr <- DR0; r$R2 <- 0
    x2_rows[[length(x2_rows) + 1]] <- r
  }
}
x2 <- do.call(rbind, lapply(x2_rows, function(r) { for (v in c("roi", "component", "Dr", "R2", "beta_from_identity", "mean_all",
  "share_constraint_violated", paste0("mean_", c("TD", "ADHDonly", "PSonly", "PSADHD")))) if (is.null(r[[v]])) r[[v]] <- NA; r }))
ident_err <- max(abs(x2$beta_PS - x2$beta_from_identity), na.rm = TRUE)
log_msg("X2 identity max |beta - (beta_chi +/- 1000 beta_R2*/Dr)/2| = ", signif(ident_err, 3))
# R2 cancellation: beta_PS identical across R2 at Dr = 274
canc <- x2[x2$block == "A_X2" & x2$Dr == DR0, ]
canc_err <- max(tapply(canc$beta_PS, paste(canc$roi, canc$component), function(b) diff(range(b))))
canc_se_err <- max(tapply(canc$se_hc3, paste(canc$roi, canc$component), function(b) diff(range(b))))
log_msg("X2 R2-cancellation: max range of beta_PS over R2 = ", signif(canc_err, 3), "; of SE = ", signif(canc_se_err, 3))
# BH-FDR within the X2 block: primary family = 6 beta_PS tests (3 ROIs x {para, dia}) at Dr = 274, R2 = 0;
# extended family = the same 6 plus their 12 simple effects
prim <- x2$block == "A_X2" & x2$Dr == DR0 & x2$R2 == 0
x2$q_BH_family6 <- NA; x2$q_BH_family6[prim] <- stats::p.adjust(x2$p_two[prim], "BH")
ext_p <- c(x2$p_two[prim], x2$PSonly_minus_TD_p[prim], x2$PSADHD_minus_ADHD_p[prim])
ext_q <- stats::p.adjust(ext_p, "BH"); k <- sum(prim)
x2$q_BH_ext18_beta <- NA; x2$q_BH_ext18_PSonly <- NA; x2$q_BH_ext18_PSADHD <- NA
x2$q_BH_ext18_beta[prim] <- ext_q[1:k]; x2$q_BH_ext18_PSonly[prim] <- ext_q[(k + 1):(2 * k)]; x2$q_BH_ext18_PSADHD[prim] <- ext_q[(2 * k + 1):(3 * k)]
wcsv(x2, "A_X2_chisep_group_effects.csv")
# BH family definition (observation O1 of VERIFY.md; POST HOC sensitivity). Q25 says "BH-FDR within block" without
# defining the block. The primary family above is X2 only (6 tests; 18 with the simple effects). Alternative: the X2 tests
# pooled with the 5 registered Tier-4 tests run in results_exploratory/tier4_exploratory (X1 x 3, X3, X6; read-only).
T4FDR <- file.path(ROOT, "results_exploratory/tier4_exploratory/fdr_table.csv")
t4 <- utils::read.csv(T4FDR, stringsAsFactors = FALSE)
regX <- t4[t4$family == "registered X", ]
stopifnot(nrow(regX) == 5)
pool_q <- function(p_x2) stats::p.adjust(c(p_x2, regX$p_two), "BH")[seq_along(p_x2)]
q11 <- pool_q(x2$p_two[prim]); q23 <- pool_q(ext_p)
bh_sens <- data.frame(model = x2$model[prim], roi = x2$roi[prim], component = x2$component[prim], p_two = x2$p_two[prim],
                      q_X2_family6 = x2$q_BH_family6[prim], q_X2_plus_regX_11 = q11,
                      q_X2_ext18 = ext_q[1:k], q_X2_ext18_plus_regX_23 = q23[1:k],
                      PSonly_minus_TD_p = x2$PSonly_minus_TD_p[prim], q_PSonly_ext18 = ext_q[(k + 1):(2 * k)],
                      q_PSonly_ext18_plus_regX_23 = q23[(k + 1):(2 * k)],
                      PSADHD_minus_ADHD_p = x2$PSADHD_minus_ADHD_p[prim], q_PSADHD_ext18 = ext_q[(2 * k + 1):(3 * k)],
                      q_PSADHD_ext18_plus_regX_23 = q23[(2 * k + 1):(3 * k)], stringsAsFactors = FALSE)
wcsv(bh_sens, "A_X2_BH_family_sensitivity.csv")
log_msg(sprintf("X2 BH family sensitivity: smallest q  X2-only 6 = %.4f; X2 + registered X (11) = %.4f; ext 18 = %.4f; ext 18 + registered X (23) = %.4f",
                min(bh_sens$q_X2_family6), min(q11), min(ext_q), min(q23)))
# descriptive values of the components at the paper constants
x2_desc <- do.call(rbind, lapply(names(X2_ROIS), function(roi) {
  cv <- X2_ROIS[[roi]][1]; rv <- X2_ROIS[[roi]][2]
  do.call(rbind, lapply(names(DR_SET), function(nm) {
    s <- sep(first, cv, rv, DR_SET[[nm]], 0)
    data.frame(roi = roi, Dr_label = nm, Dr = DR_SET[[nm]], mean_chi = mean(first[[cv]]), mean_R2s = mean(first[[rv]]),
               mean_R2s_over_Dr_ppb = mean(1000 * first[[rv]] / DR_SET[[nm]]),
               mean_para = mean(s$para), sd_para = stats::sd(s$para), mean_dia = mean(s$dia), sd_dia = stats::sd(s$dia),
               n_para_lt0 = sum(s$para < 0), n_dia_gt0 = sum(s$dia > 0), r_para_dia = stats::cor(s$para, s$dia))
  }))
}))
wcsv(x2_desc, "A_X2_component_descriptives.csv")

# ------------------------------------------------------------------------------------------------------------
# B. Joint nuisance model and Dice-quality analyses (POST HOC)
# ------------------------------------------------------------------------------------------------------------
NUIS <- c(dice = "dice_min", motion = "motion_gsr", csf_offset = "ref2_minus_wb_ppb", icv = "icv_seg_ml")
nb <- list()
add_nb <- function(m, label, note = "") {
  r <- ps_row(m, label, "B_nuisance", bound = TRUE, note = note)
  r$pB_frozen_style <- frozen_style_pB(m); nb[[length(nb) + 1]] <<- r; invisible(m)
}
add_nb(m1_ref, "M1 (frozen primary)")
single <- list()
for (k in names(NUIS)) single[[k]] <- add_nb(fit_m1y(first, "chi_SNVTA", paste0(" + scale(", NUIS[[k]], ")"), need = NUIS[[k]]),
                                             paste0("M1 + ", NUIS[[k]], " (alone)"))
joint <- add_nb(fit_m1y(first, "chi_SNVTA", paste0(" + scale(", NUIS, ")", collapse = ""), need = NUIS),
                "M1 + Dice + motion + CSF offset + ICV (joint)")
ver_rows[["B_joint"]] <- export_design(joint, "B_joint")
nbt <- do.call(rbind, nb)
chk("M1 + Dice beta (frozen 3.749)", single$dice$b["PS_e"], 3.749, 5e-4)
chk("M1 + motion beta (frozen 4.227)", single$motion$b["PS_e"], 4.227, 5e-4)
chk("M1 + ICV beta (frozen 5.337)", single$icv$b["PS_e"], 5.337, 5e-4)
chk("M1 + Dice frozen-style pB (0.0019)", nbt$pB_frozen_style[nbt$model == "M1 + dice_min (alone)"], 0.0019, 5e-5)
chk("M1 + motion frozen-style pB (0.0004)", nbt$pB_frozen_style[nbt$model == "M1 + motion_gsr (alone)"], 0.0004, 5e-5)
chk("M1 + ICV frozen-style pB (0.0003)", nbt$pB_frozen_style[nbt$model == "M1 + icv_seg_ml (alone)"], 0.0003, 5e-5)
jc <- coef_table(joint); jc <- jc[grepl("scale\\((dice_min|motion_gsr|ref2_minus_wb_ppb|icv_seg_ml)\\)|^PS_e$|^ADHD_e$", jc$term), ]
jc$ci_lo <- jc$est - stats::qt(.975, joint$df) * jc$se_hc3; jc$ci_hi <- jc$est + stats::qt(.975, joint$df) * jc$se_hc3
wcsv(nbt, "B_nuisance_models.csv"); wcsv(jc, "B_joint_model_coefficients.csv")
# collinearity of the nuisance block with PS (VIF of PS_e in the joint design, via R2 of PS_e on the other columns)
Xj <- stats::model.matrix(joint$fit); jps <- which(colnames(Xj) == "PS_e")
vif_ps <- 1 / (1 - summary(stats::lm(Xj[, jps] ~ Xj[, -c(1, jps)]))$r.squared)
Xm <- stats::model.matrix(core$m$fit); mps <- which(colnames(Xm) == "PS_e")
vif_ps_m1 <- 1 / (1 - summary(stats::lm(Xm[, mps] ~ Xm[, -c(1, mps)]))$r.squared)

# R2* in SN-VTA with and without Dice (and joint)
r2_rows <- list(
  ps_row(fit_m1y(first, "r2s_SNVTA"), "R2* M1", "B_r2s"),
  ps_row(fit_m1y(first, "r2s_SNVTA", " + scale(dice_min)", need = "dice_min"), "R2* M1 + Dice", "B_r2s"),
  ps_row(fit_m1y(first, "r2s_SNVTA", paste0(" + scale(", NUIS, ")", collapse = ""), need = NUIS), "R2* M1 + joint nuisance", "B_r2s"))
r2t <- do.call(rbind, r2_rows)
chk("K2 R2* beta (frozen 1.098)", r2t$beta_PS[1], 1.098, 5e-4)
m_r2d <- fit_m1y(first, "r2s_SNVTA", " + scale(dice_min)", need = "dice_min")
r2t$dice_coef <- c(NA, term_test(m_r2d, "scale(dice_min)")$est, NA)
r2t$dice_coef_p <- c(NA, term_test(m_r2d, "scale(dice_min)")$p_two, NA)
wcsv(r2t, "B_r2s_models.csv")

# Dice tertiles (label-free cut points over the 80 first-usable sessions)
cuts <- stats::quantile(first$dice_min, c(1, 2) / 3, type = 7)
first$dice_t <- factor(cut(first$dice_min, c(-Inf, cuts, Inf), labels = c("T1_low", "T2_mid", "T3_high")))
tert_tab <- as.data.frame.matrix(table(first$dice_t, factor(first$cell, levels = CELLS)))
tert_tab$tertile <- rownames(tert_tab)
tert_tab$dice_range <- tapply(first$dice_min, first$dice_t, function(x) sprintf("%.4f-%.4f", min(x), max(x)))[rownames(tert_tab)]
tert_tab$mean_chi <- tapply(first$chi_SNVTA, first$dice_t, mean)[rownames(tert_tab)]
tert_fisher_cell <- stats::fisher.test(table(first$dice_t, factor(first$cell, levels = CELLS)))$p.value
tert_fisher_PS <- stats::fisher.test(table(first$dice_t, first$PS_e))$p.value
TEMPLATE <- ps_row(m1_ref, "", "B_dice_tertile", bound = TRUE); TEMPLATE[1, ] <- NA; TEMPLATE$block <- "B_dice_tertile"; TEMPLATE$outcome <- "chi_SNVTA"
tr <- list()
for (tl in levels(first$dice_t)) {
  dd <- first[first$dice_t == tl, ]
  m <- try_fit(fit_m1y(dd, "chi_SNVTA"), paste("M1 in", tl))
  if (!is.null(m) && is.finite(m$se["PS_e"])) tr[[length(tr) + 1]] <- ps_row(m, paste("within", tl, ": full M1"), "B_dice_tertile", bound = TRUE)
  else { r0 <- TEMPLATE; r0$model <- paste("within", tl, ": full M1"); r0$n <- nrow(dd); r0$note <- "fit failed or HC3 undefined"
         tr[[length(tr) + 1]] <- r0 }
  dd <- set_race3(dd)
  mr <- fit_ols_hc3(chi_SNVTA ~ PS_e * ADHD_e + age + sex_e, dd); mr$data <- dd; mr$y <- "chi_SNVTA"
  tr[[length(tr) + 1]] <- ps_row(mr, paste("within", tl, ": reduced (PS*ADHD + age + sex)"), "B_dice_tertile", bound = TRUE,
                                 note = "reduced model because n per tertile is about 27")
}
mfe <- fit_m1y(first, "chi_SNVTA", " + dice_t"); tr[[length(tr) + 1]] <- ps_row(mfe, "pooled M1 + Dice-tertile fixed effects", "B_dice_tertile", bound = TRUE)
mfe2 <- fit_m1y(first, "chi_SNVTA", " + scale(dice_min) + dice_t")
tr[[length(tr) + 1]] <- ps_row(mfe2, "pooled M1 + Dice tertile FE + linear Dice", "B_dice_tertile", bound = TRUE)
mqm <- fit_m1y(first[first$dice_t != "T1_low", ], "chi_SNVTA")
tr[[length(tr) + 1]] <- ps_row(mqm, "quality-matched: drop lowest Dice tertile (M1)", "B_dice_tertile", bound = TRUE)
mqm2 <- fit_m1y(first[first$dice_t != "T1_low", ], "chi_SNVTA", " + scale(dice_min)")
tr[[length(tr) + 1]] <- ps_row(mqm2, "quality-matched: drop lowest Dice tertile (M1 + Dice)", "B_dice_tertile", bound = TRUE)
# heterogeneity of beta_PS across tertiles: M1 + dice_t + PS_e:dice_t, HC3 Wald (2 df), F reference
dh <- first; dh$dice_t <- factor(dh$dice_t); stats::contrasts(dh$dice_t) <- stats::contr.sum(3)
mh <- fit_m1y(dh, "chi_SNVTA", " + dice_t + PS_e:dice_t")
jj <- grep("^PS_e:dice_t", names(mh$b)); W <- drop(t(mh$b[jj]) %*% solve(mh$V[jj, jj]) %*% mh$b[jj])
het <- list(terms = names(mh$b)[jj], wald_chi2 = W, df1 = length(jj), F = W / length(jj), df2 = mh$df,
            p_F = stats::pf(W / length(jj), length(jj), mh$df, lower.tail = FALSE))
trt <- do.call(rbind, tr)
wcsv(trt, "B_dice_tertile_models.csv"); wcsv(tert_tab, "B_dice_tertile_cells.csv")
cor_dice_chi <- stats::cor.test(first$dice_min, first$chi_SNVTA, method = "spearman", exact = FALSE)

# ------------------------------------------------------------------------------------------------------------
# C. Low-TD-tail check (POST HOC). Values only; no image is opened.
# ------------------------------------------------------------------------------------------------------------
TAIL <- 35
first$rstud_M1 <- stats::rstudent(core$m$fit)
tail_i <- which(first$cell == "TD" & first$chi_SNVTA < TAIL)
QCV <- c("dice_min", "motion_gsr", "fimp_SNVTA", "inmask_SNVTA", "Q_phase", "roi_vol", "r2s_SNVTA", "ref2_minus_wb_ppb", "icv_seg_ml", "age")
tail_tab <- first[tail_i, c("code", "session_id", "ses_order", "chi_SNVTA", "qc_grade", "dice_min", "motion_gsr", "fimp_SNVTA",
                             "inmask_SNVTA", "Q_phase", "r2s_SNVTA", "roi_vol", "chi_SNVTA_median", "chi_SNonly", "chi_VTAPBP",
                             "ref2_minus_wb_ppb", "age", "rstud_M1")]
tail_tab <- tail_tab[order(tail_tab$chi_SNVTA), ]
td <- first[first$cell == "TD", ]; td_rest <- td[!td$sess_key %in% first$sess_key[tail_i], ]
pct <- function(x, ref) 100 * mean(ref <= x)
for (v in QCV) {
  tail_tab[[paste0("pctTDrest_", v)]] <- vapply(first[[v]][tail_i][order(first$chi_SNVTA[tail_i])], pct, 0, ref = td_rest[[v]])
  tail_tab[[paste0("pctAll80_", v)]] <- vapply(first[[v]][tail_i][order(first$chi_SNVTA[tail_i])], pct, 0, ref = first[[v]])
}
qc_cmp <- do.call(rbind, lapply(QCV, function(v) {
  a <- first[[v]][tail_i]; b <- td_rest[[v]]; w <- if (length(a) >= 1) suppressWarnings(stats::wilcox.test(a, b, exact = TRUE)) else NULL
  data.frame(metric = v, tail_mean = mean(a), tail_min = min(a), tail_max = max(a),
             TDrest_n = length(b), TDrest_median = stats::median(b), TDrest_q25 = stats::quantile(b, .25), TDrest_q75 = stats::quantile(b, .75),
             TDrest_min = min(b), TDrest_max = max(b), all80_median = stats::median(first[[v]]),
             wilcox_p_tail_vs_TDrest = if (is.null(w)) NA else w$p.value, row.names = NULL)
}))
qc_grade_tab <- list(tail = as.list(table(factor(first$qc_grade[tail_i], levels = 0:1))),
                     TD_rest = as.list(table(factor(td_rest$qc_grade, levels = 0:1))),
                     all80 = as.list(table(factor(first$qc_grade, levels = 0:1))))
wcsv(tail_tab, "C_low_TD_tail_sessions.csv"); wcsv(qc_cmp, "C_low_TD_tail_qc_comparison.csv")
cr <- list(ps_row(m1_ref, "M1 (all 80)", "C_tail", bound = TRUE))
m_notail <- fit_m1y(first[-tail_i, ], "chi_SNVTA")
cr[[2]] <- ps_row(m_notail, sprintf("M1 without the %d TD sessions with chi < %g ppb", length(tail_i), TAIL), "C_tail", bound = TRUE)
lf_i <- which(first$chi_SNVTA < TAIL)
m_lf <- fit_m1y(first[-lf_i, ], "chi_SNVTA")
cr[[3]] <- ps_row(m_lf, sprintf("M1 without every first-usable session with chi < %g ppb (label-free; %d sessions)", TAIL, length(lf_i)),
                  "C_tail", bound = TRUE, note = paste("cells of dropped:", paste(table(factor(first$cell[lf_i], levels = CELLS)), collapse = "/")))
crt <- do.call(rbind, cr)
crt$huber_beta <- c(core$huber$est,
                    huber_fit(m_notail$formula, m_notail$data, "PS_e", stats::residuals(m_notail$fit))$est,
                    huber_fit(m_lf$formula, m_lf$data, "PS_e", stats::residuals(m_lf$fit))$est)
ver_rows[["C_notail"]] <- export_design(m_notail, "C_notail")
wcsv(crt, "C_low_TD_tail_refits.csv")

# ------------------------------------------------------------------------------------------------------------
# D. Putamen / GP positive-control diagnostics (POST HOC)
# ------------------------------------------------------------------------------------------------------------
add_derived <- function(dd) {
  dd$chi_PUT_unref <- dd$chi_PUT + dd$ref_ppb                 # before whole-brain referencing
  dd$chi_PUT_csfref <- dd$chi_PUT - dd$ref2_minus_wb_ppb      # relative to lateral-ventricle CSF
  dd$chi_GP_unref <- dd$chi_GP + dd$ref_ppb
  dd$chi_GP_csfref <- dd$chi_GP - dd$ref2_minus_wb_ppb
  s <- sep(dd, "chi_PUT", "r2s_PUT", DR0, 0); dd$chi_para_PUT <- s$para; dd$chi_dia_PUT <- s$dia
  s <- sep(dd, "chi_GP", "r2s_GP", DR0, 0); dd$chi_para_GP <- s$para; dd$chi_dia_GP <- s$dia
  dd
}
first <- add_derived(first); usable <- add_derived(usable)
D_OUT <- c(r2s_PUT = "primary", chi_PUT = "primary", chi_PUT_raw = "primary", r2s_GP = "primary", chi_GP = "primary",
           chi_GP_raw = "co-report", chi_PUT_unref = "co-report (referencing)", chi_PUT_csfref = "co-report (referencing)",
           chi_GP_unref = "co-report (referencing)", chi_GP_csfref = "co-report (referencing)", ref_ppb = "co-report (referencing)",
           ref2_minus_wb_ppb = "co-report (referencing)", chi_para_PUT = "co-report (X2 model-driven)", chi_dia_PUT = "co-report (X2 model-driven)",
           chi_para_GP = "co-report (X2 model-driven)", chi_dia_GP = "co-report (X2 model-driven)", r2s_SNVTA = "co-report", chi_SNVTA = "co-report")
COVS <- c(`age only` = "", `+ sex` = " + sex_e", `+ sex + cell` = " + sex_e + PS_e * ADHD_e")
dr <- list()
for (y in names(D_OUT)) {
  sdy <- stats::sd(first[[y]], na.rm = TRUE)
  for (cn in names(COVS)) {
    dd <- first[!is.na(first[[y]]), ]
    fo <- stats::as.formula(paste(y, "~ age", COVS[[cn]]))
    m <- fit_ols_hc3(fo, dd); tt <- term_test(m, "age")
    dr[[length(dr) + 1]] <- data.frame(outcome = y, role = D_OUT[[y]], sample = "first-usable (OLS, HC3)", covariates = cn, n = m$n,
      n_subjects = length(unique(dd$sub)), slope_per_y = tt$est, se = tt$se, ci_lo = tt$ci95[1], ci_hi = tt$ci95[2], df = tt$df,
      p_two = tt$p_two, p_one_sided_pos = tt$p_upper, slope_per_y_in_SD = tt$est / sdy, stringsAsFactors = FALSE)
    du <- usable[!is.na(usable[[y]]), ]
    # nlme default optimiser (nlminb); the opt = "optim" control was dropped on 2026-10-07 (VERIFY.md sec. 3.1)
    lm_ <- try_fit(nlme::lme(stats::as.formula(paste(y, "~ age", COVS[[cn]])), random = ~ 1 | sub, data = du, method = "REML"),
                   paste("lme", y, cn))
    if (!is.null(lm_)) {
      tb <- summary(lm_)$tTable
      q <- stats::qt(.975, tb["age", "DF"])
      dr[[length(dr) + 1]] <- data.frame(outcome = y, role = D_OUT[[y]], sample = "all usable sessions (lme REML, random intercept)",
        covariates = cn, n = nrow(du), n_subjects = length(unique(du$sub)), slope_per_y = tb["age", "Value"], se = tb["age", "Std.Error"],
        ci_lo = tb["age", "Value"] - q * tb["age", "Std.Error"], ci_hi = tb["age", "Value"] + q * tb["age", "Std.Error"],
        df = tb["age", "DF"], p_two = tb["age", "p-value"], p_one_sided_pos = stats::pt(tb["age", "t-value"], tb["age", "DF"], lower.tail = FALSE),
        slope_per_y_in_SD = tb["age", "Value"] / sdy, stringsAsFactors = FALSE)
    }
  }
  # within- vs between-person age (all usable; + sex + cell)
  du <- usable[!is.na(usable[[y]]), ]
  du$age_b <- stats::ave(du$age, du$sub); du$age_w <- du$age - du$age_b
  lw <- try_fit(nlme::lme(stats::as.formula(paste(y, "~ age_w + age_b + sex_e + PS_e * ADHD_e")), random = ~ 1 | sub, data = du,
                          method = "REML"), paste("lme wb", y))          # nlme default optimiser (see above)
  if (!is.null(lw)) for (tm in c("age_w", "age_b")) {
    tb <- summary(lw)$tTable; q <- stats::qt(.975, tb[tm, "DF"])
    dr[[length(dr) + 1]] <- data.frame(outcome = y, role = D_OUT[[y]], sample = paste("all usable (lme REML),",
      if (tm == "age_w") "WITHIN-person age" else "BETWEEN-person mean age"), covariates = "+ sex + cell", n = nrow(du),
      n_subjects = length(unique(du$sub)), slope_per_y = tb[tm, "Value"], se = tb[tm, "Std.Error"],
      ci_lo = tb[tm, "Value"] - q * tb[tm, "Std.Error"], ci_hi = tb[tm, "Value"] + q * tb[tm, "Std.Error"], df = tb[tm, "DF"],
      p_two = tb[tm, "p-value"], p_one_sided_pos = stats::pt(tb[tm, "t-value"], tb[tm, "DF"], lower.tail = FALSE),
      slope_per_y_in_SD = tb[tm, "Value"] / sdy, stringsAsFactors = FALSE)
  }
}
drt <- do.call(rbind, dr)
pc2x_chk <- drt[drt$sample == "first-usable (OLS, HC3)" & drt$covariates == "age only", ]
chk("PC2x chi_PUT slope (frozen 0.382)", pc2x_chk$slope_per_y[pc2x_chk$outcome == "chi_PUT"], 0.382, 5e-4)
chk("PC2x chi_GP slope (frozen 2.397)", pc2x_chk$slope_per_y[pc2x_chk$outcome == "chi_GP"], 2.397, 5e-4)
# registered expectation (prereg Q21 / PC2w: 2.6-3.4 ppb (putamen) and 4.4-6.9 ppb (GP) over 1.95 y)
EXPECT <- list(chi_PUT = c(2.6, 3.4) / 1.95, chi_GP = c(4.4, 6.9) / 1.95)
drt$expected_lo <- NA; drt$expected_hi <- NA
for (y in names(EXPECT)) for (yy in paste0(y, c("", "_raw", "_unref", "_csfref"))) {
  i <- drt$outcome == yy; drt$expected_lo[i] <- EXPECT[[y]][1]; drt$expected_hi[i] <- EXPECT[[y]][2]
}
drt$ci_excludes_expected_range <- ifelse(is.na(drt$expected_lo), NA, drt$ci_hi < drt$expected_lo | drt$ci_lo > drt$expected_hi)
# chi slope implied if the whole R2* slope were paramagnetic (chi_dia constant): 1000 * slope(R2*) / Dr ppb/y
impl <- do.call(rbind, lapply(c("PUT", "GP"), function(r) {
  a <- drt[drt$outcome == paste0("r2s_", r), c("sample", "covariates", "slope_per_y", "ci_lo", "ci_hi")]
  b <- drt[drt$outcome == paste0("chi_", r), c("sample", "covariates", "slope_per_y", "ci_lo", "ci_hi")]
  mm <- merge(a, b, by = c("sample", "covariates"), suffixes = c("_r2s", "_chi"))
  data.frame(roi = r, mm, implied_chi_slope_ppb_per_y = 1000 * mm$slope_per_y_r2s / DR0,
             implied_lo = 1000 * mm$ci_lo_r2s / DR0, implied_hi = 1000 * mm$ci_hi_r2s / DR0)
}))
wcsv(drt, "D_age_slopes.csv"); wcsv(impl, "D_implied_chi_slope_from_R2s.csv")
# putamen "shortfall" estimate (age slope of chi_dia_PUT = (chi_PUT - 1000 * r2s_PUT / Dr) / 2) by Dr, sample and
# covariates (review item W11; POST HOC). Shows how far the estimate depends on Dr and on the covariate set.
sf <- list()
for (Dr in c(137, DR0, 321, 548)) for (smp in c("first", "usable")) for (cn in names(COVS)) {
  dd <- if (smp == "first") first else usable
  dd$chi_dia_PUT_Dr <- (dd$chi_PUT - 1000 * dd$r2s_PUT / Dr) / 2
  dd <- dd[!is.na(dd$chi_dia_PUT_Dr), ]
  fo <- stats::as.formula(paste("chi_dia_PUT_Dr ~ age", COVS[[cn]]))
  if (smp == "first") {
    tt <- term_test(fit_ols_hc3(fo, dd), "age")
    est <- tt$est; se <- tt$se; lo <- tt$ci95[1]; hi <- tt$ci95[2]; dfa <- tt$df; p <- tt$p_two
    lab <- "first-usable (OLS, HC3)"
  } else {
    ml <- try_fit(nlme::lme(fo, random = ~ 1 | sub, data = dd, method = "REML"), paste("lme shortfall", Dr, cn))
    if (is.null(ml)) next
    tb <- summary(ml)$tTable; q <- stats::qt(.975, tb["age", "DF"])
    est <- tb["age", "Value"]; se <- tb["age", "Std.Error"]; lo <- est - q * se; hi <- est + q * se; dfa <- tb["age", "DF"]
    p <- tb["age", "p-value"]; lab <- "all usable sessions (lme REML, random intercept)"
  }
  sf[[length(sf) + 1]] <- data.frame(Dr = Dr, sample = lab, covariates = cn, n = nrow(dd), slope = est, se = se, ci_lo = lo,
                                     ci_hi = hi, df = dfa, p_two = p, stringsAsFactors = FALSE)
}
sft <- do.call(rbind, sf)
stopifnot(nrow(sft) == 24)
wcsv(sft, "D_putamen_shortfall_by_Dr.csv")
cor_put <- list(r_chi_r2s_first = stats::cor(first$chi_PUT, first$r2s_PUT),
                r_chi_r2s_first_age_partial = stats::cor(stats::residuals(stats::lm(chi_PUT ~ age, first)), stats::residuals(stats::lm(r2s_PUT ~ age, first))),
                r_chi_r2s_GP_first = stats::cor(first$chi_GP, first$r2s_GP),
                age_range_first = range(first$age), age_range_usable = range(usable$age),
                n_subjects_with_2plus_usable = sum(table(usable$sub) >= 2))

# ------------------------------------------------------------------------------------------------------------
# verification files and results.json
# ------------------------------------------------------------------------------------------------------------
wcsv(do.call(rbind, ver_rows), "verify/R_reference_estimates.csv")
df2list <- function(x) lapply(split(x, seq_len(nrow(x))), as.list)
res <- list(
  meta = list(label = "POST HOC / registered-exploratory (X2 only); never a claim; confirmatory outcome B is not re-decided",
              run_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), R = R.version.string,
              input = INPUT, input_sha256 = sha256_file(INPUT), gates_sha256 = sha256_file(GATES),
              script = SELF, script_sha256 = if (file.exists(SELF)) sha256_file(SELF) else NA,
              frozen_helpers_sha256 = setNames(lapply(c("constants.R", "utils.R", "design.R", "plausibility.R", "m1.R"),
                                                      function(f) sha256_file(file.path(CODE, "R", f))),
                                               c("constants.R", "utils.R", "design.R", "plausibility.R", "m1.R")),
              samples = list(first_usable = as.list(cell_counts(first)), usable_sessions = nrow(usable),
                             usable_subjects = length(unique(usable$sub))),
              bound_ppb = BOUND, rL = rL,
              revision = "2026-10-07: fixes from verification/VERIFY.md applied (CHANGES_2026-10-07.md); still POST HOC, nothing is a finding",
              lme_optimizer = "nlme default (nlminb); the opt = 'optim' control was dropped (VERIFY.md sec. 3.1)",
              declared_deviations = list(
                X2_constants = paste("Q25 required the X2 constants in the amendment-2 module before any chi existed; they were not there.",
                                     "Dr = 274 Hz/ppm and R2 = 0 were fixed after unlock from Dimov et al. 2022, not tuned; Dr range reported."),
                R2s_voxel_set = paste("Prereg (LEADchi_prereg_v1.0_FINAL.md line 253) defines R2*_ROI on the same voxel sets as the matching chi ROI,",
                                      "restricted to Mask 3. extract_rois.py (lines 289-295) uses Mask 3 & finite R2* within the ROI, without the",
                                      "chi support S or the implausible-voxel exclusion. Departure from the plan text; also affects the frozen K2 and PC3."),
                BH_family = paste("Q25 says 'BH-FDR within block' without defining the block. X2-only families (6 and 18 tests) are used here;",
                                  "pooling with the 5 registered Tier-4 tests is reported as a sensitivity (A_X2_BH_family_sensitivity.csv)."))),
  checks = checks,
  A_X2 = list(constants = list(Dr_Hz_per_ppm = DR0, Dr_theoretical = 321, R2_assumed = 0, Dr_sensitivity = as.list(DR_SET),
                               R2_illustrative = R2_SET[-1], source = "Dimov et al. 2022 Tomography 8:1544-1551, sec. 2.1 Eq. 1 (page within the article not verified)"),
              identity_max_abs_err = ident_err, R2_cancellation_max_beta_range = canc_err, R2_cancellation_max_se_range = canc_se_err,
              group_effects = df2list(x2), component_descriptives = df2list(x2_desc),
              bh_family_sensitivity = df2list(bh_sens), tier4_fdr_table = T4FDR, tier4_fdr_table_sha256 = sha256_file(T4FDR),
              tier4_registered_X_p = setNames(as.list(regX$p_two), regX$test)),
  B_nuisance = list(models = df2list(nbt), joint_coefficients = df2list(jc), vif_PS_joint = vif_ps, vif_PS_M1 = vif_ps_m1,
                    r2s = df2list(r2t), dice_tertile_cuts = unname(cuts), dice_tertile_cells = df2list(tert_tab),
                    dice_tertile_fisher_p_cell = tert_fisher_cell, dice_tertile_fisher_p_PS = tert_fisher_PS,
                    dice_tertile_models = df2list(trt), dice_heterogeneity = het,
                    spearman_dice_chi = list(rho = unname(cor_dice_chi$estimate), p = cor_dice_chi$p.value)),
  C_low_TD_tail = list(threshold_ppb = TAIL, n_tail = length(tail_i), sessions = df2list(tail_tab), qc_comparison = df2list(qc_cmp),
                       qc_grade = qc_grade_tab, refits = df2list(crt)),
  D_positive_control = list(expected_ppb_per_y = EXPECT, slopes = df2list(drt), implied_chi_from_R2s = df2list(impl), correlations = cor_put,
                            putamen_shortfall_by_Dr = df2list(sft)))
jsonlite::write_json(res, file.path(OUT, "results.json"), auto_unbox = TRUE, digits = 10, pretty = TRUE, na = "null")
writeLines(.LOG$lines, file.path(OUT, "run.log"))
utils::capture.output(sessionInfo(), file = file.path(OUT, "sessionInfo.txt"))
log_msg("done; all checks pass: ", all(vapply(checks, `[[`, TRUE, "pass")))
