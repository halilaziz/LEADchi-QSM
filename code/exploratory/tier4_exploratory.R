#!/usr/bin/env Rscript
# LEAD-chi -- TASK T4: Tier-4 EXPLORATORY modules X1, X3, X6 (prereg Q25 / sec. 4.5 tier 4) and descriptive
# direction diagnostics of the unexpected positive beta_PS. EVERYTHING here is EXPLORATORY / post hoc
# (Annex C.4): hypothesis-generating only, BH-FDR within block, never a claim. The confirmatory outcome is frozen
# in mirror/results/confirmatory/report.md (H1 final outcome B) and is NOT re-decided here.
#
# Conventions: the frozen analysis helpers are sourced READ-ONLY from mirror/scripts/analysis/R (same effect coding
# PS_e / ADHD_e / sex_e = +/-0.5, race3 contr.sum over the levels present, ns(age, 3), scale(roi_vol), scale(Q_phase),
# OLS + HC3, t reference with df = n - p, the same plausibility edits and the same first-usable sample).
#
#   Rscript tier4_exploratory.R          (paths are fixed below; no server access; nothing is written outside OUT)

suppressPackageStartupMessages({
  library(splines); library(sandwich); library(MASS); library(lme4); library(nlme); library(boot); library(jsonlite)
  library(digest); library(ggplot2)
})
Sys.setenv(LANGUAGE = "en")
ROOT    <- "/path/to/Penn_LEAD/LEADchi"
CODE    <- file.path(ROOT, "mirror/scripts/analysis")
INPUT   <- file.path(ROOT, "mirror/unlocked/merged_qsmxt_primary_SNVTA.csv")
GATES   <- file.path(ROOT, "mirror/gates/gates_preunlock.json")
QC_IDS  <- file.path(ROOT, "mirror/qc/qc_ids.csv")
QC_P1   <- file.path(ROOT, "mirror/qc/rating_pass1.csv")
OUT     <- file.path(ROOT, "results_exploratory/tier4_exploratory")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
for (f in c("constants.R", "utils.R", "design.R", "plausibility.R", "m1.R", "t1.R")) source(file.path(CODE, "R", f))
LABEL <- "EXPLORATORY (Tier 4 / post hoc; hypothesis-generating only; never a claim)"
log_msg("T4 tier-4 exploratory start; input ", INPUT)

# ------------------------------------------------------------------------------------------------------------
# 0. Data: identical preparation to run_all.R (validate -> code -> plausibility -> samples)
# ------------------------------------------------------------------------------------------------------------
d0 <- utils::read.csv(INPUT, stringsAsFactors = FALSE, na.strings = c("NA", ""))
schema <- read_schema(file.path(CODE, "schema/roi_table_schema.csv"))
validate_table(d0, schema)
d0 <- code_vars(d0)
gates <- jsonlite::read_json(GATES, simplifyVector = TRUE)
pl <- apply_plausibility(d0); d <- pl$d
# forced-review flag (automated plausibility / slab-coverage trigger used by the QC rater; amendment 1c) joined
# through the coded QC-ID map: code + ses_order -> qc_id -> rating sheet (pass 1 carries all 181 sessions)
qids <- utils::read.csv(QC_IDS, stringsAsFactors = FALSE); p1 <- utils::read.csv(QC_P1, stringsAsFactors = FALSE)
qids$forced_review <- p1$forced_review[match(qids$qc_id, p1$qc_id)]
d$forced_review <- qids$forced_review[match(paste(d$code, d$ses_order), paste(qids$code, qids$ses_order))]
S <- build_samples(d)
first <- set_race3(S$first)
first$cell <- factor(first$cell, levels = CELLS)
stopifnot(nrow(first) == 80)
log_msg("first-usable ", cell_str(first), "; forced_review joined for ", sum(!is.na(first$forced_review)), " of ", nrow(first))
CELL_COL <- c(TD = "#2a78d6", `ADHD-only` = "#eb6834", `PS-only` = "#1baf7a", `PS+ADHD` = "#eda100")

# ------------------------------------------------------------------------------------------------------------
# helpers
# ------------------------------------------------------------------------------------------------------------
ci_str <- function(est, lo, hi, d = 2) sprintf("%s [%s, %s]", f3(est, d), f3(lo, d), f3(hi, d))
# generic M1-form fit for any outcome y (registered RHS; `vol` lets a matched ROI volume replace roi_vol)
m1_rhs_y <- function(dd, vol = "roi_vol", extra = "")
  paste0("PS_e * ADHD_e + ns(age, df = 3) + sex_e", race_term(dd), " + scale(", vol, ") + scale(Q_phase)", extra)
fit_y <- function(dd, y, vol = "roi_vol", extra = "") {
  dd <- dd[!is.na(dd[[y]]) & !is.na(dd[[vol]]), ]
  dd <- set_race3(dd)
  fo <- stats::as.formula(paste(y, "~", m1_rhs_y(dd, vol, extra)))
  m <- fit_ols_hc3(fo, dd); m$data <- dd; m$y <- y; m
}
# one row of group-effect summaries from an M1-form fit (ppb and standardised units)
ps_row <- function(m, label, block, module, note = "") {
  tt <- term_test(m, "PS_e"); ta <- term_test(m, "ADHD_e"); ti <- term_test(m, "PS_e:ADHD_e")
  s1 <- lin_test(m, c(PS_e = 1, `PS_e:ADHD_e` = -0.5), "PS-only - TD")
  s2 <- lin_test(m, c(PS_e = 1, `PS_e:ADHD_e` = 0.5), "PS+ADHD - ADHD-only")
  sdy <- stats::sd(m$data[[m$y]])
  data.frame(module = module, block = block, test = label, outcome = m$y, n = m$n,
             cells = paste(cell_counts(m$data)[CELLS], collapse = "/"),
             beta_PS = tt$est, se_hc3 = tt$se, ci_lo = tt$ci95[1], ci_hi = tt$ci95[2], t = tt$t, df = tt$df, p_two = tt$p_two,
             delta_resid = tt$est / m$sigma, beta_per_sd_raw = tt$est / sdy, sd_raw = sdy, sigma_resid = m$sigma,
             PSonly_minus_TD = s1$est, PSonly_minus_TD_lo = s1$ci95[1], PSonly_minus_TD_hi = s1$ci95[2], PSonly_minus_TD_p = s1$p_two,
             PSADHD_minus_ADHD = s2$est, PSADHD_minus_ADHD_lo = s2$ci95[1], PSADHD_minus_ADHD_hi = s2$ci95[2], PSADHD_minus_ADHD_p = s2$p_two,
             beta_ADHD = ta$est, beta_ADHD_lo = ta$ci95[1], beta_ADHD_hi = ta$ci95[2], p_ADHD = ta$p_two,
             beta_int = ti$est, beta_int_lo = ti$ci95[1], beta_int_hi = ti$ci95[2], p_int = ti$p_two,
             formula = gsub("\\s+", " ", paste(deparse(m$formula), collapse = "")), note = note, stringsAsFactors = FALSE)
}
# covariate-adjusted cell means: model-matrix column means (proportional weights) with the cell's design values
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
# covariate-adjusted PS main effect of a nuisance / QC variable (no roi_vol / Q_phase on the RHS; HC3)
nuis_ps <- function(dd, v, label, block = "diagnostic: artefact channels", module = "D5") {
  dd <- set_race3(dd[!is.na(dd[[v]]), ])
  fo <- stats::as.formula(paste0(v, " ~ PS_e * ADHD_e + ns(age, df = 3) + sex_e", race_term(dd)))
  m <- fit_ols_hc3(fo, dd); tt <- term_test(m, "PS_e")
  kw <- stats::kruskal.test(dd[[v]], factor(dd$cell, levels = CELLS))$p.value
  rho <- suppressWarnings(stats::cor.test(dd[[v]], dd$chi_SNVTA, method = "spearman"))
  data.frame(module = module, block = block, test = label, outcome = v, n = m$n, beta_PS = tt$est, se_hc3 = tt$se,
             ci_lo = tt$ci95[1], ci_hi = tt$ci95[2], p_two = tt$p_two, beta_per_sd_raw = tt$est / stats::sd(dd[[v]]),
             kruskal_p_4cells = kw, spearman_rho_with_chi_SNVTA = unname(rho$estimate), spearman_p = rho$p.value,
             stringsAsFactors = FALSE)
}
cell_plot <- function(dd, v, ylab, title, adj = NULL, file, h = 4.2, w = 6.5) {
  dd <- dd[!is.na(dd[[v]]), ]; ns <- table(dd$cell)
  lab <- paste0(names(ns), "\n(n = ", as.integer(ns), ")")
  p <- ggplot(dd, aes(x = cell, y = .data[[v]], colour = cell)) +
    geom_boxplot(outlier.shape = NA, width = 0.5, linewidth = 0.4, fill = NA) +
    geom_jitter(width = 0.14, height = 0, size = 1.7, alpha = 0.75, shape = 16) +
    scale_colour_manual(values = CELL_COL, guide = "none") + scale_x_discrete(labels = lab) +
    labs(x = NULL, y = ylab, title = title, subtitle = LABEL) +
    theme_minimal(base_size = 11) + theme(panel.grid.major.x = element_blank(), plot.subtitle = element_text(size = 7.5, colour = "grey35"))
  if (!is.null(adj)) p <- p + geom_pointrange(data = adj, aes(x = cell, y = adj_mean, ymin = lo, ymax = hi), inherit.aes = FALSE,
                                              colour = "black", shape = 18, size = 0.55, linewidth = 0.6,
                                              position = position_nudge(x = 0.36))
  ggsave(file.path(OUT, file), p, width = w, height = h, dpi = 160, bg = "white"); invisible(p)
}

RES <- list(label = LABEL, code = "tier4_exploratory.R", run_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
            input = INPUT, input_sha256 = sha256_file(INPUT), gates_sha256 = sha256_file(GATES),
            confirmatory_frozen = "mirror/results/confirmatory/report.md: H1 final outcome B (beta_PS = +4.963 ppb, p .1169; p_B .0004); not re-decided here",
            sample = list(n = nrow(first), cells = as.list(cell_counts(first))))
tests <- list()   # rows for the BH-FDR table
add_test <- function(x) tests[[length(tests) + 1]] <<- x

# ------------------------------------------------------------------------------------------------------------
# 1. X1: dentate (R, L; G5-DN pass only = non-missing chi_DN) and red-nucleus group effects, registered M1 form
# ------------------------------------------------------------------------------------------------------------
log_msg("X1")
x1 <- list()
mDNR <- fit_y(first, "chi_DN_R"); mDNL <- fit_y(first, "chi_DN_L"); mRN <- fit_y(first, "chi_RN")
add_test(ps_row(mDNR, "X1 right dentate chi: beta_PS (M1 form)", "registered X (X1/X3/X6)", "X1"))
add_test(ps_row(mDNL, "X1 left dentate chi: beta_PS (M1 form)", "registered X (X1/X3/X6)", "X1"))
add_test(ps_row(mRN,  "X1 red nucleus chi: beta_PS (M1 form)",  "registered X (X1/X3/X6)", "X1"))
# co-reports (not in the FDR block): matched dentate volume; bilateral dentate mean; SN-VTA restricted to the DN subsample
mDNRv <- fit_y(first, "chi_DN_R", vol = "DN_R_vol")
first$chi_DN_bil <- rowMeans(first[, c("chi_DN_R", "chi_DN_L")])
mDNb <- fit_y(first, "chi_DN_bil")
dn_sub <- first[!is.na(first$chi_DN_R), ]
mSN_dnsub <- fit_y(dn_sub, "chi_SNVTA")
x1$coreports <- rbind(ps_row(mDNRv, "X1 co-report: right dentate with scale(DN_R_vol) in place of scale(roi_vol)", "co-report", "X1"),
                      ps_row(mDNb,  "X1 co-report: bilateral dentate mean", "co-report", "X1"),
                      ps_row(mSN_dnsub, "X1 co-report: SN-VTA chi, M1 refitted in the G5-DN-pass subsample", "co-report", "X1",
                             note = "descriptive comparison only; the H1 estimate is frozen"))
x1$raw_by_cell <- do.call(rbind, lapply(c("chi_DN_R", "chi_DN_L", "chi_RN", "r2s_DN_R", "r2s_RN"), function(v) raw_by_cell(first, v)))
x1$adj_cell_means <- rbind(cbind(outcome = "chi_DN_R", adj_cell_means(mDNR)), cbind(outcome = "chi_DN_L", adj_cell_means(mDNL)),
                           cbind(outcome = "chi_RN", adj_cell_means(mRN)))
x1$dn_availability_by_cell <- as.data.frame.matrix(table(first$cell, ifelse(is.na(first$chi_DN_R), "DN missing (G5-DN fail)", "DN present")))
x1$dn_availability_fisher_p <- stats::fisher.test(table(first$cell, is.na(first$chi_DN_R)))$p.value
x1$coefficients <- list(chi_DN_R = coef_table(mDNR), chi_DN_L = coef_table(mDNL), chi_RN = coef_table(mRN))
# R2* of the same ROIs (iron vs myelin reading for X1; descriptive)
x1$r2s <- rbind(ps_row(fit_y(first, "r2s_DN_R"), "X1 descriptive: right dentate R2*", "co-report", "X1"),
                ps_row(fit_y(first, "r2s_RN"), "X1 descriptive: red nucleus R2*", "co-report", "X1"))
cell_plot(first, "chi_DN_R", "right dentate chi (ppb)", "X1: right dentate chi by cell (G5-DN pass)", adj_cell_means(mDNR), "fig_X1_dentate_R_by_cell.png")
cell_plot(first, "chi_RN", "red nucleus chi (ppb)", "X1: red nucleus chi by cell", adj_cell_means(mRN), "fig_X1_red_nucleus_by_cell.png")
RES$X1 <- x1

# ------------------------------------------------------------------------------------------------------------
# 2. X3: sex x PS interaction on SN-VTA chi (M1 + PS_e:sex_e)
# ------------------------------------------------------------------------------------------------------------
log_msg("X3")
mX3 <- fit_y(first, "chi_SNVTA", extra = " + PS_e:sex_e")
ti <- term_test(mX3, "PS_e:sex_e")
ps_girls <- lin_test(mX3, c(PS_e = 1, `PS_e:sex_e` = -0.5), "PS effect in girls (sex_e = -0.5)")
ps_boys  <- lin_test(mX3, c(PS_e = 1, `PS_e:sex_e` = 0.5), "PS effect in boys (sex_e = +0.5)")
add_test(data.frame(module = "X3", block = "registered X (X1/X3/X6)", test = "X3 PS x sex interaction on SN-VTA chi (M1 + PS_e:sex_e)",
                    outcome = "chi_SNVTA", n = mX3$n, beta_PS = ti$est, se_hc3 = ti$se, ci_lo = ti$ci95[1], ci_hi = ti$ci95[2],
                    p_two = ti$p_two, beta_per_sd_raw = ti$est / stats::sd(first$chi_SNVTA),
                    note = "beta = difference in the PS effect, boys minus girls (ppb)", stringsAsFactors = FALSE))
# co-report: also sex x ADHD, so the PS x sex term is not carrying a sex x ADHD difference
mX3b <- fit_y(first, "chi_SNVTA", extra = " + PS_e:sex_e + ADHD_e:sex_e")
tib <- term_test(mX3b, "PS_e:sex_e")
x3 <- list(interaction = ti[c("est", "se", "t", "df", "p_two", "ci95")],
           simple_effects = list(girls = ps_girls[c("est", "se", "p_two", "ci95")], boys = ps_boys[c("est", "se", "p_two", "ci95")]),
           with_ADHD_x_sex = list(interaction = tib[c("est", "se", "p_two", "ci95")],
                                  girls = lin_test(mX3b, c(PS_e = 1, `PS_e:sex_e` = -0.5), "girls")[c("est", "se", "p_two", "ci95")],
                                  boys = lin_test(mX3b, c(PS_e = 1, `PS_e:sex_e` = 0.5), "boys")[c("est", "se", "p_two", "ci95")]),
           n_by_cell_and_sex = as.data.frame.matrix(table(first$cell, first$sex)),
           raw_by_cell_and_sex = do.call(rbind, lapply(c("F", "M"), function(s) cbind(sex = s, raw_by_cell(first[first$sex == s, ], "chi_SNVTA")))),
           coefficients = coef_table(mX3))
first$cell_sex <- interaction(first$cell, first$sex, sep = " / ")
pX3 <- ggplot(first, aes(x = cell, y = chi_SNVTA, colour = cell, shape = sex)) +
  geom_boxplot(aes(group = cell_sex), outlier.shape = NA, width = 0.6, linewidth = 0.4, fill = NA, position = position_dodge(0.7)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.7, seed = 1), size = 1.7, alpha = 0.8) +
  scale_colour_manual(values = CELL_COL, guide = "none") + scale_shape_manual(values = c(F = 16, M = 17), name = "sex") +
  labs(x = NULL, y = "SN-VTA chi (ppb)", title = "X3: SN-VTA chi by cell and sex",
       subtitle = paste0(LABEL, "\nPS x sex interaction = ", ci_str(ti$est, ti$ci95[1], ti$ci95[2]), " ppb (boys - girls), p = ", fp(ti$p_two))) +
  theme_minimal(base_size = 11) + theme(panel.grid.major.x = element_blank(), plot.subtitle = element_text(size = 7.5, colour = "grey35"), legend.position = "top")
ggsave(file.path(OUT, "fig_X3_sex_by_PS.png"), pX3, width = 7, height = 4.4, dpi = 160, bg = "white")
RES$X3 <- x3

# ------------------------------------------------------------------------------------------------------------
# 3. X6: RSAS social anhedonia (cross-sectional, same session) vs SN-VTA chi
#    primary: S2 form with RSAS in place of PRIME, without the BDI/ARI terms (not registered for X6);
#    co-reports: + PS_e; + BDI + ARI; within PS; rank-inverse-normal; Spearman; all-session between/within
# ------------------------------------------------------------------------------------------------------------
log_msg("X6")
x6 <- list()
d6 <- first[!is.na(first$RSAS_total), ]
x6$n_by_cell <- as.list(cell_counts(d6)); x6$rsas_by_cell <- raw_by_cell(d6, "RSAS_total")
x6$rsas_missing_by_cell <- as.data.frame.matrix(table(first$cell, ifelse(is.na(first$RSAS_total), "RSAS missing", "RSAS present")))
x6_rhs <- function(dd, extra = "") paste0("scale(RSAS_total) + ADHD_e + ns(age, df = 3) + sex_e", race_term(dd), " + scale(roi_vol) + scale(Q_phase)", extra)
d6 <- set_race3(d6)
m6 <- fit_ols_hc3(stats::as.formula(paste("chi_SNVTA ~", x6_rhs(d6))), d6); t6 <- term_test(m6, "scale(RSAS_total)")
add_test(data.frame(module = "X6", block = "registered X (X1/X3/X6)", test = "X6 RSAS social anhedonia (per SD) vs SN-VTA chi, cross-sectional",
                    outcome = "chi_SNVTA", n = m6$n, beta_PS = t6$est, se_hc3 = t6$se, ci_lo = t6$ci95[1], ci_hi = t6$ci95[2],
                    p_two = t6$p_two, beta_per_sd_raw = t6$est / stats::sd(d6$chi_SNVTA),
                    note = "beta = ppb per SD of RSAS (ADHD_e, age, sex, race, roi_vol, Q_phase adjusted)", stringsAsFactors = FALSE))
m6p <- fit_ols_hc3(stats::as.formula(paste("chi_SNVTA ~", x6_rhs(d6, " + PS_e"))), d6)
d6b <- set_race3(d6[!is.na(d6$BDI_total) & !is.na(d6$ARI_total), ])
m6b <- fit_ols_hc3(stats::as.formula(paste("chi_SNVTA ~", x6_rhs(d6b, " + scale(BDI_total) + scale(ARI_total)"))), d6b)
d6w <- set_race3(d6[d6$PS_e > 0, ])
m6w <- try_fit(fit_ols_hc3(stats::as.formula(paste("chi_SNVTA ~", x6_rhs(d6w))), d6w), "X6 within PS")
d6$RSAS_rin <- stats::qnorm((rank(d6$RSAS_total) - 0.5) / nrow(d6))
m6r <- fit_ols_hc3(stats::as.formula(paste("chi_SNVTA ~", sub("scale(RSAS_total)", "RSAS_rin", x6_rhs(d6), fixed = TRUE))), d6)
sp <- suppressWarnings(stats::cor.test(d6$RSAS_total, d6$chi_SNVTA, method = "spearman"))
x6$primary <- c(t6[c("est", "se", "t", "df", "p_two", "ci95")], list(n = m6$n, sd_RSAS = stats::sd(d6$RSAS_total)))
x6$coreports <- list(
  plus_PS_e = c(term_test(m6p, "scale(RSAS_total)")[c("est", "se", "p_two", "ci95")], list(beta_PS_in_this_model = term_test(m6p, "PS_e")[c("est", "se", "p_two")])),
  plus_BDI_ARI = c(term_test(m6b, "scale(RSAS_total)")[c("est", "se", "p_two", "ci95")], list(n = m6b$n)),
  within_PS = if (is.null(m6w)) "fit failed" else c(term_test(m6w, "scale(RSAS_total)")[c("est", "se", "p_two", "ci95")], list(n = m6w$n)),
  rank_inverse_normal = term_test(m6r, "RSAS_rin")[c("est", "se", "p_two", "ci95")],
  spearman_raw = list(rho = unname(sp$estimate), p = sp$p.value, n = nrow(d6)))
# all-session between / within-person RSAS (REML maximum rule; descriptive)
da <- S$usable[!is.na(S$usable$RSAS_total), ]
if (nrow(da) >= 40 && length(unique(da$sub)) >= 30) {
  da <- set_race3(da); da$RSAS_bar <- stats::ave(da$RSAS_total, da$sub); da$RSAS_w <- da$RSAS_total - da$RSAS_bar
  da$age_bar <- stats::ave(da$age, da$sub); da$age_w <- da$age - da$age_bar
  sdr <- stats::sd(da$RSAS_total); da$RSAS_bar_z <- (da$RSAS_bar - mean(da$RSAS_bar)) / sdr; da$RSAS_w_z <- da$RSAS_w / sdr
  fx <- stats::as.formula(paste0("chi_SNVTA ~ RSAS_bar_z + RSAS_w_z + ADHD_e + ns(age_bar, 3) + age_w + sex_e", race_term(da), " + scale(roi_vol) + scale(Q_phase)"))
  rm6 <- try_fit(reml_max(fx, da, profile = FALSE), "X6 all-session")
  x6$all_session <- if (is.null(rm6)) "fit failed" else list(
    n_sessions = nrow(da), n_subjects = length(unique(da$sub)), reml_best = rm6$best, icc = rm6$icc,
    between = list(est = unname(rm6$b["RSAS_bar_z"]), se = unname(rm6$se["RSAS_bar_z"]), df = unname(rm6$containment_df["RSAS_bar_z"]),
                   p = 2 * stats::pt(-abs(rm6$b["RSAS_bar_z"] / rm6$se["RSAS_bar_z"]), rm6$containment_df["RSAS_bar_z"])),
    within = list(est = unname(rm6$b["RSAS_w_z"]), se = unname(rm6$se["RSAS_w_z"]), df = unname(rm6$containment_df["RSAS_w_z"]),
                  p = 2 * stats::pt(-abs(rm6$b["RSAS_w_z"] / rm6$se["RSAS_w_z"]), rm6$containment_df["RSAS_w_z"])))
} else x6$all_session <- "too few sessions"
x6$coefficients <- coef_table(m6)
pX6 <- ggplot(d6, aes(x = RSAS_total, y = chi_SNVTA, colour = cell)) + geom_point(size = 1.9, alpha = 0.8, shape = 16) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, se = TRUE, colour = "grey25", linewidth = 0.7, fill = "grey80") +
  scale_colour_manual(values = CELL_COL, name = NULL) +
  labs(x = "RSAS total (same session)", y = "SN-VTA chi (ppb)", title = paste0("X6: RSAS vs SN-VTA chi (n = ", nrow(d6), ")"),
       subtitle = paste0(LABEL, "\nadjusted slope ", ci_str(t6$est, t6$ci95[1], t6$ci95[2]), " ppb per SD RSAS, p = ", fp(t6$p_two), "; Spearman rho = ", f3(sp$estimate, 2), " (raw)")) +
  theme_minimal(base_size = 11) + theme(plot.subtitle = element_text(size = 7.5, colour = "grey35"), legend.position = "top")
ggsave(file.path(OUT, "fig_X6_RSAS_vs_chi.png"), pX6, width = 7, height = 4.6, dpi = 160, bg = "white")
RES$X6 <- x6

# ------------------------------------------------------------------------------------------------------------
# 4. Direction diagnostics (descriptive, post hoc)
# ------------------------------------------------------------------------------------------------------------
log_msg("direction diagnostics")
DG <- list()
# 4a. raw and adjusted SN-VTA chi by cell
mSN <- fit_y(first, "chi_SNVTA")
DG$snvta_raw_by_cell <- raw_by_cell(first, "chi_SNVTA"); DG$snvta_adj_cell_means <- adj_cell_means(mSN)
DG$snvta_m1_check <- list(beta_PS = unname(mSN$b["PS_e"]), se = unname(mSN$se["PS_e"]),
                          note = "must equal the frozen confirmatory M1 (4.963, 3.125); printed as a reproduction check only")
cell_plot(first, "chi_SNVTA", "SN-VTA chi (ppb)", "SN-VTA chi by cell: raw (box/points); M1-adjusted mean, 95% CI (black)", adj_cell_means(mSN), "fig_D1_chi_SNVTA_by_cell.png")

# 4b. regional: is the PS elevation SN-specific or global?
rois <- c(chi_SNVTA = "SN-VTA (primary)", chi_SNonly = "SN only (SNc+SNr)", chi_VTAPBP = "VTA + PBP", chi_GP = "globus pallidus",
          chi_PUT = "putamen", chi_CAU = "caudate", chi_STR = "striatum (CAU+PUT+NAcc)", chi_THAL = "thalamus", chi_RN = "red nucleus",
          chi_STN = "subthalamic nucleus", chi_DN_R = "right dentate", chi_DN_L = "left dentate", chi_DGM = "deep GM composite (CAU+GP)")
DG$regional_raw_by_cell <- do.call(rbind, lapply(names(rois), function(v) raw_by_cell(first, v)))
reg_rows <- do.call(rbind, lapply(names(rois), function(v) ps_row(fit_y(first, v), paste0("regional PS effect: ", rois[[v]]), "diagnostic: regional", "D2")))
reg_rows$roi_label <- rois[reg_rows$outcome]
DG$regional_ps_effects <- reg_rows
# post hoc diagnostics enter the FDR table (regional block), except the primary ROI (frozen) and the X1 ROIs (already in the X block)
for (i in which(!reg_rows$outcome %in% c("chi_SNVTA", "chi_DN_R", "chi_DN_L", "chi_RN"))) add_test(reg_rows[i, ])
# global-iron adjustment: SN-VTA PS effect with the deep-GM composite / GP as a covariate; and a 6-ROI global z composite
mSN_dgm <- fit_y(first, "chi_SNVTA", extra = " + scale(chi_DGM)"); mSN_gp <- fit_y(first, "chi_SNVTA", extra = " + scale(chi_GP)")
zs <- scale(first[, c("chi_GP", "chi_PUT", "chi_CAU", "chi_THAL", "chi_RN", "chi_STN")]); first$chi_global_z <- rowMeans(zs)
mGlob <- fit_y(first, "chi_global_z")
DG$specificity <- list(
  SNVTA_adjusted_for_DGM = c(term_test(mSN_dgm, "PS_e")[c("est", "se", "p_two", "ci95")], list(DGM_slope = term_test(mSN_dgm, "scale(chi_DGM)")[c("est", "p_two")])),
  SNVTA_adjusted_for_GP = c(term_test(mSN_gp, "PS_e")[c("est", "se", "p_two", "ci95")], list(GP_slope = term_test(mSN_gp, "scale(chi_GP)")[c("est", "p_two")])),
  global_z_composite_6ROI = c(term_test(mGlob, "PS_e")[c("est", "se", "p_two", "ci95")], list(unit = "mean of within-sample z of GP, PUT, CAU, THAL, RN, STN")),
  correlation_SNVTA_with = as.list(round(stats::cor(first$chi_SNVTA, first[, c("chi_GP", "chi_PUT", "chi_CAU", "chi_THAL", "chi_RN", "chi_STN", "chi_DGM")], use = "pairwise")[1, ], 3)))
add_test(data.frame(module = "D2", block = "diagnostic: regional", test = "SN-VTA PS effect adjusted for the deep-GM composite (global-iron covariate)", outcome = "chi_SNVTA",
                    n = mSN_dgm$n, beta_PS = term_test(mSN_dgm, "PS_e")$est, se_hc3 = term_test(mSN_dgm, "PS_e")$se, ci_lo = term_test(mSN_dgm, "PS_e")$ci95[1],
                    ci_hi = term_test(mSN_dgm, "PS_e")$ci95[2], p_two = term_test(mSN_dgm, "PS_e")$p_two, beta_per_sd_raw = term_test(mSN_dgm, "PS_e")$est / stats::sd(first$chi_SNVTA),
                    note = "post hoc", stringsAsFactors = FALSE))
add_test(data.frame(module = "D2", block = "diagnostic: regional", test = "regional PS effect: 6-ROI global z composite (GP, PUT, CAU, THAL, RN, STN)", outcome = "chi_global_z",
                    n = mGlob$n, beta_PS = term_test(mGlob, "PS_e")$est, se_hc3 = term_test(mGlob, "PS_e")$se, ci_lo = term_test(mGlob, "PS_e")$ci95[1],
                    ci_hi = term_test(mGlob, "PS_e")$ci95[2], p_two = term_test(mGlob, "PS_e")$p_two, beta_per_sd_raw = term_test(mGlob, "PS_e")$est / stats::sd(first$chi_global_z),
                    note = "z units", stringsAsFactors = FALSE))
# forest plot of standardised PS effects by ROI (beta / raw SD of the ROI), 95% CI
fr <- reg_rows; fr$lo <- fr$ci_lo / fr$sd_raw; fr$hi <- fr$ci_hi / fr$sd_raw
fr$roi_label <- factor(fr$roi_label, levels = rev(rois))
pF <- ggplot(fr, aes(x = beta_per_sd_raw, y = roi_label)) + geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.6, colour = "#2a78d6") + geom_point(size = 2.4, colour = "#2a78d6") +
  geom_text(aes(label = paste0("n = ", n)), x = max(fr$hi) + 0.08, size = 2.7, colour = "grey30", hjust = 0) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(x = "beta_PS in raw-SD units of the ROI (M1 form, HC3 95% CI)", y = NULL, title = "PS main effect by ROI: SN-specific or global?",
       subtitle = paste0(LABEL, "; positive = higher chi in PS youth")) +
  theme_minimal(base_size = 11) + theme(plot.subtitle = element_text(size = 7.5, colour = "grey35"))
ggsave(file.path(OUT, "fig_D2_forest_PS_effect_by_ROI.png"), pF, width = 7.5, height = 5.2, dpi = 160, bg = "white")
# faceted regional box/strip (seven ROIs named in the task, plus STR and the deep-GM composite)
lv <- c("chi_SNVTA", "chi_SNonly", "chi_VTAPBP", "chi_GP", "chi_PUT", "chi_CAU", "chi_THAL", "chi_RN", "chi_STR")
long <- do.call(rbind, lapply(lv, function(v) data.frame(cell = first$cell, roi = rois[[v]], chi = first[[v]])))
long$roi <- factor(long$roi, levels = rois[lv])
pR <- ggplot(long, aes(x = cell, y = chi, colour = cell)) + geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.35, fill = NA) +
  geom_jitter(width = 0.13, height = 0, size = 1.1, alpha = 0.7, shape = 16) + facet_wrap(~ roi, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = CELL_COL, guide = "none") + scale_x_discrete(labels = c("TD", "ADHD", "PS", "PS+ADHD")) +
  labs(x = NULL, y = "chi (ppb)", title = "Regional chi by cell (first-usable, n = 80)", subtitle = LABEL) +
  theme_minimal(base_size = 10) + theme(panel.grid.major.x = element_blank(), plot.subtitle = element_text(size = 7.5, colour = "grey35"))
ggsave(file.path(OUT, "fig_D2_regional_by_cell.png"), pR, width = 8.5, height = 7.5, dpi = 160, bg = "white")

# 4c. referencing: whole-brain reference value and the CSF reference by cell (is the effect a referencing artefact?)
first$chi_SNVTA_unref <- first$chi_SNVTA + first$ref_ppb   # chi before subtracting the whole-brain reference
DG$reference_raw_by_cell <- do.call(rbind, lapply(c("ref_ppb", "chi_CSF_ref", "chi_SNVTA_ref2", "chi_SNVTA_unref", "ref_frac_wm", "ref_frac_gm", "ref_frac_csf", "fov_brain_frac"),
                                                  function(v) raw_by_cell(first, v)))
r_ref  <- nuis_ps(first, "ref_ppb", "whole-brain reference value (ref_ppb) by PS", "diagnostic: artefact channels", "D3")
r_csf  <- ps_row(fit_y(first, "chi_CSF_ref"), "CSF reference minus whole-brain reference (M1 form = frozen K3 dref, reproduced) by PS", "diagnostic: artefact channels", "D3",
                 note = "must equal the frozen K3 beta_dref 1.486 (p .2742); K3 trigger was FALSE and is not re-decided")
r_csf_red <- nuis_ps(first, "chi_CSF_ref", "CSF reference offset by PS (reduced RHS: no roi_vol / Q_phase; descriptive)", "co-report", "D3")
r_unr  <- ps_row(fit_y(first, "chi_SNVTA_unref"), "SN-VTA chi before whole-brain referencing (chi + ref_ppb): beta_PS", "diagnostic: artefact channels", "D3")
r_ref2 <- ps_row(fit_y(first, "chi_SNVTA_ref2"), "SN-VTA chi CSF-referenced: beta_PS (frozen K3 co-report reproduced)", "co-report", "D3",
                 note = "identical to the frozen K3 CSF-referenced value (3.477, p .3507); shown for the reader, not re-decided")
r_wm   <- nuis_ps(first, "ref_frac_wm", "reference-domain WM fraction by PS", "diagnostic: artefact channels", "D3")
r_csff <- nuis_ps(first, "ref_frac_csf", "reference-domain CSF fraction by PS", "diagnostic: artefact channels", "D3")
r_fov  <- nuis_ps(first, "fov_brain_frac", "brain-in-FOV fraction by PS", "diagnostic: artefact channels", "D3")
DG$reference_tests <- list(ref_ppb = r_ref, csf_ref_m1_form = r_csf, csf_ref_reduced = r_csf_red, unreferenced = r_unr, csf_referenced = r_ref2, wm_frac = r_wm, csf_frac = r_csff, fov = r_fov)
for (x in list(r_ref, r_csf, r_unr, r_wm, r_csff, r_fov)) add_test(x)
cell_plot(first, "ref_ppb", "whole-brain reference (ppb)", "Whole-brain reference value by cell", NULL, "fig_D3_ref_ppb_by_cell.png")
cell_plot(first, "chi_CSF_ref", "CSF reference - whole-brain reference (ppb)", "CSF reference offset by cell", NULL, "fig_D3_csf_ref_by_cell.png")

# 4d. R2* SN-VTA by cell (iron vs myelin direction; K2 frozen one-sided p .9825)
DG$r2s_raw_by_cell <- do.call(rbind, lapply(c("r2s_SNVTA", "r2s_SNonly", "r2s_GP", "r2s_PUT", "r2s_RN", "r2s_DN_R"), function(v) raw_by_cell(first, v)))
mR2 <- fit_y(first, "r2s_SNVTA")
r_r2 <- ps_row(mR2, "R2* SN-VTA: beta_PS (two-sided; frozen K2 is one-sided)", "diagnostic: artefact channels", "D4")
add_test(r_r2)
DG$r2s_tests <- list(SNVTA = r_r2, by_roi = do.call(rbind, lapply(c("r2s_SNonly", "r2s_GP", "r2s_PUT"), function(v) ps_row(fit_y(first, v), paste0("R2* ", v), "co-report", "D4"))),
                     cor_chi_r2s_by_cell = do.call(rbind, lapply(CELLS, function(cl) { x <- first[first$cell == cl, ]
                       ct <- stats::cor.test(x$chi_SNVTA, x$r2s_SNVTA); data.frame(cell = cl, n = nrow(x), r = unname(ct$estimate), lo = ct$conf.int[1], hi = ct$conf.int[2]) })),
                     cor_chi_r2s_all = unname(stats::cor(first$chi_SNVTA, first$r2s_SNVTA)))
cell_plot(first, "r2s_SNVTA", "SN-VTA R2* (s^-1)", "SN-VTA R2* by cell", adj_cell_means(mR2), "fig_D4_r2s_SNVTA_by_cell.png")
pS <- ggplot(first, aes(x = r2s_SNVTA, y = chi_SNVTA, colour = cell)) + geom_point(size = 1.9, alpha = 0.8, shape = 16) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.7) + scale_colour_manual(values = CELL_COL, name = NULL) +
  labs(x = "SN-VTA R2* (s^-1)", y = "SN-VTA chi (ppb)", title = "SN-VTA chi vs R2* by cell", subtitle = paste0(LABEL, "; overall r = ", f3(DG$r2s_tests$cor_chi_r2s_all, 2))) +
  theme_minimal(base_size = 11) + theme(plot.subtitle = element_text(size = 7.5, colour = "grey35"), legend.position = "top")
ggsave(file.path(OUT, "fig_D4_chi_vs_r2s_by_cell.png"), pS, width = 7, height = 4.6, dpi = 160, bg = "white")

# 4e. QC grade / forced review / G5 Dice / motion index / Q_phase / roi_vol / tilt by cell
qc_cont <- c(dice_min = "G5 Dice (minimum of four)", motion_gsr = "magnitude ghost-to-signal motion index", Q_phase = "Q_phase (log phase-residual index)",
             roi_vol = "SN-VTA atlas volume (mm^3)", b0_tilt_deg = "B0 tilt (deg)", pitch_deg = "head pitch (deg)", fimp_SNVTA = "implausible-voxel share, SN-VTA",
             inmask_SNVTA = "in-support share, SN-VTA", slabcov_SNVTA = "slab coverage, SN-VTA", v4_vol_ml = "4th-ventricle volume (ml)", icv_seg_ml = "intracranial volume (ml)")
qc_cont <- qc_cont[names(qc_cont) %in% names(first)]
DG$qc_raw_by_cell <- do.call(rbind, lapply(names(qc_cont), function(v) raw_by_cell(first, v)))
qc_rows <- do.call(rbind, lapply(names(qc_cont), function(v) nuis_ps(first, v, paste0("QC/nuisance by PS: ", qc_cont[[v]]), "diagnostic: artefact channels", "D5")))
DG$qc_tests <- qc_rows
for (v in c("dice_min", "motion_gsr", "Q_phase", "roi_vol")) add_test(qc_rows[qc_rows$outcome == v, ])
bin_tab <- function(v, lab) { t <- table(first$cell, first[[v]], useNA = "no"); p <- if (ncol(t) > 1) stats::fisher.test(t)$p.value else NA
  list(variable = v, label = lab, table = as.data.frame.matrix(t), fisher_p = p) }
DG$qc_binary <- list(qc_grade = bin_tab("qc_grade", "final visual QC grade (0/1; grade 2 never usable)"),
                     forced_review = bin_tab("forced_review", "forced-review flag (automated plausibility / slab coverage)"),
                     dental_D = bin_tab("dental_D", "rubric code D (dental/metal)"), t1_artifact = bin_tab("t1_artifact", "PennLINC T1w 'artifact'"),
                     epoch = bin_tab("epoch", "protocol epoch"), f0_cluster = bin_tab("f0_cluster", "ImagingFrequency cluster"), ref2_type = bin_tab("ref2_type", "CSF reference type"))
for (v in c("qc_grade", "forced_review")) {
  t <- DG$qc_binary[[v]]
  add_test(data.frame(module = "D5", block = "diagnostic: artefact channels", test = paste0("QC by cell (Fisher exact): ", t$label), outcome = v, n = sum(!is.na(first[[v]])),
                      beta_PS = NA, se_hc3 = NA, ci_lo = NA, ci_hi = NA, p_two = t$fisher_p, beta_per_sd_raw = NA, note = "4-cell Fisher test; no estimate", stringsAsFactors = FALSE))
}
qlong <- do.call(rbind, lapply(c("dice_min", "motion_gsr", "Q_phase", "roi_vol"), function(v) data.frame(cell = first$cell, metric = qc_cont[[v]], value = first[[v]])))
qlong$metric <- factor(qlong$metric, levels = qc_cont[c("dice_min", "motion_gsr", "Q_phase", "roi_vol")])
pQ <- ggplot(qlong, aes(x = cell, y = value, colour = cell)) + geom_boxplot(outlier.shape = NA, width = 0.55, linewidth = 0.35, fill = NA) +
  geom_jitter(width = 0.13, height = 0, size = 1.1, alpha = 0.7, shape = 16) + facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = CELL_COL, guide = "none") + scale_x_discrete(labels = c("TD", "ADHD", "PS", "PS+ADHD")) +
  labs(x = NULL, y = NULL, title = "QC / nuisance metrics by cell", subtitle = LABEL) +
  theme_minimal(base_size = 10) + theme(panel.grid.major.x = element_blank(), plot.subtitle = element_text(size = 7.5, colour = "grey35"))
ggsave(file.path(OUT, "fig_D5_qc_by_cell.png"), pQ, width = 7.5, height = 6, dpi = 160, bg = "white")

# 4f. age, sex, race by cell
DG$age_by_cell <- raw_by_cell(first, "age")
DG$age_anova_p <- summary(stats::aov(age ~ cell, data = first))[[1]][["Pr(>F)"]][1]
DG$age_ps_welch_p <- stats::t.test(age ~ PS_e > 0, data = first)$p.value
DG$sex_by_cell <- as.data.frame.matrix(table(first$cell, first$sex)); DG$sex_fisher_p <- stats::fisher.test(table(first$cell, first$sex))$p.value
DG$race3_by_cell <- as.data.frame.matrix(table(first$cell, first$race3_chr)); DG$race3_fisher_p <- stats::fisher.test(table(first$cell, first$race3_chr))$p.value
DG$race_detail_by_cell <- as.data.frame.matrix(table(first$cell, first$race))
DG$ethnicity_by_cell <- as.data.frame.matrix(table(first$cell, ifelse(is.na(first$ethnicity), "missing", first$ethnicity)))
DG$threshold_psychosis_by_cell <- as.data.frame.matrix(table(first$cell, first$dx_psychosis))
# race and chi: does race3 relate to chi, and does dropping race change beta_PS? (frozen sensitivity: 'without race' beta 4.893)
mNoRace <- fit_ols_hc3(stats::as.formula("chi_SNVTA ~ PS_e * ADHD_e + ns(age, df = 3) + sex_e + scale(roi_vol) + scale(Q_phase)"), first)
DG$race_effect_on_chi <- list(race3_coefficients = coef_table(mSN)[grepl("race3", coef_table(mSN)$term), ],
                              chi_by_race3 = do.call(rbind, lapply(c("Black", "White", "Other"), function(r) { x <- first$chi_SNVTA[first$race3_chr == r]; data.frame(race3 = r, n = length(x), mean = mean(x), sd = stats::sd(x)) })),
                              beta_PS_without_race = unname(mNoRace$b["PS_e"]))
add_test(data.frame(module = "D6", block = "diagnostic: artefact channels", test = "age by PS (Welch t, raw)", outcome = "age", n = nrow(first), beta_PS = mean(first$age[first$PS_e > 0]) - mean(first$age[first$PS_e < 0]),
                    se_hc3 = NA, ci_lo = NA, ci_hi = NA, p_two = DG$age_ps_welch_p, beta_per_sd_raw = NA, note = "mean difference in years, PS minus non-PS", stringsAsFactors = FALSE))
add_test(data.frame(module = "D6", block = "diagnostic: artefact channels", test = "race3 by cell (Fisher exact)", outcome = "race3", n = nrow(first), beta_PS = NA, se_hc3 = NA, ci_lo = NA, ci_hi = NA,
                    p_two = DG$race3_fisher_p, beta_per_sd_raw = NA, note = "4-cell Fisher test; no estimate", stringsAsFactors = FALSE))

# 4g. chi_SNVTA vs age by cell
DG$age_slopes_by_cell <- do.call(rbind, lapply(CELLS, function(cl) { x <- first[first$cell == cl, ]; f <- stats::lm(chi_SNVTA ~ age, data = x); ci <- stats::confint(f)["age", ]
  data.frame(cell = cl, n = nrow(x), slope_ppb_per_y = unname(stats::coef(f)["age"]), lo = ci[1], hi = ci[2], p = summary(f)$coefficients["age", 4], age_mean = mean(x$age), age_range = paste(range(x$age), collapse = "-")) }))
mAge <- fit_ols_hc3(stats::as.formula(paste0("chi_SNVTA ~ PS_e * ADHD_e + age + PS_e:age + sex_e", race_term(first), " + scale(roi_vol) + scale(Q_phase)")), first)
mAge2 <- fit_y(first, "chi_SNVTA", extra = " + PS_e:age")
DG$age_x_PS <- list(linear_age_model = term_test(mAge, "PS_e:age")[c("est", "se", "p_two", "ci95")], age_slope_pooled = term_test(mAge, "age")[c("est", "se", "p_two", "ci95")],
                    m1_plus_PSxage = term_test(mAge2, "PS_e:age")[c("est", "se", "p_two", "ci95")],
                    td_only_slope = { x <- first[first$cell == "TD", ]; f <- stats::lm(chi_SNVTA ~ age + sex_e, data = x); c(slope = unname(stats::coef(f)["age"]), p = summary(f)$coefficients["age", 4]) })
tA <- term_test(mAge, "PS_e:age")
add_test(data.frame(module = "D7", block = "diagnostic: artefact channels", test = "PS x age (linear) interaction on SN-VTA chi", outcome = "chi_SNVTA", n = mAge$n, beta_PS = tA$est, se_hc3 = tA$se,
                    ci_lo = tA$ci95[1], ci_hi = tA$ci95[2], p_two = tA$p_two, beta_per_sd_raw = NA, note = "ppb per year, PS minus non-PS", stringsAsFactors = FALSE))
pA <- ggplot(first, aes(x = age, y = chi_SNVTA, colour = cell)) + geom_point(size = 1.9, alpha = 0.8, shape = 16) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.7) + facet_wrap(~ cell, ncol = 4) + scale_colour_manual(values = CELL_COL, guide = "none") +
  labs(x = "age (y)", y = "SN-VTA chi (ppb)", title = "SN-VTA chi vs age by cell (per-cell OLS lines)",
       subtitle = paste0(LABEL, "; PS x age (linear) = ", ci_str(tA$est, tA$ci95[1], tA$ci95[2]), " ppb/y, p = ", fp(tA$p_two))) +
  theme_minimal(base_size = 10) + theme(plot.subtitle = element_text(size = 7.5, colour = "grey35"))
ggsave(file.path(OUT, "fig_D7_chi_vs_age_by_cell.png"), pA, width = 9, height = 3.6, dpi = 160, bg = "white")
RES$diagnostics <- DG

# ------------------------------------------------------------------------------------------------------------
# 5. BH-FDR table across the exploratory tests of this block
#    q_block: BH within (i) the registered X block (X1 x3, X3, X6) and (ii) the post hoc diagnostic block;
#    q_all:   BH across every test in this module (most conservative, for the reader).
# ------------------------------------------------------------------------------------------------------------
log_msg("FDR table")
keep <- c("module", "block", "test", "outcome", "n", "beta_PS", "se_hc3", "ci_lo", "ci_hi", "p_two", "beta_per_sd_raw", "note")
fdr <- do.call(rbind, lapply(tests, function(x) { for (k in setdiff(keep, names(x))) x[[k]] <- NA; x[, keep] }))
names(fdr)[names(fdr) == "beta_PS"] <- "estimate"
fdr$family <- ifelse(fdr$block == "registered X (X1/X3/X6)", "registered X", "post hoc diagnostics")
fdr$q_within_family <- NA_real_
for (fam in unique(fdr$family)) { i <- fdr$family == fam; fdr$q_within_family[i] <- stats::p.adjust(fdr$p_two[i], method = "BH") }
fdr$q_all_module <- stats::p.adjust(fdr$p_two, method = "BH")
fdr$sig_q05_within_family <- !is.na(fdr$q_within_family) & fdr$q_within_family < 0.05
fdr$m_within_family <- ave(fdr$p_two, fdr$family, FUN = length)
fdr <- fdr[order(fdr$family != "registered X", fdr$p_two), ]
utils::write.csv(fdr, file.path(OUT, "fdr_table.csv"), row.names = FALSE)
RES$fdr <- list(n_tests = nrow(fdr), n_registered_X = sum(fdr$family == "registered X"), n_post_hoc = sum(fdr$family != "registered X"),
                any_q05_within_family = any(fdr$sig_q05_within_family), min_q_registered = min(fdr$q_within_family[fdr$family == "registered X"]),
                min_q_post_hoc = min(fdr$q_within_family[fdr$family != "registered X"]), min_q_all = min(fdr$q_all_module))

# ------------------------------------------------------------------------------------------------------------
# 6. outputs
# ------------------------------------------------------------------------------------------------------------
wcsv <- function(x, f) utils::write.csv(x, file.path(OUT, f), row.names = FALSE)
wcsv(rbind(x1$raw_by_cell, DG$snvta_raw_by_cell, DG$regional_raw_by_cell, DG$reference_raw_by_cell, DG$r2s_raw_by_cell, DG$qc_raw_by_cell, DG$age_by_cell), "raw_by_cell.csv")
wcsv(rbind(cbind(outcome = "chi_SNVTA", DG$snvta_adj_cell_means), x1$adj_cell_means, cbind(outcome = "r2s_SNVTA", adj_cell_means(mR2))), "adjusted_cell_means.csv")
wcsv(rbind(qc_rows, r_ref, r_csf_red, r_wm, r_csff, r_fov), "qc_and_reference_by_PS.csv")
wcsv(rbind(reg_rows[, setdiff(names(reg_rows), "roi_label")], x1$coreports, r_unr, r_ref2, r_csf, r_r2), "group_effects_by_outcome.csv")
wcsv(DG$age_slopes_by_cell, "age_slopes_by_cell.csv")
wcsv(x6$rsas_by_cell, "X6_rsas_by_cell.csv")
jsonlite::write_json(RES, file.path(OUT, "results_tier4.json"), auto_unbox = TRUE, pretty = TRUE, digits = 8, na = "null")
writeLines(.LOG$lines, file.path(OUT, "run.log"))
log_msg("done")
