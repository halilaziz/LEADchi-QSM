#!/usr/bin/env Rscript
# LEAD-chi TASK R: adversarial re-derivation of the frozen confirmatory quantities from the unlocked table.
# EXPLORATORY / REVIEW ONLY (Tier 3-4 labelling applies to every number produced here). Nothing here re-decides H1.
# Written independently of the frozen analysis code: own coding, own first-usable rule, HC3 by explicit matrix algebra
# (sandwich::vcovHC is used only as a cross-check), own contrast algebra, own permutation and bootstrap code.
suppressPackageStartupMessages({ library(splines); library(MASS); library(jsonlite) })
set.seed(20261005)

IN   <- "/path/to/Penn_LEAD/LEADchi/mirror/unlocked/merged_qsmxt_primary_SNVTA.csv"
GATE <- "/path/to/Penn_LEAD/LEADchi/mirror/gates/gates_preunlock.json"
RES  <- "/path/to/Penn_LEAD/LEADchi/mirror/results/confirmatory/results.json"
OUT  <- "/path/to/Penn_LEAD/LEADchi/results_exploratory/review"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(IN, stringsAsFactors = FALSE, na.strings = c("NA", ""))
g <- jsonlite::read_json(GATE, simplifyVector = TRUE)
R <- jsonlite::read_json(RES, simplifyVector = TRUE)
out <- list()

# ---------------------------------------------------------------- coding and samples (own implementation)
d$PS_e <- ifelse(d$study_group == "PRO/CHR", 0.5, -0.5)
d$ADHD_e <- ifelse(d$dx_adhd == 1, 0.5, -0.5)
d$sex_e <- ifelse(d$sex == "M", 0.5, -0.5)
d$race3c <- ifelse(d$race == "Black or African American", "Black", ifelse(d$race == "White", "White", "Other"))
d$cell <- ifelse(d$PS_e > 0, ifelse(d$ADHD_e > 0, "PS+ADHD", "PS-only"), ifelse(d$ADHD_e > 0, "ADHD-only", "TD"))
d$cell <- factor(d$cell, levels = c("TD", "ADHD-only", "PS-only", "PS+ADHD"))
d$sess_key <- paste0(d$participant_id, "_", d$session_id)
# plausibility re-verification (own): envelopes on the ROI columns used below
env <- list(chi_SNVTA = c(-50, 300), chi_STR = c(-80, 200), chi_GP = c(-50, 300), chi_PUT = c(-80, 200), chi_DN_L = c(-50, 300))
plaus <- do.call(rbind, lapply(names(env), function(v) {
  i <- which(!is.na(d[[v]]) & (d[[v]] < env[[v]][1] | d[[v]] > env[[v]][2]))
  if (!length(i)) return(NULL)
  data.frame(sess_key = d$sess_key[i], column = v, value = d[[v]][i], usable = d$usable[i]) }))
out$plausibility_hits <- plaus
stopifnot(!any(plaus$usable[plaus$column == "chi_SNVTA"] == 1))
d <- d[order(d$participant_id, d$ses_order), ]
u <- d[d$usable == 1, ]
first <- u[!duplicated(u$participant_id), ]
flag_first <- d$sess_key[d$first_usable %in% 1]
out$samples <- list(n_rows = nrow(d), n_subjects = length(unique(d$participant_id)), n_usable_sessions = nrow(u),
                    n_first_usable = nrow(first), first_usable_flag_agrees = setequal(flag_first, first$sess_key),
                    cells_first = as.list(table(first$cell)), cells_all_sessions = as.list(table(d$cell)),
                    n_subjects_by_cell = as.list(table(d$cell[!duplicated(d$participant_id)])))
cat("first-usable n =", nrow(first), "; flag agrees:", out$samples$first_usable_flag_agrees, "\n")

set_race <- function(x) { lev <- intersect(c("Black", "White", "Other"), unique(x$race3c))
  x$race3 <- factor(x$race3c, levels = lev); if (length(lev) >= 2) contrasts(x$race3) <- contr.sum(length(lev)); x }
first <- set_race(first)

# ---------------------------------------------------------------- HC3 by explicit algebra
hc <- function(fit, type = "HC3") {
  X <- model.matrix(fit); e <- resid(fit); h <- hatvalues(fit); n <- nrow(X); p <- ncol(X)
  w <- switch(type, HC0 = e^2, HC1 = e^2 * n / (n - p), HC2 = e^2 / (1 - h), HC3 = e^2 / (1 - h)^2)
  XtXi <- solve(crossprod(X)); XtXi %*% crossprod(X * w, X) %*% XtXi }
ctest <- function(fit, cvec, V = NULL, type = "HC3") {
  if (is.null(V)) V <- hc(fit, type); b <- coef(fit); cv <- setNames(rep(0, length(b)), names(b)); cv[names(cvec)] <- cvec
  est <- sum(cv * b); se <- sqrt(drop(t(cv) %*% V %*% cv)); df <- df.residual(fit); t <- est / se
  list(est = est, se = se, t = t, df = df, p_two = 2 * pt(-abs(t), df), p_lower = pt(t, df), p_upper = pt(t, df, lower.tail = FALSE),
       ci95 = est + c(-1, 1) * qt(.975, df) * se, ci90 = est + c(-1, 1) * qt(.95, df) * se) }

fo <- chi_SNVTA ~ PS_e * ADHD_e + ns(age, df = 3) + sex_e + race3 + scale(roi_vol) + scale(Q_phase)
fit <- lm(fo, data = first)
V3 <- hc(fit, "HC3")
V3_sandwich <- sandwich::vcovHC(fit, type = "HC3")
out$hc3_algebra_vs_sandwich_maxabsdiff <- max(abs(V3 - V3_sandwich))
m1 <- ctest(fit, c(PS_e = 1), V3)
sig <- sigma(fit); dfr <- df.residual(fit)
rob <- MASS::rlm(fo, data = first, psi = MASS::psi.huber, k = 1.345, scale.est = "proposal 2", maxit = 200)
sig_hub <- rob$s; sig_mad <- 1.4826 * median(abs(resid(fit) - median(resid(fit))))
sigB <- min(sig, sig_hub)
rL <- g$rL_primary; dB <- 0.6 * sqrt(min(1, rL))
tcrit <- qt(.95, dfr)
tB <- (m1$est + dB * sigB) / m1$se; pB <- pt(tB, dfr, lower.tail = FALSE)
out$M1 <- list(n = nobs(fit), p_params = length(coef(fit)), df = dfr, beta = m1$est, se_hc3 = m1$se, t = m1$t, p_two = m1$p_two,
               ci95 = m1$ci95, ci90 = m1$ci90, sigma_resid = sig, sigma_huber = sig_hub, huber_converged = rob$converged,
               sigma_mad = sig_mad, sigma_B = sigB, delta = m1$est / sig, delta_B = dB, bound_ppb = -dB * sigB,
               tB = tB, pB = pB, tcrit_one_sided = tcrit, A = (m1$est < 0 && m1$p_two < .05), B = (pB < .05),
               outcome = if (m1$est < 0 && m1$p_two < .05) { if (pB < .05) "A+B" else "A" } else if (pB < .05) "B" else "C",
               huber_beta = unname(coef(rob)["PS_e"]), huber_se = unname(summary(rob)$coefficients["PS_e", 2]),
               R2 = summary(fit)$r.squared, se_classical = unname(sqrt(diag(vcov(fit)))["PS_e"]),
               se_HC0 = sqrt(hc(fit, "HC0")["PS_e", "PS_e"]), se_HC1 = sqrt(hc(fit, "HC1")["PS_e", "PS_e"]), se_HC2 = sqrt(hc(fit, "HC2")["PS_e", "PS_e"]),
               # the largest deficit (in sigma_B units) excluded at one-sided 5 %: the lower one-sided 95 % bound of delta
               delta_lower_1sided95_sigmaB = (m1$est - tcrit * m1$se) / sigB,
               delta_lower_1sided95_sigmaresid = (m1$est - tcrit * m1$se) / sig,
               pB_at_deltaB_zero = pt(m1$t, dfr, lower.tail = FALSE))
# compare with the frozen report
frz <- R$H1$decision
out$M1_vs_frozen <- data.frame(quantity = c("beta", "se_hc3", "p_two", "sigma_resid", "sigma_rob", "sigma_B", "delta", "delta_B", "tB", "pB"),
  review = c(m1$est, m1$se, m1$p_two, sig, sig_hub, sigB, m1$est / sig, dB, tB, pB),
  frozen = c(frz$beta_ppb, frz$se_hc3, frz$p_two, frz$sigma_resid, frz$sigma_rob, frz$sigma_B, frz$delta, frz$B_registered$delta_B, frz$B_registered$tB, frz$B_registered$pB))
out$M1_vs_frozen$rel_diff <- abs(out$M1_vs_frozen$review - out$M1_vs_frozen$frozen) / pmax(abs(out$M1_vs_frozen$frozen), 1e-12)
print(out$M1_vs_frozen, digits = 10)

# ---------------------------------------------------------------- Bound-test sensitivity grid (scale x SE x bound)
scales <- c(sigma_resid = sig, huber = sig_hub, mad_ols = sig_mad, sd_total_pooled = g$pooled$sd_total_SNVTA,
            sd_TD_raw = sd(first$chi_SNVTA[first$cell == "TD"]), sd_raw_all = sd(first$chi_SNVTA))
ses <- c(classical = out$M1$se_classical, HC0 = out$M1$se_HC0, HC1 = out$M1$se_HC1, HC2 = out$M1$se_HC2, HC3 = m1$se)
bounds <- c(registered_0.581 = dB, unadjusted_0.6 = 0.6, half_0.5 = 0.5, cross_session = 0.6 * sqrt(g$rL_cross),
            resid_split = 0.6 * sqrt(R$H1$B_coreport_residualised_split$rL), d0.3 = 0.3, d0.2 = 0.2, d0.1 = 0.1, d0.05 = 0.05, zero = 0)
grid <- expand.grid(scale = names(scales), se = names(ses), bound = names(bounds), stringsAsFactors = FALSE)
grid$sigma <- scales[grid$scale]; grid$SE <- ses[grid$se]; grid$delta_B <- bounds[grid$bound]
grid$tB <- (m1$est + grid$delta_B * grid$sigma) / grid$SE; grid$pB <- pt(grid$tB, dfr, lower.tail = FALSE); grid$B <- grid$pB < .05
write.csv(grid, file.path(OUT, "bound_test_sensitivity_grid.csv"), row.names = FALSE)
out$bound_grid_summary <- list(n_cells = nrow(grid), share_B = mean(grid$B),
  non_B_cells = grid[!grid$B, c("scale", "se", "bound", "pB")],
  sigmaB_needed_for_B_at_registered_bound = (tcrit * m1$se - m1$est) / dB,
  note = "sigma_B needed for B is the residual scale at which p_B would equal .05 under the registered bound and HC3 SE")

# ---------------------------------------------------------------- K1 simple effects, multiplicity, max-T permutation
k1 <- list(PSonly_minus_TD = ctest(fit, c(PS_e = 1, `PS_e:ADHD_e` = -0.5), V3),
           PSADHD_minus_ADHDonly = ctest(fit, c(PS_e = 1, `PS_e:ADHD_e` = 0.5), V3),
           ADHD_main = ctest(fit, c(ADHD_e = 1), V3), interaction = ctest(fit, c(`PS_e:ADHD_e` = 1), V3))
k1tab <- data.frame(contrast = names(k1), est = sapply(k1, `[[`, "est"), se = sapply(k1, `[[`, "se"), p_two = sapply(k1, `[[`, "p_two"))
k1tab$p_holm_2simple <- NA; k1tab$p_holm_2simple[1:2] <- p.adjust(k1tab$p_two[1:2], "holm")
k1tab$p_holm_4contrasts <- p.adjust(k1tab$p_two, "holm")
k1tab$p_holm_with_betaPS <- NA
# the family actually interpreted in the report: beta_PS (registered) + 2 simple effects
k1tab$p_holm_with_betaPS[1:2] <- p.adjust(c(m1$p_two, k1tab$p_two[1:2]), "holm")[2:3]
# TOST on ADHD at +/- 0.5 sigma_B (registered K1b component; here only reproduced)
bnd <- 0.5 * sigB; ta <- k1$ADHD_main
out$tost <- list(bound_ppb = bnd, p_lower = pt((ta$est + bnd) / ta$se, dfr, lower.tail = FALSE), p_upper = pt((ta$est - bnd) / ta$se, dfr))
# max-T Freedman-Lane permutation (PS labels exchangeable within ADHD strata): reduced model drops PS_e and PS_e:ADHD_e
X <- model.matrix(fit); y <- first$chi_SNVTA; n <- nrow(X)
jj <- match(c("PS_e", "PS_e:ADHD_e"), colnames(X)); Xr <- X[, -jj, drop = FALSE]
fr <- drop(Xr %*% solve(crossprod(Xr), crossprod(Xr, y))); er <- y - fr
XtXi <- solve(crossprod(X)); A <- XtXi %*% t(X); h <- rowSums((X %*% XtXi) * X)
cmat <- rbind(PSonly_TD = c(PS_e = 1, `PS_e:ADHD_e` = -0.5), PSADHD_ADHD = c(PS_e = 1, `PS_e:ADHD_e` = 0.5))
Cfull <- matrix(0, 2, ncol(X), dimnames = list(rownames(cmat), colnames(X))); Cfull[, colnames(cmat)] <- cmat
Acon <- Cfull %*% A                                   # 2 x n: contrast estimates = Acon %*% Y
tcon <- function(Y) { B <- A %*% Y; Rm <- Y - X %*% B; E <- Rm / (1 - h)
  est <- Acon %*% Y; se <- sqrt(sapply(1:2, function(k) colSums((Acon[k, ] * E)^2))); if (is.null(dim(se))) se <- matrix(se, ncol = 2)
  est / t(se) }
t_obs <- tcon(matrix(y, ncol = 1))
stopifnot(abs(t_obs[1] - k1$PSonly_minus_TD$t) < 1e-8, abs(t_obs[2] - k1$PSADHD_minus_ADHDonly$t) < 1e-8)
nperm <- 10000; strata <- first$ADHD_e; cnt_max <- 0; cnt_single <- 0; cnt_beta <- 0
ab <- A[match("PS_e", colnames(X)), ]
tbeta <- function(Y) { B <- A %*% Y; Rm <- Y - X %*% B; B[match("PS_e", colnames(X)), ] / sqrt(colSums((ab * Rm / (1 - h))^2)) }
tb_obs <- tbeta(matrix(y, ncol = 1))
for (chunk in split(seq_len(nperm), ceiling(seq_len(nperm) / 1000))) {
  m <- length(chunk); idx <- matrix(rep(seq_len(n), m), n, m)
  for (s in unique(strata)) { r <- which(strata == s); idx[r, ] <- apply(idx[r, , drop = FALSE], 2, function(v) v[sample.int(length(v))]) }
  Ys <- fr + matrix(er[idx], n, m); Ts <- tcon(Ys)
  cnt_max <- cnt_max + sum(apply(abs(Ts), 2, max) >= abs(t_obs[1]) - 1e-12)   # max over the 2 contrasts, per permutation
  cnt_single <- cnt_single + sum(abs(Ts[1, ]) >= abs(t_obs[1]) - 1e-12)
  cnt_beta <- cnt_beta + sum(abs(tbeta(Ys)) >= abs(tb_obs) - 1e-12) }
out$K1 <- list(table = k1tab, PSonly_TD_perm_p_unadjusted = (cnt_single + 1) / (nperm + 1),
               PSonly_TD_perm_p_maxT_over_2_simple_effects = (cnt_max + 1) / (nperm + 1),
               betaPS_freedman_lane_p_two_sided = (cnt_beta + 1) / (nperm + 1), n_perm = nperm,
               frozen_PSonly_TD = R$H1$K1$simple_PSonly_minus_TD[c("est", "se", "p_two")],
               frozen_freedman_lane_p = R$H1$freedman_lane$p)
print(k1tab, digits = 4); cat("max-T perm p (PS-only - TD):", out$K1$PSonly_TD_perm_p_maxT_over_2_simple_effects, "\n")

# ---------------------------------------------------------------- bootstrap delta CI and BF (reproduction, own code)
nboot <- 5000; cells <- split(seq_len(n), first$cell); bd <- numeric(nboot)
for (b in seq_len(nboot)) { i <- unlist(lapply(cells, function(ix) ix[sample.int(length(ix), replace = TRUE)]))
  f <- lm(fo, data = first[i, ]); bd[b] <- coef(f)["PS_e"] / sigma(f) }
out$boot_delta_ci95 <- unname(quantile(bd, c(.025, .975), na.rm = TRUE)); out$boot_delta_ci95_frozen <- R$H1$boot_delta$ci95
sed <- m1$se / sig; dl <- m1$est / sig
bf <- function(mu, s) dnorm(dl, mu, sqrt(sed^2 + s^2)) / dnorm(dl, 0, sed)
out$bayes <- list(BF10_attenuated = bf(-0.35, 0.25), BF10_sceptical = bf(0, 0.5), frozen = R$H1$bayes$BF)

# ---------------------------------------------------------------- K2, K3
fit2 <- lm(update(fo, r2s_SNVTA ~ .), data = first); k2 <- ctest(fit2, c(PS_e = 1))
first$dref <- first$chi_SNVTA - first$chi_SNVTA_ref2
fit3 <- lm(update(fo, dref ~ .), data = first); k3 <- ctest(fit3, c(PS_e = 1))
fitc <- lm(update(fo, chi_SNVTA_ref2 ~ .), data = first); kc <- ctest(fitc, c(PS_e = 1))
robc <- MASS::rlm(update(fo, chi_SNVTA_ref2 ~ .), data = first, psi = MASS::psi.huber, k = 1.345, scale.est = "proposal 2", maxit = 200)
sigBc <- min(sigma(fitc), robc$s)
out$K2 <- list(beta = k2$est, se = k2$se, p_one_sided_lower = k2$p_lower, p_two = k2$p_two, frozen_est = R$H1$K2$est)
out$K3 <- list(beta_dref = k3$est, p = k3$p_two, threshold = max(.25 * abs(m1$est), .10 * sig), trigger = abs(k3$est) > max(.25 * abs(m1$est), .10 * sig) && k3$p_two < .10,
               n_no_csf = sum(is.na(first$dref)), csf_beta = kc$est, csf_p = kc$p_two,
               csf_pB = pt((kc$est + dB * sigBc) / kc$se, df.residual(fitc), lower.tail = FALSE), frozen = R$H1$K3[c("beta_dref", "p_two", "trigger", "csf_beta", "csf_p", "csf_pB")])

# ---------------------------------------------------------------- S1, S2, Holm
s1d <- first[!is.na(first$chi_STR) & !is.na(first$str_vol), ]
s1d$D <- (s1d$chi_SNVTA - g$pooled$mean_SNVTA) / g$pooled$sd_SNVTA - (s1d$chi_STR - g$pooled$mean_STR) / g$pooled$sd_STR
fs1 <- lm(D ~ PS_e * ADHD_e + ns(age, 3) + sex_e + race3 + scale(roi_vol) + scale(str_vol) + scale(Q_phase), data = s1d)
L <- ctest(fs1, c(PS_e = 1, ADHD_e = -1))
s2d <- first[first$age >= 11 & !is.na(first$PRIME_total) & !is.na(first$BDI_total) & !is.na(first$ARI_total), ]
s2d <- set_race(s2d)
fs2 <- lm(chi_SNVTA ~ scale(PRIME_total) + scale(BDI_total) + scale(ARI_total) + ADHD_e + ns(age, 3) + sex_e + race3 + scale(roi_vol) + scale(Q_phase), data = s2d)
s2 <- ctest(fs2, c(`scale(PRIME_total)` = 1))
holm <- p.adjust(c(S1 = L$p_two, S2 = s2$p_two), "holm")
out$S1 <- list(n = nobs(fs1), L = L$est, se = L$se, p_two = L$p_two, df = L$df, beta_PS_D = ctest(fs1, c(PS_e = 1))$est, beta_ADHD_D = ctest(fs1, c(ADHD_e = 1))$est,
               frozen = R$S1$L[c("est", "se", "p_two")])
out$S2 <- list(n = nobs(fs2), cells = as.list(table(s2d$cell)), beta_PRIME = s2$est, se = s2$se, p_two = s2$p_two, df = s2$df, frozen = R$S2$beta_PRIME[c("est", "se", "p_two")])
out$holm <- list(review = as.list(holm), frozen = R$holm$p_holm)
cat(sprintf("S1 L = %.4f (SE %.4f, p %.4f); S2 beta_PRIME = %.4f (SE %.4f, p %.4f); Holm %.4f / %.4f\n", L$est, L$se, L$p_two, s2$est, s2$se, s2$p_two, holm[1], holm[2]))

# ---------------------------------------------------------------- PC2x (registered) and its variants (exploratory)
pc2 <- function(v, dd, extra = NULL) { f <- lm(as.formula(paste(v, "~ age", if (!is.null(extra)) paste("+", extra) else "")), data = dd)
  s <- summary(f)$coefficients["age", ]; c(n = nobs(f), slope = unname(s[1]), t = unname(s[3]), p_one = unname(pt(s[3], df.residual(f), lower.tail = FALSE))) }
out$PC2x <- list(registered = list(PUT = pc2("chi_PUT", first), GP = pc2("chi_GP", first)),
                 with_sex_cell = list(PUT = pc2("chi_PUT", first, "sex_e + cell"), GP = pc2("chi_GP", first, "sex_e + cell")),
                 all_usable_sessions = list(PUT = pc2("chi_PUT", u), GP = pc2("chi_GP", u)),
                 TD_only = list(PUT = pc2("chi_PUT", first[first$cell == "TD", ]), GP = pc2("chi_GP", first[first$cell == "TD", ])),
                 spearman = list(PUT = cor.test(first$age, first$chi_PUT, method = "spearman", alternative = "greater", exact = FALSE)[c("estimate", "p.value")],
                                 GP = cor.test(first$age, first$chi_GP, method = "spearman", alternative = "greater", exact = FALSE)[c("estimate", "p.value")],
                                 SNVTA = cor.test(first$age, first$chi_SNVTA, method = "spearman", alternative = "greater", exact = FALSE)[c("estimate", "p.value")]),
                 age_range = range(first$age), age_by_cell = tapply(first$age, first$cell, function(a) c(mean = mean(a), sd = sd(a))),
                 r_age_PS = cor(first$age, first$PS_e), frozen = R$PC2x)
# exploratory: how much does beta_PS depend on the age adjustment? (never decisive)
out$age_adjust_sensitivity <- data.frame(model = c("registered ns(age,3)", "linear age", "no age"),
  beta_PS = c(m1$est, ctest(lm(chi_SNVTA ~ PS_e * ADHD_e + age + sex_e + race3 + scale(roi_vol) + scale(Q_phase), data = first), c(PS_e = 1))$est,
              ctest(lm(chi_SNVTA ~ PS_e * ADHD_e + sex_e + race3 + scale(roi_vol) + scale(Q_phase), data = first), c(PS_e = 1))$est),
  p_two = c(m1$p_two, ctest(lm(chi_SNVTA ~ PS_e * ADHD_e + age + sex_e + race3 + scale(roi_vol) + scale(Q_phase), data = first), c(PS_e = 1))$p_two,
            ctest(lm(chi_SNVTA ~ PS_e * ADHD_e + sex_e + race3 + scale(roi_vol) + scale(Q_phase), data = first), c(PS_e = 1))$p_two))

# ---------------------------------------------------------------- split-echo reliability: raw vs covariate-residualised (own)
sp <- first[!is.na(first$chi_SNVTA_split13) & !is.na(first$chi_SNVTA_split24), ]
sb <- function(r) 2 * r / (1 + r); fz_lo <- function(r, n) tanh(atanh(r) - 1.959964 / sqrt(n - 3))
r_raw <- cor(sp$chi_SNVTA_split13, sp$chi_SNVTA_split24)
r13 <- resid(lm(update(fo, chi_SNVTA_split13 ~ .), data = sp)); r24 <- resid(lm(update(fo, chi_SNVTA_split24 ~ .), data = sp))
r_res <- cor(r13, r24)
# which covariates carry the shared variance? R2 of each half on the covariates
out$split_echo <- list(n = nrow(sp), r_raw = r_raw, R_SB_raw = sb(r_raw), rL_raw = sb(fz_lo(r_raw, nrow(sp))), gates_rL = g$rL_primary, gates_n = g$PC4_detail$n,
                       r_resid = r_res, R_SB_resid = sb(r_res), rL_resid = sb(fz_lo(r_res, nrow(sp))),
                       R2_split13_on_covariates = summary(lm(update(fo, chi_SNVTA_split13 ~ .), data = sp))$r.squared,
                       R2_split24_on_covariates = summary(lm(update(fo, chi_SNVTA_split24 ~ .), data = sp))$r.squared,
                       R2_M1 = summary(fit)$r.squared, sd_split13 = sd(sp$chi_SNVTA_split13), sd_split24 = sd(sp$chi_SNVTA_split24), sd_full = sd(sp$chi_SNVTA),
                       r_split13_full = cor(sp$chi_SNVTA_split13, sp$chi_SNVTA), r_split24_full = cor(sp$chi_SNVTA_split24, sp$chi_SNVTA),
                       r_resid_age_only = cor(resid(lm(chi_SNVTA_split13 ~ ns(age, 3), sp)), resid(lm(chi_SNVTA_split24 ~ ns(age, 3), sp))),
                       r_resid_roivol_Qphase_only = cor(resid(lm(chi_SNVTA_split13 ~ scale(roi_vol) + scale(Q_phase), sp)), resid(lm(chi_SNVTA_split24 ~ scale(roi_vol) + scale(Q_phase), sp))),
                       cor_half_with_roi_vol = c(s13 = cor(sp$chi_SNVTA_split13, sp$roi_vol), s24 = cor(sp$chi_SNVTA_split24, sp$roi_vol), full = cor(sp$chi_SNVTA, sp$roi_vol)),
                       cor_half_with_Q_phase = c(s13 = cor(sp$chi_SNVTA_split13, sp$Q_phase), s24 = cor(sp$chi_SNVTA_split24, sp$Q_phase), full = cor(sp$chi_SNVTA, sp$Q_phase)))

# ---------------------------------------------------------------- realised power / sensitivity (post hoc, design-based on the realised SE)
se_d <- m1$se / sig; se_dB <- m1$se / sigB
pw <- function(delta, se, two = TRUE) { if (two) { q <- qt(.975, dfr); pt(-q, dfr, ncp = delta / se) + pt(q, dfr, ncp = delta / se, lower.tail = FALSE) } else pt(qt(.95, dfr), dfr, ncp = delta / se, lower.tail = FALSE) }
out$power <- list(se_delta_realised_resid = se_d, se_delta_realised_sigmaB = se_dB, se_delta_registered_n112 = 0.206, se_delta_amendment2_n80 = 0.24,
  power_A = sapply(c(d0.3 = .3, d0.4 = .4, d0.5 = .5, d0.6 = .6), function(x) pw(x, se_d)),
  MDE80_A_sigma = uniroot(function(x) pw(x, se_d) - .8, c(0.1, 3))$root,
  power_B_registered_bound = sapply(c(true0 = 0, true_m0.1 = -.1, true_m0.2 = -.2, true_m0.3 = -.3), function(x) pw(x + dB, se_dB, two = FALSE)),
  false_B_at_bound = pw(0, se_dB, two = FALSE),
  interaction_se_sigma = k1$interaction$se / sig, MDE80_interaction_sigma = uniroot(function(x) pw(x, k1$interaction$se / sig) - .8, c(0.1, 6))$root,
  note = "power computed from the realised HC3 SE with a noncentral t (df 68); post hoc, descriptive only")

# ---------------------------------------------------------------- QC selection effects by cell (session and subject level)
d$grade2 <- as.integer(d$qc_grade %in% 2); d$g5fail <- as.integer(d$g5_pass %in% 0)
tok <- strsplit(ifelse(is.na(d$fail_reason), "", d$fail_reason), ";")
all_tok <- sort(unique(unlist(tok))); all_tok <- all_tok[all_tok != ""]
tokmat <- sapply(all_tok, function(t) vapply(tok, function(z) t %in% z, logical(1)))
fr_by_cell <- aggregate(tokmat, by = list(cell = d$cell), FUN = sum); fr_by_cell$n_sessions <- as.vector(table(d$cell))
write.csv(fr_by_cell, file.path(OUT, "fail_reasons_by_cell_sessions.csv"), row.names = FALSE)
ft <- function(tab) tryCatch(fisher.test(tab, simulate.p.value = nrow(tab) * ncol(tab) > 4, B = 20000)$p.value, error = function(e) NA)
qc_sess <- list(
  qc_grade_by_cell = table(d$cell, d$qc_grade, useNA = "ifany"),
  grade2_share_by_cell = tapply(d$grade2, d$cell, mean), fisher_grade2_by_cell = ft(table(d$cell, d$grade2)),
  g5fail_share_by_cell = tapply(d$g5fail, d$cell, mean), fisher_g5fail_by_cell = ft(table(d$cell, d$g5fail)),
  usable_share_by_cell = tapply(d$usable, d$cell, mean), fisher_usable_by_cell = ft(table(d$cell, d$usable)),
  grade2_share_by_ADHD = tapply(d$grade2, d$ADHD_e, mean), fisher_grade2_by_ADHD = ft(table(d$ADHD_e, d$grade2)),
  g5fail_share_by_ADHD = tapply(d$g5fail, d$ADHD_e, mean), fisher_g5fail_by_ADHD = ft(table(d$ADHD_e, d$g5fail)),
  usable_share_by_ADHD = tapply(d$usable, d$ADHD_e, mean), fisher_usable_by_ADHD = ft(table(d$ADHD_e, d$usable)),
  usable_share_by_PS = tapply(d$usable, d$PS_e, mean), fisher_usable_by_PS = ft(table(d$PS_e, d$usable)),
  dental_D_by_cell = tapply(d$dental_D %in% 1, d$cell, sum))
# subject level (reproduce the report's Fisher tests)
subj <- d[!duplicated(d$participant_id), ]; subj$in_H1 <- subj$participant_id %in% first$participant_id
qc_subj <- list(in_H1_by_cell = table(subj$cell, subj$in_H1), fisher_in_H1 = fisher.test(table(subj$cell, subj$in_H1))$p.value,
                earliest_usable_by_cell = table(subj$cell, subj$usable), fisher_earliest = fisher.test(table(subj$cell, subj$usable))$p.value,
                earliest_usable_by_ADHD = table(subj$ADHD_e, subj$usable), fisher_earliest_by_ADHD = fisher.test(table(subj$ADHD_e, subj$usable))$p.value,
                frozen = list(p_in_H1 = R$qc$fisher_p_in_H1, p_earliest = R$qc$fisher_p_earliest))
# does the visual grade / G5 relate to chi itself (where chi exists), overall and by cell?
dc <- d[!is.na(d$chi_SNVTA) & d$chi_SNVTA > -50, ]
chi_by_grade <- list(n = nrow(dc), mean_by_grade = tapply(dc$chi_SNVTA, dc$qc_grade, mean), sd_by_grade = tapply(dc$chi_SNVTA, dc$qc_grade, sd),
                     kruskal_p = kruskal.test(chi_SNVTA ~ factor(qc_grade), data = dc)$p.value,
                     mean_by_g5 = tapply(dc$chi_SNVTA, dc$g5fail, mean), wilcox_g5_p = wilcox.test(chi_SNVTA ~ g5fail, data = dc)$p.value,
                     mean_by_usable = tapply(dc$chi_SNVTA, dc$usable, mean), wilcox_usable_p = wilcox.test(chi_SNVTA ~ usable, data = dc)$p.value,
                     mean_by_cell_x_usable = tapply(dc$chi_SNVTA, list(dc$cell, dc$usable), mean),
                     # beta_PS on ALL sessions with a plausible chi (one per subject, earliest), ignoring every QC rule: descriptive only
                     n_cells_grade2 = table(dc$cell, dc$grade2))
dcf <- dc[order(dc$participant_id, dc$ses_order), ]; dcf <- dcf[!duplicated(dcf$participant_id), ]; dcf <- set_race(dcf)
fit_noqc <- lm(fo, data = dcf); nq <- ctest(fit_noqc, c(PS_e = 1))
chi_by_grade$M1_no_QC_at_all <- list(n = nobs(fit_noqc), beta = nq$est, se = nq$se, p_two = nq$p_two, note = "earliest session with any plausible chi per subject; ignores visual QC, placement, G5 and T1w rules; descriptive")
# automated-exclusions-only reproduction (fail reasons only among the visual ones)
vis_tok <- all_tok[grepl("qc_grade|placement", all_tok)]
relax <- vapply(tok, function(z) length(z) > 0 && all(z %in% vis_tok), logical(1))
da <- d; da$usable[relax] <- 1; ua <- da[da$usable == 1 & !is.na(da$chi_SNVTA), ]; fa <- set_race(ua[!duplicated(ua$participant_id), ])
fit_auto <- lm(fo, data = fa); na_ <- ctest(fit_auto, c(PS_e = 1))
chi_by_grade$M1_automated_only <- list(n = nobs(fit_auto), beta = na_$est, se = na_$se, p_two = na_$p_two, visual_tokens_relaxed = vis_tok,
                                        frozen = R$sensitivity$table[grepl("automated", R$sensitivity$table$sensitivity), c("n", "beta_ppb", "p_two")])
chi_by_grade$strict_grade0_first_usable_n <- sum(first$qc_grade %in% 0); chi_by_grade$grade0_sessions_all <- sum(d$qc_grade %in% 0)
# registered Q13 / amendment-1b by-cell tabulations the report omitted (first-usable sample)
cov_by_cell <- function(v) { x <- first[[v]]; if (all(is.na(x))) return(NULL)
  data.frame(variable = v, t(tapply(x, first$cell, function(z) round(median(z, na.rm = TRUE), 4))), kruskal_p = kruskal.test(x ~ first$cell)$p.value,
             wilcox_PS_p = wilcox.test(x ~ first$PS_e)$p.value, wilcox_ADHD_p = wilcox.test(x ~ first$ADHD_e)$p.value, check.names = FALSE) }
cbc <- do.call(rbind, lapply(c("dice_min", "motion_gsr", "Q_phase", "fov_brain_frac", "inmask_SNVTA", "fimp_SNVTA", "slabcov_SNVTA", "pitch_deg", "roll_deg", "yaw_deg", "abs_pitch",
                               "ref_frac_wm", "ref_frac_csf", "roi_vol", "age", "ref2_minus_wb_ppb"), cov_by_cell))
write.csv(cbc, file.path(OUT, "covariates_by_cell_first_usable.csv"), row.names = FALSE)
out$qc_session <- lapply(qc_sess, function(x) if (is.table(x)) as.data.frame.matrix(x) else x)
out$qc_subject <- lapply(qc_subj, function(x) if (is.table(x)) as.data.frame.matrix(x) else x)
out$chi_vs_qc <- lapply(chi_by_grade, function(x) if (is.table(x)) as.data.frame.matrix(x) else x)
out$covariates_by_cell <- cbc

# ---------------------------------------------------------------- misc report checks
out$misc <- list(
  n_params_M1 = length(coef(fit)), df_check = nobs(fit) - length(coef(fit)),
  T1_occ0_in_total_SD = m1$est / g$pooled$sd_total_SNVTA,
  T1_wording_frozen = R$T1$Delta_vs_occasion0,
  R2s_PS_two_sided_p = k2$p_two,
  spec_curve = list(n_evaluable = R$spec_curve$n_evaluable, grade0_level_n = sum(first$qc_grade %in% 0),
                    note = "162 grade-0 specifications are non-evaluable by construction (grade-0 n < 30), so 324 = 486 - 162 and 'robust' can never be assessable"),
  PC4_gate_n = g$PC4_detail$n, resid_split_n = R$H1$B_coreport_residualised_split$n,
  threshold_psychosis = as.list(table(first$cell[first$dx_psychosis %in% 1])),
  stim_prescribed_by_cell = as.list(table(first$cell[first$stim_prescribed %in% 1])),
  n_ses1_first_usable = sum(first$session_id == "ses-1"), first_usable_session_by_cell = as.data.frame.matrix(table(first$cell, first$session_id)))

jsonlite::write_json(out, file.path(OUT, "review_rederivation.json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null", null = "null")
cat("written:", file.path(OUT, "review_rederivation.json"), "\n")
