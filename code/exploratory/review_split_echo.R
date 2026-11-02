# TASK R: split-echo reliability discrepancy check (EXPLORATORY / REVIEW ONLY)
suppressPackageStartupMessages({library(splines); library(jsonlite)})
d <- read.csv("/path/to/Penn_LEAD/LEADchi/mirror/unlocked/merged_qsmxt_primary_SNVTA.csv", stringsAsFactors = FALSE, na.strings = c("NA", ""))
d <- d[order(d$participant_id, d$ses_order), ]; u <- d[d$usable == 1, ]; f <- u[!duplicated(u$participant_id), ]
f$PS_e <- ifelse(f$study_group == "PRO/CHR", .5, -.5); f$ADHD_e <- ifelse(f$dx_adhd == 1, .5, -.5); f$sex_e <- ifelse(f$sex == "M", .5, -.5)
f$race3 <- factor(ifelse(f$race == "Black or African American", "Black", ifelse(f$race == "White", "White", "Other"))); contrasts(f$race3) <- contr.sum(3)
s13 <- f$chi_SNVTA_split13; s24 <- f$chi_SNVTA_split24; full <- f$chi_SNVTA
cat("summary split13:\n"); print(summary(s13)); cat("summary split24:\n"); print(summary(s24)); cat("summary full:\n"); print(summary(full))
o <- order(-abs(s24 - median(s24))); ex <- data.frame(sess = paste0(f$participant_id, "_", f$session_id)[o[1:8]], cell = paste(f$study_group, f$dx_adhd)[o[1:8]],
  full = full[o[1:8]], s13 = s13[o[1:8]], s24 = s24[o[1:8]], qc = f$qc_grade[o[1:8]], fimp = f$fimp_SNVTA[o[1:8]], Qphase = f$Q_phase[o[1:8]])
print(ex, digits = 5)
z24 <- (s24 - median(s24)) / mad(s24); z13 <- (s13 - median(s13)) / mad(s13)
cat("n |robust z| > 5: s24 =", sum(abs(z24) > 5), "; s13 =", sum(abs(z13) > 5), "\n")
cat("n outside envelope [-50,300]: s24 =", sum(s24 < -50 | s24 > 300), "; s13 =", sum(s13 < -50 | s13 > 300), "\n")
sb <- function(r) 2 * r / (1 + r); lo <- function(r, n) tanh(atanh(r) - 1.959964 / sqrt(n - 3))
rp <- function(keep, lab) { r <- cor(s13[keep], s24[keep]); n <- sum(keep)
  cat(sprintf("%-55s n=%3d  r=%.4f  R_SB=%.4f  rL=%.4f  spearman=%.4f\n", lab, n, r, sb(r), sb(lo(r, n)), cor(s13[keep], s24[keep], method = "spearman"))) }
rp(rep(TRUE, 80), "all 80 first-usable (as the frozen resid co-report used)")
rp(abs(z24) <= 5 & abs(z13) <= 5, "excluding |robust z| > 5 in either half")
rp(s24 >= -50 & s24 <= 300 & s13 >= -50 & s13 <= 300, "excluding envelope failures in either half")
k2 <- order(-abs(z24))[-(1:2)]; keep2 <- rep(FALSE, 80); keep2[k2] <- TRUE; rp(keep2, "dropping the 2 most extreme s24 sessions (n = 78)")
# residualised reliability after excluding the extreme sessions (what the registered co-report should have been)
fo <- ~ PS_e * ADHD_e + ns(age, df = 3) + sex_e + race3 + scale(roi_vol) + scale(Q_phase)
rr <- function(keep, lab) { g <- f[keep, ]; r1 <- resid(lm(update(fo, chi_SNVTA_split13 ~ .), data = g)); r2 <- resid(lm(update(fo, chi_SNVTA_split24 ~ .), data = g))
  r <- cor(r1, r2); n <- nrow(g); rL <- sb(lo(r, n)); dB <- 0.6 * sqrt(max(rL, 0))
  cat(sprintf("%-55s n=%3d  r_resid=%.4f  R_SB=%.4f  rL=%.4f  -> delta_B(resid)=%.3f\n", lab, n, r, sb(r), rL, dB)); c(r = r, rL = rL, dB = dB) }
a <- rr(rep(TRUE, 80), "residualised, all 80 (frozen co-report)")
b <- rr(abs(z24) <= 5 & abs(z13) <= 5, "residualised, excluding |z| > 5")
cc <- rr(keep2, "residualised, n = 78 (2 extreme dropped)")
# bound-test co-report with the cleaned residualised reliability
R <- jsonlite::read_json("/path/to/Penn_LEAD/LEADchi/mirror/results/confirmatory/results.json", simplifyVector = TRUE)
dec <- R$H1$decision
for (x in list(b, cc)) { tB <- (dec$beta_ppb + x["dB"] * dec$sigma_B) / dec$se_hc3; cat(sprintf("  bound at delta_B=%.3f: tB=%.3f pB=%.5f\n", x["dB"], tB, pt(tB, dec$df, lower.tail = FALSE))) }
# does the full chi of the extreme sessions look fine? (they are usable: visual QC passed)
cat("extreme sessions' full-chi robust z:", round(((full - median(full)) / mad(full))[o[1:2]], 2), "\n")
# correlation of each half with the full reconstruction, excluding extremes
cat("r(s13, full) / r(s24, full) excluding extremes:", round(cor(s13[keep2], full[keep2]), 3), "/", round(cor(s24[keep2], full[keep2]), 3), "\n")
# multiverse variants for the two extreme sessions (are the 9 reconstructions consistent for them?)
mv <- read.csv("/path/to/Penn_LEAD/LEADchi/mirror/unlocked/multiverse_qsmxt_SNVTA.csv", stringsAsFactors = FALSE)
for (s in o[1:2]) { m <- mv[mv$participant_id == f$participant_id[s] & mv$session_id == f$session_id[s], ]
  cat(f$participant_id[s], f$session_id[s], ": multiverse chi_wb_mean over 9 nodes: ", paste(round(m$chi_wb_mean, 1), collapse = " "), "\n") }
write.csv(ex, "/path/to/Penn_LEAD/LEADchi/results_exploratory/review/split_echo_extreme_sessions.csv", row.names = FALSE)

# ---- the six unsigned-12-bit phase sessions (amendment 1 #14): primary chi, split halves, multiverse, usable status
uns <- c("sub-20303_ses-3","sub-20712_ses-1","sub-20818_ses-1","sub-20823_ses-1","sub-20866_ses-1","sub-21038_ses-1")
d$sk <- paste0(d$participant_id, "_", d$session_id)
w <- d[d$sk %in% uns, c("sk","study_group","dx_adhd","usable","first_usable","fail_reason","qc_grade","chi_SNVTA","chi_SNVTA_split13","chi_SNVTA_split24","chi_GP","chi_STR","r2s_SNVTA")]
print(w, digits = 4)
for (s in uns) { m <- mv[paste0(mv$participant_id, "_", mv$session_id) == s, ]
  if (nrow(m)) cat(s, ": multiverse chi_wb_mean (9 nodes): ", paste(round(m$chi_wb_mean, 1), collapse = " "), "\n") }
# are any split24 envelope failures NOT in the unsigned list?
bad24 <- d$sk[!is.na(d$chi_SNVTA_split24) & (d$chi_SNVTA_split24 < -50 | d$chi_SNVTA_split24 > 300)]
bad13 <- d$sk[!is.na(d$chi_SNVTA_split13) & (d$chi_SNVTA_split13 < -50 | d$chi_SNVTA_split13 > 300)]
cat("split24 envelope failures (all 181 rows):", paste(bad24, collapse = ", "), "\n")
cat("split13 envelope failures (all 181 rows):", paste(bad13, collapse = ", "), "\n")
cat("of which unsigned-phase sessions:", paste(intersect(c(bad24, bad13), uns), collapse = ", "), "\n")
# exploratory M1 without the unsigned-phase sessions that are in first-usable (descriptive only)
fx <- f[!paste0(f$participant_id, "_", f$session_id) %in% uns, ]
fx$race3 <- factor(ifelse(fx$race == "Black or African American", "Black", ifelse(fx$race == "White", "White", "Other"))); contrasts(fx$race3) <- contr.sum(3)
fit <- lm(update(fo, chi_SNVTA ~ .), data = fx); V <- sandwich::vcovHC(fit, type = "HC3")
b <- coef(fit)["PS_e"]; se <- sqrt(V["PS_e","PS_e"]); cat(sprintf("M1 without the unsigned-phase sessions in H1: n=%d beta=%.3f se=%.3f p=%.4f\n", nobs(fit), b, se, 2*pt(-abs(b/se), df.residual(fit))))
cat("unsigned sessions in first-usable:", paste(intersect(paste0(f$participant_id, "_", f$session_id), uns), collapse = ", "), "\n")
