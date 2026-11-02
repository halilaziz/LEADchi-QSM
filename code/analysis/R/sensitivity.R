# LEAD-chi analysis: registered H1 sensitivity battery (sec. 24.4; Q13 sensitivity exclusions; Q14 puberty rule).
# Every row refits M1 (same decision function, same r_L / PC4) with ONE registered change and reports beta, SE, p
# and the M1 A/B/C outcome. None changes the primary decision. A row whose input column is absent from the merged
# table is reported as "not run (column absent)", never silently dropped.

# PDS-like puberty score (sec. 12.5): mean of the available Tanner self-report items at the SAME session, "I don't
# know" (code 5) -> NA; boys items 1-5 (height, skin, facial hair, body hair, voice); girls items 1, 2, 5, 6 (height,
# skin, body hair, breasts) plus menarche (item 3: yes -> 4, no -> 1); at least 3 items required.
pds_score <- function(d) {
  g <- function(v) if (v %in% names(d)) suppressWarnings(as.numeric(d[[v]])) else rep(NA_real_, nrow(d))
  dk <- function(x) ifelse(x %in% 1:4, x, NA_real_)
  boys <- cbind(dk(g("tanner_boy_1")), dk(g("tanner_boy_2")), dk(g("tanner_boy_3")), dk(g("tanner_boy_4")), dk(g("tanner_boy_5")))
  men <- g("tanner_girl_3"); menc <- ifelse(men %in% 1, 4, ifelse(men %in% 0, 1, NA_real_))
  girls <- cbind(dk(g("tanner_girl_1")), dk(g("tanner_girl_2")), dk(g("tanner_girl_5")), dk(g("tanner_girl_6")), menc)
  M <- ifelse(matrix(d$sex == "M", nrow(d), 5), boys, girls)
  k <- rowSums(!is.na(M))
  list(pds = ifelse(k >= 3, rowMeans(M, na.rm = TRUE), NA_real_),
       menarche_e = ifelse(d$sex == "F", ifelse(men %in% 1, 0.5, ifelse(men %in% 0, -0.5, NA_real_)), 0))
}

# Registered fallback for the puberty MI (mice is not installed; Q14): complete case, plus ONE stochastic regression
# imputation of the PDS from age, sex_e, PS_e, ADHD_e and the subject's nearest other-session PDS with its age gap
# (when available), seed-fixed.
impute_pds_once <- function(first, all, seed) {
  set.seed(seed)
  pa <- pds_score(all); all$pds <- pa$pds
  first$pds <- pds_score(first)$pds
  near <- vapply(seq_len(nrow(first)), function(i) {
    o <- all[all$sub == first$sub[i] & all$sess_key != first$sess_key[i] & !is.na(all$pds), ]
    if (!nrow(o)) return(NA_real_); o$pds[which.min(abs(o$age - first$age[i]))] }, numeric(1))
  gap <- vapply(seq_len(nrow(first)), function(i) {
    o <- all[all$sub == first$sub[i] & all$sess_key != first$sess_key[i] & !is.na(all$pds), ]
    if (!nrow(o)) return(NA_real_); (o$age - first$age[i])[which.min(abs(o$age - first$age[i]))] }, numeric(1))
  first$pds_near <- near; first$pds_gap <- gap
  obs <- !is.na(first$pds)
  use_near <- sum(obs & !is.na(near)) >= 20
  fo <- if (use_near) pds ~ age + sex_e + PS_e + ADHD_e + pds_near + pds_gap else pds ~ age + sex_e + PS_e + ADHD_e
  fit <- stats::lm(fo, data = first[obs, ])
  mis <- which(!obs & (if (use_near) !is.na(near) else TRUE))
  imp <- first$pds
  if (length(mis)) imp[mis] <- stats::predict(fit, newdata = first[mis, ]) + stats::rnorm(length(mis), 0, stats::sigma(fit))
  list(pds_imputed = imp, n_observed = sum(obs), n_imputed = length(mis), model = deparse(fo))
}

sens_row <- function(name, d, rL, pc4, y = "chi_SNVTA", extra = "", note = "") {
  r <- try_fit({
    d <- set_race3(d); fo <- m1_formula(y, d, extra); m <- fit_ols_hc3(fo, d)
    hb <- huber_fit(fo, d, "PS_e", ols_resid = stats::residuals(m$fit)); dec <- h1_decide(m, hb$scale, rL, pc4)
    data.frame(sensitivity = name, n = m$n, beta_ppb = dec$beta_ppb, se_hc3 = dec$se_hc3, p_two = dec$p_two,
               pB = dec$B_registered$pB, outcome = dec$outcome, note = note, stringsAsFactors = FALSE)
  }, paste("sensitivity", name))
  if (is.null(r)) data.frame(sensitivity = name, n = nrow(d), beta_ppb = NA, se_hc3 = NA, p_two = NA, pB = NA,
                             outcome = "fit failed", note = note, stringsAsFactors = FALSE) else r
}
not_run <- function(name, why) data.frame(sensitivity = name, n = NA, beta_ppb = NA, se_hc3 = NA, p_two = NA, pB = NA,
                                          outcome = "not run", note = why, stringsAsFactors = FALSE)
has <- function(d, v) v %in% names(d) && any(!is.na(d[[v]]))

# First-usable sample re-derived when some exclusion reasons are ignored: sessions whose ONLY failure reasons lie in
# `ignore` become usable (e.g. (c') plausibility; or the visual criteria c_qc_grade and d_placement).
first_relaxed <- function(all, ignore) {
  if (!"fail_reason" %in% names(all)) return(NULL)
  relax <- vapply(strsplit(ifelse(is.na(all$fail_reason), "", all$fail_reason), ";"),
                  function(z) length(z) > 0 && all(z %in% ignore), logical(1))
  a <- all; a$usable[relax] <- 1
  a <- a[order(a$participant_id, a$ses_order), ]; u <- a[a$usable == 1 & !is.na(a$chi_SNVTA), ]
  u[!duplicated(u$participant_id), ]
}
first_without_cprime <- function(all) first_relaxed(all, "cprime_plausibility")

run_h1_sensitivities <- function(S, gates, seed, other_engine = NULL) {
  f <- S$first; a <- S$all; rL <- gates$rL_primary; pc4 <- isTRUE(gates$PC4)
  rows <- list(); add <- function(x) rows[[length(rows) + 1]] <<- x
  # ---- sample exclusions (Q13) ----
  add(if (has(f, "dx_psychosis")) sens_row("without the 6 threshold-psychosis youths", f[!(f$dx_psychosis %in% 1), ], rL, pc4) else
        not_run("without threshold psychosis", "dx_psychosis absent"))
  add(sens_row("without 4 TD youths (dx-coding gap / prescribed stimulants)", f[!f$participant_id %in% REG$td_coding_gap, ], rL, pc4))
  add(if (has(f, "stim_prescribed")) sens_row("without prescribed-stimulant endorsers (same session)", f[!(f$stim_prescribed %in% 1), ], rL, pc4) else
        not_run("without stimulant endorsers", "stim_prescribed absent"))
  add(if (has(f, "t1_artifact")) sens_row("without T1w 'artifact' sessions", f[!(f$t1_artifact %in% 1), ], rL, pc4) else
        not_run("without T1w artifact sessions", "t1_artifact absent"))
  add(sens_row("without sub-21568 (no T1-QC entry)", f[f$participant_id != REG$t1qc_missing_sub, ], rL, pc4))
  add(if (has(f, "epoch")) sens_row("without early-epoch sessions", f[f$epoch != "early", ], rL, pc4) else
        not_run("without early-epoch sessions", "epoch absent"))
  add(if (has(f, "age_ses1")) sens_row("baseline age >= 11 only", f[!is.na(f$age_ses1) & f$age_ses1 >= 11, ], rL, pc4) else
        not_run("baseline age >= 11", "age_ses1 absent"))
  add(sens_row("exposure-sensitivity set (4 inspected sessions dropped; Q18)", S$exposure, rL, pc4))
  fwc <- first_without_cprime(a)
  add(if (!is.null(fwc)) sens_row("without the (c') plausibility exclusions", fwc, rL, pc4) else
        not_run("without (c') exclusions", "fail_reason absent"))
  add(if (has(f, "chi_SNVTA_impkeep")) sens_row("implausible voxels retained", f[!is.na(f$chi_SNVTA_impkeep), ], rL, pc4, y = "chi_SNVTA_impkeep") else
        not_run("implausible voxels retained", "chi_SNVTA_impkeep absent"))
  add(if (has(f, "dental_D")) sens_row("without rubric code D (dental/metal) sessions", f[!(f$dental_D %in% 1), ], rL, pc4) else
        not_run("without code-D sessions", "dental_D absent"))
  add(if (has(f, "qc_grade")) sens_row("strict QC (visual grade 0 only)", f[f$qc_grade %in% 0, ], rL, pc4) else
        not_run("strict QC (grade 0)", "qc_grade absent"))
  fau <- first_relaxed(a, c("c_qc_grade", "d_placement"))
  add(if (!is.null(fau)) sens_row("automated exclusions only (no visual grade or placement rating)", fau, rL, pc4) else
        not_run("automated exclusions only", "fail_reason absent"))
  # ---- model and covariates ----
  add(try_fit({ d <- set_race3(f); fo <- stats::as.formula(sub(" + race3", "", paste("chi_SNVTA ~", m1_rhs(d)), fixed = TRUE))
    m <- fit_ols_hc3(fo, d); hb <- huber_fit(fo, d, "PS_e", stats::residuals(m$fit)); dec <- h1_decide(m, hb$scale, rL, pc4)
    data.frame(sensitivity = "without race", n = m$n, beta_ppb = dec$beta_ppb, se_hc3 = dec$se_hc3, p_two = dec$p_two,
               pB = dec$B_registered$pB, outcome = dec$outcome, note = "race3 removed", stringsAsFactors = FALSE) }, "sens no race"))
  if (has(f, "ethnicity")) { f$hisp_e <- ifelse(grepl("^Hispanic", f$ethnicity), 0.5, -0.5)
    add(sens_row("+ Hispanic ethnicity", f, rL, pc4, extra = " + hisp_e")) } else add(not_run("+ Hispanic ethnicity", "ethnicity absent"))
  ps <- pds_score(f); f$pds <- ps$pds; f$menarche_e <- ps$menarche_e
  if (any(!is.na(f$pds))) {
    add(sens_row("+ puberty (PDS), complete case", f[!is.na(f$pds), ], rL, pc4, extra = " + scale(pds)"))
    imp <- impute_pds_once(f, a, seed + 71); f$pds_imp <- imp$pds_imputed
    add(sens_row("+ puberty (PDS), one stochastic regression imputation (mice fallback)", f[!is.na(f$pds_imp), ], rL, pc4,
                 extra = " + scale(pds_imp)", note = paste0("observed ", imp$n_observed, ", imputed ", imp$n_imputed)))
    add(sens_row("+ menarche (girls; 0 for boys)", f[!is.na(f$menarche_e), ], rL, pc4, extra = " + menarche_e"))
  } else add(not_run("+ puberty / menarche", "Tanner items absent"))
  add(if (has(f, "icv_seg_ml")) sens_row("+ intracranial volume (SynthSeg, coded side)", f[!is.na(f$icv_seg_ml), ], rL, pc4, extra = " + scale(icv_seg_ml)") else
        not_run("+ ICV", "icv_seg_ml absent"))
  if (has(f, "dx_mdd") || has(f, "dx_moodnos")) { f$mood_e <- ifelse(f$dx_mdd %in% 1 | f$dx_moodnos %in% 1, 0.5, -0.5)
    add(sens_row("+ mood diagnosis (MDD or mood NOS)", f, rL, pc4, extra = " + mood_e")) } else add(not_run("+ mood diagnosis", "dx flags absent"))
  if (has(f, "epoch")) { f$epoch_e <- ifelse(f$epoch == "early", 0.5, -0.5)
    add(if (length(unique(f$epoch_e)) > 1) sens_row("+ protocol epoch", f, rL, pc4, extra = " + epoch_e") else not_run("+ epoch", "one epoch only")) } else
    add(not_run("+ protocol epoch", "epoch absent"))
  if (has(f, "f0_cluster")) { f$f0c <- factor(f$f0_cluster)
    add(if (nlevels(f$f0c) > 1) sens_row("+ f0 cluster", f, rL, pc4, extra = " + f0c") else not_run("+ f0 cluster", "one cluster only")) } else
    add(not_run("+ f0 cluster", "f0_cluster absent"))
  if (has(f, "age_ses1")) { f$yrs_dx <- f$age - f$age_ses1
    add(sens_row("+ years since the diagnostic assessment", f[!is.na(f$yrs_dx), ], rL, pc4, extra = " + yrs_dx")) } else add(not_run("+ years since diagnosis", "age_ses1 absent"))
  add(if (has(f, "abs_pitch")) sens_row("+ |pitch|", f[!is.na(f$abs_pitch), ], rL, pc4, extra = " + scale(abs_pitch)") else not_run("+ |pitch|", "abs_pitch absent"))
  add(if (has(f, "dwi_FD")) sens_row("+ DWI mean FD (same session)", f[!is.na(f$dwi_FD), ], rL, pc4, extra = " + scale(dwi_FD)") else not_run("+ DWI FD", "dwi_FD absent"))
  add(if (has(f, "dice_min")) sens_row("+ registration Dice (minimum of the G5 Dice values)", f[!is.na(f$dice_min), ], rL, pc4, extra = " + scale(dice_min)") else
        not_run("+ registration Dice", "dice_min absent"))
  add(if (has(f, "motion_gsr")) sens_row("+ magnitude ghost-to-signal motion index", f[!is.na(f$motion_gsr), ], rL, pc4, extra = " + scale(motion_gsr)") else
        not_run("+ motion index", "motion_gsr absent"))
  frac <- intersect(c("ref_frac_wm", "ref_frac_csf", "fov_brain_frac"), names(f))
  add(if (length(frac) == 3 && all(vapply(frac, function(v) any(!is.na(f[[v]])), TRUE)))
    sens_row("+ reference-domain WM and CSF fractions and brain-in-FOV fraction", f[stats::complete.cases(f[, frac]), ], rL, pc4,
             extra = " + scale(ref_frac_wm) + scale(ref_frac_csf) + scale(fov_brain_frac)") else not_run("+ reference-domain fractions", "columns absent"))
  add(if (has(f, "chi_SNVTA_harm")) sens_row("coverage-harmonised whole-brain reference", f[!is.na(f$chi_SNVTA_harm), ], rL, pc4, y = "chi_SNVTA_harm") else
        not_run("coverage-harmonised reference", "chi_SNVTA_harm absent"))
  # sub-regions (descriptive; the medial SN-VTA is exposed to flow and venous partial volume without flow compensation)
  add(if (has(f, "chi_SNonly")) sens_row("SN-only complex (descriptive)", f[!is.na(f$chi_SNonly), ], rL, pc4, y = "chi_SNonly") else
        not_run("SN-only complex", "chi_SNonly absent"))
  add(if (has(f, "chi_VTAPBP")) sens_row("VTA + PBP only (descriptive)", f[!is.na(f$chi_VTAPBP), ], rL, pc4, y = "chi_VTAPBP") else
        not_run("VTA + PBP only", "chi_VTAPBP absent"))
  # rank-based ANCOVA (outcome ranks; HC3); its B test is not meaningful and is not reported
  rk <- try_fit({ d <- set_race3(f); d$chi_rank <- rank(d$chi_SNVTA); m <- fit_ols_hc3(m1_formula("chi_rank", d), d); t <- term_test(m, "PS_e")
    data.frame(sensitivity = "rank-based ANCOVA (ranks of chi_SNVTA)", n = m$n, beta_ppb = NA, se_hc3 = NA, p_two = t$p_two, pB = NA,
               outcome = if (t$est < 0 && t$p_two < REG$alpha) "rank test p < .05 (lower in PS)" else "rank test not significant",
               note = "beta in rank units not shown", stringsAsFactors = FALSE) }, "rank ANCOVA")
  add(rk)
  # other engine (only if A2 was approved and the other engine's merged table is supplied)
  if (!is.null(other_engine)) {
    o <- other_engine; o <- o[order(o$participant_id, o$ses_order), ]; ou <- o[o$usable == 1, ]; of <- ou[!duplicated(ou$participant_id), ]
    add(sens_row(paste0("other engine (", unique(o$engine), ")"), of, rL, pc4))
  } else add(not_run("other engine", "no other-engine table supplied (A2 not approved or not run)"))
  out <- do.call(rbind, rows[!vapply(rows, is.null, TRUE)])
  list(table = out)
}
