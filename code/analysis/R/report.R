# LEAD-chi analysis: human-readable report (Markdown). Every number also goes to results.json.
# The "Registered reading" section prints ONLY sentences registered verbatim in sec. 4.6 of the registration.

md_table <- function(df) {
  if (is.null(df) || !nrow(df)) return("(none)\n")
  h <- paste0("| ", paste(names(df), collapse = " | "), " |")
  s <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  r <- apply(df, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))
  paste(c(h, s, r, ""), collapse = "\n")
}

REG_SENT2 <- list(
  dilution = paste("If a trait deficit existed only in the minority who will convert to psychosis, the group-level effect",
                   "would be about 0.09-0.18 SD; such a deficit can be neither detected nor excluded here."),
  T1_not_declined = paste("Within the same youths, the psychosis-spectrum SN-VTA susceptibility difference did not decline by",
                          "0.4 total SD or more over about two years. This is within-person stability, not a replication."),
  T1_declined = paste("Within the same youths, the psychosis-spectrum SN-VTA susceptibility difference declined over about",
                      "two years. This is within-person change, not a failed replication."),
  T1_indeterminate = "The temporal-stability test (T1) was indeterminate.",
  T1_not_evaluated = "T1 was not evaluated because outcome A does not hold; follow-up estimates are descriptive.",
  S1_dissociation = "Dissociation: the PS effect is relatively nigral and the ADHD effect relatively striatal.",
  S1_ps = "The PS effect is relatively more nigral than striatal.",
  S1_adhd = "The ADHD effect is relatively more striatal than nigral.",
  S1_contrast = "The topographic contrast L is significant; neither component is significant on its own.",
  S1_null = "No topographic dissociation was established (S1 not significant under Holm); S1 was powered only for large contrasts.",
  S2_sig = "SN-VTA chi covaries with PRIME, a symptom dimension that partly indexes PS status.",
  S2_beyond = "The PRIME association holds beyond PS diagnosis.",
  S2_null = "No association between SN-VTA chi and PRIME was established (S2 not significant under Holm).",
  S3_coupled = "Right dentate chi is coupled to right proximal superior cerebellar peduncle tissue-weighted neurite density (coupled; not causal).",
  S3_specific = "The coupling is specific to the dentate-to-SCP pair (both negative controls excluded).",
  S3_null = "S3 is adequately powered only for strong latent coupling; a null S3 does not exclude moderate coupling.",
  S3_exploratory = "S3 is exploratory because a dentate gate or PC5 failed.",
  opposite = "A Holm-significant result in the direction opposite to the prediction is reported as such and not claimed.")

write_report <- function(R, path) {
  L <- character(0); a <- function(...) L <<- c(L, paste0(...))
  h <- R$H1; d <- h$decision
  if (R$meta$data_origin == "synthetic")
    a("> **SYNTHETIC OUTCOMES. NOT RESULTS.** chi, R2* and all imaging covariates in this table were simulated on the real",
      " Penn LEAD design (scenario: ", R$meta$scenario, "). Nothing here is an estimate of any real quantity.\n")
  a("# LEAD-chi confirmatory analysis report\n")
  a("- Code: ", REG$code_version, "; run ", R$meta$started_utc, " (UTC); ", R$meta$r_version)
  a("- Environment pin: ", R$meta$version_check)
  a("- Input: `", R$meta$input, "` (SHA-256 `", R$meta$input_sha256, "`)")
  a("- Gates: `", R$meta$gates, "` (SHA-256 `", R$meta$gates_sha256, "`)")
  a("- Branch from gates: engine **", R$branch$engine, "**, pipeline **", R$branch$pipeline, "**, ROI **", R$branch$roi,
    "**; primary interpretable: **", R$branch$primary_interpretable, "**", if (length(R$branch$caveats)) paste0("; caveats: ", paste(R$branch$caveats, collapse = "; ")) else "")
  a("- Resampling counts: ", paste(names(R$meta$draws), unlist(R$meta$draws), sep = " = ", collapse = ", "),
    if (!R$meta$draws_registered) " **(NON-REGISTERED counts: demonstration only)**" else " (registered)")
  a("")
  a("## 1. Samples (TD / ADHD-only / PS-only / PS+ADHD)\n"); a(md_table(R$samples))
  a("QC-by-group (subject level; Fisher's exact): not in H1 p = ", fp(R$qc$fisher_p_in_H1), "; earliest MEGRE session failed p = ",
    fp(R$qc$fisher_p_earliest), if (!is.null(R$qc$ipw)) paste0("; IPW-M1 beta_PS = ", f3(R$qc$ipw$est), " (p = ", fp(R$qc$ipw$p_two), ")") else "; IPW-M1 not triggered", "\n")
  a("## 2. Plausibility re-verification (sec. 12.4)\n")
  a(md_table(R$plausibility$log)); a("Label-free robust z (first-usable; |z| > 5 flagged, retained):\n"); a(md_table(R$plausibility$robust_z))
  a("## 3. H1 (M1; HC3; two-sided alpha .05)\n")
  a("`", h$formula, "`, n = ", h$n, "\n")
  a("| Quantity | Value |\n|---|---|")
  a("| beta_PS (ppb) [95% CI]; [90% CI] | ", f3(d$beta_ppb), " [", f3(d$ci95_ppb[1]), ", ", f3(d$ci95_ppb[2]), "]; [", f3(d$ci90_ppb[1]), ", ", f3(d$ci90_ppb[2]), "] |")
  a("| SE_HC3; t; df; two-sided p | ", f3(d$se_hc3), "; ", f3(d$t), "; ", d$df, "; ", fp(d$p_two), " |")
  a("| sigma_resid; sigma_rob (", h$huber$scale_source, "); sigma_B = min | ", f3(d$sigma_resid), "; ", f3(d$sigma_rob), "; ", f3(d$sigma_B), " |")
  a("| delta = beta / sigma_resid [cell-stratified bootstrap 95% CI] | ", f3(d$delta), " [", f3(h$boot_delta$ci95[1]), ", ", f3(h$boot_delta$ci95[2]), "] |")
  a("| noncentral-t 95% CI of delta; raw Hedges g | [", f3(h$nct_ci_delta[1]), ", ", f3(h$nct_ci_delta[2]), "]; ", f3(h$hedges_g_raw), " |")
  a("| B: r_L; B evaluable; delta_B; bound in ppb (delta_B x sigma_B) | ", f3(d$rL), "; ", d$B_evaluable, "; ", f3(d$B_registered$delta_B), "; ", f3(h$ppb$bound_ppb), " |")
  a("| t_B; one-sided p_B | ", f3(d$B_registered$tB), "; ", fp(d$B_registered$pB), " |")
  a("| B co-reports (never decisive): bound -0.6; -0.5; cross-session bound; residualised split-echo bound | ", fp(d$B_coreport[[1]]$pB), "; ", fp(d$B_coreport[[2]]$pB), "; ",
    fp(h$B_coreport_cross_session$pB), "; ", fp(h$B_coreport_residualised_split$pB), " |")
  a("| Freedman-Lane p (", h$freedman_lane$N, " permutations within ADHD strata) [MC 95% CI] | ", fp(h$freedman_lane$p), " [",
    f3(h$freedman_lane$mc_ci95[1], 4), ", ", f3(h$freedman_lane$mc_ci95[2], 4), "] |")
  a("| Huber M-estimate of beta_PS (SE) | ", f3(h$huber$est), " (", f3(h$huber$se), ") |")
  for (k in names(h$bayes$BF)) a("| BF10 under ", h$bayes$BF[[k]]$prior, " (", k, "); BF01 | ", f3(h$bayes$BF[[k]]$BF10), "; ", f3(h$bayes$BF[[k]]$BF01), " |")
  a("| sceptical-prior shrunken delta [95% CrI] | ", f3(h$bayes$shrunken_sceptical$mean), " [", f3(h$bayes$shrunken_sceptical$ci95[1]), ", ",
    f3(h$bayes$shrunken_sceptical$ci95[2]), "] |")
  a("| adult reference (Vano 2025; descriptive; reference, resolution, field and pipeline differ) | ", REG$adult_ppb, " ppb |")
  a("")
  a("**M1 outcome (whole-brain reference): ", d$outcome, "**", if (d$disconfirming) " (**disconfirming: beta > 0, p < .05**)" else "", "\n")
  a("### K3 reference agreement (CSF-only secondary reference)\n")
  a("beta_PS(dref) = ", f3(h$K3$beta_dref), " ppb (p ", fp(h$K3$p_two), "); threshold max(0.25|beta|, 0.10 sigma) = ", f3(h$K3$threshold),
    " and p < ", h$K3$p_threshold, "; trigger: **", h$K3$trigger, "**. Sessions without a valid CSF reference: ", h$K3$n_no_csf,
    ". CSF-referenced outcome: ", h$K3$csf_outcome, " (beta ", f3(h$K3$csf_beta), ", p ", fp(h$K3$csf_p), ", p_B ", fp(h$K3$csf_pB), "). Rule: ", h$K3$rule, ".\n")
  a("**Final H1 outcome (M1 -> K3): ", h$final_outcome, ".**\n")
  a("### K1 specificity labels (precedence K1c > K1b > K1a)\n")
  q <- h$K1
  a(md_table(data.frame(contrast = c("PS-only - TD (mandatory co-report)", "PS+ADHD - ADHD-only", "ADHD main effect", "PS x ADHD interaction"),
                        est_ppb = f3(c(q$simple_PSonly_minus_TD$est, q$simple_PSADHD_minus_ADHDonly$est, q$ADHD$est, q$interaction$est)),
                        se = f3(c(q$simple_PSonly_minus_TD$se, q$simple_PSADHD_minus_ADHDonly$se, q$ADHD$se, q$interaction$se)),
                        p_two = sapply(c(q$simple_PSonly_minus_TD$p_two, q$simple_PSADHD_minus_ADHDonly$p_two, q$ADHD$p_two, q$interaction$p_two), fp))))
  a("ADHD TOST at +/-", REG$tost_bound, " sigma_B (+/-", f3(q$tost$bound_ppb), " ppb): p = ", fp(q$tost$p_lower), " / ", fp(q$tost$p_upper),
    "; success: ", q$tost$success, ". Interaction 95% CI includes 0: ", q$interaction_ci_includes_0, " (non-rejection only).\n")
  a("**K1 label: ", q$label, "**. Mandatory caveat: all 6 threshold-psychosis youths sit in PS+ADHD.\n")
  a("### K2 iron vs myelin\n")
  a("R2* beta_PS = ", f3(h$K2$est), " s^-1 (n = ", h$K2$n, "), one-sided p = ", fp(h$K2$p_one_sided_lower), ".\n")
  a("### K4 precision, influence, fragility and robustness\n")
  a("- PC4: ", h$precision$PC4_pass, "; r_L ", f3(h$precision$rL), "; B evaluable: ", h$precision$B_evaluable, "; label: ", h$precision$label,
    "; attenuation-corrected delta = ", f3(h$precision$delta_attenuation_corrected), ".")
  di <- h$influence$diagnostics
  a("- Dominance: top-1 share ", f3(di$top1_share), ", top-2 share ", f3(di$top2_share), ", max |t_i| ", f3(di$max_abs_rstudent),
    ", max Cook's D ", f3(di$max_cooks_d), "; triggered: ", di$triggered,
    if (di$triggered) paste0(" (refit outcome ", h$influence$refit$outcome, ")") else "", "; **influence-dependent: ", h$influence$influence_dependent, "**.")
  a("- Leave-one-out (", h$loo$compared, "): ", h$loo$n_category_change, " of ", h$loo$n_refits, " deletions change the outcome; **fragile: ", h$loo$fragile, "**.")
  a("- ses-1-only (mandatory co-report): n = ", h$ses1_only$n, ", beta ", f3(h$ses1_only$beta), " ppb, p ", fp(h$ses1_only$p), ", outcome ", h$ses1_only$outcome, ".")
  a("- Exposure sensitivity (Q18(4); ", h$exposure_sensitivity$n_dropped, " of the 4 listed sessions in this sample",
    if (!is.null(h$exposure_sensitivity$dropped)) paste0(": ", paste(h$exposure_sensitivity$dropped, collapse = ", ")) else "",
    "): n = ", h$exposure_sensitivity$n, ", beta ", f3(h$exposure_sensitivity$beta), ", p ", fp(h$exposure_sensitivity$p), ", outcome ",
    h$exposure_sensitivity$outcome, "; changes conclusion: ", h$exposure_sensitivity$changes_conclusion, " (flagged, never overrides).")
  sc <- R$spec_curve
  if (!is.null(sc$n_evaluable))
    a("- Specification curve: ", sc$n_evaluable, " of ", sc$n_specifications, " evaluable; median delta ", f3(sc$median_delta), "; share delta < 0 ",
      f3(sc$share_delta_negative), "; share significant ", f3(sc$share_significant_predicted), "; Stouffer Z ", f3(sc$stouffer_Z),
      "; joint permutation p (median t, ", sc$joint$N, ") = ", fp(sc$joint$p), "; **robust: ", sc$robust, "** (", sc$robust_reason, ").")
  else a("- Specification curve: ", sc$note)
  a("- PC2x (age -> putamen / GP chi, one-sided): ", paste(vapply(R$PC2x$tests, function(t) paste0(t$roi, " slope ", f3(t$slope_ppb_per_y), " ppb/y p ", fp(t$p_one_sided)), ""), collapse = "; "),
    "; pass: ", R$PC2x$pass, ".\n")
  a("### H1 sensitivity battery (sec. 24.4; reported separately, never decisive)\n")
  st <- R$sensitivity$table
  a(md_table(transform(st, beta_ppb = f3(beta_ppb), se_hc3 = f3(se_hc3), p_two = sapply(p_two, fp), pB = sapply(pB, fp))))
  a("## 4. Key secondary (Holm; family: ", paste(R$holm$members, collapse = ", "), "; S3 status: ", R$holm$S3_status, ")\n")
  a(md_table(data.frame(test = c("S1 L (SN-VTA vs STR x PS vs ADHD)", "S2 beta_PRIME (per SD)", "S3a beta(chi_DN_R) (standardised)"),
                        n = c(R$S1$n, R$S2$n, R$S3$n),
                        est = f3(c(R$S1$L$est, R$S2$beta_PRIME$est, R$S3$beta_DN$est)),
                        se_hc3 = f3(c(R$S1$L$se, R$S2$beta_PRIME$se, R$S3$beta_DN$se)),
                        p_two = sapply(c(R$S1$L$p_two, R$S2$beta_PRIME$p_two, R$S3$beta_DN$p_two), fp),
                        p_holm = sapply(c("S1", "S2", "S3a"), function(k) if (k %in% R$holm$members) fp(R$holm$p_holm[[k]]) else "not in family"),
                        supported = sapply(c("S1", "S2", "S3a"), function(k) if (k %in% R$holm$members) R$holm$supported[[k]] else NA),
                        influence_dependent = sapply(c("S1", "S2", "S3a"), function(k) R$sec_influence[[k]]$influence_dependent))))
  a("- **S1** pooled constants: ", R$S1$pooled_source, ". beta_PS(D) = ", f3(R$S1$beta_PS_D$est), " (one-sided p ", fp(R$S1$beta_PS_D$p_lower),
    "); beta_ADHD(D) = ", f3(R$S1$beta_ADHD_D$est), " (one-sided p ", fp(R$S1$beta_ADHD_D$p_upper), "). Wording: **", R$S1_labels$wording, "**.",
    " Stacked GLS L = ", f3(R$S1$stacked_gls$L), " (p ", fp(R$S1$stacked_gls$p_two), ")",
    if (!is.null(R$S1$robust_scaled_L)) paste0("; robust-scaled L = ", f3(R$S1$robust_scaled_L$est), " (p ", fp(R$S1$robust_scaled_L$p_two), ")") else "", ".")
  a("- **S2** + PS_e: beta_PRIME = ", f3(R$S2$plus_PS_e$est), " (p ", fp(R$S2$plus_PS_e$p_two), "); 'beyond diagnosis' allowed: ",
    R$S2$beyond_diagnosis_wording_allowed, ". Within PS: ", f3(R$S2$within_PS$est), " (p ", fp(R$S2$within_PS$p_two), "); rank-INT: ",
    f3(R$S2$rank_inverse_normal$est), " (p ", fp(R$S2$rank_inverse_normal$p_two), ").")
  for (k in names(R$S2$sensitivities)) { s <- R$S2$sensitivities[[k]]; a("  - S2 sensitivity ", k, ": beta ", f3(s$est), " (p ", fp(s$p_two), ")") }
  if (!is.null(R$S2_allsession$between)) a("  - S2 all-session between-person PRIME: ", f3(R$S2_allsession$between$est), " (SE ", f3(R$S2_allsession$between$se),
                                          "); within-person ", f3(R$S2_allsession$within$est), " (SE ", f3(R$S2_allsession$within$se), ")")
  nc <- R$S3$negative_controls
  a("- **S3** confirmatory gates (PC4-DN, G2-DN, PC5): ", R$S3$confirmatory, ". Freedman-Lane p = ", fp(R$S3$freedman_lane$p),
    ". NC1 (DN - PUT) = ", f3(nc$NC1_DN_minus_PUT$obs), ", one-sided 95% lower bound ", f3(nc$NC1_DN_minus_PUT$lower95_one_sided),
    "; NC2 (pSCP - MCP matched) = ", f3(nc$NC2_pSCP_minus_MCPmatched$obs), ", lower bound ", f3(nc$NC2_pSCP_minus_MCPmatched$lower95_one_sided),
    ". Wording: **", R$S3_labels$wording, "**.")
  for (k in names(R$S3$sensitivities)) { s <- R$S3$sensitivities[[k]]
    a("  - S3 sensitivity ", k, ": beta ", f3(s$est), " (p ", fp(s$p_two), ")", if (!is.null(s$label)) paste0("; reduction ", f3(s$reduction), "; label: ", s$label) else "") }
  a("")
  t1 <- R$T1
  a("## 5. T1 temporal stability (held-out follow-ups; Delta = ", REG$t1_margin, " total SD = ", f3(t1$margin$Delta_ppb), " ppb)\n")
  a("REML maximum: fit `", t1$reml$best, "`; ICC ", f3(t1$reml$icc), "; boundary: ", t1$reml$boundary,
    "; profile confirms the global maximum on the grid: ", t1$reml$profile$best_is_global_on_grid, ".\n")
  a(md_table(transform(t1$reml$table, theta = f3(theta, 5), reml_deviance = f3(reml_deviance, 6), gap_to_best = f3(gap_to_best, 8))))
  a("| Quantity | Value |\n|---|---|")
  a("| rows / subjects / returners | ", t1$n_rows, " / ", t1$n_subjects, " / ", t1$n_returners, " |")
  a("| beta(PS x occ) (total-SD units), SE, containment df | ", f3(t1$beta_PSxocc), ", ", f3(t1$se), ", ", t1$containment_df, " |")
  a("| margin test H0: beta >= Delta: decision p (containment df); PB p (co-report, ", t1$margin$pb$N, ") | ", fp(t1$margin$p_decision), "; ", fp(t1$margin$pb$p), " |")
  a("| decline test H0: beta = 0: decision p (containment df); PB p (co-report) | ", fp(t1$decline$p_decision), "; ", fp(t1$decline$pb$p), " |")
  a("| OLS cluster-robust (HC2; not CR2 for the LMM): beta; SE; wild cluster p (margin, ", t1$cluster_robust$wild_margin$N, ") | ", f3(t1$cluster_robust$ols_beta), "; ",
    f3(t1$cluster_robust$se_vcovCL_HC2), "; ", fp(t1$cluster_robust$wild_margin$p), " |")
  a("| occasion-0 beta_PS (total-SD units) vs Delta | ", f3(t1$occasion0_beta_tsd), ": ", t1$Delta_vs_occasion0, " |")
  a("| follow-up-only M1 beta_PS (ppb), p, n | ", f3(t1$followup_only$beta_PS), ", ", fp(t1$followup_only$p_two), ", ", t1$followup_only$n, " |")
  a("| exposure sensitivity beta; p margin | ", f3(t1$exposure_sensitivity$beta), "; ", fp(t1$exposure_sensitivity$p_margin_containment), " |")
  a("")
  a("**T1 label: ", t1$label, "**. Assumption: ", t1$assumption, ".\n")
  a("## 6. Pooled all-session sensitivity (not a test)\n")
  a("beta_PS = ", f3(R$pooled$beta_PS), " ppb (SE ", f3(R$pooled$se), ", containment df ", R$pooled$containment_df, ", p ", fp(R$pooled$p_containment),
    "); REML fit `", R$pooled$reml$best, "`; OLS cluster-robust (HC2) SE ", f3(R$pooled$ols_cluster_robust$se_vcovCL_HC2), "; wild cluster p ", fp(R$pooled$ols_cluster_robust$wild$p), ".\n")
  a("## 7. Registered reading (verbatim registered sentences only)\n"); a(paste0("- ", R$reading, collapse = "\n"), "\n")
  a("## 8. Log\n"); a("```"); a(.LOG$lines); a("```")
  writeLines(L, path)
}

# Registered reading (sec. 4.6): a list of registered sentences selected by the registered rules.
registered_reading <- function(R) {
  h <- R$H1; o <- h$final_outcome; out <- character(0); add <- function(s) out <<- c(out, s)
  if (!R$branch$primary_interpretable) add("Pipeline not validated (positive-control or registration gate failure): H1 is reported descriptively only.")
  add(outcome_sentence(o, h$decision$disconfirming))
  if (o %in% c("B", "C")) add(REG_SENT2$dilution)
  if (!h$precision$B_evaluable) add(REGISTERED_SENTENCES$low_precision)
  if (o %in% c("A", "A+B") && !h$decision$disconfirming) {
    lab <- h$K1$label
    if (startsWith(lab, "K1c")) add(REGISTERED_SENTENCES$K1c) else if (startsWith(lab, "K1b")) add(REGISTERED_SENTENCES$K1b) else
      if (startsWith(lab, "K1a")) add(REGISTERED_SENTENCES$K1a)
    add(h$K2$wording)
    t1 <- R$T1$label
    add(if (t1 == "did not decline by >= Delta") REG_SENT2$T1_not_declined else if (t1 == "declined") REG_SENT2$T1_declined else REG_SENT2$T1_indeterminate)
    if (isTRUE(R$spec_curve$robust)) add(REGISTERED_SENTENCES$robust)
  } else add(REG_SENT2$T1_not_evaluated)
  if (isTRUE(h$loo$fragile)) add(REGISTERED_SENTENCES$fragile)
  if (isTRUE(h$influence$influence_dependent)) add(REGISTERED_SENTENCES$influence)
  # key secondary
  sl <- R$S1_labels$wording
  add(switch(sl, "dissociation" = REG_SENT2$S1_dissociation, "the PS effect is relatively more nigral" = REG_SENT2$S1_ps,
             "the ADHD effect is relatively more striatal" = REG_SENT2$S1_adhd,
             "L significant; neither component significant (worded as the contrast only)" = REG_SENT2$S1_contrast, REG_SENT2$S1_null))
  if (isTRUE(R$holm$supported$S2)) { add(REG_SENT2$S2_sig); if (isTRUE(R$S2$beyond_diagnosis_wording_allowed)) add(REG_SENT2$S2_beyond) } else add(REG_SENT2$S2_null)
  if (!R$holm$S3_in_family) add(REG_SENT2$S3_exploratory) else if (isTRUE(R$S3_labels$S3a_supported)) {
    add(REG_SENT2$S3_coupled); if (isTRUE(R$S3_labels$S3b_specific)) add(REG_SENT2$S3_specific) } else add(REG_SENT2$S3_null)
  if (any(unlist(R$holm$opposite_direction_significant))) add(REG_SENT2$opposite)
  if (length(R$branch$caveats)) add(paste0("Caveats (registered branch outcomes): ", paste(R$branch$caveats, collapse = "; "), "."))
  out
}
