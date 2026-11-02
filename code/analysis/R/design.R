# LEAD-chi analysis: schema validation, variable coding (Q11, sec. 12.5) and analysis samples (Q13).

read_schema <- function(path) utils::read.csv(path, stringsAsFactors = FALSE)

validate_table <- function(d, schema) {
  req <- schema$column[schema$required %in% c("yes", "s3")]
  miss <- setdiff(req, names(d))
  if (length(miss)) stop("ROI table is missing required columns: ", paste(miss, collapse = ", "))
  for (v in c("data_origin", "engine", "pipeline", "roi_def")) {
    u <- unique(d[[v]]); if (length(u) != 1) stop("column ", v, " must hold a single value; found: ", paste(u, collapse = ", "))
  }
  if (!unique(d$data_origin) %in% c("real", "synthetic")) stop("data_origin must be 'real' or 'synthetic'")
  key <- paste(d$participant_id, d$session_id)
  if (anyDuplicated(key)) stop("duplicate (participant_id, session_id) rows: ", paste(head(key[duplicated(key)]), collapse = "; "))
  if (!all(d$study_group %in% c("PRO/CHR", "ADHD", "TD/NC"))) stop("unexpected study_group values")
  if (!all(d$dx_adhd %in% c(0, 1))) stop("dx_adhd must be 0/1")
  if (!all(d$sex %in% c("M", "F"))) stop("sex must be M/F")
  if (!all(d$usable %in% c(0, 1))) stop("usable must be 0/1")
  if (!all(d$s3_eligible %in% c(0, 1))) stop("s3_eligible must be 0/1")
  # sec. 12.2: the secondary reference is lateral-ventricle CSF only; no valid CSF -> 'none' and a missing value
  u <- d[d$usable == 1, ]
  if (!all(u$ref2_type %in% c("CSF", "CSF_erode1", "none"))) stop("ref2_type must be CSF / CSF_erode1 / none")
  if (any(u$ref2_type == "none" & !is.na(u$chi_SNVTA_ref2))) stop("ref2_type 'none' requires a missing chi_SNVTA_ref2")
  if (any(is.na(d$age))) stop("age missing in some rows")
  for (v in c("chi_SNVTA", "roi_vol", "Q_phase"))
    if (any(is.na(u[[v]]))) stop("usable rows must have non-missing ", v)
  o <- d[order(d$participant_id, d$ses_order), ]
  bad <- tapply(o$age, o$participant_id, function(a) any(diff(a) <= 0))
  if (any(bad, na.rm = TRUE)) stop("ses_order not chronological for: ", paste(names(bad)[which(bad)], collapse = ", "))
  types <- setNames(schema$type, schema$column)
  for (v in intersect(names(types)[types %in% c("numeric", "integer")], names(d)))
    if (!is.numeric(d[[v]]) && !all(is.na(d[[v]]))) stop("column ", v, " must be numeric")
  invisible(TRUE)
}

# Registered coding (Q11; sec. 12.5).
code_vars <- function(d) {
  d$PS_e   <- ifelse(d$study_group == "PRO/CHR", 0.5, -0.5)
  d$ADHD_e <- ifelse(d$dx_adhd == 1, 0.5, -0.5)
  d$sex_e  <- ifelse(d$sex == "M", 0.5, -0.5)
  d$race3_chr <- ifelse(d$race == "Black or African American", "Black", ifelse(d$race == "White", "White", "Other"))
  d$cell <- ifelse(d$PS_e > 0, ifelse(d$ADHD_e > 0, "PS+ADHD", "PS-only"), ifelse(d$ADHD_e > 0, "ADHD-only", "TD"))
  d$sub <- d$participant_id
  d$sess_key <- paste0(d$participant_id, "_", d$session_id)
  d
}

# race3 as a deviation-coded factor over the levels present in the analysed sample (contr.sum).
set_race3 <- function(d) {
  lev <- intersect(c("Black", "White", "Other"), unique(d$race3_chr))
  d$race3 <- factor(d$race3_chr, levels = lev)
  if (length(lev) >= 2) stats::contrasts(d$race3) <- stats::contr.sum(length(lev))
  d
}
race_term <- function(d) if (nlevels(d$race3) >= 2) " + race3" else ""

# Samples (Q13). First usable = earliest usable session; occasion 1 = the earliest usable session after it.
# Exposure sensitivity (Q18(4)): the four listed session rows are dropped, never substituted by a later session.
build_samples <- function(d) {
  d <- d[order(d$participant_id, d$ses_order), ]
  u <- d[d$usable == 1, ]
  u$k_usable <- stats::ave(u$ses_order, u$participant_id, FUN = seq_along)
  first <- u[u$k_usable == 1, ]
  if ("first_usable" %in% names(d)) {
    fl <- d$sess_key[d$first_usable %in% 1]
    if (!setequal(fl, first$sess_key))
      stop("frozen first_usable flags disagree with the earliest usable session for: ",
           paste(head(setdiff(union(fl, first$sess_key), intersect(fl, first$sess_key)), 10), collapse = ", "))
  }
  occ1 <- u[u$k_usable == 2, ]
  first$occ <- 0; occ1$occ <- 1
  cn <- intersect(names(first), names(occ1))
  two_occ <- rbind(first[, cn], occ1[, cn])
  # S2: first usable session only, age >= 11 at that session, PRIME/BDI/ARI present; no substitution (Q13)
  s2 <- first[first$age >= REG$s2_age_min & !is.na(first$PRIME_total) & !is.na(first$BDI_total) & !is.na(first$ARI_total), ]
  s3_vars <- c("twICVF_pSCP_R", "b0loc_pSCP_R", "tract_len", "nSL", "FD", "NDC", "chi_DN_R", "chi_DGM",
               "DN_R_vol", "abs_pitch")
  s3_all <- first[first$s3_eligible == 1, ]
  s3 <- s3_all[stats::complete.cases(s3_all[, s3_vars]), ]
  list(all = d, usable = u, first = first, occ1 = occ1, two_occ = two_occ, s2 = s2, s3 = s3,
       s3_dropped_incomplete = setdiff(s3_all$sess_key, s3$sess_key),
       ses1_only = first[first$session_id == "ses-1", ],
       exposure = first[!first$sess_key %in% REG$exposure_sessions, ],
       s3_exposure = s3[!s3$sess_key %in% REG$exposure_sessions, ])
}

sample_table <- function(S) {
  rows <- list(
    `All MEGRE sessions in table` = S$all, `Usable sessions` = S$usable, `H1 / S1 first-usable` = S$first,
    `S2 (age >= 11, PRIME+BDI+ARI)` = S$s2, `S3 (same-session DWI, >= n_min nodes)` = S$s3,
    `T1 occasion 1 (held-out follow-ups)` = S$occ1, `ses-1-only first-usable` = S$ses1_only,
    `Exposure sensitivity, H1 (4 inspected sessions dropped)` = S$exposure,
    `Exposure sensitivity, S3` = S$s3_exposure)
  do.call(rbind, lapply(names(rows), function(n) data.frame(sample = n, t(cell_counts(rows[[n]])), check.names = FALSE)))
}
