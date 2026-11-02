#!/usr/bin/env Rscript
# LEAD-chi confirmatory analysis -- single entry point (registered attachment; frozen and hashed at posting).
#
#   Rscript run_all.R --input <merged_roi_table.csv> --gates <gates_preunlock.json> --multiverse <multiverse.csv>
#                     --out <output_dir> [--other-engine <merged table of the other engine>] [--label <text>]
#
# Reads ONLY the merged session-level ROI table and the long multiverse table written by unlock_and_merge.py
# (schemas: schema/roi_table_schema.csv, schema/multiverse_schema.csv) and the aggregate gate JSON written by
# gates_preunlock.py. Registered order (Annex A.6 step 8): QC-by-group -> PC2x -> H1 (M1, A/B/C, K1-K3, co-reports,
# influence, fragility) -> S1-S3 + Holm -> T1 -> robustness (sensitivity battery, specification curve / K4 "robust",
# pooled all-session model, S2 all-session model). The registered text governs; Annex H maps rules to functions.
# --n-perm / --n-perm-spec / --n-boot / --n-pb / --n-wild override resampling counts for TESTING ONLY; the report
# and results.json then carry a "NON-REGISTERED counts" flag. Nothing here reads imaging data.

t_start <- Sys.time()
Sys.setenv(LANGUAGE = "en"); invisible(Sys.setlocale("LC_MESSAGES", "C"))
all_args <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("^--file=", "", all_args[grep("^--file=", all_args)][1]))
CODE_DIR <- dirname(script_path)
suppressPackageStartupMessages({
  library(splines); library(sandwich); library(lmtest); library(lme4); library(nlme); library(MASS); library(boot)
  library(jsonlite); library(digest)
})
CODE_FILES <- c(file.path("R", c("constants.R", "utils.R", "design.R", "plausibility.R", "m1.R", "secondary.R", "t1.R",
                                 "sensitivity.R", "spec_curve.R", "report.R")),
                "schema/roi_table_schema.csv", "schema/multiverse_schema.csv", "run_all.R")
for (f in CODE_FILES[grepl("\\.R$", CODE_FILES) & CODE_FILES != "run_all.R"]) source(file.path(CODE_DIR, f))

opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(opt$input) || is.null(opt$gates) || is.null(opt$out))
  stop("usage: Rscript run_all.R --input <roi_table.csv> --gates <gates.json> --multiverse <mv.csv> --out <dir>")
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
draws <- list(n_perm = REG$n_perm, n_perm_spec = REG$n_perm_spec, n_boot_delta = REG$n_boot_delta,
              n_boot_s3b = REG$n_boot_s3b, n_pb = REG$n_pb, n_wild = REG$n_wild)
if (!is.null(opt$n_perm)) draws$n_perm <- as.integer(opt$n_perm)
if (!is.null(opt$n_perm_spec)) draws$n_perm_spec <- as.integer(opt$n_perm_spec)
if (!is.null(opt$n_boot)) { draws$n_boot_delta <- as.integer(opt$n_boot); draws$n_boot_s3b <- as.integer(opt$n_boot) }
if (!is.null(opt$n_pb)) draws$n_pb <- as.integer(opt$n_pb)
if (!is.null(opt$n_wild)) draws$n_wild <- as.integer(opt$n_wild)
draws_registered <- all(unlist(draws) == c(REG$n_perm, REG$n_perm_spec, REG$n_boot_delta, REG$n_boot_s3b, REG$n_pb, REG$n_wild))
REG$n_wild <- draws$n_wild

log_msg("LEAD-chi run_all.R start; input ", opt$input)
d <- utils::read.csv(opt$input, stringsAsFactors = FALSE, na.strings = c("NA", ""))
schema <- read_schema(file.path(CODE_DIR, "schema/roi_table_schema.csv"))
validate_table(d, schema)
d <- code_vars(d)
origin <- unique(d$data_origin)

# ---- environment pin (Q23): real data -> any version difference stops the run; synthetic -> recorded ----
inst <- c(R = paste(R.version$major, R.version$minor, sep = "."),
          vapply(setdiff(names(REG$pinned), "R"), function(p) as.character(utils::packageVersion(p)), ""))
inst <- gsub("-", ".", inst); pinned <- unlist(REG$pinned)
vdiff <- names(pinned)[pinned != inst[names(pinned)]]
version_check <- if (!length(vdiff)) "all pinned versions match" else paste0("DIFFERS: ", paste0(vdiff, " ", inst[vdiff], " (pinned ", pinned[vdiff], ")", collapse = "; "))
if (length(vdiff) && origin == "real") stop("environment differs from the registered pin: ", version_check)

# ---- gate JSON validation (Annex A.6; Q23): missing or invalid decision fields -> hard stop ----
gates <- jsonlite::read_json(opt$gates, simplifyVector = TRUE)
gmiss <- setdiff(REG$gate_required, names(gates))
if (length(gmiss)) stop("gate JSON lacks required fields: ", paste(gmiss, collapse = ", "))
if (!identical(gates$data_origin, origin)) stop("gates data_origin (", gates$data_origin, ") != table data_origin")
for (k in c("PC1", "PC3", "PC4", "PC4_DN", "PC5", "G2", "G2_DN", "G3", "G4", "G5", "A2_approved"))
  if (!is.logical(gates[[k]]) || length(gates[[k]]) != 1 || is.na(gates[[k]])) stop("gate field ", k, " must be TRUE or FALSE")
if (!gates$PC2w %in% c("PASS", "INDETERMINATE", "FAIL")) stop("gate field PC2w must be PASS / INDETERMINATE / FAIL")
chk_num <- function(x, lo, hi, nm) if (!is.numeric(x) || length(x) != 1 || !is.finite(x) || x < lo || x > hi)
  stop("gate field ", nm, " must be a finite number in [", lo, ", ", hi, "]; found: ", paste(x, collapse = ","))
chk_num(gates$rL_primary, -1, 1, "rL_primary"); chk_num(gates$r_point_primary, -1, 1, "r_point_primary")
for (k in c("mean_SNVTA", "mean_STR")) chk_num(gates$pooled[[k]], -1e4, 1e4, paste0("pooled.", k))
for (k in c("sd_SNVTA", "sd_STR", "sd_total_SNVTA")) chk_num(gates$pooled[[k]], 1e-9, 1e4, paste0("pooled.", k))

# ---- branch from the gate JSON (Q24 branch order) ----
pc_fail <- !isTRUE(gates$PC1) || identical(gates$PC2w, "FAIL") || !isTRUE(gates$PC3)
branch <- list(engine = if (isTRUE(gates$A2_approved) && !isTRUE(gates$G2) && !pc_fail) "sepia" else "qsmxt",   # fallback = QSMxT PDF+iLSQR
               pipeline = if (!pc_fail) "primary" else "fallback_PDF_iLSQR",
               roi = if (isTRUE(gates$G3)) "SNVTA" else "SN_only", caveats = character(0))
if (pc_fail && is.null(gates$fallback_PC_pass)) stop("gate JSON lacks fallback_PC_pass although a positive control failed")
branch$primary_interpretable <- (!pc_fail || isTRUE(gates$fallback_PC_pass)) && isTRUE(gates$G5)
for (k in c("engine", "pipeline")) if (unique(d[[k]]) != branch[[k]])
  stop("table ", k, " = ", unique(d[[k]]), " but the gates select ", branch[[k]], ": wrong input table for this branch")
if (unique(d$roi_def) != branch$roi) stop("table roi_def = ", unique(d$roi_def), " but G3 selects ", branch$roi)
if (!identical(gates$constants_roi, branch$roi)) stop("gate constants were computed for ", gates$constants_roi, " but G3 selects ", branch$roi)
if (!branch$primary_interpretable) branch$caveats <- c(branch$caveats, "pipeline not validated (positive-control or registration gate failure): methods-report path")
if (!isTRUE(gates$A2_approved) && !isTRUE(gates$G2)) branch$caveats <- c(branch$caveats, "inversion-algorithm agreement (within-QSMxT G2) not established")
if (!isTRUE(gates$G3) && identical(gates$G3_SNonly, FALSE)) branch$caveats <- c(branch$caveats, "method-dependent ROI measure (G3 failed for SN-VTA and for the SN-only fallback)")
if (!isTRUE(gates$G4)) branch$caveats <- c(branch$caveats, "QC reliability not established (G4 intra-rater kappa < 0.60 after the registered re-rating)")
log_msg("branch: engine=", branch$engine, " pipeline=", branch$pipeline, " roi=", branch$roi)

# ---- plausibility re-verification and samples ----
pl <- apply_plausibility(d); d <- pl$d
S <- build_samples(d)
log_msg("samples: first-usable ", cell_str(S$first), "; S2 ", cell_str(S$s2), "; S3 ", cell_str(S$s3), "; T1 occ1 ", cell_str(S$occ1))
seed <- REG$seed

# ---- 1. QC-by-group and PC2x ----
qc <- qc_by_group(S, gates$rL_primary, isTRUE(gates$PC4))
PC2x <- pc2x(S$first)
if (identical(gates$PC2w, "INDETERMINATE")) {
  branch$caveats <- c(branch$caveats, "PC2w indeterminate: PC2x evaluated after unlock")
  if (!isTRUE(PC2x$pass)) branch$caveats <- c(branch$caveats, "developmental positive control not confirmed (applies to every H1 statement)")
}

# ---- 2. H1 and K1-K3 ----
log_msg("H1")
H1 <- run_h1(S, gates, draws$n_perm, draws$n_boot_delta, seed)

# ---- 3. S1-S3 and Holm ----
log_msg("S1-S3")
S1 <- run_s1(S, gates); S2 <- run_s2(S); S3 <- run_s3(S, gates, draws$n_perm, draws$n_boot_s3b, seed)
holm_res <- holm_family(S1, S2, S3)
S1_labels <- s1_labels(S1, holm_res$reject$S1)
S3_labels <- s3_labels(S3, if (holm_res$S3_in_family) holm_res$reject$S3a else FALSE)
sec_infl <- secondary_influence(S1, S2, S3, holm_res, gates)

# ---- 4. T1 ----
log_msg("T1")
T1 <- run_t1(S, H1, gates, draws$n_pb, seed)

# ---- 5. robustness ----
log_msg("robustness: sensitivity battery")
other <- if (!is.null(opt$other_engine)) code_vars(utils::read.csv(opt$other_engine, stringsAsFactors = FALSE, na.strings = c("NA", ""))) else NULL
SENS <- run_h1_sensitivities(S, gates, seed, other)
log_msg("robustness: specification curve")
if (is.null(opt$multiverse)) {
  if (origin == "real") stop("--multiverse is required for real data (the K4 'robust' label depends on it)")
  SC <- list(note = "not run: no multiverse table supplied (synthetic run)", robust = FALSE)
} else {
  mv <- utils::read.csv(opt$multiverse, stringsAsFactors = FALSE, na.strings = c("NA", ""))
  SC <- run_spec_curve(S, mv, H1, draws$n_perm_spec, seed + 41, branch$roi)
}
log_msg("robustness: pooled all-session and S2 all-session models")
pooled <- run_pooled(S, draws$n_wild, seed)
S2_all <- run_s2_allsession(S)

strip <- function(x) { x$model <- NULL; x$data <- NULL; x }
R <- list(meta = list(code_version = REG$code_version, data_origin = origin, scenario = if (is.null(opt$label)) "" else opt$label,
                      input = normalizePath(opt$input), input_sha256 = sha256_file(opt$input),
                      gates = normalizePath(opt$gates), gates_sha256 = sha256_file(opt$gates),
                      multiverse = if (!is.null(opt$multiverse)) normalizePath(opt$multiverse) else NULL,
                      multiverse_sha256 = if (!is.null(opt$multiverse)) sha256_file(opt$multiverse) else NULL,
                      started_utc = format(t_start, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"), r_version = R.version.string,
                      packages = as.list(inst), version_check = version_check,
                      seed = REG$seed, draws = draws, draws_registered = draws_registered),
          branch = branch, gates = gates, samples = sample_table(S),
          plausibility = list(log = pl$log, robust_z = robust_z_counts(S$first)),
          qc = qc, PC2x = PC2x, H1 = H1, S1 = strip(S1), S1_labels = S1_labels, S2 = strip(S2), S3 = strip(S3), S3_labels = S3_labels,
          holm = holm_res, sec_influence = sec_infl, T1 = T1, sensitivity = SENS, spec_curve = SC, pooled = pooled,
          S2_allsession = S2_all)
R$reading <- registered_reading(R)
# Pre-registration (synthetic) outputs must not display which of the four quantitatively inspected sessions fall in
# which cell or sample; after unlock (data_origin "real") they are reported in full.
if (origin == "synthetic") {
  R$H1$exposure_sensitivity$dropped <- NULL; R$T1$exposure_sensitivity$dropped <- NULL
  ix <- grepl("^Exposure sensitivity", R$samples$sample); R$samples[ix, CELLS] <- NA
  R$samples$sample[ix] <- paste0(R$samples$sample[ix], " [cells masked in synthetic runs]")
}
R$meta$finished_utc <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
R$meta$runtime_min <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))

jsonify <- function(x) {
  if (is.data.frame(x)) return(x)
  if (is.list(x)) return(lapply(x, jsonify))
  if (is.table(x)) return(as.list(x))
  if (is.atomic(x) && !is.null(names(x)) && length(x) > 1) return(as.list(x))
  if (is.atomic(x) && !is.null(names(x)) && length(x) == 1) return(unname(x))
  x
}
res_json <- file.path(opt$out, "results.json"); rep_md <- file.path(opt$out, "report.md")
jsonlite::write_json(jsonify(R), res_json, auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null", null = "null")
write_report(R, rep_md)
writeLines(capture.output(utils::sessionInfo()), file.path(opt$out, "sessionInfo.txt"))
man <- list(written_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
            inputs = list(roi_table = list(path = normalizePath(opt$input), sha256 = sha256_file(opt$input)),
                          gates = list(path = normalizePath(opt$gates), sha256 = sha256_file(opt$gates)),
                          multiverse = if (!is.null(opt$multiverse)) list(path = normalizePath(opt$multiverse), sha256 = sha256_file(opt$multiverse)) else NULL),
            code = setNames(lapply(CODE_FILES, function(f) sha256_file(file.path(CODE_DIR, f))), CODE_FILES),
            outputs = list(results.json = sha256_file(res_json), report.md = sha256_file(rep_md),
                           sessionInfo.txt = sha256_file(file.path(opt$out, "sessionInfo.txt"))))
jsonlite::write_json(man, file.path(opt$out, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
log_msg("done: H1 final outcome ", H1$final_outcome, "; runtime ", f3(R$meta$runtime_min, 2), " min")
