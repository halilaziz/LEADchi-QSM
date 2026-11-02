#!/usr/bin/env Rscript
# LEAD-chi human visual-QC validation -- POST HOC refit of M1 on re-derived first-usable samples.
# NOT part of the registered confirmatory analysis; nothing here re-decides H1. Called by analyze_human_qc.py.
#
#   Rscript refit_m1_human_qc.R --merged <unlocked merged table> --scen <scenario_usable.csv> --gates <gates json>
#                               --rdir <frozen R code dir> --out <json> [--bound -6.13180464348368]
#
# The frozen helpers (constants.R, utils.R, design.R, plausibility.R, m1.R) are SOURCED read-only from the frozen
# mirror, so the design matrix, coding (effect codes, contr.sum race3, ns(age, df = 3)), OLS + HC3 (t, df = n - p),
# Huber sigma_rob and the decision-rule arithmetic are exactly the registered ones. The only change per scenario is
# the `usable` / `first_usable` columns (re-derived by analyze_human_qc.py with the registered first-usable rule).
# The same steps as run_all.R are applied: validate_table -> code_vars -> apply_plausibility -> build_samples ->
# run_m1_core(first, rL_primary, PC4).
suppressPackageStartupMessages({ library(splines); library(sandwich); library(MASS); library(jsonlite) })
Sys.setenv(LANGUAGE = "en"); invisible(Sys.setlocale("LC_MESSAGES", "C"))
args <- commandArgs(trailingOnly = TRUE)
opt <- list(); i <- 1
while (i <= length(args)) { opt[[sub("^--", "", args[i])]] <- args[i + 1]; i <- i + 2 }
for (k in c("merged", "scen", "gates", "rdir", "out")) if (is.null(opt[[k]])) stop("missing --", k)
BOUND_FIXED <- if (is.null(opt$bound)) -6.13180464348368 else as.numeric(opt$bound)

for (f in c("constants.R", "utils.R", "design.R", "plausibility.R", "m1.R")) source(file.path(opt$rdir, "R", f))
schema <- read_schema(file.path(opt$rdir, "schema", "roi_table_schema.csv"))
RNGkind("Mersenne-Twister", "Inversion", "Rejection")

d0 <- utils::read.csv(opt$merged, stringsAsFactors = FALSE, na.strings = c("NA", ""))
sc <- utils::read.csv(opt$scen, stringsAsFactors = FALSE, na.strings = c("NA", ""))
gates <- jsonlite::read_json(opt$gates, simplifyVector = TRUE)
rL <- gates$rL_primary; pc4 <- isTRUE(gates$PC4)
key0 <- paste0(d0$participant_id, "_", d0$session_id)
m <- match(key0, sc$sess_key)
if (anyNA(m) || nrow(sc) != nrow(d0)) stop("scenario table does not cover the merged table one-to-one")
scen_names <- sub("^usable__", "", grep("^usable__", names(sc), value = TRUE))

fit_one <- function(s) {
  d <- d0
  d$usable <- as.integer(sc[[paste0("usable__", s)]][m])
  d$first_usable <- as.integer(sc[[paste0("first_usable__", s)]][m])
  validate_table(d, schema)
  d <- code_vars(d)
  d <- apply_plausibility(d)$d
  S <- build_samples(d)                      # stops if first_usable disagrees with the earliest usable session
  core <- run_m1_core(S$first, rL, pc4)
  tt <- term_test(core$m, "PS_e"); dec <- core$dec
  tB <- (tt$est - BOUND_FIXED) / tt$se
  ct <- coef_table(core$m)
  list(scenario = s, n = core$m$n, df = core$m$df, cells = as.list(cell_counts(S$first)),
       n_usable_sessions = nrow(S$usable), first_usable_keys = sort(S$first$sess_key),
       beta_PS_ppb = tt$est, se_hc3 = tt$se, t = tt$t, p_two = tt$p_two, ci95_ppb = tt$ci95, ci90_ppb = tt$ci90,
       sigma_resid = core$m$sigma, sigma_rob = core$huber$scale, huber_converged = core$huber$converged,
       sigma_B = dec$sigma_B, delta = dec$delta,
       bound_fixed = list(bound_ppb = BOUND_FIXED, note = "registered bound frozen at the confirmatory run (-6.1318 ppb)",
                          tB = tB, df = tt$df, pB_one_sided = stats::pt(tB, tt$df, lower.tail = FALSE),
                          reject_at_05 = stats::pt(tB, tt$df, lower.tail = FALSE) < REG$alpha),
       bound_rederived = c(dec$B_registered, list(note = "delta_B = 0.6*sqrt(rL_primary) times sigma_B recomputed in this sample")),
       A_beta_negative_p05 = dec$A, B_rederived = dec$B, outcome_label_descriptive = dec$outcome,
       disconfirming = dec$disconfirming, formula = paste(deparse(core$formula), collapse = " "),
       coefficients = ct)
}
res <- lapply(scen_names, function(s) tryCatch(fit_one(s), error = function(e) list(scenario = s, error = conditionMessage(e))))
names(res) <- scen_names
out <- list(label = "POST HOC sensitivity (human visual QC replaces/complements the AI grade); not confirmatory",
            R = R.version.string, sandwich = as.character(utils::packageVersion("sandwich")),
            MASS = as.character(utils::packageVersion("MASS")), rL_primary = rL, PC4 = pc4,
            frozen_code_sha256 = vapply(c("constants.R", "utils.R", "design.R", "plausibility.R", "m1.R"),
                                        function(f) tryCatch(sha256_file(file.path(opt$rdir, "R", f)),
                                                             error = function(e) NA_character_), ""),
            scenarios = res)
jsonlite::write_json(out, opt$out, auto_unbox = TRUE, digits = NA, pretty = TRUE, na = "null")
cat("wrote", opt$out, "\n")
