#!/usr/bin/env Rscript
# Data helper for figs_pcn.py (PCN submission figures). READ-ONLY: sources the frozen helpers in
# mirror/scripts/analysis/R unchanged, rebuilds the first-usable sample exactly as tier4_exploratory.R and
# t3a_tier3_estimation.R do (validate -> code -> plausibility -> build_samples), refits the E-z TD-only model
# chi_SNVTA ~ ns(age, 3) + sex_e (as in t3a), and prints one JSON object to STDOUT. Writes no file.
suppressPackageStartupMessages({ library(splines); library(sandwich); library(jsonlite); library(digest) })
Sys.setenv(LANGUAGE = "en")
ROOT  <- "/path/to/Penn_LEAD/LEADchi"
CODE  <- file.path(ROOT, "mirror/scripts/analysis")
INPUT <- file.path(ROOT, "mirror/unlocked/merged_qsmxt_primary_SNVTA.csv")
for (f in c("constants.R", "utils.R", "design.R", "plausibility.R", "m1.R", "t1.R")) source(file.path(CODE, "R", f))
d <- utils::read.csv(INPUT, stringsAsFactors = FALSE, na.strings = c("NA", ""))
schema <- read_schema(file.path(CODE, "schema/roi_table_schema.csv")); validate_table(d, schema)
d <- code_vars(d); d <- apply_plausibility(d)$d
S <- build_samples(d); fu <- S$first
stopifnot(nrow(fu) == 80)
# E-z TD-only model (identical to t3a_tier3_estimation.R)
td <- fu[fu$cell == "TD", ]
ns_td <- ns(td$age, df = 3); K <- attr(ns_td, "knots"); Bk <- attr(ns_td, "Boundary.knots")
m_td <- stats::lm(chi_SNVTA ~ ns(age, knots = K, Boundary.knots = Bk) + sex_e, data = td)
ag <- seq(Bk[1], Bk[2], length.out = 100)
mu <- as.numeric(stats::predict(m_td, newdata = data.frame(age = ag, sex_e = 0)))
nonTD <- fu[fu$cell != "TD", ]
out <- list(
  input_sha256 = sha256_file(INPUT),
  first = fu[, c("sub", "session_id", "cell", "age", "sex", "dx_psychosis", "chi_SNVTA", "r2s_SNVTA")],
  ez = list(knots = K, boundary = Bk, sigma_td = stats::sigma(m_td), age_grid = ag, mu_grid = mu,
            check_nonTD = data.frame(sub = nonTD$sub, mu_TD = as.numeric(stats::predict(m_td, newdata = nonTD)))))
cat(jsonlite::toJSON(out, digits = NA, auto_unbox = TRUE, dataframe = "columns"))
