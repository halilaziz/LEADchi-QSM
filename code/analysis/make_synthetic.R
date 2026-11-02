#!/usr/bin/env Rscript
# Write SYNTHETIC ROI tables (schema/roi_table_schema.csv) on the real Penn LEAD design, one per scenario, plus a
# synthetic gate JSON for each. Outcomes are simulated; design, availability and phenotype scores are real metadata.
#   Rscript make_synthetic.R [--out synthetic] [--seed 20260924]
all_args <- commandArgs(trailingOnly = FALSE)
CODE_DIR <- dirname(normalizePath(sub("^--file=", "", all_args[grep("^--file=", all_args)][1])))
suppressPackageStartupMessages(library(jsonlite))
for (f in c("constants.R", "utils.R", "synthetic.R")) source(file.path(CODE_DIR, "R", f))
opt <- parse_args(commandArgs(trailingOnly = TRUE))
out <- if (is.null(opt$out)) file.path(CODE_DIR, "synthetic") else opt$out
seed <- if (is.null(opt$seed)) REG$seed else as.integer(opt$seed)
dir.create(out, showWarnings = FALSE, recursive = TRUE)
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
des <- build_real_design()
cat(sprintf("real design: %d MEGRE sessions, %d subjects, %d T1-eligible sessions\n", nrow(des),
            length(unique(des$participant_id)), sum(des$t1_ok)))
for (sc in names(SCENARIOS)) {
  tab <- simulate_table(des, sc, seed)
  fn <- file.path(out, sprintf("SYNTHETIC_roi_table_%s.csv", sc))
  utils::write.csv(tab, fn, row.names = FALSE, na = "NA")
  jsonlite::write_json(synthetic_gates(tab, sc), file.path(out, sprintf("SYNTHETIC_gates_%s.json", sc)),
                       auto_unbox = TRUE, pretty = TRUE, digits = NA)
  utils::write.csv(simulate_multiverse(tab, seed), file.path(out, sprintf("SYNTHETIC_multiverse_%s.csv", sc)), row.names = FALSE, na = "NA")
  cat(sprintf("%-14s rows %d, usable %d, first-usable %d -> %s\n", sc, nrow(tab), sum(tab$usable), sum(tab$first_usable), fn))
}
