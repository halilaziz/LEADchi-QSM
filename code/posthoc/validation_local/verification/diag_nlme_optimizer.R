#!/usr/bin/env Rscript
# DIAGNOSTIC ONLY (POST HOC verification): adjudicates the small lme discrepancies in D_age_slopes.csv.
# Fits one model three ways with nlme: lmeControl(opt = "optim") as in validation_local.R, the nlme default (nlminb),
# and nlminb with tight tolerances. Reads the merged table; writes nothing (prints to stdout only).
suppressPackageStartupMessages(library(nlme))
d <- read.csv("/path/to/Penn_LEAD/LEADchi/mirror/unlocked/merged_qsmxt_primary_SNVTA.csv", stringsAsFactors = FALSE)
u <- d[d$usable == 1, ]
u$PS_e <- ifelse(u$study_group == "PRO/CHR", .5, -.5); u$ADHD_e <- ifelse(u$dx_adhd == 1, .5, -.5); u$sex_e <- ifelse(u$sex == "M", .5, -.5)
for (f in list(chi_PUT ~ age + sex_e + PS_e * ADHD_e, chi_PUT ~ age + sex_e)) {
  fits <- list(optim = lme(f, random = ~1 | participant_id, data = u, method = "REML", control = lmeControl(opt = "optim")),
               nlminb_default = lme(f, random = ~1 | participant_id, data = u, method = "REML"),
               nlminb_tight = lme(f, random = ~1 | participant_id, data = u, method = "REML",
                                  control = lmeControl(msTol = 1e-14, tolerance = 1e-12, niterEM = 0, msMaxIter = 1000)))
  cat(deparse(f), "\n")
  for (nm in names(fits)) { m <- fits[[nm]]; tt <- summary(m)$tTable
    cat(sprintf("  %-15s age %.8f SE %.8f p %.6f REML logLik %.10f\n", nm, tt["age", "Value"], tt["age", "Std.Error"], tt["age", "p-value"], logLik(m))) }
}
