suppressPackageStartupMessages({ library(splines); library(sandwich); library(MASS); library(jsonlite); library(digest) })
ROOT <- "/path/to/Penn_LEAD/LEADchi"
CODE <- file.path(ROOT, "mirror/scripts/analysis")
GATES <- file.path(ROOT, "mirror/gates/gates_preunlock.json")
FROZEN_BOUND <- -6.13180464348368
for (f in c("constants.R", "utils.R", "design.R", "plausibility.R", "m1.R", "sensitivity.R")) source(file.path(CODE, "R", f))
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
d0 <- utils::read.csv(file.path(ROOT,"results_posthoc/humanqc_primary/input/merged_H_all180.csv"), stringsAsFactors = FALSE, na.strings = c("NA", ""))
validate_table(d0, read_schema(file.path(CODE, "schema/roi_table_schema.csv")))
d0 <- code_vars(d0)
gates <- jsonlite::read_json(GATES, simplifyVector = TRUE); rL <- gates$rL_primary; pc4 <- isTRUE(gates$PC4)
d <- apply_plausibility(d0)$d
S <- build_samples(d)
first <- set_race3(S$first); first$cell <- factor(first$cell, levels = CELLS)
fixedB <- function(est, se, df, bound = FROZEN_BOUND) c(bound=bound, tB=(est-bound)/se, pB=stats::pt((est - bound) / se, df, lower.tail = FALSE))
rep1 <- function(dd, lab){ dd <- set_race3(dd); core <- run_m1_core(dd, rL, pc4); tt <- term_test(core$m, "PS_e")
  cat(sprintf("%s: n=%d df=%d beta=%.4f se=%.4f CI[%.3f,%.3f] p=%.4f | frozen bound pB=%.6f tB=%.3f | recomputed bound %.3f pB=%.6f | outcome %s\n",
   lab, core$m$n, core$m$df, tt$est, tt$se, tt$ci95[1], tt$ci95[2], tt$p_two, fixedB(tt$est,tt$se,core$m$df)["pB"], fixedB(tt$est,tt$se,core$m$df)["tB"],
   core$dec$B_registered$bound_ppb, core$dec$B_registered$pB, core$dec$outcome)) }
rep1(first, "full")
o <- first[order(first$chi_SNVTA), c("participant_id","session_id","cell","chi_SNVTA")][1:5,]; print(o)
rep1(first[first$participant_id != "sub-22045", ], "drop sub-22045 youth")
# alternative: drop the session and let the youth's next usable session enter
a <- S$usable; cat("sub-22045 usable sessions:\n"); print(a[a$participant_id=="sub-22045", c("session_id","chi_SNVTA")])
alt <- first; alt <- alt[alt$participant_id != "sub-22045",]
row2 <- set_race3(S$usable[S$usable$participant_id=="sub-22045" & S$usable$session_id=="ses-2",])
cm <- intersect(names(alt), names(row2)); alt2 <- rbind(alt[,cm], row2[,cm])
rep1(alt2, "swap sub-22045 to ses-2")
