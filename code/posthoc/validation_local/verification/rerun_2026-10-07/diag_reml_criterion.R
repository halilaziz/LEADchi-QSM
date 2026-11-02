#!/usr/bin/env Rscript
# DIAGNOSTIC ONLY (POST HOC, 2026-10-07): profiled REML criterion at nlme's default solution vs at the verifier's theta
# (verification/mine/D_age_slopes.csv) for three lme models. Lower = closer to the optimum. Reads only; prints to stdout.
suppressPackageStartupMessages({library(splines); library(sandwich); library(MASS); library(nlme); library(jsonlite); library(digest)})
ROOT <- "/path/to/Penn_LEAD/LEADchi"; CODE <- file.path(ROOT, "mirror/scripts/analysis")
for (f in c("constants.R", "utils.R", "design.R", "plausibility.R", "m1.R")) source(file.path(CODE, "R", f))
d0 <- utils::read.csv(file.path(ROOT, "mirror/unlocked/merged_qsmxt_primary_SNVTA.csv"), stringsAsFactors = FALSE, na.strings = c("NA", ""))
d0 <- code_vars(d0); pl <- apply_plausibility(d0); S <- build_samples(pl$d); u <- S$usable
u$chi_GP_csfref <- u$chi_GP - u$ref2_minus_wb_ppb
mine <- read.csv(file.path(ROOT, "results_posthoc/validation_local/verification/mine/D_age_slopes.csv"), stringsAsFactors = FALSE)
nll <- function(theta, X, y, g) {   # profiled REML criterion (same algebra as the verifier)
  N <- nrow(X); p <- ncol(X); ng <- as.vector(table(g)); Sx <- rowsum(X, g); Sy <- as.vector(rowsum(y, g))
  c_ <- theta / (1 + ng * theta)
  XHX <- crossprod(X) - t(Sx * c_) %*% Sx; XHy <- crossprod(X, y) - colSums(Sx * (c_ * Sy)); yHy <- sum(y^2) - sum(c_ * Sy^2)
  b <- solve(XHX, XHy); s2 <- drop(yHy - crossprod(XHy, b)) / (N - p)
  0.5 * ((N - p) * log(s2) + sum(log1p(ng * theta)) + as.numeric(determinant(XHX)$modulus))
}
for (case in list(list(y = "chi_GP_csfref", f = "chi_GP_csfref ~ age + sex_e", cv = "+ sex"),
                  list(y = "chi_GP", f = "chi_GP ~ age + sex_e", cv = "+ sex"),
                  list(y = "chi_PUT", f = "chi_PUT ~ age + sex_e + PS_e * ADHD_e", cv = "+ sex + cell"))) {
  du <- u[!is.na(u[[case$y]]), ]; du <- du[order(du$sub), ]
  m <- lme(as.formula(case$f), random = ~1|sub, data = du, method = "REML")
  vc <- as.numeric(VarCorr(m)[, "Variance"]); th_nlme <- vc[1] / vc[2]
  X <- model.matrix(as.formula(case$f), du); y <- du[[case$y]]; g <- factor(du$sub, levels = unique(du$sub))
  th_ver <- mine$theta[mine$outcome == case$y & mine$sample == "all usable sessions (lme REML, random intercept)" & mine$covariates == case$cv]
  opt <- optimize(nll, c(1e-6, 50), X = X, y = y, g = g, tol = 1e-14)
  cat(sprintf("%-40s theta nlme %.10f verifier %.10f R-optimize %.10f | nll nlme %.12f verifier %.12f\n", case$f, th_nlme, th_ver, opt$minimum,
              nll(th_nlme, X, y, g), nll(th_ver, X, y, g)))
}
